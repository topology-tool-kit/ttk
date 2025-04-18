#include <ttkRipsPersistenceGenerators.h>
#include <ttkRipsPersistenceDiagram.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkRipsPersistenceGenerators);

ttkRipsPersistenceGenerators::ttkRipsPersistenceGenerators() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

int ttkRipsPersistenceGenerators::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceGenerators::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceGenerators::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  ttk::Timer tm{};

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkPointSet *pointSet = vtkPointSet::GetData(inputVector[1]);
  vtkUnstructuredGrid *outputGenerators
    = vtkUnstructuredGrid::GetData(outputVector, 0);
  vtkUnstructuredGrid *outputPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 1);

  if(!input || !pointSet)
    return 0;

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = input->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      }
    }
  }

  if(input->GetNumberOfRows() <= 0 || ScalarFields.size() <= 0) {
    this->printErr("Input matrix has invalid dimensions (rows: "
                   + std::to_string(input->GetNumberOfRows())
                   + ", columns: " + std::to_string(ScalarFields.size()) + ")");
    return 0;
  }
  else if (input->GetNumberOfRows() != pointSet->GetNumberOfPoints()) {
    this->printErr("Input and 3D representation have different");
    this->printErr("numbers of points: resp. "
                   + std::to_string(input->GetNumberOfRows())
                   + " and "
                   + std::to_string(pointSet->GetNumberOfPoints()));
    return 0;
  }

  std::vector<vtkAbstractArray *> arrays;
  arrays.reserve(ScalarFields.size());
  for(const auto &s : ScalarFields)
    arrays.push_back(input->GetColumnByName(s.data()));

  const int numberOfPoints = input->GetNumberOfRows();
  const int dimension = ScalarFields.size();
  std::vector<std::vector<double>> points(numberOfPoints);
  for(int i = 0; i < numberOfPoints; ++i) {
    for(int j = 0; j < dimension; ++j)
      points[i].push_back(arrays[j]->GetVariantValue(i).ToDouble());
  }

  this->printMsg("Computing Rips pers. generators",0.0, tm.getElapsedTime(), 1);
  this->printMsg("#dimensions: " + std::to_string(dimension)
                   + ", #points: " + std::to_string(numberOfPoints),
                 0.0, tm.getElapsedTime(), 1);

  std::vector<ttk::rpd::Diagram> diagram(0);
  std::vector<ttk::rpd::Generator> generators(0);
  this->execute(points, diagram, generators);

  GeneratorsToVTU(outputGenerators, pointSet->GetPoints(), generators, !OutputCascade);
  DiagramToVTU(outputPersistenceDiagram, diagram, SimplexMaximumDiameter);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}
