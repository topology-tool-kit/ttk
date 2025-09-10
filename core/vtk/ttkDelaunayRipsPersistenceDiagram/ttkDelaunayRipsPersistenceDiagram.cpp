#include <ttkDelaunayRipsPersistenceDiagram.h>
#include <ttkRipsPersistenceDiagram.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDelaunayRipsPersistenceDiagram);

ttkDelaunayRipsPersistenceDiagram::ttkDelaunayRipsPersistenceDiagram() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkDelaunayRipsPersistenceDiagram::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceDiagram::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceDiagram::RequestData(vtkInformation *ttkNotUsed(request),
                                                   vtkInformationVector **inputVector,
                                                   vtkInformationVector *outputVector) {

  ttk::Timer tm{};

  vtkTable *input = vtkTable::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector);

  if(!input)
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

  std::vector<vtkAbstractArray *> arrays;
  arrays.reserve(ScalarFields.size());
  for(const auto &s : ScalarFields)
    arrays.push_back(input->GetColumnByName(s.data()));

  const int numberOfPoints = input->GetNumberOfRows();
  const int dimension = ScalarFields.size();
  const bool doGenerators = OutputGenerators && (dimension == 2 || dimension == 3);

  PointCloud points(numberOfPoints);
  for(int i = 0; i < numberOfPoints; ++i) {
    for(int j = 0; j < dimension; ++j)
      points[i].push_back(arrays[j]->GetVariantValue(i).ToDouble());
  }
  this->printMsg(
    "Computing Delaunay-Rips persistence diagram", 1.0, tm.getElapsedTime(), 1);
  this->printMsg("#dimensions: " + std::to_string(dimension)
                   + ", #points: " + std::to_string(numberOfPoints),
                 0.0, tm.getElapsedTime(), 1);

  MultidimensionalDiagram diagram;
  std::vector<Generator> generators;

  if(this->execute(points, diagram, generators) != 0)
    return 0;

  if(doGenerators) { // todo
    /*vtkNew<vtkPoints> vtkPoints{};
    MakeVtkPoints(vtkPoints, points);
    GeneratorsToVTU(outputPersistenceDiagram, vtkPoints, generators, true);*/
  } else
    DiagramToVTU(outputPersistenceDiagram, diagram, inf);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}