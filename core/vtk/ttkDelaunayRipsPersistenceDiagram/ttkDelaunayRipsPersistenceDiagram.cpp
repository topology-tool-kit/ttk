#include <ttkDelaunayRipsPersistenceDiagram.h>
#include <ttkRipsPersistenceDiagram.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDelaunayRipsPersistenceDiagram);

using namespace ttk::rpd;

ttkDelaunayRipsPersistenceDiagram::ttkDelaunayRipsPersistenceDiagram() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkDelaunayRipsPersistenceDiagram::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceDiagram::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceDiagram::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer tm{};

  vtkInformation *info = inputVector[0]->GetInformationObject(0);
  vtkDataObject *input = info->Get(vtkDataObject::DATA_OBJECT());
  vtkUnstructuredGrid *outputPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector);

  if(!input)
    return 0;

  PointCloud points;
  int numberOfPoints = 0;
  int dimension = 0;

  if(vtkTable *table = vtkTable::SafeDownCast(input)) {
    if(SelectFieldsWithRegexp) {
      // select all input columns whose name is matching the regexp
      ScalarFields.clear();
      const auto n = table->GetNumberOfColumns();
      for(int i = 0; i < n; ++i) {
        const auto &name = table->GetColumnName(i);
        if(std::regex_match(name, std::regex(RegexpString))) {
          ScalarFields.emplace_back(name);
        }
      }
    }

    if(table->GetNumberOfRows() <= 0 || ScalarFields.size() <= 1) {
      this->printErr("Input matrix has invalid dimensions (rows: "
                     + std::to_string(table->GetNumberOfRows()) + ", columns: "
                     + std::to_string(ScalarFields.size()) + ")");
      return 0;
    }

    std::vector<vtkAbstractArray *> arrays;
    arrays.reserve(ScalarFields.size());
    for(const auto &s : ScalarFields)
      arrays.push_back(table->GetColumnByName(s.data()));

    numberOfPoints = table->GetNumberOfRows();
    dimension = ScalarFields.size();

    points.resize(numberOfPoints);
    for(int i = 0; i < numberOfPoints; ++i) {
      for(int j = 0; j < dimension; ++j)
        points[i].push_back(arrays[j]->GetVariantValue(i).ToDouble());
    }
  }

  else if(vtkPointSet *pointset = vtkPointSet::SafeDownCast(input)) {
    numberOfPoints = pointset->GetNumberOfPoints();
    dimension = 3;
    points.resize(numberOfPoints, std::vector<double>(3));
    for(int i = 0; i < numberOfPoints; ++i)
      pointset->GetPoint(i, points[i].data());
  }

  this->printMsg("Computing Delaunay-Rips persistence diagram", 1.0,
                 tm.getElapsedTime(), getThreadNumber());
  this->printMsg("#dimensions: " + std::to_string(dimension)
                   + ", #points: " + std::to_string(numberOfPoints),
                 0.0, tm.getElapsedTime(), getThreadNumber());

  MultidimensionalDiagram diagram;
  if(this->execute(points, diagram) != 0)
    return 0;

  DiagramToVTU(outputPersistenceDiagram, diagram, inf);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), getThreadNumber());

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}