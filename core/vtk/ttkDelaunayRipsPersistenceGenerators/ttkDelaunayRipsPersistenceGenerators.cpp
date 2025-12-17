#include <ttkDelaunayRipsPersistenceGenerators.h>
#include <ttkRipsPersistenceDiagram.h>
#include <ttkRipsPersistenceGenerators.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkDelaunayRipsPersistenceGenerators);

using namespace ttk::rpd;

static void MakeVtkPoints(vtkPoints *vtkPoints,
                          const std::vector<std::vector<double>> &pointsData) {

  const int dimension = pointsData[0].size();
  vtkPoints->SetNumberOfPoints(pointsData.size());

  for(unsigned i = 0; i < pointsData.size(); ++i) {
    if(dimension >= 3)
      vtkPoints->SetPoint(
        i, pointsData[i][0], pointsData[i][1], pointsData[i][2]);
    else
      vtkPoints->SetPoint(i, pointsData[i][0], pointsData[i][1], 0.);
  }
}

void GeneratorsToVTU(vtkUnstructuredGrid *vtu,
                     vtkPoints *inputPoints,
                     const std::vector<Generator2> &generators) {
  const auto cd = vtu->GetCellData();

  int n_triangles = 0;
  for(auto const &g : generators)
    n_triangles += g.first.size();

  // cell data arrays
  vtkNew<vtkIntArray> triangleId{};
  triangleId->SetName("TriangleIdentifier");
  triangleId->SetNumberOfTuples(n_triangles);
  cd->AddArray(triangleId);

  vtkNew<vtkIntArray> classId{};
  classId->SetName("ClassIdentifier");
  classId->SetNumberOfTuples(n_triangles);
  cd->AddArray(classId);

  vtkNew<vtkDoubleArray> classBirth{};
  classBirth->SetName("ClassBirth");
  classBirth->SetNumberOfTuples(n_triangles);
  cd->AddArray(classBirth);

  vtkNew<vtkDoubleArray> classDeath{};
  classDeath->SetName("ClassDeath");
  classDeath->SetNumberOfTuples(n_triangles);
  cd->AddArray(classDeath);

  vtkNew<vtkDoubleArray> classPersistence{};
  classPersistence->SetName("ClassPersistence");
  classPersistence->SetNumberOfTuples(n_triangles);
  cd->AddArray(classPersistence);

  vtkNew<vtkIntArray> classDimension{};
  classDimension->SetName("ClassDimension");
  classDimension->SetNumberOfTuples(n_triangles);
  cd->AddArray(classDimension);

  // grid
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(n_triangles + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(3 * n_triangles);

  unsigned i = 0;
  for(unsigned j = 0; j < generators.size(); ++j) {
    const Generator2 &g = generators[j];
    for(auto const &t : g.first) {
      const unsigned i0 = 3 * i, i1 = 3 * i + 1, i2 = 3 * i + 2;
      triangleId->SetTuple1(i, i);
      classId->SetTuple1(i, j);
      classBirth->SetTuple1(i, g.second.first);
      classDeath->SetTuple1(i, g.second.second);
      classPersistence->SetTuple1(i, g.second.second - g.second.first);
      classDimension->SetTuple1(i, 2);

      connectivity->SetTuple1(i0, t[0]);
      connectivity->SetTuple1(i1, t[1]);
      connectivity->SetTuple1(i2, t[2]);
      offsets->SetTuple1(i, 3 * i);

      ++i;
    }
  }
  offsets->SetTuple1(n_triangles, connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  vtu->SetPoints(inputPoints);
  vtu->SetCells(VTK_TRIANGLE, cells);
}

ttkDelaunayRipsPersistenceGenerators::ttkDelaunayRipsPersistenceGenerators() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

int ttkDelaunayRipsPersistenceGenerators::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceGenerators::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkDelaunayRipsPersistenceGenerators::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer tm{};

  vtkInformation *info = inputVector[0]->GetInformationObject(0);
  vtkDataObject *input = info->Get(vtkDataObject::DATA_OBJECT());
  vtkUnstructuredGrid *outputPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);
  vtkUnstructuredGrid *outputGenerators1
    = vtkUnstructuredGrid::GetData(outputVector, 1);
  vtkUnstructuredGrid *outputGenerators2
    = vtkUnstructuredGrid::GetData(outputVector, 2);

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
  std::vector<Generator1> generators1;
  std::vector<Generator2> generators2;

  if(this->execute(points, diagram, generators1, generators2) != 0)
    return 0;

  DiagramToVTU(outputPersistenceDiagram, diagram, inf);
  const vtkNew<vtkPoints> vtkPoints{};
  MakeVtkPoints(vtkPoints, points);
  GeneratorsToVTU(outputGenerators1, vtkPoints, generators1, true);
  GeneratorsToVTU(outputGenerators2, vtkPoints, generators2);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), getThreadNumber());

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());
  outputGenerators1->GetFieldData()->ShallowCopy(input->GetFieldData());
  outputGenerators2->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}