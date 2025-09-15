#include <ttkRipsPersistenceDiagram.h>
#include <ttkRipsPersistenceGenerators.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkRipsPersistenceGenerators);

using EdgeParametrization
  = std::unordered_map<ttk::rpd::Edge, double, boost::hash<ttk::rpd::Edge>>;
static void ParametrizeGenerator(EdgeParametrization &parametrization,
                                 const ttk::rpd::Generator1 &generator) {

  const int n = generator.first.size();
  int id_a = generator.first[0].first;
  int id_b = generator.first[0].second;
  parametrization[generator.first[0]] = 0.;
  for(int i = 1; i < n; ++i) {
    for(const ttk::rpd::Edge &e : generator.first) {
      if(e.first == id_b && e.second != id_a) {
        parametrization[e] = double(i) / n;
        id_a = id_b;
        id_b = e.second;
        break;
      }
      if(e.second == id_b && e.first != id_a) {
        parametrization[e] = double(i) / n;
        id_a = id_b;
        id_b = e.first;
        break;
      }
    }
  }
}

void GeneratorsToVTU(vtkUnstructuredGrid *vtu,
                     vtkPoints *inputPoints,
                     const std::vector<ttk::rpd::Generator1> &generators1,
                     bool parametrize) {

  const auto cd = vtu->GetCellData();

  int n_edges = 0;
  for(auto const &g : generators1)
    n_edges += g.first.size();

  // cell data arrays
  vtkNew<vtkIntArray> simplexId{};
  simplexId->SetName("simplexIdentifier");
  simplexId->SetNumberOfTuples(n_edges);
  cd->AddArray(simplexId);

  vtkNew<vtkIntArray> classId{};
  classId->SetName("ClassIdentifier");
  classId->SetNumberOfTuples(n_edges);
  cd->AddArray(classId);

  vtkNew<vtkDoubleArray> classBirth{};
  classBirth->SetName("ClassBirth");
  classBirth->SetNumberOfTuples(n_edges);
  cd->AddArray(classBirth);

  vtkNew<vtkDoubleArray> classDeath{};
  classDeath->SetName("ClassDeath");
  classDeath->SetNumberOfTuples(n_edges);
  cd->AddArray(classDeath);

  vtkNew<vtkDoubleArray> classPersistence{};
  classPersistence->SetName("ClassPersistence");
  classPersistence->SetNumberOfTuples(n_edges);
  cd->AddArray(classPersistence);

  vtkNew<vtkIntArray> classDimension{};
  classDimension->SetName("ClassDimension");
  classDimension->SetNumberOfTuples(n_edges);
  cd->AddArray(classDimension);

  vtkNew<vtkDoubleArray> generatorParametrization{};
  generatorParametrization->SetName("GeneratorParametrization");
  generatorParametrization->SetNumberOfTuples(n_edges);
  cd->AddArray(generatorParametrization);

  // grid
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(n_edges + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * n_edges);

  unsigned i = 0;
  for(unsigned j = 0; j < generators1.size(); ++j) {
    const ttk::rpd::Generator1 &g = generators1[j];
    EdgeParametrization parametrization;
    if(parametrize)
      ParametrizeGenerator(parametrization, g);
    for(auto const &e : g.first) {
      const unsigned i0 = 2 * i, i1 = 2 * i + 1;
      simplexId->SetTuple1(i, i);
      classId->SetTuple1(i, j);
      classBirth->SetTuple1(i, g.second.first);
      classDeath->SetTuple1(i, g.second.second);
      classPersistence->SetTuple1(i, g.second.second - g.second.first);
      classDimension->SetTuple1(i, 1);
      if(parametrize)
        generatorParametrization->SetTuple1(i, parametrization[e]);
      else
        generatorParametrization->SetTuple1(i, 0.);

      connectivity->SetTuple1(i0, e.first);
      connectivity->SetTuple1(i1, e.second);
      offsets->SetTuple1(i, 2 * i);

      ++i;
    }
  }
  offsets->SetTuple1(n_edges, connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  vtu->SetPoints(inputPoints);
  vtu->SetCells(VTK_LINE, cells);
}

ttkRipsPersistenceGenerators::ttkRipsPersistenceGenerators() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

int ttkRipsPersistenceGenerators::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceGenerators::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceGenerators::RequestData(
  vtkInformation *ttkNotUsed(request),
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
  } else if(input->GetNumberOfRows() != pointSet->GetNumberOfPoints()) {
    this->printErr("Input and 3D representation have different");
    this->printErr("numbers of points: resp. "
                   + std::to_string(input->GetNumberOfRows()) + " and "
                   + std::to_string(pointSet->GetNumberOfPoints()));
    return 0;
  }

  std::vector<vtkAbstractArray *> arrays;
  arrays.reserve(ScalarFields.size());
  for(const auto &s : ScalarFields)
    arrays.push_back(input->GetColumnByName(s.data()));

  std::vector<std::vector<double>> points;

  if(!InputIsDistanceMatrix) {
    const int numberOfPoints = input->GetNumberOfRows();
    const int dimension = ScalarFields.size();
    points.resize(numberOfPoints);
    for(int i = 0; i < numberOfPoints; ++i) {
      for(int j = 0; j < dimension; ++j)
        points[i].push_back(arrays[j]->GetVariantValue(i).ToDouble());
    }
    this->printMsg(
      "Computing Rips pers. generators", 0.0, tm.getElapsedTime(), 1);
    this->printMsg("#dimensions: " + std::to_string(dimension)
                     + ", #points: " + std::to_string(numberOfPoints),
                   0.0, tm.getElapsedTime(), 1);
  }

  else {
    const unsigned n = input->GetNumberOfRows();
    if(n != ScalarFields.size()) {
      this->printErr("Input distance matrix is not squared.");
      this->printErr("(rows: " + std::to_string(input->GetNumberOfRows())
                     + ", columns: " + std::to_string(ScalarFields.size())
                     + ")");
      return 0;
    }

    points = {std::vector<double>(n * (n - 1) / 2)};
    for(unsigned i = 1; i < n; ++i) {
      for(unsigned j = 0; j < i; ++j)
        points[0][i * (i - 1) / 2 + j]
          = arrays[j]->GetVariantValue(i).ToDouble();
    }
    this->printMsg(
      "Computing Rips pers. generators", 0.0, tm.getElapsedTime(), 1);
    this->printMsg(
      "(" + std::to_string(n) + "x" + std::to_string(n) + " distance matrix)",
      0.0, tm.getElapsedTime(), 1);
  }

  std::vector<ttk::rpd::Diagram> diagram(0);
  std::vector<ttk::rpd::Generator1> generators(0);
  this->execute(points, diagram, generators);

  GeneratorsToVTU(
    outputGenerators, pointSet->GetPoints(), generators, !OutputCascade);
  DiagramToVTU(outputPersistenceDiagram, diagram, SimplexMaximumDiameter);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}
