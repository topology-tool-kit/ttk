#include <ttkRipsPersistenceDiagram.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <regex>

vtkStandardNewMacro(ttkRipsPersistenceDiagram);

void MakeVtkPoints(vtkPoints *vtkPoints,
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

using EdgeParametrization
  = std::unordered_map<ttk::rpd::Edge, double, boost::hash<ttk::rpd::Edge>>;
static void ParametrizeGenerator(EdgeParametrization &parametrization,
                                 const ttk::rpd::Generator &generator) {

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

void DiagramToVTU(vtkUnstructuredGrid *vtu,
                  const ttk::rpd::MultidimensionalDiagram &diagram,
                  double SimplexMaximumDiameter) {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  int n_pairs = 0;
  for(auto const &diagram_d : diagram)
    n_pairs += diagram_d.size();
  if(SimplexMaximumDiameter == ttk::rpd::inf)
    n_pairs--;

  // point data arrays
  vtkNew<ttkSimplexIdTypeArray> vertsId{};
  vertsId->SetName(ttk::VertexScalarFieldName);
  vertsId->SetNumberOfTuples(2 * n_pairs);
  pd->AddArray(vertsId);

  vtkNew<vtkIntArray> critType{};
  critType->SetName(ttk::PersistenceCriticalTypeName);
  critType->SetNumberOfTuples(2 * n_pairs);
  pd->AddArray(critType);

  // cell data arrays
  vtkNew<ttkSimplexIdTypeArray> pairsId{};
  pairsId->SetName(ttk::PersistencePairIdentifierName);
  pairsId->SetNumberOfTuples(n_pairs);
  cd->AddArray(pairsId);

  vtkNew<vtkIntArray> pairsDim{};
  pairsDim->SetName(ttk::PersistencePairTypeName);
  pairsDim->SetNumberOfTuples(n_pairs);
  cd->AddArray(pairsDim);

  vtkNew<vtkDoubleArray> persistence{};
  persistence->SetName(ttk::PersistenceName);
  persistence->SetNumberOfTuples(n_pairs);
  cd->AddArray(persistence);

  vtkNew<vtkDoubleArray> birthScalars{};
  birthScalars->SetName(ttk::PersistenceBirthName);
  birthScalars->SetNumberOfTuples(n_pairs);
  cd->AddArray(birthScalars);

  vtkNew<vtkUnsignedCharArray> isFinite{};
  isFinite->SetName(ttk::PersistenceIsFinite);
  isFinite->SetNumberOfTuples(n_pairs);
  cd->AddArray(isFinite);

  // grid
  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * n_pairs);
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(n_pairs + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * n_pairs);

  unsigned i = 0;
  unsigned i_max = 0;
  double birth_max = 0.;
  for(unsigned d = 0; d < diagram.size(); ++d) {
    for(auto const &pair : diagram[d]) {
      if(d == 0 && pair.second.second == ttk::rpd::inf
         && SimplexMaximumDiameter == ttk::rpd::inf)
        continue;

      const unsigned i0 = 2 * i, i1 = 2 * i + 1;
      pairsId->SetTuple1(i, i);
      pairsDim->SetTuple1(i, d);

      const double death = std::min(SimplexMaximumDiameter, pair.second.second);
      isFinite->SetTuple1(i, pair.second.second < ttk::rpd::inf);
      persistence->SetTuple1(i, death - pair.first.second);
      birthScalars->SetTuple1(i, pair.first.second);
      points->SetPoint(i0, pair.first.second, pair.first.second, 0);
      points->SetPoint(i1, pair.first.second, death, 0);

      if(pair.first.second > birth_max) {
        birth_max = pair.first.second;
        i_max = i;
      }

      connectivity->SetTuple1(i0, i0);
      connectivity->SetTuple1(i1, i1);
      offsets->SetTuple1(i, 2 * i);

      critType->SetTuple1(i0, d);
      critType->SetTuple1(i1, d + 1);

      vertsId->SetTuple1(i0, *std::max_element(pair.first.first.begin(),
                                               pair.first.first.end()));
      vertsId->SetTuple1(i1, *std::max_element(pair.second.first.begin(),
                                               pair.second.first.end()));

      ++i;
    }
  }

  offsets->SetTuple1(n_pairs, connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  vtu->SetPoints(points);
  vtu->SetCells(VTK_LINE, cells);

  // add diagonal
  std::array<vtkIdType, 2> diag{0, 2 * i_max};
  vtu->InsertNextCell(VTK_LINE, 2, diag.data());
  pairsId->InsertTuple1(n_pairs, -1);
  pairsDim->InsertTuple1(n_pairs, -1);
  isFinite->InsertTuple1(n_pairs, false);
  persistence->InsertTuple1(n_pairs, 0.);
  birthScalars->InsertTuple1(n_pairs, 0.);
}

void GeneratorsToVTU(vtkUnstructuredGrid *vtu,
                     vtkPoints *inputPoints,
                     const std::vector<ttk::rpd::Generator> &generators,
                     bool parametrize) {

  const auto cd = vtu->GetCellData();

  int n_edges = 0;
  for(auto const &g : generators)
    n_edges += g.first.size();

  // cell data arrays
  vtkNew<vtkIntArray> edgesId{};
  edgesId->SetName("EdgeIdentifier");
  edgesId->SetNumberOfTuples(n_edges);
  cd->AddArray(edgesId);

  vtkNew<vtkIntArray> polygonId{};
  polygonId->SetName("ClassIdentifier");
  polygonId->SetNumberOfTuples(n_edges);
  cd->AddArray(polygonId);

  vtkNew<vtkDoubleArray> polygonBirth{};
  polygonBirth->SetName("ClassBirth");
  polygonBirth->SetNumberOfTuples(n_edges);
  cd->AddArray(polygonBirth);

  vtkNew<vtkDoubleArray> polygonDeath{};
  polygonDeath->SetName("ClassDeath");
  polygonDeath->SetNumberOfTuples(n_edges);
  cd->AddArray(polygonDeath);

  vtkNew<vtkDoubleArray> polygonPersistence{};
  polygonPersistence->SetName("ClassPersistence");
  polygonPersistence->SetNumberOfTuples(n_edges);
  cd->AddArray(polygonPersistence);

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
  for(unsigned j = 0; j < generators.size(); ++j) {
    const ttk::rpd::Generator &g = generators[j];
    EdgeParametrization parametrization;
    if(parametrize)
      ParametrizeGenerator(parametrization, g);
    for(auto const &e : g.first) {
      const unsigned i0 = 2 * i, i1 = 2 * i + 1;
      edgesId->SetTuple1(i, i);
      polygonId->SetTuple1(i, j);
      polygonBirth->SetTuple1(i, g.second.first);
      polygonDeath->SetTuple1(i, g.second.second);
      polygonPersistence->SetTuple1(i, g.second.second - g.second.first);
      if(parametrize)
        generatorParametrization->SetTuple1(i, parametrization[e]);

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

ttkRipsPersistenceDiagram::ttkRipsPersistenceDiagram() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkRipsPersistenceDiagram::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceDiagram::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRipsPersistenceDiagram::RequestData(vtkInformation *ttkNotUsed(request),
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

  std::vector<std::vector<double>> points;
  bool doGenerators = false;
  if(!InputIsDistanceMatrix || BackEnd == BACKEND::GEOMETRY) {
    const int numberOfPoints = input->GetNumberOfRows();
    const int dimension = ScalarFields.size();
    doGenerators
      = OutputGenerators && (dimension == 2) && BackEnd == BACKEND::GEOMETRY;

    points = std::vector<std::vector<double>>(numberOfPoints);
    for(int i = 0; i < numberOfPoints; ++i) {
      for(int j = 0; j < dimension; ++j)
        points[i].push_back(arrays[j]->GetVariantValue(i).ToDouble());
    }
    this->printMsg(
      "Computing Rips persistence diagram", 1.0, tm.getElapsedTime(), 1);
    this->printMsg("#dimensions: " + std::to_string(dimension)
                     + ", #points: " + std::to_string(numberOfPoints),
                   0.0, tm.getElapsedTime(), 1);
  } else {
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
      "Computing Rips persistence diagram", 1.0, tm.getElapsedTime(), 1);
    this->printMsg(
      "(" + std::to_string(n) + "x" + std::to_string(n) + " distance matrix)",
      0.0, tm.getElapsedTime(), 1);
  }

  this->printMsg(
    "Simplex maximum dimension: " + std::to_string(SimplexMaximumDimension),
    0.0, tm.getElapsedTime(), 1);
  this->printMsg(
    "Simplex maximum diameter: " + std::to_string(SimplexMaximumDiameter), 0.0,
    tm.getElapsedTime(), 1);
  if(BackEnd == BACKEND::RIPSER)
    this->printMsg("Backend: Ripser", 0.0, tm.getElapsedTime(), 1);
  else if(BackEnd == BACKEND::GEOMETRY)
    this->printMsg("Backend: Geometric", 0.0, tm.getElapsedTime(), 1);

  ttk::rpd::MultidimensionalDiagram diagram(0);
  std::vector<ttk::rpd::Generator> generators(0);

  if(this->execute(points, diagram, generators) != 0)
    return 0;

  if(doGenerators) {
    vtkNew<vtkPoints> vtkPoints{};
    MakeVtkPoints(vtkPoints, points);
    GeneratorsToVTU(outputPersistenceDiagram, vtkPoints, generators, true);
  } else
    DiagramToVTU(
      outputPersistenceDiagram, diagram,
      (BackEnd == BACKEND::GEOMETRY) ? ttk::rpd::inf : SimplexMaximumDiameter);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}