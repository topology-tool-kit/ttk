#include <ttkRipsPersistenceDiagram.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkTable.h>
#include <boost/functional/hash.hpp>

#include <regex>

vtkStandardNewMacro(ttkRipsPersistenceDiagram);

void DiagramToVTU(vtkUnstructuredGrid *vtu,
                  const ttk::rpd::MultidimensionalDiagram &diagram,
                  double SimplexMaximumDiameter) {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  int n_pairs = 0;
  for(auto const &diagram_d : diagram)
    n_pairs += diagram_d.size();

  if(SimplexMaximumDiameter == ttk::rpd::inf) {
    double maxFiniteValue = 0.;
    for(auto const &diag : diagram) {
      for(auto const &[b, d] : diag) {
        if(d.second < ttk::rpd::inf)
          maxFiniteValue = std::max(maxFiniteValue, d.second);
      }
    }
    SimplexMaximumDiameter = 1.5 * maxFiniteValue;
  }

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
  if(!InputIsDistanceMatrix || BackEnd == BACKEND::GEOMETRY) {
    const int numberOfPoints = input->GetNumberOfRows();
    const int dimension = ScalarFields.size();

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
    "Homology maximum dimension: " + std::to_string(HomologyMaximumDimension),
    0.0, tm.getElapsedTime(), 1);
  this->printMsg(
    "Simplex maximum diameter: " + std::to_string(SimplexMaximumDiameter), 0.0,
    tm.getElapsedTime(), 1);
  if(BackEnd == BACKEND::RIPSER)
    this->printMsg("Backend: Ripser", 0.0, tm.getElapsedTime(), 1);
  else if(BackEnd == BACKEND::GEOMETRY)
    this->printMsg("Backend: Geometric", 0.0, tm.getElapsedTime(), 1);

  ttk::rpd::MultidimensionalDiagram diagram(0);

  if(this->execute(points, diagram) != 0)
    return 0;

  DiagramToVTU(
    outputPersistenceDiagram, diagram,
    (BackEnd == BACKEND::GEOMETRY) ? ttk::rpd::inf : SimplexMaximumDiameter);

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);

  // shallow copy input Field Data
  outputPersistenceDiagram->GetFieldData()->ShallowCopy(input->GetFieldData());

  // return success
  return 1;
}
