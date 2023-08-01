#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include <ttkPersistentGenerators.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistentGenerators);

ttkPersistentGenerators::ttkPersistentGenerators() {
  this->setDebugMsgPrefix("PersistentGenerators");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistentGenerators::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistentGenerators::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename triangulationType>
int ttkPersistentGenerators::dispatch(vtkPolyData *output,
                                      vtkDataArray *const inputScalarsArray,
                                      const SimplexId *const inputOrder,
                                      const triangulationType &triangulation) {

  this->buildGradient(ttkUtils::GetVoidPointer(inputScalarsArray),
                      inputScalarsArray->GetMTime(), inputOrder, triangulation);
  std::vector<GeneratorType> cycles{};
  std::vector<std::vector<SimplexId>> connComps{};
  this->computePersistentGenerators(
    cycles, connComps, inputOrder, triangulation);

  const SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  std::vector<SimplexId> isVisited(numberOfVertices, -1);
  std::vector<SimplexId> visitedIds{};

  vtkNew<vtkPoints> points{};
  vtkNew<vtkCellArray> cells{};

  vtkNew<ttkSimplexIdTypeArray> cycleId{};
  cycleId->SetName("CycleId");
  vtkNew<ttkSimplexIdTypeArray> edgeId{};
  edgeId->SetName("EdgeId");
  vtkNew<ttkSimplexIdTypeArray> birthId{};
  birthId->SetName("BirthId");
  vtkNew<ttkSimplexIdTypeArray> deathId{};
  deathId->SetName("DeathId");
  vtkNew<ttkSimplexIdTypeArray> ccId{};
  ccId->SetName("ComponentId");
  vtkNew<vtkUnsignedCharArray> iFin{};
  iFin->SetName(ttk::PersistenceIsFinite);
  vtkSmartPointer<vtkDataArray> const pers{inputScalarsArray->NewInstance()};
  pers->SetName(ttk::PersistenceName);
  vtkNew<vtkSignedCharArray> mask{};
  mask->SetName(ttk::MaskScalarFieldName);
  vtkNew<ttkSimplexIdTypeArray> vertsId{};
  vertsId->SetName(ttk::VertexScalarFieldName);
  vtkSmartPointer<vtkDataArray> const sf{inputScalarsArray->NewInstance()};
  sf->SetName(inputScalarsArray->GetName());
  vtkNew<vtkSignedCharArray> nbOnBoundary{};
  nbOnBoundary->SetNumberOfComponents(1);
  nbOnBoundary->SetName(ttk::MorseSmaleCriticalPointsOnBoundaryName);

  const auto addVertex = [&](const SimplexId v) {
    if(v == -1) {
      return vtkIdType(-1);
    }
    std::array<float, 3> p{};
    triangulation.getVertexPoint(v, p[0], p[1], p[2]);
    return points->InsertNextPoint(p.data());
  };

  const auto addEdge = [&](const SimplexId e) {
    std::array<vtkIdType, 2> pts{};
    for(int i = 0; i < 2; ++i) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, i, v);
      if(isVisited[v] == -1) {
        pts[i] = addVertex(v);
        isVisited[v] = pts[i];
        visitedIds.emplace_back(v);
      } else {
        pts[i] = isVisited[v];
      }
    }
    cells->InsertNextCell(2, pts.data());
  };

  for(size_t i = 0; i < cycles.size(); ++i) {
    if(cycles[i].boundary.empty()) {
      continue;
    }

    const auto &cycle{cycles[i]};
    const auto cbirth{cycle.boundary[0]};
    const auto cdeath{cycle.critTriangleId};
    const auto cpers{inputScalarsArray->GetTuple1(cycle.critVertsIds[0])
                     - inputScalarsArray->GetTuple1(cycle.critVertsIds[1])};

    for(const auto e : cycle.boundary) {
      addEdge(e);
      cycleId->InsertNextTuple1(i);
      edgeId->InsertNextTuple1(e);
      birthId->InsertNextTuple1(cbirth);
      deathId->InsertNextTuple1(cdeath);
      pers->InsertNextTuple1(cpers);
      iFin->InsertNextTuple1(cdeath != -1);
      nbOnBoundary->InsertNextTuple1(
        triangulation.isEdgeOnBoundary(cbirth)
        + (cdeath != -1 && triangulation.getDimensionality() == 3
             ? triangulation.isTriangleOnBoundary(cdeath)
             : 0));
    }
    for(const auto cc : connComps[i]) {
      ccId->InsertNextTuple1(cc);
    }

    // copy input scalar field
    for(const auto v : visitedIds) {
      vertsId->InsertNextTuple1(v);
      sf->InsertNextTuple1(inputScalarsArray->GetTuple1(v));
      mask->InsertNextTuple1(1);
    }

    // fill mask (0 on critical edge greater vertex, 1 on other vertices)
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(cbirth, 0, v0);
    triangulation.getEdgeVertex(cbirth, 1, v1);
    mask->SetTuple1(isVisited[v0], 0);
    mask->SetTuple1(isVisited[v1], 0);

    // ensure cycles don't share points
    for(const auto v : visitedIds) {
      isVisited[v] = -1;
    }
    visitedIds.clear();
  }

  output->SetPoints(points);
  output->GetPointData()->AddArray(sf);
  output->GetPointData()->AddArray(mask);
  output->GetPointData()->AddArray(vertsId);
  output->SetLines(cells);
  output->GetCellData()->AddArray(cycleId);
  output->GetCellData()->AddArray(edgeId);
  output->GetCellData()->AddArray(birthId);
  output->GetCellData()->AddArray(deathId);
  output->GetCellData()->AddArray(ccId);
  output->GetCellData()->AddArray(iFin);
  output->GetCellData()->AddArray(pers);
  output->GetCellData()->AddArray(nbOnBoundary);

  return 1;
}

int ttkPersistentGenerators::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto *input = vtkDataSet::GetData(inputVector[0]);
  auto *output = vtkPolyData::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Wrong triangulation");
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
  if(inputScalars == nullptr) {
    this->printErr("Wrong input scalars");
    return 0;
  }

  vtkDataArray *offsetField
    = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);
  if(offsetField == nullptr) {
    this->printErr("Wrong input offsets");
    return 0;
  }

  ttkTemplateMacro(
    triangulation->getType(),
    this->dispatch(output, inputScalars,
                   ttkUtils::GetPointer<SimplexId>(offsetField),
                   *static_cast<TTK_TT *>(triangulation->getData())));

  return 1;
}
