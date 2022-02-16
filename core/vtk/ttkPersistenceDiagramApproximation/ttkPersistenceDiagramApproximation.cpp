#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkPersistenceDiagramApproximation.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagramApproximation);

ttkPersistenceDiagramApproximation::ttkPersistenceDiagramApproximation() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkPersistenceDiagramApproximation::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramApproximation::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 1) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  } else if(port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkPersistenceDiagramApproximation::setPersistenceDiagram(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  const std::vector<ttk::PersistencePair> &diagram,
  vtkDataArray *inputScalarsArray,
  const scalarType *const inputScalars,
  const triangulationType *triangulation) const {

  if(diagram.empty()) {
    return 1;
  }

  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  // point data arrays

  vtkNew<ttkSimplexIdTypeArray> vertexIdentifierScalars{};
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);
  vertexIdentifierScalars->SetNumberOfTuples(2 * diagram.size());

  vtkNew<vtkIntArray> nodeTypeScalars{};
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName(ttk::PersistenceCriticalTypeName);
  nodeTypeScalars->SetNumberOfTuples(2 * diagram.size());

  vtkNew<vtkFloatArray> coordsScalars{};
  vtkSmartPointer<vtkDataArray> birthScalars{inputScalarsArray->NewInstance()};
  vtkSmartPointer<vtkDataArray> deathScalars{inputScalarsArray->NewInstance()};

  if(this->ShowInsideDomain) {
    birthScalars->SetNumberOfComponents(1);
    birthScalars->SetName(ttk::PersistenceBirthName);
    birthScalars->SetNumberOfTuples(2 * diagram.size());

    deathScalars->SetNumberOfComponents(1);
    deathScalars->SetName(ttk::PersistenceDeathName);
    deathScalars->SetNumberOfTuples(2 * diagram.size());
  } else {
    coordsScalars->SetNumberOfComponents(3);
    coordsScalars->SetName(ttk::PersistenceCoordinatesName);
    coordsScalars->SetNumberOfTuples(2 * diagram.size());
  }

  // cell data arrays

  vtkNew<ttkSimplexIdTypeArray> pairIdentifierScalars{};
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName(ttk::PersistencePairIdentifierName);
  pairIdentifierScalars->SetNumberOfTuples(diagram.size());

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName(ttk::PersistenceName);
  persistenceScalars->SetNumberOfTuples(diagram.size());

  vtkNew<vtkIntArray> extremumIndexScalars{};
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName(ttk::PersistencePairTypeName);
  extremumIndexScalars->SetNumberOfTuples(diagram.size());

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * diagram.size());
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(diagram.size() + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * diagram.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto a = diagram[i].birth;
    const auto b = diagram[i].death;
    const auto ta = diagram[i].birthType;
    const auto tb = diagram[i].deathType;
    // inputScalarsArray->GetTuple is not thread safe...
    const auto sa = inputScalars[a];
    const auto sb = inputScalars[b];

    if(this->ShowInsideDomain) {
      std::array<float, 3> coords{};
      triangulation->getVertexPoint(a, coords[0], coords[1], coords[2]);
      points->SetPoint(2 * i, coords[0], coords[1], coords[2]);
      triangulation->getVertexPoint(b, coords[0], coords[1], coords[2]);
      points->SetPoint(2 * i + 1, coords[0], coords[1], coords[2]);
    } else {
      points->SetPoint(2 * i, sa, sa, 0);
      points->SetPoint(2 * i + 1, sa, sb, 0);
    }
    connectivity->SetTuple1(2 * i, 2 * i);
    connectivity->SetTuple1(2 * i + 1, 2 * i + 1);
    offsets->SetTuple1(i, 2 * i);

    // point data
    vertexIdentifierScalars->SetTuple1(2 * i, a);
    vertexIdentifierScalars->SetTuple1(2 * i + 1, b);
    nodeTypeScalars->SetTuple1(2 * i, static_cast<ttk::SimplexId>(ta));
    nodeTypeScalars->SetTuple1(2 * i + 1, static_cast<ttk::SimplexId>(tb));

    if(this->ShowInsideDomain) {
      birthScalars->SetTuple1(2 * i, sa);
      birthScalars->SetTuple1(2 * i + 1, sa);
      deathScalars->SetTuple1(2 * i, sa);
      deathScalars->SetTuple1(2 * i + 1, sb);
    } else {
      std::array<float, 3> coords{};
      triangulation->getVertexPoint(a, coords[0], coords[1], coords[2]);
      coordsScalars->SetTuple3(2 * i, coords[0], coords[1], coords[2]);
      triangulation->getVertexPoint(b, coords[0], coords[1], coords[2]);
      coordsScalars->SetTuple3(2 * i + 1, coords[0], coords[1], coords[2]);
    }

    // cell data
    pairIdentifierScalars->SetTuple1(i, i);
    persistenceScalars->SetTuple1(i, diagram[i].persistence);
    if(i == 0) {
      extremumIndexScalars->SetTuple1(i, -1);
    } else {
      const auto type = diagram[i].pairType;
      if(type == 0) {
        extremumIndexScalars->SetTuple1(i, minIndex);
      } else if(type == 1) {
        extremumIndexScalars->SetTuple1(i, saddleSaddleIndex);
      } else if(type == 2) {
        extremumIndexScalars->SetTuple1(i, maxIndex);
      }
    }
  }
  offsets->SetTuple1(diagram.size(), connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  persistenceDiagram->SetPoints(points);
  persistenceDiagram->SetCells(VTK_LINE, cells);

  if(!this->ShowInsideDomain) {
    // add diagonal (first point -> last birth/penultimate point)
    std::array<vtkIdType, 2> diag{0, 2 * (cells->GetNumberOfCells() - 1)};
    persistenceDiagram->InsertNextCell(VTK_LINE, 2, diag.data());
    pairIdentifierScalars->InsertTuple1(diagram.size(), -1);
    extremumIndexScalars->InsertTuple1(diagram.size(), -1);
    // persistence of min-max pair
    const auto maxPersistence = diagram[0].persistence;
    persistenceScalars->InsertTuple1(diagram.size(), 2 * maxPersistence);
  }

  // add data arrays
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  if(this->ShowInsideDomain) {
    persistenceDiagram->GetPointData()->AddArray(birthScalars);
    persistenceDiagram->GetPointData()->AddArray(deathScalars);
  } else {
    persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  }
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);

  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkPersistenceDiagramApproximation::dispatch(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  vtkUnstructuredGrid *outputBounds,
  vtkDataArray *const inputScalarsArray,
  const scalarType *const inputScalars,
  scalarType *outputScalars,
  SimplexId *outputOffsets,
  int *outputMonotonyOffsets,
  const SimplexId *const inputOrder,
  const triangulationType *triangulation) {

  int status{};
  std::vector<ttk::PersistencePair> CTDiagram{};

  BackEnd = BACKEND::APPROXIMATE_TOPOLOGY;

  double *range = inputScalarsArray->GetRange(0);
  this->setDeltaApproximate(range[1] - range[0]);
  this->setOutputScalars(outputScalars);
  this->setOutputOffsets(outputOffsets);
  this->setOutputMonotonyOffsets(outputMonotonyOffsets);

  status = this->execute(CTDiagram, inputScalars, inputOrder, triangulation);

  // something wrong in baseCode
  if(status != 0) {
    this->printErr("PersistenceDiagram::execute() error code : "
                   + std::to_string(status));
    return 0;
  }

  setPersistenceDiagram(outputCTPersistenceDiagram, CTDiagram,
                        inputScalarsArray, outputScalars, triangulation);

  if(BackEnd == BACKEND::APPROXIMATE_TOPOLOGY) {
    drawBottleneckBounds(outputBounds, CTDiagram, inputScalarsArray,
                         outputScalars, inputScalars, triangulation);
  }

  return 1;
}

int ttkPersistenceDiagramApproximation::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputCTPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);
  vtkDataSet *outputField = vtkDataSet::GetData(outputVector, 1);
  vtkUnstructuredGrid *outputBounds
    = vtkUnstructuredGrid::GetData(outputVector, 2);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars) {
    this->printErr("Wrong input scalars");
    return 0;
  }
#endif

  vtkDataArray *offsetField
    = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!offsetField) {
    this->printErr("Wrong input offsets");
    return 0;
  }
  if(offsetField->GetDataType() != VTK_INT
     and offsetField->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported");
    return 0;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  outputOffsets->SetNumberOfComponents(1);
  outputOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputOffsets->SetName("outputOffsets");

  vtkNew<vtkIntArray> outputMonotonyOffsets{};
  outputMonotonyOffsets->SetNumberOfComponents(1);
  outputMonotonyOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputMonotonyOffsets->SetName("outputMonotonyffsets");
  outputMonotonyOffsets->FillComponent(0, 0);

  vtkSmartPointer<vtkDataArray> outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->SetNumberOfComponents(1);
  outputScalars->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputScalars->DeepCopy(inputScalars);

  std::stringstream ss;
  ss << inputScalars->GetName() << "_Approximated";
  outputScalars->SetName(ss.str().c_str());

  int status{};
  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    status = this->dispatch(
      outputCTPersistenceDiagram, outputBounds, inputScalars,
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
      static_cast<int *>(ttkUtils::GetVoidPointer(outputMonotonyOffsets)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsetField)),
      static_cast<TTK_TT *>(triangulation->getData())));

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  if(BackEnd == BACKEND::APPROXIMATE_TOPOLOGY) {
    outputField->ShallowCopy(input);
    outputField->GetPointData()->AddArray(outputScalars);
    outputField->GetPointData()->AddArray(outputOffsets);
    outputField->GetPointData()->AddArray(outputMonotonyOffsets);
  } else {
    printWrn("The exact Persistence Diagram was computed");
    printWrn("Other outputs are empty");
  }

  return status;
}

template <typename scalarType, typename triangulationType>
int ttkPersistenceDiagramApproximation::drawBottleneckBounds(
  vtkUnstructuredGrid *outputBounds,
  const std::vector<ttk::PersistencePair> &diagram,
  vtkDataArray *inputScalarsArray,
  const scalarType *const outputScalars,
  const scalarType *const inputScalars,
  const triangulationType *triangulation) const {

  if(diagram.empty()) {
    return 1;
  }

  vtkNew<vtkUnstructuredGrid> bounds{};
  double *range = inputScalarsArray->GetRange(0);
  double delta = (range[1] - range[0]) * Epsilon;
  // std::cout << "DELTA for BOUNDS " << delta << " " << getEpsilon() <<
  // std::endl;

  vtkNew<ttkSimplexIdTypeArray> BirthId{};
  BirthId->SetNumberOfComponents(1);
  BirthId->SetName("BirthId");
  vtkNew<ttkSimplexIdTypeArray> DeathId{};
  DeathId->SetNumberOfComponents(1);
  DeathId->SetName("DeathId");

  vtkNew<ttkSimplexIdTypeArray> pairIdentifierScalars{};
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> unfolded{};
  unfolded->SetNumberOfComponents(1);
  unfolded->SetName("TrueValues");

  vtkNew<vtkIntArray> extremumIndexScalars{};
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  vtkNew<vtkPoints> points{};

  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto a = diagram[i].birth;
    const auto b = diagram[i].death;

    const auto sa = outputScalars[a];
    const auto sb = outputScalars[b];
    if(sb - sa >= 2 * delta) {
      vtkIdType p0 = points->InsertNextPoint(sa - delta, sb - delta, 0);
      vtkIdType p1 = points->InsertNextPoint(sa + delta, sb - delta, 0);
      vtkIdType p2 = points->InsertNextPoint(sa - delta, sb + delta, 0);
      vtkIdType p3 = points->InsertNextPoint(sa + delta, sb + delta, 0);

      vtkNew<vtkIdList> quad{};
      quad->InsertNextId(p0);
      quad->InsertNextId(p1);
      quad->InsertNextId(p3);
      quad->InsertNextId(p2);
      bounds->InsertNextCell(VTK_QUAD, quad);

      pairIdentifierScalars->InsertNextTuple1(i);
      BirthId->InsertNextTuple1(a);
      DeathId->InsertNextTuple1(b);
      persistenceScalars->InsertNextTuple1(sb - sa);
      if(i == 0) {
        extremumIndexScalars->InsertNextTuple1(-1);
      } else {
        const auto type = diagram[i].pairType;
        if(type == 0) {
          extremumIndexScalars->InsertNextTuple1(minIndex);
        } else if(type == 1) {
          extremumIndexScalars->InsertNextTuple1(saddleSaddleIndex);
        } else if(type == 2) {
          extremumIndexScalars->InsertNextTuple1(maxIndex);
        }
      }
      if((abs((double)(outputScalars[a] - inputScalars[a])) > 1e-6)
         or (abs((double)(outputScalars[b] - inputScalars[b])) > 1e-6)) {
        unfolded->InsertNextTuple1(0);
      } else {
        unfolded->InsertNextTuple1(1);
      }
    }
  }

  // ____________________________________________
  // Case of the uncertain zone, near the diagonal
  // ____________________________________________
  const auto s0 = inputScalars[diagram[0].birth];
  const auto s2 = inputScalars[diagram[0].death];
  auto s1 = inputScalars[diagram.back().birth];
  s1 = s1 > s2 / 2 ? s1 : s2 / 2;
  vtkIdType p0 = points->InsertNextPoint(s0, s0, 0);
  vtkIdType p1 = points->InsertNextPoint(s0, s0 + 2 * delta, 0);
  vtkIdType p2 = points->InsertNextPoint(s1, s1, 0);
  vtkIdType p3 = points->InsertNextPoint(s1, s1 + 2 * delta, 0);

  vtkNew<vtkIdList> quad{};
  quad->InsertNextId(p0);
  quad->InsertNextId(p1);
  quad->InsertNextId(p3);
  quad->InsertNextId(p2);
  bounds->InsertNextCell(VTK_QUAD, quad);
  pairIdentifierScalars->InsertNextTuple1(-1);
  BirthId->InsertNextTuple1(-1);
  DeathId->InsertNextTuple1(-1);
  persistenceScalars->InsertNextTuple1(Epsilon);
  extremumIndexScalars->InsertNextTuple1(-1);
  unfolded->InsertNextTuple1(-1);
  // ____________________________________________

  bounds->SetPoints(points);

  bounds->GetCellData()->AddArray(BirthId);
  bounds->GetCellData()->AddArray(DeathId);
  bounds->GetCellData()->AddArray(persistenceScalars);
  bounds->GetCellData()->AddArray(pairIdentifierScalars);
  bounds->GetCellData()->AddArray(extremumIndexScalars);
  bounds->GetCellData()->AddArray(unfolded);

  outputBounds->ShallowCopy(bounds);

  return 0;
}
