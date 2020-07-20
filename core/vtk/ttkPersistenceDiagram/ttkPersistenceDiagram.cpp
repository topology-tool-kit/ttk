#include <ttkMacros.h>
#include <ttkPersistenceDiagram.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;
using namespace dcg;

vtkStandardNewMacro(ttkPersistenceDiagram)

  ttkPersistenceDiagram::ttkPersistenceDiagram() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkPersistenceDiagram::~ttkPersistenceDiagram() {
  if(CTDiagram_) {
    switch(scalarDataType) { vtkTemplateMacro(deleteDiagram<VTK_TT>()); }
  }
}

int ttkPersistenceDiagram::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
  }
  return 1;
}

template <typename VTK_TT>
int ttkPersistenceDiagram::deleteDiagram() {
  using tuple_t = tuple<SimplexId, CriticalType, SimplexId, CriticalType,
                        VTK_TT, SimplexId>;
  vector<tuple_t> *CTDiagram = (vector<tuple_t> *)CTDiagram_;
  delete CTDiagram;
  return 0;
}

template <typename scalarType,
          typename vtkSimplexArray,
          class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagramInfo(
  ttk::SimplexId id,
  vtkSmartPointer<vtkSimplexArray> vertexIdentifierScalars,
  vtkSmartPointer<vtkIntArray> nodeTypeScalars,
  vtkSmartPointer<vtkFloatArray> coordsScalars,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkSmartPointer<vtkPoints> points,
  vtkIdType ids[3],
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  double p[3] = {0, 0, 0};
  const ttk::SimplexId a = std::get<0>(diagram[id]);
  const ttk::SimplexId na
    = static_cast<ttk::SimplexId>(std::get<1>(diagram[id]));
  const ttk::SimplexId b = std::get<2>(diagram[id]);
  const ttk::SimplexId nb
    = static_cast<ttk::SimplexId>(std::get<3>(diagram[id]));

  nodeTypeScalars->InsertTuple1(2 * id, na);
  nodeTypeScalars->InsertTuple1(2 * id + 1, nb);

  vertexIdentifierScalars->InsertTuple1(2 * id, a);
  vertexIdentifierScalars->InsertTuple1(2 * id + 1, b);

  float coords[3];
  triangulation->getVertexPoint(a, coords[0], coords[1], coords[2]);
  coordsScalars->InsertTuple3(2 * id, coords[0], coords[1], coords[2]);

  triangulation->getVertexPoint(b, coords[0], coords[1], coords[2]);
  coordsScalars->InsertTuple3(2 * id + 1, coords[0], coords[1], coords[2]);

  p[0] = inputScalars->GetTuple1(a);
  p[1] = inputScalars->GetTuple1(a);
  ids[0] = points->InsertNextPoint(p);

  p[0] = inputScalars->GetTuple1(a);
  p[1] = inputScalars->GetTuple1(b);
  ids[1] = points->InsertNextPoint(p);

  return 0;
}

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::getPersistenceDiagram(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkSmartPointer<vtkIntArray> nodeTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkSmartPointer<ttkSimplexIdTypeArray> pairIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> extremumIndexScalars
    = vtkSmartPointer<vtkIntArray>::New();
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkSmartPointer<vtkFloatArray> coordsScalars
    = vtkSmartPointer<vtkFloatArray>::New();
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];
    vtkIdType oldIds[2];

    scalarType maxPersistenceValue = std::numeric_limits<scalarType>::min();
    oldIds[0] = 0;
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const scalarType persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      setPersistenceDiagramInfo(i, vertexIdentifierScalars, nodeTypeScalars,
                                coordsScalars, diagram, points, ids,
                                inputScalars, triangulation);

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
    oldIds[1] = ids[0];

    // add diag
    persistenceDiagram->InsertNextCell(VTK_LINE, 2, oldIds);
    pairIdentifierScalars->InsertTuple1(diagramSize, -1);
    extremumIndexScalars->InsertTuple1(diagramSize, -1);
    persistenceScalars->InsertTuple1(diagramSize, 2 * maxPersistenceValue);
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);

  return 0;
}

template <typename scalarType,
          typename vtkSimplexArray,
          class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagramInfoInsideDomain(
  ttk::SimplexId id,
  vtkSmartPointer<vtkSimplexArray> vertexIdentifierScalars,
  vtkSmartPointer<vtkIntArray> nodeTypeScalars,
  vtkDataArray *birthScalars,
  vtkDataArray *deathScalars,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkSmartPointer<vtkPoints> points,
  vtkIdType ids[3],
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  float p[3];
  const ttk::SimplexId a = std::get<0>(diagram[id]);
  const ttk::SimplexId na
    = static_cast<ttk::SimplexId>(std::get<1>(diagram[id]));
  const ttk::SimplexId b = std::get<2>(diagram[id]);
  const ttk::SimplexId nb
    = static_cast<ttk::SimplexId>(std::get<3>(diagram[id]));
  const double sa = inputScalars->GetTuple1(a);
  const double sb = inputScalars->GetTuple1(b);

  nodeTypeScalars->InsertTuple1(2 * id, na);
  nodeTypeScalars->InsertTuple1(2 * id + 1, nb);
  vertexIdentifierScalars->InsertTuple1(2 * id, a);
  vertexIdentifierScalars->InsertTuple1(2 * id + 1, b);
  birthScalars->InsertTuple1(2 * id, sa);
  birthScalars->InsertTuple1(2 * id + 1, sa);
  deathScalars->InsertTuple1(2 * id, sa);
  deathScalars->InsertTuple1(2 * id + 1, sb);

  triangulation->getVertexPoint(a, p[0], p[1], p[2]);
  ids[0] = points->InsertNextPoint(p);

  triangulation->getVertexPoint(b, p[0], p[1], p[2]);
  ids[1] = points->InsertNextPoint(p);

  return 0;
}

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::getPersistenceDiagramInsideDomain(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkSmartPointer<vtkIntArray> nodeTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkSmartPointer<ttkSimplexIdTypeArray> pairIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> extremumIndexScalars
    = vtkSmartPointer<vtkIntArray>::New();
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkDataArray *birthScalars = inputScalars->NewInstance();
  birthScalars->SetNumberOfComponents(1);
  birthScalars->SetName("Birth");

  vtkDataArray *deathScalars = inputScalars->NewInstance();
  deathScalars->SetNumberOfComponents(1);
  deathScalars->SetName("Death");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];

    scalarType maxPersistenceValue = std::numeric_limits<scalarType>::min();
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const scalarType persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      setPersistenceDiagramInfoInsideDomain(
        i, vertexIdentifierScalars, nodeTypeScalars, birthScalars, deathScalars,
        diagram, points, ids, inputScalars, triangulation);

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(birthScalars);
  persistenceDiagram->GetPointData()->AddArray(deathScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);
  return 0;
}

template <typename VTK_TT, typename TTK_TT>
int ttkPersistenceDiagram::dispatch(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  vtkDataArray *inputScalarDataArray,
  const VTK_TT *inputScalars,
  int inputOffsetsDataType,
  const void *inputOffsets,
  const TTK_TT *triangulation) {

  int ret = 0;

  using tuple_t = tuple<SimplexId, CriticalType, SimplexId, CriticalType,
                        VTK_TT, SimplexId>;

  if(CTDiagram_ && computeDiagram_) {
    vector<tuple_t> *tmpDiagram = (vector<tuple_t> *)CTDiagram_;
    delete tmpDiagram;
    CTDiagram_ = new vector<tuple_t>();
  } else if(!CTDiagram_) {
    CTDiagram_ = new vector<tuple_t>();
    computeDiagram_ = true;
  }

  vector<tuple_t> *CTDiagram = (vector<tuple_t> *)CTDiagram_;

  if(computeDiagram_) {
    if(inputOffsetsDataType == VTK_INT)
      ret = this->execute<VTK_TT, int, TTK_TT>(
        *CTDiagram, inputScalars, (int *)inputOffsets, triangulation);
    if(inputOffsetsDataType == VTK_ID_TYPE)
      ret = this->execute<VTK_TT, vtkIdType, TTK_TT>(
        *CTDiagram, inputScalars, (vtkIdType *)inputOffsets, triangulation);
#ifndef TTK_ENABLE_KAMIKAZE
    if(ret) {
      std::stringstream msg;
      msg << "PersistenceDiagram::execute() error code : " << ret;
      this->printErr(msg.str());
      return -4;
    }
#endif
  }

  if(ShowInsideDomain)
    ret = getPersistenceDiagramInsideDomain<VTK_TT>(
      outputCTPersistenceDiagram, ftm::TreeType::Contour, *CTDiagram,
      inputScalarDataArray, triangulation);
  else
    ret = getPersistenceDiagram<VTK_TT>(outputCTPersistenceDiagram,
                                        ftm::TreeType::Contour, *CTDiagram,
                                        inputScalarDataArray, triangulation);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("Build of contour tree persistence diagram has failed.");
    return -5;
  }
#endif

  return ret;
}

int ttkPersistenceDiagram::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputCTPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);

  vtkPointData *pointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    this->printErr("Input has no point data.");
    return 0;
  }
#endif

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif
  // TODO: Remove when FTM and MSC are migrated
  setupTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars) {
    this->printErr("Wrong input scalars");
    return 0;
  }
#endif

  vtkDataArray *offsetField = ttkAlgorithm::GetOptionalArray(
    ForceInputOffsetScalarField, 1, ttk::OffsetScalarFieldName, inputVector);

  if(!offsetField) {
    offsetField = pointData->GetArray(ttk::OffsetScalarFieldName);
  }

  if(!offsetField) {
    const SimplexId numberOfVertices = input->GetNumberOfPoints();

    offsetField = ttkSimplexIdTypeArray::New();
    offsetField->SetNumberOfComponents(1);
    offsetField->SetNumberOfTuples(numberOfVertices);
    offsetField->SetName(ttk::OffsetScalarFieldName);
    for(SimplexId i = 0; i < numberOfVertices; ++i)
      offsetField->SetTuple1(i, i);
  }
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

  if(this->GetMTime() < inputScalars->GetMTime())
    computeDiagram_ = true;

  vector<tuple<Cell, Cell>> dmt_pairs_temp;
  setDMTPairs(&dmt_pairs_temp);

  scalarDataType = inputScalars->GetDataType();

  int status = 0;
  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    (status = this->dispatch<VTK_TT, TTK_TT>(
       outputCTPersistenceDiagram, inputScalars,
       (VTK_TT *)ttkUtils::GetVoidPointer(inputScalars),
       offsetField->GetDataType(), ttkUtils::GetVoidPointer(offsetField),
       (TTK_TT *)(triangulation->getData()))))

    // something wrong in baseCode
    if(status) {
    std::stringstream msg;
    msg << "PersistenceDiagram::execute() error code : " << status;
    this->printErr(msg.str());
    return 0;
  }

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  computeDiagram_ = false;

  return 1;
}
