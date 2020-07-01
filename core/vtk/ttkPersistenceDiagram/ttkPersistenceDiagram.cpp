#include <ttkPersistenceDiagram.h>

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
    switch(scalarDataType) {
      vtkTemplateMacro(deleteDiagram<VTK_TT>());
    }
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

template <typename VTK_TT, typename TTK_TT>
int ttkPersistenceDiagram::dispatch(vtkUnstructuredGrid *outputCTPersistenceDiagram,
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
      ret = this->execute<VTK_TT, int, TTK_TT>(*CTDiagram, 
        inputScalars, (int *)inputOffsets,
        triangulation);
    if(inputOffsetsDataType == VTK_ID_TYPE)
      ret = this->execute<VTK_TT, vtkIdType, TTK_TT>(
        *CTDiagram, inputScalars, (vtkIdType *)inputOffsets,
        triangulation);
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
    ret = getPersistenceDiagramInsideDomain<VTK_TT>(outputCTPersistenceDiagram,
      ftm::TreeType::Contour, *CTDiagram, inputScalarDataArray, triangulation);
  else
    ret = getPersistenceDiagram<VTK_TT>(outputCTPersistenceDiagram, ftm::TreeType::Contour, *CTDiagram, inputScalarDataArray, triangulation);
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

  Timer timer;

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
#ifndef TTK_ENABLE_KAMIKAZE
    // something wrong in baseCode
    if(status) {
    std::stringstream msg;
    msg << "PersistenceDiagram::execute() error code : " << status;
    this->printErr(msg.str());
    return 0;
  }
#endif

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  computeDiagram_ = false;

  printMsg("Completed", 1, timer.getElapsedTime(), threadNumber_);
  printMsg(ttk::debug::Separator::L1);

  return 1;
}
