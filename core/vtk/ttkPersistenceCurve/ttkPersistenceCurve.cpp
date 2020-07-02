#include <ttkPersistenceCurve.h>

using namespace std;
using namespace ttk;

using namespace ftm;

vtkStandardNewMacro(ttkPersistenceCurve)

  ttkPersistenceCurve::ttkPersistenceCurve() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(4);
}

ttkPersistenceCurve::~ttkPersistenceCurve() {
}

vtkTable *ttkPersistenceCurve::GetOutput() {
  return this->GetOutput(0);
}

vtkTable *ttkPersistenceCurve::GetOutput(int port) {
  return vtkTable::SafeDownCast(this->GetOutputDataObject(port));
}

int ttkPersistenceCurve::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceCurve::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  switch(port) {
    case 0:
    case 1:
    case 2:
    case 3:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
      return 1;
      break;
  }
  return 0;
}

template <typename VTK_TT, typename TTK_TT>
int ttkPersistenceCurve::dispatch(vtkTable *outputJTPersistenceCurve,
                                  vtkTable *outputMSCPersistenceCurve,
                                  vtkTable *outputSTPersistenceCurve,
                                  vtkTable *outputCTPersistenceCurve,
                                  const VTK_TT *inputScalars,
                                  int inputOffsetsDataType,
                                  const void *inputOffsets,
                                  const TTK_TT *triangulation) {

  int ret = 0;
  std::vector<std::pair<VTK_TT, SimplexId>> JTPlot{};
  std::vector<std::pair<VTK_TT, SimplexId>> STPlot{};
  std::vector<std::pair<VTK_TT, SimplexId>> MSCPlot{};
  std::vector<std::pair<VTK_TT, SimplexId>> CTPlot{};

  if(inputOffsetsDataType == VTK_INT) {
    ret = this->execute<VTK_TT, int, TTK_TT>(JTPlot, STPlot, MSCPlot, CTPlot,
                                             inputScalars, (int *)inputOffsets,
                                             triangulation);
  }
  if(inputOffsetsDataType == VTK_ID_TYPE) {
    ret = this->execute<VTK_TT, vtkIdType, TTK_TT>(
      JTPlot, STPlot, MSCPlot, CTPlot, inputScalars, (vtkIdType *)inputOffsets,
      triangulation);
  }

  ret = getPersistenceCurve<vtkDoubleArray, VTK_TT>(
    outputJTPersistenceCurve, TreeType::Join, JTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("build of join tree persistence curve has failed.");
    return -1;
  }
#endif

  ret = getMSCPersistenceCurve<vtkDoubleArray, VTK_TT>(
    outputMSCPersistenceCurve, MSCPlot);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("Build of saddle-saddle persistence curve has failed.");
    return -1;
  }
#endif

  ret = getPersistenceCurve<vtkDoubleArray, VTK_TT>(
    outputSTPersistenceCurve, TreeType::Split, STPlot);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("Build of split tree persistence curve has failed.");
    return -1;
  }
#endif

  ret = getPersistenceCurve<vtkDoubleArray, VTK_TT>(
    outputCTPersistenceCurve, TreeType::Contour, CTPlot);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("Build of contour tree persistence curve has failed.");
    return -1;
  }
#endif

  return ret;
}

int ttkPersistenceCurve::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  Timer timer;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkTable *outputJTPersistenceCurve = vtkTable::GetData(outputVector, 0);
  vtkTable *outputMSCPersistenceCurve = vtkTable::GetData(outputVector, 1);
  vtkTable *outputSTPersistenceCurve = vtkTable::GetData(outputVector, 2);
  vtkTable *outputCTPersistenceCurve = vtkTable::GetData(outputVector, 3);

  vtkPointData *pointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkPersistenceCurve] Error : input has no point data." << endl;
    return -1;
  }
#endif

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif
  preconditionTriangulation(triangulation);
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

  int status = 0;
  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(), 
    (status = this->dispatch<VTK_TT, TTK_TT>(
       outputJTPersistenceCurve, outputMSCPersistenceCurve,
       outputSTPersistenceCurve, outputCTPersistenceCurve,
       (VTK_TT *)ttkUtils::GetVoidPointer(inputScalars),
       offsetField->GetDataType(), ttkUtils::GetVoidPointer(offsetField),
       (TTK_TT *)(triangulation->getData()))))
#ifndef TTK_ENABLE_KAMIKAZE
    // something wrong in baseCode
    if(status) {
    std::stringstream msg;
    msg << "PersistenceCurve::execute() error code : " << status;
    this->printErr(msg.str());
    return 0;
  }
#endif

  printMsg("Completed", 1, timer.getElapsedTime(), threadNumber_);
  printMsg(ttk::debug::Separator::L1);

  return 1;
}
