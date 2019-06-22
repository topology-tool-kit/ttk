#include <ttkTrackingFromPersistenceDiagrams.h>

vtkStandardNewMacro(ttkTrackingFromPersistenceDiagrams)

  using dataType = double;

ttkTrackingFromPersistenceDiagrams::ttkTrackingFromPersistenceDiagrams() {
  outputMesh_ = nullptr;
  UseAllCores = false;

  DistanceAlgorithm = "ttk";
  PVAlgorithm = -1;
  Alpha = 1.0;
  Tolerance = 1.0;
  DoPostProc = false;
  PostProcThresh = 0.0;
  PX = 1;
  PY = 1;
  PZ = 1;
  PE = 1;
  PS = 1;

  WassersteinMetric = "1";
  UseGeometricSpacing = false;
  Is3D = true;
  Spacing = 1.0;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkTrackingFromPersistenceDiagrams::~ttkTrackingFromPersistenceDiagrams() {
  if(outputMesh_)
    outputMesh_->Delete();
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkTrackingFromPersistenceDiagrams::needsToAbort() {
  return GetAbortExecute() > 0;
}

// transmit progress status -- to copy paste in other wrappers
int ttkTrackingFromPersistenceDiagrams::updateProgress(const float &progress) {
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] " << progress * 100
        << "% processed...." << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkTrackingFromPersistenceDiagrams::FillInputPortInformation(
  int port, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  // Output pointers informations
  vtkInformation *outInfo;

  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkUnstructuredGrid::DATA_OBJECT()));

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] Number of inputs: "
        << numInputs << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  // Get input data
  std::vector<vtkDataSet *> input(numInputs);
  for(int i = 0; i < numInputs; i++) {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  doIt<double>(input, mesh, numInputs);

  {
    std::stringstream msg;
    msg << "[ttkTrackingFromPersistenceDiagrams] Memory usage: "
        << m.getElapsedUsage() << " MB." << std::endl;
    dMsg(std::cout, msg.str(), memoryMsg);
  }

  return 1;
}
