#include <ttkSurfaceQuadrangulation.h>

#define MODULE_S "[ttkSurfaceQuadrangulation] "
#define MODULE_ERROR_S MODULE_S "Error: "

vtkStandardNewMacro(ttkSurfaceQuadrangulation);

ttkSurfaceQuadrangulation::ttkSurfaceQuadrangulation()
  : UseAllCores{true}, ThreadNumber{}, ForceInputIdentifiersField{false},
    ForceInputOffsetIdentifiersField{false}, SubdivisionLevel{5},
    RelaxationIterations{100}, criticalPoints_{}, criticalPointsIdentifiers_{},
    separatrices_{} {

  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));

  // critical points + 1-separatrices
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkSurfaceQuadrangulation::FillInputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0 || port == 1) {
    // Morse-Smale Complex output
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 0;
}

int ttkSurfaceQuadrangulation::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    // surface
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  }
  return 0;
}

int ttkSurfaceQuadrangulation::getCriticalPoints(vtkUnstructuredGrid *input) {

  // store handler to points
  criticalPoints_ = input->GetPoints();

  // store handle to points identifiers
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);
  auto pointData = input->GetPointData();

  if(pointData->GetArray(vsfn) != nullptr) {
    criticalPointsIdentifiers_ = pointData->GetArray(vsfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(criticalPoints_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points" << endl;
    return -1;
  }
  if(criticalPointsIdentifiers_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale critical points identifiers"
         << endl;
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  return 0;
}

int ttkSurfaceQuadrangulation::getSeparatrices(vtkUnstructuredGrid *input) {

  // get separatrices points
  separatrices_ = input->GetPoints();

  auto cellData = input->GetCellData();
  separatrixId_ = cellData->GetArray("SeparatrixId");
  separatrixSourceId_ = cellData->GetArray("SourceId");
  separatrixDestinationId_ = cellData->GetArray("DestinationId");

#ifndef TTK_ENABLE_KAMIKAZE
  if(separatrices_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong Morse-Smale separatrices points." << endl;
    return -1;
  }
  if(separatrixId_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices id." << endl;
    return -1;
  }
  if(separatrixSourceId_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices source id." << endl;
    return -1;
  }
  if(separatrixDestinationId_ == nullptr) {
    cerr << MODULE_ERROR_S "wrong separatrices destination id." << endl;
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  return 0;
}

int ttkSurfaceQuadrangulation::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  auto criticalPoints = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto separatrices = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkPolyData::SafeDownCast(outputs[0]);

  int res = 0;

  res += getCriticalPoints(criticalPoints);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << MODULE_ERROR_S "wrong points." << endl;
    return -1;
  }
#endif

  res += getSeparatrices(separatrices);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << MODULE_ERROR_S "wrong separatrices." << endl;
    return -1;
  }
#endif

  surfaceQuadrangulation_.setSubdivisionLevel(SubdivisionLevel);
  surfaceQuadrangulation_.setRelaxationIterations(RelaxationIterations);

  std::vector<ttk::SimplexId> outputVertices;
  std::vector<std::vector<ttk::SimplexId>> outputEdges;

  res += surfaceQuadrangulation_.execute();

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/SurfaceQuadrangulation
  if(res != 0) {
    cerr << MODULE_S "SurfaceQuadrangulation.execute() error code: " << res
         << endl;
    return -9;
  }
#endif

  // update result
  auto outputQuadrangulation = vtkSmartPointer<vtkIntArray>::New();
  outputQuadrangulation->SetNumberOfComponents(1);
  outputQuadrangulation->SetVoidArray(
    static_cast<void *>(outputVertices.data()), outputVertices.size(), 1);

  // output->GetPointData()->AddArray(outputQuadrangulation);

  std::stringstream msg;
  msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
