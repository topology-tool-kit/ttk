#include <ttkSurfaceQuadrangulation.h>

#define MODULE_S "[ttkSurfaceQuadrangulation] "
#define MODULE_ERROR_S MODULE_S "Error: "
#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)    \
  if(COND) {                            \
    cerr << MODULE_ERROR_S MSG << endl; \
    return RET;                         \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro(ttkSurfaceQuadrangulation);

ttkSurfaceQuadrangulation::ttkSurfaceQuadrangulation()
  : UseAllCores{true}, ThreadNumber{} {

  // critical points + 1-separatrices + segmentation
  SetNumberOfInputPorts(3);
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

int ttkSurfaceQuadrangulation::getCriticalPoints(vtkUnstructuredGrid *input) {

  // store handler to points
  auto cp = input->GetPoints();

  // store handle to points identifiers
  auto pointData = input->GetPointData();
  auto cpci = pointData->GetArray("CellId");

  TTK_ABORT_KK(cp == nullptr, "wrong Morse-Smale critical points", -1);
  TTK_ABORT_KK(cpci == nullptr, "wrong critical points cell identifiers", -2);

  surfaceQuadrangulation_.setCriticalPointsNumber(cp->GetNumberOfPoints());
  surfaceQuadrangulation_.setCriticalPointsCellIds(cpci->GetVoidPointer(0));

  return 0;
}

int ttkSurfaceQuadrangulation::getSeparatrices(vtkUnstructuredGrid *input) {

  // get separatrices points
  auto separatrices = input->GetPoints();

  auto cellData = input->GetCellData();
  auto sepSourceId = cellData->GetArray("SourceId");
  auto sepDestId = cellData->GetArray("DestinationId");

  TTK_ABORT_KK(
    separatrices == nullptr, "wrong Morse-Smale separatrices points", -1);
  TTK_ABORT_KK(sepSourceId == nullptr, "wrong separatrices source id", -2);
  TTK_ABORT_KK(sepDestId == nullptr, "wrong separatrices desination id", -3);

  surfaceQuadrangulation_.setSeparatrices(sepSourceId->GetNumberOfValues(),
                                          sepSourceId->GetVoidPointer(0),
                                          sepDestId->GetVoidPointer(0));
  return 0;
}

int ttkSurfaceQuadrangulation::getSegmentation(vtkUnstructuredGrid *input) {

  auto segmentation = input->GetPointData();
  auto segmf = segmentation->GetArray("MorseSmaleManifold");

  TTK_ABORT_KK(segmentation == nullptr, "wrong Morse-Smale segmentation", -1);
  TTK_ABORT_KK(segmf == nullptr, "wrong segmentation manifold data", -1);

  surfaceQuadrangulation_.setSegmentation(
    segmf->GetNumberOfValues(), segmf->GetVoidPointer(0));

  return 0;
}

int ttkSurfaceQuadrangulation::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  auto cp = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto spr = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto seg = vtkUnstructuredGrid::SafeDownCast(inputs[2]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getCriticalPoints(cp);

  TTK_ABORT_KK(res != 0, "wrong points", -1);

  res += getSeparatrices(spr);

  TTK_ABORT_KK(res != 0, "wrong separatrices", -1);

  surfaceQuadrangulation_.setOutputCells(&outQuadrangles_);

  res += surfaceQuadrangulation_.execute();

  TTK_ABORT_KK(
    res != 0, "SurfaceQuadrangulation.execute() error code: " << res, -9);

  // update result: get critical points from input
  output->ShallowCopy(cp);

  // vtkCellArray of quadrangle values containing outArray
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < outQuadrangles_.size() / 5; i++) {
    cells->InsertNextCell(4, &outQuadrangles_[5 * i + 1]);
  }

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
