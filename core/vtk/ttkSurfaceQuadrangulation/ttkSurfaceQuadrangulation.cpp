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
  auto cpcd = pointData->GetArray("CellDimension");
  auto cpid = pointData->GetArray(ttk::VertexScalarFieldName);

  TTK_ABORT_KK(cp == nullptr, "wrong Morse-Smale critical points", -1);
  TTK_ABORT_KK(cpci == nullptr, "wrong critical points cell identifiers", -2);
  TTK_ABORT_KK(cpcd == nullptr, "wrong critical points cell dimension", -3);
  TTK_ABORT_KK(cpci == nullptr, "wrong critical points identifiers", -4);

  surfaceQuadrangulation_.setCriticalPointsNumber(cp->GetNumberOfPoints());
  surfaceQuadrangulation_.setCriticalPointsCellIds(cpci->GetVoidPointer(0));
  surfaceQuadrangulation_.setCriticalPointsType(cpcd->GetVoidPointer(0));
  surfaceQuadrangulation_.setCriticalPointsIdentifiers(cpid->GetVoidPointer(0));

  return 0;
}

int ttkSurfaceQuadrangulation::getSeparatrices(vtkUnstructuredGrid *input) {

  // check if separatrices have points
  auto separatrices = input->GetPoints();
  TTK_ABORT_KK(
    separatrices == nullptr, "no points in Morse-Smale separatrices", -1);

  // get separatrices point data
  auto pointData = input->GetPointData();
  auto id = pointData->GetArray("CellId");
  auto mask = pointData->GetArray(ttk::MaskScalarFieldName);

  TTK_ABORT_KK(id == nullptr, "wrong separatrices cell id", -2);
  TTK_ABORT_KK(mask == nullptr, "wrong separatrices mask", -3);

  surfaceQuadrangulation_.setSeparatrices(
    id->GetNumberOfValues(), id->GetVoidPointer(0), mask->GetVoidPointer(0),
    separatrices->GetVoidPointer(0));
  return 0;
}

int ttkSurfaceQuadrangulation::getSegmentation(vtkUnstructuredGrid *input) {

  auto segmentation = input->GetPointData();
  auto segmf = segmentation->GetArray("MorseSmaleManifold");

  TTK_ABORT_KK(segmentation == nullptr, "wrong Morse-Smale segmentation", -1);
  TTK_ABORT_KK(segmf == nullptr, "wrong segmentation manifold data", -2);

  surfaceQuadrangulation_.setSegmentation(
    segmf->GetNumberOfValues(), segmf->GetVoidPointer(0));

  triangulation_ = ttkTriangulation::getTriangulation(input);
  TTK_ABORT_KK(triangulation_ == nullptr, "invalid triangulation", -3);
  triangulation_->setWrapper(this);
  surfaceQuadrangulation_.setWrapper(this);
  surfaceQuadrangulation_.setupTriangulation(triangulation_);

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

  res += getSegmentation(seg);

  TTK_ABORT_KK(res != 0, "wrong segmentation", -1);

  surfaceQuadrangulation_.setOutputQuads(&outQuadrangles_, &outQuadPoints_);
  surfaceQuadrangulation_.setDualQuadrangulation(DualQuadrangulation);

  res += surfaceQuadrangulation_.execute();

  TTK_ABORT_KK(
    res != 0, "SurfaceQuadrangulation.execute() error code: " << res, -9);


  // output points: critical points + generated separatrices middles
  auto points = vtkSmartPointer<vtkPoints>::New();

  // copy critical points
  for(size_t i = 0; i < cp->GetNumberOfPoints(); i++) {
    points->InsertNextPoint(cp->GetPoints()->GetPoint(i));
  }

  // copy separatrices middles
  for(size_t i = 0; i < outQuadPoints_.size() / 3; i++) {
    points->InsertNextPoint(&outQuadPoints_[3 * i]);
  }

  output->SetPoints(points);

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
