#include <ttkMorseSmaleQuadrangulation.h>

#define MODULE_S "[ttkMorseSmaleQuadrangulation] "
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

vtkStandardNewMacro(ttkMorseSmaleQuadrangulation);

ttkMorseSmaleQuadrangulation::ttkMorseSmaleQuadrangulation()
  : UseAllCores{true}, ThreadNumber{} {

  // critical points + 1-separatrices + segmentation
  SetNumberOfInputPorts(3);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkMorseSmaleQuadrangulation::FillInputPortInformation(
  int port, vtkInformation *info) {

  if(port == 0 || port == 1) {
    // Morse-Smale Complex output
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 0;
}

int ttkMorseSmaleQuadrangulation::getCriticalPoints(
  vtkUnstructuredGrid *input) {

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

  baseWorker_.setCriticalPoints(
    cp->GetNumberOfPoints(), cp->GetVoidPointer(0), cpid->GetVoidPointer(0),
    cpci->GetVoidPointer(0), cpcd->GetVoidPointer(0));

  return 0;
}

int ttkMorseSmaleQuadrangulation::getSeparatrices(vtkUnstructuredGrid *input) {

  // check if separatrices have points
  auto separatrices = input->GetPoints();
  TTK_ABORT_KK(
    separatrices == nullptr, "no points in Morse-Smale separatrices", -1);

  // get separatrices point data
  auto pointData = input->GetPointData();
  auto id = pointData->GetArray("CellId");
  auto dim = pointData->GetArray("CellDimension");
  auto mask = pointData->GetArray(ttk::MaskScalarFieldName);

  TTK_ABORT_KK(id == nullptr, "wrong separatrices cell id", -2);
  TTK_ABORT_KK(dim == nullptr, "wrong separatrices cell dimension", -3);
  TTK_ABORT_KK(mask == nullptr, "wrong separatrices mask", -4);

  baseWorker_.setSeparatrices(id->GetNumberOfTuples(), id->GetVoidPointer(0),
                              dim->GetVoidPointer(0), mask->GetVoidPointer(0),
                              separatrices->GetVoidPointer(0));
  return 0;
}

int ttkMorseSmaleQuadrangulation::getTriangulation(vtkUnstructuredGrid *input) {
  triangulation_ = ttkTriangulation::getTriangulation(input);
  TTK_ABORT_KK(triangulation_ == nullptr, "invalid triangulation", -1);

  auto points = input->GetPoints();
  TTK_ABORT_KK(points == nullptr, "wrong points", -2);

  triangulation_->setWrapper(this);
  baseWorker_.setWrapper(this);
  baseWorker_.setupTriangulation(triangulation_);
  baseWorker_.setInputPoints(points->GetVoidPointer(0));

  return 0;
}

int ttkMorseSmaleQuadrangulation::doIt(std::vector<vtkDataSet *> &inputs,
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

  res += getTriangulation(seg);

  TTK_ABORT_KK(res != 0, "wrong segmentation", -1);

  baseWorker_.setDualQuadrangulation(DualQuadrangulation);
  baseWorker_.setShowResError(ShowResError);

  res += baseWorker_.execute();

  if(res != 0) {
    vtkWarningMacro(MODULE_ERROR_S
                    "Consider another (eigen) function, persistence threshold "
                    "or refine your input triangulation");
    if(!ShowResError) {
      return res;
    }
  }

  auto &outQuadrangles = baseWorker_.outputCells_;
  auto &outQuadPoints = baseWorker_.outputPoints_;
  auto &outPointsIds = baseWorker_.outputPointsIds_;
  auto &outPointsType = baseWorker_.outputPointsTypes_;
  auto &outPointsCells = baseWorker_.outputPointsCells_;

  // output points: critical points + generated separatrices middles
  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < outQuadPoints.size() / 3; i++) {
    points->InsertNextPoint(&outQuadPoints[3 * i]);
  }
  output->SetPoints(points);

  // quad vertices identifiers
  auto identifiers = vtkSmartPointer<vtkIntArray>::New();
  identifiers->SetName(ttk::VertexScalarFieldName);
  identifiers->SetVoidArray(outPointsIds.data(), outPointsIds.size(), 1);
  output->GetPointData()->AddArray(identifiers);

  // quad vertices type
  auto type = vtkSmartPointer<vtkIntArray>::New();
  type->SetName("QuadVertType");
  type->SetVoidArray(outPointsType.data(), outPointsType.size(), 1);
  output->GetPointData()->AddArray(type);

  // quad vertices cells
  auto cellid = vtkSmartPointer<vtkIntArray>::New();
  cellid->SetName("QuadCellId");
  cellid->SetVoidArray(outPointsCells.data(), outPointsCells.size(), 1);
  output->GetPointData()->AddArray(cellid);

  // vtkCellArray of quadrangle values containing outArray
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < outQuadrangles.size() / 5; i++) {
    cells->InsertNextCell(4, &outQuadrangles[5 * i + 1]);
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
