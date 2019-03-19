#include <ttkQuadrangulationSubdivision.h>

#define MODULE_S "[ttkQuadrangulationSubdivision] "
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

vtkStandardNewMacro(ttkQuadrangulationSubdivision);

ttkQuadrangulationSubdivision::ttkQuadrangulationSubdivision()
  : UseAllCores{true}, ThreadNumber{}, ForceInputIdentifiersField{false},
    ForceInputOffsetIdentifiersField{false}, SubdivisionLevel{3},
    RelaxationIterations{100} {

  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));

  // MSC quadrangulation + initial 2D mesh
  SetNumberOfInputPorts(2);
  // quad mesh (containing ttkVertexIdentifiers of critical points)
  SetNumberOfOutputPorts(1);
}

int ttkQuadrangulationSubdivision::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  }
  return 1;
}

int ttkQuadrangulationSubdivision::getTriangulation(
  vtkUnstructuredGrid *const input) {

  auto trg = ttkTriangulation::getTriangulation(input);

  TTK_ABORT_KK(trg == nullptr, "input triangulation pointer is NULL.", -1);

  trg->setWrapper(this);
  baseWorker_.setWrapper(this);
  baseWorker_.setupTriangulation(trg);
  Modified();

  TTK_ABORT_KK(trg->isEmpty(), "ttkTriangulation allocation problem.", -2);

  return 0;
}

int ttkQuadrangulationSubdivision::getQuadVertices(
  vtkUnstructuredGrid *const input) {
  auto cells = input->GetCells();

  TTK_ABORT_KK(cells == nullptr, "invalid input quadrangle cells", -3);
  TTK_ABORT_KK(
    cells->GetData() == nullptr, "invalid input quadrangle cell data", -4);

  auto points = input->GetPoints();

  TTK_ABORT_KK(points == nullptr, "invalid input critical points", -5);
  TTK_ABORT_KK(
    points->GetData() == nullptr, "invalid input quadrangle cell data", -6);

  auto pointData = input->GetPointData();
  auto identifiers = pointData->GetArray(
    static_cast<const char *>(ttk::VertexScalarFieldName));

  TTK_ABORT_KK(pointData == nullptr, "invalid input quadrangle point data", -7);
  TTK_ABORT_KK(identifiers == nullptr,
               "invalid input quadrangle vertices identifiers", -8);

  auto address = static_cast<vtkIdType *>(cells->GetData()->GetVoidPointer(0));

  baseWorker_.setInputQuads(
    cells->GetData()->GetVoidPointer(0), cells->GetNumberOfCells());
  baseWorker_.setInputVertices(
    points->GetData()->GetVoidPointer(0), points->GetNumberOfPoints());
  baseWorker_.setInputVertexIdentifiers(identifiers->GetVoidPointer(0));

  return 0;
}

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  baseWorker_.setSubdivisionLevel(SubdivisionLevel);
  baseWorker_.setRelaxationIterations(RelaxationIterations);
  baseWorker_.setOutputPoints(&outVertices_);
  baseWorker_.setOutputQuads(&outQuadrangles_);

  auto quads = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto mesh = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(mesh);

  TTK_ABORT_KK(res != 0, "Cannot get mesh triangulation", -1);

  res += getQuadVertices(quads);

  TTK_ABORT_KK(res != 0, "Cannot get quad vertices", -2);

  res += baseWorker_.execute();

  TTK_ABORT_KK(
    res != 0, "QuadrangulationSubdivision.execute() error code: " << res, -3);

  auto outArrayQuads = vtkSmartPointer<vtkIdTypeArray>::New();
  auto cells = vtkSmartPointer<vtkCellArray>::New();

  // outArray points to outQuadrangles_ data, but does not deallocate it
  outArrayQuads->SetArray(outQuadrangles_.data(), outQuadrangles_.size(), 1);

  // vtkCellArray of quadrangle values containing outArray
  cells->SetCells(outQuadrangles_.size(), outArrayQuads);

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < outVertices_.size(); i+=3) {
    points->InsertNextPoint(&outVertices_[i]);
  }

  // update output: get quadrangle vertices
  output->SetPoints(points);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
