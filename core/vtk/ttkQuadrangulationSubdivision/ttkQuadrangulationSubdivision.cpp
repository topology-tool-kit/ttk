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
    ForceInputOffsetIdentifiersField{false}, SubdivisionLevel{5},
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

  auto pointData = input->GetPointData();
  auto identifiers = pointData->GetArray(
    static_cast<const char *>(ttk::VertexScalarFieldName));

  TTK_ABORT_KK(pointData == nullptr, "invalid input quadrangle point data", -4);
  TTK_ABORT_KK(identifiers == nullptr,
               "invalid input quadrangle vertices identifiers", -4);

  baseWorker_.setInputQuadranglesNumber(cells->GetNumberOfCells());
  baseWorker_.setInputQuadrangles(cells->GetData()->GetVoidPointer(0));
  baseWorker_.setInputQuadIdentifiers(identifiers->GetVoidPointer(0));

  return 0;
}

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  baseWorker_.setSubdivisionLevel(SubdivisionLevel);
  baseWorker_.setRelaxationIterations(RelaxationIterations);
  baseWorker_.setOutputQuads(&outQuadrangles_);

  auto quads = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto mesh = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(mesh);

  TTK_ABORT_KK(res != 0, "Cannot get mesh triangulation", -1);

  res += getQuadVertices(quads);

  TTK_ABORT_KK(res != 0, "Cannot get quad vertices", -2);

  // output points
  auto points = vtkSmartPointer<vtkPoints>::New();
  // total number of points after subdivision:
  // critical points +  MSC edges * (5 ^ SubdivisionLevel),
  // given that MSC edges = 5 values / quad
  // and that each subdivision generates 5 new points
  vtkIdType numPoints = quads->GetNumberOfPoints()
                        + quads->GetNumberOfCells() / 5
                            * static_cast<vtkIdType>(pow(5, SubdivisionLevel));
  points->SetNumberOfPoints(numPoints);
  // pass to base worker
  baseWorker_.setOutputPointNumber(numPoints);
  baseWorker_.setOutputPoints(points->GetVoidPointer(0));

  res += baseWorker_.execute();

  TTK_ABORT_KK(
    res != 0, "QuadrangulationSubdivision.execute() error code: " << res, -3);

  auto outArray = vtkSmartPointer<vtkIdTypeArray>::New();
  auto cells = vtkSmartPointer<vtkCellArray>::New();

  // outArray points to outQuadrangles_ data, but does not deallocate it
  outArray->SetArray(outQuadrangles_.data(), outQuadrangles_.size(), 1);

  // vtkCellArray of quadrangle values containing outArray
  cells->SetCells(outQuadrangles_.size(), outArray);

  // update output: get quadrangle values
  output->SetCells(VTK_QUAD, cells);

  // update output: get quadrangle vertices
  output->SetPoints(points);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
