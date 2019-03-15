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

int ttkQuadrangulationSubdivision::getTriangulation(vtkDataSet *const input) {

  auto trg = ttkTriangulation::getTriangulation(input);

  TTK_ABORT_KK(trg == nullptr, "input triangulation pointer is NULL.", -1);

  trg->setWrapper(this);
  quadrangulationSubdivision_.setWrapper(this);
  quadrangulationSubdivision_.setupTriangulation(trg);
  Modified();

  TTK_ABORT_KK(trg->isEmpty(), "ttkTriangulation allocation problem.", -2);

  return 0;
}

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  quadrangulationSubdivision_.setSubdivisionLevel(SubdivisionLevel);
  quadrangulationSubdivision_.setRelaxationIterations(RelaxationIterations);

  auto quads = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto mesh = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(mesh);

  TTK_ABORT_KK(
    res != 0, "Cannot get mesh triangulation", -3);

  res += quadrangulationSubdivision_.execute();

  TTK_ABORT_KK(
    res != 0, "QuadrangulationSubdivision.execute() error code: " << res, -4);

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
