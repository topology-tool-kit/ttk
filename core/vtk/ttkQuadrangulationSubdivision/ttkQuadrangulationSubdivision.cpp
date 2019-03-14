#include <ttkQuadrangulationSubdivision.h>

#define MODULE_S "[ttkQuadrangulationSubdivision] "
#define MODULE_ERROR_S MODULE_S "Error: "

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

int ttkQuadrangulationSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  quadrangulationSubdivision_.setSubdivisionLevel(SubdivisionLevel);
  quadrangulationSubdivision_.setRelaxationIterations(RelaxationIterations);

  auto quads = vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  auto mesh = vtkUnstructuredGrid::SafeDownCast(inputs[1]);
  auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int res = 0;

  quadrangulationSubdivision_.setOutputCells(&outQuadrangles_);

  res += quadrangulationSubdivision_.execute();

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/QuadrangulationSubdivision
  if(res != 0) {
    cerr << MODULE_S "QuadrangulationSubdivision.execute() error code: " << res
         << endl;
    return -9;
  }
#endif

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
