/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkFTRGraph.h>
#include <ttkProgramBase.h>

int main(int argc, char **argv) {
  vtkProgram<ttkFTRGraph> program;

  int scalarFieldId = 0;
  bool singleSweep = false;

  program.parser_.setArgument("f", &scalarFieldId, "Scalar field id", true);
  program.parser_.setOption("s", &singleSweep, "Single sweep");

  auto ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetScalarFieldId(scalarFieldId);
  program.ttkObject_->SetSingleSweep(singleSweep);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  ret = program.save();

  return ret;
}
