/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Continuous scatterplot computation program.

// include the local headers
#include <ttkContinuousScatterPlot.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkContinuousScatterPlot> program;

  // specify local parameters to the TTK module with default values.
  int uId = 0, vId = 1;
  int xRes = 1920, yRes = 1080;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "u", &uId, "Identifier of the u-component field", true);
  program.parser_.setArgument(
    "v", &vId, "Identifier of the v-component field", true);
  program.parser_.setArgument("x", &xRes, "Width of the scatterplot", true);
  program.parser_.setArgument("y", &yRes, "Height of the scatterplot", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetUcomponentId(uId);
  program.ttkObject_->SetVcomponentId(vId);
  program.ttkObject_->SetScatterplotResolution(xRes, yRes);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
