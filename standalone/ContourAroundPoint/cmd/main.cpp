/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkContourAroundPoint.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkContourAroundPoint> program;

  // specify local parameters to the TTK module with default values.
  double ui_sizeFilter = 7.;
  double ui_extension = 67.;
  //  string ui_scalars = "tas";

  // register these arguments to the command line parser
  program.parser_.setArgument("S", &ui_sizeFilter, "Size filter", true);
  program.parser_.setArgument("E", &ui_extension, "Extension", true);
  //  program.parser_.setOption("D", &ui_scalars,
  //    "Scalar variable");

  int ret = 0;
  ret = program.init(argc, argv);
  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->Setui_sizeFilter(ui_sizeFilter);
  program.ttkObject_->Setui_extension(ui_extension);
  //  program.ttkObject_->Setui_scalars(ui_scalars);

  // execute data processing
  ret = program.run();
  if(ret != 0)
    return ret;

  // save the output
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();

  return ret;
}
