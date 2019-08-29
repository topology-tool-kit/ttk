/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkBlank.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkBlank> program;

  // TODO-1:
  // specify local parameters to the TTK module with default values.
  bool someOption = false;
  int someIntegerArgument = -1;
  double someDoubleArgument = -1.0;
  // end of TODO-1

  // TODO-2:
  // register these arguments to the command line parser
  program.parser_.setArgument(
    "D", &someDoubleArgument, "Some optional double argument", true);
  program.parser_.setArgument(
    "I", &someIntegerArgument, "Some optional integer argument", true);
  program.parser_.setOption(
    "O", &someOption, "Some option to enable or disable");
  // end of TODO-2

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // TODO-3:
  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetSomeIntegerArgument(someIntegerArgument);
  program.ttkObject_->SetSomeDoubleArgument(someDoubleArgument);
  program.ttkObject_->SetSomeOption(someOption);
  // end of TODO-3

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  // optional TODO-4:
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();
  /// end of optional TODO-4

  return ret;
}
