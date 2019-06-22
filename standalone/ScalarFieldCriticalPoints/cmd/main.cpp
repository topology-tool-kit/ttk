/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for critical point computation.

// include the local headers
#include <ttkProgramBase.h>
#include <ttkScalarFieldCriticalPoints.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkScalarFieldCriticalPoints> program;

  // specify local parameters to the TTK module with default values.
  int scalarFieldId = 0, offsetFieldId = -1;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "F", &scalarFieldId, "Input scalar field identifier", true);
  program.parser_.setArgument(
    "O", &offsetFieldId, "Input vertex offset field identifier", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetScalarFieldId(scalarFieldId);
  program.ttkObject_->SetOffsetFieldId(offsetFieldId);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
