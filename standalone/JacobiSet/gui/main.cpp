/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief GUI program for Jacobi set computation.

// include the local headers
#include <ttkJacobiSet.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkJacobiSet> program;

int main(int argc, char **argv) {

  // specify local parameters to the TTK module with default values.
  int uComponentId = 0, vComponentId = 1;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "u", &uComponentId, "Identifier of the u-component field", true);
  program.parser_.setArgument(
    "v", &vComponentId, "Identifier of the v-component field", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetUcomponentId(uComponentId);
  program.ttkObject_->SetVcomponentId(vComponentId);

  program.run();

  return 0;
}
