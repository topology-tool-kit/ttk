/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for Jacobi set computation.

// include the local headers
#include <ttkJacobiSet.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkJacobiSet> program;

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

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
