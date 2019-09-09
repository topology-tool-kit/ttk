/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for bivariate Reeb space computation.

// include the local headers
#include <ttkReebSpace.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkReebSpace> program;

int main(int argc, char **argv) {

  // specify local parameters to the TTK module with default values.
  double threshold = 0;
  int criterion = 0;
  int uId = 0, vId = 1;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "s", &threshold, "Simplification threshold, in [0, 1] (default: 0)", true);
  program.parser_.setArgument(
    "c", &criterion, "Simplification criterion, in [0, 2] (default: 1)", true);

  program.parser_.setArgument(
    "u", &uId, "Identifier of the u-component field (default: 0)", true);
  program.parser_.setArgument(
    "v", &vId, "Identifier of the v-component field (default: 1)", true);

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetUcomponentId(uId);
  program.ttkObject_->SetVcomponentId(vId);
  program.ttkObject_->SetSimplificationCriterion(criterion);
  program.ttkObject_->SetSimplificationThreshold(threshold);

  vector<int> hiddenOutputs = {0, 2, 3};
  program.hideOutputs(hiddenOutputs);

  // execute data processing
  ret = program.run();

  return ret;
}
