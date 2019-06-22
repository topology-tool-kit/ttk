/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program example for scalar field smoothing.

// include the local headers
#include <ttkProgramBase.h>
#include <ttkScalarFieldSmoother.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkScalarFieldSmoother> program;

  // specify local parameters to the TTK module with default values.
  int iterationNumber = 1, scalarFieldId = 0, maskId = -1;

  // register these arguments to the command line parser
  program.parser_.setArgument("I", &iterationNumber, "Iteration number", true);
  program.parser_.setArgument(
    "F", &scalarFieldId, "Scalar field identifier", true);
  program.parser_.setArgument("M", &maskId, "Mask field identifier", true);

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetNumberOfIterations(iterationNumber);
  program.ttkObject_->SetScalarFieldIdentifier(scalarFieldId);
  if(maskId != -1) {
    program.ttkObject_->SetForceInputMaskScalarField(true);
    program.ttkObject_->SetMaskIdentifier(maskId);
  }

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
