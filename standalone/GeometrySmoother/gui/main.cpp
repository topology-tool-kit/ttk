/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief dummy GUI program example for geometry smoothing.

// include the local headers
#include <ttkGeometrySmoother.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkGeometrySmoother> program;

int main(int argc, char **argv) {

  // specify local parameters to the TTK module with default values.
  int iterationNumber = 1, maskId = -1;

  program.parser_.setArgument(
    "I", &iterationNumber, "Number of smoothing iterations", true);
  program.parser_.setArgument("M", &maskId, "Mask field identifier", true);

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetNumberOfIterations(iterationNumber);
  if(maskId != -1) {
    program.ttkObject_->SetForceInputMaskScalarField(true);
    program.ttkObject_->SetMaskIdentifier(maskId);
  }

  program.run();

  return 0;
}
