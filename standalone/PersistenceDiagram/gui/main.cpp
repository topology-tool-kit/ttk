/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief GUI program for persistence diagram computation.

// include the local headers
#include <ttkPersistenceDiagram.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkPersistenceDiagram> program;

int main(int argc, char **argv) {

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

  program.run();

  return 0;
}
