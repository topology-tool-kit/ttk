/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2018.
///
/// \brief Manifold check GUI program example.

// include the local headers
#include <ttkManifoldCheck.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkManifoldCheck> program;

class myKeyHandler : public ttkKeyHandler {

public:
  int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key) {

    stringstream msg;
    msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
    dMsg(cout, msg.str(), infoMsg);

    // TODO-4
    // depending on the value of "key", trigger the right functions on the
    // program object (or its contained ttkObject_).
    // end of TODO-4

    return 0;
  }
};

int main(int argc, char **argv) {

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler);

  // TODO-5
  // specify if the display of certain outputs should be disabled
  // vector<int> hiddenOutputs = {0, 2};
  // program.hideOutputs(hiddenOutputs);
  // end of TODO-5

  program.run();

  return 0;
}
