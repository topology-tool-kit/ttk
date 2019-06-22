/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2018.
///
/// \brief point merger program.

// include the local headers
#include <ttkPointMerger.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkPointMerger> program;

class myKeyHandler : public ttkKeyHandler {

public:
  int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key) {

    stringstream msg;
    msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
    dMsg(cout, msg.str(), infoMsg);

    return 0;
  }
};

int main(int argc, char **argv) {

  bool boundaryOnly = true;
  double distanceThreshold = 0.001;

  program.parser_.setOption("b", &boundaryOnly, "Boundary only");
  program.parser_.setArgument(
    "t", &distanceThreshold, "Distance threshold", true);

  int ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetBoundaryOnly(boundaryOnly);
  program.ttkObject_->SetDistanceThreshold(distanceThreshold);

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler);

  program.run();

  return 0;
}
