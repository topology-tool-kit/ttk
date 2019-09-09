/// \author Maxime Soler <soler.maxime@total.com>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include <ttkBottleneckDistance.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkBottleneckDistance> program;

class myKeyHandler : public ttkKeyHandler {

public:
  int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key) {

    stringstream msg;
    msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
    dMsg(cout, msg.str(), infoMsg);

    // depending on the value of "key", trigger the right functions on the
    // program object (or its contained ttkObject_).

    return 0;
  }
};

int main(int argc, char **argv) {

  string wasserstein;
  double alpha = 1.0;
  double spacing = 5.0;
  bool useMatchingMesh = false;

  bool persistence = false;

  program.parser_.setArgument(
    "w", &wasserstein, "Wasserstein type (inf = bottleneck, 0, 1, 2, ...)");

  program.parser_.setArgument("a", &alpha, "Geometric factor.", true);

  program.parser_.setArgument("s", &spacing, "Spacing (absolute).", true);

  program.parser_.setOption("m", &useMatchingMesh, "Use matching mesh.");
  program.parser_.setOption("p", &persistence, "Use persistence metric");

  int ret = program.init(argc, argv);
  if(ret != 0)
    return ret;

  program.ttkObject_->SetWassersteinMetric(wasserstein);
  program.ttkObject_->SetUsePersistenceMetric(persistence);
  program.ttkObject_->SetUseOutputMatching(useMatchingMesh);
  program.ttkObject_->SetAlpha(alpha);
  program.ttkObject_->SetSpacing(spacing);

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler);

  // specify if the display of certain outputs should be disabled
  // vector<int> hiddenOutputs = {0, 2};
  // program.hideOutputs(hiddenOutputs);

  program.run();

  return ret;
}
