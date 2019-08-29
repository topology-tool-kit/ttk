/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include <ttkFTRGraph.h>
#include <ttkUserInterfaceBase.h>

vtkUserInterface<ttkFTRGraph> program;

class myKeyHandler : public ttkKeyHandler {
public:
  int OnKeyPress(vtkRenderWindowInteractor *interactor, std::string &key) {
    std::stringstream msg;
    msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
    dMsg(cout, msg.str(), infoMsg);

    return 0;
  }
};

int main(int argc, char **argv) {
  int scalarFieldId = 0;
  bool singleSweep = false;

  program.parser_.setArgument("f", &scalarFieldId, "Scalar field id", true);
  program.parser_.setOption("s", &singleSweep, "Single sweep");

  const auto ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetScalarFieldId(scalarFieldId);
  program.ttkObject_->SetSingleSweep(singleSweep);

  myKeyHandler myHandler;
  program.setKeyHandler(&myHandler);

  program.run();

  return 0;
}
