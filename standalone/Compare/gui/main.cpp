/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy GUI program example.

// include the local headers
#include <ttkCompare.h>
#include <ttkUserInterfaceBase.h>

using namespace std;
using namespace ttk;

vtkUserInterface<ttkCompare> program;

class myKeyHandler : public ttkKeyHandler
{
  public:
   int OnKeyPress(vtkRenderWindowInteractor *interactor, string &key)
   {
      stringstream msg;
      msg << "[myKeyHandler] The user pressed the key `" << key << "'." << endl;
      dMsg(cout, msg.str(), infoMsg);

      return 0;
   }
};

int main(int argc, char **argv)
{
   bool meshOnly = false;
   program.parser_.setOption("m", &meshOnly, "Some option to enable or disable");

   int ret = program.init(argc, argv);

   if (ret != 0)
      return ret;

   program.ttkObject_->SetmeshOnly(meshOnly);

   myKeyHandler myHandler;
   program.setKeyHandler(&myHandler);

   program.run();

   return 0;
}
