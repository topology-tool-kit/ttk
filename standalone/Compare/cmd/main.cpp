/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkCompare.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv)
{
   vtkProgram<ttkCompare> program;

   bool meshOnly = false;
   program.parser_.setOption("m", &meshOnly, "Some option to enable or disable");

   int ret = 0;
   ret     = program.init(argc, argv);

   if (ret != 0)
      return ret;

   program.ttkObject_->SetmeshOnly(meshOnly);

   ret = program.run();

   if (ret != 0)
      return ret;

   ret = program.save();

   return ret;
}
