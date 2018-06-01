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
   bool diffCode = false;
   program.parser_.setOption("m", &meshOnly, "Only compare mesh, no scalars");
   program.parser_.setOption("d", &diffCode, "simply return 0 if data sets are equals");

   int ret = 0;
   ret     = program.init(argc, argv);

   if (ret != 0)
      return ret;

   program.ttkObject_->SetmeshOnly(meshOnly);

   ret = program.run();

   if (ret != 0)
      return ret;

   if (diffCode) {
      ret = program.ttkObject_->getDiffCode();
   } else {
      ret = program.save();
   }

   return ret;
}
