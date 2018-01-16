/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <ttkFTRGraph.h>
#include <ttkProgramBase.h>

int main(int argc, char **argv)
{
   vtkProgram<ttkFTRGraph> program;

   int ret = 0;
   ret     = program.init(argc, argv);

   if (ret != 0)
      return ret;

   // execute data processing
   ret = program.run();

   if (ret != 0)
      return ret;

   ret = program.save();

   return ret;
}
