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

   int scalarFieldId = 0;

   program.parser_.setArgument("f", &scalarFieldId, "Scalar field id", true);

   int ret = 0;
   ret     = program.init(argc, argv);

   if (ret != 0)
      return ret;

   program.ttkObject_->SetScalarFieldId(scalarFieldId);

   // execute data processing
   ret = program.run();

   if (ret != 0)
      return ret;

   ret = program.save();

   return ret;
}
