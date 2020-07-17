/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for mesh subdivision.

// include the local headers
#include <ttkMeshSubdivision.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkMeshSubdivision> program;

  // specify local parameters to the TTK module with default values.
  int iterationNumber = 1;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "I", &iterationNumber, "Number of subdivision iterations", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetIterationNumber(iterationNumber);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
