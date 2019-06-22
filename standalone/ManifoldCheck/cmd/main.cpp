/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2018.
///
/// \brief Manifold check program example.

// include the local headers
#include <ttkManifoldCheck.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkManifoldCheck> program;

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  // optional TODO-4:
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();
  /// end of optional TODO-4

  return ret;
}
