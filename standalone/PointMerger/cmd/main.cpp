/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2018.
///
/// \brief point merger program.

// include the local headers
#include <ttkPointMerger.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkPointMerger> program;

  bool boundaryOnly = true;
  double distanceThreshold = 0.001;

  program.parser_.setOption("b", &boundaryOnly, "Boundary only");
  program.parser_.setArgument(
    "t", &distanceThreshold, "Distance threshold", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetBoundaryOnly(boundaryOnly);
  program.ttkObject_->SetDistanceThreshold(distanceThreshold);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  ret = program.save();

  return ret;
}
