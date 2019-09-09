/// \author Maxime Soler <soler.maxime@total.com>.
/// \date 16/02/2017
////
/// \brief nothing

// include local headers
#include <ttkLDistance.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkLDistance> program;

  string distanceType;

  program.parser_.setArgument(
    "d", &distanceType, "Distance type (1, 2, ..., inf)", true);

  int ret = program.init(argc, argv);
  if(ret != 0)
    return ret;

  // Change here arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetDistanceType(distanceType);

  // Perform data processing.
  ret = program.run();
  if(ret != 0)
    return ret;

  // save the output
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();

  return ret;
}
