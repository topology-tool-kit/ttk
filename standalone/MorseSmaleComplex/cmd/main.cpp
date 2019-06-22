/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for Morse-Smale complex computation.

// include the local headers
#include <ttkMorseSmaleComplex.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkMorseSmaleComplex> program;

  // specify local parameters to the TTK module with default values.
  int scalarFieldId = 0, offsetFieldId = -1;
  bool plCompliantExtrema = true, plCompliantSaddles = false;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "F", &scalarFieldId, "Input scalar field identifier", true);
  program.parser_.setArgument(
    "O", &offsetFieldId, "Input vertex offset field identifier", true);
  program.parser_.setOption("plE", &plCompliantExtrema, "PL-compliant extrema");
  program.parser_.setOption("plS", &plCompliantSaddles, "PL-compliant saddles");

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetScalarFieldId(scalarFieldId);
  program.ttkObject_->SetOffsetFieldId(offsetFieldId);
  program.ttkObject_->SetReverseSaddleMaximumConnection(plCompliantExtrema);
  program.ttkObject_->SetReverseSaddleSaddleConnection(plCompliantSaddles);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
