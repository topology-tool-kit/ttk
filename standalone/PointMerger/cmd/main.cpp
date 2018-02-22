/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2018.
///
/// \brief point merger program.

// include the local headers
#include                  <ttkPointMerger.h>
#include                  <ttkProgramBase.h>

int main(int argc, char **argv) {

  vtkProgram<ttkPointMerger> program;
  
  double distanceThreshold = 0.000001;

  program.parser_.setArgument("D", &distanceThreshold,
    "Distance threshold", true);
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  program.ttkObject_->SetDistanceThreshold(distanceThreshold);
  
  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  ret = program.save();
  
  return ret;
}
