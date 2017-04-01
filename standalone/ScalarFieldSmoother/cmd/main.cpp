/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program example for scalar field smoothing.

// include the local headers
#include                  <vtkScalarFieldSmoother.h>
#include                  <vtkProgramBase.h>

int main(int argc, char **argv) {

  vtkProgram<vtkScalarFieldSmoother> program;
  
  // specify local parameters to the TTK module with default values.
  int iterationNumber = 1, scalarFieldId = 0;

  // register these arguments to the command line parser
  program.parser_.setArgument("I", &iterationNumber,
    "Iteration number", true);
  program.parser_.setArgument("F", &scalarFieldId,
    "Scalar field identifier", true);
  
  int ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetNumberOfIterations(iterationNumber);
  program.ttkObject_->SetScalarFieldIdentifier(scalarFieldId);
  
  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
  
  // save the output
  ret = program.save();
  
  return ret;
}
