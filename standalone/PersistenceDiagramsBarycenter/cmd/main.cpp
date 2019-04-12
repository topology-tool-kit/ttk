/// \brief Command line program for persistence diagram barycenters computation.

// include the local headers
#include                  <ttkPersistenceDiagramsBarycenter.h>
#include                  <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkPersistenceDiagramsBarycenter> program;
  
  // specify local parameters to the TTK module with default values.
  double timeLimit = -1;
  int numberOfInputs = 1;
  int method=0;
  double geom=1.0;
  int prog=1;
  // register these arguments to the command line parser

  // program.parser_.setArgument("G", &geom,
  //       "geometry", true);
  program.parser_.setArgument("M", &method,
        "algorithm : 0= prog PB, 2=Auction", true);
  // program.parser_.setArgument("N", &numberOfInputs,
  //       "Number of Input Diagrams", true);
  // program.parser_.setArgument("P", &prog,
  //   "use progressive computation", true);
  program.parser_.setArgument("T", &timeLimit,
    "Time Limit for Computation, in seconds. No time limit by default", true);

  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  // program.ttkObject_->SetAlpha(geom);
  // program.ttkObject_->SetUseProgressive(prog);
  if(timeLimit==-1){
    timeLimit = 999999999999;
  }
  program.ttkObject_->SetMethod(method);
  program.ttkObject_->SetTimeLimit(timeLimit);
  program.ttkObject_->setNumberOfInputsFromCommandLine(program.getNumberOfInputs());
// program.ttkObject_->AddInputConnection(0,program.ttkObject_->GetInputAlgorithm(0,0));
  
  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  // save the output
  ret = program.save();
  
  return ret;
}

