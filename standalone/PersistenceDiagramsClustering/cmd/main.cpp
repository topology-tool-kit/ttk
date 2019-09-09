/// \brief Command line program for persistence diagram barycenters computation.

// include the local headers
#include <ttkPersistenceDiagramsClustering.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkPersistenceDiagramsClustering> program;

  // specify local parameters to the TTK module with default values.
  double timeLimit = -1;
  int numberOfInputs = 1;
  int numberOfClusters = 1;
  double geometry_penalization = 1.;
  int kmeanspp = 1;
  int accelerated = 1;
  int pairType = -1;
  int randomness = 0;
  int method = 0;
  int use_prog = 1;
  int write_distances = 0;

  // register these arguments to the command line parser
  program.parser_.setArgument("M", &method, "Select algorithm", true);
  program.parser_.setArgument("U", &use_prog, "use progressivity", true);
  program.parser_.setArgument(
    "P", &pairType,
    "Set the type of critical pairs considered in the clustering. \n\t\t\t-1 : "
    "all critical pairs\n\t\t\t 0 : only min-saddle pairs \n\t\t\t 1 : "
    "onlysaddle-saddle pairs\n\t\t\t 2 : only saddle-max pairs\n\t\t",
    true);
  program.parser_.setArgument(
    "A", &accelerated, "Use the kmean acceleration", true);
  program.parser_.setArgument("I", &kmeanspp, "Use the kmeanspp init.", true);
  program.parser_.setArgument(
    "K", &numberOfClusters, "Number of Clusters", true);
  program.parser_.setArgument(
    "T", &timeLimit,
    "Time Limit for Computation, in seconds. No time limit by default", true);
  program.parser_.setArgument("R", &randomness, "Randomness", true);
  program.parser_.setArgument(
    "G", &geometry_penalization,
    "Geometry Penalization, bet 0. and 1., 1. means no lifting", true);
  program.parser_.setArgument(
    "write-distances", &write_distances,
    "0 : don't write - 1 : write accelerated KMeans approximations - 2 : "
    "compute and write actual distances",
    true);

  // program.parser_.printArgs();
  // std::cout<<"number of args "<<program.parser_.getNumberOfArgs()<<std::endl;

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.

  if(timeLimit == -1) {
    timeLimit = 999999999999;
  }
  int deterministic = 1 - randomness;
  program.ttkObject_->SetMethod(method);
  program.ttkObject_->SetUseProgressive(use_prog);
  program.ttkObject_->SetDeterministic(deterministic);
  program.ttkObject_->SetUseAccelerated(accelerated);
  program.ttkObject_->SetPairTypeClustering(pairType);
  program.ttkObject_->SetUseKmeansppInit(kmeanspp);
  program.ttkObject_->SetTimeLimit(timeLimit);
  program.ttkObject_->SetAlpha(geometry_penalization);
  program.ttkObject_->SetNumberOfClusters(numberOfClusters);
  program.ttkObject_->SetDistanceWritingOptions(write_distances);

  program.ttkObject_->setNumberOfInputsFromCommandLine(
    program.getNumberOfInputs());
  // program.ttkObject_->AddInputConnection(0,program.ttkObject_->GetInputAlgorithm(0,0));

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
