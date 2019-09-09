/// \author Julien Tierny <julien.tierny@lip6.fr>.
/// \date February 2017.
///
/// \brief Command line program for contour tree computation.

// include the local headers
#include <ttkContourForests.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkContourForests> program;

  // specify local parameters to the TTK module with default values.
  int scalarFieldId = 0, offsetFieldId = -1, treeType = 2, arcSampling = 20,
      arcSmoothing = 15;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "F", &scalarFieldId, "Scalar field identifier", true);
  program.parser_.setArgument(
    "O", &offsetFieldId, "Field identifier for vertex offsets", true);
  program.parser_.setArgument(
    "T", &treeType, "Tree type (0: join, 1: split, 2: contour)", true);
  program.parser_.setArgument("Sa", &arcSampling, "Arc sampling", true);
  program.parser_.setArgument(
    "So", &arcSmoothing, "Iteration number for arc smoothing", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetFieldId(scalarFieldId);
  program.ttkObject_->SetInputOffsetFieldId(offsetFieldId);
  program.ttkObject_->SetTreeType(treeType);
  program.ttkObject_->SetArcResolution(arcSampling);
  program.ttkObject_->SetSkeletonSmoothing(arcSmoothing);
  // globalDebugLevel_ = 2;

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  // save the output
  ret = program.save();

  return ret;
}
