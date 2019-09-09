/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2017-06-30
///
/// \brief GUI program for contour tree computation.

// include the local headers
#include <ttkFTMTree.h>
#include <ttkUserInterfaceBase.h>

vtkUserInterface<ttkFTMTree> program;

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {
  // specify local parameters to the TTK module with default values.
  int scalarFieldId = 0, offsetFieldId = -1, treeType = 2, arcSampling = 20;

  // register these arguments to the command line parser
  program.parser_.setArgument(
    "F", &scalarFieldId, "Scalar field identifier", true);
  program.parser_.setArgument(
    "O", &offsetFieldId, "Field identifier for vertex offsets", true);
  program.parser_.setArgument(
    "T", &treeType, "Tree type (0: join, 1: split, 2: contour)", true);
  program.parser_.setArgument("Sa", &arcSampling, "Arc sampling", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
  program.ttkObject_->SetScalarFieldId(scalarFieldId);
  program.ttkObject_->SetOffsetFieldId(offsetFieldId);
  program.ttkObject_->SetTreeType(treeType);
  program.ttkObject_->SetSuperArcSamplingLevel(arcSampling);
  globalDebugLevel_ = 2;

  vector<int> hiddenOutputs = {0, 2};
  program.hideOutputs(hiddenOutputs);
  program.run();

  return 0;
}
