/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>.
/// \date January 2019.
///
/// \brief FTM algorithm standalone binary.

// include the local headers
#include <ttkFTMTree.h>
#include <ttkProgramBase.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  vtkProgram<ttkFTMTree> program;

  int fieldId = 0;
  int treeType = 0;

  program.parser_.setArgument("f", &fieldId, "Field identifier", true);
  program.parser_.setArgument(
    "T", &treeType, "Tree type {0: JT, 1: ST, 2: CT}", true);

  int ret = 0;
  ret = program.init(argc, argv);

  if(ret != 0)
    return ret;

  program.ttkObject_->SetScalarFieldId(fieldId);
  program.ttkObject_->SetTreeType(treeType);

  // execute data processing
  ret = program.run();

  if(ret != 0)
    return ret;

  ret = program.save();

  return ret;
}
