#include <VPath.h>

using namespace std;
using namespace ttk;
using namespace vp;

VPath::VPath(){
  this->setDebugMsgPrefix("VPath");
}

VPath::~VPath() = default;

int VPath::execute(vector<dcg::Cell> &output, const bool &isForward){

  printMsg("Computing VPath...");

  return 0;
}
