#include <ManifoldCheck.h>

using namespace std;
using namespace ttk;

ManifoldCheck::ManifoldCheck() {

  vertexLinkComponentNumber_ = NULL;
  edgeLinkComponentNumber_ = NULL;
  triangleLinkComponentNumber_ = NULL;
  this->setDebugMsgPrefix("ManifoldCheck");
}

ManifoldCheck::~ManifoldCheck() {
}
