#include <ManifoldCheck.h>

using namespace std;
using namespace ttk;

ManifoldCheck::ManifoldCheck() {

  vertexLinkComponentNumber_ = nullptr;
  edgeLinkComponentNumber_ = nullptr;
  triangleLinkComponentNumber_ = nullptr;
  this->setDebugMsgPrefix("ManifoldCheck");
}

ManifoldCheck::~ManifoldCheck() = default;
