#include <CompactTriangulation.h>

using namespace std;
using namespace ttk;

CompactTriangulation::CompactTriangulation() {
  setDebugMsgPrefix("CompactTriangulation");
  clear();
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
}

CompactTriangulation::~CompactTriangulation() {
}

int CompactTriangulation::clear() {

  vertexNumber_ = 0;
  cellNumber_ = 0;
  nodeNumber_ = 0;
  doublePrecision_ = false;

  return AbstractTriangulation::clear();
}