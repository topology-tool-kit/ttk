#include <CompactTriangulation.h>

using namespace std;
using namespace ttk;

CompactTriangulation::CompactTriangulation() {
  setDebugMsgPrefix("CompactTriangulation");
  clear();

#ifdef TTK_ENABLE_OPENMP
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
#else
  caches_.resize(1);
  cacheMaps_.resize(1);
#endif
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