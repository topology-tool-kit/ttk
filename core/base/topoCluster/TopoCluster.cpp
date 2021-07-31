#include <TopoCluster.h>

using namespace std;
using namespace ttk;

TopoCluster::TopoCluster() {
  setDebugMsgPrefix("TopoCluster");
  clear();

#ifdef TTK_ENABLE_OPENMP
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
#else
  caches_.resize(1);
  cacheMaps_.resize(1);
#endif
  initCache(0.2f);
}

TopoCluster::~TopoCluster() {
}

int TopoCluster::clear() {

  AbstractTriangulation::clear();

  vertexNumber_ = 0;
  cellNumber_ = 0;
  nodeNumber_ = 0;
  doublePrecision_ = false;

  // {
  //   stringstream msg;
  //   msg << "[TopoCluster] Triangulation cleared." << endl;
  //   dMsg(cout, msg.str(), detailedInfoMsg);
  // }

  return AbstractTriangulation::clear();
}