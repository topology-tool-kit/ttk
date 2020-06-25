#include <ScalarFieldCriticalPoints.h>

ttk::ScalarFieldCriticalPoints::ScalarFieldCriticalPoints() {

  dimension_ = 0;
  vertexNumber_ = 0;
  vertexLinkEdgeLists_ = NULL;
  criticalPoints_ = NULL;
  sosOffsets_ = NULL;

  forceNonManifoldCheck = false;

  setDebugMsgPrefix("ScalarFieldCriticalPoints");

  //   threadNumber_ = 1;
}

ttk::ScalarFieldCriticalPoints::~ScalarFieldCriticalPoints() {
}
