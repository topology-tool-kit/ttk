#include <ScalarFieldSmoother.h>

ttk::ScalarFieldSmoother::ScalarFieldSmoother() {
  this->setDebugMsgPrefix("ScalarFieldSmoother");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
#ifdef TTK_ENABLE_MPI
void ttk::ScalarFieldSmoother::updateDebugPrefix() {
  this->setDebugMsgPrefix("ScalarFieldSmoother");
}
#endif
ttk::ScalarFieldSmoother::~ScalarFieldSmoother() = default;
