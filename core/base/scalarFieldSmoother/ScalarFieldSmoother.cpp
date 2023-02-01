#include <ScalarFieldSmoother.h>

ttk::ScalarFieldSmoother::ScalarFieldSmoother() {
  this->setDebugMsgPrefix("ScalarFieldSmoother");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}

ttk::ScalarFieldSmoother::~ScalarFieldSmoother() = default;
