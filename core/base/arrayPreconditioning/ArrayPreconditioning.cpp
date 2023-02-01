#include <ArrayPreconditioning.h>

ttk::ArrayPreconditioning::ArrayPreconditioning() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ArrayPreconditioning");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
