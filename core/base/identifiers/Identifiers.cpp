#include <Identifiers.h>

ttk::Identifiers::Identifiers() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("Identifiers");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
#ifdef TTK_ENABLE_MPI
void ttk::Identifiers::updateDebugPrefix() {
  this->setDebugMsgPrefix("Identifiers");
}
#endif
