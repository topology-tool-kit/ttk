#include <BaseClass.h>

#ifdef TTK_ENABLE_OPENMP
COMMON_EXPORTS int ttk::globalThreadNumber_ = omp_get_max_threads();
#else
COMMON_EXPORTS int ttk::globalThreadNumber_ = 1;
#endif

using namespace ttk;

BaseClass::BaseClass() : lastObject_{false}, wrapper_{nullptr} {
  threadNumber_ = ttk::globalThreadNumber_;
}

int BaseClass::setWrapper(const Wrapper *wrapper) {
  wrapper_ = const_cast<Wrapper *>(wrapper);

  // dynamic_cast is impossible because the Wrapper class is incomplete
  setThreadNumber(reinterpret_cast<BaseClass *>(wrapper_)->threadNumber_);
  return 0;
}
