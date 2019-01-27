#include<BaseClass.h>

using namespace ttk;

BaseClass::BaseClass():
  lastObject_{false},
  threadNumber_{1},
  wrapper_{nullptr}
{
#ifdef TTK_ENABLE_OPENMP
  threadNumber_=omp_get_num_procs();
#else
  threadNumber_ = 1;
#endif
}

int BaseClass::setWrapper(const Wrapper *wrapper){
  wrapper_=const_cast<Wrapper*>(wrapper);

  // dynamic_cast is impossible because the Wrapper class is incomplete
  setThreadNumber(reinterpret_cast<BaseClass*>(wrapper_)->threadNumber_);
  return 0;
}
