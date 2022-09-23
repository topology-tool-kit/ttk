#include <SurfaceGeometrySmoother.h>

ttk::SurfaceGeometrySmoother::SurfaceGeometrySmoother() {
  this->setDebugMsgPrefix("SurfaceGeometrySmoother");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
