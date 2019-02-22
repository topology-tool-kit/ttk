#include <HarmonicFieldComputation.h>

ttk::HarmonicFieldComputation::HarmonicFieldComputation()
    : vertexNumber_{}, constraintNumber_{}, useCotanMethod_{false},
      triangulation_{}, constraints_{}, outputScalarFieldPointer_{},
      vertexIdentifierScalarFieldPointer_{} {
#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_OPENMP)
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP && TTK_ENABLE_EIGEN
}
