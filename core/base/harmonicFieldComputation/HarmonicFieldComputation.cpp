#include <HarmonicFieldComputation.h>

ttk::HarmonicFieldComputation::HarmonicFieldComputation()
    : triangulation_{}, vertexNumber_{}, constraintNumber_{},
      inputScalarFieldPointer_{}, vertexIdentifierScalarFieldPointer_{},
      inputOffsetScalarFieldPointer_{}, considerIdentifierAsBlackList_{},
      addPerturbation_{}, outputScalarFieldPointer_{},
      outputOffsetScalarFieldPointer_{} {
  considerIdentifierAsBlackList_ = false;
  addPerturbation_ = false;
#ifdef TTK_ENABLE_OPENMP
  Eigen::SetNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP
}

ttk::HarmonicFieldComputation::~HarmonicFieldComputation() {}
