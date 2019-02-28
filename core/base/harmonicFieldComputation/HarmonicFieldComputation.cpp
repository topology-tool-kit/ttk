#include <HarmonicFieldComputation.h>

ttk::HarmonicFieldComputation::HarmonicFieldComputation()
    : vertexNumber_{}, constraintNumber_{}, useCotanWeights_{false},
      triangulation_{}, sources_{}, constraints_{},
      outputScalarFieldPointer_{}, solverType_{Cholesky} {
#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_OPENMP)
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP && TTK_ENABLE_EIGEN
}

ttk::SolverType ttk::HarmonicFieldComputation::findBestSolver() const {
  // TODO
  return ttk::SolverType::Cholesky;
}
