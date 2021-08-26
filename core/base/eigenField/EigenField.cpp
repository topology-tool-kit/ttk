#include <EigenField.h>
#include <Laplacian.h>

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

// main routine
template <typename T, class TriangulationType>
int ttk::EigenField::execute(const TriangulationType &triangulation,
                             T *const outputFieldPointer,
                             const unsigned int eigenNumber,
                             bool computeStatistics,
                             T *const outputStatistics) const {

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)

  Timer tm;
  Memory mem;

  this->printMsg("Beginning computation...");

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  using SpMat = Eigen::SparseMatrix<T>;
  using DMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  // number of vertices
  const auto vertexNumber = triangulation.getNumberOfVertices();

  // graph laplacian of current mesh
  SpMat lap;
  // compute graph laplacian using cotangent weights
  Laplacian::cotanWeights<T>(lap, triangulation);
  // lap is square
  eigen_plain_assert(lap.cols() == lap.rows());

  auto n = lap.cols();
  auto m = eigenNumber;
  // threshold: minimal number of eigenpairs to get a converging solution
  const size_t minEigenNumber = 20;

  if(eigenNumber == 0) {
    // default value
    m = n / 1000;
  } else if(eigenNumber < minEigenNumber) {
    m = minEigenNumber;
  }

  Spectra::SparseSymMatProd<T> op(lap);
  Spectra::SymEigsSolver<decltype(op)> solver(op, m, 2 * m);

  solver.init();

  // number of eigenpairs correctly computed
  int nconv = solver.compute(Spectra::SortRule::LargestAlge);

  switch(solver.info()) {
    case Spectra::CompInfo::NumericalIssue:
      this->printMsg("Numerical Issue!", ttk::debug::Priority::ERROR);
      break;
    case Spectra::CompInfo::NotConverging:
      this->printMsg("No Convergence! (" + std::to_string(nconv) + " out of "
                       + std::to_string(eigenNumber) + " values computed)",
                     ttk::debug::Priority::ERROR);
      break;
    case Spectra::CompInfo::NotComputed:
      this->printMsg("Invalid Input!", ttk::debug::Priority::ERROR);
      break;
    default:
      break;
  }

  DMat eigenvectors = solver.eigenvectors();

  auto outputEigenFunctions = static_cast<T *>(outputFieldPointer);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    for(size_t j = 0; j < eigenNumber; ++j) {
      // cannot avoid copy here...
      outputEigenFunctions[i * eigenNumber + j] = eigenvectors(i, j);
    }
  }

  if(computeStatistics && outputStatistics != nullptr) {

    // number of statistics components
    const int statsComp = 4;

    auto outputStats = static_cast<T *>(outputStatistics);

    // compute statistics on eigenfunctions
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber; ++i) {
      const auto k = i * statsComp;
      // init current tuple computation
      outputStats[k] = outputEigenFunctions[i * eigenNumber];
      outputStats[k + 1] = outputEigenFunctions[i * eigenNumber];
      outputStats[k + 2] = outputEigenFunctions[i * eigenNumber];
      outputStats[k + 3] = outputEigenFunctions[i * eigenNumber]
                           * outputEigenFunctions[i * eigenNumber];
      // loop from 1
      for(size_t j = 1; j < eigenNumber; ++j) {
        outputStats[k] = std::min<T>(
          outputEigenFunctions[i * eigenNumber + j], outputStats[k]);
        outputStats[k + 1] = std::max<T>(
          outputEigenFunctions[i * eigenNumber + j], outputStats[k]);
        outputStats[k + 2] += outputEigenFunctions[i * eigenNumber + j];
        outputStats[k + 3] += outputEigenFunctions[i * eigenNumber + j]
                              * outputEigenFunctions[i * eigenNumber + j];
        ;
      }
    }
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_,
                 mem.getElapsedUsage());

#else
  TTK_FORCE_USE(triangulation);
  TTK_FORCE_USE(outputFieldPointer);
  TTK_FORCE_USE(eigenNumber);
  TTK_FORCE_USE(computeStatistics);
  TTK_FORCE_USE(outputStatistics);

  this->printMsg(
    std::vector<std::string>{"Spectra support disabled, computation skipped!",
                             "Please re-compile TTK with Eigen AND Spectra "
                             "support to enable this feature."},
    debug::Priority::ERROR);
#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  return 0;
}

// explicit template specializations for double and float types
#define EIGENFIELD_SPECIALIZE(TYPE)                                            \
  template int ttk::EigenField::execute<TYPE>(                                 \
    const Triangulation &, TYPE *const, const unsigned int, bool, TYPE *const) \
    const

EIGENFIELD_SPECIALIZE(double);
EIGENFIELD_SPECIALIZE(float);
