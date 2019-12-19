#include <EigenField.h>
#include <Laplacian.h>

#define MODULE_S "[EigenField] "

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif // __GNUC__
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif // __GNUC__

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

// main routine
template <typename T>
int ttk::EigenField::execute(Triangulation *triangulation,
                             T *const outputFieldPointer,
                             const unsigned int eigenNumber,
                             bool computeStatistics,
                             T *const outputStatistics) const {

  Timer t;
  const auto vertexNumber = triangulation->getNumberOfVertices();

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  using SpMat = Eigen::SparseMatrix<T>;
  using DMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  {
    std::stringstream msg;
    msg << MODULE_S "Beginnning computation..." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  // graph laplacian of current mesh
  SpMat lap;
  // compute graph laplacian using cotangent weights
  Laplacian::cotanWeights<T>(lap, *triangulation);
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
  Spectra::SymEigsSolver<T, Spectra::LARGEST_ALGE, decltype(op)> solver(
    &op, m, 2 * m);

  solver.init();

  // number of eigenpairs correctly computed
  int nconv = solver.compute();

  {
    std::stringstream msg;
    switch(solver.info()) {
      case Spectra::COMPUTATION_INFO::NUMERICAL_ISSUE:
        msg << MODULE_S "Numerical Issue!" << std::endl;
        break;
      case Spectra::COMPUTATION_INFO::NOT_CONVERGING:
        msg << MODULE_S "No Convergence! (" << nconv << " out of "
            << eigenNumber << " values computed)" << std::endl;
        break;
      case Spectra::COMPUTATION_INFO::NOT_COMPUTED:
        msg << MODULE_S "Invalid Input!" << std::endl;
        break;
      default:
        break;
    }
    dMsg(std::cout, msg.str(), infoMsg);
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
      auto k = i * statsComp;
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

  {
    std::stringstream msg;
    msg << MODULE_S "Ending computation after " << t.getElapsedTime() << "s ("
        << threadNumber_ << " thread(s))" << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

#else

  {
    std::stringstream msg;
    msg << MODULE_S << std::endl;
    msg << MODULE_S << std::endl;
    msg << MODULE_S "Spectra support disabled, computation skipped!"
        << std::endl;
    msg << MODULE_S "Please re-compile TTK with Eigen AND Spectra support to "
                    "enable this feature."
        << std::endl;
    msg << MODULE_S << std::endl;
    msg << MODULE_S << std::endl;
    dMsg(std::cerr, msg.str(), infoMsg);
  }

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

  return 0;
}

// explicit instantiations for double and float types
template int ttk::EigenField::execute<double>(Triangulation *,
                                              double *const,
                                              const unsigned int,
                                              bool,
                                              double *const) const;
template int ttk::EigenField::execute<float>(
  Triangulation *, float *const, const unsigned int, bool, float *const) const;
