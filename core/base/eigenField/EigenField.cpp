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
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif // __GNUC__

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

// main routine
template <typename T>
int ttk::EigenField::execute() const {

  Timer t;

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  using SpMat = Eigen::SparseMatrix<T>;
  using DMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  {
    std::stringstream msg;
    msg << MODULE_S "Beginnning computation" << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  // graph laplacian of current mesh
  SpMat lap;
  // compute graph laplacian using cotangent weights
  Laplacian::cotanWeights<T>(lap, *triangulation_);
  // lap is square
  eigen_plain_assert(lap.cols() == lap.rows());

  auto n = lap.cols();
  auto m = eigenNumber_;
  // threshold: minimal number of eigenpairs to get a converging solution
  const size_t minEigenNumber = 20;

  if(eigenNumber_ == 0) {
    // default value
    m = n / 1000;
  } else if(eigenNumber_ < minEigenNumber) {
    m = minEigenNumber;
  }

  Spectra::SparseGenMatProd<T> op(lap);
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
            << eigenNumber_ << " values computed)" << std::endl;
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

  auto outputScalarField = static_cast<T *>(outputScalarFieldPointer_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    // cannot avoid copy here...
    outputScalarField[i] = eigenvectors(i, 0);
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
template int ttk::EigenField::execute<double>() const;
template int ttk::EigenField::execute<float>() const;
