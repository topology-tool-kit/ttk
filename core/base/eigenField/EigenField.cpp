#include <EigenField.h>
#include <Laplacian.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>
#endif // TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Eigenvalues>
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
#endif // TTK_ENABLE_EIGEN

template <typename scalarFieldType,
          typename SparseMatrixType = Eigen::SparseMatrix<scalarFieldType>,
          typename VectorType
          = Eigen::Matrix<scalarFieldType, Eigen::Dynamic, 1>>
int ttk::EigenField::eigenfunctions(const SparseMatrixType A,
                                    VectorType &eigenVector) const {

  auto n = A.cols();
  auto m = eigenNumber_;
  const size_t minEigenNumber = 20;

  if(eigenNumber_ == 0) {
    // default value
    m = n / 1000;
  } else if(eigenNumber_ < minEigenNumber) {
    m = minEigenNumber;
  }

  // A is square
  eigen_plain_assert(n == A.rows());

  Spectra::SparseGenMatProd<scalarFieldType> op(A);
  Spectra::SymEigsSolver<scalarFieldType, Spectra::LARGEST_ALGE, decltype(op)>
    eigs(&op, m, 2 * m);

  eigs.init();
  int nconv = eigs.compute();

  int ret = eigs.info();

  // Retrieve results
  if(ret == Spectra::SUCCESSFUL) {
    eigenVector = eigs.eigenvectors().col(m - 1);
  } else if(ret == Spectra::NOT_CONVERGING) {
    std::cout << "not converging" << std::endl;
    if(nconv > 0) {
      eigenVector = eigs.eigenvectors().col(nconv - 1);
    }
  } else if(ret == Spectra::NUMERICAL_ISSUE) {
    std::cout << "numerical issue" << std::endl;
  }

  return ret;
}

// main routine
template <typename scalarFieldType>
int ttk::EigenField::execute() const {

  using std::cerr;
  using std::cout;
  using std::endl;
  using std::stringstream;

#ifdef TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  Timer t;

  {
    stringstream msg;
    msg << "[EigenField] Beginnning computation" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // graph laplacian of current mesh
  using SpMat = Eigen::SparseMatrix<scalarFieldType>;
  SpMat lap;
  Laplacian::cotanWeights<scalarFieldType>(lap, *triangulation_);

  using Vect = Eigen::Matrix<scalarFieldType, Eigen::Dynamic, 1>;
  Vect out(lap.rows());

  auto res = eigenfunctions<scalarFieldType>(lap, out);

  auto outputScalarField
    = static_cast<scalarFieldType *>(outputScalarFieldPointer_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    // cannot avoid copy here...
    outputScalarField[i] = out(i, 0);
  }

  {
    stringstream msg;
    auto info = static_cast<Eigen::ComputationInfo>(res);
    switch(info) {
      case Eigen::ComputationInfo::Success:
        msg << "[EigenField] Success!" << endl;
        break;
      case Eigen::ComputationInfo::NumericalIssue:
        msg << "[EigenField] Numerical Issue!" << endl;
        break;
      case Eigen::ComputationInfo::NoConvergence:
        msg << "[EigenField] No Convergence!" << endl;
        break;
      case Eigen::ComputationInfo::InvalidInput:
        msg << "[EigenField] Invalid Input!" << endl;
        break;
    }
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  {
    stringstream msg;
    msg << "[EigenField] Ending computation after " << t.getElapsedTime()
        << "s (" << threadNumber_ << " thread(s))" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

#else
  {
    stringstream msg;
    msg << "[EigenField]" << endl;
    msg << "[EigenField]" << endl;
    msg << "[EigenField] Eigen support disabled, computation "
           "skipped!"
        << endl;
    msg << "[EigenField] Please re-compile TTK with Eigen support to enable"
        << " this feature." << endl;
    msg << "[EigenField]" << endl;
    msg << "[EigenField]" << endl;
    dMsg(cerr, msg.str(), infoMsg);
  }
#endif // TTK_ENABLE_EIGEN

  return 0;
}

// explicit instantiations for double and float types
template int ttk::EigenField::execute<double>() const;
template int ttk::EigenField::execute<float>() const;
