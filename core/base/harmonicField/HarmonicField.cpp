#include <HarmonicField.h>
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

ttk::SolvingMethodType ttk::HarmonicField::findBestSolver() const {

  // for switching between Cholesky factorization and Iterate
  // (conjugate gradients) method
  const SimplexId threshold = 500000;

  // compare threshold to number of non-zero values in laplacian matrix
  if(2 * edgeNumber_ + vertexNumber_ > threshold) {
    return ttk::SolvingMethodType::Iterative;
  }
  return ttk::SolvingMethodType::Cholesky;
}

template <typename SparseMatrixType,
          typename SparseVectorType,
          typename SolverType>
int ttk::HarmonicField::solve(SparseMatrixType const &lap,
                              SparseMatrixType const &penalty,
                              SparseVectorType const &constraints,
                              SparseMatrixType &sol) const {
  SolverType solver(lap - penalty);
  sol = solver.solve(penalty * constraints);
  return solver.info();
}

template <typename scalarFieldType,
          typename SparseMatrixType = Eigen::SparseMatrix<scalarFieldType>,
          typename VectorType
          = Eigen::Matrix<scalarFieldType, Eigen::Dynamic, 1>>
int ttk::HarmonicField::eigenfunctions(const SparseMatrixType A,
                                       VectorType &eigenVector,
                                       const size_t eigenNumber) const {

  auto n = A.cols();
  auto m = eigenNumber;

  if(eigenNumber == 0) {
    // default value
    m = n / 1000;
  } else if(eigenNumber < 20) {
    m = 20;
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
int ttk::HarmonicField::execute() const {

  using std::cerr;
  using std::cout;
  using std::endl;
  using std::stringstream;

#ifdef TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  using SpMat = Eigen::SparseMatrix<scalarFieldType>;
  using SpVec = Eigen::SparseVector<scalarFieldType>;
  using TripletType = Eigen::Triplet<scalarFieldType>;

  Timer t;

  // scalar field constraints vertices
  auto identifiers = static_cast<SimplexId *>(sources_);
  // scalar field: 0 everywhere except on constraint vertices
  auto sf = static_cast<scalarFieldType *>(constraints_);

  {
    stringstream msg;
    msg << "[HarmonicField] Beginnning computation" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // filter unique constraint identifiers
  std::set<SimplexId> uniqueIdentifiersSet;
  for(SimplexId i = 0; i < constraintNumber_; ++i) {
    uniqueIdentifiersSet.insert(identifiers[i]);
  }

  // vector of unique constraint identifiers
  std::vector<SimplexId> uniqueIdentifiers(
    uniqueIdentifiersSet.begin(), uniqueIdentifiersSet.end());
  // vector of unique constraint values
  std::vector<scalarFieldType> uniqueValues(uniqueIdentifiers.size());

  // put identifier corresponding constraints in vector
  for(size_t i = 0; i < uniqueIdentifiers.size(); ++i) {
    for(SimplexId j = 0; j < constraintNumber_; ++j) {
      if(uniqueIdentifiers[i] == identifiers[j]) {
        uniqueValues[i] = sf[j];
        break;
      }
    }
  }

  // unique constraint number
  size_t uniqueConstraintNumber = uniqueValues.size();

  // graph laplacian of current mesh
  SpMat lap;
  if(useCotanWeights_) {
    Laplacian::cotanWeights<scalarFieldType>(lap, *triangulation_);
  } else {
    Laplacian::discreteLaplacian<scalarFieldType>(lap, *triangulation_);
  }

  {
    Timer te;

    using Vect = Eigen::Matrix<scalarFieldType, Eigen::Dynamic, 1>;
    Vect out(lap.rows());

    eigenfunctions<scalarFieldType>(lap, out, size_t(logAlpha_));

    auto outputScalarField
      = static_cast<scalarFieldType *>(outputScalarFieldPointer_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      // cannot avoid copy here...
      outputScalarField[i] = out(i, 0);
    }

    stringstream msg;
    msg << "[HarmonicField] Laplacian eigenvalues " << te.getElapsedTime()
        << "s (discrete laplacian, " << threadNumber_ << " thread(s))" << endl;
    dMsg(cout, msg.str(), infoMsg);

    return 0;
  }

  // constraints vector
  SpVec constraints(vertexNumber_);
  for(size_t i = 0; i < uniqueConstraintNumber; ++i) {
    // put constraint at identifier index
    constraints.coeffRef(uniqueIdentifiers[i]) = uniqueValues[i];
  }

  auto sm = ttk::SolvingMethodType::Cholesky;

  switch(solvingMethod_) {
    case ttk::SolvingMethodUserType::Auto:
      sm = findBestSolver();
      break;
    case ttk::SolvingMethodUserType::Cholesky:
      sm = ttk::SolvingMethodType::Cholesky;
      break;
    case ttk::SolvingMethodUserType::Iterative:
      sm = ttk::SolvingMethodType::Iterative;
      break;
  }

  // penalty matrix
  SpMat penalty(vertexNumber_, vertexNumber_);
  // penalty value
  const scalarFieldType alpha = pow10(logAlpha_);

  std::vector<TripletType> triplets;
  triplets.reserve(uniqueConstraintNumber);
  for(size_t i = 0; i < uniqueConstraintNumber; ++i) {
    triplets.emplace_back(
      TripletType(uniqueIdentifiers[i], uniqueIdentifiers[i], alpha));
  }
  penalty.setFromTriplets(triplets.begin(), triplets.end());

  int res = 0;
  SpMat sol;

  switch(sm) {
    case ttk::SolvingMethodType::Cholesky:
      res = solve<SpMat, SpVec, Eigen::SimplicialCholesky<SpMat>>(
        lap, penalty, constraints, sol);
      break;
    case ttk::SolvingMethodType::Iterative:
      res = solve<SpMat, SpVec,
                  Eigen::ConjugateGradient<SpMat, Eigen::Upper | Eigen::Lower>>(
        lap, penalty, constraints, sol);
      break;
  }

  {
    stringstream msg;
    auto info = static_cast<Eigen::ComputationInfo>(res);
    switch(info) {
      case Eigen::ComputationInfo::Success:
        msg << "[HarmonicField] Success!" << endl;
        break;
      case Eigen::ComputationInfo::NumericalIssue:
        msg << "[HarmonicField] Numerical Issue!" << endl;
        break;
      case Eigen::ComputationInfo::NoConvergence:
        msg << "[HarmonicField] No Convergence!" << endl;
        break;
      case Eigen::ComputationInfo::InvalidInput:
        msg << "[HarmonicField] Invalid Input!" << endl;
        break;
    }
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  auto outputScalarField
    = static_cast<scalarFieldType *>(outputScalarFieldPointer_);

  // convert to dense Eigen matrix
  Eigen::Matrix<scalarFieldType, Eigen::Dynamic, 1> solDense(sol);

  // copy solver solution into output array
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    // cannot avoid copy here...
    outputScalarField[i] = -solDense(i, 0);
  }

  {
    stringstream msg;
    msg << "[HarmonicField] Ending computation after " << t.getElapsedTime()
        << "s (";
    if(useCotanWeights_) {
      msg << "cotan weights, ";
    } else {
      msg << "discrete laplacian, ";
    }
    if(sm == ttk::SolvingMethodType::Iterative) {
      msg << "iterative solver, ";
    } else {
      msg << "Cholesky, ";
    }
    msg << threadNumber_ << " thread(s))" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

#else
  {
    stringstream msg;
    msg << "[HarmonicField]" << endl;
    msg << "[HarmonicField]" << endl;
    msg << "[HarmonicField] Eigen support disabled, computation "
           "skipped!"
        << endl;
    msg << "[HarmonicField] Please re-compile TTK with Eigen support to enable"
        << " this feature." << endl;
    msg << "[HarmonicField]" << endl;
    msg << "[HarmonicField]" << endl;
    dMsg(cerr, msg.str(), infoMsg);
  }
#endif // TTK_ENABLE_EIGEN

  return 0;
}

// explicit instantiations for double and float types
template int ttk::HarmonicField::execute<double>() const;
template int ttk::HarmonicField::execute<float>() const;
