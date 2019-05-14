#include <HarmonicField.h>
#include <Laplacian.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>
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
    uniqueIdentifiersSet.insert(sources_[i]);
  }

  // vector of unique constraint identifiers
  std::vector<SimplexId> uniqueIdentifiers(
    uniqueIdentifiersSet.begin(), uniqueIdentifiersSet.end());
  // vector of unique constraint values
  std::vector<scalarFieldType> uniqueValues(uniqueIdentifiers.size());

  // put identifier corresponding constraints in vector
  for(size_t i = 0; i < uniqueIdentifiers.size(); ++i) {
    for(SimplexId j = 0; j < constraintNumber_; ++j) {
      if(uniqueIdentifiers[i] == sources_[j]) {
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
