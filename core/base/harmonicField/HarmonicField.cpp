#include <HarmonicField.h>
#include <Laplacian.h>
#include <set>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>
#endif // TTK_ENABLE_EIGEN

ttk::HarmonicField::SolvingMethodType
  ttk::HarmonicField::findBestSolver(const SimplexId vertexNumber,
                                     const SimplexId edgeNumber) const {

  // for switching between Cholesky factorization and Iterate
  // (conjugate gradients) method
  const SimplexId threshold = 500000;

  // compare threshold to number of non-zero values in laplacian matrix
  if(2 * edgeNumber + vertexNumber > threshold) {
    return SolvingMethodType::ITERATIVE;
  }
  return SolvingMethodType::CHOLESKY;
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
template <class T, class TriangulationType>
int ttk::HarmonicField::execute(const TriangulationType &triangulation,
                                const SimplexId constraintNumber,
                                const SimplexId *const sources,
                                const T *const constraints,
                                T *const outputScalarField,
                                const bool useCotanWeights,
                                const SolvingMethodUserType solvingMethod,
                                const double logAlpha) const {

#ifdef TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  if(constraintNumber < 1) {
    this->printErr("Cannot solve Laplace problem with no boundary constraints");
    return 1;
  }

  using SpMat = Eigen::SparseMatrix<T>;
  using SpVec = Eigen::SparseVector<T>;
  using TripletType = Eigen::Triplet<T>;

  Timer tm;

  const auto vertexNumber = triangulation.getNumberOfVertices();
  const auto edgeNumber = triangulation.getNumberOfEdges();

  // find the right solving method
  auto findSolvingMethod = [&]() -> SolvingMethodType {
    SolvingMethodType res{};
    switch(solvingMethod) {
      case SolvingMethodUserType::AUTO:
        res = findBestSolver(vertexNumber, edgeNumber);
        break;
      case SolvingMethodUserType::CHOLESKY:
        res = SolvingMethodType::CHOLESKY;
        break;
      case SolvingMethodUserType::ITERATIVE:
        res = SolvingMethodType::ITERATIVE;
        break;
    }
    return res;
  };

  auto sm = findSolvingMethod();
  std::string begMsg{"Beginning computation... ("};
  if(useCotanWeights) {
    begMsg.append("cotan weights, ");
  } else {
    begMsg.append("discrete laplacian, ");
  }
  if(sm == SolvingMethodType::ITERATIVE) {
    begMsg.append("iterative method)");
  } else {
    begMsg.append("Cholesky method)");
  }

  this->printMsg(begMsg);

  // filter unique constraint identifiers
  std::set<SimplexId> uniqueIdentifiersSet;
  for(SimplexId i = 0; i < constraintNumber; ++i) {
    uniqueIdentifiersSet.insert(sources[i]);
  }

  // (constraint id, scalar field value)
  std::vector<std::pair<SimplexId, T>> idValues{};

  for(const auto id : uniqueIdentifiersSet) {
    for(SimplexId j = 0; j < constraintNumber; ++j) {
      if(id == sources[j]) {
        idValues.emplace_back(id, constraints[j]);
        break;
      }
    }
  }

  // unique constraint number
  size_t uniqueConstraintNumber = idValues.size();

  // graph laplacian of current mesh
  SpMat lap;
  if(useCotanWeights) {
    Laplacian::cotanWeights<T>(lap, triangulation);
  } else {
    Laplacian::discreteLaplacian<T>(lap, triangulation);
  }

  // constraints vector
  SpVec constraintsMat(vertexNumber);
  for(const auto &pair : idValues) {
    // put constraint at identifier index
    constraintsMat.coeffRef(pair.first) = pair.second;
  }

  // penalty matrix
  SpMat penalty(vertexNumber, vertexNumber);
  // penalty value
  const T alpha = Geometry::powIntTen(logAlpha);

  std::vector<TripletType> triplets;
  triplets.reserve(uniqueConstraintNumber);
  for(const auto &pair : idValues) {
    triplets.emplace_back(TripletType(pair.first, pair.first, alpha));
  }
#ifndef __clang_analyzer__
  penalty.setFromTriplets(triplets.begin(), triplets.end());
#endif // __clang_analyzer__

  int res = 0;
  SpMat sol;

  switch(sm) {
    case SolvingMethodType::CHOLESKY:
      res = solve<SpMat, SpVec, Eigen::SimplicialCholesky<SpMat>>(
        lap, penalty, constraintsMat, sol);
      break;
    case SolvingMethodType::ITERATIVE:
      res = solve<SpMat, SpVec,
                  Eigen::ConjugateGradient<SpMat, Eigen::Upper | Eigen::Lower>>(
        lap, penalty, constraintsMat, sol);
      break;
  }

  auto info = static_cast<Eigen::ComputationInfo>(res);
  switch(info) {
    case Eigen::ComputationInfo::NumericalIssue:
      this->printMsg("Numerical Issue!", ttk::debug::Priority::ERROR);
      break;
    case Eigen::ComputationInfo::NoConvergence:
      this->printMsg("No Convergence!", ttk::debug::Priority::ERROR);
      break;
    case Eigen::ComputationInfo::InvalidInput:
      this->printMsg("Invalid Input!", ttk::debug::Priority::ERROR);
      break;
    default:
      break;
  }

  // convert to dense Eigen matrix
  Eigen::Matrix<T, Eigen::Dynamic, 1> solDense(sol);

  // copy solver solution into output array
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    // cannot avoid copy here...
    outputScalarField[i] = -solDense(i, 0);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

#else
  TTK_FORCE_USE(triangulation);
  TTK_FORCE_USE(constraintNumber);
  TTK_FORCE_USE(sources);
  TTK_FORCE_USE(constraints);
  TTK_FORCE_USE(outputScalarField);
  TTK_FORCE_USE(useCotanWeights);
  TTK_FORCE_USE(solvingMethod);
  TTK_FORCE_USE(logAlpha);

  this->printMsg(
    std::vector<std::string>{
      "Eigen support disabled, computation skipped!",
      "Please re-compile TTK with Eigen support to enable this feature."},
    debug::Priority::ERROR);
#endif // TTK_ENABLE_EIGEN

  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  return 0;
}

// explicit template specializations for double and float types
#define HARMONICFIELD_SPECIALIZE(TYPE)                                       \
  template int ttk::HarmonicField::execute<TYPE>(                            \
    const Triangulation &, const SimplexId, const SimplexId *const,          \
    const TYPE *const, TYPE *const, const bool, const SolvingMethodUserType, \
    const double) const

HARMONICFIELD_SPECIALIZE(float);
HARMONICFIELD_SPECIALIZE(double);
