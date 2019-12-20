#include <Geometry.h>
#include <HarmonicField.h>
#include <Laplacian.h>
#include <cmath>
#include <set>
#include <type_traits>

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
template <typename ScalarFieldType>
int ttk::HarmonicField::execute(Triangulation *triangulation,
                                SimplexId constraintNumber,
                                SimplexId *sources,
                                ScalarFieldType *constraints,
                                ScalarFieldType *outputScalarField,
                                bool useCotanWeights,
                                SolvingMethodUserType solvingMethod,
                                double logAlpha) const {

#ifdef TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_OPENMP
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP

  using SpMat = Eigen::SparseMatrix<ScalarFieldType>;
  using SpVec = Eigen::SparseVector<ScalarFieldType>;
  using TripletType = Eigen::Triplet<ScalarFieldType>;

  Timer tm;
  Memory mem;

  const auto vertexNumber = triangulation->getNumberOfVertices();
  const auto edgeNumber = triangulation->getNumberOfEdges();

  // filter unique constraint identifiers
  std::set<SimplexId> uniqueIdentifiersSet;
  for(SimplexId i = 0; i < constraintNumber; ++i) {
    uniqueIdentifiersSet.insert(sources[i]);
  }

  // vector of unique constraint identifiers
  std::vector<SimplexId> uniqueIdentifiers(
    uniqueIdentifiersSet.begin(), uniqueIdentifiersSet.end());
  // vector of unique constraint values
  std::vector<ScalarFieldType> uniqueValues(uniqueIdentifiers.size());

  // put identifier corresponding constraints in vector
  for(size_t i = 0; i < uniqueIdentifiers.size(); ++i) {
    for(SimplexId j = 0; j < constraintNumber; ++j) {
      if(uniqueIdentifiers[i] == sources[j]) {
        uniqueValues[i] = constraints[j];
        break;
      }
    }
  }

  // unique constraint number
  size_t uniqueConstraintNumber = uniqueValues.size();

  // graph laplacian of current mesh
  SpMat lap;
  if(useCotanWeights) {
    Laplacian::cotanWeights<ScalarFieldType>(lap, *triangulation);
  } else {
    Laplacian::discreteLaplacian<ScalarFieldType>(lap, *triangulation);
  }

  // constraints vector
  SpVec constraintsMat(vertexNumber);
  for(size_t i = 0; i < uniqueConstraintNumber; ++i) {
    // put constraint at identifier index
    constraintsMat.coeffRef(uniqueIdentifiers[i]) = uniqueValues[i];
  }

  auto sm = SolvingMethodType::CHOLESKY;

  switch(solvingMethod) {
    case SolvingMethodUserType::AUTO:
      sm = findBestSolver(vertexNumber, edgeNumber);
      break;
    case SolvingMethodUserType::CHOLESKY:
      sm = SolvingMethodType::CHOLESKY;
      break;
    case SolvingMethodUserType::ITERATIVE:
      sm = SolvingMethodType::ITERATIVE;
      break;
  }

  // penalty matrix
  SpMat penalty(vertexNumber, vertexNumber);
  // penalty value
  const ScalarFieldType alpha = pow10(logAlpha);

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
      this->PrintMsg("Numerical Issue!", ttk::debug::Priority::ERROR);
      break;
    case Eigen::ComputationInfo::NoConvergence:
      this->PrintMsg("No Convergence!", ttk::debug::Priority::ERROR);
      break;
    case Eigen::ComputationInfo::InvalidInput:
      this->PrintMsg("Invalid Input!", ttk::debug::Priority::ERROR);
      break;
    default:
      break;
  }

  // convert to dense Eigen matrix
  Eigen::Matrix<ScalarFieldType, Eigen::Dynamic, 1> solDense(sol);

  // copy solver solution into output array
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    // cannot avoid copy here...
    outputScalarField[i] = -solDense(i, 0);
  }

  std::string endMsg{"Complete ("};
  if(useCotanWeights) {
    endMsg.append("cotan weights, ");
  } else {
    endMsg.append("discrete laplacian, ");
  }
  if(sm == SolvingMethodType::ITERATIVE) {
    endMsg.append("iterative)");
  } else {
    endMsg.append("Cholesky)");
  }

  this->PrintMsg(endMsg, 1.0, tm.getElapsedTime(), mem.getElapsedUsage());

#else
  this->PrintMsg(
    std::vector<std::string>{
      "Eigen support disabled, computation skipped!",
      "Please re-compile TTK with Eigen support to enable this feature."},
    ttk::debug::Priority::ERROR);
#endif // TTK_ENABLE_EIGEN

  return 0;
}

// explicit template specializations for double and float types
#define HARMONICFIELD_SPECIALIZE(TYPE)                                 \
  template int ttk::HarmonicField::execute<TYPE>(                      \
    Triangulation * triangulation, SimplexId constraintNumber,         \
    SimplexId * sources, TYPE * constraints, TYPE * outputScalarField, \
    bool useCotanWeights = true,                                       \
    SolvingMethodUserType solvingMethod = SolvingMethodUserType::AUTO, \
    double logAlpha = 5.0) const

HARMONICFIELD_SPECIALIZE(float);
HARMONICFIELD_SPECIALIZE(double);
