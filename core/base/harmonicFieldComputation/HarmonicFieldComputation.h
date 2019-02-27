/// \ingroup base
/// \class ttk::HarmonicFieldComputation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkHarmonicFieldComputation.cpp % for a usage example.

#pragma once

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>
#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Sparse>
#endif // TTK_ENABLE_EIGEN

namespace ttk {

class HarmonicFieldComputation : public Debug {

public:
  HarmonicFieldComputation();

  // default destructor
  ~HarmonicFieldComputation() override = default;
  // default copy constructor
  HarmonicFieldComputation(const HarmonicFieldComputation &) = default;
  // default move constructor
  HarmonicFieldComputation(HarmonicFieldComputation &&) = default;
  // default copy assignment operator
  HarmonicFieldComputation &
  operator=(const HarmonicFieldComputation &) = default;
  // default move assignment operator
  HarmonicFieldComputation &operator=(HarmonicFieldComputation &&) = default;

  inline int setVertexNumber(SimplexId vertexNumber) {
    vertexNumber_ = vertexNumber;
    return 0;
  }
  inline int setConstraintNumber(SimplexId constraintNumber) {
    constraintNumber_ = constraintNumber;
    return 0;
  }
  inline int setUseCotanWeights(bool useCotanWeights) {
    useCotanWeights_ = useCotanWeights;
    return 0;
  }
  inline int setupTriangulation(Triangulation *triangulation) {
    triangulation_ = triangulation;
    if (triangulation_ != nullptr) {
      vertexNumber_ = triangulation_->getNumberOfVertices();
      triangulation_->preprocessVertexNeighbors();
    }
    if (useCotanWeights_) {
      // cotan weights method needs more pre-processing
      triangulation_->preprocessEdgeTriangles();
    }
    return 0;
  }
  inline int setSources(void *data) {
    sources_ = data;
    return 0;
  }
  inline int setConstraints(void *data) {
    constraints_ = data;
    return 0;
  }
  inline int setOutputScalarFieldPointer(void *data) {
    outputScalarFieldPointer_ = data;
    return 0;
  }

  template <typename scalarFieldType> int execute() const;

  template <typename SparseMatrixType, typename TripletType,
            typename scalarFieldType>
  SparseMatrixType compute_laplacian() const;

  template <typename SparseMatrixType, typename TripletType,
            typename scalarFieldType>
  SparseMatrixType compute_laplacian_with_cotan_weights() const;

protected:
  // number of vertices in the mesh
  SimplexId vertexNumber_;
  // number of constraints
  SimplexId constraintNumber_;
  // cotan weights vs simple laplacian resolution
  bool useCotanWeights_;
  // the mesh
  Triangulation *triangulation_;
  // array of mesh points with scalar constraints
  // should be of constraintNumber_ size
  void *sources_;
  // array of scalar constraints on sources_
  // should be of constraintNumber_ size
  void *constraints_;
  // output of harmonic field computation
  void *outputScalarFieldPointer_;
};
} // namespace ttk

template <typename SparseMatrixType, typename TripletType,
          typename scalarFieldType>
SparseMatrixType ttk::HarmonicFieldComputation::compute_laplacian() const {

  SparseMatrixType lap(vertexNumber_, vertexNumber_);
  std::vector<TripletType> triplets;
#define USE_SYMMETRIC_LAPLACIAN
#ifdef USE_SYMMETRIC_LAPLACIAN
  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    SimplexId nneigh = triangulation_->getVertexNeighborNumber(SimplexId(i));
    triplets.emplace_back(TripletType(i, i, scalarFieldType(nneigh)));
    // rest: neighbors mapping
    for (SimplexId j = 0; j < nneigh; ++j) {
      SimplexId neighid = -1;
      triangulation_->getVertexNeighbor(i, j, neighid);
      triplets.emplace_back(TripletType(i, neighid, -1.0));
    }
  }
#else
  for (SimplexId i = 0; i < vertexNumber_; ++j) {
    triplets.emplace_back(TripletType(i, i, 1.0));
    SimplexId nneigh = triangulation_->getVertexNeighborNumber(SimplexId(i));
    // rest: neighbors mapping
    for (SimplexId j = 0; j < nneigh; ++j) {
      SimplexId neighid = -1;
      triangulation_->getVertexNeighbor(i, j, neighid);
      triplets.emplace_back(TripletType(i, j, -1.0 / scalarFieldType(nneigh)));
    }
  }
#endif // USE_SYMMETRIC_LAPLACIAN
  lap.setFromTriplets(triplets.begin(), triplets.end());
  return lap;
}

template <typename SparseMatrixType, typename TripletType,
          typename scalarFieldType>
SparseMatrixType
ttk::HarmonicFieldComputation::compute_laplacian_with_cotan_weights() const {
  using std::cout;
  using std::endl;

  {
    std::stringstream msg;
    msg << "[HarmonicFieldComputation] Beginning graph laplacian computation"
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // laplacian matrix with cotan weights
  SparseMatrixType lap(vertexNumber_, vertexNumber_);
  // used to fill more efficiently the matrix
  std::vector<TripletType> triplets;

  // iterate over all edges
  SimplexId edgesNumber = triangulation_->getNumberOfEdges();

  // hold sum of cotan weights for every vertex
  std::vector<scalarFieldType> vertexWeightSum(vertexNumber_, 0.0);

  for (SimplexId i = 0; i < edgesNumber; ++i) {

    // the two vertices of the current edge (+ a third)
    std::vector<SimplexId> edgeVertices(3);
    for (SimplexId j = 0; j < 2; ++j) {
      triangulation_->getEdgeVertex(i, j, edgeVertices[j]);
    }

    // get the triangles that share the current edge
    SimplexId trianglesNumber = triangulation_->getEdgeTriangleNumber(i);
    std::vector<SimplexId> edgeTriangles(trianglesNumber);
    for (SimplexId j = 0; j < trianglesNumber; ++j) {
      triangulation_->getEdgeTriangle(i, j, edgeTriangles[j]);
    }

    // iterate over current edge triangles
    std::vector<scalarFieldType> angles;

    for (const auto &j : edgeTriangles) {

      // get the third vertex of the triangle
      SimplexId thirdNeigh;
      for (SimplexId k = 0; k < 3; ++k) {
        triangulation_->getTriangleVertex(j, k, thirdNeigh);
        if (thirdNeigh != edgeVertices[0] && thirdNeigh != edgeVertices[1]) {
          edgeVertices[2] = thirdNeigh;
          break;
        }
      }
      // compute the 3D coords of the three vertices
      float coords[9];
      for (SimplexId k = 0; k < 3; ++k) {
        triangulation_->getVertexPoint(edgeVertices[k], coords[3 * k],
                                       coords[3 * k + 1], coords[3 * k + 2]);
      }
      angles.emplace_back(ttk::Geometry::angle(&coords[6], // edgeVertices[2]
                                               &coords[0], // edgeVertices[0]
                                               &coords[6], // edgeVertices[2]
                                               &coords[3]) // edgeVertices[1]
      );
    }

    scalarFieldType cotan_weight =
        0.5 * (std::tan(1.0 / angles[0]) + std::tan(1.0 / angles[1]));

    // since we iterate over the edges, fill the laplacian matrix
    // symmetrically for the two vertices
    triplets.emplace_back(
        TripletType(edgeVertices[0], edgeVertices[1], -cotan_weight));
    triplets.emplace_back(
        TripletType(edgeVertices[1], edgeVertices[0], -cotan_weight));

    // store the cotan weight sum for the two vertices of the current edge
    vertexWeightSum[edgeVertices[0]] += cotan_weight;
    vertexWeightSum[edgeVertices[1]] += cotan_weight;
  }

  // on the diagonal: sum of cotan weights for every vertex
  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    triplets.emplace_back(TripletType(i, i, vertexWeightSum[i]));
  }

  // put triplets into lap
  lap.setFromTriplets(triplets.begin(), triplets.end());

  {
    std::stringstream msg;
    msg << "[HarmonicFieldComputation] Graph laplacian computed" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  return lap;
}

// if the package is a pure template typename, uncomment the following line
// #include                  <HarmonicFieldComputation.cpp>

// main routine
template <typename scalarFieldType>
int ttk::HarmonicFieldComputation::execute() const {

  using std::cout;
  using std::endl;
  using std::stringstream;

  // scalar field constraints vertices
  auto identifiers = static_cast<SimplexId *>(sources_);
  // scalar field: 0 everywhere except on constraint vertices
  auto sf = static_cast<scalarFieldType *>(constraints_);

  Timer t;

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Beginnning computation" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // get unique constraint vertices
  std::set<SimplexId> identifiersSet;
  for (SimplexId i = 0; i < constraintNumber_; ++i) {
    identifiersSet.insert(identifiers[i]);
  }
  // contains vertices with constraints
  std::vector<SimplexId> identifiersVec(identifiersSet.begin(),
                                        identifiersSet.end());

#ifdef TTK_ENABLE_EIGEN
  using SpMat = Eigen::SparseMatrix<scalarFieldType>;
  using SpVec = Eigen::SparseVector<scalarFieldType>;

  // graph laplacian of current mesh
  SpMat lap;
  if (useCotanWeights_) {
    lap = compute_laplacian_with_cotan_weights<SpMat, TripletType,
                                               scalarFieldType>();
  } else {
    lap = compute_laplacian<SpMat, TripletType, scalarFieldType>();
  }

  // constraints vector
  SpVec constraints(vertexNumber_);
  for (size_t i = 0; i < identifiersVec.size(); i++) {
    // put constraint at identifier index
    constraints.coeffRef(identifiersVec[i]) = sf[i];
  }

  // penalty matrix
  SpMat penalty(vertexNumber_, vertexNumber_);
  const scalarFieldType alpha = 1.0e8;
  for (auto i : identifiersVec) {
    penalty.coeffRef(i, i) = alpha;
  }

  Eigen::SimplicialCholesky<SpMat> solver(lap - penalty);
  SpMat sol = solver.solve(penalty * constraints);

  {
    stringstream msg;
    switch (solver.info()) {
    case Eigen::ComputationInfo::Success:
      msg << "[HarmonicFieldComputation] Success!" << endl;
      break;
    case Eigen::ComputationInfo::NumericalIssue:
      msg << "[HarmonicFieldComputation] Numerical Issue!" << endl;
      break;
    case Eigen::ComputationInfo::NoConvergence:
      msg << "[HarmonicFieldComputation] No Convergence!" << endl;
      break;
    case Eigen::ComputationInfo::InvalidInput:
      msg << "[HarmonicFieldComputation] Invalid Input!" << endl;
      break;
    }
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  auto outputScalarField =
      static_cast<scalarFieldType *>(outputScalarFieldPointer_);
  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    // TODO avoid copy here
    outputScalarField[i] = std::abs(sol.coeffRef(i, 0));
  }

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

#else
  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Eigen support disabled, computation "
           "skipped "
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
#endif // TTK_ENABLE_EIGEN

  return 0;
}
