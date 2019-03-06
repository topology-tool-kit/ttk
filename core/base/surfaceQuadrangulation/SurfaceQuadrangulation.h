/// \ingroup base
/// \class ttk::SurfaceQuadrangulation
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkSurfaceQuadrangulation.cpp % for a usage example.

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

  enum struct SolvingMethodUserType { Auto = 0, Cholesky = 1, Iterative = 2 };
  enum struct SolvingMethodType { Cholesky, Iterative };

  class SurfaceQuadrangulation : public Debug {

  public:
    SurfaceQuadrangulation();

    // default destructor
    ~SurfaceQuadrangulation() override = default;
    // default copy constructor
    SurfaceQuadrangulation(const SurfaceQuadrangulation &) = default;
    // default move constructor
    SurfaceQuadrangulation(SurfaceQuadrangulation &&) = default;
    // default copy assignment operator
    SurfaceQuadrangulation &operator=(const SurfaceQuadrangulation &) = default;
    // default move assignment operator
    SurfaceQuadrangulation &operator=(SurfaceQuadrangulation &&) = default;

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
      if(triangulation_ != nullptr) {
        vertexNumber_ = triangulation_->getNumberOfVertices();
        triangulation_->preprocessVertexNeighbors();
        edgeNumber_ = triangulation_->getNumberOfEdges();
      }
      if(useCotanWeights_) {
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
    inline int setSolvingMethod(int solvingMethod) {
      solvingMethod_ = static_cast<SolvingMethodUserType>(solvingMethod);
      return 0;
    }
    inline int setLogAlpha(double logAlpha) {
      logAlpha_ = logAlpha;
      return 0;
    }

    SolvingMethodType findBestSolver() const;

    template <typename SparseMatrixType,
              typename SparseVectorType,
              typename SolverType>
    int solve(SparseMatrixType const &lap,
              SparseMatrixType const &penalty,
              SparseVectorType const &constraints,
              SparseMatrixType &sol) const;

    template <typename scalarFieldType>
    int execute() const;

    template <typename SparseMatrixType,
              typename TripletType,
              typename scalarFieldType>
    SparseMatrixType compute_laplacian() const;

    template <typename SparseMatrixType,
              typename TripletType,
              typename scalarFieldType>
    SparseMatrixType compute_laplacian_with_cotan_weights() const;

  protected:
    // number of vertices in the mesh
    SimplexId vertexNumber_;
    // number of edges in the mesh
    SimplexId edgeNumber_;
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
    // user-selected solver
    SolvingMethodUserType solvingMethod_;
    // log10 of penalty value
    double logAlpha_;
  };
} // namespace ttk

template <typename SparseMatrixType,
          typename TripletType,
          typename scalarFieldType>
SparseMatrixType ttk::SurfaceQuadrangulation::compute_laplacian() const {

  // laplacian matrix
  SparseMatrixType lap(vertexNumber_, vertexNumber_);

  // number of triplets to insert into laplacian matrix: vertexNumber_
  // values on the diagonal + 2 values per edge
  std::vector<TripletType> triplets(vertexNumber_ + 2 * edgeNumber_);

#define USE_SYMMETRIC_LAPLACIAN

#ifdef USE_SYMMETRIC_LAPLACIAN

  // on the diagonal: number of neighbors
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    SimplexId nneigh = triangulation_->getVertexNeighborNumber(SimplexId(i));
    triplets[i] = TripletType(i, i, scalarFieldType(nneigh));
  }

  // neighbors mapping: loop over edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber_; ++i) {
    // the two vertices of the current edge
    std::vector<SimplexId> edgeVertices(2);
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation_->getEdgeVertex(i, j, edgeVertices[j]);
    }
    // fill triplets for both vertices of current edge
    triplets[vertexNumber_ + 2 * i]
      = TripletType(edgeVertices[0], edgeVertices[1], -1.0);
    triplets[vertexNumber_ + 2 * i + 1]
      = TripletType(edgeVertices[1], edgeVertices[0], -1.0);
  }

#else // USE_SYMMETRIC_LAPLACIAN

  // on the diagonal: 1.0
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    triplets[i] = TripletType(i, i, 1.0);
  }

  // neighbors mapping: loop over edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber; ++i) {
    // the two vertices of the current edge
    std::vector<SimplexId> edgeVertices(2);
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation_->getEdgeVertex(i, j, edgeVertices[j]);
    }
    SimplexId nneigh0
      = triangulation_->getVertexNeighborNumber(edgeVertices[0]);
    SimplexId nneigh1
      = triangulation_->getVertexNeighborNumber(edgeVertices[1]);

    // fill triplets for both vertices of current edge
    triplets[vertexNumber_ + 2 * i]
      = TripletType(edgeVertices[0], edgeVertices[1], -1.0 / nneigh0);
    triplets[vertexNumber_ + 2 * i + 1]
      = TripletType(edgeVertices[1], edgeVertices[0], -1.0 / nneigh1);
  }

#endif // USE_SYMMETRIC_LAPLACIAN

  lap.setFromTriplets(triplets.begin(), triplets.end());
  return lap;
}

template <typename SparseMatrixType,
          typename TripletType,
          typename scalarFieldType>
SparseMatrixType
  ttk::SurfaceQuadrangulation::compute_laplacian_with_cotan_weights() const {
  using std::cout;
  using std::endl;

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Beginning graph laplacian computation"
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // laplacian matrix with cotan weights
  SparseMatrixType lap(vertexNumber_, vertexNumber_);

  // used to fill more efficiently the matrix
  std::vector<TripletType> triplets(vertexNumber_ + 2 * edgeNumber_);

  // hold sum of cotan weights for every vertex
  std::vector<scalarFieldType> vertexWeightSum(vertexNumber_, 0.0);

  // iterate over all edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < edgeNumber_; ++i) {

    // the two vertices of the current edge (+ a third)
    std::vector<SimplexId> edgeVertices(3);
    for(SimplexId j = 0; j < 2; ++j) {
      triangulation_->getEdgeVertex(i, j, edgeVertices[j]);
    }

    // get the triangles that share the current edge
    // in 2D only 2, in 3D, maybe more...
    SimplexId trianglesNumber = triangulation_->getEdgeTriangleNumber(i);
    // stores the triangles ID for every triangle around the current edge
    std::vector<SimplexId> edgeTriangles(trianglesNumber);
    for(SimplexId j = 0; j < trianglesNumber; ++j) {
      triangulation_->getEdgeTriangle(i, j, edgeTriangles[j]);
    }

    // iterate over current edge triangles
    std::vector<scalarFieldType> angles;

    for(const auto &j : edgeTriangles) {

      // get the third vertex of the triangle
      SimplexId thirdNeigh;
      // a triangle has only three vertices
      for(SimplexId k = 0; k < 3; ++k) {
        triangulation_->getTriangleVertex(j, k, thirdNeigh);
        if(thirdNeigh != edgeVertices[0] && thirdNeigh != edgeVertices[1]) {
          // store the third vertex ID into the edgeVertices array to
          // be more easily handled
          edgeVertices[2] = thirdNeigh;
          break;
        }
      }
      // compute the 3D coords of the three vertices
      float coords[9];
      for(SimplexId k = 0; k < 3; ++k) {
        triangulation_->getVertexPoint(
          edgeVertices[k], coords[3 * k], coords[3 * k + 1], coords[3 * k + 2]);
      }
      angles.emplace_back(ttk::Geometry::angle(&coords[6], // edgeVertices[2]
                                               &coords[0], // edgeVertices[0]
                                               &coords[6], // edgeVertices[2]
                                               &coords[3]) // edgeVertices[1]
      );
    }

    // cotan weights for every triangle around the current edge
    scalarFieldType cotan_weight = 0.0;
    // C++ has no map statement until C++17 (std::transform)
    for(auto &angle : angles) {
      cotan_weight += 1.0 / std::tan(angle);
    }

    // since we iterate over the edges, fill the laplacian matrix
    // symmetrically for the two vertices
    triplets[2 * i]
      = TripletType(edgeVertices[0], edgeVertices[1], -cotan_weight);
    triplets[2 * i + 1]
      = TripletType(edgeVertices[1], edgeVertices[0], -cotan_weight);

    // store the cotan weight sum for the two vertices of the current edge
    vertexWeightSum[edgeVertices[0]] += cotan_weight;
    vertexWeightSum[edgeVertices[1]] += cotan_weight;
  }

  // on the diagonal: sum of cotan weights for every vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber_; ++i) {
    triplets[2 * edgeNumber_ + i] = TripletType(i, i, vertexWeightSum[i]);
  }

  // put triplets into lap
  lap.setFromTriplets(triplets.begin(), triplets.end());

  {
    std::stringstream msg;
    msg << "[SurfaceQuadrangulation] Graph laplacian computed" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  return lap;
}

// if the package is a pure template typename, uncomment the following line
// #include                  <SurfaceQuadrangulation.cpp>

template <typename SparseMatrixType,
          typename SparseVectorType,
          typename SolverType>
int ttk::SurfaceQuadrangulation::solve(SparseMatrixType const &lap,
                                       SparseMatrixType const &penalty,
                                       SparseVectorType const &constraints,
                                       SparseMatrixType &sol) const {
  SolverType solver(lap - penalty);
  sol = solver.solve(penalty * constraints);
  return solver.info();
}
