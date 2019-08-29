/// \ingroup base
/// \class ttk::HarmonicField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2019
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
///
///
/// \sa ttkHarmonicField.cpp % for a usage example.

#pragma once

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>
#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

namespace ttk {

  enum struct SolvingMethodUserType { Auto = 0, Cholesky = 1, Iterative = 2 };
  enum struct SolvingMethodType { Cholesky, Iterative };

  class HarmonicField : public Debug {

  public:
    // default constructor
    HarmonicField() = default;
    // default destructor
    ~HarmonicField() override = default;
    // default copy constructor
    HarmonicField(const HarmonicField &) = default;
    // default move constructor
    HarmonicField(HarmonicField &&) = default;
    // default copy assignment operator
    HarmonicField &operator=(const HarmonicField &) = default;
    // default move assignment operator
    HarmonicField &operator=(HarmonicField &&) = default;

    inline int setVertexNumber(const SimplexId vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }
    inline int setConstraintNumber(const SimplexId constraintNumber) {
      constraintNumber_ = constraintNumber;
      return 0;
    }
    inline int setUseCotanWeights(const bool useCotanWeights) {
      useCotanWeights_ = useCotanWeights;
      return 0;
    }
    inline int setupTriangulation(Triangulation *const triangulation) {
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
    inline int setSources(void *const data) {
      sources_ = static_cast<SimplexId *>(data);
      return 0;
    }
    inline int setConstraints(void *const data) {
      constraints_ = data;
      return 0;
    }
    inline int setOutputScalarFieldPointer(void *const data) {
      outputScalarFieldPointer_ = data;
      return 0;
    }
    inline int setSolvingMethod(const int solvingMethod) {
      solvingMethod_ = static_cast<SolvingMethodUserType>(solvingMethod);
      return 0;
    }
    inline int setLogAlpha(const double logAlpha) {
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

  private:
    // number of vertices in the mesh
    SimplexId vertexNumber_{};
    // number of edges in the mesh
    SimplexId edgeNumber_{};
    // number of constraints
    SimplexId constraintNumber_{};
    // cotan weights vs simple laplacian resolution
    bool useCotanWeights_{true};
    // the mesh
    Triangulation *triangulation_{};
    // array of mesh points with scalar constraints
    // should be of constraintNumber_ size
    SimplexId *sources_{};
    // array of scalar constraints on sources_
    // should be of constraintNumber_ size
    void *constraints_{};
    // output of harmonic field computation
    void *outputScalarFieldPointer_{};
    // user-selected solver
    SolvingMethodUserType solvingMethod_{ttk::SolvingMethodUserType::Auto};
    // log10 of penalty value
    double logAlpha_{5};
  };
} // namespace ttk
