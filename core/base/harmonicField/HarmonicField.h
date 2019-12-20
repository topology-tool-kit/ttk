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
#include <Triangulation.h>

namespace ttk {

  class HarmonicField : virtual public Debug {

  protected:
    enum class SolvingMethodUserType { AUTO, CHOLESKY, ITERATIVE };
    enum class SolvingMethodType { CHOLESKY, ITERATIVE };

    inline void setupTriangulation(Triangulation *const triangulation,
                                   bool cotanWeights) {
      if(triangulation == nullptr) {
        return;
      }
      triangulation->preconditionVertexNeighbors();
      triangulation->preconditionEdges();
      if(cotanWeights) {
        // cotan weights method needs more pre-processing
        triangulation->preconditionEdgeTriangles();
      }
    }

    SolvingMethodType findBestSolver(const SimplexId vertexNumber,
                                     const SimplexId edgeNumber) const;

    template <typename SparseMatrixType,
              typename SparseVectorType,
              typename SolverType>
    int solve(SparseMatrixType const &lap,
              SparseMatrixType const &penalty,
              SparseVectorType const &constraints,
              SparseMatrixType &sol) const;

    template <typename ScalarFieldType>
    int execute(Triangulation *triangulation,
                SimplexId constraintNumber,
                SimplexId *sources,
                ScalarFieldType *constraints,
                ScalarFieldType *outputScalarField,
                bool useCotanWeights = true,
                SolvingMethodUserType solvingMethod
                = SolvingMethodUserType::AUTO,
                double logAlpha = 5.0) const;
  };
} // namespace ttk
