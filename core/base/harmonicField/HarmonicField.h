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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n

#pragma once

// base code includes
#include <Laplacian.h>
#include <Triangulation.h>

namespace ttk {

  class HarmonicField : virtual public Debug {

  protected:
    enum class SolvingMethodUserType { AUTO, CHOLESKY, ITERATIVE };
    enum class SolvingMethodType { CHOLESKY, ITERATIVE };

    HarmonicField() {
      this->setDebugMsgPrefix("HarmonicField");
    }

    inline void
      preconditionTriangulation(AbstractTriangulation &triangulation) const {
      Laplacian::preconditionTriangulation(triangulation);
    }

    template <class T, class TriangulationType = AbstractTriangulation>
    int execute(const TriangulationType &triangulation,
                const SimplexId constraintNumber,
                const SimplexId *const sources,
                const T *const constraints,
                T *const outputScalarField,
                const bool useCotanWeights = true,
                const SolvingMethodUserType solvingMethod
                = SolvingMethodUserType::AUTO,
                const double logAlpha = 5.0) const;

  private:
    SolvingMethodType findBestSolver(const SimplexId vertexNumber,
                                     const SimplexId edgeNumber) const;

    template <typename SparseMatrixType,
              typename SparseVectorType,
              typename SolverType>
    int solve(SparseMatrixType const &lap,
              SparseMatrixType const &penalty,
              SparseVectorType const &constraints,
              SparseMatrixType &sol) const;
  };
} // namespace ttk
