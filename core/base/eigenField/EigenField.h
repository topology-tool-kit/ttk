/// \ingroup base
/// \class ttk::EigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date April 2019
///
/// \brief TTK processing package for computing eigenfunctions of a
/// triangular mesh.
///
/// \sa ttkEigenField.cpp % for a usage example.

#pragma once

// base code includes
#include <Debug.h>
#include <Geometry.h>
#include <Laplacian.h>
#include <Triangulation.h>

#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_SPECTRA)
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#endif // TTK_ENABLE_EIGEN && TTK_ENABLE_SPECTRA

namespace ttk {

  class EigenField : virtual public Debug {
  public:
    EigenField() {
      this->setDebugMsgPrefix("EigenField");
    }

    inline void preconditionTriangulation(Triangulation *triangulation) const {
      if(triangulation != nullptr) {
        triangulation->preconditionVertexNeighbors();
        // cotan weights method needs more pre-processing
        triangulation->preconditionEdgeTriangles();
      }
    }

    template <typename T>
    int execute(Triangulation *triangulation,
                T *const outputFieldPointer,
                const unsigned int eigenNumber = 500,
                bool computeStatistics = false,
                T *const outputStatistics = nullptr) const;
  };

} // namespace ttk
