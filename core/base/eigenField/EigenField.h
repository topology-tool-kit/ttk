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
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  class EigenField : public Debug {
  public:
    inline void setupTriangulation(Triangulation *triangulation) {
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
