#pragma once

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {
  namespace Laplacian {
    /**
     * @brief Compute the Laplacian matrix of the graph
     *
     * @param[out] output Laplacian matrix
     * @param[in] triangulation Access to neighbor vertices, should be already
     * preprocessed
     *
     * @return 0 in case of success
     */
    template <typename T, typename SparseMatrixType>
    int discreteLaplacian(SparseMatrixType &output,
                          const Triangulation &triangulation);

    /**
     * @brief Compute the Laplacian matrix of the graph using the
     * cotangente weights method
     *
     * @param[out] output Laplacian matrix
     * @param[in] triangulation Access to neighbor vertices, should be already
     * preprocessed
     *
     * @return 0 in case of success
     */
    template <typename T, typename SparseMatrixType>
    int cotanWeights(SparseMatrixType &output,
                     const Triangulation &triangulation);

  } // namespace Laplacian
} // namespace ttk
