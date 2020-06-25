#pragma once

#include <Triangulation.h>

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
    template <typename T,
              class TriangulationType = AbstractTriangulation,
              typename SparseMatrixType>
    int discreteLaplacian(SparseMatrixType &output,
                          const TriangulationType &triangulation);

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
    template <typename T,
              class TriangulationType = AbstractTriangulation,
              typename SparseMatrixType>
    int cotanWeights(SparseMatrixType &output,
                     const TriangulationType &triangulation);

  } // namespace Laplacian
} // namespace ttk
