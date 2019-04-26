#pragma once

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {
  namespace Laplacian {
    /**
     * @brief Compute the Laplacian shortest path from source
     *
     * @param[in] source Source vertex from which to compute the Laplacian
     * algorithm
     * @param[in] triangulation Access to neighbor vertices, should be already
     * preprocessed
     * @param[out] outputDists Vector of distances to source for every vertex in
     * the mesh
     * @param[in] bounds Pointer to a vector of vertices that will
     * stop the algorithm once all are reached, set it to nullptr for
     * processing the whole mesh
     *
     * @return 0 in case of success
     */
    template <typename T>
    int shortestPath(const SimplexId source,
                     Triangulation &triangulation,
                     std::vector<T> &outputDists,
                     const std::vector<SimplexId> &bounds = std::vector<SimplexId>());

  } // namespace Laplacian
} // namespace ttk
