#pragma once

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {
  namespace Dijkstra {
    /**
     * @brief Compute the Dijkstra shortest path from source
     *
     * @param[in] source Source vertex for the Dijkstra algorithm
     * @param[in] triangulation Access to neighbor vertices
     * @param[out] outputDists Distances to source for every mesh vertex
     * @param[in] bounds Stop the algorithim if all vertices are reached
     * @param[in] mask Vector masking the triangulation
     *
     * @return 0 in case of success
     */
    template <typename T>
    int shortestPath(const SimplexId source,
                     Triangulation &triangulation,
                     std::vector<T> &outputDists,
                     const std::vector<SimplexId> &bounds
                     = std::vector<SimplexId>(),
                     const std::vector<bool> &mask = std::vector<bool>());

  } // namespace Dijkstra
} // namespace ttk
