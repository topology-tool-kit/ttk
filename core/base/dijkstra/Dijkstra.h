#pragma once

#include <Geometry.h>
#include <Triangulation.h>

#include <array>
#include <functional>
#include <limits>
#include <queue>

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
    template <typename T,
              typename triangulationType = ttk::AbstractTriangulation>
    int shortestPath(const SimplexId source,
                     const triangulationType &triangulation,
                     std::vector<T> &outputDists,
                     const std::vector<SimplexId> &bounds
                     = std::vector<SimplexId>(),
                     const std::vector<bool> &mask = std::vector<bool>()) {

      // should we process the whole mesh or stop at some point?
      bool processAllVertices = bounds.empty();
      // total number of vertices in the mesh
      size_t vertexNumber = triangulation.getNumberOfVertices();
      // is there a mask?
      bool isMask = !mask.empty();

      // check mask size
      if(isMask && mask.size() != vertexNumber) {
        return 1;
      }

      // list all reached bounds
      std::vector<bool> reachedBounds;

      // alloc and fill reachedBounds
      if(!processAllVertices) {
        reachedBounds.resize(bounds.size(), false);
      }

      // preprocess output vector
      outputDists.clear();
      outputDists.resize(vertexNumber, std::numeric_limits<T>::infinity());

      // link vertex and current distance to source
      using pq_t = std::pair<T, SimplexId>;

      // priority queue storing pairs of (distance, vertices TTK id)
      std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> pq;

      // init pipeline
      pq.push(std::make_pair(T(0.0F), source));
      outputDists[source] = T(0.0F);

      while(!pq.empty()) {
        auto elem = pq.top();
        pq.pop();
        auto vert = elem.second;
        std::array<float, 3> vCoords{};
        triangulation.getVertexPoint(vert, vCoords[0], vCoords[1], vCoords[2]);

        auto nneigh = triangulation.getVertexNeighborNumber(vert);

        for(SimplexId i = 0; i < nneigh; i++) {
          // neighbor Id
          SimplexId neigh{};
          triangulation.getVertexNeighbor(vert, i, neigh);

          // limit to masked vertices
          if(isMask && !mask[neigh]) {
            continue;
          }

          // neighbor coordinates
          std::array<float, 3> nCoords{};
          triangulation.getVertexPoint(
            neigh, nCoords[0], nCoords[1], nCoords[2]);
          // (square) distance between vertex and neighbor
          T distVN = Geometry::distance(vCoords.data(), nCoords.data());
          if(outputDists[neigh] > outputDists[vert] + distVN) {
            outputDists[neigh] = outputDists[vert] + distVN;
            if(!processAllVertices) {
              // check if neigh in bounds
              auto it = std::find(bounds.begin(), bounds.end(), neigh);
              if(it != bounds.end()) {
                // mark it as found
                reachedBounds[it - bounds.begin()] = true;
              }
              // break if all are found
              if(std::all_of(reachedBounds.begin(), reachedBounds.end(),
                             [](const bool v) { return v; })) {
                break;
              }
            }
            pq.push(std::make_pair(outputDists[neigh], neigh));
          }
        }
      }

      return 0;
    }

  } // namespace Dijkstra
} // namespace ttk
