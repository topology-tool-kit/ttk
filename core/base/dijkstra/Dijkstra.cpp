#include <Dijkstra.h>
#include <Geometry.h>

#include <array>
#include <functional>
#include <queue>

template <typename T, typename triangulationType>
int ttk::Dijkstra::shortestPath(const ttk::SimplexId source,
                                const triangulationType &triangulation,
                                std::vector<T> &outputDists,
                                const std::vector<ttk::SimplexId> &bounds,
                                const std::vector<bool> &mask) {

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
      triangulation.getVertexPoint(neigh, nCoords[0], nCoords[1], nCoords[2]);
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

// explicit template specializations for floating-point types
#define DIJKSTRA_SPECIALIZE(DATATYPE, TRIANGLTYPE)                 \
  template int ttk::Dijkstra::shortestPath<DATATYPE, TRIANGLTYPE>( \
    const ttk::SimplexId source, const TRIANGLTYPE &triangulation, \
    std::vector<DATATYPE> &outputDists,                            \
    const std::vector<ttk::SimplexId> &bounds, const std::vector<bool> &mask);

DIJKSTRA_SPECIALIZE(float, ttk::ExplicitTriangulation)
DIJKSTRA_SPECIALIZE(float, ttk::ImplicitTriangulation)
DIJKSTRA_SPECIALIZE(float, ttk::PeriodicImplicitTriangulation)
DIJKSTRA_SPECIALIZE(double, ttk::ExplicitTriangulation)
DIJKSTRA_SPECIALIZE(double, ttk::ImplicitTriangulation)
DIJKSTRA_SPECIALIZE(double, ttk::PeriodicImplicitTriangulation)
