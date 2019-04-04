#include <Dijkstra.h>
#include <queue>

template <typename T>
int ttk::Dijkstra::shortestPath(const ttk::SimplexId source,
                                ttk::Triangulation &triangulation,
                                std::vector<T> &outputDists,
                                const std::vector<ttk::SimplexId> *bounds) {

  // should we process the whole mesh or stop at some point?
  bool processAllVertices = (bounds == nullptr);
  // total number of vertices in the mesh
  size_t vertexNumber = triangulation.getNumberOfVertices();

  // list all reached bounds
  std::vector<bool> reachedBounds;

  if(!processAllVertices) {
    reachedBounds.resize(bounds->size(), false);
  }

  // map vertex TTK id -> if already visited
  std::vector<bool> visited(vertexNumber);

  // link vertex and current distance to source
  using pq_t = std::pair<T, SimplexId>;

  // priority queue storing pairs of (distance, vertices TTK id)
  std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> pq;

  // map TTK id to current distance to source
  std::vector<T> dists(vertexNumber, std::numeric_limits<T>::infinity());

  // init pipeline
  pq.push(std::make_pair(T(0.0f), source));
  dists[source] = 0.0f;

  while(!pq.empty()) {
    auto elem = pq.top();
    pq.pop();
    auto vert = elem.second;
    float vCoords[3];
    triangulation.getVertexPoint(vert, vCoords[0], vCoords[1], vCoords[2]);

    auto nneigh = triangulation.getVertexNeighborNumber(vert);

    for(SimplexId i = 0; i < nneigh; i++) {
      // neighbor Id
      SimplexId neigh{};
      triangulation.getVertexNeighbor(vert, i, neigh);
      // neighbor coordinates
      float nCoords[3];
      triangulation.getVertexPoint(neigh, nCoords[0], nCoords[1], nCoords[2]);
      // (square) distance between vertex and neighbor
      T distVN = Geometry::distance(&vCoords[0], &nCoords[0]);
      if(dists[neigh] > dists[vert] + distVN) {
        dists[neigh] = dists[vert] + distVN;
        if(!processAllVertices) {
          // check if neigh in bounds
          auto it = std::find(bounds->begin(), bounds->end(), neigh);
          if(it != bounds->end()) {
            // mark it as found
            reachedBounds[it - bounds->begin()] = true;
          }
          // break if all are found
          if(std::all_of(reachedBounds.begin(), reachedBounds.end(),
                         [](const bool v) { return v; })) {
            break;
          }
        }
        pq.push(std::make_pair(dists[neigh], neigh));
      }
    }
  }

  return 0;
}

// explicit intantiations for floating-point types
template int
  ttk::Dijkstra::shortestPath<float>(const ttk::SimplexId source,
                                     ttk::Triangulation &triangulation,
                                     std::vector<float> &outputDists,
                                     const std::vector<ttk::SimplexId> *bounds
                                     = nullptr);
template int
  ttk::Dijkstra::shortestPath<double>(const ttk::SimplexId source,
                                      ttk::Triangulation &triangulation,
                                      std::vector<double> &outputDists,
                                      const std::vector<ttk::SimplexId> *bounds
                                      = nullptr);
