#include <ZeroSkeleton.h>

using namespace ttk;

ZeroSkeleton::ZeroSkeleton() {
  setDebugMsgPrefix("ZeroSkeleton");
}

int ZeroSkeleton::buildVertexEdges(
  const SimplexId &vertexNumber,
  const std::vector<std::array<SimplexId, 2>> &edgeList,
  FlatJaggedArray &vertexEdges) const {

  std::vector<SimplexId> offsets(vertexNumber + 1);
  // number of edges processed per vertex
  std::vector<SimplexId> edgesId(vertexNumber);

  Timer t;

  printMsg("Building vertex edges", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // store number of edges per vertex
  for(const auto &e : edgeList) {
    offsets[e[0] + 1]++;
    offsets[e[1] + 1]++;
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> data(offsets.back());

  // fill flat data vector using offsets and edges count vectors
  for(size_t i = 0; i < edgeList.size(); ++i) {
    const auto &e{edgeList[i]};
    data[offsets[e[0]] + edgesId[e[0]]] = i;
    edgesId[e[0]]++;
    data[offsets[e[1]] + edgesId[e[1]]] = i;
    edgesId[e[1]]++;
  }

  // fill FlatJaggedArray struct
  vertexEdges.setData(std::move(data), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex edges", 1,
           t.getElapsedTime(), 1);

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 11.85 s
  // 24 threads: 20.93 s [not efficient]

  return 0;
}

// 2D cells (triangles)
int ZeroSkeleton::buildVertexLinks(
  const FlatJaggedArray &vertexStars,
  const std::vector<std::array<SimplexId, 3>> &cellEdges,
  const std::vector<std::array<SimplexId, 2>> &edgeList,
  FlatJaggedArray &vertexLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexStars.empty())
    return -1;
  if(cellEdges.empty())
    return -2;
  if(edgeList.empty())
    return -3;
#endif

  Timer t;

  const SimplexId vertexNumber = vertexStars.subvectorsNumber();
  std::vector<SimplexId> offsets(vertexNumber + 1);
  // one edge per star
  std::vector<SimplexId> links(vertexStars.dataSize());

  printMsg("Building vertex links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; i++) {

    // copy the vertexStars offsets array
    offsets[i] = vertexStars.offset(i);

    for(SimplexId j = 0; j < vertexStars.size(i); j++) {
      // for each cell/triangle in vertex i's star, get the opposite edge
      for(size_t k = 0; k < cellEdges[vertexStars.get(i, j)].size(); k++) {
        const auto e = cellEdges[vertexStars.get(i, j)][k];
        if(i != edgeList[e][0] && i != edgeList[e][1]) {
          links[offsets[i] + j] = e;
          break;
        }
      }
    }
  }

  // don't forget the last offset
  offsets[vertexNumber] = vertexStars.offset(vertexNumber);

  vertexLinks.setData(std::move(links), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex links", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

// 3D cells (tetrahedron)
int ZeroSkeleton::buildVertexLinks(
  const FlatJaggedArray &vertexStars,
  const std::vector<std::array<SimplexId, 4>> &cellTriangles,
  const std::vector<std::array<SimplexId, 3>> &triangleList,
  FlatJaggedArray &vertexLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexStars.empty())
    return -1;
  if(cellTriangles.empty())
    return -2;
  if(triangleList.empty())
    return -3;
#endif

  Timer tm;

  const SimplexId vertexNumber = vertexStars.subvectorsNumber();
  std::vector<SimplexId> offsets(vertexNumber + 1);
  // one triangle per star
  std::vector<SimplexId> links(vertexStars.dataSize());

  printMsg("Building vertex links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; i++) {

    // copy the vertexStars offsets array
    offsets[i] = vertexStars.offset(i);

    for(SimplexId j = 0; j < vertexStars.size(i); j++) {
      // for each cell/tetra in vertex i's star, get the opposite triangle
      for(size_t k = 0; k < cellTriangles[vertexStars.get(i, j)].size(); k++) {
        const auto t = cellTriangles[vertexStars.get(i, j)][k];
        if(i != triangleList[t][0] && i != triangleList[t][1]
           && i != triangleList[t][2]) {
          links[offsets[i] + j] = t;
          break;
        }
      }
    }
  }

  // don't forget the last offset
  offsets[vertexNumber] = vertexStars.offset(vertexNumber);

  vertexLinks.setData(std::move(links), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex links", 1,
           tm.getElapsedTime(), threadNumber_);

  return 0;
}

int ZeroSkeleton::buildVertexNeighbors(
  const SimplexId &vertexNumber,
  FlatJaggedArray &vertexNeighbors,
  const std::vector<std::array<SimplexId, 2>> &edgeList) const {

  std::vector<SimplexId> offsets(vertexNumber + 1);
  // number of neighbors processed per vertex
  std::vector<SimplexId> neighborsId(vertexNumber);

  Timer t;

  printMsg("Building vertex neighbors", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // store number of neighbors per vertex
  for(const auto &e : edgeList) {
    offsets[e[0] + 1]++;
    offsets[e[1] + 1]++;
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> neighbors(offsets.back());

  // fill flat neighbors vector using offsets and neighbors count vectors
  for(const auto &e : edgeList) {
    neighbors[offsets[e[0]] + neighborsId[e[0]]] = e[1];
    neighborsId[e[0]]++;
    neighbors[offsets[e[1]] + neighborsId[e[1]]] = e[0];
    neighborsId[e[1]]++;
  }

  // fill FlatJaggedArray struct
  vertexNeighbors.setData(std::move(neighbors), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex neighbors", 1,
           t.getElapsedTime(), 1);

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // (only merging step, after edge list creation)
  // 1 thread: 9.16 s
  // 24 threads: 13.21 s [not efficient in parallel]

  return 0;
}

int ZeroSkeleton::buildVertexStars(const SimplexId &vertexNumber,
                                   const CellArray &cellArray,
                                   FlatJaggedArray &vertexStars) const {

  Timer t;

  printMsg("Building vertex stars", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::vector<SimplexId> offsets(vertexNumber + 1);
  // number of cells processed per vertex
  std::vector<SimplexId> cellIds(vertexNumber);

  const auto cellNumber = cellArray.getNbCells();

  // store number of stars per vertex
  for(SimplexId i = 0; i < cellNumber; ++i) {
    const auto nbVertCell = cellArray.getCellVertexNumber(i);
    for(SimplexId j = 0; j < nbVertCell; ++j) {
      offsets[cellArray.getCellVertex(i, j) + 1]++;
    }
  }

  // compute partial sum of number of stars per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat data vector
  std::vector<SimplexId> data(offsets.back());

  // fill flat data vector using offsets and edges count vectors
  for(SimplexId i = 0; i < cellNumber; ++i) {
    const auto nbVertCell = cellArray.getCellVertexNumber(i);
    for(SimplexId j = 0; j < nbVertCell; ++j) {
      const auto v = cellArray.getCellVertex(i, j);
      data[offsets[v] + cellIds[v]] = i;
      cellIds[v]++;
    }
  }

  // fill FlatJaggedArray struct
  vertexStars.setData(std::move(data), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex stars", 1,
           t.getElapsedTime(), 1);

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 0.53 s
  // 24 threads: 7.99 s

  return 0;
}
