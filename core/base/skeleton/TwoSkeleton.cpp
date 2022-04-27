#include <TwoSkeleton.h>
#include <boost/container/small_vector.hpp>

using namespace ttk;

TwoSkeleton::TwoSkeleton() {
  setDebugMsgPrefix("TwoSkeleton");
}

int TwoSkeleton::buildCellNeighborsFromEdges(
  const CellArray &cellArray,
  FlatJaggedArray &cellNeighbors,
  const FlatJaggedArray &edgeStars) const {

  Timer t;

  printMsg("Building cell neighbors", 0, 0, 1, debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  const SimplexId edgeNumber = edgeStars.subvectorsNumber();
  std::vector<SimplexId> offsets(cellNumber + 1);
  // number of neigbhors processed per cell
  std::vector<SimplexId> neighborsId(cellNumber);

  for(SimplexId i = 0; i < edgeNumber; i++) {
    if(edgeStars.size(i) == 2) {
      // tetra cells in edge i's star
      const auto cs0 = edgeStars.get(i, 0);
      const auto cs1 = edgeStars.get(i, 1);
      offsets[cs0 + 1]++;
      offsets[cs1 + 1]++;
    }
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> neighbors(offsets.back());

  // fill flat neighbors vector using offsets and neighbors count vectors
  for(SimplexId i = 0; i < edgeNumber; i++) {
    if(edgeStars.size(i) == 2) {
      // tetra cells in edge i's star
      const auto cs0 = edgeStars.get(i, 0);
      const auto cs1 = edgeStars.get(i, 1);
      neighbors[offsets[cs0] + neighborsId[cs0]] = cs1;
      neighborsId[cs0]++;
      neighbors[offsets[cs1] + neighborsId[cs1]] = cs0;
      neighborsId[cs1]++;
    }
  }

  // fill FlatJaggedArray struct
  cellNeighbors.setData(std::move(neighbors), std::move(offsets));

  printMsg("Built " + std::to_string(cellNumber) + " cell neighbors", 1,
           t.getElapsedTime(), 1);

  return 0;
}

int TwoSkeleton::buildCellNeighborsFromVertices(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  FlatJaggedArray &cellNeighbors,
  FlatJaggedArray *vertexStars) const {

  auto localVertexStars = vertexStars;
  FlatJaggedArray defaultVertexStars{};

  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if(localVertexStars->empty()) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.buildVertexStars(vertexNumber, cellArray, *localVertexStars);
  }

  Timer t;

  printMsg(
    "Building cell neighbors", 0, 0, threadNumber_, debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  using boost::container::small_vector;
  // for each cell/triangle, a vector of neighbors
  std::vector<small_vector<SimplexId, 3>> neighbors(cellNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);

    for(SimplexId j = 0; j < nbVertCell; j++) {

      SimplexId v0 = cellArray.getCellVertex(cid, j);
      SimplexId v1 = cellArray.getCellVertex(cid, (j + 1) % nbVertCell);

      // perform an intersection of the 2 sorted star lists
      SimplexId pos0 = 0, pos1 = 0;
      SimplexId intersection = -1;

      while(pos0 < localVertexStars->size(v0)
            && pos1 < localVertexStars->size(v1)) {

        SimplexId biggest = localVertexStars->get(v0, pos0);
        if(localVertexStars->get(v1, pos1) > biggest) {
          biggest = localVertexStars->get(v1, pos1);
        }

        for(SimplexId l = pos0; l < localVertexStars->size(v0); l++) {
          if(localVertexStars->get(v0, l) < biggest) {
            pos0++;
          } else {
            break;
          }
        }
        for(SimplexId l = pos1; l < localVertexStars->size(v1); l++) {
          if(localVertexStars->get(v1, l) < biggest) {
            pos1++;
          } else {
            break;
          }
        }

        if(localVertexStars->get(v0, pos0) == localVertexStars->get(v1, pos1)) {
          if(localVertexStars->get(v0, pos0) != cid) {
            intersection = localVertexStars->get(v0, pos0);
            break;
          }

          pos0++;
          pos1++;
        }
      }

      if(intersection != -1) {
        neighbors[cid].emplace_back(intersection);
      }
    }
  }

  // convert to a FlatJaggedArray
  cellNeighbors.fillFrom(neighbors);

  printMsg("Built " + std::to_string(cellNumber) + " cell neighbors", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildEdgeTriangles(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  FlatJaggedArray &edgeTriangleList,
  const std::vector<std::array<SimplexId, 2>> &edgeList,
  std::vector<std::array<SimplexId, 3>> *triangleEdgeList) const {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexNumber <= 0)
    return -1;
#endif

  auto localTriangleEdgeList = triangleEdgeList;
  std::vector<std::array<SimplexId, 3>> defaultTriangleEdgeList{};
  if(!localTriangleEdgeList) {
    localTriangleEdgeList = &defaultTriangleEdgeList;
  }

  if(localTriangleEdgeList->empty()) {
    buildTriangleEdgeList(
      vertexNumber, cellArray, *localTriangleEdgeList, edgeList);
  }

  const auto edgeNumber{edgeList.size()};

  std::vector<SimplexId> offsets(edgeNumber + 1);
  // number of neighbors processed per vertex
  std::vector<SimplexId> trianglesId(edgeNumber);

  Timer t;

  printMsg("Building edge triangles", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  // store number of triangles per edge
  for(const auto &te : *localTriangleEdgeList) {
    offsets[te[0] + 1]++;
    offsets[te[1] + 1]++;
    offsets[te[2] + 1]++;
  }

  // compute partial sum of number of triangles per edge
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat edge triangle vector
  std::vector<SimplexId> edgeTriangles(offsets.back());

  // fill flat neighbors vector using offsets and neighbors count vectors
  for(size_t i = 0; i < localTriangleEdgeList->size(); ++i) {
    const auto &te{(*localTriangleEdgeList)[i]};
    edgeTriangles[offsets[te[0]] + trianglesId[te[0]]] = i;
    trianglesId[te[0]]++;
    edgeTriangles[offsets[te[1]] + trianglesId[te[1]]] = i;
    trianglesId[te[1]]++;
    edgeTriangles[offsets[te[2]] + trianglesId[te[2]]] = i;
    trianglesId[te[2]]++;
  }

  // fill FlatJaggedArray struct
  edgeTriangleList.setData(std::move(edgeTriangles), std::move(offsets));

  printMsg("Built " + std::to_string(edgeNumber) + " edge triangles", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildTriangleList(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  std::vector<std::array<SimplexId, 3>> *triangleList,
  FlatJaggedArray *triangleStars,
  std::vector<std::array<SimplexId, 4>> *cellTriangleList) const {

  Timer t;

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexNumber <= 0)
    return -1;
  if((!triangleList) && (!triangleStars) && (!cellTriangleList)) {
    // we've got nothing to do here.
    return -2;
  }
#endif

  // check parameters consistency (this method is useless in 2D)
  const auto dim = cellArray.getCellVertexNumber(0) - 1;
  if(dim == 2) {
    this->printWrn("Calling buildTriangleList is useless in 2D, skipping...");
    return -1;
  }

  printMsg(
    "Building triangles", 0, 0, threadNumber_, ttk::debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();

  // we need cellTriangleList to compute triangleStars
  std::vector<std::array<SimplexId, 4>> defaultCellTriangleList{};
  if(triangleStars != nullptr && cellTriangleList == nullptr) {
    cellTriangleList = &defaultCellTriangleList;
  }

  if(cellTriangleList) {
    cellTriangleList->resize(cellNumber, {-1, -1, -1, -1});
  }

  struct TriangleData {
    // the two higher vertices id of the triangle
    std::array<SimplexId, 2> highVerts{};
    // the triangle id
    SimplexId id{-1};
    TriangleData(std::array<SimplexId, 2> hVerts, SimplexId i)
      : highVerts{hVerts}, id{i} {
    }
  };

  using boost::container::small_vector;
  // for each vertex, a vector of TriangleData
  std::vector<small_vector<TriangleData, 8>> triangleTable(vertexNumber);

  SimplexId nTriangles{};

  printMsg("Building triangles", 0.25, t.getElapsedTime(), 1,
           debug::LineMode::REPLACE);

  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    // a tetra cell has 4 faces
    for(size_t j = 0; j < 4; j++) {
      std::array<SimplexId, 3> triangle{};
      for(size_t k = 0; k < 3; k++) {
        // TODO: ASSUME Regular Mesh Here!
        triangle[k] = cellArray.getCellVertex(cid, (j + k) % 4);
      }
      std::sort(triangle.begin(), triangle.end());
      auto &ttable = triangleTable[triangle[0]];

      // check if current triangle already registered in triangleTable
      // via another tetra in its star
      bool found{false};
      for(auto &d : ttable) {
        if(d.highVerts[0] == triangle[1] && d.highVerts[1] == triangle[2]) {
          found = true;
          if(cellTriangleList != nullptr) {
            (*cellTriangleList)[cid][j] = d.id;
          }
          break;
        }
      }
      if(!found) {
        // new triangle added
        ttable.emplace_back(
          TriangleData{{triangle[1], triangle[2]}, nTriangles});
        if(cellTriangleList != nullptr) {
          (*cellTriangleList)[cid][j] = nTriangles;
        }
        nTriangles++;
      }
    }
  }

  printMsg(
    "Building triangles", 0.5, t.getElapsedTime(), 1, debug::LineMode::REPLACE);

  // resize vectors to the correct size
  if(triangleList) {
    triangleList->resize(nTriangles);
  }

  // fill data buffers in parallel

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < vertexNumber; ++i) {
    const auto &ttable = triangleTable[i];
    for(const auto &data : ttable) {
      if(triangleList != nullptr) {
        (*triangleList)[data.id] = {i, data.highVerts[0], data.highVerts[1]};
      }
    }
  }

  printMsg("Building triangles", 0.75, t.getElapsedTime(), 1,
           debug::LineMode::REPLACE);

  if(cellTriangleList != nullptr && triangleStars != nullptr) {
    std::vector<SimplexId> offsets(nTriangles + 1);
    // number of cells processed per vertex
    std::vector<SimplexId> starIds(nTriangles);

    // store number of cells per triangle
    for(const auto &c : *cellTriangleList) {
      offsets[c[0] + 1]++;
      offsets[c[1] + 1]++;
      offsets[c[2] + 1]++;
      offsets[c[3] + 1]++;
    }

    // compute partial sum of number of cells per triangle
    for(size_t i = 1; i < offsets.size(); ++i) {
      offsets[i] += offsets[i - 1];
    }

    // allocate flat triangle stars vector
    std::vector<SimplexId> triangleSt(offsets.back());

    // fill flat neighbors vector using offsets and neighbors count vectors
    for(size_t i = 0; i < cellTriangleList->size(); ++i) {
      const auto &ct{(*cellTriangleList)[i]};
      triangleSt[offsets[ct[0]] + starIds[ct[0]]] = i;
      starIds[ct[0]]++;
      triangleSt[offsets[ct[1]] + starIds[ct[1]]] = i;
      starIds[ct[1]]++;
      triangleSt[offsets[ct[2]] + starIds[ct[2]]] = i;
      starIds[ct[2]]++;
      triangleSt[offsets[ct[3]] + starIds[ct[3]]] = i;
      starIds[ct[3]]++;
    }

    // fill FlatJaggedArray struct
    triangleStars->setData(std::move(triangleSt), std::move(offsets));
  }

  printMsg("Built " + std::to_string(nTriangles) + " triangles", 1,
           t.getElapsedTime(), 1);
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)// 1 thread: 58.5631 s//
  // 24 threads: 87.5816 s (~) ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - no
  // tet adj// 1 thread: 5.7427 s// 4 threads: 9.14764 s (~)
  // ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - tet adjacency// 1
  // thread: 8.59854 s// 4 threads: 15.837 s
  // (~) ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - tet adjacency only// 1
  // thread: 6.85897 s// 4 threads: 14.78 s (~)
  return 0;
}

int TwoSkeleton::buildTriangleEdgeList(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  std::vector<std::array<SimplexId, 3>> &triangleEdgeList,
  const std::vector<std::array<SimplexId, 2>> &edgeList,
  FlatJaggedArray *vertexEdgeList,
  std::vector<std::array<SimplexId, 3>> *triangleList,
  FlatJaggedArray *triangleStarList,
  std::vector<std::array<SimplexId, 4>> *cellTriangleList) const {

  auto localVertexEdgeList = vertexEdgeList;
  FlatJaggedArray defaultVertexEdgeList{};
  if(!localVertexEdgeList) {
    localVertexEdgeList = &defaultVertexEdgeList;
  }

  if(localVertexEdgeList->empty()) {

    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.setThreadNumber(threadNumber_);

    zeroSkeleton.buildVertexEdges(
      vertexNumber, edgeList, (*localVertexEdgeList));
  }

  // NOTE:
  // triangleStarList and cellTriangleList: we do not need these guys but we
  // can compute them for free optionally.

  auto localTriangleList = triangleList;
  std::vector<std::array<SimplexId, 3>> defaultTriangleList{};
  if(!localTriangleList) {
    localTriangleList = &defaultTriangleList;
  }
  if(!localTriangleList->size()) {

    buildTriangleList(vertexNumber, cellArray, localTriangleList,
                      triangleStarList, cellTriangleList);
  }

  triangleEdgeList.resize(localTriangleList->size());

  Timer tm;

  printMsg("Building triangle edges", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  // now for each triangle, grab its vertices, add the edges in the triangle
  // with no duplicate
  // let's do the real stuff
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < localTriangleList->size(); i++) {
    const auto &t = (*localTriangleList)[i];
    const auto beg0 = localVertexEdgeList->get_ptr(t[0], 0);
    const auto end0 = beg0 + localVertexEdgeList->size(t[0]);
    const auto beg1 = localVertexEdgeList->get_ptr(t[1], 0);
    const auto end1 = beg1 + localVertexEdgeList->size(t[1]);
    const auto beg2 = localVertexEdgeList->get_ptr(t[2], 0);
    const auto end2 = beg2 + localVertexEdgeList->size(t[2]);
    std::set_intersection(beg0, end0, beg1, end1, &triangleEdgeList[i][0]);
    std::set_intersection(beg0, end0, beg2, end2, &triangleEdgeList[i][1]);
    std::set_intersection(beg1, end1, beg2, end2, &triangleEdgeList[i][2]);
    std::sort(triangleEdgeList[i].begin(), triangleEdgeList[i].end());
  }

  SimplexId triangleNumber = localTriangleList->size();

  printMsg("Built " + std::to_string(triangleNumber) + " triangle edges", 1,
           tm.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildTriangleLinks(
  const std::vector<std::array<SimplexId, 3>> &triangleList,
  const FlatJaggedArray &triangleStars,
  const CellArray &cellArray,
  FlatJaggedArray &triangleLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleList.empty())
    return -1;
  if((triangleStars.empty())
     || (triangleStars.subvectorsNumber() != triangleList.size()))
    return -2;
#endif

  Timer tm;

  // check parameters consistency (this method is useless in 2D)
  const auto dim = cellArray.getCellVertexNumber(0) - 1;
  if(dim == 2) {
    this->printWrn("Calling buildTriangleLinks is useless in 2D, skipping...");
    return -1;
  }

  const SimplexId triangleNumber = triangleList.size();
  std::vector<SimplexId> offsets(triangleNumber + 1);
  // one vertex per star
  std::vector<SimplexId> links(triangleStars.dataSize());

  printMsg("Building triangle links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < triangleNumber; i++) {

    // copy the triangleStars offsets array
    offsets[i] = triangleStars.offset(i);

    const auto &t = triangleList[i];
    for(SimplexId j = 0; j < triangleStars.size(i); j++) {
      // for each tetra in triangle i's star, get the opposite vertex
      for(size_t k = 0; k < 4; k++) {
        const auto v = cellArray.getCellVertex(triangleStars.get(i, j), k);
        if(v != t[0] && v != t[1] && v != t[2]) {
          links[offsets[i] + j] = v;
          break;
        }
      }
    }
  }

  // don't forget the last offset
  offsets[triangleNumber] = triangleStars.offset(triangleNumber);

  triangleLinks.setData(std::move(links), std::move(offsets));

  printMsg("Built " + std::to_string(triangleNumber) + " triangle links", 1,
           tm.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildVertexTriangles(
  const SimplexId &vertexNumber,
  const std::vector<std::array<SimplexId, 3>> &triangleList,
  FlatJaggedArray &vertexTriangles) const {

  Timer tm;

  printMsg("Building vertex triangles", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::vector<SimplexId> offsets(vertexNumber + 1);
  // number of triangles processed per vertex
  std::vector<SimplexId> trianglesId(vertexNumber);

  // store number of triangles per vertex
  for(const auto &t : triangleList) {
    offsets[t[0] + 1]++;
    offsets[t[1] + 1]++;
    offsets[t[2] + 1]++;
  }

  // compute partial sum of number of neighbors per vertex
  for(size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  // allocate flat neighbors vector
  std::vector<SimplexId> data(offsets.back());

  // fill flat data vector using offsets and triangles count vectors
  for(size_t i = 0; i < triangleList.size(); ++i) {
    const auto &t{triangleList[i]};
    data[offsets[t[0]] + trianglesId[t[0]]] = i;
    trianglesId[t[0]]++;
    data[offsets[t[1]] + trianglesId[t[1]]] = i;
    trianglesId[t[1]]++;
    data[offsets[t[2]] + trianglesId[t[2]]] = i;
    trianglesId[t[2]]++;
  }

  // fill FlatJaggedArray struct
  vertexTriangles.setData(std::move(data), std::move(offsets));

  printMsg("Built " + std::to_string(vertexNumber) + " vertex triangles", 1,
           tm.getElapsedTime(), 1);

  return 0;
}
