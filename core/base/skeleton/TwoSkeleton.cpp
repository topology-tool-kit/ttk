#include <TwoSkeleton.h>

using namespace std;
using namespace ttk;

TwoSkeleton::TwoSkeleton() {

  setDebugMsgPrefix("TwoSkeleton");
}

TwoSkeleton::~TwoSkeleton() {
}

int TwoSkeleton::buildCellNeighborsFromVertices(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &cellNeighbors,
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
  cellNeighbors.resize(cellNumber);
  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);
    cellNeighbors[cid].reserve(nbVertCell);
  }

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
        cellNeighbors[cid].emplace_back(intersection);
      }
    }
  }

  printMsg("Built " + to_string(cellNumber) + " cell neighbors", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildEdgeTriangles(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &edgeTriangleList,
  vector<std::array<SimplexId, 2>> *edgeList,
  vector<std::array<SimplexId, 3>> *triangleList,
  FlatJaggedArray *vertexTriangleList) const {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexNumber <= 0)
    return -1;
#endif

  auto localEdgeList = edgeList;
  vector<std::array<SimplexId, 2>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  auto localVertexTriangleList = vertexTriangleList;
  FlatJaggedArray defaultVertexTriangleList{};
  if(!localVertexTriangleList) {
    localVertexTriangleList = &defaultVertexTriangleList;
  }

  auto localTriangleList = triangleList;
  vector<std::array<SimplexId, 3>> defaultTriangleList{};
  if(!localTriangleList) {
    localTriangleList = &defaultTriangleList;
  }

  OneSkeleton oneSkeleton;
  oneSkeleton.setDebugLevel(debugLevel_);
  oneSkeleton.setThreadNumber(threadNumber_);

  // now do the pre-computation
  if(localEdgeList->empty()) {
    oneSkeleton.buildEdgeList(vertexNumber, cellArray, (*localEdgeList));
  }

  if(localTriangleList->empty()) {
    buildTriangleList(vertexNumber, cellArray, localTriangleList);
  }

  if(localVertexTriangleList->empty()) {
    buildVertexTriangles(
      vertexNumber, (*localTriangleList), (*localVertexTriangleList));
  }

  edgeTriangleList.resize(localEdgeList->size());

  Timer t;

  printMsg("Building edge triangles", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < localEdgeList->size(); i++) {
    const auto &e = (*localEdgeList)[i];
    const auto beg0 = localVertexTriangleList->get_ptr(e[0], 0);
    const auto end0 = beg0 + localVertexTriangleList->size(e[0]);
    const auto beg1 = localVertexTriangleList->get_ptr(e[1], 0);
    const auto end1 = beg1 + localVertexTriangleList->size(e[1]);
    // merge the two vertex triangle lists
    std::set_intersection(
      beg0, end0, beg1, end1, std::back_inserter(edgeTriangleList[i]));
  }

  SimplexId edgeNumber = localEdgeList->size();

  printMsg("Built " + to_string(edgeNumber) + " edge triangles", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildTriangleList(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<std::array<SimplexId, 3>> *triangleList,
  vector<vector<SimplexId>> *triangleStars,
  vector<std::array<SimplexId, 4>> *cellTriangleList) const {

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

  printMsg(
    "Building triangles", 0, 0, threadNumber_, ttk::debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  if(cellTriangleList) {
    cellTriangleList->resize(cellNumber, {-1, -1, -1, -1});
  }

  struct TriangleData {
    // the two higher vertices id of the triangle
    std::array<SimplexId, 2> highVerts{};
    // the triangle id
    SimplexId id{-1};
    // tetra ids that have current triangle as face
    std::array<SimplexId, 2> cellIds{-1, -1};
    TriangleData(std::array<SimplexId, 2> hVerts, SimplexId i, SimplexId cellId)
      : highVerts{hVerts}, id{i}, cellIds{cellId, -1} {
    }
  };

  // for each vertex, a vector of TriangleData
  std::vector<std::vector<TriangleData>> triangleTable(vertexNumber);

  printMsg("Building triangles", 0.12, t.getElapsedTime(), 1,
           ttk::debug::LineMode::REPLACE);

  for(SimplexId i = 0; i < vertexNumber; i++) {
    triangleTable[i].reserve(16);
  }

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
          // hypothesis: at most two tetras sharing the same triangle
          d.cellIds[1] = cid;
          if(cellTriangleList != nullptr) {
            (*cellTriangleList)[cid][j] = d.id;
          }
          break;
        }
      }
      if(!found) {
        // new triangle added
        ttable.emplace_back(
          TriangleData{{triangle[1], triangle[2]}, nTriangles, cid});
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

  printMsg("Building triangles", 0.67, t.getElapsedTime(), 1,
           debug::LineMode::REPLACE);

  if(triangleStars) {
    triangleStars->resize(nTriangles);
    for(auto &v : *triangleStars) {
      v.reserve(2);
    }
  }

  printMsg("Building triangles", 0.75, t.getElapsedTime(), 1,
           debug::LineMode::REPLACE);

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
      for(size_t j = 0; j < 2; ++j) {
        if(data.cellIds[j] == -1) {
          break;
        }
        if(triangleStars != nullptr) {
          (*triangleStars)[data.id].emplace_back(data.cellIds[j]);
        }
      }
    }
  }

  printMsg(
    "Built " + to_string(nTriangles) + " triangles", 1, t.getElapsedTime(), 1);
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
  vector<std::array<SimplexId, 3>> &triangleEdgeList,
  FlatJaggedArray *vertexEdgeList,
  vector<std::array<SimplexId, 2>> *edgeList,
  vector<std::array<SimplexId, 3>> *triangleList,
  vector<vector<SimplexId>> *triangleStarList,
  vector<std::array<SimplexId, 4>> *cellTriangleList) const {

  auto localEdgeList = edgeList;
  vector<std::array<SimplexId, 2>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {
    OneSkeleton oneSkeleton;
    oneSkeleton.setDebugLevel(debugLevel_);
    oneSkeleton.setThreadNumber(threadNumber_);

    oneSkeleton.buildEdgeList(vertexNumber, cellArray, (*localEdgeList));
  }

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
      vertexNumber, (*localEdgeList), (*localVertexEdgeList));
  }

  // NOTE:
  // triangleStarList and cellTriangleList: we do not need these guys but we
  // can compute them for free optionally.

  auto localTriangleList = triangleList;
  vector<std::array<SimplexId, 3>> defaultTriangleList{};
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

  printMsg("Built " + to_string(triangleNumber) + " triangle edges", 1,
           tm.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildTriangleLinks(
  const vector<std::array<SimplexId, 3>> &triangleList,
  const vector<vector<SimplexId>> &triangleStars,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &triangleLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleList.empty())
    return -1;
  if((triangleStars.empty()) || (triangleStars.size() != triangleList.size()))
    return -2;
#endif

  Timer t;

  printMsg("Building triangle links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  triangleLinks.resize(triangleList.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < triangleList.size(); i++) {

    for(size_t j = 0; j < triangleStars[i].size(); j++) {

      for(size_t k = 0; k < 4; k++) {
        SimplexId vertexId = cellArray.getCellVertex(triangleStars[i][j], k);

        if((vertexId != triangleList[i][0]) && (vertexId != triangleList[i][1])
           && (vertexId != triangleList[i][2])) {
          triangleLinks[i].push_back(vertexId);
          break;
        }
      }
    }
  }

  printMsg("Built " + std::to_string(triangleList.size()) + " triangle links",
           1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int TwoSkeleton::buildVertexTriangles(
  const SimplexId &vertexNumber,
  const vector<std::array<SimplexId, 3>> &triangleList,
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
