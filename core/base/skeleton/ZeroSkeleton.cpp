#include <OneSkeleton.h>
#include <ZeroSkeleton.h>

using namespace std;
using namespace ttk;

ZeroSkeleton::ZeroSkeleton() {
  setDebugMsgPrefix("ZeroSkeleton");
}

int ZeroSkeleton::buildVertexEdges(
  const SimplexId &vertexNumber,
  const vector<std::array<SimplexId, 2>> &edgeList,
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

int ZeroSkeleton::buildVertexLink(const SimplexId &vertexId,
                                  const CellArray &cellArray,
                                  vector<LongSimplexId> &vertexLink) const {

  const LongSimplexId cellNumber = cellArray.getNbCells();
  LongSimplexId totalLinkSize = 0;
  vector<SimplexId> vertexStar;
  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);
    for(SimplexId j = 1; j < nbVertCell; j++) {
      if(cellArray.getCellVertex(cid, j) == vertexId) {
        vertexStar.emplace_back(cid);
        totalLinkSize += cellArray.getCellVertexNumber(cid);
        break;
      }
    }
  }

  vertexLink.resize(totalLinkSize);
  vector<SimplexId> faceIds;
  const SimplexId nbCellInStar = vertexStar.size();

  // construct the link
  for(SimplexId i = 0; i < nbCellInStar; i++) {
    const SimplexId cellId = vertexStar[i];
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cellId);
    faceIds.resize(nbVertCell);
    faceIds[0] = nbVertCell; // first entry is number of vertices

    bool hasPivotVertex = false;

    // iterate on the cell's faces
    for(int k = 0; k < 2; k++) {
      faceIds[1] = cellArray.getCellVertex(cellId, k);

      if(faceIds[1] != vertexId) {

        if(nbVertCell > 2) {

          for(SimplexId l = k + 1; l <= nbVertCell - 1; l++) {
            faceIds[2] = cellArray.getCellVertex(cellId, l);

            if(faceIds[2] != vertexId) {

              if(nbVertCell == 4) {
                // tet case, faceIds has 4 entries to fill
                for(SimplexId m = l + 1; m < nbVertCell; m++) {
                  faceIds[3] = cellArray.getCellVertex(cellId, m);

                  if(faceIds[3] != vertexId) {

                    // all the vertices of the face are different from our
                    // vertex. let's add that face to the link
                    for(size_t n = 0; n < faceIds.size(); n++) {
                      vertexLink[i * (nbVertCell) + n] = faceIds[n];
                    }
                    hasPivotVertex = true;
                    break;
                  }
                }
              } else if(nbVertCell == 3) {
                // triangle case
                // we're holding to an edge that does not contain our vertex
                for(size_t n = 0; n < faceIds.size(); n++) {
                  vertexLink[i * (nbVertCell) + n] = faceIds[n];
                }
                hasPivotVertex = true;
                break;
              }
              if(hasPivotVertex)
                break;
            }
          }
        } else if(nbVertCell == 2) {
          // edge-mesh case
          // we're holding a neighbor different from our vertex
          for(size_t n = 0; n < faceIds.size(); n++) {
            vertexLink[i * (nbVertCell) + n] = faceIds[n];
          }
          break;
        }
      }
      if(hasPivotVertex)
        break;
    }
  }

  return 0;
}

int ZeroSkeleton::buildVertexLinks(const SimplexId &vertexNumber,
                                   const CellArray &cellArray,
                                   vector<vector<LongSimplexId>> &vertexLinks,
                                   FlatJaggedArray *vertexStars) const {

  auto localVertexStars = vertexStars;
  FlatJaggedArray defaultVertexStars{};
  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if((SimplexId)localVertexStars->subvectorsNumber() != vertexNumber) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.buildVertexStars(vertexNumber, cellArray, *localVertexStars);
  }

  Timer t;

  printMsg("Building vertex links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  const SimplexId nbVertLinks = vertexLinks.size();
  if(nbVertLinks != vertexNumber) {
    vertexLinks.resize(vertexNumber);
    for(SimplexId i = 0; i < nbVertLinks; i++) {
      SimplexId nbVertsInLinkI = 0;
      for(SimplexId j = 0; j < localVertexStars->size(i); ++j) {
        const auto cid = localVertexStars->get(i, j);
        nbVertsInLinkI += cellArray.getCellVertexNumber(cid);
      }
      vertexLinks[i].resize(nbVertsInLinkI);
    }
  }

  vector<vector<SimplexId>> faceIds(threadNumber_);
  for(ThreadId i = 0; i < threadNumber_; i++) {
    faceIds[i].reserve(cellArray.getCellVertexNumber(i));
  }
  // NOTE:
  // 1-thread:
  // memory allocation: 2.15 s.
  // processing: 1 s. [< 33%]
  // 8-thread (4 cores): 0.38 s. [< 18%] total 2.53
  // processing:

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber; i++) {

    ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    const SimplexId nbCellInStarI = localVertexStars->size(i);
    for(SimplexId j = 0; j < nbCellInStarI; j++) {

      const SimplexId cellId = localVertexStars->get(i, j);
      const SimplexId nbVertCell = cellArray.getCellVertexNumber(cellId);

      faceIds[threadId].resize(nbVertCell); // no realloc here
      // the first entry should be the number of vertices
      // tetrahedra (4) => link made of triangles (3)
      // triangle (3) => link made of edges (2)
      faceIds[threadId][0] = nbVertCell - 1;

      // tet case (4)
      // 0 - 1 - 2
      // 0 - 1 - 3
      // 0 - 2 - 3
      // 1 - 2 - 3
      // triangle case (3)
      // 0 - 1
      // 0 - 2
      // 1 - 2
      // edge case (2)
      // 0
      // 1
      bool hasPivotVertex = false;

      // iterate on the cell's faces
      for(int k = 0; k < 2; k++) {

        faceIds[threadId][1] = cellArray.getCellVertex(cellId, k);

        if(faceIds[threadId][1] != i) {

          if(nbVertCell > 2) {

            for(SimplexId l = k + 1; l <= nbVertCell - 1; l++) {
              faceIds[threadId][2] = cellArray.getCellVertex(cellId, l);

              if(faceIds[threadId][2] != i) {

                if(nbVertCell == 4) {
                  // tet case, faceIds[threadId] has 4 entries to fill
                  for(SimplexId m = l + 1; m < nbVertCell; m++) {
                    faceIds[threadId][3] = cellArray.getCellVertex(cellId, m);

                    // now test if this face contains our vertex or not
                    // there's should be only one face
                    if(faceIds[threadId][3] != i) {
                      // all the vertices of the face are different from our
                      // pivot vertex.
                      // let's add that face to the link
                      for(SimplexId n = 0; n < 4; n++) {
                        // TODO: ASSUME Regular Mesh Here !
                        vertexLinks[i][j * nbVertCell + n]
                          = faceIds[threadId][n];
                      }
                      hasPivotVertex = true;
                      break;
                    }
                  }
                } else if(nbVertCell == 3) {
                  // triangle case
                  // we're holding to an edge that does not contain our vertex
                  for(SimplexId n = 0; n < 3; n++) {
                    vertexLinks[i][j * nbVertCell + n] = faceIds[threadId][n];
                  }
                  break;
                }
                if(hasPivotVertex)
                  break;
              }
            }
          } else if(nbVertCell == 2) {
            // edge-mesh case
            // we holding a neighbor different from us
            for(SimplexId n = 0; n < 2; n++) {
              vertexLinks[i][j * nbVertCell + n] = faceIds[threadId][n];
            }
            break;
          }
        }
        if(hasPivotVertex)
          break;
      }
    }
  }

  if(debugLevel_ >= (int)(debug::Priority::DETAIL)) {
    for(size_t i = 0; i < vertexLinks.size(); i++) {
      stringstream msg;
      // TODO: ASSUME Regular Mesh Here!
      msg << "Vertex #" << i << " ("
          << vertexLinks[i].size() / (vertexLinks[i][0] + 1)
          << " items, length: " << vertexLinks[i].size() << "): ";
      printMsg(msg.str(), debug::Priority::DETAIL);
      for(size_t j = 0; j < vertexLinks[i].size() / (vertexLinks[i][0] + 1);
          j++) {
        stringstream msgLine;
        msgLine << "  - " << j << ":";
        const SimplexId nbVertCell = cellArray.getCellVertexNumber(j);
        for(SimplexId k = 0; k < nbVertCell; k++) {
          // TODO: ASSUME Regular Mesh Here!
          msgLine << " " << vertexLinks[i][j * nbVertCell + k];
        }
        printMsg(msgLine.str(), debug::Priority::DETAIL);
      }
    }
  }

  printMsg("Built " + std::to_string(vertexNumber) + " vertex links", 1,
           t.getElapsedTime(), threadNumber_);

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 10.47 s
  // 24 threads: 6.25 s

  return 0;
}

// 2D cells (triangles)
int ZeroSkeleton::buildVertexLinks(
  const FlatJaggedArray &vertexStars,
  const vector<std::array<SimplexId, 3>> &cellEdges,
  const vector<std::array<SimplexId, 2>> &edgeList,
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
  const vector<std::array<SimplexId, 4>> &cellTriangles,
  const vector<std::array<SimplexId, 3>> &triangleList,
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
  const CellArray &cellArray,
  FlatJaggedArray &vertexNeighbors,
  vector<std::array<SimplexId, 2>> *edgeList) const {

  std::vector<SimplexId> offsets(vertexNumber + 1);
  // number of neighbors processed per vertex
  std::vector<SimplexId> neighborsId(vertexNumber);

  auto localEdgeList = edgeList;
  vector<std::array<SimplexId, 2>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {
    OneSkeleton osk;
    osk.setDebugLevel(debugLevel_);
    osk.setThreadNumber(threadNumber_);
    osk.buildEdgeList(vertexNumber, cellArray, localEdgeList);
  }

  Timer t;

  printMsg("Building vertex neighbors", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // store number of neighbors per vertex
  for(const auto &e : *localEdgeList) {
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
  for(const auto &e : *localEdgeList) {
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

  if(debugLevel_ >= static_cast<int>(debug::Priority::VERBOSE)) {
    for(size_t i = 0; i < vertexStars.subvectorsNumber(); i++) {
      stringstream msg;
      msg << "Vertex #" << i << " (" << vertexStars.size(i) << " cell(s)): ";
      for(SimplexId j = 0; j < vertexStars.size(i); j++) {
        msg << " " << vertexStars.get(i, j);
      }
      printMsg(msg.str(), debug::Priority::VERBOSE);
    }
  }

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 0.53 s
  // 24 threads: 7.99 s

  return 0;
}
