#include <OneSkeleton.h>
#include <ZeroSkeleton.h>

using namespace std;
using namespace ttk;

ZeroSkeleton::ZeroSkeleton() {
  setDebugMsgPrefix("ZeroSkeleton");
}

ZeroSkeleton::~ZeroSkeleton() {
}

int ZeroSkeleton::buildVertexEdges(
  const SimplexId &vertexNumber,
  const vector<pair<SimplexId, SimplexId>> &edgeList,
  vector<vector<SimplexId>> &vertexEdges) const {

  Timer t;

  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  ThreadId oldThreadNumber = threadNumber_;
  threadNumber_ = 1;

  vertexEdges.resize(vertexNumber);
  for(SimplexId i = 0; i < (SimplexId)vertexEdges.size(); i++) {
    vertexEdges[i].clear();
  }

  if(threadNumber_ == 1) {

    int timeBuckets = 10;
    if(timeBuckets > (int)edgeList.size()) {
      timeBuckets = edgeList.size();
    }

    printMsg("Building vertex edges", 0, 0, threadNumber_,
             ttk::debug::LineMode::REPLACE);

    for(SimplexId i = 0; i < (SimplexId)edgeList.size(); i++) {
      vertexEdges[edgeList[i].first].push_back(i);
      vertexEdges[edgeList[i].second].push_back(i);

      if(debugLevel_ >= (int)(debug::Priority::INFO)) {
        if(!(i % ((edgeList.size()) / timeBuckets))) {
          printMsg("Building vertex edges", (i / (float)edgeList.size()),
                   t.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);
        }
      }
    }
  } else {
    vector<vector<vector<SimplexId>>> threadedVertexEdges(threadNumber_);
    for(ThreadId i = 0; i < threadNumber_; i++) {
      threadedVertexEdges[i].resize(vertexNumber);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)edgeList.size(); i++) {

      ThreadId threadId = 0;

#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif

      threadedVertexEdges[threadId][edgeList[i].first].push_back(i);
      threadedVertexEdges[threadId][edgeList[i].second].push_back(i);
    }

    // now merge the thing
    for(ThreadId i = 0; i < threadNumber_; i++) {
      for(SimplexId j = 0; j < (SimplexId)threadedVertexEdges[i].size(); j++) {
        for(SimplexId k = 0; k < (SimplexId)threadedVertexEdges[i][j].size();
            k++) {

          bool hasFound = false;
          for(SimplexId l = 0; l < (SimplexId)vertexEdges[j].size(); l++) {
            if(vertexEdges[j][l] == threadedVertexEdges[i][j][k]) {
              hasFound = true;
              break;
            }
          }
          if(!hasFound)
            vertexEdges[j].push_back(threadedVertexEdges[i][j][k]);
        }
      }
    }
  }

  printMsg("Built " + std::to_string(vertexNumber) + " vertex edges", 1,
           t.getElapsedTime(), threadNumber_);

  threadNumber_ = oldThreadNumber;

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
                    for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
                      vertexLink[i * (nbVertCell) + n] = faceIds[n];
                    }
                    hasPivotVertex = true;
                    break;
                  }
                }
              } else if(nbVertCell == 3) {
                // triangle case
                // we're holding to an edge that does not contain our vertex
                for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
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
          for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
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

int ZeroSkeleton::buildVertexLinks(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<vector<LongSimplexId>> &vertexLinks,
  vector<vector<SimplexId>> *vertexStars) const {

  auto localVertexStars = vertexStars;
  vector<vector<SimplexId>> defaultVertexStars{};
  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if((SimplexId)localVertexStars->size() != vertexNumber) {
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
      for(auto cid : (*localVertexStars)[i]) {
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

    const SimplexId nbCellInStarI = (*localVertexStars)[i].size();
    for(SimplexId j = 0; j < nbCellInStarI; j++) {

      const SimplexId cellId = (*localVertexStars)[i][j];
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
    for(SimplexId i = 0; i < (SimplexId)vertexLinks.size(); i++) {
      stringstream msg;
      // TODO: ASSUME Regular Mesh Here!
      msg << "Vertex #" << i << " ("
          << vertexLinks[i].size() / (vertexLinks[i][0] + 1)
          << " items, length: " << vertexLinks[i].size() << "): ";
      printMsg(msg.str(), debug::Priority::DETAIL);
      for(SimplexId j = 0;
          j < (SimplexId)vertexLinks[i].size() / (vertexLinks[i][0] + 1); j++) {
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

int ZeroSkeleton::buildVertexLinks(
  const vector<vector<SimplexId>> &vertexStars,
  const vector<vector<SimplexId>> &cellEdges,
  const vector<pair<SimplexId, SimplexId>> &edgeList,
  vector<vector<SimplexId>> &vertexLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexStars.empty())
    return -1;
  if(cellEdges.empty())
    return -2;
  if(edgeList.empty())
    return -3;
#endif

  Timer t;

  vertexLinks.resize(vertexStars.size());

  printMsg("Building vertex links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)vertexLinks.size(); i++) {

    for(SimplexId j = 0; j < (SimplexId)vertexStars[i].size(); j++) {
      for(SimplexId k = 0; k < (SimplexId)cellEdges[vertexStars[i][j]].size();
          k++) {
        SimplexId edgeId = cellEdges[vertexStars[i][j]][k];

        SimplexId vertexId0 = edgeList[edgeId].first;
        SimplexId vertexId1 = edgeList[edgeId].second;

        if((vertexId0 != i) && (vertexId1 != i)) {
          vertexLinks[i].push_back(edgeId);
        }
      }
    }
  }

  printMsg("Built " + std::to_string(vertexLinks.size()) + " vertex links", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int ZeroSkeleton::buildVertexLinks(
  const vector<vector<SimplexId>> &vertexStars,
  const vector<vector<SimplexId>> &cellTriangles,
  const vector<vector<SimplexId>> &triangleList,
  vector<vector<SimplexId>> &vertexLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexStars.empty())
    return -1;
  if(cellTriangles.empty())
    return -2;
  if(triangleList.empty())
    return -3;
#endif

  Timer t;

  vertexLinks.resize(vertexStars.size());

  printMsg("Building vertex links", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)vertexLinks.size(); i++) {

    for(SimplexId j = 0; j < (SimplexId)vertexStars[i].size(); j++) {
      for(SimplexId k = 0;
          k < (SimplexId)cellTriangles[vertexStars[i][j]].size(); k++) {
        SimplexId triangleId = cellTriangles[vertexStars[i][j]][k];

        bool hasVertex = false;
        for(int l = 0; l < 3; l++) {
          if(i == triangleList[triangleId][l]) {
            hasVertex = true;
            break;
          }
        }

        if(!hasVertex) {
          vertexLinks[i].push_back(triangleId);
        }
      }
    }
  }

  printMsg("Built " + std::to_string(vertexLinks.size()) + " vertex links", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int ZeroSkeleton::buildVertexNeighbors(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &oneSkeleton,
  vector<pair<SimplexId, SimplexId>> *edgeList) const {

  oneSkeleton.resize(vertexNumber);

  auto localEdgeList = edgeList;
  vector<pair<SimplexId, SimplexId>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {
    OneSkeleton osk;
    osk.setDebugLevel(debugLevel_);
    osk.setThreadNumber(threadNumber_);
    osk.buildEdgeList(vertexNumber, cellArray, *localEdgeList);
  }

  Timer t;

  printMsg("Building vertex neighbors", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  const SimplexId nbEdge = localEdgeList->size();
  for(SimplexId i = 0; i < nbEdge; i++) {
    oneSkeleton[(*localEdgeList)[i].first].emplace_back(
      (*localEdgeList)[i].second);
    oneSkeleton[(*localEdgeList)[i].second].emplace_back(
      (*localEdgeList)[i].first);
  }

  printMsg("Built " + std::to_string(vertexNumber) + " vertex neighbors", 1,
           t.getElapsedTime(), threadNumber_);

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // (only merging step, after edge list creation)
  // 1 thread: 9.16 s
  // 24 threads: 13.21 s [not efficient in parallel]

  return 0;
}

int ZeroSkeleton::buildVertexStars(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &vertexStars) const {

  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  ThreadId oldThreadNumber = threadNumber_;
  threadNumber_ = 1;

  Timer t;

  printMsg("Building vertex stars", 0, 0, threadNumber_,
           ttk::debug::LineMode::REPLACE);

  std::vector<std::vector<std::vector<SimplexId>>> zeroSkelThr(threadNumber_);
  for(auto &vi : zeroSkelThr) {
    vi.resize(vertexNumber);
    for(auto &vij : vi) {
      vij.reserve(32);
    }
  }

  if(threadNumber_ > 1) {
    vertexStars.resize(vertexNumber);
    for(auto &vs : vertexStars) {
      vs.reserve(32);
    }
  }

  const SimplexId cellNumber = cellArray.getNbCells();
  const SimplexId timeBuckets = std::min<ttk::SimplexId>(10, cellNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId cid = 0; cid < cellNumber; cid++) {

    ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);
    for(SimplexId j = 0; j < nbVertCell; j++) {
      zeroSkelThr[threadId][cellArray.getCellVertex(cid, j)].emplace_back(cid);
    }

    if(debugLevel_ >= (int)(debug::Priority::INFO) && threadNumber_ == 1) {
      if(!(cid % ((cellNumber) / timeBuckets))) {
        printMsg("Building vertex stars", (cid / (float)cellNumber),
                 t.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);
      }
    }
  }

  if(threadNumber_ > 1) {
    // now merge the thing
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(size_t j = 0; j < vertexStars.size(); j++) {
      auto &dst = vertexStars[j];

      for(size_t i = 0; i < zeroSkelThr.size(); ++i) {
        const auto &vij = zeroSkelThr[i][j];

        // looping on the vertex star
        for(size_t k = 0; k < vij.size(); k++) {
          if(std::find(dst.begin(), dst.end(), vij[k]) == dst.end())
            dst.emplace_back(vij[k]);
        }
      }
    }
  } else {
    vertexStars = std::move(zeroSkelThr[0]);
  }

  printMsg("Built " + std::to_string(vertexNumber) + " vertex stars", 1,
           t.getElapsedTime(), threadNumber_);

  if(debugLevel_ >= static_cast<int>(debug::Priority::DETAIL)) {
    for(SimplexId i = 0; i < (SimplexId)vertexStars.size(); i++) {
      stringstream msg;
      msg << "Vertex #" << i << " (" << vertexStars[i].size() << " cell(s)): ";
      for(SimplexId j = 0; j < (SimplexId)vertexStars[i].size(); j++) {
        msg << " " << vertexStars[i][j];
      }
      printMsg(msg.str(), debug::Priority::INFO);
    }
  }

  threadNumber_ = oldThreadNumber;

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 0.53 s
  // 24 threads: 7.99 s

  return 0;
}
