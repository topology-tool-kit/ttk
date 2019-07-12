#include <OneSkeleton.h>
#include <ZeroSkeleton.h>

using namespace std;
using namespace ttk;

ZeroSkeleton::ZeroSkeleton() {
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
    for(SimplexId i = 0; i < (SimplexId)edgeList.size(); i++) {
      vertexEdges[edgeList[i].first].push_back(i);
      vertexEdges[edgeList[i].second].push_back(i);
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

  {
    stringstream msg;
    msg << "[ZeroSkeleton] Vertex edges built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  threadNumber_ = oldThreadNumber;

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 11.85 s
  // 24 threads: 20.93 s [not efficient]

  return 0;
}

int ZeroSkeleton::buildVertexLink(const SimplexId &vertexId,
                                  const SimplexId &cellNumber,
                                  const LongSimplexId *cellArray,
                                  vector<LongSimplexId> &vertexLink) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray[0];

  vector<SimplexId> vertexStar;
  for(SimplexId i = 0; i < cellNumber; i++) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      if(cellArray[(verticesPerCell + 1) * i + j] == vertexId) {
        vertexStar.push_back(i);
        break;
      }
    }
  }

  vertexLink.resize(vertexStar.size() * verticesPerCell);
  vector<SimplexId> faceIds(verticesPerCell);
  faceIds[0] = verticesPerCell - 1;

  // construct the link
  for(SimplexId i = 0; i < (SimplexId)vertexStar.size(); i++) {
    SimplexId cellId = vertexStar[i];

    bool hasPivotVertex = false;

    // iterate on the cell's faces
    for(int k = 0; k < 2; k++) {
      faceIds[1] = cellArray[(verticesPerCell + 1) * cellId + 1 + k];

      if(faceIds[1] != vertexId) {

        if(verticesPerCell > 2) {

          for(SimplexId l = k + 1; l <= verticesPerCell - 1; l++) {
            faceIds[2] = cellArray[(verticesPerCell + 1) * cellId + 1 + l];

            if(faceIds[2] != vertexId) {

              if(verticesPerCell == 4) {
                // tet case, faceIds has 4 entries to fill
                for(SimplexId m = l + 1; m < verticesPerCell; m++) {
                  faceIds[3] = cellArray[(verticesPerCell + 1) * cellId + m];

                  if(faceIds[3] != vertexId) {

                    // all the vertices of the face are different from our
                    // vertex. let's add that face to the link
                    for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
                      vertexLink[i * (verticesPerCell) + n] = faceIds[n];
                    }
                    hasPivotVertex = true;
                    break;
                  }
                }
              } else if(verticesPerCell == 3) {
                // triangle case
                // we're holding to an edge that does not contain our vertex
                for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
                  vertexLink[i * (verticesPerCell) + n] = faceIds[n];
                }
                hasPivotVertex = true;
                break;
              }
              if(hasPivotVertex)
                break;
            }
          }
        } else if(verticesPerCell == 2) {
          // edge-mesh case
          // we're holding a neighbor different from our vertex
          for(SimplexId n = 0; n < (SimplexId)faceIds.size(); n++) {
            vertexLink[i * (verticesPerCell) + n] = faceIds[n];
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
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<LongSimplexId>> &vertexLinks,
  vector<vector<SimplexId>> *vertexStars) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  Timer t;

  auto localVertexStars = vertexStars;
  vector<vector<SimplexId>> defaultVertexStars{};
  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if((SimplexId)localVertexStars->size() != vertexNumber) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.buildVertexStars(
      vertexNumber, cellNumber, cellArray, *localVertexStars);
  }

  // WARNING
  // assuming triangulations
  int verticesPerCell = cellArray[0];

  if((SimplexId)vertexLinks.size() != vertexNumber) {
    vertexLinks.resize(vertexNumber);
    for(SimplexId i = 0; i < (SimplexId)vertexLinks.size(); i++) {
      vertexLinks[i].resize((*localVertexStars)[i].size() * verticesPerCell);
    }
  }

  vector<vector<SimplexId>> faceIds(threadNumber_);
  for(ThreadId i = 0; i < threadNumber_; i++) {
    faceIds[i].resize(verticesPerCell);
    // the first entry should be the number of vertices
    // tet-mesh (4) => link made of triangles (3)
    // triangle-mesh (3) => link made of edges (2)
    faceIds[i][0] = verticesPerCell - 1;
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

    for(SimplexId j = 0; j < (SimplexId)(*localVertexStars)[i].size(); j++) {

      SimplexId cellId = (*localVertexStars)[i][j];

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

        faceIds[threadId][1]
          = cellArray[(verticesPerCell + 1) * cellId + 1 + k];

        if(faceIds[threadId][1] != i) {

          if(verticesPerCell > 2) {

            for(SimplexId l = k + 1; l <= verticesPerCell - 1; l++) {
              faceIds[threadId][2]
                = cellArray[(verticesPerCell + 1) * cellId + 1 + l];

              if(faceIds[threadId][2] != i) {

                if(verticesPerCell == 4) {
                  // tet case, faceIds[threadId] has 4 entries to fill
                  for(SimplexId m = l + 1; m < verticesPerCell; m++) {
                    faceIds[threadId][3]
                      = cellArray[(verticesPerCell + 1) * cellId + 1 + m];

                    // now test if this face contains our vertex or not
                    // there's should be only one face
                    if(faceIds[threadId][3] != i) {
                      // all the vertices of the face are different from our
                      // pivot vertex.
                      // let's add that face to the link
                      for(SimplexId n = 0;
                          n < (SimplexId)faceIds[threadId].size(); n++) {
                        vertexLinks[i][j * (verticesPerCell) + n]
                          = faceIds[threadId][n];
                      }
                      hasPivotVertex = true;
                      break;
                    }
                  }
                } else if(verticesPerCell == 3) {
                  // triangle case
                  // we're holding to an edge that does not contain our vertex
                  for(SimplexId n = 0; n < (SimplexId)faceIds[threadId].size();
                      n++) {
                    vertexLinks[i][j * (verticesPerCell) + n]
                      = faceIds[threadId][n];
                  }
                  hasPivotVertex = true;
                  break;
                }
                if(hasPivotVertex)
                  break;
              }
            }
          } else if(verticesPerCell == 2) {
            // edge-mesh case
            // we holding a neighbor different from us
            for(SimplexId n = 0; n < (SimplexId)faceIds[threadId].size(); n++) {
              vertexLinks[i][j * (verticesPerCell) + n] = faceIds[threadId][n];
            }
            hasPivotVertex = true;
            break;
          }
        }
        if(hasPivotVertex)
          break;
      }
    }
  }

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    stringstream msg;
    for(SimplexId i = 0; i < (SimplexId)vertexLinks.size(); i++) {
      msg << "[ZeroSkeleton] Vertex #" << i << " ("
          << vertexLinks[i].size() / (vertexLinks[i][0] + 1)
          << " items, length: " << vertexLinks[i].size() << "): " << endl;
      for(SimplexId j = 0;
          j < (SimplexId)vertexLinks[i].size() / (vertexLinks[i][0] + 1); j++) {
        msg << "[ZeroSkeleton]   - " << j << ":";
        for(SimplexId k = 0; k < verticesPerCell; k++) {
          msg << " " << vertexLinks[i][j * (verticesPerCell) + k];
        }
        msg << endl;
      }
    }
    dMsg(cout, msg.str(), Debug::advancedInfoMsg);
  }

  {
    stringstream msg;
    msg << "[ZeroSkeleton] Vertex links built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

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

  {
    stringstream msg;
    msg << "[ZeroSkeleton] Vertex links built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

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

  {
    stringstream msg;
    msg << "[ZeroSkeleton] Vertex links built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ZeroSkeleton::buildVertexNeighbors(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &oneSkeleton,
  vector<pair<SimplexId, SimplexId>> *edgeList) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  Timer t;

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
    osk.buildEdgeList(vertexNumber, cellNumber, cellArray, *localEdgeList);
  }

  for(SimplexId i = 0; i < (SimplexId)localEdgeList->size(); i++) {
    oneSkeleton[(*localEdgeList)[i].first].push_back(
      (*localEdgeList)[i].second);
    oneSkeleton[(*localEdgeList)[i].second].push_back(
      (*localEdgeList)[i].first);
  }

  {
    stringstream msg;
    msg << "[ZeroSkeleton] One-skeleton built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // (only merging step, after edge list creation)
  // 1 thread: 9.16 s
  // 24 threads: 13.21 s [not efficient in parallel]

  return 0;
}

int ZeroSkeleton::buildVertexStars(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &vertexStars) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  ThreadId oldThreadNumber = threadNumber_;
  threadNumber_ = 1;

  Timer t;

  vertexStars.resize(vertexNumber);
  for(SimplexId i = 0; i < vertexNumber; i++)
    vertexStars[i].reserve(32);

  vector<vector<vector<SimplexId>> *> threadedZeroSkeleton(threadNumber_);
  if(threadNumber_ == 1) {
    threadedZeroSkeleton[0] = &vertexStars;
  } else {
    for(ThreadId i = 0; i < threadNumber_; i++) {
      threadedZeroSkeleton[i] = new vector<vector<SimplexId>>();
      threadedZeroSkeleton[i]->resize(vertexNumber);
      for(SimplexId j = 0; j < (SimplexId)threadedZeroSkeleton[i]->size();
          j++) {
        (*threadedZeroSkeleton[i])[j].reserve(32);
      }
    }
  }

  SimplexId vertexNumberPerCell = cellArray[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber; i++) {

    ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    for(SimplexId j = 0; j < vertexNumberPerCell; j++) {
      (*threadedZeroSkeleton[threadId])[cellArray[(vertexNumberPerCell + 1) * i
                                                  + 1 + j]]
        .push_back(i);
    }
  }

  if(threadNumber_ > 1) {
    // now merge the thing
    for(ThreadId i = 0; i < threadNumber_; i++) {

      // looping on vertices
      for(SimplexId j = 0; j < (SimplexId)threadedZeroSkeleton[i]->size();
          j++) {

        // looping on the vertex star
        for(SimplexId k = 0;
            k < (SimplexId)(*threadedZeroSkeleton[i])[j].size(); k++) {

          bool hasFound = false;
          if(threadNumber_) {
            for(SimplexId l = 0; l < (SimplexId)vertexStars[j].size(); l++) {
              if(vertexStars[j][l] == (*threadedZeroSkeleton[i])[j][k]) {
                hasFound = true;
                break;
              }
            }
          }
          if(!hasFound)
            vertexStars[j].push_back((*threadedZeroSkeleton[i])[j][k]);
        }
      }
    }

    for(ThreadId i = 0; i < threadNumber_; i++) {
      delete threadedZeroSkeleton[i];
    }
  }

  {
    stringstream msg;
    msg << "[ZeroSkeleton] Vertex stars built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    stringstream msg;

    for(SimplexId i = 0; i < (SimplexId)vertexStars.size(); i++) {
      msg << "[ZeroSkeleton] Vertex #" << i << " (" << vertexStars[i].size()
          << " cell(s)): ";
      for(SimplexId j = 0; j < (SimplexId)vertexStars[i].size(); j++) {
        msg << " " << vertexStars[i][j];
      }
      msg << endl;
    }
    dMsg(cout, msg.str(), Debug::advancedInfoMsg);
  }

  threadNumber_ = oldThreadNumber;

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 0.53 s
  // 24 threads: 7.99 s

  return 0;
}
