#include <ThreeSkeleton.h>

using namespace std;
using namespace ttk;

ThreeSkeleton::ThreeSkeleton() {
}

ThreeSkeleton::~ThreeSkeleton() {
}

int ThreeSkeleton::buildCellEdges(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &cellEdges,
  vector<pair<SimplexId, SimplexId>> *edgeList,
  vector<vector<SimplexId>> *vertexEdges) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexNumber <= 0)
    return -1;
  if(cellNumber <= 0)
    return -2;
  if(!cellArray)
    return -3;
#endif

  Timer t;

  auto localEdgeList = edgeList;
  auto localVertexEdges = vertexEdges;
  vector<pair<SimplexId, SimplexId>> defaultEdgeList{};
  vector<vector<SimplexId>> defaultVertexEdges{};

  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {

    OneSkeleton oneSkeleton;
    oneSkeleton.setDebugLevel(debugLevel_);
    oneSkeleton.setThreadNumber(threadNumber_);
    oneSkeleton.buildEdgeList(
      vertexNumber, cellNumber, cellArray, *localEdgeList);
  }

  if(!localVertexEdges) {
    localVertexEdges = &defaultVertexEdges;
  }

  if(!localVertexEdges->size()) {

    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.buildVertexEdges(
      vertexNumber, *localEdgeList, *localVertexEdges);
  }

  cellEdges.resize(cellNumber);
  for(SimplexId i = 0; i < (SimplexId)cellEdges.size(); i++) {
    // optimized for tet meshes
    cellEdges[i].reserve(6);
  }

  int vertexPerCell = cellArray[0];

  // for each cell, for each pair of vertices, find the edge
  // TODO: check for parallel efficiency here
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber; i++) {

    SimplexId cellId = (vertexPerCell + 1) * i;

    for(SimplexId j = 0; j < vertexPerCell; j++) {

      for(SimplexId k = j + 1; k < vertexPerCell; k++) {

        SimplexId vertexId0 = cellArray[cellId + 1 + j];
        SimplexId vertexId1 = cellArray[cellId + 1 + k];

        // loop around the edges of vertexId0 in search of vertexId1
        SimplexId edgeId = -1;
        for(SimplexId l = 0;
            l < (SimplexId)(*localVertexEdges)[vertexId0].size(); l++) {

          SimplexId localEdgeId = (*localVertexEdges)[vertexId0][l];
          if(((*localEdgeList)[localEdgeId].first == vertexId1)
             || ((*localEdgeList)[localEdgeId].second == vertexId1)) {
            edgeId = localEdgeId;
            break;
          }
        }

        cellEdges[i].push_back(edgeId);
      }
    }
  }

  {
    stringstream msg;
    msg << "[ThreeSkeleton] Cell edges built in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ThreeSkeleton::buildCellNeighborsFromTriangles(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &cellNeighbors,
  vector<vector<SimplexId>> *triangleStars) const {

  Timer t;

  auto localTriangleStars = triangleStars;
  vector<vector<SimplexId>> defaultTriangleStars{};
  if(!localTriangleStars) {
    localTriangleStars = &defaultTriangleStars;
  }

  if(!localTriangleStars->size()) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setThreadNumber(threadNumber_);
    twoSkeleton.setDebugLevel(debugLevel_);
    twoSkeleton.buildTriangleList(
      vertexNumber, cellNumber, cellArray, NULL, localTriangleStars);
  }

  SimplexId vertexPerCell = cellArray[0];

  cellNeighbors.resize(cellNumber);
  for(SimplexId i = 0; i < (SimplexId)cellNeighbors.size(); i++) {
    cellNeighbors[i].reserve(vertexPerCell);
  }

  // NOTE: not efficient so far in parallel
  ThreadId oldThreadNumber = threadNumber_;
  threadNumber_ = 1;

  if(threadNumber_ == 1) {

    for(SimplexId i = 0; i < (SimplexId)localTriangleStars->size(); i++) {

      if((*localTriangleStars)[i].size() == 2) {

        // interior triangle
        cellNeighbors[(*localTriangleStars)[i][0]].push_back(
          (*localTriangleStars)[i][1]);

        cellNeighbors[(*localTriangleStars)[i][1]].push_back(
          (*localTriangleStars)[i][0]);
      }
    }
  } else {

    vector<vector<vector<SimplexId>>> threadedCellNeighbors(threadNumber_);

    for(ThreadId i = 0; i < threadNumber_; i++) {
      threadedCellNeighbors[i].resize(cellNumber);
      for(SimplexId j = 0; j < (SimplexId)threadedCellNeighbors[i].size(); j++)
        threadedCellNeighbors[i][j].reserve(4);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)(*localTriangleStars).size(); i++) {

      ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif

      if((*localTriangleStars)[i].size() == 2) {

        // interior triangles

        threadedCellNeighbors[threadId][(*localTriangleStars)[i][0]].push_back(
          (*localTriangleStars)[i][1]);

        threadedCellNeighbors[threadId][(*localTriangleStars)[i][1]].push_back(
          (*localTriangleStars)[i][0]);
      }
    }

    // now merge things
    for(ThreadId i = 0; i < threadNumber_; i++) {
      for(SimplexId j = 0; j < (SimplexId)threadedCellNeighbors[i].size();
          j++) {

        for(SimplexId k = 0; k < (SimplexId)threadedCellNeighbors[i][j].size();
            k++) {
          cellNeighbors[j].push_back(threadedCellNeighbors[i][j][k]);
        }
      }
    }
  }

  {
    stringstream msg;
    msg << "[ThreeSkeleton] Cell neighbors (" << cellNumber
        << " cells) computed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  threadNumber_ = oldThreadNumber;

  // ethaneDiol.vtu, 8.7Mtets, vger (4coresHT)
  // 1 thread: 9.80 s
  // 4 threads: 14.18157 s

  // ethaneDiol.vtu, 8.7Mtets, richard (8coresHT)
  // 1 thread: 9.6763 s
  // 4 threads: 16.7489 s

  // ethaneDiol.vtu, 8.7Mtets, vger (4coresHT), VTK implementation
  // 1 thread: 14.9023 s
  // 4 threads: 8.58048 s

  // ethaneDiol.vtu, 8.7Mtets, richard (8coresHT), VTK implementation
  // 1 thread: 15.1387  s
  // 8 threads: 5.38286 s

  return 0;
}

int ThreeSkeleton::buildCellNeighborsFromVertices(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &cellNeighbors,
  vector<vector<SimplexId>> *vertexStars) const {

  if(cellArray[0] == 3) {

    TwoSkeleton twoSkeleton;
    twoSkeleton.setDebugLevel(debugLevel_);
    twoSkeleton.setThreadNumber(threadNumber_);
    return twoSkeleton.buildCellNeighborsFromVertices(
      vertexNumber, cellNumber, cellArray, cellNeighbors, vertexStars);
  }

  if(cellArray[0] == 2) {
    // 1D
    stringstream msg;
    msg << "[ThreeSkeleton] buildCellNeighborsFromVertices in 1D:" << endl;
    msg << "[ThreeSkeleton] Not implemented! TODO!" << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
    return -1;
  }

  Timer t;

  auto localVertexStars = vertexStars;
  vector<vector<SimplexId>> defaultVertexStars{};

  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }

  if(!localVertexStars->size()) {

    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.buildVertexStars(
      vertexNumber, cellNumber, cellArray, *localVertexStars);
  }

  int vertexPerCell = cellArray[0];

  cellNeighbors.resize(cellNumber);
  for(SimplexId i = 0; i < (SimplexId)cellNeighbors.size(); i++)
    cellNeighbors[i].reserve(vertexPerCell);

    // pre-sort vertex stars
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber; i++)
    sort((*localVertexStars)[i].begin(), (*localVertexStars)[i].end());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber; i++) {

    // go triangle by triangle
    for(SimplexId j = 0; j < vertexPerCell; j++) {

      SimplexId v0
        = cellArray[(vertexPerCell + 1) * i + 1 + (j) % vertexPerCell];
      SimplexId v1
        = cellArray[(vertexPerCell + 1) * i + 1 + (j + 1) % vertexPerCell];
      SimplexId v2
        = cellArray[(vertexPerCell + 1) * i + 1 + (j + 2) % vertexPerCell];

      // perform an intersection of the 3 (sorted) star lists
      SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
      SimplexId intersection = -1;

      while((pos0 < (SimplexId)(*localVertexStars)[v0].size())
            && (pos1 < (SimplexId)(*localVertexStars)[v1].size())
            && (pos2 < (SimplexId)(*localVertexStars)[v2].size())) {

        SimplexId biggest = (*localVertexStars)[v0][pos0];
        if((*localVertexStars)[v1][pos1] > biggest) {
          biggest = (*localVertexStars)[v1][pos1];
        }
        if((*localVertexStars)[v2][pos2] > biggest) {
          biggest = (*localVertexStars)[v2][pos2];
        }

        for(SimplexId l = pos0; l < (SimplexId)(*localVertexStars)[v0].size();
            l++) {
          if((*localVertexStars)[v0][l] < biggest) {
            pos0++;
          } else {
            break;
          }
        }
        for(SimplexId l = pos1; l < (SimplexId)(*localVertexStars)[v1].size();
            l++) {
          if((*localVertexStars)[v1][l] < biggest) {
            pos1++;
          } else {
            break;
          }
        }
        for(SimplexId l = pos2; l < (SimplexId)(*localVertexStars)[v2].size();
            l++) {
          if((*localVertexStars)[v2][l] < biggest) {
            pos2++;
          } else {
            break;
          }
        }

        if((pos0 < (SimplexId)(*localVertexStars)[v0].size())
           && (pos1 < (SimplexId)(*localVertexStars)[v1].size())
           && (pos2 < (SimplexId)(*localVertexStars)[v2].size())) {

          if(((*localVertexStars)[v0][pos0] == (*localVertexStars)[v1][pos1])
             && ((*localVertexStars)[v0][pos0]
                 == (*localVertexStars)[v2][pos2])) {

            if((*localVertexStars)[v0][pos0] != i) {
              intersection = (*localVertexStars)[v0][pos0];
              break;
            }

            pos0++;
            pos1++;
            pos2++;
          }
        }
      }

      if(intersection != -1) {
        cellNeighbors[i].push_back(intersection);
      }
    }
  }

  {
    stringstream msg;
    msg << "[ThreeSkeleton] Cell neighbors (" << cellNumber
        << " cells) computed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  // ethaneDiol.vtu, 8.7Mtets, richard (4coresHT)
  // 1 thread: 9.39488 s
  // 2 threads: 5.93132 s [faster than any other implementation]
  // 4 threads: 3.3445 s [~x3, not too bad]

  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 4.4475 s [speedup on VTK ~2.5]
  // 24 threads: 1.4488 s

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 29.4951 s [speedup on VTK ~3.3]
  // 24 threads: 12.7058 s (parallel speed up ~3, speed up on VTK: 2)

  // NOTE:
  // when dealing with simple types such as int, it's better to have heap
  // variable (instead of a vector<int> indexed by threadId).
  // it is MUCH MUCH MUCH better
  // especially with write access (cache misses)

  return 0;
}
