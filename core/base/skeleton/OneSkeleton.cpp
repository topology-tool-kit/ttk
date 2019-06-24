#include <OneSkeleton.h>

using namespace std;
using namespace ttk;

OneSkeleton::OneSkeleton() {
}

OneSkeleton::~OneSkeleton() {
}

int OneSkeleton::buildEdgeLinks(
  const vector<pair<SimplexId, SimplexId>> &edgeList,
  const vector<vector<SimplexId>> &edgeStars,
  const LongSimplexId *cellArray,
  vector<vector<SimplexId>> &edgeLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeList.empty())
    return -1;
  if((edgeStars.empty()) || (edgeStars.size() != edgeList.size()))
    return -2;
  if(!cellArray)
    return -3;
#endif

  Timer t;

  edgeLinks.resize(edgeList.size());

  SimplexId verticesPerCell = cellArray[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)edgeLinks.size(); i++) {
    for(SimplexId j = 0; j < (SimplexId)edgeStars[i].size(); j++) {

      SimplexId vertexId = -1;
      for(int k = 0; k < 3; k++) {
        if((cellArray[(verticesPerCell + 1) * edgeStars[i][j] + 1 + k]
            != edgeList[i].first)
           && (cellArray[(verticesPerCell + 1) * edgeStars[i][j] + 1 + k]
               != edgeList[i].second)) {
          vertexId = cellArray[(verticesPerCell + 1) * edgeStars[i][j] + 1 + k];
          break;
        }
      }
      if(vertexId != -1) {
        edgeLinks[i].push_back(vertexId);
      }
    }
  }

  {
    stringstream msg;
    msg << "[OneSkeleton] Edge links built in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int OneSkeleton::buildEdgeLinks(
  const vector<pair<SimplexId, SimplexId>> &edgeList,
  const vector<vector<SimplexId>> &edgeStars,
  const vector<vector<SimplexId>> &cellEdges,
  vector<vector<SimplexId>> &edgeLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeList.empty())
    return -1;
  if((edgeStars.empty()) || (edgeStars.size() != edgeList.size()))
    return -2;
  if(cellEdges.empty())
    return -3;
#endif

  Timer t;

  edgeLinks.resize(edgeList.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)edgeLinks.size(); i++) {

    SimplexId otherEdgeId = -1;

    for(SimplexId j = 0; j < (SimplexId)edgeStars[i].size(); j++) {

      SimplexId linkEdgeId = -1;

      for(SimplexId k = 0; k < (SimplexId)cellEdges[edgeStars[i][j]].size();
          k++) {
        otherEdgeId = cellEdges[edgeStars[i][j]][k];

        if((edgeList[otherEdgeId].first != edgeList[i].first)
           && (edgeList[otherEdgeId].first != edgeList[i].second)
           && (edgeList[otherEdgeId].second != edgeList[i].first)
           && (edgeList[otherEdgeId].second != edgeList[i].second)) {
          linkEdgeId = otherEdgeId;
          break;
        }
      }

      edgeLinks[i].push_back(linkEdgeId);
    }
  }

  {
    stringstream msg;
    msg << "[OneSkeleton] Edge links built in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int OneSkeleton::buildEdgeList(
  const SimplexId &vertexNumber,
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<pair<SimplexId, SimplexId>> &edgeList) const {

  ThreadId oldThreadNumber = threadNumber_;

  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  threadNumber_ = 1;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  Timer t;

  vector<vector<vector<SimplexId>>> threadedEdgeTable(threadNumber_);

  for(SimplexId i = 0; i < (SimplexId)threadedEdgeTable.size(); i++) {
    threadedEdgeTable[i].resize(vertexNumber);
  }

  // WARNING!
  // assuming triangulations here
  SimplexId verticesPerCell = cellArray[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < cellNumber; i++) {

    ThreadId threadId = 0, tmpVertexId = 0;
    pair<SimplexId, SimplexId> edgeIds;
#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    // tet case
    // 0 - 1
    // 0 - 2
    // 0 - 3
    // 1 - 2
    // 1 - 3
    // 2 - 3
    for(SimplexId j = 0; j <= verticesPerCell - 2; j++) {
      for(SimplexId k = j + 1; k <= verticesPerCell - 1; k++) {
        // edge processing
        edgeIds.first = cellArray[(verticesPerCell + 1) * i + 1 + j];
        edgeIds.second = cellArray[(verticesPerCell + 1) * i + 1 + k];

        if(edgeIds.first > edgeIds.second) {
          tmpVertexId = edgeIds.first;
          edgeIds.first = edgeIds.second;
          edgeIds.second = tmpVertexId;
        }

        bool hasFound = false;
        for(SimplexId l = 0;
            l < (SimplexId)threadedEdgeTable[threadId][edgeIds.first].size();
            l++) {
          if(edgeIds.second == threadedEdgeTable[threadId][edgeIds.first][l]) {
            hasFound = true;
            break;
          }
        }
        if(!hasFound) {
          threadedEdgeTable[threadId][edgeIds.first].push_back(edgeIds.second);
        }
        // end of edge processing
      }
    }
  }

  // now merge the thing
  SimplexId edgeCount = 0;
  vector<vector<SimplexId>> edgeTable;

  if(threadNumber_ > 1) {
    edgeTable.resize(vertexNumber);
    for(SimplexId i = 0; i < (SimplexId)threadedEdgeTable.size(); i++) {

      for(SimplexId j = 0; j < (SimplexId)threadedEdgeTable[i].size(); j++) {

        for(SimplexId k = 0; k < (SimplexId)threadedEdgeTable[i][j].size();
            k++) {

          // search if it already exists
          bool hasFound = false;

          for(SimplexId l = 0; l < (SimplexId)edgeTable[j].size(); l++) {

            if(edgeTable[j][l] == threadedEdgeTable[i][j][k]) {
              hasFound = true;
              break;
            }
          }
          if(!hasFound) {
            edgeTable[j].push_back(threadedEdgeTable[i][j][k]);
            edgeCount++;
          }
        }
      }
    }
  } else {
    for(SimplexId i = 0; i < (SimplexId)threadedEdgeTable[0].size(); i++)
      edgeCount += threadedEdgeTable[0][i].size();
  }

  vector<vector<SimplexId>> *masterEdgeTable = &edgeTable;
  if(threadNumber_ == 1) {
    masterEdgeTable = &(threadedEdgeTable[0]);
  }

  edgeList.resize(edgeCount);
  edgeCount = 0;
  for(SimplexId i = 0; i < (SimplexId)masterEdgeTable->size(); i++) {

    for(SimplexId j = 0; j < (SimplexId)(*masterEdgeTable)[i].size(); j++) {

      edgeList[edgeCount].first = i;
      edgeList[edgeCount].second = (*masterEdgeTable)[i][j];
      edgeCount++;
    }
  }

  {
    stringstream msg;
    msg << "[OneSkeleton] Edge-list built in " << t.getElapsedTime() << " s. ("
        << edgeList.size() << " edges, " << threadNumber_ << " thread(s))."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  threadNumber_ = oldThreadNumber;

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 10.4979 s
  // 24 threads: 12.3994 s [not efficient in parallel]

  return 0;
}

int OneSkeleton::buildEdgeLists(
  const vector<vector<LongSimplexId>> &cellArrays,
  vector<vector<pair<SimplexId, SimplexId>>> &edgeLists) const {

  Timer t;

  edgeLists.resize(cellArrays.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)cellArrays.size(); i++) {
    buildEdgeSubList(cellArrays[i].size() / (cellArrays[i][0] + 1),
                     cellArrays[i].data(), edgeLists[i]);
  }

  {
    stringstream msg;
    msg << "[OneSkeleton] Multiple edge-lists built in " << t.getElapsedTime()
        << " s. (" << edgeLists.size() << " meshes)" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    stringstream msg;
    for(SimplexId i = 0; i < (SimplexId)edgeLists.size(); i++) {
      msg << "[OneSkeleton] Surface #" << i << " (" << edgeLists[i].size()
          << " edges):" << endl;
      for(SimplexId j = 0; j < (SimplexId)edgeLists[i].size(); j++) {
        msg << "[OneSkeleton] - [" << edgeLists[i][j].first << " - "
            << edgeLists[i][j].second << "]" << endl;
      }
    }
    dMsg(cout, msg.str(), Debug::advancedInfoMsg);
  }

  // computing the edge list of each vertex link:
  // 24 threads (12 cores): 1.69s.
  // 1 thread: 7.2 (> x4)

  return 0;
}

int OneSkeleton::buildEdgeStars(const SimplexId &vertexNumber,
                                const SimplexId &cellNumber,
                                const LongSimplexId *cellArray,
                                vector<vector<SimplexId>> &starList,
                                vector<pair<SimplexId, SimplexId>> *edgeList,
                                vector<vector<SimplexId>> *vertexStars) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellArray)
    return -1;
#endif

  Timer t;

  auto localEdgeList = edgeList;
  vector<pair<SimplexId, SimplexId>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {
    buildEdgeList(vertexNumber, cellNumber, cellArray, *localEdgeList);
  }

  starList.resize(localEdgeList->size());
  for(SimplexId i = 0; i < (SimplexId)starList.size(); i++)
    starList[i].reserve(16);

  auto localVertexStars = vertexStars;
  vector<vector<SimplexId>> defaultVertexStars{};
  if(!localVertexStars) {
    localVertexStars = &defaultVertexStars;
  }
  if((SimplexId)localVertexStars->size() != vertexNumber) {
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.buildVertexStars(
      vertexNumber, cellNumber, cellArray, *localVertexStars);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)localEdgeList->size(); i++) {

    SimplexId vertex0 = (*localEdgeList)[i].first;
    SimplexId vertex1 = (*localEdgeList)[i].second;

    // merge the two vertex stars
    for(SimplexId j = 0; j < (SimplexId)(*localVertexStars)[vertex0].size();
        j++) {

      bool hasFound = false;
      for(SimplexId k = 0; k < (SimplexId)(*localVertexStars)[vertex1].size();
          k++) {
        if((*localVertexStars)[vertex0][j] == (*localVertexStars)[vertex1][k]) {
          hasFound = true;
          break;
        }
      }
      if(hasFound) {
        // common to the two vertex stars
        starList[i].push_back((*localVertexStars)[vertex0][j]);
      }
    }
  }

  {
    stringstream msg;
    msg << "[OneSkeleton] Edge stars built in " << t.getElapsedTime() << " s. ("
        << starList.size() << " edges, " << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // with edge list and vertex stars
  // 1 thread: 13 s
  // 24 threads: 48 s (~ x4)

  return 0;
}

int OneSkeleton::buildEdgeSubList(
  const SimplexId &cellNumber,
  const LongSimplexId *cellArray,
  vector<pair<SimplexId, SimplexId>> &edgeList) const {

  // NOTE: here we're dealing with a subportion of the mesh.
  // hence our lookup strategy (based on the number of total vertices) is no
  // longer efficient. let's use a standard map instead
  // NOTE: when dealing with the entire mesh (case above), our vertex based
  // look up strategy is about 7 times faster than the standard map.
  // For mesh portions, the standard map is orders of magnitude faster

  map<pair<SimplexId, SimplexId>, bool> edgeMap;
  edgeList.clear();

  SimplexId verticesPerCell = cellArray[0];
  for(SimplexId i = 0; i < cellNumber; i++) {

    pair<SimplexId, SimplexId> edgeIds;
    int tmpVertexId;
    // tet case
    // 0 - 1
    // 0 - 2
    // 0 - 3
    // 1 - 2
    // 1 - 3
    // 2 - 3
    for(SimplexId j = 0; j <= verticesPerCell - 2; j++) {
      for(SimplexId k = j + 1; k <= verticesPerCell - 1; k++) {

        edgeIds.first = cellArray[(verticesPerCell + 1) * i + 1 + j];
        edgeIds.second = cellArray[(verticesPerCell + 1) * i + 1 + k];

        if(edgeIds.first > edgeIds.second) {
          tmpVertexId = edgeIds.first;
          edgeIds.first = edgeIds.second;
          edgeIds.second = tmpVertexId;
        }

        map<pair<SimplexId, SimplexId>, bool>::iterator it
          = edgeMap.find(edgeIds);

        if(it == edgeMap.end()) {
          // not found, let's add this edge
          edgeList.push_back(edgeIds);
          edgeMap[edgeIds] = true;
        }
      }
    }
  }

  return 0;
}
