#include <OneSkeleton.h>

using namespace std;
using namespace ttk;

OneSkeleton::OneSkeleton() {
  setDebugMsgPrefix("OneSkeleton");
}

OneSkeleton::~OneSkeleton() {
}

int OneSkeleton::buildEdgeLinks(
  const vector<pair<SimplexId, SimplexId>> &edgeList,
  const vector<vector<SimplexId>> &edgeStars,
  const CellArray &cellArray,
  vector<vector<SimplexId>> &edgeLinks) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeList.empty())
    return -1;
  if(edgeStars.size() != edgeList.size())
    return -2;
#endif

  Timer t;

  printMsg(
    "Building edge links", 0, 0, threadNumber_, ttk::debug::LineMode::REPLACE);

  edgeLinks.resize(edgeList.size());

  const SimplexId nbEdges = edgeLinks.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < nbEdges; i++) {
    const SimplexId localNbTriangle = edgeStars[i].size();
    for(SimplexId j = 0; j < localNbTriangle; j++) {

      // Look for the right vertex in each triangle
      for(int k = 0; k < 3; k++) {
        const SimplexId tmpVertexId
          = cellArray.getCellVertex(edgeStars[i][j], k);
        if(tmpVertexId != edgeList[i].first
           && tmpVertexId != edgeList[i].second) {
          // found the vertex in the triangle
          edgeLinks[i].push_back(tmpVertexId);
          break;
        }
      }
    }
  }

  printMsg("Built " + to_string(edgeLinks.size()) + " edge links", 1,
           t.getElapsedTime(), threadNumber_);

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

  printMsg(
    "Building edge links", 0, 0, threadNumber_, debug::LineMode::REPLACE);

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

  printMsg("Built " + to_string(edgeLinks.size()) + " edge links", 1,
           t.getElapsedTime(), threadNumber_);

  return 0;
}

int OneSkeleton::buildEdgeList(
  const SimplexId &vertexNumber,
  const CellArray &cellArray,
  vector<pair<SimplexId, SimplexId>> &edgeList) const {

  ThreadId oldThreadNumber = threadNumber_;

  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  // TODO: we should remove parallel part, only dead code.
  threadNumber_ = 1;

  Timer t;

  // TODO: describe the content of this complex container
  // for each thread
  // - a vector of size nb vertices containing
  // - -
  vector<vector<vector<SimplexId>>> threadedEdgeTable(threadNumber_);

  for(auto &vect : threadedEdgeTable) {
    vect.resize(vertexNumber);
  }

  printMsg("Building edges", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  const SimplexId cellNumber = cellArray.getNbCells();
  const int timeBuckets = std::min(10, cellNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId cid = 0; cid < cellNumber; cid++) {

    ThreadId threadId = 0;
    pair<SimplexId, SimplexId> edgeIds;
#ifdef TTK_ENABLE_OPENMP
    threadId = omp_get_thread_num();
#endif

    const SimplexId nbVertsInCell = cellArray.getCellVertexNumber(cid);

    // tet case
    // 0 - 1
    // 0 - 2
    // 0 - 3
    // 1 - 2
    // 1 - 3
    // 2 - 3
    for(SimplexId j = 0; j <= nbVertsInCell - 2; j++) {
      for(SimplexId k = j + 1; k <= nbVertsInCell - 1; k++) {
        // edge processing
        edgeIds.first = cellArray.getCellVertex(cid, j);
        edgeIds.second = cellArray.getCellVertex(cid, k);

        if(edgeIds.first > edgeIds.second) {
          std::swap(edgeIds.first, edgeIds.second);
        }

        bool hasFound = false;
        // TODO Traversing the list each time is suboptimal. Sort + uniq
        for(const auto &l : threadedEdgeTable[threadId][edgeIds.first]) {
          if(edgeIds.second == l) {
            hasFound = true;
            break;
          }
        }
        if(!hasFound) {
          threadedEdgeTable[threadId][edgeIds.first].emplace_back(
            edgeIds.second);
        }
        // end of edge processing
      }
    }
    if(debugLevel_ >= (int)(debug::Priority::INFO)) {
      if(!(cid % ((cellNumber) / timeBuckets)))
        printMsg("Building edges", (cid / (float)cellNumber),
                 t.getElapsedTime(), 1, debug::LineMode::REPLACE);
    }
  }

  // now merge the thing
  SimplexId edgeCount = 0;
  vector<vector<SimplexId>> edgeTable;

  if(threadNumber_ > 1) {
    edgeTable.resize(vertexNumber);
    // All these .size are recomputed each time. Us
    const SimplexId nbThreads = threadedEdgeTable.size();
    for(SimplexId i = 0; i < nbThreads; i++) {

      const LongSimplexId nbVertThreadI = threadedEdgeTable[i].size();
      for(LongSimplexId j = 0; j < nbVertThreadI; j++) {

        const LongSimplexId nbAttachedVertsJ = threadedEdgeTable[i][j].size();
        for(LongSimplexId k = 0; k < nbAttachedVertsJ; k++) {

          // search if it already exists
          bool hasFound = false;

          // TODO here again, linear search is not optimal.
          const LongSimplexId nbVertEdgeJ = edgeTable[j].size();
          for(SimplexId l = 0; l < nbVertEdgeJ; l++) {
            if(edgeTable[j][l] == threadedEdgeTable[i][j][k]) {
              hasFound = true;
              break;
            }
          }
          if(!hasFound) {
            edgeTable[j].emplace_back(threadedEdgeTable[i][j][k]);
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

  printMsg(
    "Built " + to_string(edgeList.size()) + " edges", 1, t.getElapsedTime(), 1);

  threadNumber_ = oldThreadNumber;

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 10.4979 s
  // 24 threads: 12.3994 s [not efficient in parallel]

  return 0;
}

// int OneSkeleton::buildEdgeLists(
//   const vector<vector<LongSimplexId>> &cellArrays,
//   vector<vector<pair<SimplexId, SimplexId>>> &edgeLists) const {
//   Timer t;
//   printMsg(
//     "Building edge lists", 0, 0, threadNumber_, debug::LineMode::REPLACE);
//   edgeLists.resize(cellArrays.size());
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
//   for(SimplexId i = 0; i < (SimplexId)cellArrays.size(); i++) {
//     buildEdgeSubList(cellArrays[i].size() / (cellArrays[i][0] + 1),
//                      cellArrays[i].data(), edgeLists[i]);
//   }
//   printMsg("Built " + to_string(edgeLists.size()) + " edge lists", 1,
//            t.getElapsedTime(), threadNumber_);
//   if(debugLevel_ >= (int)(debug::Priority::DETAIL)) {
//     for(SimplexId i = 0; i < (SimplexId)edgeLists.size(); i++) {
//       {
//         stringstream stringStream;
//         stringStream << "Surface #" << i << " (" << edgeLists[i].size()
//                      << " edges):";
//         printMsg(stringStream.str(), debug::Priority::DETAIL);
//       }
//       for(SimplexId j = 0; j < (SimplexId)edgeLists[i].size(); j++) {
//         stringstream stringStream;
//         stringStream << "- [" << edgeLists[i][j].first << " - "
//                      << edgeLists[i][j].second << "]";
//         printMsg(stringStream.str(), debug::Priority::DETAIL);
//       }
//     }
//   }
//   // computing the edge list of each vertex link:
//   // 24 threads (12 cores): 1.69s.
//   // 1 thread: 7.2 (> x4)
//   return 0;
// }

int OneSkeleton::buildEdgeStars(const SimplexId &vertexNumber,
                                const CellArray &cellArray,
                                vector<vector<SimplexId>> &starList,
                                vector<pair<SimplexId, SimplexId>> *edgeList,
                                vector<vector<SimplexId>> *vertexStars) const {

  auto localEdgeList = edgeList;
  vector<pair<SimplexId, SimplexId>> defaultEdgeList{};
  if(!localEdgeList) {
    localEdgeList = &defaultEdgeList;
  }

  if(!localEdgeList->size()) {
    buildEdgeList(vertexNumber, cellArray, *localEdgeList);
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
      vertexNumber, cellArray, *localVertexStars);
  }

  Timer t;

  printMsg(
    "Building edge stars", 0, 0, threadNumber_, debug::LineMode::REPLACE);

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

  printMsg("Built " + to_string(starList.size()) + " edge stars", 1,
           t.getElapsedTime(), threadNumber_);

  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // with edge list and vertex stars
  // 1 thread: 13 s
  // 24 threads: 48 s (~ x4)

  return 0;
}

int OneSkeleton::buildEdgeSubList(
  const CellArray &cellArray,
  vector<pair<SimplexId, SimplexId>> &edgeList) const {

  // NOTE: here we're dealing with a subportion of the mesh.
  // hence our lookup strategy (based on the number of total vertices) is no
  // longer efficient. let's use a standard map instead
  // NOTE: when dealing with the entire mesh (case above), our vertex based
  // look up strategy is about 7 times faster than the standard map.
  // For mesh portions, the standard map is orders of magnitude faster

  map<pair<SimplexId, SimplexId>, bool> edgeMap;
  edgeList.clear();

  const SimplexId cellNumber = cellArray.getNbCells();
  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    const SimplexId nbVertCell = cellArray.getCellVertexNumber(cid);

    pair<SimplexId, SimplexId> edgeIds;
    int tmpVertexId;
    // tet case
    // 0 - 1
    // 0 - 2
    // 0 - 3
    // 1 - 2
    // 1 - 3
    // 2 - 3
    for(SimplexId j = 0; j <= nbVertCell - 2; j++) {
      for(SimplexId k = j + 1; k <= nbVertCell - 1; k++) {

        edgeIds.first = cellArray.getCellVertex(cid, j);
        edgeIds.second = cellArray.getCellVertex(cid, k);

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
