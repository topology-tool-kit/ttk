#include <ContourTree.h>

#include <cmath>
#include <iomanip>
#include <queue>
#include <set>

using namespace std;
using namespace ttk;

struct MyCmp {
  const vector<double> *vertexScalars_;
  const vector<int> *vertexOffsets_;
  const vector<Node> *nodeList_;
  bool isAscendingOrder_;

  MyCmp(const vector<double> *vertexScalars,
        const vector<int> *vertexOffsets,
        const vector<Node> *nodeList,
        bool isAscendingOrder)
    : vertexScalars_(vertexScalars), vertexOffsets_(vertexOffsets),
      nodeList_(nodeList), isAscendingOrder_(isAscendingOrder) {
  }

  bool operator()(int node1, int node2) {
    int vertex1 = (*nodeList_)[node1].getVertexId();
    int vertex2 = (*nodeList_)[node2].getVertexId();
    bool cmp = ((*vertexScalars_)[vertex1] < (*vertexScalars_)[vertex2])
               || ((*vertexScalars_)[vertex1] == (*vertexScalars_)[vertex2]
                   && (*vertexOffsets_)[vertex1] < (*vertexOffsets_)[vertex2]);
    if(isAscendingOrder_)
      return cmp;
    else
      return !cmp;
  }
};

struct filtrationCtCmp {
  bool operator()(const pair<bool, pair<double, pair<int, int>>> &v0,
                  const pair<bool, pair<double, pair<int, int>>> &v1) const {

    if(v0.first) {
      return ((v0.second.first < v1.second.first)
              || ((v0.second.first == v1.second.first)
                  && (v0.second.second.first < v1.second.second.first)));
    } else {
      return ((v0.second.first > v1.second.first)
              || ((v0.second.first == v1.second.first)
                  && (v0.second.second.first > v1.second.second.first)));
    }
  }
} filtrationCmp;

struct _persistenceCmp {
  bool operator()(const pair<pair<int, int>, double> &p0,
                  const pair<pair<int, int>, double> &p1) {

    return p0.second < p1.second;
  }

} _pCmp;

struct _persistenceCmp3 {
  _persistenceCmp3(const vector<double> *vertexScalars)
    : vertexScalars_{vertexScalars} {
  }

  bool operator()(const pair<pair<int, int>, double> &p0,
                  const pair<pair<int, int>, double> &p1) {
    return (*vertexScalars_)[p0.first.first]
           > (*vertexScalars_)[p1.first.first];
  }

  const vector<double> *vertexScalars_;
};

struct _persistencePairCmp {
  bool operator()(const pair<double, double> &p0,
                  const pair<double, double> &p1) {
    return p0.first < p1.first;
  }
} _pPairCmp;

void SuperArc::smooth(const vector<Node> &nodeList,
                      const vector<vector<double>> *vertexPositions,
                      bool order) {
  int N = getNumberOfBarycenters();
  if(N) {
    /// init ///
    vector<vector<double>> barycenterList(barycenterList_.size());
    for(unsigned int i = 0; i < barycenterList.size(); ++i)
      barycenterList[i].resize(3);

    int up_vId;
    int down_vId;
    if(order) {
      up_vId = nodeList[getUpNodeId()].getVertexId();
      down_vId = nodeList[getDownNodeId()].getVertexId();
    } else {
      up_vId = nodeList[getDownNodeId()].getVertexId();
      down_vId = nodeList[getUpNodeId()].getVertexId();
    }

    double p0[3];
    double p1[3];
    for(unsigned int k = 0; k < 3; ++k) {
      p0[k] = (*vertexPositions)[down_vId][k];
      p1[k] = (*vertexPositions)[up_vId][k];
    }

    /// filtering ///
    if(N > 1) {
      // first
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k] = (p0[k] + barycenterList_[1][k]) * 0.5;

      // main
      for(int i = 1; i < N - 1; ++i) {
        for(unsigned int k = 0; k < 3; ++k)
          barycenterList[i][k]
            = (barycenterList_[i - 1][k] + barycenterList_[i + 1][k]) * 0.5;
      }
      // last
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[N - 1][k] = (p1[k] + barycenterList_[N - 1][k]) * 0.5;
    } else {
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k] = (p0[k] + p1[k]) * 0.5;
    }

    /// copy ///
    for(int i = 0; i < N; ++i) {
      for(unsigned int k = 0; k < 3; ++k)
        barycenterList_[i][k] = barycenterList[i][k];
    }
  }
}

void SuperArc::sortRegularNodes(const vector<double> *vertexScalars,
                                const vector<int> *vertexOffsets,
                                const vector<Node> *nodeList,
                                bool order) {
  MyCmp cmp(vertexScalars, vertexOffsets, nodeList, order);

  if(order)
    std::sort(regularNodeList_.begin(), regularNodeList_.end(), cmp);
  else
    std::sort(regularNodeList_.rbegin(), regularNodeList_.rend(), cmp);
}

double PersistenceMetric::computeSuperArcMetric(
  const int &downVertexId,
  const int &upVertexId,
  const vector<int> &ttkNotUsed(interiorNodeIds)) {

  if(!tree_)
    return -DBL_MAX;

  double downScalar = 0, upScalar = 0;

  tree_->getVertexScalar(downVertexId, downScalar);
  tree_->getVertexScalar(upVertexId, upScalar);

  return fabs(upScalar - downScalar);
}

SubLevelSetTree::SubLevelSetTree() {
  this->setDebugMsgPrefix("SubLevelSetTree");
}

int SubLevelSetTree::appendRegularNode(const int &superArcId,
                                       const int &nodeId) {

  superArcList_[superArcId].appendRegularNode(nodeId);
  vertex2superArc_[nodeList_[nodeId].getVertexId()] = superArcId;
  vertex2superArcNode_[nodeList_[nodeId].getVertexId()]
    = superArcList_[superArcId].getNumberOfRegularNodes() - 1;

  // create a low level arc with the previous node
  int nId = -1;
  if(superArcList_[superArcId].getNumberOfRegularNodes() == 1) {
    // the down neighbor is the arc extremity
    nId = superArcList_[superArcId].getDownNodeId();
  } else {
    nId = superArcList_[superArcId].getRegularNodeId(
      superArcList_[superArcId].getNumberOfRegularNodes() - 2);
  }

  makeArc(nId, nodeId);

  return 0;
}

int SubLevelSetTree::build() {

  Timer timer;

  if(!vertexNumber_)
    return -1;
  if((!vertexScalars_) || ((int)vertexScalars_->size() != vertexNumber_))
    return -2;
  if((!vertexSoSoffsets_) || ((int)vertexSoSoffsets_->size() != vertexNumber_))
    return -3;
  // TODO uncomment after triangulation templatization
  // if(triangulation_->isEmpty())
  //   return -4;
  if((!minimumList_) && (!maximumList_))
    return -5;
  if((minimumList_) && (maximumList_))
    return -6;

  vector<UnionFind> seeds;
  vector<vector<int>> seedSuperArcs;
  vector<UnionFind *> vertexSeeds(vertexNumber_, (UnionFind *)NULL);
  vector<UnionFind *> starSets;
  vector<bool> visitedVertices(vertexNumber_, false);

  SimplexId vertexId = -1, nId = -1;
  UnionFind *seed = NULL, *firstUf = NULL;

  const vector<int> *extremumList = NULL;

  bool isMergeTree = true;

  if(minimumList_) {
    extremumList = minimumList_;
  } else {
    isMergeTree = false;
    extremumList = maximumList_;
  }

  set<
    // is merging
    pair<bool,
         // function value
         pair<double,
              // vertex offset
              pair<int,
                   // vertex id
                   int>>>,
    filtrationCtCmp>
    filtrationFront;

  seeds.resize(extremumList->size());
  seedSuperArcs.resize(seeds.size());

  for(int i = 0; i < (int)extremumList->size(); i++) {
    // link each minimum to a union find seed
    vertexSeeds[(*extremumList)[i]] = &(seeds[i]);

    // open an arc
    seedSuperArcs[i].push_back(openSuperArc(makeNode((*extremumList)[i])));

    // add each minimum to the filtration front
    pair<bool, pair<double, pair<int, int>>> v;
    v.first = isMergeTree;
    v.second.first = (*vertexScalars_)[(*extremumList)[i]];
    v.second.second.first = (*vertexSoSoffsets_)[(*extremumList)[i]];
    v.second.second.second = (*extremumList)[i];

    filtrationFront.insert(v);
  }

  bool merge = false;

  // filtration loop
  do {

    vertexId = filtrationFront.begin()->second.second.second;
    filtrationFront.erase(filtrationFront.begin());

    starSets.clear();

    merge = false;
    firstUf = NULL;

    SimplexId neighborNumber
      = triangulation_->getVertexNeighborNumber(vertexId);
    for(SimplexId i = 0; i < neighborNumber; i++) {
      triangulation_->getVertexNeighbor(vertexId, i, nId);

      if(vertexSeeds[nId]) {
        seed = vertexSeeds[nId]->find();
        starSets.push_back(seed);

        // is it merging things?
        if(!firstUf)
          firstUf = seed;
        else if(seed != firstUf)
          merge = true;
      }

      if(!visitedVertices[nId]) {

        pair<bool, pair<double, pair<int, int>>> v;
        v.first = isMergeTree;
        v.second.first = (*vertexScalars_)[nId];
        v.second.second.first = (*vertexSoSoffsets_)[nId];
        v.second.second.second = nId;

        filtrationFront.insert(v);
        visitedVertices[nId] = true;
      }
    }

    if(!vertexSeeds[vertexId]) {

      vertexSeeds[vertexId] = UnionFind::makeUnion(starSets);

      int newNodeId = makeNode(vertexId);

      if(merge) {

        vector<int> seedIds;
        for(int i = 0; i < (int)starSets.size(); i++) {
          int seedId = starSets[i] - &(seeds[0]);
          bool found = false;
          for(int j = 0; j < (int)seedIds.size(); j++) {
            if(seedIds[j] == seedId) {
              found = true;
              break;
            }
          }
          if(!found)
            seedIds.push_back(seedId);
        }

        for(int i = 0; i < (int)seedIds.size(); i++) {
          int superArcId
            = seedSuperArcs[seedIds[i]][seedSuperArcs[seedIds[i]].size() - 1];
          closeSuperArc(superArcId, newNodeId);
        }

        int seedId = vertexSeeds[vertexId] - &(seeds[0]);
        if(!filtrationFront.empty())
          seedSuperArcs[seedId].push_back(openSuperArc(newNodeId));
      } else if(starSets.size()) {
        // we're dealing with a degree-2 node
        int seedId = starSets[0] - &(seeds[0]);
        int superArcId
          = seedSuperArcs[seedId][seedSuperArcs[seedId].size() - 1];

        if(filtrationFront.empty()) {
          // last vertex
          closeSuperArc(superArcId, newNodeId);
        } else {
          if(maintainRegularVertices_)
            appendRegularNode(superArcId, newNodeId);
        }
      }
    }

    visitedVertices[vertexId] = true;

  } while(!filtrationFront.empty());

  const std::string tree = isMergeTree ? "JoinTree" : "SplitTree";
  this->printMsg(
    tree + " computed", 1.0, timer.getElapsedTime(), this->threadNumber_);
  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Nodes", std::to_string(getNumberOfNodes())},
    {"#Arcs", std::to_string(getNumberOfArcs())}});

  print();

  return 0;
}

int SubLevelSetTree::clearArc(const int &vertexId0, const int &vertexId1) {

  if((vertexId0 < 0) || (vertexId0 >= vertexNumber_))
    return -1;
  if((vertexId1 < 0) || (vertexId1 >= vertexNumber_))
    return -2;

  int nodeId0 = vertex2node_[vertexId0], nodeId1 = vertex2node_[vertexId1];

  for(int i = 0; i < nodeList_[nodeId0].getNumberOfUpArcs(); i++) {
    if(nodeList_[arcList_[nodeList_[nodeId0].getUpArcId(i)].getUpNodeId()]
         .getVertexId()
       == vertexId1) {
      nodeList_[nodeId0].removeUpArcId(i);
      break;
    }
  }

  for(int i = 0; i < nodeList_[nodeId1].getNumberOfDownArcs(); i++) {
    if(nodeList_[arcList_[nodeList_[nodeId1].getDownArcId(i)].getDownNodeId()]
         .getVertexId()
       == vertexId0) {
      nodeList_[nodeId1].removeDownArcId(i);
      break;
    }
  }

  return 0;
}

int SubLevelSetTree::clearRegularNode(const int &vertexId) {

  if((vertexId < 0) || (vertexId >= vertexNumber_))
    return -1;

  int nodeId = vertex2node_[vertexId];

  if(!((nodeList_[nodeId].getNumberOfDownArcs() == 1)
       && (nodeList_[nodeId].getNumberOfUpArcs() == 1)))
    return -2;

  int downArcId = nodeList_[nodeId].getDownArcId(0),
      upArcId = nodeList_[nodeId].getUpArcId(0);

  int upNodeId = arcList_[upArcId].getUpNodeId();

  // removes the arcs from the current node
  nodeList_[nodeId].removeDownArcId(0);
  nodeList_[nodeId].removeUpArcId(0);

  // update the arc from the down node
  arcList_[downArcId].setUpNodeId(upNodeId);

  // remove the arc from the up node
  for(int i = 0; i < nodeList_[upNodeId].getNumberOfDownArcs(); i++) {
    if(nodeList_[upNodeId].getDownArcId(i) == upArcId) {
      nodeList_[upNodeId].removeDownArcId(i);
      break;
    }
  }
  // and add the updated arc
  nodeList_[upNodeId].addDownArcId(downArcId);

  return 0;
}

int SubLevelSetTree::clearRoot(const int &vertexId) {

  if((vertexId < 0) || (vertexId >= vertexNumber_))
    return -1;

  int nodeId = vertex2node_[vertexId];

  if((nodeList_[nodeId].getNumberOfDownArcs() == 1)
     && (!nodeList_[nodeId].getNumberOfUpArcs())) {
    clearArc(
      nodeList_[arcList_[nodeList_[nodeId].getDownArcId(0)].getDownNodeId()]
        .getVertexId(),
      vertexId);
  }

  return 0;
}

int SubLevelSetTree::exportArcPosToVtk(const int &arcId,
                                       const int &pointId,
                                       vector<int> &vertexIds,
                                       const vector<float> *origin,
                                       const vector<float> *voxelSize,
                                       ofstream &o) {

  vector<float> myOrigin(3), myVoxelSize(3);

  if(origin) {
    myOrigin = *(origin);
  } else {
    myOrigin[0] = myOrigin[1] = myOrigin[2] = 0;
  }

  if(voxelSize) {
    myVoxelSize = *(voxelSize);
  } else {
    myVoxelSize[0] = myVoxelSize[1] = myVoxelSize[2] = 1;
  }

  double offset = myVoxelSize[0];

  vector<double> v;

  int downNodeId = superArcList_[arcId].downNodeId_;
  int upNodeId = superArcList_[arcId].upNodeId_;

  v = (*vertexPositions_)[nodeList_[downNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId);

  v = (*vertexPositions_)[nodeList_[downNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 1);

  v = (*vertexPositions_)[nodeList_[downNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 2);

  v = (*vertexPositions_)[nodeList_[downNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 3);

  v = (*vertexPositions_)[nodeList_[upNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 4);

  v = (*vertexPositions_)[nodeList_[upNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 5);

  v = (*vertexPositions_)[nodeList_[upNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 6);

  v = (*vertexPositions_)[nodeList_[upNodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 7);

  return 0;
}

int SubLevelSetTree::exportNodeColorToVtk(const int &nodeId, ofstream &o) {

  for(int i = 0; i < 8; i++) {

    if(!nodeList_[nodeId].downSuperArcList_.size()) {
      // extremum
      if(minimumList_) {
        // minimum
        o << "0 ";
      } else {
        o << "1 ";
      }
    } else if(!nodeList_[nodeId].upSuperArcList_.size()) {
      if(minimumList_) {
        // maximum
        o << "1 ";
      } else {
        o << "0 ";
      }
    } else {
      o << "0.5 ";
    }
  }

  return 0;
}

int SubLevelSetTree::exportNodePosToVtk(const int &nodeId,
                                        const int &pointId,
                                        vector<int> &vertexIds,
                                        const vector<float> *origin,
                                        const vector<float> *voxelSize,
                                        ofstream &o) {

  vector<float> myOrigin(3), myVoxelSize(3);

  if(origin) {
    myOrigin = *(origin);
  } else {
    myOrigin[0] = myOrigin[1] = myOrigin[2] = 0;
  }

  if(voxelSize) {
    myVoxelSize = *(voxelSize);
  } else {
    myVoxelSize[0] = myVoxelSize[1] = myVoxelSize[2] = 1;
  }

  double offset = myVoxelSize[0];

  vector<double> v;

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }

  v[0] -= offset;
  v[1] -= offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[1] -= offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 1);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[1] += offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 2);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[1] += offset;
  v[2] -= offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 3);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[1] -= offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 4);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[1] -= offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 5);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] += offset;
  v[1] += offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 6);

  v = (*vertexPositions_)[nodeList_[nodeId].vertexId_];
  // fix a bug in paraview
  for(int i = 0; i < 3; i++) {
    v[i] -= myOrigin[i];
    v[i] /= myVoxelSize[i];
  }
  v[0] -= offset;
  v[1] += offset;
  v[2] += offset;
  o << v[0] << " " << v[1] << " " << v[2] << endl;
  vertexIds.push_back(pointId + 7);

  return 0;
}

int SubLevelSetTree::exportPersistenceCurve(const string &fileName) const {
  vector<pair<double, int>> persistencePlot;

  getPersistencePlot(persistencePlot);

  ofstream file(fileName.data(), ios::out);

  if(!file) {
    return -1;
  }

  for(int i = 0; i < (int)persistencePlot.size(); i++) {

    file << setprecision(REAL_SIGNIFICANT_DIGITS) << persistencePlot[i].first
         << " " << persistencePlot[i].second << endl;
  }

  file.close();

  return 0;
}

int SubLevelSetTree::exportPersistenceDiagram(const string &fileName) const {

  vector<pair<double, double>> diagram;

  getPersistenceDiagram(diagram);

  ofstream file(fileName.data(), ios::out);

  if(!file) {
    return -1;
  }

  for(int i = 0; i < (int)diagram.size(); i++) {
    file << setprecision(REAL_SIGNIFICANT_DIGITS) << diagram[i].first << " "
         << diagram[i].first << endl;
    file << setprecision(REAL_SIGNIFICANT_DIGITS) << diagram[i].first << " "
         << diagram[i].second << endl;
    file << setprecision(REAL_SIGNIFICANT_DIGITS) << diagram[i].first << " "
         << diagram[i].first << endl;
  }

  file.close();

  return 0;
}

int SubLevelSetTree::exportToSvg(const string &fileName,
                                 const double &scaleX,
                                 const double &scaleY) {

  Timer t;

  if(!vertexScalars_)
    return -1;

  bool hasLayout = buildPlanarLayout(scaleX, scaleY);

  // put the root right above the highest saddle
  double minY = -1, maxY = -1;
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    const SuperArc *a = &(superArcList_[i]);
    const Node *n = &(nodeList_[a->getDownNodeId()]);
    if((!i) || (minY > n->layoutY_)) {
      minY = n->layoutY_;
    }
    if((!i) || (maxY < n->layoutY_)) {
      maxY = n->layoutY_;
    }

    n = &(nodeList_[a->getUpNodeId()]);
    if((!i) || (minY > n->layoutY_)) {
      minY = n->layoutY_;
    }
    if((!i) || (maxY < n->layoutY_)) {
      maxY = n->layoutY_;
    }
  }

  for(int i = 0; i < (int)superArcList_.size(); i++) {
    const SuperArc *a = &(superArcList_[i]);
    Node *upNode = &(nodeList_[a->getUpNodeId()]);
    Node *downNode = &(nodeList_[a->getDownNodeId()]);

    if(!upNode->getNumberOfUpSuperArcs()) {
      // root
      if(minimumList_) {
        // join tree
        upNode->layoutY_
          = 0.4 * (downNode->layoutY_ - minY) + downNode->layoutY_;
        break;
      }
      if(maximumList_) {
        // split tree
        upNode->layoutY_
          = downNode->layoutY_ - 0.4 * (maxY - downNode->layoutY_);
        break;
      }
    }
  }

  string dotFileName = fileName + ".dot";

  ofstream dotFile(dotFileName.data(), ios::out);

  if(!dotFile) {
    this->printErr("Could not open file `" + dotFileName + "'!");
    return -2;
  }

  dotFile << "digraph \"";
  if(minimumList_)
    dotFile << "Join";
  else
    dotFile << "Split";
  dotFile << " Tree\"{" << endl;

  double minValue = 0, maxValue = 0;
  for(int i = 0; i < (int)vertexScalars_->size(); i++) {
    if((!i) || (minValue > (*vertexScalars_)[i])) {
      minValue = (*vertexScalars_)[i];
    }
    if((!i) || (maxValue < (*vertexScalars_)[i])) {
      maxValue = (*vertexScalars_)[i];
    }
  }

  // super-arc first!
  for(int i = 0; i < (int)superArcList_.size(); i++) {

    if(!superArcList_[i].pruned_) {

      int downNodeId = superArcList_[i].getDownNodeId();
      int upNodeId = superArcList_[i].getUpNodeId();

      dotFile << "  \"f = ";
      dotFile << (*vertexScalars_)[nodeList_[downNodeId].vertexId_];
      if(vertexPositions_) {
        dotFile << "\\n p = ("
                << (*vertexPositions_)[nodeList_[downNodeId].vertexId_][0]
                << " "
                << (*vertexPositions_)[nodeList_[downNodeId].vertexId_][1]
                << " "
                << (*vertexPositions_)[nodeList_[downNodeId].vertexId_][2]
                << ")";
      }
      dotFile << "\" -> \"f = ";
      dotFile << (*vertexScalars_)[nodeList_[upNodeId].vertexId_];
      if(vertexPositions_) {
        dotFile << "\\n p = ("
                << (*vertexPositions_)[nodeList_[upNodeId].vertexId_][0] << " "
                << (*vertexPositions_)[nodeList_[upNodeId].vertexId_][1] << " "
                << (*vertexPositions_)[nodeList_[upNodeId].vertexId_][2] << ")";
      }
      dotFile << "\"" << endl;
    }
  }

  for(int i = 0; i < (int)nodeList_.size(); i++) {

    if((nodeList_[i].downSuperArcList_.size())
       || (nodeList_[i].upSuperArcList_.size())) {
      // not a regular node
      if(!nodeList_[i].pruned_) {
        dotFile << "  \"f = ";
        dotFile << (*vertexScalars_)[nodeList_[i].vertexId_];
        if(vertexPositions_) {
          dotFile << "\\n p = ("
                  << (*vertexPositions_)[nodeList_[i].vertexId_][0] << " "
                  << (*vertexPositions_)[nodeList_[i].vertexId_][1] << " "
                  << (*vertexPositions_)[nodeList_[i].vertexId_][2] << ")";
        }
        dotFile << "\" [";

        // properties
        if(minimumList_) {
          // join tree

          // color
          if(!nodeList_[i].upSuperArcList_.size()) {
            // global maximum --> green
            dotFile << "fillcolor=green,";
          } else if(!nodeList_[i].downSuperArcList_.size()) {
            // minimum --> blue
            dotFile << "fillcolor=blue,fontcolor=white,";
          } else {
            dotFile << "fillcolor=lightyellow,shape=diamond,";
          }

        } else {
          // split tree

          // color
          if(!nodeList_[i].upSuperArcList_.size()) {
            // global minimum --> blue
            dotFile << "fillcolor=blue,fontcolor=white,";
          } else if(!nodeList_[i].downSuperArcList_.size()) {
            // maximum --> green
            dotFile << "fillcolor=green,";
          } else {
            dotFile << "fillcolor=lightyellow,shape=diamond,";
          }
        }

        dotFile << "style=filled,";

        if(hasLayout) {
          dotFile << "pos=\"" << nodeList_[i].layoutX_ << ", "
                  << nodeList_[i].layoutY_ << "!\"";
        }

        dotFile << "]";

        dotFile << endl;
      }
    }
  }

  dotFile << "}" << endl;

  dotFile.close();

  stringstream commandLine;
  commandLine << "dot -Kneato -Tsvg " << dotFileName << " -o " << fileName
              << " &> /dev/null";
  this->printMsg("Calling GraphViz to generate the SVG file.");
  this->printMsg("This may take a long time...");

  int cmdRet = system(commandLine.str().data());

  if((!cmdRet) || (cmdRet == 256)) {
    this->printMsg(
      "Output file `" + fileName + "' generated", 1.0, t.getElapsedTime(), 1);
  } else {
    this->printErr("Could not find the `GraphViz' package!");
    this->printErr("Please install it and re-run this program.");
  }

  //   OsCall::rmFile(dotFileName);

  return 0;
}

int SubLevelSetTree::exportToVtk(const string &fileName,
                                 const vector<float> *origin,
                                 const vector<float> *voxelSize) {

  if(!vertexScalars_)
    return -1;

  if(!vertexPositions_)
    return -2;

  ofstream o(fileName.data(), ios::out);

  if(!o) {
    this->printErr("Could not open file `" + fileName + "'!");
    return -3;
  }

  int superArcNumber = 0;
  int nodeNumber = 0;
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      superArcNumber++;
    }
  }
  nodeNumber = superArcNumber + 1;

  int pointNumber = 0, cellNumber = 0;
  pointNumber = 8 * nodeNumber + 8 * superArcNumber;
  cellNumber = nodeNumber + superArcNumber;

  vector<bool> nodeColorOut(nodeList_.size(), false);
  vector<bool> nodePosOut(nodeList_.size(), false);
  vector<bool> nodeMeshOut(nodeList_.size(), false);
  vector<vector<int>> nodeIds(nodeList_.size());
  vector<vector<int>> arcIds(superArcList_.size());

  o << "<?xml version=\"1.0\"?>" << endl;
  o << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
    << " byte_order=\"LittleEndian\">" << endl;
  o << "  <UnstructuredGrid>" << endl;
  o << "    <Piece NumberOfPoints=\"" << pointNumber << "\" NumberOfCells=\""
    << cellNumber << "\">";
  o << "      <PointData Scalars=\"Color Code\">" << endl;
  o << "        <DataArray type=\"Float32\" Name=\"Critical Color\""
    << " format=\"ascii\" NumberOfComponents=\"1\">" << endl;

  // node color
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      int downNodeId = superArcList_[i].downNodeId_;
      int upNodeId = superArcList_[i].upNodeId_;

      if(!nodeColorOut[downNodeId]) {
        o << "          ";
        exportNodeColorToVtk(downNodeId, o);
        o << endl;
        nodeColorOut[downNodeId] = true;
      }
      if(!nodeColorOut[upNodeId]) {
        o << "          ";
        exportNodeColorToVtk(upNodeId, o);
        o << endl;
        nodeColorOut[upNodeId] = true;
      }
    }
  }

  // arc color
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      o << "\t\t\t\t\t";
      for(int j = 0; j < 8; j++) {
        o << "2 ";
      }
      o << endl;
    }
  }
  o << "          </DataArray>" << endl;
  o << "        </PointData>" << endl;

  // node position

  o << "        <Points>" << endl;
  o << "          <DataArray type=\"Float32\" NumberOfComponents=\"3\""
    << " format=\"ascii\">" << endl;

  int pointId = 0;
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      int downNodeId = superArcList_[i].downNodeId_;
      int upNodeId = superArcList_[i].upNodeId_;

      if(!nodePosOut[downNodeId]) {
        o << "            ";
        exportNodePosToVtk(
          downNodeId, pointId, nodeIds[downNodeId], origin, voxelSize, o);
        nodePosOut[downNodeId] = true;
        o << endl;
        pointId += 8;
      }
      if(!nodePosOut[upNodeId]) {
        o << "            ";
        exportNodePosToVtk(
          upNodeId, pointId, nodeIds[upNodeId], origin, voxelSize, o);
        nodePosOut[upNodeId] = true;
        o << endl;
        pointId += 8;
      }
    }
  }
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      o << "            ";
      exportArcPosToVtk(i, pointId, arcIds[i], origin, voxelSize, o);
      o << endl;
      pointId += 8;
    }
  }
  o << "          </DataArray>" << endl;
  o << "        </Points>" << endl;

  // cells now
  o << "      <Cells>" << endl;

  o << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    << " format=\"ascii\">" << endl;

  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      int downNodeId = superArcList_[i].downNodeId_;
      int upNodeId = superArcList_[i].upNodeId_;

      if(!nodeMeshOut[downNodeId]) {
        o << "        " << nodeIds[downNodeId][0] << " "
          << nodeIds[downNodeId][1] << " " << nodeIds[downNodeId][2] << " "
          << nodeIds[downNodeId][3] << " " << nodeIds[downNodeId][4] << " "
          << nodeIds[downNodeId][5] << " " << nodeIds[downNodeId][6] << " "
          << nodeIds[downNodeId][7] << endl;
        nodeMeshOut[downNodeId] = true;
      }
      if(!nodeMeshOut[upNodeId]) {
        o << "        " << nodeIds[upNodeId][0] << " " << nodeIds[upNodeId][1]
          << " " << nodeIds[upNodeId][2] << " " << nodeIds[upNodeId][3] << " "
          << nodeIds[upNodeId][4] << " " << nodeIds[upNodeId][5] << " "
          << nodeIds[upNodeId][6] << " " << nodeIds[upNodeId][7] << endl;
        nodeMeshOut[upNodeId] = true;
      }
    }
  }
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      o << "        " << arcIds[i][0] << " " << arcIds[i][1] << " "
        << arcIds[i][2] << " " << arcIds[i][3] << " " << arcIds[i][4] << " "
        << arcIds[i][5] << " " << arcIds[i][6] << " " << arcIds[i][7] << endl;
    }
  }

  o << "        </DataArray>" << endl;

  o << "        <DataArray type=\"Int32\" Name=\"offsets\""
    << "  format=\"ascii\">" << endl;
  pointId = 8;
  for(int i = 0; i < nodeNumber; i++) {
    o << "        " << pointId << endl;
    pointId += 8;
  }
  for(int i = 0; i < superArcNumber; i++) {
    o << "        " << pointId << endl;
    pointId += 8;
  }
  o << "        </DataArray>" << endl;

  o << "        <DataArray type=\"Int32\" Name=\"types\""
    << " format=\"ascii\">" << endl;
  for(int i = 0; i < nodeNumber; i++) {
    o << "          12" << endl;
  }
  for(int i = 0; i < superArcNumber; i++) {
    o << "          12" << endl;
  }
  o << "        </DataArray>" << endl;

  o << "      </Cells>" << endl;
  o << "    </Piece>" << endl;
  o << "  </UnstructuredGrid>" << endl;
  o << "</VTKFile>" << endl;
  o.close();

  return 0;
}

int SubLevelSetTree::flush() {

  nodeList_.clear();
  arcList_.clear();
  superArcList_.clear();
  vertex2node_.clear();
  vertex2superArc_.clear();
  vertex2superArcNode_.clear();

  vertex2node_.resize(vertexNumber_, -1);
  vertex2superArc_.resize(vertexNumber_, -1);
  vertex2superArcNode_.resize(vertexNumber_, -1);

  return 0;
}

int SubLevelSetTree::buildExtremumList(vector<int> &extremumList,
                                       const bool &isSubLevelSet) {

  // TODO uncomment after triangulation templatization
  if((!triangulation_) /*|| triangulation_->isEmpty()*/)
    return -1;

  if((!vertexScalars_) || (!vertexScalars_->size()))
    return -2;

  if((!vertexSoSoffsets_) || (!vertexSoSoffsets_->size()))
    return -3;

  vector<pair<bool, pair<double, pair<int, int>>>> tmpList;

  for(SimplexId i = 0; i < triangulation_->getNumberOfVertices(); i++) {

    bool isExtremum = true;
    SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(i);
    for(SimplexId j = 0; j < neighborNumber; j++) {
      SimplexId otherId;
      triangulation_->getVertexNeighbor(i, j, otherId);

      if(isSubLevelSet) {
        // looking for minima
        if(isSosHigherThan(i, otherId)) {
          // not a minimum then
          isExtremum = false;
          break;
        }
      } else {
        // looking for maxima
        if(!isSosHigherThan(i, otherId)) {
          // not a maximum then
          isExtremum = false;
          break;
        }
      }
    }
    if(isExtremum) {
      pair<bool, pair<double, pair<int, int>>> entry;
      entry.first = isSubLevelSet;
      entry.second.first = (*vertexScalars_)[i];
      entry.second.second.first = (*vertexSoSoffsets_)[i];
      entry.second.second.second = i;
      tmpList.push_back(entry);
    }
  }

  sort(tmpList.begin(), tmpList.end(), filtrationCmp);

  extremumList.resize(tmpList.size());
  for(int i = 0; i < (int)extremumList.size(); i++) {
    extremumList[i] = tmpList[i].second.second.second;
  }

  if(isSubLevelSet) {
    minimumList_ = &(extremumList);
  } else {
    maximumList_ = &(extremumList);
  }

  return 0;
}

int SubLevelSetTree::closeSuperArc(const int &superArcId, const int &nodeId) {

  if((superArcId < 0) || (superArcId >= (int)superArcList_.size()))
    return -1;

  if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
    return -2;

  superArcList_[superArcId].setUpNodeId(nodeId);

  nodeList_[nodeId].addDownSuperArcId(superArcId);

  // create a low level arc with the previous node
  int nId = -1;
  if(!superArcList_[superArcId].getNumberOfRegularNodes()) {
    // the down neighbor is the arc extremity
    nId = superArcList_[superArcId].getDownNodeId();
  } else {
    nId = superArcList_[superArcId].getRegularNodeId(
      superArcList_[superArcId].getNumberOfRegularNodes() - 1);
  }

  makeArc(nId, nodeId);

  return 0;
}

int SubLevelSetTree::getPersistenceDiagram(
  vector<pair<double, double>> &diagram,
  vector<pair<pair<int, int>, double>> *pairs) const {

  vector<pair<pair<int, int>, double>> defaultPairs{};
  if(!pairs) {
    pairs = &defaultPairs;
  }

  if(!pairs->size())
    getPersistencePairs(*pairs);

  // fast fix :(
  diagram.resize(pairs->size());

  for(int i = 0; i < (int)pairs->size(); i++) {
    if((maximumList_) && (!minimumList_)) {
      // split tree
      diagram[i].second = (*vertexScalars_)[(*pairs)[i].first.first];
      diagram[i].first = (*vertexScalars_)[(*pairs)[i].first.second];
    } else {
      // join tree or contour tree
      diagram[i].first = (*vertexScalars_)[(*pairs)[i].first.first];
      diagram[i].second = (*vertexScalars_)[(*pairs)[i].first.second];
    }
  }

  std::sort(diagram.begin(), diagram.end(), _pPairCmp);

  return 0;
}

int SubLevelSetTree::getPersistencePairs(
  vector<pair<pair<int, int>, double>> &pairs,
  std::vector<std::pair<std::pair<int, int>, double>> *ttkNotUsed(mergePairs),
  std::vector<std::pair<std::pair<int, int>, double>> *ttkNotUsed(
    splitPairs)) const {

  Timer t;

  if(!superArcList_.size())
    return -1;

  bool isMergeTree = true;

  vector<int> *extremumList = minimumList_;
  if(!extremumList) {
    extremumList = maximumList_;
    isMergeTree = false;
  }
  if(!extremumList)
    return -2;

  pairs.resize(extremumList->size());

  vector<int> min2arc(vertexNumber_, -1);

  // compute the initial persistence
  int leafId = 0;
  for(int i = 0; i < getNumberOfSuperArcs(); ++i) {
    const SuperArc *a = getSuperArc(i);
    const Node *n = getNode(a->getDownNodeId());
    int nodeId = n->getVertexId();

    if(!n->getNumberOfDownSuperArcs()) {
      // internal, split trees and merges trees have their leaves at the bottom
      double extremumScalar = 0;
      double saddleScalar = 0;
      int saddleId = getNode(a->getUpNodeId())->getVertexId();

      extremumScalar = (*vertexScalars_)[nodeId];
      saddleScalar = (*vertexScalars_)[saddleId];

      pairs[leafId].first.first = nodeId;
      pairs[leafId].first.second = saddleId;
      pairs[leafId].second = fabs(saddleScalar - extremumScalar);
      float persistence = pairs[leafId].second;
      if(isnan(persistence))
        pairs[leafId].second = 0;

      // save that in the map
      min2arc[nodeId] = i;
      leafId++;
    }
  }

  // sort the pairs and start the processing
  sort(pairs.begin(), pairs.end(), _pCmp);

  vector<int> pairedSaddles(vertexNumber_, 0);
  // now start the pairing...
  for(unsigned int i = 0; i < pairs.size(); ++i) {
    pairedSaddles[pairs[i].first.second]++;

    // now update the persistence of all others
    for(unsigned int j = i + 1; j < pairs.size(); ++j) {
      if((pairs[j].first.second == pairs[i].first.second)
         && (pairedSaddles[pairs[i].first.second] + 1
             == (getVertexNode(pairs[i].first.second))
                  ->getNumberOfDownSuperArcs())) {
        // that pair contains a minimum competing with our guy
        // and the saddle is fully paired

        const SuperArc *a = getSuperArc(min2arc[pairs[j].first.first]);
        do {
          if(!(getNode(a->getUpNodeId())->getNumberOfUpSuperArcs()))
            break;

          int nextArcId = getNode(a->getUpNodeId())->getUpSuperArcId(0);
          a = getSuperArc(nextArcId);

          if(a) {
            int saddleId = getNode(a->getUpNodeId())->getVertexId();
            if(!getNode(a->getUpNodeId())->getNumberOfDownSuperArcs())
              break;

            if(pairedSaddles[saddleId]
               < (getVertexNode(saddleId))->getNumberOfDownSuperArcs() - 1) {

              double extremumScalar = 0;
              double saddleScalar = 0;

              extremumScalar = (*vertexScalars_)[pairs[j].first.first];
              saddleScalar = (*vertexScalars_)[saddleId];
              pairs[j].first.second = saddleId;
              pairs[j].second = fabs(saddleScalar - extremumScalar);
              float persistence = pairs[j].second;
              if(isnan(persistence))
                pairs[j].second = 0;
              break;
            }
          }
        } while(a);

        if(a) {
          // we indeed update the persistence of the competitor
          sort(pairs.begin() + j, pairs.end(), _pCmp);
        }

        break;
      }
    }
  }

  // the last entry should have infinite persistence
  // pairs.erase(pairs.end() - 1);

  // let's pair the global maximum with the global min
  for(int i = 0; i < getNumberOfSuperArcs(); i++) {
    const SuperArc *a = getSuperArc(i);
    const Node *up = getNode(a->getUpNodeId());
    if(!up->getNumberOfUpSuperArcs()) {
      // global min
      pairs.back().first.second = up->getVertexId();
      pairs.back().second
        = fabs((*vertexScalars_)[pairs.back().first.first]
               - (*vertexScalars_)[pairs.back().first.second]);
      float persistence = pairs.back().second;
      if(isnan(persistence))
        pairs.back().second = 0;
      break;
    }
  }

  std::string pairtype = isMergeTree ? "(0-1)" : "(1-2)";
  std::string pairextr = isMergeTree ? "min" : "max";

  this->printMsg(
    std::vector<std::vector<std::string>>{
      {"#" + pairtype + " pairs", std::to_string(pairs.size())}},
    debug::Priority::DETAIL);
  for(const auto &p : pairs) {
    this->printMsg(
      std::vector<std::vector<std::string>>{
        {"#" + pairextr, std::to_string(p.first.first)},
        {"#saddles", std::to_string(p.first.second)},
        {"Persistence", std::to_string(p.second)}},
      debug::Priority::DETAIL);
  }

  this->printMsg(std::to_string(pairs.size()) + " persistence pairs computed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

int SubLevelSetTree::getPersistencePlot(
  vector<pair<double, int>> &plot,
  vector<pair<pair<int, int>, double>> *persistencePairs) const {

  vector<pair<pair<int, int>, double>> defaultPersistencePairs{};
  if(!persistencePairs) {
    persistencePairs = &defaultPersistencePairs;
  }

  if(!persistencePairs->size())
    getPersistencePairs(*persistencePairs);

  plot.resize(persistencePairs->size());

  for(int i = 0; i < (int)plot.size(); i++) {
    plot[i].first = (*persistencePairs)[i].second;
    if(plot[i].first < Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS)) {
      plot[i].first = Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS);
    }
    plot[i].second = persistencePairs->size() - i;
  }

  return 0;
}

bool SubLevelSetTree::buildPlanarLayout(const double &scaleX,
                                        const double &scaleY) {

  if((minimumList_) && (maximumList_)) {
    this->printErr("Contour tree planar layout not implemented.");
    this->printErr("Planar layout is only implemented for merge-trees.");
    return false;
  }

  if((!vertexScalars_) || (!vertexScalars_->size()))
    return false;

  // global dimansions [0,1]x[0,1]

  // 1) set the y coordinate for everyone
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    int downNodeId = superArcList_[i].downNodeId_;
    int upNodeId = superArcList_[i].upNodeId_;

    nodeList_[downNodeId].layoutY_
      = ((*vertexScalars_)[nodeList_[downNodeId].vertexId_] - minScalar_)
        / (maxScalar_ - minScalar_);

    nodeList_[upNodeId].layoutY_
      = ((*vertexScalars_)[nodeList_[upNodeId].vertexId_] - minScalar_)
        / (maxScalar_ - minScalar_);
  }

  // 2) for each node count the number of leaves under it
  //    --> start a front from the leaves and increment a table
  vector<bool> inFront(nodeList_.size(), false);
  vector<int> node2LeafNumber(nodeList_.size(), 0);
  // use the same data-structure for front propagation
  set<pair<bool, pair<double, pair<int, int>>>, filtrationCtCmp> front;

  // 2a) start from the leaves
  for(int i = 0; i < (int)superArcList_.size(); i++) {

    if(!superArcList_[i].pruned_) {

      int downNodeId = superArcList_[i].downNodeId_;
      if(!nodeList_[downNodeId].downSuperArcList_.size()) {
        pair<bool, pair<double, pair<int, int>>> n;
        n.first = (minimumList_ ? true : false);
        n.second.first = (*vertexScalars_)[nodeList_[downNodeId].vertexId_];
        n.second.second.first
          = (*vertexSoSoffsets_)[nodeList_[downNodeId].vertexId_];
        n.second.second.second = nodeList_[downNodeId].vertexId_;
        front.insert(n);

        node2LeafNumber[downNodeId] = 1;
        inFront[downNodeId] = true;
      }
    }
  }

  // filtration loop
  do {

    int vertexId = front.begin()->second.second.second;
    front.erase(front.begin());

    int nodeId = vertex2node_[vertexId];

    // retrieve the number of leaves given by its children
    for(int i = 0; i < (int)nodeList_[nodeId].downSuperArcList_.size(); i++) {
      int arcId = nodeList_[nodeId].downSuperArcList_[i];
      int childId = superArcList_[arcId].downNodeId_;

      node2LeafNumber[nodeId] += node2LeafNumber[childId];
    }

    // add its parents to the front
    for(int i = 0; i < (int)nodeList_[nodeId].upSuperArcList_.size(); i++) {
      int arcId = nodeList_[nodeId].upSuperArcList_[i];
      int parentId = superArcList_[arcId].upNodeId_;

      if(!inFront[parentId]) {
        pair<bool, pair<double, pair<int, int>>> n;
        n.first = (minimumList_ ? true : false);
        n.second.first = (*vertexScalars_)[nodeList_[parentId].vertexId_];
        n.second.second.first
          = (*vertexSoSoffsets_)[nodeList_[parentId].vertexId_];
        n.second.second.second = nodeList_[parentId].vertexId_;
        front.insert(n);

        node2LeafNumber[parentId] = 0;
        inFront[parentId] = true;
      }
    }

  } while(!front.empty());

  // 3) start from the node below the root, these two guys have the same x
  //   each node has its X coordinate given by a ratio

  // let's find the root
  int rootId = -1;
  int subRootId = -1;

  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      int upNodeId = superArcList_[i].upNodeId_;
      if(!nodeList_[upNodeId].upSuperArcList_.size()) {
        // that's it
        rootId = upNodeId;
        break;
      }
    }
  }

  // in theory, the root has exactly one arc below it
  subRootId = superArcList_[nodeList_[rootId].downSuperArcList_[0]].downNodeId_;

  // now start the recursion
  vector<pair<double, double>> nodeInterval(
    nodeList_.size(), pair<double, double>(0, 1));
  queue<int> nodeQueue;

  //   printf("pushing root: %f\n",
  //     (*vertexScalars_)[nodeList_[rootId].vertexId_]);
  nodeQueue.push(rootId);

  do {

    int nodeId = nodeQueue.front();
    nodeQueue.pop();

    // get the number of leaves of the first child
    double ratio = 0.5;

    if((nodeList_[nodeId].downSuperArcList_.size()) && (nodeId != subRootId)
       && (nodeId != rootId)) {

      int arcId = nodeList_[nodeId].downSuperArcList_[0];
      int childId = superArcList_[arcId].downNodeId_;

      ratio = node2LeafNumber[childId] / ((double)node2LeafNumber[nodeId]);
      if(ratio > 1)
        ratio = 1;
    }

    //     printf("node %f: [%f-%f] r=%f l=%d\n",
    //       (*vertexScalars_)[nodeList_[nodeId].vertexId_],
    //       nodeInterval[nodeId].first, nodeInterval[nodeId].second, ratio,
    //       node2LeafNumber[nodeId]);
    nodeList_[nodeId].layoutX_
      = nodeInterval[nodeId].first
        + ratio * (nodeInterval[nodeId].second - nodeInterval[nodeId].first);
    if(nodeList_[nodeId].layoutX_ > 1) {
      nodeList_[nodeId].layoutX_ = 1;
    }

    // pass the right interval to the children and add them to the queue
    double xOffset = 0;
    for(int i = 0; i < (int)nodeList_[nodeId].downSuperArcList_.size(); i++) {
      int arcId = nodeList_[nodeId].downSuperArcList_[i];
      int childId = superArcList_[arcId].downNodeId_;

      ratio = node2LeafNumber[childId] / ((double)node2LeafNumber[nodeId]);
      if(ratio > 1)
        ratio = 1;

      nodeInterval[childId].first = nodeInterval[nodeId].first + xOffset;
      nodeInterval[childId].second
        = nodeInterval[childId].first
          + ratio * (nodeInterval[nodeId].second - nodeInterval[nodeId].first);
      xOffset = nodeInterval[childId].second;

      nodeQueue.push(childId);
    }

  } while(!nodeQueue.empty());

  // scale a bit
  int maintainedArcNumber = 0;
  for(int i = 0; i < (int)superArcList_.size(); i++) {
    if(!superArcList_[i].pruned_) {
      maintainedArcNumber++;
    }
  }

  // test
  //   for(int i = 0; i < (int) superArcList_.size(); i++){
  //     if(!superArcList_[i].pruned_){
  //       Node *downNode = &(nodeList_[superArcList_[i].getDownNodeId()]);
  //       printf("node %f --> %fx%f\n",
  //         (*vertexScalars_)[downNode->vertexId_],
  //         downNode->layoutX_, downNode->layoutY_);
  //     }
  //   }
  //

  for(int i = 0; i < (int)nodeList_.size(); i++) {
    nodeList_[i].layoutX_ *= 30 * scaleX;
    nodeList_[i].layoutY_ *= 1000 * scaleY;
  }

  return true;
}

int SubLevelSetTree::buildSaddleList(vector<int> &vertexList) const {

  vertexList.clear();

  for(int i = 0; i < (int)superArcList_.size(); i++) {
    const Node *down = &(nodeList_[superArcList_[i].getDownNodeId()]);

    if(down->getNumberOfDownSuperArcs()) {
      vertexList.push_back(down->getVertexId());
    }
  }

  return 0;
}

int SubLevelSetTree::makeArc(const int &nodeId0, const int &nodeId1) {

  if((nodeId0 < 0) || (nodeId0 >= (int)nodeList_.size()))
    return -1;
  if((nodeId1 < 0) || (nodeId1 >= (int)nodeList_.size()))
    return -2;

  arcList_.resize(arcList_.size() + 1);
  arcList_[arcList_.size() - 1].setDownNodeId(nodeId0);
  arcList_[arcList_.size() - 1].setUpNodeId(nodeId1);

  nodeList_[nodeId0].addUpArcId((int)arcList_.size() - 1);
  nodeList_[nodeId1].addDownArcId((int)arcList_.size() - 1);

  return (int)arcList_.size() - 1;
}

int SubLevelSetTree::makeNode(const int &vertexId) {

  if((vertexId < 0) || (vertexId >= vertexNumber_))
    return -1;

  if(vertex2node_[vertexId] == -1) {
    nodeList_.resize(nodeList_.size() + 1);
    nodeList_[nodeList_.size() - 1].setVertexId(vertexId);

    vertex2node_[vertexId] = (int)nodeList_.size() - 1;

    return (int)nodeList_.size() - 1;
  }

  return vertex2node_[vertexId];
}

int SubLevelSetTree::openSuperArc(const int &nodeId) {

  if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
    return -1;

  superArcList_.resize(superArcList_.size() + 1);
  superArcList_[superArcList_.size() - 1].setDownNodeId(nodeId);

  nodeList_[nodeId].addUpSuperArcId((int)superArcList_.size() - 1);

  return (int)superArcList_.size() - 1;
}

bool SubLevelSetTree::isSosHigherThan(const int &vertexId0,
                                      const int &vertexId1) const {

  return (((*vertexScalars_)[vertexId0] > (*vertexScalars_)[vertexId1])
          || (((*vertexScalars_)[vertexId0] == (*vertexScalars_)[vertexId1])
              && ((*vertexSoSoffsets_)[vertexId0]
                  > (*vertexSoSoffsets_)[vertexId1])));
}

bool SubLevelSetTree::isSosLowerThan(const int &vertexId0,
                                     const int &vertexId1) const {

  return (((*vertexScalars_)[vertexId0] < (*vertexScalars_)[vertexId1])
          || (((*vertexScalars_)[vertexId0] == (*vertexScalars_)[vertexId1])
              && ((*vertexSoSoffsets_)[vertexId0]
                  < (*vertexSoSoffsets_)[vertexId1])));
}

int SubLevelSetTree::moveRegularNode(const Node *n,
                                     const Node *oldDown,
                                     const Node *oldUp,
                                     const Node *newDown,
                                     const Node *newUp) {

  int arcId = -1;
  int nodeId = n - &(nodeList_[0]);

  // remove n from oldDown and oldUp and update n
  for(int i = 0; i < oldUp->getNumberOfDownArcs(); i++) {
    arcId = oldUp->getDownArcId(i);
    if(arcList_[arcId].getDownNodeId() == nodeId) {
      arcList_[arcId].setDownNodeId(oldDown - &(nodeList_[0]));
      break;
    }
  }
  for(int i = 0; i < n->getNumberOfUpArcs(); i++) {
    if(arcId == n->getUpArcId(i)) {
      nodeList_[nodeId].removeUpArcId(i);
      break;
    }
  }

  for(int i = 0; i < oldDown->getNumberOfUpArcs(); i++) {
    if(arcList_[oldDown->getUpArcId(i)].getUpNodeId() == nodeId) {
      nodeList_[oldDown - &(nodeList_[0])].removeUpArcId(i);
      nodeList_[oldDown - &(nodeList_[0])].addUpArcId(arcId);
      break;
    }
  }
  for(int i = 0; i < n->getNumberOfDownArcs(); i++) {
    if(arcList_[n->getDownArcId(i)].getDownNodeId()
       == (oldDown - &(nodeList_[0]))) {
      nodeList_[nodeId].removeDownArcId(i);
      break;
    }
  }

  // disconnect newDown and newUp
  if(newUp) {
    for(int i = 0; i < newDown->getNumberOfUpArcs(); i++) {
      arcId = newDown->getUpArcId(i);
      if(arcList_[arcId].getUpNodeId() == (newUp - &(nodeList_[0]))) {
        nodeList_[newDown - &(nodeList_[0])].removeUpArcId(i);
        break;
      }
    }
    for(int i = 0; i < newUp->getNumberOfDownArcs(); i++) {
      if(newUp->getDownArcId(i) == arcId) {
        nodeList_[newUp - &(nodeList_[0])].removeDownArcId(i);
        break;
      }
    }
  }

  makeArc(newDown - &(nodeList_[0]), n - &(nodeList_[0]));
  if(newUp)
    makeArc(n - &(nodeList_[0]), newUp - &(nodeList_[0]));

  return 0;
}

int SubLevelSetTree::print() const {

  stringstream msg;

  this->printMsg(
    "Node list (" + std::to_string(getNumberOfNodes()) + " nodes):",
    debug::Priority::DETAIL);

  int minCount = 0, saddleCount = 0, maxCount = 0, regularCount = 0;

  for(int i = 0; i < getNumberOfNodes(); i++) {
    const Node *n = getNode(i);

    if(vertex2superArc_[n->getVertexId()] == -1) {
      // regular node

      this->printMsg(
        std::vector<std::vector<std::string>>{
          {"Id", std::to_string(i), "VertId", std::to_string(n->getVertexId()),
           "D", std::to_string(n->getNumberOfDownSuperArcs()), "U",
           std::to_string(n->getNumberOfUpSuperArcs())}},
        debug::Priority::DETAIL);

      if(!n->getNumberOfDownSuperArcs())
        minCount++;
      else if(!n->getNumberOfUpSuperArcs())
        maxCount++;
      else
        saddleCount++;
    }
  }

  this->printMsg(
    "Arc list (" + std::to_string(getNumberOfSuperArcs()) + " arcs):",
    debug::Priority::DETAIL);

  for(int i = 0; i < getNumberOfSuperArcs(); i++) {
    const SuperArc *a = getSuperArc(i);

    const Node *down = getNode(a->getDownNodeId()),
               *up = getNode(a->getUpNodeId());
    if((up) && (down)) {

      this->printMsg(
        std::vector<std::vector<std::string>>{
          {"Id", std::to_string(i), "D",
           std::to_string(down->getNumberOfDownSuperArcs()), "U",
           std::to_string(up->getNumberOfUpSuperArcs()), "V",
           std::to_string(a->getNumberOfRegularNodes())}},
        debug::Priority::DETAIL);
    } else {
      this->printErr("Arc inconsistency! " + std::to_string(a->getDownNodeId())
                     + "->" + std::to_string(a->getUpNodeId()));
    }

    for(int j = 0; j < a->getNumberOfRegularNodes(); j++) {
      this->printMsg(
        std::vector<std::vector<std::string>>{
          {"Regular vertex",
           std::to_string(nodeList_[a->getRegularNodeId(j)].getVertexId())}},
        debug::Priority::VERBOSE);
    }

    regularCount += a->getNumberOfRegularNodes();
  }

  this->printMsg(
    std::vector<std::vector<std::string>>{
      {"#Minima", std::to_string(minCount)},
      {"#Saddles", std::to_string(saddleCount)},
      {"#Maxima", std::to_string(maxCount)},
      {"#Regular", std::to_string(regularCount)},
      {"#Sum",
       std::to_string(minCount + saddleCount + maxCount + regularCount)},
      {"#Input vertices", std::to_string(vertexNumber_)},
    },
    debug::Priority::DETAIL);

  return 0;
}

int SubLevelSetTree::simplify(const double &simplificationThreshold,
                              ContourTreeSimplificationMetric *metric) {

  Timer t;

  //   if((simplificationThreshold < 0)||(simplificationThreshold > 1))
  //     return -1;

  PersistenceMetric defaultMetric;

  if(metric == nullptr) {
    metric = &defaultMetric;
  }

  metric->tree_ = this;

  //   double unNormalizedThreshold = simplificationThreshold
  //     * (maxScalar_ - minScalar_);

  if(!nodeList_.size())
    return -2;

  if(!superArcList_.size())
    return -3;

  if((minimumList_) && (maximumList_)) {
    this->printErr("Contour tree simplification not implemented.");
    this->printErr("Simplification is only implemented for merge-trees.");
    return -4;
  }

  int simplifiedArcNumber = 0;
  double maximumMetricScore = 0;

  if(!originalNodeList_.size()) {
    // first time we simplify
    originalNodeList_ = nodeList_;
    originalSuperArcList_ = superArcList_;
  } else {
    // not the first time
    nodeList_ = originalNodeList_;
    superArcList_ = originalSuperArcList_;
  }

  if(!simplificationThreshold)
    return 0;

  //   if(!unNormalizedThreshold) return -4;

  vector<pair<real, int>> regularNodeList;

  int arcToSimplify = -1;
  double currentScore = 0, minScore;

  do {

    arcToSimplify = -1;
    minScore = 1.1 * (maxScalar_ - minScalar_);

    // search for an elligible arc.
    for(int i = 0; i < (int)superArcList_.size(); i++) {

      if(!superArcList_[i].pruned_) {

        int downNodeId = superArcList_[i].downNodeId_;
        int upNodeId = superArcList_[i].upNodeId_;

        if(!nodeList_[downNodeId].downSuperArcList_.size()) {

          currentScore = metric->computeSuperArcMetric(
            nodeList_[downNodeId].vertexId_, nodeList_[upNodeId].vertexId_,
            superArcList_[i].regularNodeList_);

          if((currentScore < 0) || (currentScore > (maxScalar_ - minScalar_))) {
            this->printWrn("Out-of-range score! (user-defined metric)");
          } else {
            if(currentScore < minScore) {
              arcToSimplify = i;
              minScore = currentScore;
            }
          }
        }
      }
    }

    if(minScore > simplificationThreshold)
      break;

    if(minScore > maximumMetricScore)
      maximumMetricScore = minScore;

    if(arcToSimplify == -1)
      break;

    // perform the simplification
    int downNodeId = superArcList_[arcToSimplify].downNodeId_;
    int upNodeId = superArcList_[arcToSimplify].upNodeId_;

    int pivotId = upNodeId;
    int leafId = downNodeId;

    if(!nodeList_[downNodeId].downSuperArcList_.size()) {
      // remove an extremum
      pivotId = upNodeId;
      leafId = downNodeId;
    }

    bool isSimpleSaddle = (nodeList_[pivotId].downSuperArcList_.size()
                             + nodeList_[pivotId].upSuperArcList_.size()
                           == 3);

    int brotherId = 0;

    for(int i = 0; i < (int)nodeList_[pivotId].downSuperArcList_.size(); i++) {
      if(nodeList_[pivotId].downSuperArcList_[i] != arcToSimplify) {
        brotherId = nodeList_[pivotId].downSuperArcList_[i];
        break;
      }
    }

    int brotherExtremityId = superArcList_[brotherId].downNodeId_;

    int parentId = nodeList_[pivotId].upSuperArcList_[0];

    // 1) update the parent extremity (only if simple saddles)
    if(isSimpleSaddle) {
      superArcList_[parentId].downNodeId_ = brotherExtremityId;
    }

    // 2) tell the brotherExtremityId that it's no longer linked to the
    // brother but to the parent
    if(isSimpleSaddle) {
      for(int i = 0;
          i < (int)nodeList_[brotherExtremityId].upSuperArcList_.size(); i++) {
        if(nodeList_[brotherExtremityId].upSuperArcList_[i] == brotherId) {
          nodeList_[brotherExtremityId].upSuperArcList_[i] = parentId;
          break;
        }
      }
    } else {
      // multi-saddle, we need to erase arcToSimplify from the pivot
      for(int i = 0; i < (int)nodeList_[pivotId].downSuperArcList_.size();
          i++) {
        if(nodeList_[pivotId].downSuperArcList_[i] == arcToSimplify) {
          nodeList_[pivotId].downSuperArcList_.erase(
            nodeList_[pivotId].downSuperArcList_.begin() + i);
          break;
        }
      }
    }

    // 3) do the merging
    if(isSimpleSaddle) {
      // copy the regular nodes of the brother into the parent
      superArcList_[parentId].regularNodeList_.insert(
        superArcList_[parentId].regularNodeList_.end(),
        superArcList_[brotherId].regularNodeList_.begin(),
        superArcList_[brotherId].regularNodeList_.end());

      // copy the regular nodes of the arcToSimplify into the parent
      superArcList_[parentId].regularNodeList_.insert(
        superArcList_[parentId].regularNodeList_.end(),
        superArcList_[arcToSimplify].regularNodeList_.begin(),
        superArcList_[arcToSimplify].regularNodeList_.end());

      // add the pivot
      superArcList_[parentId].regularNodeList_.push_back(pivotId);

      // add the leaf
      superArcList_[parentId].regularNodeList_.push_back(leafId);

      // NOTE: optionally, the list of regular nodes could be sorted

      // mark the brother as pruned
      superArcList_[brotherId].pruned_ = true;
      nodeList_[pivotId].pruned_ = true;
      simplifiedArcNumber++;
    } else {
      // just merge with the first brother we found....
      superArcList_[brotherId].regularNodeList_.insert(
        superArcList_[brotherId].regularNodeList_.end(),
        superArcList_[arcToSimplify].regularNodeList_.begin(),
        superArcList_[arcToSimplify].regularNodeList_.end());

      // add the leaf
      superArcList_[brotherId].regularNodeList_.push_back(leafId);
    }
    superArcList_[arcToSimplify].pruned_ = true;
    nodeList_[leafId].pruned_ = true;
    simplifiedArcNumber++;

  } while(minScore < simplificationThreshold);

  // sort the regular nodes
  Timer sortTimer;
  for(int i = 0; i < getNumberOfSuperArcs(); ++i) {
    if(!superArcList_[i].pruned_) {
      superArcList_[i].sortRegularNodes(
        vertexScalars_, vertexSoSoffsets_, &nodeList_);
    }
  }

  this->printMsg("Regular node sorting", 1.0, sortTimer.getElapsedTime(),
                 this->threadNumber_);

  this->printMsg(std::to_string(simplifiedArcNumber)
                   + " super-arcs pruned (threshold="
                   + std::to_string(simplificationThreshold) + ")",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  this->printMsg(
    "Biggest simplification metric score: " + std::to_string(maximumMetricScore)
    + " ("
    + std::to_string(maximumMetricScore * 100 / (maxScalar_ - minScalar_))
    + "%)");

  return 0;
}

int SubLevelSetTree::computeBarycenters() {
  vector<int> sample;
  vector<double> barycenter(3);
  int vertexId;

  const SuperArc *a;
  for(int i = 0; i < getNumberOfSuperArcs(); ++i) {
    a = getSuperArc(i);
    if(!a->isPruned()) {
      for(int j = 0; j < a->getNumberOfSamples(); ++j) {
        a->getSample(j, sample);

        for(unsigned int k = 0; k < 3; ++k)
          barycenter[k] = 0;

        for(unsigned int k = 0; k < sample.size(); ++k) {
          vertexId = getNode(sample[k])->getVertexId();

          for(unsigned int l = 0; l < 3; ++l)
            barycenter[l] += (*vertexPositions_)[vertexId][l];
        }
        if(sample.size()) {
          for(unsigned int k = 0; k < 3; ++k)
            barycenter[k] /= sample.size();

          // update the arc
          superArcList_[i].appendBarycenter(barycenter);
        }
      }
    }
  }

  return 0;
}

int SubLevelSetTree::getSkeletonScalars(
  const vector<double> &scalars,
  vector<vector<double>> &skeletonScalars) const {
  skeletonScalars.clear();
  skeletonScalars.resize(getNumberOfSuperArcs());

  vector<int> sample;
  int nodeId;
  int vertexId;

  double f;
  double f0;
  double f1;
  double fmin;
  double fmax;
  int nodeMinId;
  int nodeMaxId;
  int nodeMinVId;
  int nodeMaxVId;
  const SuperArc *a;
  for(int i = 0; i < getNumberOfSuperArcs(); ++i) {
    a = getSuperArc(i);

    if(!a->isPruned()) {
      if(maximumList_ && !minimumList_) {
        nodeMinId = a->getUpNodeId();
        nodeMaxId = a->getDownNodeId();
      } else {
        nodeMaxId = a->getUpNodeId();
        nodeMinId = a->getDownNodeId();
      }

      nodeMaxVId = getNode(nodeMaxId)->getVertexId();
      nodeMinVId = getNode(nodeMinId)->getVertexId();

      fmax = scalars[nodeMaxVId];
      fmin = scalars[nodeMinVId];

      // init: min
      f0 = fmin;

      // iteration
      for(int j = 0; j < a->getNumberOfSamples(); ++j) {
        a->getSample(j, sample);

        f = 0;
        for(unsigned int k = 0; k < sample.size(); ++k) {
          nodeId = sample[k];
          vertexId = getNode(nodeId)->getVertexId();
          f += scalars[vertexId];
        }
        if(sample.size()) {
          f /= sample.size();

          f1 = f;
          // update the arc
          skeletonScalars[i].push_back((f0 + f1) / 2);
          f0 = f1;
        }
      }

      // end: max
      f1 = fmax;

      // update the arc
      skeletonScalars[i].push_back((f0 + f1) / 2);
    }
  }

  return 0;
}

int SubLevelSetTree::computeSkeleton(unsigned int arcResolution) {
  sample(arcResolution);
  computeBarycenters();

  isSkeletonComputed_ = true;
  return 0;
}

int SubLevelSetTree::smoothSkeleton(unsigned int skeletonSmoothing) {
  for(unsigned int i = 0; i < skeletonSmoothing; ++i) {
    for(int j = 0; j < getNumberOfSuperArcs(); ++j) {
      if(!superArcList_[j].isPruned()) {
        if(minimumList_)
          superArcList_[j].smooth(nodeList_, vertexPositions_, true);
        else
          superArcList_[j].smooth(nodeList_, vertexPositions_, false);
      }
    }
  }

  return 0;
}

int SubLevelSetTree::sample(unsigned int samplingLevel) {
  vector<vector<int>> sampleList(samplingLevel);

  const SuperArc *a;
  for(int i = 0; i < getNumberOfSuperArcs(); ++i) {
    a = getSuperArc(i);

    if(!a->isPruned()) {

      for(unsigned int j = 0; j < samplingLevel; ++j)
        sampleList[j].clear();

      double fmax, fmin;
      int nodeMaxId, nodeMinId;
      int nodeMaxVId, nodeMinVId;
      double delta;
      if(a->getNumberOfRegularNodes()) {
        if(minimumList_) {
          nodeMaxId = a->getUpNodeId();
          nodeMinId = a->getDownNodeId();
        } else {
          nodeMaxId = a->getDownNodeId();
          nodeMinId = a->getUpNodeId();
        }

        nodeMaxVId = getNode(nodeMaxId)->getVertexId();
        nodeMinVId = getNode(nodeMinId)->getVertexId();

        fmax = (*vertexScalars_)[nodeMaxVId];
        fmin = (*vertexScalars_)[nodeMinVId];

        delta = (fmax - fmin) / samplingLevel;

        double f;
        int nodeId;
        int vertexId;
        for(int j = 0; j < a->getNumberOfRegularNodes(); ++j) {
          nodeId = a->getRegularNodeId(j);
          vertexId = getNode(nodeId)->getVertexId();
          f = (*vertexScalars_)[vertexId];

          for(unsigned int k = 0; k < samplingLevel; ++k) {
            if(f <= (k + 1) * delta + fmin) {
              sampleList[k].push_back(nodeId);
              break;
            }
          }
        }

        // update the arc
        for(unsigned int j = 0; j < sampleList.size(); ++j)
          superArcList_[i].appendSample(sampleList[j]);
      }
    }
  }

  return 0;
}

int SubLevelSetTree::clearSkeleton() {
  for(int j = 0; j < getNumberOfSuperArcs(); ++j) {
    superArcList_[j].clearBarycenters();
    superArcList_[j].clearSamples();
  }

  return 0;
}

ContourTree::ContourTree() {
  this->setDebugMsgPrefix("ContourTree");
}

int ContourTree::build() {

  Timer timer;

  if(!vertexNumber_)
    return -1;
  if((!vertexScalars_) || ((int)vertexScalars_->size() != vertexNumber_))
    return -2;
  if(triangulation_->getNumberOfVertices() != vertexNumber_)
    return -3;

  mergeTree_.setDebugLevel(debugLevel_);
  splitTree_.setDebugLevel(debugLevel_);

  // 0) init data structures
  if(!vertexSoSoffsets_) {
    vertexSoSoffsets_ = new vector<int>;
    vertexSoSoffsets_->resize(vertexNumber_);
    for(int i = 0; i < (int)vertexSoSoffsets_->size(); i++)
      (*vertexSoSoffsets_)[i] = i;
  }

  // build the actual extrema list

  minimumList_ = new vector<int>;
  maximumList_ = new vector<int>;

  for(int i = 0; i < vertexNumber_; i++) {

    bool isMin = true, isMax = true;
    SimplexId neighborNumber = triangulation_->getVertexNeighborNumber(i);
    for(SimplexId j = 0; j < neighborNumber; j++) {
      SimplexId nId;
      triangulation_->getVertexNeighbor(i, j, nId);

      if(((*vertexScalars_)[nId] > (*vertexScalars_)[i])
         || (((*vertexScalars_)[nId] == (*vertexScalars_)[i])
             && ((*vertexSoSoffsets_)[nId] > (*vertexSoSoffsets_)[i])))
        isMax = false;

      if(((*vertexScalars_)[nId] < (*vertexScalars_)[i])
         || (((*vertexScalars_)[nId] == (*vertexScalars_)[i])
             && ((*vertexSoSoffsets_)[nId] < (*vertexSoSoffsets_)[i])))
        isMin = false;

      if((!isMin) && (!isMax))
        break;
    }

    if((isMin) && (!isMax))
      minimumList_->push_back(i);
    if((isMax) && (!isMin))
      maximumList_->push_back(i);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      // 1) build the merge tree
      mergeTree_.setMinimumList(*minimumList_);
      mergeTree_.setNumberOfVertices(vertexNumber_);
      mergeTree_.setVertexScalars(vertexScalars_);
      mergeTree_.setVertexPositions(vertexPositions_);
      mergeTree_.setTriangulation(triangulation_);
      mergeTree_.setVertexSoSoffsets(vertexSoSoffsets_);
      mergeTree_.build();
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    {
      // 2) build the split tree
      splitTree_.setMaximumList(*maximumList_);
      splitTree_.setNumberOfVertices(vertexNumber_);
      splitTree_.setVertexScalars(vertexScalars_);
      splitTree_.setVertexPositions(vertexPositions_);
      splitTree_.setVertexSoSoffsets(vertexSoSoffsets_);
      splitTree_.setTriangulation(triangulation_);
      splitTree_.build();
    }
  }

  // note: at this point, the split tree is layed out upside down.

  // 3) merge the two trees into the contour tree
  combineTrees();

  // 4) update the high level super arc structure
  finalize();

  this->printMsg("ContourTree computed", 1.0, timer.getElapsedTime());

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Nodes", std::to_string(getNumberOfNodes())},
    {"#Arcs", std::to_string(getNumberOfArcs())}});

  print();

  return 0;
}

int ContourTree::combineTrees() {

  Timer timer;

  if((!mergeTree_.getNumberOfNodes()) || (!splitTree_.getNumberOfNodes()))
    return -1;

  queue<const Node *> nodeQueue;
  const Node *mergeNode = NULL, *splitNode = NULL;
  int initNumber = 0;

  do {

    int initQueueSize = (int)nodeQueue.size();

    for(int i = 0; i < mergeTree_.getNumberOfNodes(); i++) {
      mergeNode = mergeTree_.getNode(i);
      if(isNodeEligible(mergeNode)) {
        nodeQueue.push(mergeNode);
      }
    }
    for(int i = 0; i < splitTree_.getNumberOfNodes(); i++) {
      splitNode = splitTree_.getNode(i);
      if(isNodeEligible(splitNode)) {
        nodeQueue.push(splitNode);
      }
    }

    // no more eligible nodes
    if((int)nodeQueue.size() == initQueueSize)
      break;

    const Node *n0 = NULL, *n1 = NULL, *other = NULL;

    do {

      n0 = nodeQueue.front();
      nodeQueue.pop();

#ifndef TTK_ENABLE_KAMIKAZE
      if(n0 == nullptr) {
        continue;
      }
#endif // TTK_ENABLE_KAMIKAZE

      if((mergeTree_.getNode(n0 - mergeTree_.getNode(0)) == n0)
         && (isNodeEligible(n0))) {
        // merge node
        // get the corresponding split node
        other = splitTree_.getVertexNode(n0->getVertexId());

        if(!((other->getNumberOfUpArcs())
             && (other->getNumberOfDownArcs() > 1))) {

          n1 = NULL;

          if(n0->getNumberOfUpArcs()) {

            n1 = mergeTree_.getNodeUpNeighbor(n0, 0);

            int newNode0 = makeNode(n0->getVertexId()),
                newNode1 = makeNode(n1->getVertexId());

            // "we move v and its incident edge from JT to the contour tree".
            makeArc(newNode0, newNode1);

            mergeTree_.clearArc(n0->getVertexId(), n1->getVertexId());
          }
          if((other->getNumberOfDownArcs() == 1)
             && (other->getNumberOfUpArcs() == 1)) {
            // "in ST, either v is a degree 2 node or..."
            // "in the former case, we suppress v in ST while maintaining the
            // connection"
            splitTree_.clearRegularNode(other->getVertexId());
          } else {
            // "or it's the root"
            // "in the latter, we delete v from ST"
            splitTree_.clearRoot(other->getVertexId());
          }

          // update n1 if it's eligible
          if((n1) && (isNodeEligible(n1))) {

            nodeQueue.push(n1);
          }
        } else {
          // unsure
          if(nodeQueue.empty())
            break;
          nodeQueue.push(n0);
        }
      }

      // symmetric case with nodes coming from the split tree
      else if((splitTree_.getNode(n0 - splitTree_.getNode(0)) == n0)
              && (isNodeEligible(n0))) {

        // merge node
        // get the corresponding split node
        other = mergeTree_.getVertexNode(n0->getVertexId());

        if(!((other->getNumberOfUpArcs())
             && (other->getNumberOfDownArcs() > 1))) {

          n1 = NULL;

          if(n0->getNumberOfUpArcs()) {

            n1 = splitTree_.getNodeUpNeighbor(n0, 0);

            int newNode0 = makeNode(n0->getVertexId()),
                newNode1 = makeNode(n1->getVertexId());

            makeArc(newNode1, newNode0);
            splitTree_.clearArc(n0->getVertexId(), n1->getVertexId());
          }
          if((other->getNumberOfDownArcs() == 1)
             && (other->getNumberOfUpArcs() == 1)) {
            // degree 2 node
            mergeTree_.clearRegularNode(other->getVertexId());
          } else {
            // other has to be the root
            mergeTree_.clearRoot(other->getVertexId());
          }

          // update n1 if it's eligible
          if((n1) && (isNodeEligible(n1))) {
            nodeQueue.push(n1);
          }
        } else {
          if(nodeQueue.empty())
            break;
          nodeQueue.push(n0);
        }
      }
    } while(nodeQueue.size());

    initNumber++;

    if((int)nodeList_.size() == vertexNumber_)
      break;

    // if one of the two trees became a line, we need to re-iterate the process
  } while((int)nodeList_.size() < vertexNumber_);

  if((int)nodeList_.size() != vertexNumber_) {
    this->printErr("Incomplete contour tree! ("
                   + std::to_string(nodeList_.size()) + " vs. "
                   + std::to_string(vertexNumber_) + ")");
  }

  this->printMsg(
    "MergeTree and SplitTree combined", 1.0, timer.getElapsedTime(), 1);

  return 0;
}

int ContourTree::finalize() {

  int minCount = 0, mergeCount = 0, splitCount = 0, maxCount = 0, regCount = 0;

  for(int i = 0; i < (int)nodeList_.size(); i++) {
    if(!nodeList_[i].getNumberOfDownArcs()) {
      minCount++;
    } else if(!nodeList_[i].getNumberOfUpArcs()) {
      maxCount++;
    } else if((nodeList_[i].getNumberOfDownArcs() == 1)
              && (nodeList_[i].getNumberOfUpArcs() == 1)) {
      regCount++;
    } else if(nodeList_[i].getNumberOfDownArcs() > 1) {
      mergeCount++;
    } else if(nodeList_[i].getNumberOfUpArcs() > 1) {
      splitCount++;
    }
  }

  vector<bool> inQueue(nodeList_.size(), false);
  queue<int> nodeIdQueue;

  for(int i = 0; i < (int)nodeList_.size(); i++) {
    if(!nodeList_[i].getNumberOfDownArcs()) {
      nodeIdQueue.push(i);
      inQueue[i] = true;
    }
  }

  while(nodeIdQueue.size()) {

    int nodeId = nodeIdQueue.front();
    nodeIdQueue.pop();

    for(int i = 0; i < nodeList_[nodeId].getNumberOfUpArcs(); i++) {
      int nextNodeId = finalizeSuperArc(nodeId, i);

      if(!inQueue[nextNodeId]) {
        nodeIdQueue.push(nextNodeId);
        inQueue[nextNodeId] = true;
      }
    }
  }

  return 0;
}

int ContourTree::finalizeSuperArc(const int &nodeId, const int &arcId) {

  if((nodeId < 0) || (nodeId >= (int)nodeList_.size()))
    return -1;
  if((arcId < 0) || (arcId >= nodeList_[nodeId].getNumberOfUpArcs()))
    return -2;

  int superArcId = openSuperArc(nodeId);

  int currentNodeId = nodeId;

  do {

    if(nodeList_[currentNodeId].getNumberOfUpArcs()) {
      if(currentNodeId != nodeId) {

        superArcList_[superArcId].appendRegularNode(currentNodeId);
        vertex2superArc_[nodeList_[currentNodeId].getVertexId()] = superArcId;
        vertex2superArcNode_[nodeList_[currentNodeId].getVertexId()]
          = superArcList_[superArcId].getNumberOfRegularNodes() - 1;

        currentNodeId
          = arcList_[nodeList_[currentNodeId].getUpArcId(0)].getUpNodeId();
      } else {
        currentNodeId
          = arcList_[nodeList_[currentNodeId].getUpArcId(arcId)].getUpNodeId();
      }
    }

  } while((currentNodeId != nodeId)
          && (nodeList_[currentNodeId].getNumberOfUpArcs() == 1)
          && (nodeList_[currentNodeId].getNumberOfDownArcs() == 1));

  if(currentNodeId != nodeId) {
    superArcList_[superArcId].setUpNodeId(currentNodeId);
    nodeList_[currentNodeId].addDownSuperArcId(superArcId);
  }

  return currentNodeId;
}

bool ContourTree::isNodeEligible(const Node *n) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(n == nullptr) {
    return false;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const Node *merge = NULL, *split = NULL;

  if(mergeTree_.getNode(n - mergeTree_.getNode(0)) == n) {
    merge = n;
    split = splitTree_.getVertexNode(merge->getVertexId());

    // "choose a non-root leaf that is not a split in ST"
    // non-root
    if((!merge->getNumberOfDownArcs())
       && (merge->getNumberOfUpArcs())
       // it is not a split in ST
       && (((split->getNumberOfDownArcs() < 2) && (split->getNumberOfUpArcs()))
           // or it's the root
           || ((split->getNumberOfDownArcs() == 1)
               && (!split->getNumberOfUpArcs())))) {

      return true;
    }
  }

  if(splitTree_.getNode(n - splitTree_.getNode(0)) == n) {
    split = n;
    merge = mergeTree_.getVertexNode(split->getVertexId());

    if((!split->getNumberOfDownArcs()) && (split->getNumberOfUpArcs())
       && (((merge->getNumberOfDownArcs() < 2) && (merge->getNumberOfUpArcs()))
           || ((merge->getNumberOfDownArcs() == 1)
               && (!merge->getNumberOfUpArcs())))) {

      return true;
    }
  }

  return false;
}

int ContourTree::computeSkeleton(unsigned int arcResolution) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { SubLevelSetTree::computeSkeleton(arcResolution); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { mergeTree_.computeSkeleton(arcResolution); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { splitTree_.computeSkeleton(arcResolution); }
  }

  return 0;
}

int ContourTree::smoothSkeleton(unsigned int skeletonSmoothing) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { SubLevelSetTree::smoothSkeleton(skeletonSmoothing); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { mergeTree_.smoothSkeleton(skeletonSmoothing); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { splitTree_.smoothSkeleton(skeletonSmoothing); }
  }

  return 0;
}

int ContourTree::clearSkeleton() {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { SubLevelSetTree::clearSkeleton(); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { mergeTree_.clearSkeleton(); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { splitTree_.clearSkeleton(); }
  }

  return 0;
}

int ContourTree::getPersistencePairs(
  vector<pair<pair<int, int>, double>> &pairs,
  vector<pair<pair<int, int>, double>> *mergePairs,
  vector<pair<pair<int, int>, double>> *splitPairs) const {
  if(pairs.size())
    return 0;

  vector<pair<pair<int, int>, double>> defaultMergePairs{};
  if(!mergePairs) {
    mergePairs = &defaultMergePairs;
  }
  vector<pair<pair<int, int>, double>> defaultSplitPairs{};
  if(!splitPairs) {
    splitPairs = &defaultSplitPairs;
  }

  if(!mergePairs->size() || !splitPairs->size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
    {

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      {
        if(!mergePairs->size())
          mergeTree_.getPersistencePairs(*mergePairs);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
      {
        if(!splitPairs->size())
          splitTree_.getPersistencePairs(*splitPairs);
      }
    }
  }

  pairs.resize(mergePairs->size() + splitPairs->size());

  for(unsigned int i = 0; i < mergePairs->size(); ++i)
    pairs[i] = (*mergePairs)[i];

  unsigned int shift = mergePairs->size();
  for(unsigned int i = 0; i < splitPairs->size(); ++i)
    pairs[shift + i] = (*splitPairs)[i];

  std::sort(pairs.begin(), pairs.end(), _pCmp);

  return 0;
}

int ContourTree::getPersistencePlot(
  vector<pair<double, int>> &plot,
  vector<pair<pair<int, int>, double>> *mergePairs,
  vector<pair<pair<int, int>, double>> *splitPairs,
  vector<pair<pair<int, int>, double>> *pairs) const {

  vector<pair<pair<int, int>, double>> defaultPairs{};
  if(!pairs) {
    pairs = &defaultPairs;
  }

  if(!pairs->size())
    getPersistencePairs(*pairs, mergePairs, splitPairs);

  plot.resize(pairs->size());

  for(int i = 0; i < (int)plot.size(); i++) {
    plot[i].first = (*pairs)[i].second;
    if(plot[i].first < Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS)) {
      plot[i].first = Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS);
    }
    plot[i].second = pairs->size() - i;
  }

  return 0;
}

int ContourTree::getPersistenceDiagram(
  vector<pair<double, double>> &diagram,
  vector<pair<pair<int, int>, double>> *mergePairs,
  vector<pair<pair<int, int>, double>> *splitPairs,
  vector<pair<pair<int, int>, double>> *pairs) const {

  vector<pair<pair<int, int>, double>> defaultPairs{};
  if(!pairs) {
    pairs = &defaultPairs;
  }

  if(!pairs->size())
    getPersistencePairs(*pairs, mergePairs, splitPairs);

  // fast fix :(
  diagram.resize(pairs->size());

  for(int i = 0; i < (int)pairs->size(); i++) {
    if((maximumList_) && (!minimumList_)) {
      // split tree
      diagram[i].second = (*vertexScalars_)[(*pairs)[i].first.first];
      diagram[i].first = (*vertexScalars_)[(*pairs)[i].first.second];
    } else {
      // join tree or contour tree
      diagram[i].first = (*vertexScalars_)[(*pairs)[i].first.first];
      diagram[i].second = (*vertexScalars_)[(*pairs)[i].first.second];
    }
  }

  std::sort(diagram.begin(), diagram.end(), _pPairCmp);

  return 0;
}

int ContourTree::simplify(const double &simplificationThreshold,
                          ContourTreeSimplificationMetric *metric) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { mergeTree_.simplify(simplificationThreshold, metric); }

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    { splitTree_.simplify(simplificationThreshold, metric); }
  }

  return 0;
}
