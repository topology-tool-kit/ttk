#include <ManifoldCheck.h>

using namespace std;
using namespace ttk;

ManifoldCheck::ManifoldCheck() {

  triangulation_ = NULL;
  vertexLinkComponentNumber_ = NULL;
  edgeLinkComponentNumber_ = NULL;
  triangleLinkComponentNumber_ = NULL;
}

ManifoldCheck::~ManifoldCheck() {
}

int ManifoldCheck::vertexManifoldCheck(const SimplexId &vertexId) const {

  SimplexId linkSize = triangulation_->getVertexLinkNumber(vertexId);

  if(triangulation_->getDimensionality() == 1)
    return linkSize;

  vector<SimplexId> linkNeighbors;

  for(SimplexId i = 0; i < linkSize; i++) {
    SimplexId linkId = -1;
    triangulation_->getVertexLink(vertexId, i, linkId);

    bool isIn = false;
    SimplexId neighborId = -1;

    if(triangulation_->getDimensionality() == 2) {
      triangulation_->getEdgeVertex(linkId, 0, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation_->getEdgeVertex(linkId, 1, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);
    }
    if(triangulation_->getDimensionality() == 3) {
      triangulation_->getTriangleVertex(linkId, 0, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation_->getTriangleVertex(linkId, 1, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);

      triangulation_->getTriangleVertex(linkId, 2, neighborId);
      isIn = false;
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId) {
          isIn = true;
          break;
        }
      }
      if(!isIn)
        linkNeighbors.push_back(neighborId);
    }
  }

  vector<UnionFind> seeds(linkNeighbors.size());
  vector<UnionFind *> seedList(linkNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)seeds.size(); i++) {
    seedList[i] = &(seeds[i]);
  }

  for(SimplexId i = 0; i < (SimplexId)linkSize; i++) {

    SimplexId linkId = -1;
    triangulation_->getVertexLink(vertexId, i, linkId);

    SimplexId neighborId0 = -1, neighborId1 = -1, neighborId2 = -1;
    SimplexId uf0 = -1, uf1 = -1, uf2 = -1;

    if(triangulation_->getDimensionality() == 2) {
      triangulation_->getEdgeVertex(linkId, 0, neighborId0);
      triangulation_->getEdgeVertex(linkId, 1, neighborId1);

      // connect the two uf together
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId0) {
          uf0 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId1) {
          uf1 = j;
          break;
        }
      }

      seedList[uf0] = makeUnion(seedList[uf0], seedList[uf1]);
      seedList[uf1] = seedList[uf0];
    }

    if(triangulation_->getDimensionality() == 3) {
      triangulation_->getTriangleVertex(linkId, 0, neighborId0);
      triangulation_->getTriangleVertex(linkId, 1, neighborId1);
      triangulation_->getTriangleVertex(linkId, 2, neighborId2);

      // connect the two uf together
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId0) {
          uf0 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId1) {
          uf1 = j;
          break;
        }
      }
      for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
        if(linkNeighbors[j] == neighborId2) {
          uf2 = j;
          break;
        }
      }

      seedList[uf0] = makeUnion(seedList[uf0], seedList[uf1]);
      seedList[uf0] = makeUnion(seedList[uf0], seedList[uf2]);
      seedList[uf1] = seedList[uf0];
      seedList[uf2] = seedList[uf0];
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)seedList.size(); i++) {
    seedList[i] = seedList[i]->find();
  }

  vector<UnionFind *>::iterator it;
  sort(seedList.begin(), seedList.end());
  it = unique(seedList.begin(), seedList.end());
  seedList.resize(distance(seedList.begin(), it));

  return (SimplexId)seedList.size();
}

int ManifoldCheck::edgeManifoldCheck(const SimplexId &edgeId) const {

  SimplexId linkSize = triangulation_->getEdgeLinkNumber(edgeId);

  if(triangulation_->getDimensionality() == 2)
    return linkSize;

  vector<SimplexId> linkNeighbors;

  for(SimplexId i = 0; i < linkSize; i++) {
    SimplexId linkId = -1;
    triangulation_->getEdgeLink(edgeId, i, linkId);

    bool isIn = false;
    SimplexId neighborId = -1;

    triangulation_->getEdgeVertex(linkId, 0, neighborId);
    isIn = false;
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId) {
        isIn = true;
        break;
      }
    }
    if(!isIn)
      linkNeighbors.push_back(neighborId);

    triangulation_->getEdgeVertex(linkId, 1, neighborId);
    isIn = false;
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId) {
        isIn = true;
        break;
      }
    }
    if(!isIn)
      linkNeighbors.push_back(neighborId);
  }

  vector<UnionFind> seeds(linkNeighbors.size());
  vector<UnionFind *> seedList(linkNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)seeds.size(); i++) {
    seedList[i] = &(seeds[i]);
  }

  for(SimplexId i = 0; i < (SimplexId)linkSize; i++) {

    SimplexId linkId = -1;
    triangulation_->getEdgeLink(edgeId, i, linkId);

    SimplexId neighborId0 = -1, neighborId1 = -1;
    SimplexId uf0 = -1, uf1 = -1;

    triangulation_->getEdgeVertex(linkId, 0, neighborId0);
    triangulation_->getEdgeVertex(linkId, 1, neighborId1);

    // connect the two uf together
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId0) {
        uf0 = j;
        break;
      }
    }
    for(SimplexId j = 0; j < (SimplexId)linkNeighbors.size(); j++) {
      if(linkNeighbors[j] == neighborId1) {
        uf1 = j;
        break;
      }
    }

    seedList[uf0] = makeUnion(seedList[uf0], seedList[uf1]);
    seedList[uf1] = seedList[uf0];
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)seedList.size(); i++) {
    seedList[i] = seedList[i]->find();
  }

  vector<UnionFind *>::iterator it;
  sort(seedList.begin(), seedList.end());
  it = unique(seedList.begin(), seedList.end());
  seedList.resize(distance(seedList.begin(), it));

  return (SimplexId)seedList.size();
}

int ManifoldCheck::execute() const {

  Timer t;

// check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_)
    return -1;
#endif

  SimplexId vertexNumber = triangulation_->getNumberOfVertices();

  if(vertexLinkComponentNumber_) {

    vertexLinkComponentNumber_->resize(vertexNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < vertexNumber; i++) {
      (*vertexLinkComponentNumber_)[i] = vertexManifoldCheck(i);
    }
  }

  if((edgeLinkComponentNumber_) && (triangulation_->getDimensionality() >= 2)) {

    SimplexId edgeNumber = triangulation_->getNumberOfEdges();
    edgeLinkComponentNumber_->resize(edgeNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < edgeNumber; i++) {
      (*edgeLinkComponentNumber_)[i] = edgeManifoldCheck(i);
    }
  }

  if((triangleLinkComponentNumber_)
     && (triangulation_->getDimensionality() == 3)) {

    SimplexId triangleNumber = triangulation_->getNumberOfTriangles();
    triangleLinkComponentNumber_->resize(triangleNumber);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < triangleNumber; i++) {
      (*triangleLinkComponentNumber_)[i]
        = triangulation_->getTriangleLinkNumber(i);
    }
  }

  {
    stringstream msg;
    msg << "[ManifoldCheck] Data-set (" << vertexNumber
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}
