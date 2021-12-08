#include <CompactTriangulation.h>

using namespace ttk;

CompactTriangulation::CompactTriangulation() {
  setDebugMsgPrefix("CompactTriangulation");
  clear();
  caches_.resize(threadNumber_);
  cacheMaps_.resize(threadNumber_);
}

CompactTriangulation::CompactTriangulation(const CompactTriangulation &rhs)
  : AbstractTriangulation(rhs), doublePrecision_(rhs.doublePrecision_),
    maxCellDim_(rhs.maxCellDim_), cellNumber_(rhs.cellNumber_),
    vertexNumber_(rhs.vertexNumber_), nodeNumber_(rhs.nodeNumber_),
    pointSet_(rhs.pointSet_), vertexIndices_(rhs.vertexIndices_),
    vertexIntervals_(rhs.vertexIntervals_), edgeIntervals_(rhs.edgeIntervals_),
    triangleIntervals_(rhs.triangleIntervals_),
    cellIntervals_(rhs.cellIntervals_), cellArray_(rhs.cellArray_),
    externalCells_(rhs.externalCells_) {
}

CompactTriangulation &
  CompactTriangulation::operator=(const CompactTriangulation &rhs) {
  if(this != &rhs) {
    doublePrecision_ = rhs.doublePrecision_;
    maxCellDim_ = rhs.maxCellDim_;
    cellNumber_ = rhs.cellNumber_;
    vertexNumber_ = rhs.vertexNumber_;
    nodeNumber_ = rhs.nodeNumber_;
    pointSet_ = rhs.pointSet_;
    vertexIndices_ = rhs.vertexIndices_;
    vertexIntervals_ = rhs.vertexIntervals_;
    edgeIntervals_ = rhs.edgeIntervals_;
    triangleIntervals_ = rhs.triangleIntervals_;
    cellIntervals_ = rhs.cellIntervals_;
    cellArray_ = rhs.cellArray_;
    externalCells_ = rhs.externalCells_;
    // cache system is not copied
  }
  return *this;
}

CompactTriangulation::~CompactTriangulation() {
}

int CompactTriangulation::reorderVertices(std::vector<SimplexId> &vertexMap) {
  // get the number of nodes (the max value in the array)
  nodeNumber_ = 0;
  for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
    if(vertexIndices_[vid] > nodeNumber_) {
      nodeNumber_ = vertexIndices_[vid];
    }
  }
  nodeNumber_++; // since the index starts from 0
  std::vector<std::vector<SimplexId>> nodeVertices(nodeNumber_);
  for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
    nodeVertices[vertexIndices_[vid]].push_back(vid);
  }

  // update the vertex intervals
  vertexIntervals_.resize(nodeNumber_ + 1);
  vertexIntervals_[0] = -1;
  SimplexId vertexCount = 0;
  for(SimplexId nid = 0; nid < nodeNumber_; nid++) {
    for(SimplexId vid : nodeVertices[nid]) {
      vertexMap[vid] = vertexCount++;
    }
    vertexIntervals_[nid + 1] = vertexCount - 1;
  }

  // rearange the vertex coordinate values
  if(doublePrecision_) {
    std::vector<double> newPointSet(3 * vertexNumber_);
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        newPointSet[3 * vertexMap[vid] + j]
          = ((double *)pointSet_)[3 * vid + j];
      }
    }
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        ((double *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
      }
    }
  } else {
    std::vector<float> newPointSet(3 * vertexNumber_);
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        newPointSet[3 * vertexMap[vid] + j] = ((float *)pointSet_)[3 * vid + j];
      }
    }
    for(SimplexId vid = 0; vid < vertexNumber_; vid++) {
      for(int j = 0; j < 3; j++) {
        ((float *)pointSet_)[3 * vid + j] = newPointSet[3 * vid + j];
      }
    }
  }

  // change the vertex indices
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId vid = vertexIntervals_[nid - 1] + 1;
        vid <= vertexIntervals_[nid]; vid++) {
      ((int *)vertexIndices_)[vid] = nid;
    }
  }

  return 0;
}

int CompactTriangulation::reorderCells(const std::vector<SimplexId> &vertexMap,
                                       const SimplexId &cellNumber,
                                       const LongSimplexId *connectivity,
                                       const LongSimplexId *offset) {
  // change the indices in cell array
  SimplexId cellCount = 0, verticesPerCell = offset[1] - offset[0];
  std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
  LongSimplexId *cellArr = ((LongSimplexId *)connectivity);

  for(SimplexId cid = 0; cid < cellNumber; cid++) {
    SimplexId cellId = offset[cid];
    for(int j = 0; j < verticesPerCell; j++) {
      cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
    }
    std::sort(cellArr + cellId, cellArr + cellId + verticesPerCell);
    nodeCells[vertexIndices_[cellArr[cellId]]].push_back(cid);
  }

  // rearange the cell array
  cellIntervals_.resize(nodeNumber_ + 1);
  externalCells_.resize(nodeNumber_ + 1);
  cellIntervals_[0] = -1;
  std::vector<LongSimplexId> newCellArray(verticesPerCell * cellNumber);
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId cid : nodeCells[nid]) {
      SimplexId cellId = verticesPerCell * cid;
      SimplexId newCellId = verticesPerCell * cellCount;
      for(int j = 0; j < verticesPerCell; j++) {
        newCellArray[newCellId + j] = connectivity[cellId + j];
        if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
          SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
          if(externalCells_[nodeNum].empty()
             || externalCells_[nodeNum].back() != cid) {
            externalCells_[nodeNum].push_back(cid);
          }
        }
      }
      cellCount++;
    }
    cellIntervals_[nid] = cellCount - 1;
  }

  // copy the new cell array back to original one
  for(SimplexId i = 0; i < verticesPerCell * cellNumber; i++) {
    ((LongSimplexId *)connectivity)[i] = newCellArray[i];
  }

  return 0;
}

int CompactTriangulation::reorderCells(const std::vector<SimplexId> &vertexMap,
                                       const LongSimplexId *cellArray) {
  // change the indices in cell array
  SimplexId cellCount = 0, verticesPerCell = cellArray[0];
  std::vector<std::vector<SimplexId>> nodeCells(nodeNumber_ + 1);
  LongSimplexId *cellArr = ((LongSimplexId *)cellArray);

  for(SimplexId cid = 0; cid < cellNumber_; cid++) {
    SimplexId cellId = (verticesPerCell + 1) * cid;
    for(int j = 1; j <= verticesPerCell; j++) {
      cellArr[cellId + j] = vertexMap[cellArr[cellId + j]];
    }
    std::sort(cellArr + cellId + 1, cellArr + cellId + 1 + verticesPerCell);
    nodeCells[vertexIndices_[cellArr[cellId + 1]]].push_back(cid);
  }

  // rearange the cell array
  cellIntervals_.resize(nodeNumber_ + 1);
  externalCells_.resize(nodeNumber_ + 1);
  cellIntervals_[0] = -1;
  std::vector<LongSimplexId> newCellArray((verticesPerCell + 1) * cellNumber_);
  for(SimplexId nid = 1; nid <= nodeNumber_; nid++) {
    for(SimplexId cid : nodeCells[nid]) {
      SimplexId cellId = (verticesPerCell + 1) * cid;
      SimplexId newCellId = (verticesPerCell + 1) * cellCount;
      newCellArray[newCellId] = verticesPerCell;
      for(int j = 1; j <= verticesPerCell; j++) {
        newCellArray[newCellId + j] = cellArray[cellId + j];
        if(newCellArray[newCellId + j] > vertexIntervals_[nid]) {
          SimplexId nodeNum = vertexIndices_[newCellArray[newCellId + j]];
          if(externalCells_[nodeNum].empty()
             || externalCells_[nodeNum].back() != cid) {
            externalCells_[nodeNum].push_back(cid);
          }
        }
      }
      cellCount++;
    }
    cellIntervals_[nid] = cellCount - 1;
  }

  // copy the new cell array back to original one
  for(SimplexId i = 0; i < (verticesPerCell + 1) * cellNumber_; i++) {
    ((LongSimplexId *)cellArray)[i] = newCellArray[i];
  }

  return 0;
}

int CompactTriangulation::buildInternalEdgeMap(
  ImplicitCluster *const nodePtr,
  bool computeInternalEdgeList,
  bool computeInternalEdgeMap) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId edgeCount = 0, verticesPerCell = cellArray_->getCellVertexNumber(0);

  if(!nodePtr->internalEdgeMap_.empty()) {
    // if the edge map has been computed and only request the edge list
    if(computeInternalEdgeList) {
      nodePtr->internalEdgeList_ = std::vector<std::array<SimplexId, 2>>(
        nodePtr->internalEdgeMap_.size());
      for(auto iter = nodePtr->internalEdgeMap_.begin();
          iter != nodePtr->internalEdgeMap_.end(); iter++) {
        nodePtr->internalEdgeList_[iter->second - 1] = iter->first;
      }
      return 0;
    }
  }

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 2> edgeIds;

    // loop through each edge of the cell
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      // the edge does not belong to the current node
      if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[1] = cellArray_->getCellVertex(cid, k);

        // not found in the edge map - assign new edge id
        if(nodePtr->internalEdgeMap_.find(edgeIds)
           == nodePtr->internalEdgeMap_.end()) {
          edgeCount++;
          nodePtr->internalEdgeMap_[edgeIds] = edgeCount;
        }
      }
    }
  }

  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 2> edgeIds;

    // loop through each edge of the cell
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        edgeIds[1] = cellArray_->getCellVertex(cid, k);

        // the edge is in the current node
        if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
           && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
          if(nodePtr->internalEdgeMap_.find(edgeIds)
             == nodePtr->internalEdgeMap_.end()) {
            edgeCount++;
            (nodePtr->internalEdgeMap_)[edgeIds] = edgeCount;
          }
        }
      }
    }
  }

  if(computeInternalEdgeList) {
    nodePtr->internalEdgeList_
      = std::vector<std::array<SimplexId, 2>>(nodePtr->internalEdgeMap_.size());
    for(auto iter = nodePtr->internalEdgeMap_.begin();
        iter != nodePtr->internalEdgeMap_.end(); iter++) {
      nodePtr->internalEdgeList_[iter->second - 1] = iter->first;
    }
  }

  if(!computeInternalEdgeMap) {
    nodePtr->internalEdgeMap_.clear();
  }

  return 0;
}

int CompactTriangulation::buildExternalEdgeMap(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  boost::unordered_map<SimplexId, std::vector<std::array<SimplexId, 2>>>
    edgeNodes;

  // loop through the external cell list
  for(size_t i = 0; i < externalCells_[nodePtr->nid].size(); i++) {
    std::array<SimplexId, 2> edgeIds;
    SimplexId cellId = externalCells_[nodePtr->nid][i];

    // loop through each edge of the cell
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cellId, j);
        edgeIds[1] = cellArray_->getCellVertex(cellId, k);

        // check if the edge is an external edge
        if(edgeIds[0] <= vertexIntervals_[nodePtr->nid - 1]
           && edgeIds[1] > vertexIntervals_[nodePtr->nid - 1]
           && edgeIds[1] <= vertexIntervals_[nodePtr->nid]) {
          SimplexId nid = vertexIndices_[edgeIds[0]];
          edgeNodes[nid].push_back(edgeIds);
        }
      }
    }
  }

  boost::unordered_map<SimplexId,
                       std::vector<std::array<SimplexId, 2>>>::iterator iter;
  for(iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
    ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
    if(exnode) {
      if(exnode->internalEdgeMap_.empty())
        buildInternalEdgeMap(exnode, false, true);
      for(std::array<SimplexId, 2> edgePair : iter->second) {
        (nodePtr->externalEdgeMap_)[edgePair]
          = exnode->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[iter->first - 1];
      }
    } else {
      ImplicitCluster tmpCluster(iter->first);
      buildInternalEdgeMap(&tmpCluster, false, true);
      for(std::array<SimplexId, 2> edgePair : iter->second) {
        (nodePtr->externalEdgeMap_)[edgePair]
          = tmpCluster.internalEdgeMap_.at(edgePair)
            + edgeIntervals_[iter->first - 1];
      }
    }
  }

  return 0;
}

int CompactTriangulation::buildInternalTriangleMap(
  ImplicitCluster *const nodePtr,
  bool computeInternalTriangleList,
  bool computeInternalTriangleMap) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId triangleCount = 0, verticesPerCell = 4;
  std::vector<std::vector<SimplexId>> localTriangleStars;

  if(!nodePtr->internalTriangleMap_.empty()) {
    if(computeInternalTriangleList) {
      nodePtr->internalTriangleList_ = std::vector<std::array<SimplexId, 3>>(
        nodePtr->internalTriangleMap_.size());
      for(auto iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        (nodePtr->internalTriangleList_)[iter->second - 1] = iter->first;
      }
      return 0;
    }
  }

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 3> triangleIds;

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      // the triangle does not belong to the current node
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
        for(SimplexId l = k + 1; l < verticesPerCell; l++) {
          triangleIds[1] = cellArray_->getCellVertex(cid, k);
          triangleIds[2] = cellArray_->getCellVertex(cid, l);

          if(nodePtr->internalTriangleMap_.find(triangleIds)
             == nodePtr->internalTriangleMap_.end()) {
            triangleCount++;
            (nodePtr->internalTriangleMap_)[triangleIds] = triangleCount;
          }
        }
      }
    }
  }

  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 3> triangleIds;

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
         && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
        for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);

            if(nodePtr->internalTriangleMap_.find(triangleIds)
               == nodePtr->internalTriangleMap_.end()) {
              triangleCount++;
              (nodePtr->internalTriangleMap_)[triangleIds] = triangleCount;
            }
          }
        }
      }
    }
  }

  if(computeInternalTriangleList) {
    nodePtr->internalTriangleList_ = std::vector<std::array<SimplexId, 3>>(
      nodePtr->internalTriangleMap_.size());
    for(auto iter = nodePtr->internalTriangleMap_.begin();
        iter != nodePtr->internalTriangleMap_.end(); iter++) {
      (nodePtr->internalTriangleList_)[iter->second - 1] = iter->first;
    }
  }

  if(!computeInternalTriangleMap) {
    nodePtr->internalTriangleMap_.clear();
  }

  return 0;
}

int CompactTriangulation::buildExternalTriangleMap(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  boost::unordered_map<SimplexId, std::vector<std::array<SimplexId, 3>>>
    nodeTriangles;

  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 3> triangleIds;

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] <= vertexIntervals_[nodePtr->nid - 1]) {
        for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);

            if(triangleIds[1] > vertexIntervals_[nodePtr->nid - 1]
               && triangleIds[1] <= vertexIntervals_[nodePtr->nid]) {
              SimplexId nodeNum = vertexIndices_[triangleIds[0]];
              nodeTriangles[nodeNum].push_back(triangleIds);
            } else if(triangleIds[2] > vertexIntervals_[nodePtr->nid - 1]
                      && triangleIds[2] <= vertexIntervals_[nodePtr->nid]) {
              SimplexId nodeNum = vertexIndices_[triangleIds[0]];
              nodeTriangles[nodeNum].push_back(triangleIds);
            }
          }
        }
      }
    }
  }

  boost::unordered_map<SimplexId,
                       std::vector<std::array<SimplexId, 3>>>::iterator iter;
  for(iter = nodeTriangles.begin(); iter != nodeTriangles.end(); iter++) {
    ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
    if(exnode) {
      if(exnode->internalTriangleMap_.empty())
        buildInternalTriangleMap(exnode, false, true);
      for(std::array<SimplexId, 3> triangleVec : iter->second) {
        (nodePtr->externalTriangleMap_)[triangleVec]
          = exnode->internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[iter->first - 1];
      }
    } else {
      ImplicitCluster tmpCluster(iter->first);
      buildInternalTriangleMap(&tmpCluster, false, true);
      for(std::array<SimplexId, 3> triangleVec : iter->second) {
        (nodePtr->externalTriangleMap_)[triangleVec]
          = tmpCluster.internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[iter->first - 1];
      }
    }
  }

  return 0;
}

SimplexId CompactTriangulation::countInternalEdges(SimplexId nodeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodeId <= 0 || nodeId > nodeNumber_)
    return -1;
#endif

  SimplexId edgeCount = 0, verticesPerCell = cellArray_->getCellVertexNumber(0);
  boost::unordered_set<std::array<SimplexId, 2>,
                       boost::hash<std::array<SimplexId, 2>>>
    edgeSet;

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
      cid <= cellIntervals_[nodeId]; cid++) {
    std::array<SimplexId, 2> edgeIds{};

    // loop through each edge of the cell
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      // the edge does not belong to the current node
      if(edgeIds[0] > vertexIntervals_[nodeId]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[1] = cellArray_->getCellVertex(cid, k);

        // not found in the edge map - assign new edge id
        if(edgeSet.find(edgeIds) == edgeSet.end()) {
          edgeCount++;
          edgeSet.insert(edgeIds);
        }
      }
    }
  }

  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodeId]) {
    std::array<SimplexId, 2> edgeIds{};

    // loop through each edge of the cell
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        edgeIds[1] = cellArray_->getCellVertex(cid, k);

        // the edge is in the current node
        if(edgeIds[0] > vertexIntervals_[nodeId - 1]
           && edgeIds[0] <= vertexIntervals_[nodeId]) {
          if(edgeSet.find(edgeIds) == edgeSet.end()) {
            edgeCount++;
            edgeSet.insert(edgeIds);
          }
        }
      }
    }
  }

  return edgeCount;
}

int CompactTriangulation::countInternalTriangles(SimplexId nodeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodeId <= 0 || nodeId > nodeNumber_)
    return -1;
#endif

  SimplexId triangleCount = 0,
            verticesPerCell = cellArray_->getCellVertexNumber(0);
  boost::unordered_set<std::array<SimplexId, 3>> triangleSet;

  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodeId - 1] + 1;
      cid <= cellIntervals_[nodeId]; cid++) {
    std::array<SimplexId, 3> triangleIds{};

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      // the triangle does not belong to the current node
      if(triangleIds[0] > vertexIntervals_[nodeId]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
        for(SimplexId l = k + 1; l < verticesPerCell; l++) {
          triangleIds[1] = cellArray_->getCellVertex(cid, k);
          triangleIds[2] = cellArray_->getCellVertex(cid, l);

          if(triangleSet.find(triangleIds) == triangleSet.end()) {
            triangleCount++;
            triangleSet.insert(triangleIds);
          }
        }
      }
    }
  }

  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodeId]) {
    std::array<SimplexId, 3> triangleIds{};

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] > vertexIntervals_[nodeId - 1]
         && triangleIds[0] <= vertexIntervals_[nodeId]) {
        for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);

            if(triangleSet.find(triangleIds) == triangleSet.end()) {
              triangleCount++;
              triangleSet.insert(triangleIds);
            }
          }
        }
      }
    }
  }

  return triangleCount;
}

int CompactTriangulation::getClusterCellNeighbors(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  std::vector<std::vector<SimplexId>> localCellNeighbors(
    cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);
  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  std::vector<std::vector<SimplexId>> localVertexStars;

  if(nodePtr->vertexStars_.empty()) {
    getClusterVertexStars(nodePtr);
  }
  // sort the vertex star vector
  nodePtr->vertexStars_.copyTo(localVertexStars);
  for(size_t i = 0; i < localVertexStars.size(); i++) {
    sort(localVertexStars[i].begin(), localVertexStars[i].end());
  }

  boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>> nodeMaps;
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      if(cellArray_->getCellVertex(cid, j) > vertexIntervals_[nodePtr->nid]) {
        SimplexId nodeId = vertexIndices_[cellArray_->getCellVertex(cid, j)];
        if(nodeMaps.find(nodeId) == nodeMaps.end()) {
          ImplicitCluster newNode(nodeId);
          getClusterVertexStars(&newNode);
          // sort the vertex star list
          std::vector<std::vector<SimplexId>> tmpVec;
          newNode.vertexStars_.copyTo(tmpVec);
          for(size_t i = 0; i < tmpVec.size(); i++) {
            sort(tmpVec[i].begin(), tmpVec[i].end());
          }
          nodeMaps[nodeId] = tmpVec;
        }
      }
    }
  }

  if(getDimensionality() == 2) {
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      for(SimplexId j = 0; j < 3; j++) {

        SimplexId v0 = cellArray_->getCellVertex(cid, j);
        SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 3);
        SimplexId localV0 = v0 - vertexIntervals_[nodePtr->nid - 1] - 1;
        SimplexId localV1 = v1 - vertexIntervals_[nodePtr->nid - 1] - 1;

        std::vector<SimplexId> stars0, stars1;
        if(v0 <= vertexIntervals_[nodePtr->nid]) {
          stars0 = localVertexStars[localV0];
        } else {
          localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
          stars0 = nodeMaps[vertexIndices_[v0]][localV0];
        }
        if(v1 <= vertexIntervals_[nodePtr->nid]) {
          stars1 = localVertexStars[localV1];
        } else {
          localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
          stars1 = nodeMaps[vertexIndices_[v1]][localV1];
        }

        // perform an intersection of the 2 sorted star lists
        SimplexId pos0 = 0, pos1 = 0;
        SimplexId intersection = -1;

        while((pos0 < (SimplexId)stars0.size())
              && (pos1 < (SimplexId)stars1.size())) {
          SimplexId biggest = stars0[pos0];
          if(stars1[pos1] > biggest) {
            biggest = stars1[pos1];
          }

          for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
            if(stars0[l] < biggest) {
              pos0++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
            if(stars1[l] < biggest) {
              pos1++;
            } else {
              break;
            }
          }

          if(pos0 >= (SimplexId)stars0.size()
             || pos1 >= (SimplexId)stars1.size())
            break;

          if(stars0[pos0] == stars1[pos1]) {
            if(stars0[pos0] != cid) {
              intersection = stars0[pos0];
              break;
            }
            pos0++;
            pos1++;
          }
        }

        if(intersection != -1) {
          localCellNeighbors[cid - cellIntervals_[nodePtr->nid - 1] - 1]
            .push_back(intersection);
        }
      }
    }
  }

  else if(getDimensionality() == 3) {
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      // go triangle by triangle
      for(SimplexId j = 0; j < 4; j++) {

        SimplexId v0 = cellArray_->getCellVertex(cid, j % 4);
        SimplexId v1 = cellArray_->getCellVertex(cid, (j + 1) % 4);
        SimplexId v2 = cellArray_->getCellVertex(cid, (j + 2) % 4);

        SimplexId localV0 = v0 - vertexIntervals_[nodePtr->nid - 1] - 1;
        SimplexId localV1 = v1 - vertexIntervals_[nodePtr->nid - 1] - 1;
        SimplexId localV2 = v2 - vertexIntervals_[nodePtr->nid - 1] - 1;

        std::vector<SimplexId> stars0, stars1, stars2;
        if(v0 <= vertexIntervals_[nodePtr->nid]) {
          stars0 = localVertexStars[localV0];
        } else {
          localV0 = v0 - vertexIntervals_[vertexIndices_[v0] - 1] - 1;
          stars0 = nodeMaps[vertexIndices_[v0]][localV0];
        }
        if(v1 <= vertexIntervals_[nodePtr->nid]) {
          stars1 = localVertexStars[localV1];
        } else {
          localV1 = v1 - vertexIntervals_[vertexIndices_[v1] - 1] - 1;
          stars1 = nodeMaps[vertexIndices_[v1]][localV1];
        }
        if(v2 <= vertexIntervals_[nodePtr->nid]) {
          stars2 = localVertexStars[localV2];
        } else {
          localV2 = v2 - vertexIntervals_[vertexIndices_[v2] - 1] - 1;
          stars2 = nodeMaps[vertexIndices_[v2]][localV2];
        }

        // perform an intersection of the 3 (sorted) star lists
        SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
        SimplexId intersection = -1;

        while((pos0 < (SimplexId)stars0.size())
              && (pos1 < (SimplexId)stars1.size())
              && (pos2 < (SimplexId)stars2.size())) {

          SimplexId biggest = stars0[pos0];
          if(stars1[pos1] > biggest) {
            biggest = stars1[pos1];
          }
          if(stars2[pos2] > biggest) {
            biggest = stars2[pos2];
          }

          for(SimplexId l = pos0; l < (SimplexId)stars0.size(); l++) {
            if(stars0[l] < biggest) {
              pos0++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos1; l < (SimplexId)stars1.size(); l++) {
            if(stars1[l] < biggest) {
              pos1++;
            } else {
              break;
            }
          }
          for(SimplexId l = pos2; l < (SimplexId)stars2.size(); l++) {
            if(stars2[l] < biggest) {
              pos2++;
            } else {
              break;
            }
          }

          if(pos0 >= (SimplexId)stars0.size()
             || pos1 >= (SimplexId)stars1.size()
             || pos2 >= (SimplexId)stars2.size())
            break;

          if((stars0[pos0] == stars1[pos1]) && (stars0[pos0] == stars2[pos2])) {
            if(stars0[pos0] != cid) {
              intersection = stars0[pos0];
              break;
            }
            pos0++;
            pos1++;
            pos2++;
          }
        }

        if(intersection != -1) {
          localCellNeighbors[cid - cellIntervals_[nodePtr->nid - 1] - 1]
            .push_back(intersection);
        }
      }
    }
  }

  nodePtr->cellNeighbors_.fillFrom(localCellNeighbors);
  return 0;
}

int CompactTriangulation::getClusterCellTriangles(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = 4;
  boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
    nodeTriangles;

  nodePtr->tetraTriangles_ = std::vector<std::array<SimplexId, 4>>(
    cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);

  if(nodePtr->internalTriangleMap_.empty()) {
    buildInternalTriangleMap(nodePtr, false, true);
  }

  for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
      i <= cellIntervals_[nodePtr->nid]; i++) {
    std::array<SimplexId, 3> triangleVec;
    // get the internal triangle from the map
    triangleVec[0] = cellArray_->getCellVertex(i, 0);
    for(SimplexId k = 1; k < verticesPerCell - 1; k++) {
      triangleVec[1] = cellArray_->getCellVertex(i, k);
      for(SimplexId l = k + 1; l < verticesPerCell; l++) {
        triangleVec[2] = cellArray_->getCellVertex(i, l);
        (nodePtr->tetraTriangles_)[i - cellIntervals_[nodePtr->nid - 1] - 1]
                                  [k + l - 3]
          = nodePtr->internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[nodePtr->nid - 1];
      }
    }
    // group the external triangles by node id
    triangleVec[0] = cellArray_->getCellVertex(i, 1);
    triangleVec[1] = cellArray_->getCellVertex(i, 2);
    triangleVec[2] = cellArray_->getCellVertex(i, 3);
    if(triangleVec[0] <= vertexIntervals_[nodePtr->nid]) {
      (nodePtr->tetraTriangles_)[i - cellIntervals_[nodePtr->nid - 1] - 1]
        .back()
        = nodePtr->internalTriangleMap_.at(triangleVec)
          + triangleIntervals_[nodePtr->nid - 1];
    } else {
      std::vector<SimplexId> triangleTuple
        = {i, triangleVec[0], triangleVec[1], triangleVec[2]};
      SimplexId nodeNum = vertexIndices_[triangleVec[0]];
      nodeTriangles[nodeNum].push_back(triangleTuple);
    }
  }

  for(auto iter = nodeTriangles.begin(); iter != nodeTriangles.end(); iter++) {
    ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
    if(exnode) {
      if(exnode->internalTriangleMap_.empty())
        buildInternalTriangleMap(exnode, false, true);
      for(std::vector<SimplexId> triangleVec : iter->second) {
        std::array<SimplexId, 3> triangle
          = {triangleVec[1], triangleVec[2], triangleVec[3]};
        (nodePtr->tetraTriangles_)[triangleVec[0]
                                   - cellIntervals_[nodePtr->nid - 1] - 1]
          .back()
          = exnode->internalTriangleMap_.at(triangle)
            + triangleIntervals_[iter->first - 1];
      }
    } else {
      ImplicitCluster tmpCluster(iter->first);
      buildInternalTriangleMap(&tmpCluster, false, true);
      for(std::vector<SimplexId> triangleVec : iter->second) {
        std::array<SimplexId, 3> triangle
          = {triangleVec[1], triangleVec[2], triangleVec[3]};
        (nodePtr->tetraTriangles_)[triangleVec[0]
                                   - cellIntervals_[nodePtr->nid - 1] - 1]
          .back()
          = tmpCluster.internalTriangleMap_.at(triangle)
            + triangleIntervals_[iter->first - 1];
      }
    }
  }

  return 0;
}

int CompactTriangulation::getClusterEdgeLinks(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localEdgeNum
    = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
    linksCount(localEdgeNum, 0);

  if(getDimensionality() == 2) {
    if(nodePtr->edgeStars_.empty()) {
      getClusterEdgeStars(nodePtr);
    }
    // set the offsets vector
    boost::unordered_map<std::array<SimplexId, 2>, SimplexId>::const_iterator
      iter;
    for(iter = nodePtr->internalEdgeMap_.begin();
        iter != nodePtr->internalEdgeMap_.end(); iter++) {
      for(SimplexId j = 0; j < nodePtr->edgeStars_.size(iter->second - 1);
          j++) {
        SimplexId cellId = nodePtr->edgeStars_.get(iter->second - 1, j);
        for(int k = 0; k < 3; k++) {
          SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
          if((vertexId != iter->first[0]) && (vertexId != iter->first[1])) {
            offsets[iter->second]++;
            break;
          }
        }
      }
    }

    // compute partial sum of number of links per edge
    for(SimplexId i = 1; i <= localEdgeNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for edge link data
    std::vector<SimplexId> edgeLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(iter = nodePtr->internalEdgeMap_.begin();
        iter != nodePtr->internalEdgeMap_.end(); iter++) {
      for(SimplexId j = 0; j < nodePtr->edgeStars_.size(iter->second - 1);
          j++) {
        SimplexId cellId = nodePtr->edgeStars_.get(iter->second - 1, j);
        for(int k = 0; k < 3; k++) {
          SimplexId vertexId = cellArray_->getCellVertex(cellId, k);
          if((vertexId != iter->first[0]) && (vertexId != iter->first[1])) {
            SimplexId localEdgeId = iter->second - 1;
            edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
              = vertexId;
            linksCount[localEdgeId]++;
            break;
          }
        }
      }
    }

    // fill FlatJaggedArray struct
    nodePtr->edgeLinks_.setData(std::move(edgeLinkData), std::move(offsets));
  } else if(getDimensionality() == 3) {
    if(nodePtr->tetraEdges_.empty()) {
      getClusterTetraEdges(nodePtr);
    }
    if(nodePtr->internalEdgeMap_.empty()) {
      buildInternalEdgeMap(nodePtr, false, true);
    }

    // set the offsets vector
    SimplexId localCellNum
      = cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1];
    for(SimplexId cid = 0; cid < localCellNum; cid++) {
      SimplexId cellId = cid + cellIntervals_[nodePtr->nid - 1] + 1;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
           (SimplexId)cellArray_->getCellVertex(cellId, 1),
           (SimplexId)cellArray_->getCellVertex(cellId, 2),
           (SimplexId)cellArray_->getCellVertex(cellId, 3)};
      std::array<SimplexId, 2> edgePair;
      edgePair[0] = vertexIds[0];
      for(SimplexId j = 1; j < 4; j++) {
        edgePair[1] = vertexIds[j];
        offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
      }
      if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        edgePair[0] = vertexIds[1];
        for(int j = 2; j < 4; j++) {
          edgePair[1] = vertexIds[j];
          offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
        }
        if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
          edgePair = {vertexIds[2], vertexIds[3]};
          offsets[nodePtr->internalEdgeMap_.at(edgePair)]++;
        }
      }
    }

    // loop through the external cell list
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      std::array<SimplexId, 2> edgeIds;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2),
           (SimplexId)cellArray_->getCellVertex(cid, 3)};

      // loop through each edge of the cell
      for(SimplexId j = 0; j < 3; j++) {
        for(SimplexId k = j + 1; k < 4; k++) {
          edgeIds[0] = vertexIds[j];
          edgeIds[1] = vertexIds[k];

          // the edge is in the current node
          if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
             && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
            offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
          }
        }
      }
    }

    // compute partial sum of number of neighbors per vertex
    for(SimplexId i = 1; i <= localEdgeNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for edge link data
    std::vector<SimplexId> edgeLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = 0; cid < localCellNum; cid++) {
      SimplexId cellId = cid + cellIntervals_[nodePtr->nid - 1] + 1;
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cellId, 0),
           (SimplexId)cellArray_->getCellVertex(cellId, 1),
           (SimplexId)cellArray_->getCellVertex(cellId, 2),
           (SimplexId)cellArray_->getCellVertex(cellId, 3)};
      std::array<SimplexId, 2> edgePair;
      edgePair[0] = vertexIds[0];
      for(SimplexId j = 1; j < 4; j++) {
        edgePair[1] = vertexIds[j];
        SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgePair) - 1;
        edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
          = nodePtr->tetraEdges_.at(cid)[6 - j];
        linksCount[localEdgeId]++;
      }
      if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        edgePair[0] = vertexIds[1];
        for(int j = 2; j < 4; j++) {
          edgePair[1] = vertexIds[j];
          SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgePair) - 1;
          edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
            = nodePtr->tetraEdges_.at(cid)[4 - j];
          linksCount[localEdgeId]++;
        }
        if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
          edgePair = {vertexIds[2], vertexIds[3]};
          SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgePair) - 1;
          edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
            = nodePtr->tetraEdges_.at(cid)[0];
          linksCount[localEdgeId]++;
        }
      }
    }

    // loop through the external cell list
    boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2),
           (SimplexId)cellArray_->getCellVertex(cid, 3)};

      std::array<SimplexId, 2> edgeIds;
      // loop through each edge of the cell
      for(SimplexId j = 0; j < 3; j++) {
        for(SimplexId k = j + 1; k < 4; k++) {
          edgeIds[0] = vertexIds[j];
          edgeIds[1] = vertexIds[k];

          // the edge is in the current node
          if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
             && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
            std::array<SimplexId, 2> otherEdge = {-1, -1};
            for(int i = 0; i < 4; i++) {
              if(vertexIds[i] != edgeIds[0] && vertexIds[i] != edgeIds[1]) {
                if(otherEdge[0] == -1) {
                  otherEdge[0] = vertexIds[i];
                } else if(otherEdge[1] == -1) {
                  otherEdge[1] = vertexIds[i];
                } else {
                  printErr("[CompactTriangulation] More than two other "
                           "vertices are "
                           "found in the edge!\n");
                }
              }
            }
            SimplexId nodeId = vertexIndices_[otherEdge[0]];
            if(nodeMaps.find(nodeId) == nodeMaps.end()) {
              nodeMaps[nodeId] = ImplicitCluster(nodeId);
              buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
            }
            SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
            edgeLinkData[offsets[localEdgeId] + linksCount[localEdgeId]]
              = nodeMaps[nodeId].internalEdgeMap_.at(otherEdge);
            linksCount[localEdgeId]++;
          }
        }
      }
    }

    // fill FlatJaggedArray struct
    nodePtr->edgeLinks_.setData(std::move(edgeLinkData), std::move(offsets));
  }

  return 0;
}

int CompactTriangulation::getClusterEdgeStars(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  SimplexId localEdgeNum
    = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localEdgeNum + 1, 0), starsCount(localEdgeNum);

  if(nodePtr->internalEdgeMap_.empty()) {
    buildInternalEdgeMap(nodePtr, false, true);
  }

  // set the offsets vector
  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      // the edge does not belong to the current node
      if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
      }
    }
  }
  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
           && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
          offsets[nodePtr->internalEdgeMap_.at(edgeIds)]++;
        }
      }
    }
  }

  // compute partial sum of number of stars per edge
  for(SimplexId i = 1; i <= localEdgeNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for edge star data
  std::vector<SimplexId> edgeStarData(offsets.back());

  // fill the flat vector using offsets and count vectors
  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      edgeIds[0] = cellArray_->getCellVertex(cid, j);
      // the edge does not belong to the current node
      if(edgeIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
        edgeStarData[offsets[localEdgeId] + starsCount[localEdgeId]] = cid;
        starsCount[localEdgeId]++;
      }
    }
  }
  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 2> edgeIds;
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        edgeIds[0] = cellArray_->getCellVertex(cid, j);
        edgeIds[1] = cellArray_->getCellVertex(cid, k);
        if(edgeIds[0] > vertexIntervals_[nodePtr->nid - 1]
           && edgeIds[0] <= vertexIntervals_[nodePtr->nid]) {
          SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edgeIds) - 1;
          edgeStarData[offsets[localEdgeId] + starsCount[localEdgeId]] = cid;
          starsCount[localEdgeId]++;
        }
      }
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->edgeStars_.setData(std::move(edgeStarData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterEdgeTriangles(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localEdgeNum
    = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localEdgeNum + 1, 0),
    trianglesCount(localEdgeNum, 0);

  if(nodePtr->internalEdgeMap_.empty()) {
    buildInternalEdgeMap(nodePtr, false, true);
  }
  if(nodePtr->internalTriangleMap_.empty()) {
    buildInternalTriangleMap(nodePtr, false, true);
  }
  if(nodePtr->externalTriangleMap_.empty()) {
    buildExternalTriangleMap(nodePtr);
  }

  // set the offsets vector
  boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::iterator iter;
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    std::array<SimplexId, 2> edge1 = {iter->first[0], iter->first[1]};
    std::array<SimplexId, 2> edge2 = {iter->first[0], iter->first[2]};
    offsets[nodePtr->internalEdgeMap_.at(edge1)]++;
    offsets[nodePtr->internalEdgeMap_.at(edge2)]++;
    if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
      edge1 = {iter->first[1], iter->first[2]};
      offsets[nodePtr->internalEdgeMap_.at(edge1)]++;
    }
  }
  for(iter = nodePtr->externalTriangleMap_.begin();
      iter != nodePtr->externalTriangleMap_.end(); iter++) {
    std::array<SimplexId, 2> edge = {iter->first.at(1), iter->first.at(2)};
    if(edge[0] > vertexIntervals_[nodePtr->nid - 1]
       && edge[0] <= vertexIntervals_[nodePtr->nid]) {
      offsets[nodePtr->internalEdgeMap_.at(edge)]++;
    }
  }

  // compute partial sum of number of triangles per edge
  for(SimplexId i = 1; i <= localEdgeNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for edge triangle data
  std::vector<SimplexId> edgeTriangleData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    std::array<SimplexId, 2> edge1 = {iter->first[0], iter->first[1]};
    std::array<SimplexId, 2> edge2 = {iter->first[0], iter->first[2]};
    SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edge1) - 1;
    edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
      = iter->second + triangleIntervals_[nodePtr->nid - 1];
    trianglesCount[localEdgeId]++;
    localEdgeId = nodePtr->internalEdgeMap_.at(edge2) - 1;
    edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
      = iter->second + triangleIntervals_[nodePtr->nid - 1];
    trianglesCount[localEdgeId]++;

    if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
      edge1 = {iter->first[1], iter->first[2]};
      localEdgeId = nodePtr->internalEdgeMap_.at(edge1) - 1;
      edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
        = iter->second + triangleIntervals_[nodePtr->nid - 1];
      trianglesCount[localEdgeId]++;
    }
  }

  // for external triangles
  for(iter = nodePtr->externalTriangleMap_.begin();
      iter != nodePtr->externalTriangleMap_.end(); iter++) {
    std::array<SimplexId, 2> edge = {iter->first.at(1), iter->first.at(2)};
    if(edge[0] > vertexIntervals_[nodePtr->nid - 1]
       && edge[0] <= vertexIntervals_[nodePtr->nid]) {
      SimplexId localEdgeId = nodePtr->internalEdgeMap_.at(edge) - 1;
      edgeTriangleData[offsets[localEdgeId] + trianglesCount[localEdgeId]]
        = iter->second;
      trianglesCount[localEdgeId]++;
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->edgeTriangles_.setData(
    std::move(edgeTriangleData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterTetraEdges(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = 4;
  nodePtr->tetraEdges_ = std::vector<std::array<SimplexId, 6>>(
    cellIntervals_[nodePtr->nid] - cellIntervals_[nodePtr->nid - 1]);
  boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
    edgeNodes;

  if(nodePtr->internalEdgeMap_.empty()) {
    buildInternalEdgeMap(nodePtr, false, true);
  }

  for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
      i <= cellIntervals_[nodePtr->nid]; i++) {
    int cnt = 0;
    // get the internal edge id from the map
    for(SimplexId k = 1; k < verticesPerCell; k++) {
      std::array<SimplexId, 2> edgePair
        = {(SimplexId)cellArray_->getCellVertex(i, 0),
           (SimplexId)cellArray_->getCellVertex(i, k)};
      (nodePtr->tetraEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][cnt++]
        = nodePtr->internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodePtr->nid - 1];
    }
    for(SimplexId j = 1; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        std::array<SimplexId, 2> edgePair
          = {(SimplexId)cellArray_->getCellVertex(i, j),
             (SimplexId)cellArray_->getCellVertex(i, k)};
        if(edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
          (nodePtr
             ->tetraEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][cnt++]
            = nodePtr->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
        }
        // group the external edges by node id
        else {
          std::vector<SimplexId> edgeTuple{i, cnt++, edgePair[0], edgePair[1]};
          SimplexId nodeNum = vertexIndices_[edgePair[0]];
          edgeNodes[nodeNum].push_back(edgeTuple);
        }
      }
    }
  }

  for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
    ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
    if(exnode) {
      if(exnode->internalEdgeMap_.empty())
        buildInternalEdgeMap(exnode, false, true);
      for(std::vector<SimplexId> edgeTuple : iter->second) {
        std::array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
        (nodePtr->tetraEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                               - 1][edgeTuple[1]]
          = exnode->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[iter->first - 1];
      }
    } else {
      ImplicitCluster tmpCluster(iter->first);
      buildInternalEdgeMap(&tmpCluster, false, true);
      for(std::vector<SimplexId> edgeTuple : iter->second) {
        std::array<SimplexId, 2> edgePair = {edgeTuple[2], edgeTuple[3]};
        (nodePtr->tetraEdges_)[edgeTuple[0] - cellIntervals_[nodePtr->nid - 1]
                               - 1][edgeTuple[1]]
          = tmpCluster.internalEdgeMap_.at(edgePair)
            + edgeIntervals_[iter->first - 1];
      }
    }
  }

  return 0;
}

int CompactTriangulation::getClusterTriangleEdges(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif
  nodePtr->triangleEdges_ = std::vector<std::array<SimplexId, 3>>(
    triangleIntervals_[nodePtr->nid] - triangleIntervals_[nodePtr->nid - 1]);
  boost::unordered_map<SimplexId, std::vector<std::vector<SimplexId>>>
    edgeNodes;

  if(nodePtr->internalEdgeMap_.empty()) {
    buildInternalEdgeMap(nodePtr, false, true);
  }

  if(getDimensionality() == 2) {
    for(SimplexId i = cellIntervals_[nodePtr->nid - 1] + 1;
        i <= cellIntervals_[nodePtr->nid]; i++) {
      // {v0, v1}
      std::array<SimplexId, 2> edgePair
        = {(SimplexId)cellArray_->getCellVertex(i, 0),
           (SimplexId)cellArray_->getCellVertex(i, 1)};
      (nodePtr->triangleEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][0]
        = nodePtr->internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodePtr->nid - 1];
      // {v0, v2};
      edgePair[1] = cellArray_->getCellVertex(i, 2);
      (nodePtr->triangleEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][1]
        = nodePtr->internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodePtr->nid - 1];
      // {v1, v2}
      edgePair[0] = cellArray_->getCellVertex(i, 1);
      if(edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
        (nodePtr->triangleEdges_)[i - cellIntervals_[nodePtr->nid - 1] - 1][2]
          = nodePtr->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
      }
      // group the external edges by node id
      else {
        std::vector<SimplexId> edgeTuple{i, edgePair[0], edgePair[1]};
        SimplexId nodeNum = vertexIndices_[edgePair[0]];
        edgeNodes[nodeNum].push_back(edgeTuple);
      }
    }

    for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
      ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
      if(exnode) {
        if(exnode->internalEdgeMap_.empty())
          buildInternalEdgeMap(exnode, false, true);
        for(std::vector<SimplexId> edgeTuple : iter->second) {
          std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
          (nodePtr->triangleEdges_)[edgeTuple[0]
                                    - cellIntervals_[nodePtr->nid - 1] - 1][2]
            = exnode->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[iter->first - 1];
        }
      } else {
        ImplicitCluster tmpCluster(iter->first);
        buildInternalEdgeMap(&tmpCluster, false, true);
        for(std::vector<SimplexId> edgeTuple : iter->second) {
          std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
          (nodePtr->triangleEdges_)[edgeTuple[0]
                                    - cellIntervals_[nodePtr->nid - 1] - 1][2]
            = tmpCluster.internalEdgeMap_.at(edgePair)
              + edgeIntervals_[iter->first - 1];
        }
      }
    }
  } else if(getDimensionality() == 3) {
    if(nodePtr->internalTriangleMap_.empty()) {
      buildInternalTriangleMap(nodePtr, false, true);
    }

    for(auto iter = nodePtr->internalTriangleMap_.begin();
        iter != nodePtr->internalTriangleMap_.end(); iter++) {
      // since the first vertex of the triangle is in the node ...
      std::array<SimplexId, 2> edgePair = {iter->first[0], iter->first[1]};
      (nodePtr->triangleEdges_)[iter->second - 1][0]
        = nodePtr->internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodePtr->nid - 1];
      edgePair[1] = iter->first[2];
      (nodePtr->triangleEdges_)[iter->second - 1][1]
        = nodePtr->internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodePtr->nid - 1];
      edgePair[0] = iter->first[1];
      if(edgePair[0] > vertexIntervals_[nodePtr->nid - 1]
         && edgePair[0] <= vertexIntervals_[nodePtr->nid]) {
        (nodePtr->triangleEdges_)[iter->second - 1][2]
          = nodePtr->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
      } else {
        std::vector<SimplexId> edgeTuple{
          iter->second - 1, edgePair[0], edgePair[1]};
        SimplexId nodeNum = vertexIndices_[edgePair[0]];
        edgeNodes[nodeNum].push_back(edgeTuple);
      }
    }

    for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++) {
      ImplicitCluster *exnode = searchCache(iter->first, nodePtr->nid);
      if(exnode) {
        if(exnode->internalEdgeMap_.empty())
          buildInternalEdgeMap(exnode, false, true);
        for(std::vector<SimplexId> edgeTuple : iter->second) {
          std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
          (nodePtr->triangleEdges_)[edgeTuple[0]][2]
            = exnode->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[iter->first - 1];
        }
      } else {
        ImplicitCluster tmpCluster(iter->first);
        buildInternalEdgeMap(&tmpCluster, false, true);
        for(std::vector<SimplexId> edgeTuple : iter->second) {
          std::array<SimplexId, 2> edgePair = {edgeTuple[1], edgeTuple[2]};
          (nodePtr->triangleEdges_)[edgeTuple[0]][2]
            = tmpCluster.internalEdgeMap_.at(edgePair)
              + edgeIntervals_[iter->first - 1];
        }
      }
    }
  }

  return 0;
}

int CompactTriangulation::getClusterTriangleLinks(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localTriangleNum
    = triangleIntervals_[nodePtr->nid] - triangleIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localTriangleNum + 1, 0),
    linksCount(localTriangleNum, 0);

  if(nodePtr->triangleStars_.empty()) {
    getClusterTriangleStars(nodePtr);
  }

  // set the offsets vector
  boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::const_iterator
    iter;
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    for(SimplexId i = 0; i < nodePtr->triangleStars_.size(iter->second - 1);
        i++) {
      SimplexId cellId = nodePtr->triangleStars_.get(iter->second - 1, i);
      for(int j = 0; j < 4; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
        if((vertexId != iter->first[0]) && (vertexId != iter->first[1])
           && (vertexId != iter->first[2])) {
          offsets[iter->second]++;
          break;
        }
      }
    }
  }

  // compute partial sum of number of links per triangle
  for(SimplexId i = 1; i <= localTriangleNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for triangle link data
  std::vector<SimplexId> triangleLinkData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    for(SimplexId i = 0; i < nodePtr->triangleStars_.size(iter->second - 1);
        i++) {
      SimplexId cellId = nodePtr->triangleStars_.get(iter->second - 1, i);
      for(int j = 0; j < 4; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cellId, j);
        if((vertexId != iter->first[0]) && (vertexId != iter->first[1])
           && (vertexId != iter->first[2])) {
          triangleLinkData[offsets[iter->second - 1]
                           + linksCount[iter->second - 1]]
            = vertexId;
          linksCount[iter->second - 1]++;
          break;
        }
      }
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->triangleLinks_.setData(
    std::move(triangleLinkData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterTriangleStars(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  SimplexId localTriangleNum
    = triangleIntervals_[nodePtr->nid] - triangleIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localTriangleNum + 1, 0),
    starsCount(localTriangleNum, 0);

  if(nodePtr->internalTriangleMap_.empty()) {
    buildInternalTriangleMap(nodePtr, false, true);
  }

  // set the offsets vector
  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 3> triangleIds;
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
        for(SimplexId l = k + 1; l < verticesPerCell; l++) {
          triangleIds[1] = cellArray_->getCellVertex(cid, k);
          triangleIds[2] = cellArray_->getCellVertex(cid, l);
          offsets[nodePtr->internalTriangleMap_.at(triangleIds)]++;
        }
      }
    }
  }
  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 3> triangleIds;
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
         && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
        for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            offsets[nodePtr->internalTriangleMap_.at(triangleIds)]++;
          }
        }
      }
    }
  }

  // compute partial sum of number of stars per triangle
  for(SimplexId i = 1; i <= localTriangleNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for triangle star data
  std::vector<SimplexId> triangleStarData(offsets.back());

  // fill the flat vector using offsets and count vectors
  // loop through the internal cell list
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    std::array<SimplexId, 3> triangleIds;
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      // the triangle does not belong to the current node
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid]) {
        break;
      }
      for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
        for(SimplexId l = k + 1; l < verticesPerCell; l++) {
          triangleIds[1] = cellArray_->getCellVertex(cid, k);
          triangleIds[2] = cellArray_->getCellVertex(cid, l);
          SimplexId localTriangleId
            = nodePtr->internalTriangleMap_.at(triangleIds) - 1;
          triangleStarData[offsets[localTriangleId]
                           + starsCount[localTriangleId]]
            = cid;
          starsCount[localTriangleId]++;
        }
      }
    }
  }
  // loop through the external cell list
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    std::array<SimplexId, 3> triangleIds;

    // loop through each triangle of the cell
    for(SimplexId j = 0; j < verticesPerCell - 2; j++) {
      triangleIds[0] = cellArray_->getCellVertex(cid, j);
      if(triangleIds[0] > vertexIntervals_[nodePtr->nid - 1]
         && triangleIds[0] <= vertexIntervals_[nodePtr->nid]) {
        for(SimplexId k = j + 1; k < verticesPerCell - 1; k++) {
          for(SimplexId l = k + 1; l < verticesPerCell; l++) {
            triangleIds[1] = cellArray_->getCellVertex(cid, k);
            triangleIds[2] = cellArray_->getCellVertex(cid, l);
            SimplexId localTriangleId
              = nodePtr->internalTriangleMap_.at(triangleIds) - 1;
            triangleStarData[offsets[localTriangleId]
                             + starsCount[localTriangleId]]
              = cid;
            starsCount[localTriangleId]++;
          }
        }
      }
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->triangleStars_.setData(
    std::move(triangleStarData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterVertexEdges(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    edgesCount(localVertexNum, 0);

  if(nodePtr->internalEdgeMap_.empty()) {
    buildInternalEdgeMap(nodePtr, false, true);
  }
  if(nodePtr->externalEdgeMap_.empty()) {
    buildExternalEdgeMap(nodePtr);
  }

  // set the offsets vector
  boost::unordered_map<std::array<SimplexId, 2>, SimplexId>::iterator iter;
  for(iter = nodePtr->internalEdgeMap_.begin();
      iter != nodePtr->internalEdgeMap_.end(); iter++) {
    offsets[iter->first[0] - vertexIntervals_[nodePtr->nid - 1]]++;
    if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
      offsets[iter->first[1] - vertexIntervals_[nodePtr->nid - 1]]++;
    }
  }
  for(iter = nodePtr->externalEdgeMap_.begin();
      iter != nodePtr->externalEdgeMap_.end(); iter++) {
    offsets[iter->first[1] - vertexIntervals_[nodePtr->nid - 1]]++;
  }
  // compute partial sum of number of edges per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex edge data
  std::vector<SimplexId> vertexEdgeData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(iter = nodePtr->internalEdgeMap_.begin();
      iter != nodePtr->internalEdgeMap_.end(); iter++) {
    SimplexId localVertexId
      = iter->first[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
    vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
      = iter->second + edgeIntervals_[nodePtr->nid - 1];
    edgesCount[localVertexId]++;
    if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
      localVertexId = iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
      vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
        = iter->second + edgeIntervals_[nodePtr->nid - 1];
      edgesCount[localVertexId]++;
    }
  }
  for(iter = nodePtr->externalEdgeMap_.begin();
      iter != nodePtr->externalEdgeMap_.end(); iter++) {
    SimplexId localVertexId
      = iter->first[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
    vertexEdgeData[offsets[localVertexId] + edgesCount[localVertexId]]
      = iter->second;
    edgesCount[localVertexId]++;
  }

  // fill FlatJaggedArray struct
  nodePtr->vertexEdges_.setData(std::move(vertexEdgeData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterVertexLinks(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    linksCount(localVertexNum, 0);
  boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
  // triangle mesh
  if(getDimensionality() == 2) {
    if(nodePtr->internalEdgeMap_.empty()) {
      buildInternalEdgeMap(nodePtr, false, true);
    }

    // set the offsets vector
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      offsets[cellArray_->getCellVertex(cid, 0)
              - vertexIntervals_[nodePtr->nid - 1]]++;
      for(SimplexId j = 1; j < 3; j++) {
        if(cellArray_->getCellVertex(cid, j)
           <= vertexIntervals_[nodePtr->nid]) {
          offsets[cellArray_->getCellVertex(cid, j)
                  - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
    }
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      for(SimplexId j = 1; j < 3; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cid, j);
        if(vertexId > vertexIntervals_[nodePtr->nid - 1]
           && vertexId <= vertexIntervals_[nodePtr->nid]) {
          offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
    }

    // compute partial sum of number of links per vertex
    for(SimplexId i = 1; i <= localVertexNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for vertex link data
    std::vector<SimplexId> vertexLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      std::array<SimplexId, 3> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2)};
      // the first vertex of the cell must be in the cluster
      std::array<SimplexId, 2> edgePair = {vertexIds[1], vertexIds[2]};
      SimplexId nodeId = vertexIndices_[vertexIds[1]];
      if(nodeMaps.find(nodeId) == nodeMaps.end()) {
        nodeMaps[nodeId] = ImplicitCluster(nodeId);
        buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
      }
      SimplexId localVertexId
        = vertexIds[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
      vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
        = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
          + edgeIntervals_[nodeId - 1];
      linksCount[localVertexId]++;

      if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        edgePair[0] = vertexIds[0];
        localVertexId = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodePtr->internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodePtr->nid - 1];
        linksCount[localVertexId]++;
        if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
          localVertexId = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
          edgePair[1] = vertexIds[1];
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = nodePtr->internalEdgeMap_.at(edgePair)
              + edgeIntervals_[nodePtr->nid - 1];
          linksCount[localVertexId]++;
        }
      }
    }
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      std::array<SimplexId, 3> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2)};
      std::array<SimplexId, 2> edgePair = {vertexIds[0], vertexIds[2]};
      SimplexId nodeId = vertexIndices_[edgePair[0]];
      if(nodeMaps.find(nodeId) == nodeMaps.end()) {
        nodeMaps[nodeId] = ImplicitCluster(nodeId);
        buildInternalEdgeMap(&nodeMaps[nodeId], false, true);
      }
      SimplexId localVertexId
        = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
      if(vertexIds[1] > vertexIntervals_[nodePtr->nid - 1]
         && vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[2] > vertexIntervals_[nodePtr->nid - 1]
         && vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
        edgePair[1] = vertexIds[1];
        localVertexId = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodeMaps[nodeId].internalEdgeMap_.at(edgePair)
            + edgeIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
    }

    // fill FlatJaggedArray struct
    nodePtr->vertexLinks_.setData(
      std::move(vertexLinkData), std::move(offsets));

    // tetrahedral mesh
  } else if(getDimensionality() == 3) {

    if(nodePtr->internalTriangleMap_.empty()) {
      buildInternalTriangleMap(nodePtr, false, true);
    }

    // set the offsets vector
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      offsets[cellArray_->getCellVertex(cid, 0)
              - vertexIntervals_[nodePtr->nid - 1]]++;
      for(SimplexId j = 1; j < 4; j++) {
        if(cellArray_->getCellVertex(cid, j)
           <= vertexIntervals_[nodePtr->nid]) {
          offsets[cellArray_->getCellVertex(cid, j)
                  - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
    }
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      for(SimplexId j = 1; j < 4; j++) {
        SimplexId vertexId = cellArray_->getCellVertex(cid, j);
        if(vertexId > vertexIntervals_[nodePtr->nid - 1]
           && vertexId <= vertexIntervals_[nodePtr->nid]) {
          offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
        }
      }
    }

    // compute partial sum of number of links per vertex
    for(SimplexId i = 1; i <= localVertexNum; i++) {
      offsets[i] += offsets[i - 1];
    }

    // allocate the flat vector for vertex link data
    std::vector<SimplexId> vertexLinkData(offsets.back());

    // fill the flat vector using offsets and count vectors
    for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
        cid <= cellIntervals_[nodePtr->nid]; cid++) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2),
           (SimplexId)cellArray_->getCellVertex(cid, 3)};

      // v1: (v2, v3, v4)
      std::array<SimplexId, 3> triangleVec
        = {vertexIds[1], vertexIds[2], vertexIds[3]};
      SimplexId nodeId = vertexIndices_[vertexIds[1]];
      if(nodeMaps.find(nodeId) == nodeMaps.end()) {
        nodeMaps[nodeId] = ImplicitCluster(nodeId);
        buildInternalTriangleMap(&nodeMaps[nodeId], false, true);
      }
      SimplexId localVertexId
        = vertexIds[0] - vertexIntervals_[nodePtr->nid - 1] - 1;
      vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
        = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
          + triangleIntervals_[nodeId - 1];
      linksCount[localVertexId]++;
      // v2: (v1, v3, v4)
      if(vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        triangleVec[0] = vertexIds[0];
        localVertexId = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodePtr->internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[nodePtr->nid - 1];
        linksCount[localVertexId]++;
        // v3: (v1, v2, v4)
        if(vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
          triangleVec[1] = vertexIds[1];
          localVertexId = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = nodePtr->internalTriangleMap_.at(triangleVec)
              + triangleIntervals_[nodePtr->nid - 1];
          linksCount[localVertexId]++;
        }
        // v4: (v1, v2, v3)
        if(vertexIds[3] <= vertexIntervals_[nodePtr->nid]) {
          triangleVec[2] = vertexIds[2];
          localVertexId = vertexIds[3] - vertexIntervals_[nodePtr->nid - 1] - 1;
          vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
            = nodePtr->internalTriangleMap_.at(triangleVec)
              + triangleIntervals_[nodePtr->nid - 1];
          linksCount[localVertexId]++;
        }
      }
    }

    // loop through the external cell list
    for(SimplexId cid : externalCells_[nodePtr->nid]) {
      std::array<SimplexId, 4> vertexIds
        = {(SimplexId)cellArray_->getCellVertex(cid, 0),
           (SimplexId)cellArray_->getCellVertex(cid, 1),
           (SimplexId)cellArray_->getCellVertex(cid, 2),
           (SimplexId)cellArray_->getCellVertex(cid, 3)};
      // start from v2
      std::array<SimplexId, 3> triangleVec
        = {vertexIds[0], vertexIds[2], vertexIds[3]};
      SimplexId nodeId = vertexIndices_[vertexIds[0]];
      if(nodeMaps.find(nodeId) == nodeMaps.end()) {
        nodeMaps[nodeId] = ImplicitCluster(nodeId);
        buildInternalTriangleMap(&nodeMaps[nodeId], false, true);
      }
      SimplexId localVertexId
        = vertexIds[1] - vertexIntervals_[nodePtr->nid - 1] - 1;
      if(vertexIds[1] > vertexIntervals_[nodePtr->nid - 1]
         && vertexIds[1] <= vertexIntervals_[nodePtr->nid]) {
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[2] > vertexIntervals_[nodePtr->nid - 1]
         && vertexIds[2] <= vertexIntervals_[nodePtr->nid]) {
        triangleVec[1] = vertexIds[1];
        localVertexId = vertexIds[2] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
      if(vertexIds[3] > vertexIntervals_[nodePtr->nid - 1]
         && vertexIds[3] <= vertexIntervals_[nodePtr->nid]) {
        triangleVec[1] = vertexIds[1];
        triangleVec[2] = vertexIds[2];
        localVertexId = vertexIds[3] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexLinkData[offsets[localVertexId] + linksCount[localVertexId]]
          = nodeMaps[nodeId].internalTriangleMap_.at(triangleVec)
            + triangleIntervals_[nodeId - 1];
        linksCount[localVertexId]++;
      }
    }

    // fill FlatJaggedArray struct
    nodePtr->vertexLinks_.setData(
      std::move(vertexLinkData), std::move(offsets));
  }

  return 0;
}

int CompactTriangulation::getClusterVertexNeighbors(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  SimplexId localVertexNum
    = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> vertexNeighborData, offsets(localVertexNum + 1, 0);

  SimplexId v1, v2;
  std::vector<boost::unordered_set<SimplexId>> vertexNeighborSet(
    localVertexNum);

  // loop through the internal cells
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      v1 = cellArray_->getCellVertex(cid, j);
      if(v1 <= vertexIntervals_[nodePtr->nid]) {
        for(SimplexId k = j + 1; k < verticesPerCell; k++) {
          v2 = cellArray_->getCellVertex(cid, k);
          vertexNeighborSet[v1 - vertexIntervals_[nodePtr->nid - 1] - 1].insert(
            v2);
          if(v2 <= vertexIntervals_[nodePtr->nid]) {
            vertexNeighborSet[v2 - vertexIntervals_[nodePtr->nid - 1] - 1]
              .insert(v1);
          }
        }
      }
    }
  }

  // loop through external cells
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    for(SimplexId j = 0; j < verticesPerCell - 1; j++) {
      for(SimplexId k = j + 1; k < verticesPerCell; k++) {
        v1 = cellArray_->getCellVertex(cid, j);
        v2 = cellArray_->getCellVertex(cid, k);
        if(v1 > vertexIntervals_[nodePtr->nid - 1]
           && v1 <= vertexIntervals_[nodePtr->nid])
          vertexNeighborSet[v1 - vertexIntervals_[nodePtr->nid - 1] - 1].insert(
            v2);
        if(v2 > vertexIntervals_[nodePtr->nid - 1]
           && v2 <= vertexIntervals_[nodePtr->nid])
          vertexNeighborSet[v2 - vertexIntervals_[nodePtr->nid - 1] - 1].insert(
            v1);
      }
    }
  }

  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] = offsets[i - 1] + vertexNeighborSet[i - 1].size();
    vertexNeighborData.insert(vertexNeighborData.end(),
                              vertexNeighborSet[i - 1].begin(),
                              vertexNeighborSet[i - 1].end());
  }

  nodePtr->vertexNeighbors_.setData(
    std::move(vertexNeighborData), std::move(offsets));
  return 0;
}

int CompactTriangulation::getClusterVertexStars(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId verticesPerCell = cellArray_->getCellVertexNumber(0);
  SimplexId localVertexNum
    = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    starsCount(localVertexNum, 0);

  // set the offsets vector
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    SimplexId vertexId = cellArray_->getCellVertex(cid, 0);
    offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      vertexId = cellArray_->getCellVertex(cid, j);
      if(vertexId > vertexIntervals_[nodePtr->nid - 1]
         && vertexId <= vertexIntervals_[nodePtr->nid]) {
        offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
      }
    }
  }
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      if(vertexId > vertexIntervals_[nodePtr->nid - 1]
         && vertexId <= vertexIntervals_[nodePtr->nid]) {
        offsets[vertexId - vertexIntervals_[nodePtr->nid - 1]]++;
      }
    }
  }

  // compute partial sum of number of stars per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex star data
  std::vector<SimplexId> vertexStarData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(SimplexId cid = cellIntervals_[nodePtr->nid - 1] + 1;
      cid <= cellIntervals_[nodePtr->nid]; cid++) {
    SimplexId localVertexId = cellArray_->getCellVertex(cid, 0)
                              - vertexIntervals_[nodePtr->nid - 1] - 1;
    vertexStarData[offsets[localVertexId] + starsCount[localVertexId]] = cid;
    starsCount[localVertexId]++;
    for(SimplexId j = 1; j < verticesPerCell; j++) {
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      // see if it is in the current node
      if(vertexId > vertexIntervals_[nodePtr->nid - 1]
         && vertexId <= vertexIntervals_[nodePtr->nid]) {
        localVertexId = vertexId - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
          = cid;
        starsCount[localVertexId]++;
      }
    }
  }
  for(SimplexId cid : externalCells_[nodePtr->nid]) {
    for(SimplexId j = 0; j < verticesPerCell; j++) {
      // see if it is in the current node
      SimplexId vertexId = cellArray_->getCellVertex(cid, j);
      if(vertexId > vertexIntervals_[nodePtr->nid - 1]
         && vertexId <= vertexIntervals_[nodePtr->nid]) {
        SimplexId localVertexId
          = vertexId - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexStarData[offsets[localVertexId] + starsCount[localVertexId]]
          = cid;
        starsCount[localVertexId]++;
      }
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->vertexStars_.setData(std::move(vertexStarData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getClusterVertexTriangles(
  ImplicitCluster *const nodePtr) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  SimplexId localVertexNum
    = vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1];
  std::vector<SimplexId> offsets(localVertexNum + 1, 0),
    trianglesCount(localVertexNum, 0);

  if(nodePtr->internalTriangleMap_.empty()) {
    buildInternalTriangleMap(nodePtr, false, true);
  }
  if(nodePtr->externalTriangleMap_.empty()) {
    buildExternalTriangleMap(nodePtr);
  }

  // set the offsets vector
  boost::unordered_map<std::array<SimplexId, 3>, SimplexId>::iterator iter;
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    for(SimplexId j = 0; j < 3; j++) {
      if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
         && iter->first[j] <= vertexIntervals_[nodePtr->nid])
        offsets[iter->first[j] - vertexIntervals_[nodePtr->nid - 1]]++;
    }
  }
  for(iter = nodePtr->externalTriangleMap_.begin();
      iter != nodePtr->externalTriangleMap_.end(); iter++) {
    for(SimplexId j = 0; j < 3; j++) {
      if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
         && iter->first[j] <= vertexIntervals_[nodePtr->nid])
        offsets[iter->first[j] - vertexIntervals_[nodePtr->nid - 1]]++;
    }
  }

  // compute partial sum of number of triangles per vertex
  for(SimplexId i = 1; i <= localVertexNum; i++) {
    offsets[i] += offsets[i - 1];
  }

  // allocate the flat vector for vertex triangle data
  std::vector<SimplexId> vertexTriangleData(offsets.back());

  // fill the flat vector using offsets and count vectors
  for(iter = nodePtr->internalTriangleMap_.begin();
      iter != nodePtr->internalTriangleMap_.end(); iter++) {
    for(SimplexId j = 0; j < 3; j++) {
      if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
         && iter->first[j] <= vertexIntervals_[nodePtr->nid]) {
        SimplexId localVertexId
          = iter->first[j] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexTriangleData[offsets[localVertexId]
                           + trianglesCount[localVertexId]]
          = iter->second + triangleIntervals_[nodePtr->nid - 1];
        trianglesCount[localVertexId]++;
      }
    }
  }
  for(iter = nodePtr->externalTriangleMap_.begin();
      iter != nodePtr->externalTriangleMap_.end(); iter++) {
    for(SimplexId j = 0; j < 3; j++) {
      if(iter->first[j] > vertexIntervals_[nodePtr->nid - 1]
         && iter->first[j] <= vertexIntervals_[nodePtr->nid]) {
        SimplexId localVertexId
          = iter->first[j] - vertexIntervals_[nodePtr->nid - 1] - 1;
        vertexTriangleData[offsets[localVertexId]
                           + trianglesCount[localVertexId]]
          = iter->second;
        trianglesCount[localVertexId]++;
      }
    }
  }

  // fill FlatJaggedArray struct
  nodePtr->vertexTriangles_.setData(
    std::move(vertexTriangleData), std::move(offsets));

  return 0;
}

int CompactTriangulation::getBoundaryCells(ImplicitCluster *const nodePtr,
                                           const SimplexId dim) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
    return -1;
#endif

  if(getDimensionality() == 2) {
    SimplexId localEdgeNum
      = edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1];
    if(nodePtr->boundaryEdges_.empty()) {
      nodePtr->boundaryEdges_ = std::vector<bool>(localEdgeNum, false);
      if(nodePtr->edgeStars_.empty()) {
        getClusterEdgeStars(nodePtr);
      }
      for(SimplexId i = 0; i < localEdgeNum; i++) {
        if(nodePtr->edgeStars_.size(i) == 1) {
          nodePtr->boundaryEdges_.at(i) = true;
        }
      }
    }
    // if boundary vertices are requested
    if(dim == 0 && nodePtr->boundaryVertices_.empty()) {
      nodePtr->boundaryVertices_ = std::vector<bool>(
        vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
        false);
      if(nodePtr->externalEdgeMap_.empty()) {
        buildExternalEdgeMap(nodePtr);
      }
      // internal edges
      for(auto iter = nodePtr->internalEdgeMap_.begin();
          iter != nodePtr->internalEdgeMap_.end(); iter++) {
        if((nodePtr->boundaryEdges_)[iter->second - 1]) {
          (nodePtr->boundaryVertices_)[iter->first[0]
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
            = true;
          if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
            (nodePtr
               ->boundaryVertices_)[iter->first[1]
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
              = true;
          }
        }
      }
      // external edges
      boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
      for(auto iter = nodePtr->externalEdgeMap_.begin();
          iter != nodePtr->externalEdgeMap_.end(); iter++) {
        SimplexId nodeId = vertexIndices_[iter->first[0]];
        if(nodeMaps.find(nodeId) == nodeMaps.end()) {
          nodeMaps[nodeId] = ImplicitCluster(nodeId);
          getBoundaryCells(&nodeMaps[nodeId]);
        }
        if((nodeMaps[nodeId]
              .boundaryEdges_)[iter->second - edgeIntervals_[nodeId - 1] - 1]) {
          (nodePtr->boundaryVertices_)[iter->first[1]
                                       - vertexIntervals_[nodePtr->nid - 1] - 1]
            = true;
        }
      }
    }
  } else if(getDimensionality() == 3) {
    // get the boundary triangles first
    SimplexId localTriangleNum
      = triangleIntervals_[nodePtr->nid] - triangleIntervals_[nodePtr->nid - 1];
    if(nodePtr->boundaryTriangles_.empty()) {
      nodePtr->boundaryTriangles_ = std::vector<bool>(localTriangleNum, false);
      if(nodePtr->triangleStars_.empty()) {
        getClusterTriangleStars(nodePtr);
      }
      for(SimplexId i = 0; i < localTriangleNum; i++) {
        if(nodePtr->triangleStars_.size(i) == 1) {
          (nodePtr->boundaryTriangles_)[i] = true;
        }
      }
    }
    // if the boundary edges are requested
    if(dim == 1 && nodePtr->boundaryEdges_.empty()) {
      nodePtr->boundaryEdges_ = std::vector<bool>(
        edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid - 1], false);
      if(nodePtr->externalTriangleMap_.empty()) {
        buildExternalTriangleMap(nodePtr);
      }
      if(nodePtr->internalEdgeMap_.empty()) {
        buildInternalEdgeMap(nodePtr, false, true);
      }
      // internal triangles
      for(auto iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        if((nodePtr->boundaryTriangles_)[iter->second - 1]) {
          std::array<SimplexId, 2> edgePair = {iter->first[0], iter->first[1]};
          (nodePtr->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
            = true;
          edgePair[1] = iter->first[2];
          (nodePtr->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
            = true;
          if(iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
            edgePair[0] = iter->first[1];
            (nodePtr
               ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
              = true;
          }
        }
      }
      // external triangles
      boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
      for(auto iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        SimplexId nodeId = vertexIndices_[iter->first[0]];
        if(nodeMaps.find(nodeId) == nodeMaps.end()) {
          nodeMaps[nodeId] = ImplicitCluster(nodeId);
          getBoundaryCells(&nodeMaps[nodeId]);
        }
        if((nodeMaps[nodeId]
              .boundaryTriangles_)[iter->second - triangleIntervals_[nodeId - 1]
                                   - 1]) {
          if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
            std::array<SimplexId, 2> edgePair
              = {iter->first[1], iter->first[2]};
            (nodePtr
               ->boundaryEdges_)[nodePtr->internalEdgeMap_.at(edgePair) - 1]
              = true;
          }
        }
      }
    }

    // if the boundary vertices are requested
    else if(dim == 0 && nodePtr->boundaryVertices_.empty()) {
      nodePtr->boundaryVertices_ = std::vector<bool>(
        vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid - 1],
        false);
      if(nodePtr->externalTriangleMap_.empty()) {
        buildExternalTriangleMap(nodePtr);
      }
      // internal triangles
      for(auto iter = nodePtr->internalTriangleMap_.begin();
          iter != nodePtr->internalTriangleMap_.end(); iter++) {
        if((nodePtr->boundaryTriangles_)[iter->second - 1]) {
          for(int j = 0; j < 3; j++) {
            SimplexId vid = iter->first[j];
            if(vid <= vertexIntervals_[nodePtr->nid]) {
              (nodePtr
                 ->boundaryVertices_)[vid - vertexIntervals_[nodePtr->nid - 1]
                                      - 1]
                = true;
            }
          }
        }
      }
      // external triangles
      boost::unordered_map<SimplexId, ImplicitCluster> nodeMaps;
      for(auto iter = nodePtr->externalTriangleMap_.begin();
          iter != nodePtr->externalTriangleMap_.end(); iter++) {
        SimplexId nodeId = vertexIndices_[iter->first[0]];
        if(nodeMaps.find(nodeId) == nodeMaps.end()) {
          nodeMaps[nodeId] = ImplicitCluster(nodeId);
          ;
          getBoundaryCells(&nodeMaps[nodeId]);
        }
        if((nodeMaps[nodeId]
              .boundaryTriangles_)[iter->second - triangleIntervals_[nodeId - 1]
                                   - 1]) {
          if(iter->first[1] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[1] <= vertexIntervals_[nodePtr->nid]) {
            (nodePtr
               ->boundaryVertices_)[iter->first[1]
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
              = true;
          }
          if(iter->first[2] > vertexIntervals_[nodePtr->nid - 1]
             && iter->first[2] <= vertexIntervals_[nodePtr->nid]) {
            (nodePtr
               ->boundaryVertices_)[iter->first[2]
                                    - vertexIntervals_[nodePtr->nid - 1] - 1]
              = true;
          }
        }
      }
    }
  } else {
    return -1;
  }

  return 0;
}
