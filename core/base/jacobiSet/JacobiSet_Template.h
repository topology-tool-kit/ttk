#include "JacobiSet.h"

template <class dataTypeU, class dataTypeV, typename triangulationType>
int ttk::JacobiSet::execute(std::vector<std::pair<SimplexId, char>> &jacobiSet,
                            const dataTypeU *const uField,
                            const dataTypeV *const vField,
                            const triangulationType &triangulation,
                            std::vector<char> *isPareto) {

  Timer t;

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((triangulation.isEmpty())) {
    if(vertexNumber_) {
      return executeLegacy(jacobiSet, uField, vField);
    }
    return -1;
  }
  if(!uField)
    return -2;
  if(!vField)
    return -3;
#endif

  jacobiSet.clear();

  SimplexId edgeNumber = triangulation.getNumberOfEdges();

  std::vector<std::vector<std::pair<SimplexId, char>>> threadedCriticalTypes(
    threadNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < edgeNumber; i++) {

    char type = getCriticalType(i, uField, vField, triangulation);

    if(type != -2) {
      // -2: regular vertex
      ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif
      threadedCriticalTypes[threadId].push_back(
        std::pair<SimplexId, char>(i, type));
    }
  }

  // now merge the threaded lists
  for(SimplexId i = 0; i < threadNumber_; i++) {
    for(SimplexId j = 0; j < (SimplexId)threadedCriticalTypes[i].size(); j++) {
      jacobiSet.push_back(threadedCriticalTypes[i][j]);
    }
  }

  SimplexId minimumNumber = 0, saddleNumber = 0, maximumNumber = 0,
            monkeySaddleNumber = 0;

  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {
    switch(jacobiSet[i].second) {
      case 0:
        minimumNumber++;
        break;
      case 1:
        saddleNumber++;
        break;
      case 2:
        maximumNumber++;
        break;
      case -1:
        monkeySaddleNumber++;
        break;
    }
  }

  if(isPareto) {
    isPareto->resize(jacobiSet.size(), 0);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int)jacobiSet.size(); i++) {
      int edgeId = jacobiSet[i].first;
      SimplexId vertexId0 = -1, vertexId1 = -1;
      triangulation.getEdgeVertex(edgeId, 0, vertexId0);
      triangulation.getEdgeVertex(edgeId, 1, vertexId1);
      double denominator = vField[vertexId1] - vField[vertexId0];
      if(fabs(denominator) < Geometry::powIntTen(-DBL_DIG))
        denominator = 1;
      if((uField[vertexId1] - uField[vertexId0]) / denominator < 0)
        (*isPareto)[i] = 1;
    }
  }

  this->printMsg(
    std::vector<std::vector<std::string>>{
      {"#Minimum edges", std::to_string(minimumNumber)},
      {"#Saddle edges", std::to_string(saddleNumber)},
      {"#Maximum edges", std::to_string(maximumNumber)},
      {"#Multi-saddle edges", std::to_string(monkeySaddleNumber)}},
    debug::Priority::DETAIL);

  this->printMsg(
    "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg("Jacobi edge rate: "
                 + std::to_string(100.0 * jacobiSet.size() / (double)edgeNumber)
                 + "%");

  return 0;
}

template <class dataTypeU, class dataTypeV>
int ttk::JacobiSet::executeLegacy(
  std::vector<std::pair<SimplexId, char>> &jacobiSet,
  const dataTypeU *const uField,
  const dataTypeV *const vField) {

  Timer t;

  this->printWrn("Using legacy implementation...");

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!vertexNumber_)
    return -1;
  if(!uField)
    return -2;
  if(!vField)
    return -3;
  if(!edgeList_)
    return -4;
  if(!edgeFanLinkEdgeLists_)
    return -5;
  if(!edgeFans_)
    return -6;
  if(!sosOffsetsU_)
    return -7;
#endif

  SimplexId count = 0;

  jacobiSet.clear();

  // distance fields (not really memory efficient)
  // for each thread
  //      for each vertex: distance field map
  std::vector<std::vector<double>> threadedDistanceField(threadNumber_);
  for(SimplexId i = 0; i < (SimplexId)threadedDistanceField.size(); i++) {
    threadedDistanceField[i].resize(vertexNumber_);
  }

  std::vector<ScalarFieldCriticalPoints> threadedCriticalPoints(threadNumber_);
  for(ThreadId i = 0; i < threadNumber_; i++) {
    threadedCriticalPoints[i].setDomainDimension(2);
    threadedCriticalPoints[i].setVertexNumber(vertexNumber_);
  }

  std::vector<std::vector<std::pair<SimplexId, char>>> threadedCriticalTypes(
    threadNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)edgeList_->size(); i++) {

    // avoid any processing if the abort signal is sent
    if((!wrapper_) || ((wrapper_) && (!wrapper_->needsToAbort()))) {

      ThreadId threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif

      // processing here!
      SimplexId pivotVertexId = (*edgeList_)[i].first;
      SimplexId otherExtremityId = (*edgeList_)[i].second;

      // A) compute the distance field
      double projectedPivotVertex[2];
      projectedPivotVertex[0] = uField[pivotVertexId];
      projectedPivotVertex[1] = vField[pivotVertexId];

      double projectedOtherVertex[2];
      projectedOtherVertex[0] = uField[otherExtremityId];
      projectedOtherVertex[1] = vField[otherExtremityId];

      double rangeEdge[2];
      rangeEdge[0] = projectedOtherVertex[0] - projectedPivotVertex[0];
      rangeEdge[1] = projectedOtherVertex[1] - projectedPivotVertex[1];

      double rangeNormal[2];
      rangeNormal[0] = -rangeEdge[1];
      rangeNormal[1] = rangeEdge[0];

      for(SimplexId j = 0; j < (SimplexId)(*edgeFans_)[i].size() / 4; j++) {
        for(int k = 0; k < 3; k++) {

          SimplexId vertexId = (*edgeFans_)[i][j * 4 + 1 + k];

          // we can compute the distance field (in the rage)
          double projectedVertex[2];
          projectedVertex[0] = uField[vertexId];
          projectedVertex[1] = vField[vertexId];

          double vertexRangeEdge[2];
          vertexRangeEdge[0] = projectedVertex[0] - projectedPivotVertex[0];
          vertexRangeEdge[1] = projectedVertex[1] - projectedPivotVertex[1];

          // signed distance: linear function of the dot product
          threadedDistanceField[threadId][vertexId]
            = vertexRangeEdge[0] * rangeNormal[0]
              + vertexRangeEdge[1] * rangeNormal[1];
        }
      }

      // B) compute critical points
      // watch out between local and global Ids
      // what I could do is to translate the ids from global to local
      // also, lots of things in there can be done out of the loop

      // in the loop
      char type = threadedCriticalPoints[threadId].getCriticalType(
        pivotVertexId, sosOffsetsU_, (*edgeFanLinkEdgeLists_)[i]);

      if(type != -2) {
        // -2: regular vertex
        threadedCriticalTypes[threadId].push_back(
          std::pair<SimplexId, char>(i, type));
      }

      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > static_cast<int>(debug::Priority::DETAIL)) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if((wrapper_) && (!(count % ((vertexNumber_) / 10)))) {
            wrapper_->updateProgress((count + 1.0) / vertexNumber_);
          }

          count++;
        }
      }
    }
  }

  // now merge the threaded lists
  for(ThreadId i = 0; i < threadNumber_; i++) {
    for(SimplexId j = 0; j < (SimplexId)threadedCriticalTypes[i].size(); j++) {
      jacobiSet.push_back(threadedCriticalTypes[i][j]);
    }
  }

  SimplexId minimumNumber = 0, saddleNumber = 0, maximumNumber = 0,
            monkeySaddleNumber = 0;

  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {
    switch(jacobiSet[i].second) {
      case 0:
        minimumNumber++;
        break;
      case 1:
        saddleNumber++;
        break;
      case 2:
        maximumNumber++;
        break;
      case -1:
        monkeySaddleNumber++;
        break;
    }
  }

  this->printMsg(
    std::vector<std::vector<std::string>>{
      {"#Minimum edges", std::to_string(minimumNumber)},
      {"#Saddle edges", std::to_string(saddleNumber)},
      {"#Maximum edges", std::to_string(maximumNumber)},
      {"#Multi-saddle edges", std::to_string(monkeySaddleNumber)}},
    debug::Priority::DETAIL);

  this->printMsg(
    "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg(
    "Jacobi edge rate: "
    + std::to_string(100.0 * (jacobiSet.size() / ((double)edgeList_->size())))
    + "%");

  return 0;
}

template <class dataTypeU, class dataTypeV, typename triangulationType>
char ttk::JacobiSet::getCriticalType(const SimplexId &edgeId,
                                     const dataTypeU *const uField,
                                     const dataTypeV *const vField,
                                     const triangulationType &triangulation) {

  SimplexId vertexId0 = -1, vertexId1 = -1;
  triangulation.getEdgeVertex(edgeId, 0, vertexId0);
  triangulation.getEdgeVertex(edgeId, 1, vertexId1);

  double projectedPivotVertex[2];
  projectedPivotVertex[0] = uField[vertexId0];
  projectedPivotVertex[1] = vField[vertexId0];

  double projectedOtherVertex[2];
  projectedOtherVertex[0] = uField[vertexId1];
  projectedOtherVertex[1] = vField[vertexId1];

  double rangeEdge[2];
  rangeEdge[0] = projectedOtherVertex[0] - projectedPivotVertex[0];
  rangeEdge[1] = projectedOtherVertex[1] - projectedPivotVertex[1];

  double rangeNormal[2];
  rangeNormal[0] = -rangeEdge[1];
  rangeNormal[1] = rangeEdge[0];

  SimplexId starNumber = triangulation.getEdgeStarNumber(edgeId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  SimplexId neighborNumber = 0;

  for(SimplexId i = 0; i < starNumber; i++) {

    SimplexId tetId = -1;
    triangulation.getEdgeStar(edgeId, i, tetId);

    SimplexId vertexNumber = triangulation.getCellVertexNumber(tetId);
    for(SimplexId j = 0; j < vertexNumber; j++) {
      SimplexId vertexId = -1;
      triangulation.getCellVertex(tetId, j, vertexId);

      if((vertexId != -1) && (vertexId != vertexId0)
         && (vertexId != vertexId1)) {
        // new neighbor
        bool isIn = false;
        for(SimplexId k = 0; k < (SimplexId)lowerNeighbors.size(); k++) {
          if(vertexId == lowerNeighbors[k]) {
            isIn = true;
            break;
          }
        }

        if(!isIn) {
          for(SimplexId k = 0; k < (SimplexId)upperNeighbors.size(); k++) {
            if(vertexId == upperNeighbors[k]) {
              isIn = true;
              break;
            }
          }
        }

        if(!isIn) {
          // compute the actual distance field
          // A) compute the distance field
          double projectedVertex[2];
          projectedVertex[0] = uField[vertexId];
          projectedVertex[1] = vField[vertexId];

          double vertexRangeEdge[2];
          vertexRangeEdge[0] = projectedVertex[0] - projectedPivotVertex[0];
          vertexRangeEdge[1] = projectedVertex[1] - projectedPivotVertex[1];

          // signed distance: linear function of the dot product
          double distance = vertexRangeEdge[0] * rangeNormal[0]
                            + vertexRangeEdge[1] * rangeNormal[1];

          neighborNumber++;

          if(distance < 0) {
            lowerNeighbors.push_back(vertexId);
          } else if(distance > 0) {
            upperNeighbors.push_back(vertexId);
          } else {
            // degenerate
            // compute the distance field out of the offset positions
            double offsetProjectedPivotVertex[2];
            offsetProjectedPivotVertex[0] = sosOffsetsU_[vertexId0];
            offsetProjectedPivotVertex[1]
              = sosOffsetsV_[vertexId0] * sosOffsetsV_[vertexId0];

            double offsetProjectedOtherVertex[2];
            offsetProjectedOtherVertex[0] = sosOffsetsU_[vertexId1];
            offsetProjectedOtherVertex[1]
              = sosOffsetsV_[vertexId1] * sosOffsetsV_[vertexId1];

            double offsetRangeEdge[2];
            offsetRangeEdge[0]
              = offsetProjectedOtherVertex[0] - offsetProjectedPivotVertex[0];
            offsetRangeEdge[1]
              = offsetProjectedOtherVertex[1] - offsetProjectedPivotVertex[1];

            double offsetRangeNormal[2];
            offsetRangeNormal[0] = -offsetRangeEdge[1];
            offsetRangeNormal[1] = offsetRangeEdge[0];

            projectedVertex[0] = sosOffsetsU_[vertexId];
            projectedVertex[1]
              = sosOffsetsV_[vertexId] * sosOffsetsV_[vertexId];

            vertexRangeEdge[0]
              = projectedVertex[0] - offsetProjectedPivotVertex[0];
            vertexRangeEdge[1]
              = projectedVertex[1] - offsetProjectedPivotVertex[1];

            distance = vertexRangeEdge[0] * offsetRangeNormal[0]
                       + vertexRangeEdge[1] * offsetRangeNormal[1];

            if(distance < 0) {
              lowerNeighbors.push_back(vertexId);
            } else if(distance > 0) {
              upperNeighbors.push_back(vertexId);
            } else {
              this->printWrn(
                "Inconsistent (non-bijective?) offsets for vertex #"
                + std::to_string(vertexId));
            }
          }
        }
      }
    }
  }

  // at this point, we know if each vertex of the edge link is higher or not.
  if((SimplexId)(lowerNeighbors.size() + upperNeighbors.size())
     != neighborNumber) {
    // Inconsistent offsets (cf above error message)
    return -2;
  }

  if(lowerNeighbors.empty()) {
    if(rangeNormal[0] + rangeNormal[1] > 0)
      // minimum
      return 0;
    // else maximum
    return triangulation.getDimensionality() - 1;
  }
  if(upperNeighbors.empty()) {
    if(rangeNormal[0] + rangeNormal[1] > 0)
      // maximum
      return triangulation.getDimensionality() - 1;
    // else minimum
    return 0;
  }

  // let's check the connectivity now
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  for(SimplexId i = 0; i < starNumber; i++) {

    SimplexId tetId = -1;
    triangulation.getEdgeStar(edgeId, i, tetId);

    SimplexId vertexNumber = triangulation.getCellVertexNumber(tetId);
    for(SimplexId j = 0; j < vertexNumber; j++) {
      SimplexId edgeVertexId0 = -1;
      triangulation.getCellVertex(tetId, j, edgeVertexId0);
      if((edgeVertexId0 != vertexId0) && (edgeVertexId0 != vertexId1)) {
        for(SimplexId k = j + 1; k < vertexNumber; k++) {
          SimplexId edgeVertexId1 = -1;
          triangulation.getCellVertex(tetId, k, edgeVertexId1);
          if((edgeVertexId1 != vertexId0) && (edgeVertexId1 != vertexId1)) {
            // processing the edge (edgeVertexId0, edgeVertexId1)

            // we need to find out if they're lower or not
            bool lower0 = false;
            for(SimplexId l = 0; l < (SimplexId)lowerNeighbors.size(); l++) {
              if(lowerNeighbors[l] == edgeVertexId0) {
                lower0 = true;
                break;
              }
            }
            bool lower1 = false;
            for(SimplexId l = 0; l < (SimplexId)lowerNeighbors.size(); l++) {
              if(lowerNeighbors[l] == edgeVertexId1) {
                lower1 = true;
                break;
              }
            }

            std::vector<SimplexId> *neighbors = &lowerNeighbors;
            std::vector<UnionFind *> *seeds = &lowerList;

            if(!lower0) {
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }

            if(lower0 == lower1) {
              // connect their union-find sets!
              SimplexId lowerId0 = -1, lowerId1 = -1;
              for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
                if((*neighbors)[l] == edgeVertexId0) {
                  lowerId0 = l;
                }
                if((*neighbors)[l] == edgeVertexId1) {
                  lowerId1 = l;
                }
              }

              if((lowerId0 != -1) && (lowerId1 != -1)) {
                (*seeds)[lowerId0] = UnionFind::makeUnion(
                  (*seeds)[lowerId0], (*seeds)[lowerId1]);
                (*seeds)[lowerId1] = (*seeds)[lowerId0];
              }
            }

            break;
          }
        }
      }
    }
  }

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if((upperList.size() == 1) && (lowerList.size() == 1))
    return -2;

  return 1;
}

template <class dataTypeU, class dataTypeV>
int ttk::JacobiSet::perturbate(const dataTypeU *const uField,
                               const dataTypeV *const vField,
                               const dataTypeU uEpsilon,
                               const dataTypeV vEpsilon) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(!uField)
    return -1;
  if(!vField)
    return -2;
  if(!vertexNumber_)
    return -3;
#endif

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    // simulation of simplicity in 2 dimensions, need to use degree 2 polynoms

    uField[i] += i * uEpsilon;
    vField[i] += (i * vEpsilon) * (i * vEpsilon);
  }

  return 0;
}
