#include <ReebSpace.h>

using namespace std;
using namespace ttk;

ReebSpace::ReebSpace() {

  vertexNumber_ = 0;
  tetNumber_ = 0;
  uField_ = NULL;
  vField_ = NULL;
  sosOffsetsU_ = NULL;
  sosOffsetsV_ = NULL;

  totalArea_ = -1;
  totalVolume_ = -1;
  totalHyperVolume_ = -1;

  hasConnectedSheets_ = false;
  expand3sheets_ = true;

  triangulation_ = NULL;
  withRangeDrivenOctree_ = true;
}

ReebSpace::~ReebSpace() {
}

int ReebSpace::compute1sheetsOnly(
  const vector<pair<SimplexId, char>> &jacobiSet,
  vector<pair<SimplexId, SimplexId>> &jacobiClassification) {

  // bfs on the saddle jacobi edges only to identify saddle 1-sheets as well as
  // saddle 0-sheets

  Timer t;

  // alloc
  jacobiClassification.reserve(jacobiSet.size());

  // markup the saddle jacobi and visisted edges
  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {
    originalData_.edgeTypes_[jacobiSet[i].first] = jacobiSet[i].second;
  }

  vector<bool> visitedEdges(triangulation_->getNumberOfEdges(), false);

  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {

    if(visitedEdges[jacobiSet[i].first] == false) {

      SimplexId sheet1Id = originalData_.sheet1List_.size();
      originalData_.sheet1List_.resize(originalData_.sheet1List_.size() + 1);
      originalData_.sheet1List_.back().hasSaddleEdges_ = false;
      originalData_.sheet1List_.back().pruned_ = false;

      queue<SimplexId> edgeQueue;
      edgeQueue.push(jacobiSet[i].first);

      do {

        SimplexId edgeId = edgeQueue.front();
        edgeQueue.pop();

        if(!visitedEdges[edgeId]) {

          jacobiClassification.push_back(
            pair<SimplexId, SimplexId>(edgeId, sheet1Id));

          originalData_.sheet1List_.back().edgeList_.push_back(edgeId);
          originalData_.edge2sheet1_[edgeId] = sheet1Id;

          if(originalData_.edgeTypes_[edgeId] == 1) {
            originalData_.sheet1List_.back().hasSaddleEdges_ = true;
          }

          SimplexId vertexId0 = -1;
          triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
          SimplexId vertexId1 = -1;
          triangulation_->getEdgeVertex(edgeId, 1, vertexId1);

          SimplexId neighborNumber = 0;

          SimplexId vertexEdgeNumber
            = triangulation_->getVertexEdgeNumber(vertexId0);

          for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

            SimplexId vertexEdgeId = -1;
            triangulation_->getVertexEdge(vertexId0, j, vertexEdgeId);

            if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
              if(vertexEdgeId != edgeId) {

                neighborNumber++;

                if(!visitedEdges[vertexEdgeId]) {
                  edgeQueue.push(vertexEdgeId);
                }
              }
            }
          }
          if((neighborNumber > 2)
             && (originalData_.vertex2sheet0_[vertexId0] == -1)) {
            // vertexId0 is a 0-sheet
            // mark up the segmentation
            originalData_.vertex2sheet0_[vertexId0]
              = originalData_.sheet0List_.size();

            originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                             + 1);
            originalData_.sheet0List_.back().vertexId_ = vertexId0;
            originalData_.sheet0List_.back().type_ = 1;
            originalData_.sheet0List_.back().pruned_ = false;
          }

          neighborNumber = 0;
          vertexEdgeNumber = triangulation_->getVertexEdgeNumber(vertexId1);
          for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

            SimplexId vertexEdgeId = -1;
            triangulation_->getVertexEdge(vertexId1, j, vertexEdgeId);

            if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
              if(vertexEdgeId != edgeId) {

                neighborNumber++;

                if(!visitedEdges[vertexEdgeId]) {
                  edgeQueue.push(vertexEdgeId);
                }
              }
            }
          }

          if((neighborNumber > 2)
             && (originalData_.vertex2sheet0_[vertexId1] == -1)) {
            // vertexId1 is a 0-sheet
            // mark up the segmentation

            originalData_.vertex2sheet0_[vertexId1]
              = originalData_.sheet0List_.size();

            originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                             + 1);
            originalData_.sheet0List_.back().vertexId_ = vertexId1;
            originalData_.sheet0List_.back().type_ = 1;
            originalData_.sheet0List_.back().pruned_ = false;
          }
        }

        visitedEdges[edgeId] = true;

      } while(edgeQueue.size());
    }
  }

  {
    stringstream msg;
    msg << "[ReebSpace] " << originalData_.sheet1List_.size()
        << " 1-sheets and " << originalData_.sheet0List_.size()
        << " 0-sheets extracted in " << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ReebSpace::compute1sheets(
  const vector<pair<SimplexId, char>> &jacobiSet,
  vector<pair<SimplexId, SimplexId>> &jacobiClassification) {

  // bfs on the saddle jacobi edges only to identify saddle 1-sheets as well as
  // saddle 0-sheets

  Timer t;

  // alloc
  jacobiClassification.reserve(jacobiSet.size());

  // markup the saddle jacobi and visisted edges
  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {
    originalData_.edgeTypes_[jacobiSet[i].first] = jacobiSet[i].second;
  }

  vector<bool> visitedEdges(triangulation_->getNumberOfEdges(), false);
  for(SimplexId i = 0; i < (SimplexId)jacobiSet.size(); i++) {

    if(/*(saddleEdge[jacobiSet[i].first] == 1)&&*/
       (visitedEdges[jacobiSet[i].first] == false)) {
      // saddle, non-visited edge

      // we have a seed here
      SimplexId sheet1Id = originalData_.sheet1List_.size();
      originalData_.sheet1List_.resize(originalData_.sheet1List_.size() + 1);
      originalData_.sheet1List_.back().hasSaddleEdges_ = false;
      originalData_.sheet1List_.back().pruned_ = false;

      queue<SimplexId> edgeQueue;
      edgeQueue.push(jacobiSet[i].first);

      do {

        SimplexId edgeId = edgeQueue.front();
        edgeQueue.pop();

        if(!visitedEdges[edgeId]) {

          jacobiClassification.push_back(
            pair<SimplexId, SimplexId>(edgeId, sheet1Id));

          if(originalData_.edgeTypes_[edgeId] == 1) {
            // saddle edge
            originalData_.sheet1List_.back().hasSaddleEdges_ = true;
          }
          originalData_.sheet1List_.back().edgeList_.push_back(edgeId);
          originalData_.edge2sheet1_[edgeId] = sheet1Id;

          // next grab its neighbors
          // if saddleEdge and not visited continue (only if one)
          SimplexId vertexId0 = -1;
          triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
          SimplexId vertexId1 = -1;
          triangulation_->getEdgeVertex(edgeId, 1, vertexId1);

          if(originalData_.vertex2sheet0_[vertexId0]
             >= (SimplexId)originalData_.sheet0List_.size()) {
            // WEIRD BUG after multiple runs
            originalData_.vertex2sheet0_[vertexId0] = -1;
          }

          // make sure these are not already visited 0-sheets
          if(originalData_.vertex2sheet0_[vertexId0] == -1) {

            // inspect neighbor edges
            SimplexId neighborNumber = 0;
            SimplexId nextEdgeId = -1;
            SimplexId vertexEdgeNumber
              = triangulation_->getVertexEdgeNumber(vertexId0);
            for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

              SimplexId vertexEdgeId = -1;
              triangulation_->getVertexEdge(vertexId0, j, vertexEdgeId);

              if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
                neighborNumber++;
                if(vertexEdgeId != edgeId
                   /*&&(saddleEdge[(*vertexEdgeList_)[vertexId0][j]] == 1)*/) {
                  nextEdgeId = vertexEdgeId;
                }
              }
            }

            if((neighborNumber == 2) && (nextEdgeId != -1)) {
              // this is a regular vertex and we can add the next edge to the
              // queue if it's not been visited before.
              if(!visitedEdges[nextEdgeId])
                edgeQueue.push(nextEdgeId);
            } else {
              // we hit a 0-sheet that hasn't been visited before.
              // let's create it

              // mark up the segmentation
              originalData_.vertex2sheet0_[vertexId0]
                = originalData_.sheet0List_.size();

              originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                               + 1);
              originalData_.sheet0List_.back().vertexId_ = vertexId0;
              originalData_.sheet0List_.back().type_ = 1;
              originalData_.sheet0List_.back().pruned_ = false;
            }
          }

          if(originalData_.vertex2sheet0_[vertexId0] != -1) {
            // we hit a 0-sheet
            // attach it to the one sheet and reciprocally

            // attach the 0-sheet to the 1-sheet
            originalData_.sheet1List_[sheet1Id].sheet0List_.push_back(
              originalData_.vertex2sheet0_[vertexId0]);

            // attach the 1-sheet to the 0-sheet
            originalData_.sheet0List_[originalData_.vertex2sheet0_[vertexId0]]
              .sheet1List_.push_back(sheet1Id);
          }

          if(originalData_.vertex2sheet0_[vertexId1]
             >= (SimplexId)originalData_.sheet0List_.size()) {
            // WEIRD BUG after multiple runs
            originalData_.vertex2sheet0_[vertexId1] = -1;
          }
          // now do the same for the other vertex
          // make sure these are not already visited 0-sheets
          if(originalData_.vertex2sheet0_[vertexId1] == -1) {

            // inspect neighbor edges
            SimplexId neighborNumber = 0;
            SimplexId nextEdgeId = -1;

            SimplexId vertexEdgeNumber
              = triangulation_->getVertexEdgeNumber(vertexId1);

            for(SimplexId j = 0; j < vertexEdgeNumber; j++) {

              SimplexId vertexEdgeId = -1;
              triangulation_->getVertexEdge(vertexId1, j, vertexEdgeId);

              if(originalData_.edgeTypes_[vertexEdgeId] != -1) {
                neighborNumber++;
                if(vertexEdgeId != edgeId
                   /*&&(saddleEdge[(*vertexEdgeList_)[vertexId1][j]] == 1)*/) {
                  nextEdgeId = vertexEdgeId;
                }
              }
            }

            if((neighborNumber == 2) && (nextEdgeId != -1)) {
              // this is a regular vertex and we can add the next edge to the
              // queue if it's not been visited before.
              if(!visitedEdges[nextEdgeId])
                edgeQueue.push(nextEdgeId);
            } else {
              // we hit a 0-sheet that hasn't been visited before.
              // let's create it

              // mark up the segmentation
              originalData_.vertex2sheet0_[vertexId1]
                = originalData_.sheet0List_.size();

              originalData_.sheet0List_.resize(originalData_.sheet0List_.size()
                                               + 1);
              originalData_.sheet0List_.back().vertexId_ = vertexId1;
              originalData_.sheet0List_.back().type_ = 1;
              originalData_.sheet0List_.back().pruned_ = false;
            }
          }

          if(originalData_.vertex2sheet0_[vertexId1] != -1) {
            // we hit a 0-sheet
            // attach it to the one sheet and reciprocally

            // attach the 0-sheet to the 1-sheet
            originalData_.sheet1List_[sheet1Id].sheet0List_.push_back(
              originalData_.vertex2sheet0_[vertexId1]);

            // attach the 1-sheet to the 0-sheet
            originalData_.sheet0List_[originalData_.vertex2sheet0_[vertexId1]]
              .sheet1List_.push_back(sheet1Id);
          }

          visitedEdges[edgeId] = true;
        }

      } while(edgeQueue.size());
    }
  }

  {
    stringstream msg;
    msg << "[ReebSpace] " << originalData_.sheet1List_.size()
        << " 1-sheets and " << originalData_.sheet0List_.size()
        << " 0-sheets extracted in " << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ReebSpace::compute3sheet(
  const SimplexId &vertexId,
  const vector<vector<vector<SimplexId>>> &tetTriangles) {

  SimplexId sheetId = originalData_.sheet3List_.size();
  originalData_.sheet3List_.resize(originalData_.sheet3List_.size() + 1);
  originalData_.sheet3List_.back().pruned_ = false;
  originalData_.sheet3List_.back().preMerger_ = -1;
  originalData_.sheet3List_.back().Id_ = sheetId;

  queue<SimplexId> vertexQueue;
  vertexQueue.push(vertexId);

  do {

    SimplexId localVertexId = vertexQueue.front();
    vertexQueue.pop();

    if(originalData_.vertex2sheet3_[localVertexId] == -1) {
      // not visited yet

      originalData_.sheet3List_.back().vertexList_.push_back(localVertexId);
      originalData_.vertex2sheet3_[localVertexId] = sheetId;

      SimplexId vertexStarNumber
        = triangulation_->getVertexStarNumber(localVertexId);

      for(SimplexId i = 0; i < vertexStarNumber; i++) {
        SimplexId tetId = -1;
        triangulation_->getVertexStar(localVertexId, i, tetId);

        if(tetTriangles[tetId].empty()) {
          //             if(originalData_.tet2sheet3_[tetId] == -1){
          //               originalData_.tet2sheet3_[tetId] = sheetId;
          //               originalData_.sheet3List_.back().tetList_.push_back(tetId);
          //             }
          for(int j = 0; j < 4; j++) {

            SimplexId tetVertexId = -1;
            triangulation_->getCellVertex(tetId, j, tetVertexId);

            if(originalData_.vertex2sheet3_[tetVertexId] == -1) {
              vertexQueue.push(tetVertexId);
            }
          }
        } else {
          // fiber surface in there
          for(int j = 0; j < 4; j++) {
            SimplexId otherVertexId = -1;
            triangulation_->getCellVertex(tetId, j, otherVertexId);
            if((otherVertexId != localVertexId)
               && (originalData_.vertex2sheet3_[otherVertexId] == -1)) {
              // we need to see if the edge <localVertexId, otherVertexId> is
              // cut by a fiber surface triangle or not.

              bool isCut = false;
              for(SimplexId k = 0; k < (SimplexId)tetTriangles[tetId].size();
                  k++) {
                SimplexId l = 0, m = 0, n = 0;
                l = tetTriangles[tetId][k][0];
                m = tetTriangles[tetId][k][1];
                n = tetTriangles[tetId][k][2];

                for(int p = 0; p < 3; p++) {
                  pair<SimplexId, SimplexId> meshEdge;

                  if(fiberSurfaceVertexList_.size()) {
                    // the fiber surfaces have been merged
                    meshEdge
                      = fiberSurfaceVertexList_[originalData_.sheet2List_[l]
                                                  .triangleList_[m][n]
                                                  .vertexIds_[p]]
                          .meshEdge_;
                  } else {
                    // the fiber surfaces have not been merged
                    meshEdge = originalData_.sheet2List_[l]
                                 .vertexList_[m][originalData_.sheet2List_[l]
                                                   .triangleList_[m][n]
                                                   .vertexIds_[p]]
                                 .meshEdge_;
                  }

                  if(((meshEdge.first == localVertexId)
                      && (meshEdge.second == otherVertexId))
                     || ((meshEdge.second == localVertexId)
                         && (meshEdge.first == otherVertexId))) {
                    isCut = true;
                    break;
                  }
                }

                if(isCut)
                  break;
              }

              if(!isCut) {
                // add the vertex to the queue
                vertexQueue.push(otherVertexId);
              }
            }
          }
        }
      }
    }

  } while(vertexQueue.size());

  return 0;
}

int ReebSpace::compute3sheets(vector<vector<vector<SimplexId>>> &tetTriangles) {

  Timer t;

  tetTriangles.resize(tetNumber_);

  for(SimplexId i = 0; i < (SimplexId)originalData_.sheet2List_.size(); i++) {
    for(SimplexId j = 0;
        j < (SimplexId)originalData_.sheet2List_[i].triangleList_.size(); j++) {
      for(SimplexId k = 0;
          k < (SimplexId)originalData_.sheet2List_[i].triangleList_[j].size();
          k++) {

        SimplexId tetId
          = originalData_.sheet2List_[i].triangleList_[j][k].tetId_;

        vector<SimplexId> triangle(3);
        triangle[0] = i;
        triangle[1] = j;
        triangle[2] = k;

        tetTriangles[tetId].push_back(triangle);
      }
    }
  }

  // mark all the jacobi edge vertices
  for(SimplexId i = 0; i < (SimplexId)originalData_.sheet1List_.size(); i++) {
    for(SimplexId j = 0;
        j < (SimplexId)originalData_.sheet1List_[i].edgeList_.size(); j++) {

      SimplexId edgeId = originalData_.sheet1List_[i].edgeList_[j];

      SimplexId vertexId0 = -1, vertexId1 = -1;
      triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
      triangulation_->getEdgeVertex(edgeId, 1, vertexId1);

      originalData_.vertex2sheet3_[vertexId0] = -2 - i;
      originalData_.vertex2sheet3_[vertexId1] = -2 - i;
    }
  }

  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(originalData_.vertex2sheet3_[i] == -1) {
      compute3sheet(i, tetTriangles);
    }
  }

  // for 3-sheet expansion
  vector<vector<pair<SimplexId, bool>>> neighborList(
    originalData_.sheet3List_.size());
  // end of 3-sheet expansion

  // add the tets in parallel
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)originalData_.sheet3List_.size(); i++) {
    for(SimplexId j = 0;
        j < (SimplexId)originalData_.sheet3List_[i].vertexList_.size(); j++) {

      SimplexId vertexId = originalData_.sheet3List_[i].vertexList_[j];
      SimplexId sheetId = originalData_.vertex2sheet3_[vertexId];

      SimplexId vertexStarNumber
        = triangulation_->getVertexStarNumber(vertexId);

      for(SimplexId k = 0; k < vertexStarNumber; k++) {
        SimplexId tetId = -1;
        triangulation_->getVertexStar(vertexId, k, tetId);
        if(tetTriangles[tetId].empty()) {
          if(originalData_.tet2sheet3_[tetId] == -1) {
            originalData_.tet2sheet3_[tetId] = i;
            originalData_.sheet3List_[i].tetList_.push_back(tetId);
          }
        } else {

          // expending here 3-sheets.
          if(expand3sheets_) {
            for(int l = 0; l < 4; l++) {
              SimplexId otherVertexId = -1;

              triangulation_->getCellVertex(tetId, l, otherVertexId);

              if(vertexId != otherVertexId) {

                SimplexId otherSheetId
                  = originalData_.vertex2sheet3_[otherVertexId];

                if((sheetId != otherSheetId) && (otherSheetId >= 0)) {

                  bool inThere = false;
                  for(SimplexId m = 0;
                      m < (SimplexId)neighborList[sheetId].size(); m++) {
                    if(neighborList[sheetId][m].first == otherSheetId) {
                      inThere = true;
                      break;
                    }
                  }

                  if(!inThere) {
                    neighborList[sheetId].push_back(
                      pair<SimplexId, bool>(otherSheetId, true));
                  }

                  for(SimplexId m = 0;
                      m < (SimplexId)tetTriangles[tetId].size(); m++) {

                    // see if this guy is a saddle
                    SimplexId x = tetTriangles[tetId][m][0];
                    SimplexId y = tetTriangles[tetId][m][1];
                    SimplexId z = tetTriangles[tetId][m][2];

                    FiberSurface::Triangle *tr
                      = &(originalData_.sheet2List_[x].triangleList_[y][z]);

                    bool cuttingTriangle = false;
                    for(int n = 0; n < 3; n++) {
                      FiberSurface::Vertex *v
                        = &(originalData_.sheet2List_[x]
                              .vertexList_[y][tr->vertexIds_[n]]);

                      if(((v->meshEdge_.first == vertexId)
                          && (v->meshEdge_.second == otherVertexId))
                         || ((v->meshEdge_.second == vertexId)
                             && (v->meshEdge_.first == otherVertexId))) {

                        cuttingTriangle = true;
                        break;
                      }
                    }

                    if(cuttingTriangle) {

                      SimplexId polygonId = tr->polygonEdgeId_;
                      SimplexId edgeId = jacobi2edges_[polygonId];
                      if(originalData_.edgeTypes_[edgeId] == 1) {
                        // this is a saddle Jacobi edge

                        inThere = false;
                        for(SimplexId n = 0;
                            n < (SimplexId)neighborList[sheetId].size(); n++) {

                          if(neighborList[sheetId][n].first == otherSheetId) {
                            if(neighborList[sheetId][n].second == true) {
                              neighborList[sheetId][n].second = false;
                            }
                            inThere = true;
                            break;
                          }
                        }

                        if(!inThere) {
                          neighborList[sheetId].push_back(
                            pair<SimplexId, bool>(otherSheetId, false));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          // end of the 3-sheet expansion
        }
      }
    }
  }

  SimplexId totalSheetNumber = 0;

  if(expand3sheets_) {

    totalSheetNumber = originalData_.sheet3List_.size();

    // expending sheets
    for(SimplexId i = 0; i < (SimplexId)originalData_.sheet3List_.size(); i++) {
      if(originalData_.sheet3List_[i].pruned_ == false) {

        for(SimplexId j = 0; j < (SimplexId)neighborList[i].size(); j++) {

          if(neighborList[i][j].second) {

            bool isForbidden = false;
            SimplexId neighborId = neighborList[i][j].first;

            // get the sheet where this guys has been merged
            while(originalData_.sheet3List_[neighborId].preMerger_ != -1) {

              neighborId = originalData_.sheet3List_[neighborId].preMerger_;
            }

            // make sure that no forbidden neighbor has been merged in the
            // candidate neighbor
            for(SimplexId k = 0;
                k < (SimplexId)originalData_.sheet3List_[neighborId]
                      .preMergedSheets_.size();
                k++) {

              SimplexId subNeighborId
                = originalData_.sheet3List_[neighborId].preMergedSheets_[k];

              for(SimplexId l = 0; l < (SimplexId)neighborList[i].size(); l++) {
                if((neighborList[i][l].first == subNeighborId)
                   && (!neighborList[i][l].second)) {
                  isForbidden = true;
                  break;
                }
              }
              if(isForbidden)
                break;
            }

            // make sure that neighborId is not a candidate for a merge with a
            // sheet that is forbidden for i
            if(!isForbidden) {
              for(SimplexId k = 0; k < (SimplexId)neighborList[i].size(); k++) {
                if(!neighborList[i][k].second) {
                  SimplexId forbiddenNeighbor = neighborList[i][k].first;

                  // make sure forbiddenNeighbor is not a valid merger for
                  // neighborId
                  for(SimplexId l = 0;
                      l < (SimplexId)neighborList[neighborId].size(); l++) {
                    if((forbiddenNeighbor == neighborList[neighborId][l].first)
                       && (neighborList[neighborId][l].second)) {
                      isForbidden = true;
                      break;
                    }
                  }
                  if(isForbidden)
                    break;
                }
              }
            }

            if((neighborId != i) && (!isForbidden)
               && (originalData_.sheet3List_[neighborId].pruned_ == false)
               && (originalData_.sheet3List_[neighborId].vertexList_.size()
                   > originalData_.sheet3List_[i].vertexList_.size())) {

              // these can be merged
              preMergeSheets(i, neighborId);
              totalSheetNumber--;
              break;
            }
          }
        }
      }
    }
  } else {
    totalSheetNumber = originalData_.sheet3List_.size();
  }
  // end of the 3-sheet expansion

  {
    stringstream msg;
    msg << "[ReebSpace] " << totalSheetNumber << " 3-sheets computed in "
        << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ReebSpace::connect3sheetTo0sheet(ReebSpaceData &data,
                                     const SimplexId &sheet3Id,
                                     const SimplexId &sheet0Id) {

  bool alreadyConnected = false;
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[sheet3Id].sheet0List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet0List_[i] == sheet0Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet0List_.push_back(sheet0Id);

  for(SimplexId i = 0;
      i < (SimplexId)data.sheet0List_[sheet0Id].sheet3List_.size(); i++) {
    if(data.sheet0List_[sheet0Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet0List_[sheet0Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ReebSpace::connect3sheetTo1sheet(ReebSpaceData &data,
                                     const SimplexId &sheet3Id,
                                     const SimplexId &sheet1Id) {

  bool alreadyConnected = false;
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[sheet3Id].sheet1List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet1List_[i] == sheet1Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet1List_.push_back(sheet1Id);

  for(SimplexId i = 0;
      i < (SimplexId)data.sheet1List_[sheet1Id].sheet3List_.size(); i++) {
    if(data.sheet1List_[sheet1Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet1List_[sheet1Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ReebSpace::connect3sheetTo2sheet(ReebSpaceData &data,
                                     const SimplexId &sheet3Id,
                                     const SimplexId &sheet2Id) {

  bool alreadyConnected = false;
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[sheet3Id].sheet2List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet2List_[i] == sheet2Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet2List_.push_back(sheet2Id);

  for(SimplexId i = 0;
      i < (SimplexId)data.sheet2List_[sheet2Id].sheet3List_.size(); i++) {
    if(data.sheet2List_[sheet2Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet2List_[sheet2Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ReebSpace::connect3sheetTo3sheet(ReebSpaceData &data,
                                     const SimplexId &sheet3Id,
                                     const SimplexId &otherSheet3Id) {

  if(sheet3Id == otherSheet3Id)
    return -1;

  bool alreadyConnected = false;
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[sheet3Id].sheet3List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet3List_[i] == otherSheet3Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet3List_.push_back(otherSheet3Id);

  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[otherSheet3Id].sheet3List_.size(); i++) {
    if(data.sheet3List_[otherSheet3Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -3;
  }
  data.sheet3List_[otherSheet3Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ReebSpace::connectSheets() {

  Timer t;

  for(SimplexId i = 0; i < (SimplexId)originalData_.sheet2List_.size(); i++) {
    for(SimplexId j = 0;
        j < (SimplexId)originalData_.sheet2List_[i].triangleList_.size(); j++) {

      for(SimplexId k = 0;
          k < (SimplexId)originalData_.sheet2List_[i].triangleList_[j].size();
          k++) {

        SimplexId tetId
          = originalData_.sheet2List_[i].triangleList_[j][k].tetId_;

        for(int l = 0; l < 4; l++) {
          SimplexId vertexId = -1;
          triangulation_->getCellVertex(tetId, l, vertexId);

          SimplexId sheet3Id = originalData_.vertex2sheet3_[vertexId];

          if(sheet3Id >= 0) {
            connect3sheetTo2sheet(originalData_, sheet3Id, i);
          }
        }
      }
    }
  }

  // connect 3-sheets together
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(originalData_.vertex2sheet3_[i] >= 0) {

      SimplexId vertexEdgeNumber = triangulation_->getVertexEdgeNumber(i);

      for(SimplexId j = 0; j < vertexEdgeNumber; j++) {
        SimplexId edgeId = -1;
        triangulation_->getVertexEdge(i, j, edgeId);
        SimplexId otherVertexId = -1;
        triangulation_->getEdgeVertex(edgeId, 0, otherVertexId);
        if(otherVertexId == i) {
          triangulation_->getEdgeVertex(edgeId, 1, otherVertexId);
        }

        if(originalData_.vertex2sheet3_[otherVertexId] >= 0) {
          if(originalData_.vertex2sheet3_[otherVertexId]
             != originalData_.vertex2sheet3_[i]) {

            connect3sheetTo3sheet(originalData_,
                                  originalData_.vertex2sheet3_[i],
                                  originalData_.vertex2sheet3_[otherVertexId]);
          }
        }

        if(originalData_.vertex2sheet0_[otherVertexId] != -1) {
          connect3sheetTo0sheet(originalData_, originalData_.vertex2sheet3_[i],
                                originalData_.vertex2sheet0_[otherVertexId]);
        }

        if(originalData_.vertex2sheet3_[otherVertexId] < -1) {
          SimplexId sheet1Id = -2 - originalData_.vertex2sheet3_[otherVertexId];
          connect3sheetTo1sheet(
            originalData_, originalData_.vertex2sheet3_[i], sheet1Id);
        }
      }
    }
  }

  {
    stringstream msg;
    msg << "[ReebSpace] Sheet connectivity established." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  printConnectivity(cout, originalData_);

  hasConnectedSheets_ = true;

  return 0;
}

int ReebSpace::disconnect1sheetFrom0sheet(ReebSpaceData &data,
                                          const SimplexId &sheet1Id,
                                          const SimplexId &sheet0Id,
                                          const SimplexId &biggerId) {

  vector<SimplexId> newList;

  newList.reserve(data.sheet0List_[sheet0Id].sheet1List_.size());

  for(SimplexId i = 0;
      i < (SimplexId)data.sheet0List_[sheet0Id].sheet1List_.size(); i++) {

    if(data.sheet0List_[sheet0Id].sheet1List_[i] != sheet1Id) {
      newList.push_back(data.sheet0List_[sheet0Id].sheet1List_[i]);
    }
  }

  if(data.sheet0List_[sheet0Id].sheet1List_.empty()) {
    data.sheet0List_[sheet0Id].pruned_ = true;
    data.vertex2sheet3_[data.sheet0List_[sheet0Id].vertexId_] = biggerId;
  }

  return 0;
}

int ReebSpace::disconnect3sheetFrom0sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet0Id) {

  vector<SimplexId> newList;

  newList.reserve(data.sheet0List_[sheet0Id].sheet3List_.size());
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet0List_[sheet0Id].sheet3List_.size(); i++) {
    if(data.sheet0List_[sheet0Id].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet0List_[sheet0Id].sheet3List_[i]);
  }

  data.sheet0List_[sheet0Id].sheet3List_ = newList;

  return 0;
}

int ReebSpace::disconnect3sheetFrom1sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet1Id,
                                          const SimplexId &biggerId) {

  vector<SimplexId> newList;

  newList.reserve(data.sheet1List_[sheet1Id].sheet3List_.size());
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet1List_[sheet1Id].sheet3List_.size(); i++) {
    SimplexId other3SheetId = data.sheet1List_[sheet1Id].sheet3List_[i];
    if((other3SheetId != sheet3Id) && (!data.sheet3List_[other3SheetId].pruned_)
       && (data.sheet3List_[other3SheetId].tetList_.size())) {
      newList.push_back(data.sheet1List_[sheet1Id].sheet3List_[i]);
    }
  }

  data.sheet1List_[sheet1Id].sheet3List_ = newList;

  if((data.sheet1List_[sheet1Id].hasSaddleEdges_)
     && (data.sheet1List_[sheet1Id].sheet3List_.size() == 1)) {
    // this guy is no longer separting any body.
    data.sheet1List_[sheet1Id].pruned_ = true;
    data.sheet2List_[sheet1Id].pruned_ = true;

    // update the segmentation
    for(SimplexId i = 0;
        i < (SimplexId)data.sheet1List_[sheet1Id].edgeList_.size(); i++) {

      SimplexId vertexId = -1;
      triangulation_->getEdgeVertex(
        data.sheet1List_[sheet1Id].edgeList_[i], 0, vertexId);
      data.vertex2sheet3_[vertexId] = biggerId;

      vertexId = -1;
      triangulation_->getEdgeVertex(
        data.sheet1List_[sheet1Id].edgeList_[i], 1, vertexId);
      data.vertex2sheet3_[vertexId] = biggerId;
    }

    for(SimplexId i = 0;
        i < (SimplexId)data.sheet1List_[sheet1Id].sheet0List_.size(); i++) {

      disconnect1sheetFrom0sheet(
        data, sheet1Id, data.sheet1List_[sheet1Id].sheet0List_[i], biggerId);
    }
  }

  return 0;
}

int ReebSpace::disconnect3sheetFrom2sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet2Id) {

  vector<SimplexId> newList;

  newList.reserve(data.sheet2List_[sheet2Id].sheet3List_.size());
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet2List_[sheet2Id].sheet3List_.size(); i++) {
    if(data.sheet2List_[sheet2Id].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet2List_[sheet2Id].sheet3List_[i]);
  }

  data.sheet2List_[sheet2Id].sheet3List_ = newList;

  return 0;
}

int ReebSpace::disconnect3sheetFrom3sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &other3SheetId) {

  vector<SimplexId> newList;

  newList.reserve(data.sheet3List_[other3SheetId].sheet3List_.size());
  for(SimplexId i = 0;
      i < (SimplexId)data.sheet3List_[other3SheetId].sheet3List_.size(); i++) {
    if(data.sheet3List_[other3SheetId].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet3List_[other3SheetId].sheet3List_[i]);
  }

  data.sheet3List_[other3SheetId].sheet3List_ = newList;

  return 0;
}

int ReebSpace::flush() {

  totalArea_ = -1;
  totalVolume_ = -1;
  totalHyperVolume_ = -1;
  hasConnectedSheets_ = false;

  // store the segmentation for later purpose
  originalData_.vertex2sheet0_.resize(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++)
    originalData_.vertex2sheet0_[i] = -1;

  originalData_.vertex2sheet3_.resize(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++)
    originalData_.vertex2sheet3_[i] = -1;

  originalData_.edge2sheet1_.resize(triangulation_->getNumberOfEdges(), -1);
  originalData_.edgeTypes_.resize(triangulation_->getNumberOfEdges(), -1);

  originalData_.tet2sheet3_.resize(tetNumber_, -1);

  jacobi2edges_.clear();
  jacobiSetEdges_.clear();

  originalData_.sheet0List_.clear();
  originalData_.sheet1List_.clear();
  originalData_.sheet2List_.clear();
  originalData_.sheet3List_.clear();

  fiberSurfaceVertexList_.clear();

  return 0;
}

int ReebSpace::mergeSheets(const SimplexId &smallerId,
                           const SimplexId &biggerId) {

  // 1. add the vertices and tets of smaller to bigger
  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].vertexList_.size();
      i++) {
    SimplexId vertexId = currentData_.sheet3List_[smallerId].vertexList_[i];
    currentData_.sheet3List_[biggerId].vertexList_.push_back(vertexId);
    currentData_.vertex2sheet3_[vertexId] = biggerId;
  }
  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].tetList_.size(); i++) {
    SimplexId tetId = currentData_.sheet3List_[smallerId].tetList_[i];
    currentData_.sheet3List_[biggerId].tetList_.push_back(tetId);
    currentData_.tet2sheet3_[tetId] = biggerId;
  }

  // 2. update bigger's score and re-insert it in the candidate list
  currentData_.sheet3List_[biggerId].domainVolume_
    += currentData_.sheet3List_[smallerId].domainVolume_;
  currentData_.sheet3List_[biggerId].rangeArea_
    += currentData_.sheet3List_[smallerId].rangeArea_;
  currentData_.sheet3List_[biggerId].hyperVolume_
    += currentData_.sheet3List_[smallerId].hyperVolume_;

  // 3. add smaller's connections to bigger's
  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet3List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet3List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo3sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet2List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet2List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo2sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet1List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet1List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo1sheet(currentData_, biggerId, otherSheetId);
    }
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet0List_.size();
      i++) {

    SimplexId otherSheetId = currentData_.sheet3List_[smallerId].sheet0List_[i];

    if(otherSheetId != biggerId) {
      connect3sheetTo0sheet(currentData_, biggerId, otherSheetId);
    }
  }

  currentData_.sheet3List_[smallerId].pruned_ = true;

  // 4. disconnect smaller from all of its connections.
  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet3List_.size();
      i++) {
    disconnect3sheetFrom3sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet3List_[i]);
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet2List_.size();
      i++) {
    disconnect3sheetFrom2sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet2List_[i]);
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet1List_.size();
      i++) {
    disconnect3sheetFrom1sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet1List_[i], biggerId);
  }

  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[smallerId].sheet0List_.size();
      i++) {
    disconnect3sheetFrom0sheet(
      currentData_, smallerId,
      currentData_.sheet3List_[smallerId].sheet0List_[i]);
  }

  return 0;
}

int ReebSpace::preMergeSheets(const SimplexId &sheetId0,
                              const SimplexId &sheetId1) {

  // 1. add the vertices and tets of 0 to 1
  for(SimplexId i = 0;
      i < (SimplexId)originalData_.sheet3List_[sheetId0].vertexList_.size();
      i++) {
    SimplexId vertexId = originalData_.sheet3List_[sheetId0].vertexList_[i];
    originalData_.sheet3List_[sheetId1].vertexList_.push_back(vertexId);
    originalData_.vertex2sheet3_[vertexId] = sheetId1;
  }
  for(SimplexId i = 0;
      i < (SimplexId)originalData_.sheet3List_[sheetId0].tetList_.size(); i++) {
    SimplexId tetId = originalData_.sheet3List_[sheetId0].tetList_[i];
    originalData_.sheet3List_[sheetId1].tetList_.push_back(tetId);
    originalData_.tet2sheet3_[tetId] = sheetId1;
  }

  // 2. update bigger's score and re-insert it in the candidate list
  originalData_.sheet3List_[sheetId1].domainVolume_
    += originalData_.sheet3List_[sheetId0].domainVolume_;
  originalData_.sheet3List_[sheetId1].rangeArea_
    += originalData_.sheet3List_[sheetId0].rangeArea_;
  originalData_.sheet3List_[sheetId1].hyperVolume_
    += originalData_.sheet3List_[sheetId0].hyperVolume_;

  originalData_.sheet3List_[sheetId0].pruned_ = true;
  originalData_.sheet3List_[sheetId0].preMerger_ = sheetId1;
  originalData_.sheet3List_[sheetId1].preMergedSheets_.push_back(sheetId0);

  return 0;
}

int ReebSpace::prepareSimplification() {

  Timer t;

  //   currentData_ = originalData_;

  // copy parts of the original data to the current one
  // here we don't want to copy the fiber surfaces, this is just too much
  currentData_.tet2sheet3_ = originalData_.tet2sheet3_;
  currentData_.vertex2sheet0_ = originalData_.vertex2sheet0_;
  currentData_.vertex2sheet3_ = originalData_.vertex2sheet3_;
  currentData_.edge2sheet1_ = originalData_.edge2sheet1_;
  currentData_.edgeTypes_ = originalData_.edgeTypes_;

  currentData_.sheet0List_ = originalData_.sheet0List_;
  currentData_.sheet1List_ = originalData_.sheet1List_;
  currentData_.sheet3List_ = originalData_.sheet3List_;

  currentData_.sheet2List_.resize(originalData_.sheet2List_.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet2List_.size(); i++) {
    currentData_.sheet2List_[i].sheet1Id_
      = originalData_.sheet2List_[i].sheet1Id_;
    currentData_.sheet2List_[i].sheet3List_
      = originalData_.sheet2List_[i].sheet3List_;
  }

  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet3List_.size(); i++) {
    currentData_.sheet3List_[i].simplificationId_
      = currentData_.sheet3List_[i].Id_;
  }

  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet1List_.size(); i++) {
    if((currentData_.sheet1List_[i].hasSaddleEdges_)
       && (currentData_.sheet1List_[i].sheet3List_.size() == 1)) {

      currentData_.sheet1List_[i].pruned_ = true;
      currentData_.sheet2List_[i].pruned_ = true;

      for(SimplexId j = 0;
          j < (SimplexId)currentData_.sheet1List_[i].sheet0List_.size(); j++) {
        SimplexId sheet0Id = currentData_.sheet1List_[i].sheet0List_[j];
        currentData_.sheet0List_[sheet0Id].pruned_ = true;
      }
    }
  }

  {
    stringstream msg;
    msg << "[ReebSpace] Data prepared for simplification in "
        << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

int ReebSpace::printConnectivity(ostream &stream,
                                 const ReebSpaceData &data) const {

  if(debugLevel_ < advancedInfoMsg)
    return -1;

  stringstream msg;

  msg << "[ReebSpace] Connectivity..." << endl;

  msg << "[ReebSpace] " << data.sheet0List_.size() << " 0-sheets:" << endl;
  for(SimplexId i = 0; i < (SimplexId)data.sheet0List_.size(); i++) {
    msg << "[ReebSpace]  3-sheets for 0-sheet #" << i
        << " [p=" << data.sheet0List_[i].pruned_ << "]"
        << ": ";
    for(SimplexId j = 0; j < (SimplexId)data.sheet0List_[i].sheet3List_.size();
        j++) {
      msg << "#" << data.sheet0List_[i].sheet3List_[j] << ", ";
    }
    msg << endl;
  }

  msg << "[ReebSpace] " << data.sheet1List_.size() << " 1-sheets:" << endl;
  for(SimplexId i = 0; i < (SimplexId)data.sheet1List_.size(); i++) {
    msg << "[ReebSpace]  3-sheets for 1-sheet #" << i
        << " [p=" << data.sheet1List_[i].pruned_ << "]"
        << ": ";
    for(SimplexId j = 0; j < (SimplexId)data.sheet1List_[i].sheet3List_.size();
        j++) {
      msg << "#" << data.sheet1List_[i].sheet3List_[j] << ", ";
    }
    msg << endl;
  }

  msg << "[ReebSpace] " << data.sheet2List_.size() << " 2-sheets:" << endl;
  for(SimplexId i = 0; i < (SimplexId)data.sheet2List_.size(); i++) {
    msg << "[ReebSpace]  3-sheets for 2-sheet #" << i
        << " [p=" << data.sheet2List_[i].pruned_ << "]"
        << ": ";
    for(SimplexId j = 0; j < (SimplexId)data.sheet2List_[i].sheet3List_.size();
        j++) {
      msg << "#" << data.sheet2List_[i].sheet3List_[j] << ", ";
    }
    msg << endl;
  }

  msg << "[ReebSpace] " << data.sheet3List_.size() << " 3-sheets:" << endl;
  for(SimplexId i = 0; i < (SimplexId)data.sheet3List_.size(); i++) {
    msg << "[ReebSpace]  3-sheets for 3-sheet #" << i
        << " [p=" << data.sheet3List_[i].pruned_ << "]"
        << ": ";
    for(SimplexId j = 0; j < (SimplexId)data.sheet3List_[i].sheet3List_.size();
        j++) {
      msg << "#" << data.sheet3List_[i].sheet3List_[j] << ", ";
    }
    msg << endl;
  }

  dMsg(stream, msg.str(), advancedInfoMsg);
  //   dMsg(stream, msg.str(), timeMsg);

  return 0;
}

int ReebSpace::simplifySheet(const SimplexId &sheetId,
                             const SimplificationCriterion &criterion) {

  SimplexId candidateId = -1;
  double maximumScore = -1;

  // see the adjacent 3-sheets
  for(SimplexId i = 0;
      i < (SimplexId)currentData_.sheet3List_[sheetId].sheet3List_.size();
      i++) {
    SimplexId otherSheetId = currentData_.sheet3List_[sheetId].sheet3List_[i];
    if((!currentData_.sheet3List_[otherSheetId].pruned_)
       && (sheetId != otherSheetId)) {

      double otherScore = 0;

      switch(criterion) {
        case domainVolume:
          otherScore = currentData_.sheet3List_[otherSheetId].domainVolume_;
          break;
        case rangeArea:
          otherScore = currentData_.sheet3List_[otherSheetId].rangeArea_;
          break;
        case hyperVolume:
          otherScore = currentData_.sheet3List_[otherSheetId].hyperVolume_;
          break;
      }

      if((maximumScore < 0) || (otherScore > maximumScore)) {
        candidateId = otherSheetId;
        maximumScore = otherScore;
      }
    }
  }

  // see adjacent 3-sheets along 1-sheets
  //   for(int i = 0;
  //     i < (int) currentData_.sheet3List_[sheetId].sheet1List_.size(); i++){
  //
  //     int sheet1Id = currentData_.sheet3List_[sheetId].sheet1List_[i];
  //
  //     for(int j = 0;
  //       j < (int) currentData_.sheet1List_[sheet1Id].sheet3List_.size();
  //       j++){
  //
  //       int otherSheetId = currentData_.sheet1List_[sheet1Id].sheet3List_[j];
  //       if((otherSheetId != sheetId)
  //         &&(!currentData_.sheet3List_[otherSheetId].pruned_)){
  //
  //         double otherScore = 0;
  //
  //         switch(criterion){
  //           case domainVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].domainVolume_;
  //             break;
  //           case rangeArea:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].rangeArea_;
  //             break;
  //           case hyperVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].hyperVolume_;
  //             break;
  //         }
  //
  //         if((maximumScore < 0)||(otherScore > maximumScore)){
  //           candidateId = otherSheetId;
  //           maximumScore = otherScore;
  //         }
  //       }
  //     }
  //   }

  // see adjacent 3-sheets on 0-sheets
  //   for(int i = 0;
  //     i < (int) currentData_.sheet3List_[sheetId].sheet0List_.size(); i++){
  //
  //     int sheet0Id = currentData_.sheet3List_[sheetId].sheet0List_[i];
  //
  //     for(int j = 0;
  //       j < (int) currentData_.sheet0List_[sheet0Id].sheet3List_.size();
  //       j++){
  //
  //       int otherSheetId = currentData_.sheet0List_[sheet0Id].sheet3List_[j];
  //       if((otherSheetId != sheetId)
  //         &&(!currentData_.sheet3List_[otherSheetId].pruned_)){
  //
  //         double otherScore = 0;
  //
  //         switch(criterion){
  //           case domainVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].domainVolume_;
  //             break;
  //           case rangeArea:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].rangeArea_;
  //             break;
  //           case hyperVolume:
  //             otherScore =
  //               currentData_.sheet3List_[otherSheetId].hyperVolume_;
  //             break;
  //         }
  //
  //         if((maximumScore < 0)||(otherScore > maximumScore)){
  //           candidateId = otherSheetId;
  //           maximumScore = otherScore;
  //         }
  //       }
  //     }
  //   }

  // mark as pruned if so
  if(candidateId != -1) {
    mergeSheets(sheetId, candidateId);
  }

  currentData_.sheet3List_[sheetId].pruned_ = true;

  return 0;
}

int ReebSpace::simplifySheets(
  const double &simplificationThreshold,
  const SimplificationCriterion &simplificationCriterion) {

  Timer t;

  if(!currentData_.sheet3List_.size())
    return -1;

  SimplexId simplifiedSheets = 0;
  double lastThreshold = -1;

  for(SimplexId it = 0; it < (SimplexId)originalData_.sheet3List_.size();
      it++) {

    // do while, avoiding infinite loop

    double minValue = -1;
    SimplexId minId = -1;

    for(SimplexId i = 0; i < (SimplexId)currentData_.sheet3List_.size(); i++) {
      if(!currentData_.sheet3List_[i].pruned_) {

        double value = 0;
        switch(simplificationCriterion) {

          case ReebSpace::domainVolume:
            value = currentData_.sheet3List_[i].domainVolume_ / totalVolume_;
            break;

          case ReebSpace::rangeArea:
            value = currentData_.sheet3List_[i].rangeArea_ / totalArea_;
            break;

          case ReebSpace::hyperVolume:
            value
              = currentData_.sheet3List_[i].hyperVolume_ / totalHyperVolume_;
            break;
        }

        if((minId == -1) || (value < minValue)) {
          minValue = value;
          minId = i;
        }
      }
    }

    if((minId != -1) && (minValue < simplificationThreshold)) {
      simplifySheet(minId, simplificationCriterion);
      simplifiedSheets++;
      lastThreshold = minValue;
    } else {
      break;
    }
  }

  currentData_.simplificationThreshold_ = simplificationThreshold;
  currentData_.simplificationCriterion_ = simplificationCriterion;

  SimplexId simplificationId = 0;
  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet3List_.size(); i++) {
    if(!currentData_.sheet3List_[i].pruned_) {
      currentData_.sheet3List_[i].simplificationId_ = simplificationId;
      simplificationId++;
    }
  }
  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet3List_.size(); i++) {
    if(currentData_.sheet3List_[i].pruned_) {
      // find where it merged
      if(currentData_.sheet3List_[i].vertexList_.size()) {
        SimplexId vertexId = currentData_.sheet3List_[i].vertexList_[0];
        SimplexId sheetId = currentData_.vertex2sheet3_[vertexId];
        if(sheetId != i) {
          currentData_.sheet3List_[i].simplificationId_
            = currentData_.sheet3List_[sheetId].simplificationId_;
        }
      }
    }
  }

  // orphans
  for(SimplexId i = 0; i < (SimplexId)currentData_.sheet1List_.size(); i++) {
    if((!currentData_.sheet1List_[i].pruned_)
       && (currentData_.sheet1List_[i].hasSaddleEdges_)) {
      SimplexId nonSimplified = 0;
      for(SimplexId j = 0;
          j < (SimplexId)currentData_.sheet1List_[i].sheet3List_.size(); j++) {
        SimplexId sheet3Id = currentData_.sheet1List_[i].sheet3List_[j];

        if(!currentData_.sheet3List_[sheet3Id].pruned_) {
          nonSimplified++;
        }
      }

      if(nonSimplified < 2) {
        currentData_.sheet1List_[i].pruned_ = true;
      }
    }
  }

  // TODO: update segmentation for 1-sheets and 0-sheets?...

  printConnectivity(cout, currentData_);

  {
    stringstream msg;
    msg << "[ReebSpace] " << simplifiedSheets << " 3-sheets simplified in "
        << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))"
        << endl;
    if(simplifiedSheets) {
      msg << "[ReebSpace] Last 3-sheet simplified at threshold "
          << lastThreshold << endl;
    }
    msg << "[ReebSpace] " << simplificationId << " 3-sheets left." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

// int ReebSpace::triangulateTetrahedron(const int &tetId,
//   const vector<vector<int> > &triangles,
//   vector<long long int> &outputTets){
//
//   // create a local mesh to avoid large memory allocations in the constrained
//   // triangulation class.
//   vector<double> localPoints;
//   map<int, int> global2local;
//   vector<int> local2global;
//
//   int localTriangleNumber = 0;
//   vector<long long int> localMarkers;
//   vector<long long int> localCells;
//
//   // add the triangles of the tet, in order:
//   // i, j, k
//   // 0, 1, 2
//   // 0, 1, 3
//   // 0, 2, 3
//   // 1, 2, 3
//   for(int i = 0; i < 2; i++){
//     for(int j = i + 1; j < 3; j++){
//       for(int k = j + 1; k < 4; k++){
//
//         int iId = 0;
//         int globalVertexId = tetList_[5*tetId + 1 + i];
//
//         map<int, int>::iterator it = global2local.find(-(globalVertexId +
//         1));
//
//         if(it == global2local.end()){
//           iId = local2global.size();
//           global2local[-(globalVertexId + 1)] = iId;
//           local2global.push_back(-(globalVertexId + 1));
//
//           localPoints.resize(localPoints.size() + 3);
//           localPoints[3*iId] = pointSet_[3*globalVertexId];
//           localPoints[3*iId + 1] = pointSet_[3*globalVertexId + 1];
//           localPoints[3*iId + 2] = pointSet_[3*globalVertexId + 2];
//         }
//         else{
//           iId = it->second;
//         }
//
//         int jId = 0;
//         globalVertexId = tetList_[5*tetId + 1 + j];
//
//         it = global2local.find(-(globalVertexId + 1));
//
//         if(it == global2local.end()){
//           jId = local2global.size();
//           global2local[-(globalVertexId + 1)] = jId;
//           local2global.push_back(-(globalVertexId + 1));
//
//           localPoints.resize(localPoints.size() + 3);
//           localPoints[3*jId] = pointSet_[3*globalVertexId];
//           localPoints[3*jId + 1] = pointSet_[3*globalVertexId + 1];
//           localPoints[3*jId + 2] = pointSet_[3*globalVertexId + 2];
//         }
//         else{
//           jId = it->second;
//         }
//
//         int kId = 0;
//         globalVertexId = tetList_[5*tetId + 1 + k];
//
//         it = global2local.find(-(globalVertexId + 1));
//
//         if(it == global2local.end()){
//           kId = local2global.size();
//           global2local[-(globalVertexId + 1)] = kId;
//           local2global.push_back(-(globalVertexId + 1));
//
//           localPoints.resize(localPoints.size() + 3);
//           localPoints[3*kId] = pointSet_[3*globalVertexId];
//           localPoints[3*kId + 1] = pointSet_[3*globalVertexId + 1];
//           localPoints[3*kId + 2] = pointSet_[3*globalVertexId + 2];
//         }
//         else{
//           kId = it->second;
//         }
//
//         int triangleId = localTriangleNumber;
//         localTriangleNumber++;
//         localCells.resize(localCells.size() + 4);
//         localCells[4*triangleId] = 3;
//         localCells[4*triangleId + 1] = iId;
//         localCells[4*triangleId + 2] = jId;
//         localCells[4*triangleId + 3] = kId;
//         localMarkers.push_back(-1);
//       }
//     }
//   }
//
//   // add the fiber surface triangles
//   for(int i = 0; i < (int) triangles.size(); i++){
//
//     int triangleId = localTriangleNumber;
//     localTriangleNumber++;
//
//     localCells.resize(localCells.size() + 4);
//     localCells[4*triangleId] = 3;
//     localMarkers.push_back(triangles[i][0]);
//
//     for(int j = 0; j < 3; j++){
//       int localId = 0;
//       int globalVertexId =
//         sheet2List_[triangles[i][0]].triangleList_[triangles[i][1]][
//           triangles[i][2]].vertexIds_[j];
//
//       map<int, int>::iterator it = global2local.find(globalVertexId);
//       if(it == global2local.end()){
//
//         // merge fiber surface vertices with tetVertices if needed
//         bool tetVertex = false;
//         for(int k = 0; k < 4; k++){
//           vector<double> tetPoint(3);
//           tetPoint[0] = localPoints[3*k];
//           tetPoint[1] = localPoints[3*k+1];
//           tetPoint[2] = localPoints[3*k+2];
//           double distance = Geometry::distance(
//             fiberSurfaceVertexList_[globalVertexId].p_,
//             tetPoint.data());
//           if(distance < pow10(-FLT_DIG)){
//             localId = k;
//             tetVertex = true;
//             break;
//           }
//         }
//
//         if(!tetVertex){
//           // not found
//           localId = local2global.size();
//           global2local[globalVertexId] = localId;
//           local2global.push_back(globalVertexId);
//
//           localPoints.resize(localPoints.size() + 3);
//           localPoints[3*localId] =
//             fiberSurfaceVertexList_[globalVertexId].p_[0];
//           localPoints[3*localId + 1] =
//             fiberSurfaceVertexList_[globalVertexId].p_[1];
//           localPoints[3*localId + 2] =
//             fiberSurfaceVertexList_[globalVertexId].p_[2];
//         }
//       }
//       else{
//         localId = it->second;
//       }
//
//       localCells[4*triangleId + 1 + j] = localId;
//     }
//   }
//
//
//
//   int ret = 0;
//   {
//     ConstrainedTriangulation cTriangulation;
//     cTriangulation.setDebugLevel(16);
// //     cTriangulation.setDebugLevel(0);
//     cTriangulation.setThreadNumber(1);
//     cTriangulation.setInputVertexNumber(localPoints.size()/3);
//     cTriangulation.setInputPoints(
//       (const double *) localPoints.data());
//     cTriangulation.setInputCellNumber(localTriangleNumber);
//     cTriangulation.setInputCells(
//       (const long long int *) localCells.data());
//     cTriangulation.setInputCellMarkers(
//       (const long long int *) localMarkers.data());
//     cTriangulation.setBoundaryMarker(-1);
//     cTriangulation.setBoundaryRemeshing(true);
//
//     vector<float> outputPoints;
//     cTriangulation.setOutputPoints(&outputPoints);
//     cTriangulation.setOutputCells(&outputTets);
//     ret = cTriangulation.execute();
//   }
//
//   // update with global vertex identifiers
//   for(int i = 0; i < (int) outputTets.size()/5; i++){
//     for(int j = 0; j < 4; j++){
//       outputTets[5*i + 1 + j] = local2global[outputTets[5*i + 1 + j]];
//     }
//   }
//
//   return ret;
// }

// int ReebSpace::triangulateThreeSheets(){
//
//   Timer t;
//
//   vector<bool> inQueue(tetNumber_, false);
//   vector<vector<long long int> > outputTets(tetNumber_);
//   vector<vector<vector<int> > > tetTriangles(tetNumber_);
//   vector<int> tetList;
//
//   for(int i = 0; i < (int) sheet2List_.size(); i++){
//     for(int j = 0; j < (int) sheet2List_[i].triangleList_.size(); j++){
//       for(int k = 0; k < (int) sheet2List_[i].triangleList_[j].size(); k++){
//
//         int tetId = sheet2List_[i].triangleList_[j][k].tetId_;
//
//         if(!inQueue[tetId]){
//           tetList.push_back(tetId);
//           inQueue[tetId] = true;
//         }
//
//         tetTriangles[tetId].resize(
//           tetTriangles[tetId].size() + 1);
//         tetTriangles[tetId].back().resize(3);
//         tetTriangles[tetId].back()[0] = i;
//         tetTriangles[tetId].back()[1] = j;
//         tetTriangles[tetId].back()[2] = k;
//       }
//     }
//   }
//
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif
//   for(int i = 0; i < (int) tetList.size(); i++){
//     int tetId = tetList[i];
//     int ret =
//       triangulateTetrahedron(tetId, tetTriangles[tetId], outputTets[tetId]);
//     if(ret < 0){
//       break;
//     }
//   }
//
//   // now merge the tet lists into 1 big list
//   sheet3points_.resize(3*(vertexNumber_ + fiberSurfaceVertexList_.size()));
//   for(int i = 0; i < 3*vertexNumber_; i++){
//     sheet3points_[i] = pointSet_[i];
//   }
//   for(int i = 0; i < (int) fiberSurfaceVertexList_.size(); i++){
//     sheet3points_[3*vertexNumber_ + 3*i] = fiberSurfaceVertexList_[i].p_[0];
//     sheet3points_[3*vertexNumber_ + 3*i + 1] =
//     fiberSurfaceVertexList_[i].p_[1]; sheet3points_[3*vertexNumber_ + 3*i +
//     2] = fiberSurfaceVertexList_[i].p_[2];
//   }
//
//   int tetOffset = 0;
//   for(int i = 0; i < (int) outputTets.size(); i++){
//
//     tetOffset = sheet3cells_.size();
//
//     if(outputTets[i].empty()){
//       // original tet which has not been remeshed
//       sheet3cells_.resize(sheet3cells_.size() + 5);
//       sheet3cells_[tetOffset] = 4;
//       for(int j = 0; j < 4; j++)
//         sheet3cells_[tetOffset + 1 + j] = tetList_[5*i + 1 + j];
//     }
//     else{
//       // concat sheet3cells_ and outputTets[i]
//       sheet3cells_.resize(sheet3cells_.size() + outputTets[i].size());
//       for(int j = 0; j < (int) outputTets[i].size()/5; j++){
//
//         sheet3cells_[tetOffset] = 4;
//
//         int vertexId = -1;
//
//         for(int k = 0; k < 4; k++){
//           if(outputTets[i][5*j + 1 + k] < 0){
//             // tetVertex --> original id
//             vertexId = -(outputTets[i][5*j + 1 + k] + 1);
//           }
//           else{
//             vertexId = vertexNumber_ + outputTets[i][5*j + 1 + k];
//           }
//           sheet3cells_[tetOffset + 1 + k] = vertexId;
//         }
//
//         tetOffset += 5;
//       }
//     }
//   }
//
//   {
//     stringstream msg;
//     msg << "[ReebSpace] 3-sheets ("
//       << sheet3cells_.size()/5
//       << " tets) triangulated in "
//       << t.getElapsedTime() << " s. ("
//       << threadNumber_
//       << " thread(s))" << endl;
//     dMsg(cout, msg.str(), timeMsg);
//   }
//
//   return 0;
// }
