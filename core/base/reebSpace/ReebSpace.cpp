#include <ReebSpace.h>

ttk::ReebSpace::ReebSpace() {
  this->setDebugMsgPrefix("ReebSpace");
}

int ttk::ReebSpace::connect3sheetTo0sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet0Id) {

  bool alreadyConnected = false;
  for(size_t i = 0; i < data.sheet3List_[sheet3Id].sheet0List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet0List_[i] == sheet0Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet0List_.push_back(sheet0Id);

  for(size_t i = 0; i < data.sheet0List_[sheet0Id].sheet3List_.size(); i++) {
    if(data.sheet0List_[sheet0Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet0List_[sheet0Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ttk::ReebSpace::connect3sheetTo1sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet1Id) {

  bool alreadyConnected = false;
  for(size_t i = 0; i < data.sheet3List_[sheet3Id].sheet1List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet1List_[i] == sheet1Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet1List_.push_back(sheet1Id);

  for(size_t i = 0; i < data.sheet1List_[sheet1Id].sheet3List_.size(); i++) {
    if(data.sheet1List_[sheet1Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet1List_[sheet1Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ttk::ReebSpace::connect3sheetTo2sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &sheet2Id) {

  bool alreadyConnected = false;
  for(size_t i = 0; i < data.sheet3List_[sheet3Id].sheet2List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet2List_[i] == sheet2Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet2List_.push_back(sheet2Id);

  for(size_t i = 0; i < data.sheet2List_[sheet2Id].sheet3List_.size(); i++) {
    if(data.sheet2List_[sheet2Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -1;
  }
  data.sheet2List_[sheet2Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ttk::ReebSpace::connect3sheetTo3sheet(ReebSpaceData &data,
                                          const SimplexId &sheet3Id,
                                          const SimplexId &otherSheet3Id) {

  if(sheet3Id == otherSheet3Id)
    return -1;

  bool alreadyConnected = false;
  for(size_t i = 0; i < data.sheet3List_[sheet3Id].sheet3List_.size(); i++) {
    if(data.sheet3List_[sheet3Id].sheet3List_[i] == otherSheet3Id) {
      // already connected
      alreadyConnected = true;
      break;
    }
  }
  if(!alreadyConnected)
    data.sheet3List_[sheet3Id].sheet3List_.push_back(otherSheet3Id);

  for(size_t i = 0; i < data.sheet3List_[otherSheet3Id].sheet3List_.size();
      i++) {
    if(data.sheet3List_[otherSheet3Id].sheet3List_[i] == sheet3Id)
      // already connected
      return -3;
  }
  data.sheet3List_[otherSheet3Id].sheet3List_.push_back(sheet3Id);

  return 0;
}

int ttk::ReebSpace::disconnect1sheetFrom0sheet(ReebSpaceData &data,
                                               const SimplexId &sheet1Id,
                                               const SimplexId &sheet0Id,
                                               const SimplexId &biggerId) {

  std::vector<SimplexId> newList;

  newList.reserve(data.sheet0List_[sheet0Id].sheet1List_.size());

  for(size_t i = 0; i < data.sheet0List_[sheet0Id].sheet1List_.size(); i++) {

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

int ttk::ReebSpace::disconnect3sheetFrom0sheet(ReebSpaceData &data,
                                               const SimplexId &sheet3Id,
                                               const SimplexId &sheet0Id) {

  std::vector<SimplexId> newList;

  newList.reserve(data.sheet0List_[sheet0Id].sheet3List_.size());
  for(size_t i = 0; i < data.sheet0List_[sheet0Id].sheet3List_.size(); i++) {
    if(data.sheet0List_[sheet0Id].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet0List_[sheet0Id].sheet3List_[i]);
  }

  data.sheet0List_[sheet0Id].sheet3List_ = newList;

  return 0;
}

int ttk::ReebSpace::disconnect3sheetFrom2sheet(ReebSpaceData &data,
                                               const SimplexId &sheet3Id,
                                               const SimplexId &sheet2Id) {

  std::vector<SimplexId> newList;

  newList.reserve(data.sheet2List_[sheet2Id].sheet3List_.size());
  for(size_t i = 0; i < data.sheet2List_[sheet2Id].sheet3List_.size(); i++) {
    if(data.sheet2List_[sheet2Id].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet2List_[sheet2Id].sheet3List_[i]);
  }

  data.sheet2List_[sheet2Id].sheet3List_ = newList;

  return 0;
}

int ttk::ReebSpace::disconnect3sheetFrom3sheet(ReebSpaceData &data,
                                               const SimplexId &sheet3Id,
                                               const SimplexId &other3SheetId) {

  std::vector<SimplexId> newList;

  newList.reserve(data.sheet3List_[other3SheetId].sheet3List_.size());
  for(size_t i = 0; i < data.sheet3List_[other3SheetId].sheet3List_.size();
      i++) {
    if(data.sheet3List_[other3SheetId].sheet3List_[i] != sheet3Id)
      newList.push_back(data.sheet3List_[other3SheetId].sheet3List_[i]);
  }

  data.sheet3List_[other3SheetId].sheet3List_ = newList;

  return 0;
}

int ttk::ReebSpace::preMergeSheets(const SimplexId &sheetId0,
                                   const SimplexId &sheetId1) {

  // 1. add the vertices and tets of 0 to 1
  for(size_t i = 0; i < originalData_.sheet3List_[sheetId0].vertexList_.size();
      i++) {
    SimplexId vertexId = originalData_.sheet3List_[sheetId0].vertexList_[i];
    originalData_.sheet3List_[sheetId1].vertexList_.push_back(vertexId);
    originalData_.vertex2sheet3_[vertexId] = sheetId1;
  }
  for(size_t i = 0; i < originalData_.sheet3List_[sheetId0].tetList_.size();
      i++) {
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

int ttk::ReebSpace::prepareSimplification() {

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
  for(size_t i = 0; i < currentData_.sheet2List_.size(); i++) {
    currentData_.sheet2List_[i].sheet1Id_
      = originalData_.sheet2List_[i].sheet1Id_;
    currentData_.sheet2List_[i].sheet3List_
      = originalData_.sheet2List_[i].sheet3List_;
  }

  for(size_t i = 0; i < currentData_.sheet3List_.size(); i++) {
    currentData_.sheet3List_[i].simplificationId_
      = currentData_.sheet3List_[i].Id_;
  }

  for(size_t i = 0; i < currentData_.sheet1List_.size(); i++) {
    if((currentData_.sheet1List_[i].hasSaddleEdges_)
       && (currentData_.sheet1List_[i].sheet3List_.size() == 1)) {

      currentData_.sheet1List_[i].pruned_ = true;
      currentData_.sheet2List_[i].pruned_ = true;

      for(size_t j = 0; j < currentData_.sheet1List_[i].sheet0List_.size();
          j++) {
        SimplexId sheet0Id = currentData_.sheet1List_[i].sheet0List_[j];
        currentData_.sheet0List_[sheet0Id].pruned_ = true;
      }
    }
  }

  this->printMsg("Data prepared for simplification.", 1.0, t.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

int ttk::ReebSpace::printConnectivity(const ReebSpaceData &data) const {

  if(debugLevel_ < static_cast<int>(debug::Priority::DETAIL))
    return -1;

  std::stringstream msg;

  msg << "Connectivity..." << std::endl;

  msg << data.sheet0List_.size() << " 0-sheets:" << std::endl;
  for(size_t i = 0; i < data.sheet0List_.size(); i++) {
    msg << "3-sheets for 0-sheet #" << i
        << " [p=" << data.sheet0List_[i].pruned_ << "]"
        << ": ";
    for(size_t j = 0; j < data.sheet0List_[i].sheet3List_.size(); j++) {
      msg << "#" << data.sheet0List_[i].sheet3List_[j] << ", ";
    }
    msg << std::endl;
  }

  msg << data.sheet1List_.size() << " 1-sheets:" << std::endl;
  for(size_t i = 0; i < data.sheet1List_.size(); i++) {
    msg << "3-sheets for 1-sheet #" << i
        << " [p=" << data.sheet1List_[i].pruned_ << "]"
        << ": ";
    for(size_t j = 0; j < data.sheet1List_[i].sheet3List_.size(); j++) {
      msg << "#" << data.sheet1List_[i].sheet3List_[j] << ", ";
    }
    msg << std::endl;
  }

  msg << data.sheet2List_.size() << " 2-sheets:" << std::endl;
  for(size_t i = 0; i < data.sheet2List_.size(); i++) {
    msg << "3-sheets for 2-sheet #" << i
        << " [p=" << data.sheet2List_[i].pruned_ << "]"
        << ": ";
    for(size_t j = 0; j < data.sheet2List_[i].sheet3List_.size(); j++) {
      msg << "#" << data.sheet2List_[i].sheet3List_[j] << ", ";
    }
    msg << std::endl;
  }

  msg << data.sheet3List_.size() << " 3-sheets:" << std::endl;
  for(size_t i = 0; i < data.sheet3List_.size(); i++) {
    msg << "3-sheets for 3-sheet #" << i
        << " [p=" << data.sheet3List_[i].pruned_ << "]"
        << ": ";
    for(size_t j = 0; j < data.sheet3List_[i].sheet3List_.size(); j++) {
      msg << "#" << data.sheet3List_[i].sheet3List_[j] << ", ";
    }
    msg << std::endl;
  }

  std::string one_line{};
  while(std::getline(msg, one_line)) {
    this->printMsg(one_line, debug::Priority::VERBOSE);
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
//           if(distance < Geometry::powIntTen(-FLT_DIG)){
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
