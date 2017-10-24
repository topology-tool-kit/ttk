#include                  <JacobiSet.h>

template <class dataTypeU, class dataTypeV> 
  JacobiSet<dataTypeU, dataTypeV>::JacobiSet(){

  vertexNumber_ = 0;
  
  uField_ = NULL;
  vField_ = NULL;
  
  tetList_ = NULL;
  
  edgeList_ = NULL;
  edgeFanLinkEdgeLists_ = NULL;
  edgeFans_ = NULL;
  sosOffsetsU_ = NULL;
  sosOffsetsV_ = NULL;

  triangulation_ = NULL;  
}

template <class dataTypeU, class dataTypeV> 
  JacobiSet<dataTypeU, dataTypeV>::~JacobiSet(){
  
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::connectivityPreprocessing(
    const vector<vector<int> > &edgeStarList,
    vector<vector<pair<int, int> > > &edgeFanLinkEdgeLists,
    vector<vector<long long int> > &edgeFans,
    vector<int> &sosOffsets) const{

  Timer t;
 
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!vertexNumber_)
    return -1;
  if(!tetList_)
    return -2;
  if(!edgeList_)
    return -3;
  if(edgeStarList.size() != edgeList_->size())
    return -4;
#endif

  edgeFanLinkEdgeLists.resize(edgeList_->size());
  
  int count = 0;

  // edge triangle fans [each thread writes in a different spot]
  //    for each edge
  //      for each triangle 
  //        list of vertices
  edgeFans.resize(edgeList_->size());
  // pre-allocate memory
  for(int i = 0; i < (int) edgeFans.size(); i++){
    // we store 4 integers per triangle per edge
    edgeFans[i].resize(edgeStarList[i].size()*4);
  }
 
  if(!sosOffsets.size()){
    sosOffsets.resize(vertexNumber_);
    for(int i = 0; i < vertexNumber_; i++){
      sosOffsets[i] = i;
    }
  }
  
  vector<ZeroSkeleton> threadedLinkers(threadNumber_);
  vector<vector<long long int> > threadedLinks(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedLinkers[i].setDebugLevel(debugLevel_);
    threadedLinkers[i].setThreadNumber(1);
  }
  
  vector<OneSkeleton> threadedEdgeListers(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedEdgeListers[i].setDebugLevel(debugLevel_);
    threadedEdgeListers[i].setThreadNumber(1);
  }
 
  vector<vector<int> > threadedTriangleIds(threadNumber_);
  for(int i = 0; i < (int) threadedTriangleIds.size(); i++){
    threadedTriangleIds[i].resize(4);
    threadedTriangleIds[i][0] = 3;
  }
  
#ifdef TTK_ENABLE_OPENMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(int i = 0; i < (int) edgeList_->size(); i++){

    // avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      int threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif
      
      // processing here!
      int pivotVertexId = (*edgeList_)[i].first;
      int otherExtremityId = (*edgeList_)[i].second;
      
      // A) compute triangle fans
      // format: #vertices, id0, id1, id2, etc.
      for(int j = 0; j < (int) edgeStarList[i].size(); j++){
        
        int tetId = edgeStarList[i][j];
        
        // loop over the tet's triangles and add that to the list
        // no need for check, only one triangle verifies this and two tets
        // can't add the same triangle
        for(int k = 0; k < 4; k++){
          
          bool hasPivotVertex = false;
          bool hasOtherExtremity = false;
          for(int l = 0; l < 3; l++){
            threadedTriangleIds[threadId][l + 1] 
              = tetList_[5*tetId + 1 + (l + k)%4];
            if(threadedTriangleIds[threadId][l + 1] == pivotVertexId){
              hasPivotVertex = true;
            }
            if(threadedTriangleIds[threadId][l + 1] == otherExtremityId){
              hasOtherExtremity = true;
            }
          }
          
          if((hasPivotVertex)&&(!hasOtherExtremity)){
            for(int l = 0; l < (int) threadedTriangleIds[threadId].size(); l++){
              edgeFans[i][j*4 + l] = 
                threadedTriangleIds[threadId][l];
            }
            break;
          }
        }
      }
      
      // set-up the link of the edge fan
      threadedLinkers[threadId].buildVertexLink(
        pivotVertexId, edgeFans[i].size()/4, 
          edgeFans[i].data(), threadedLinks[threadId]);
      
      // now compute the edge list of the link
      threadedEdgeListers[threadId].buildEdgeSubList(
        edgeFans[i].size()/4, edgeFans[i].data(),
        edgeFanLinkEdgeLists[i]);
      
      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
#ifdef TTK_ENABLE_OPENMP
        omp_set_lock(&writeLock);
#endif
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }

        count++;
#ifdef TTK_ENABLE_OPENMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }
   
#ifdef TTK_ENABLE_OPENMP
  omp_destroy_lock(&writeLock);
#endif
  
  {
    stringstream msg;
    msg << "[JacobiSet] Edge-fans computed in "
      << t.getElapsedTime() << " s. ("
      << edgeList_->size()
      << " edges)" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::execute(
    vector<pair<int, char> > &jacobiSet){

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!triangulation_)||(triangulation_->isEmpty())){
    if(vertexNumber_){
      return executeLegacy(jacobiSet);
    }
    return -1;
  }
  if(!uField_)
    return -2;
  if(!vField_)
    return -3;
#endif

  int vertexNumber = triangulation_->getNumberOfVertices();
  
  if(!sosOffsetsU_){
    // let's use our own local copy
    sosOffsetsU_ = &localSosOffsetsU_;
  }

  if(vertexNumber != (int) sosOffsetsU_->size()){
    
    sosOffsetsU_->resize(vertexNumber);
    for(int i = 0; i < vertexNumber; i++){
      (*sosOffsetsU_)[i] = i;
    }
  }
  
  if(!sosOffsetsV_){
    // let's use our own local copy
    sosOffsetsV_ = &localSosOffsetsV_;
  }

  if(vertexNumber != (int) sosOffsetsV_->size()){
    
    sosOffsetsV_->resize(vertexNumber);
    for(int i = 0; i < vertexNumber; i++){
      (*sosOffsetsV_)[i] = vertexNumber - i;
    }
  }
  
  jacobiSet.clear();
 
  int edgeNumber = triangulation_->getNumberOfEdges();
  
  vector<vector<pair<int, char> > > threadedCriticalTypes(threadNumber_);
  
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < edgeNumber; i++){
    
    char type = getCriticalType(i);
    
    if(type != -2){
      // -2: regular vertex
      int threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif
      threadedCriticalTypes[threadId].push_back(pair<int, char>(i, type));
    }
  }
  
  // now merge the threaded lists
  for(int i = 0; i < threadNumber_; i++){
    for(int j = 0; j < (int) threadedCriticalTypes[i].size(); j++){
      jacobiSet.push_back(threadedCriticalTypes[i][j]);
    }
  }
  
  if(debugLevel_ >= Debug::infoMsg){
    int minimumNumber = 0, saddleNumber = 0, maximumNumber = 0, 
      monkeySaddleNumber = 0;
      
    for(int i = 0; i < (int) jacobiSet.size(); i++){
      switch(jacobiSet[i].second){
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
    
    {
      stringstream msg;
      msg << "[JacobiSet] Minimum edges: " << minimumNumber << endl;
      msg << "[JacobiSet] Saddle edges: " << saddleNumber << endl;
      msg << "[JacobiSet] Maximum edges: " << maximumNumber << endl;
      msg << "[JacobiSet] Multi-saddle edges: " << monkeySaddleNumber << endl;
      dMsg(cout, msg.str(), Debug::infoMsg);
    }
  }
  
  {
    stringstream msg;
    msg << "[JacobiSet] Data-set (" << edgeNumber
      << " edges) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    msg << "[JacobiSet] Jacobi edge rate: "
      << 100*(jacobiSet.size()/((double) edgeNumber))
      << "%" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::executeLegacy(
    vector<pair<int, char> > &jacobiSet){

  Timer t;
 
  {
    stringstream msg;
    msg << "[JacobiSet] Using legacy implementation..." << endl;
    dMsg(cout, msg.str(), Debug::infoMsg);
  }
  
  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if(!vertexNumber_)
    return -1;
  if(!uField_)
    return -2;
  if(!vField_)
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

  int count = 0;
  
  jacobiSet.clear();
  
  dataTypeU *uField = (dataTypeU *) uField_;
  dataTypeV *vField = (dataTypeV *) vField_;
 
  // distance fields (not really memory efficient)
  // for each thread
  //      for each vertex: distance field map
  vector<vector<double> > threadedDistanceField(threadNumber_);
  for(int i = 0; i < (int) threadedDistanceField.size(); i++){
    threadedDistanceField[i].resize(vertexNumber_);
  }
  
  vector<ScalarFieldCriticalPoints<double> > 
    threadedCriticalPoints(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedCriticalPoints[i].setDomainDimension(2);
    threadedCriticalPoints[i].setScalarValues(
      threadedDistanceField[i].data());
    threadedCriticalPoints[i].setVertexNumber(vertexNumber_);
    threadedCriticalPoints[i].setSosOffsets(sosOffsetsU_);
  }

  vector<vector<pair<int, char> > > threadedCriticalTypes(threadNumber_);

#ifdef TTK_ENABLE_OPENMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(int i = 0; i < (int) edgeList_->size(); i++){

    // avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      int threadId = 0;
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif
      
      // processing here!
      int pivotVertexId = (*edgeList_)[i].first;
      int otherExtremityId = (*edgeList_)[i].second;
     
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
      
      for(int j = 0; j < (int) (*edgeFans_)[i].size()/4; j++){
        for(int k = 0; k < 3; k++){
          
          int vertexId = (*edgeFans_)[i][j*4 + 1 + k];
          
          // we can compute the distance field (in the rage)
          double projectedVertex[2];
          projectedVertex[0] = uField[vertexId];
          projectedVertex[1] = vField[vertexId];
        
          double vertexRangeEdge[2];
          vertexRangeEdge[0] = projectedVertex[0] - projectedPivotVertex[0];
          vertexRangeEdge[1] = projectedVertex[1] - projectedPivotVertex[1];
          
          // signed distance: linear function of the dot product
          threadedDistanceField[threadId][vertexId] = 
            vertexRangeEdge[0]*rangeNormal[0] 
              + vertexRangeEdge[1]*rangeNormal[1];
        }
      }
      
      // B) compute critical points
      // watch out between local and global Ids
      // what I could do is to translate the ids from global to local
      // also, lots of things in there can be done out of the loop
      
      // in the loop
      char type = 
        threadedCriticalPoints[threadId].getCriticalType(pivotVertexId,
          (*edgeFanLinkEdgeLists_)[i]);
        
      if(type != -2){
        // -2: regular vertex
        threadedCriticalTypes[threadId].push_back(pair<int, char>(i, type));
      }
      
      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
#ifdef TTK_ENABLE_OPENMP
        omp_set_lock(&writeLock);
#endif
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }

        count++;
#ifdef TTK_ENABLE_OPENMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }
  
  // now merge the threaded lists
  for(int i = 0; i < threadNumber_; i++){
    for(int j = 0; j < (int) threadedCriticalTypes[i].size(); j++){
      jacobiSet.push_back(threadedCriticalTypes[i][j]);
    }
  }
  
#ifdef TTK_ENABLE_OPENMP
  omp_destroy_lock(&writeLock);
#endif
 
  if(debugLevel_ >= Debug::infoMsg){
    int minimumNumber = 0, saddleNumber = 0, maximumNumber = 0, 
      monkeySaddleNumber = 0;
      
    for(int i = 0; i < (int) jacobiSet.size(); i++){
      switch(jacobiSet[i].second){
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
    
    {
      stringstream msg;
      msg << "[JacobiSet] Minimum edges: " << minimumNumber << endl;
      msg << "[JacobiSet] Saddle edges: " << saddleNumber << endl;
      msg << "[JacobiSet] Maximum edges: " << maximumNumber << endl;
      msg << "[JacobiSet] Multi-saddle edges: " << monkeySaddleNumber << endl;
      dMsg(cout, msg.str(), Debug::infoMsg);
    }
  }
  
  {
    stringstream msg;
    msg << "[JacobiSet] Data-set (" << edgeList_->size()
      << " edges) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    msg << "[JacobiSet] Jacobi edge rate: "
      << 100*(jacobiSet.size()/((double) edgeList_->size()))
      << "%" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV> 
  char JacobiSet<dataTypeU, dataTypeV>::getCriticalType(const int &edgeId){
  
  dataTypeU *uField = (dataTypeU *) uField_;
  dataTypeV *vField = (dataTypeV *) vField_;
  
  int vertexId0 = -1, vertexId1 = -1;
  triangulation_->getEdgeVertex(edgeId, 0, vertexId0);
  triangulation_->getEdgeVertex(edgeId, 1, vertexId1);
  
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
  
  int starNumber = triangulation_->getEdgeStarNumber(edgeId);
  vector<int> lowerNeighbors, upperNeighbors;
  
  int neighborNumber = 0;
  
  for(int i = 0; i < starNumber; i++){
    
    int tetId = -1;
    triangulation_->getEdgeStar(edgeId, i, tetId);
    
    int vertexNumber = triangulation_->getCellVertexNumber(tetId);
    for(int j = 0; j < vertexNumber; j++){
      int vertexId = -1;
      triangulation_->getCellVertex(tetId, j, vertexId);
      
      if((vertexId != -1)&&(vertexId != vertexId0)&&(vertexId != vertexId1)){
        // new neighbor
        bool isIn = false;
        for(int k = 0; k < (int) lowerNeighbors.size(); k++){
          if(vertexId == lowerNeighbors[k]){
            isIn = true;
            break;
          }
        }
        
        if(!isIn){
          for(int k = 0; k < (int) upperNeighbors.size(); k++){
            if(vertexId == upperNeighbors[k]){
              isIn = true;
              break;
            }
          }
        }
        
        if(!isIn){
          // compute the actual distance field
          // A) compute the distance field
          double projectedVertex[2];
          projectedVertex[0] = uField[vertexId];
          projectedVertex[1] = vField[vertexId];
        
          double vertexRangeEdge[2];
          vertexRangeEdge[0] = projectedVertex[0] - projectedPivotVertex[0];
          vertexRangeEdge[1] = projectedVertex[1] - projectedPivotVertex[1];
          
          // signed distance: linear function of the dot product
          double distance = 
            vertexRangeEdge[0]*rangeNormal[0] 
              + vertexRangeEdge[1]*rangeNormal[1];
          
          neighborNumber++;
              
          if(distance < 0){
            lowerNeighbors.push_back(vertexId);
          }
          else if(distance > 0){
            upperNeighbors.push_back(vertexId);
          }
          else{
            // degenerate
            // compute the distance field out of the offset positions
            double offsetProjectedPivotVertex[2];
            offsetProjectedPivotVertex[0] = (*sosOffsetsU_)[vertexId0];
            offsetProjectedPivotVertex[1] = 
              (*sosOffsetsV_)[vertexId0]*(*sosOffsetsV_)[vertexId0];
              
            double offsetProjectedOtherVertex[2];
            offsetProjectedOtherVertex[0] = (*sosOffsetsU_)[vertexId1];
            offsetProjectedOtherVertex[1] = 
              (*sosOffsetsV_)[vertexId1]*(*sosOffsetsV_)[vertexId1];
  
            double offsetRangeEdge[2];
            offsetRangeEdge[0] = 
              offsetProjectedOtherVertex[0] - offsetProjectedPivotVertex[0];
            offsetRangeEdge[1] =
              offsetProjectedOtherVertex[1] - offsetProjectedPivotVertex[1];
              
            double offsetRangeNormal[2];
            offsetRangeNormal[0] = -offsetRangeEdge[1];
            offsetRangeNormal[1] = offsetRangeEdge[0];
            
            projectedVertex[0] = (*sosOffsetsU_)[vertexId];
            projectedVertex[1] = 
              (*sosOffsetsV_)[vertexId]*(*sosOffsetsV_)[vertexId];
              
            vertexRangeEdge[0] = 
              projectedVertex[0] - offsetProjectedPivotVertex[0];
            vertexRangeEdge[1] = 
              projectedVertex[1] - offsetProjectedPivotVertex[1];
                
            distance = 
              vertexRangeEdge[0]*offsetRangeNormal[0] 
              + vertexRangeEdge[1]*offsetRangeNormal[1];
              
            if(distance < 0){
              lowerNeighbors.push_back(vertexId);
            }
            else if (distance > 0){
              upperNeighbors.push_back(vertexId);
            }
            else{
              stringstream msg;
              msg << 
                "[JacobiSet] Inconsistent (non-bijective?) offsets for vertex #"
                << vertexId << endl;
              dMsg(cerr, msg.str(), Debug::infoMsg);
            }
          }
        }
      }
    }
  }
  
  // at this point, we know if each vertex of the edge link is higher or not.
  if((int) (lowerNeighbors.size() + upperNeighbors.size()) != neighborNumber){
    // Inconsistent offsets (cf above error message)
    return -2;
  }
  
  if(lowerNeighbors.empty()){
    // minimum
    return 0;
  }
  if(upperNeighbors.empty()){
    // maximum
    return 2;
  }
  
  // let's check the connectivity now
  vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  vector<UnionFind *> lowerList(lowerNeighbors.size());
  vector<UnionFind> upperSeeds(upperNeighbors.size());
  vector<UnionFind *> upperList(upperNeighbors.size());
  
  for(int i = 0; i < (int) lowerSeeds.size(); i++){
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(int i = 0; i < (int) upperSeeds.size(); i++){
    upperList[i] = &(upperSeeds[i]);
  }
  
  for(int i = 0; i < starNumber; i++){
    
    int tetId = -1;
    triangulation_->getEdgeStar(edgeId, i, tetId);
    
    int vertexNumber = triangulation_->getCellVertexNumber(tetId);
    for(int j = 0; j < vertexNumber; j++){
      int edgeVertexId0 = -1;
      triangulation_->getCellVertex(tetId, j, edgeVertexId0);
      if((edgeVertexId0 != vertexId0)&&(edgeVertexId0 != vertexId1)){
        for(int k = j + 1; k < vertexNumber; k++){
          int edgeVertexId1 = -1;
          triangulation_->getCellVertex(tetId, k, edgeVertexId1);
          if((edgeVertexId1 != vertexId0)&&(edgeVertexId1 != vertexId1)){
            // processing the edge (edgeVertexId0, edgeVertexId1)
            
            // we need to find out if they're lower or not
            bool lower0 = false;
            for(int l = 0; l < (int) lowerNeighbors.size(); l++){
              if(lowerNeighbors[l] == edgeVertexId0){
                lower0 = true;
                break;
              }
            }
            bool lower1 = false;
            for(int l = 0; l < (int) lowerNeighbors.size(); l++){
              if(lowerNeighbors[l] == edgeVertexId1){
                lower1 = true;
                break;
              }
            }
            
            vector<int> *neighbors = &lowerNeighbors;
            vector<UnionFind *> *seeds = &lowerList;
            
            if(!lower0){
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }
            
            if(lower0 == lower1){
              // connect their union-find sets!
              int lowerId0 = -1, lowerId1 = -1;
              for(int l = 0; l < (int) neighbors->size(); l++){
                if((*neighbors)[l] == edgeVertexId0){
                  lowerId0 = l;
                }
                if((*neighbors)[l] == edgeVertexId1){
                  lowerId1 = l;
                }
              }
              
              if((lowerId0 != -1)&&(lowerId1 != -1)){
                (*seeds)[lowerId0] = 
                  makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
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
  for(int i = 0; i < (int) lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(int i = 0; i < (int) upperList.size(); i++)
    upperList[i] = upperList[i]->find();
  
  vector<UnionFind *>::iterator it;
  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));
  
  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));
  
  if((upperList.size() == 1)&&(lowerList.size() == 1))
    return -2;
  
  return 1;
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::perturbate(
    const dataTypeU &uEpsilon, const dataTypeV &vEpsilon) const{
      
#ifndef TTK_ENABLE_KAMIKAZE
  if(!uField_)
    return -1;
  if(!vField_)
    return -2;
  if(!vertexNumber_)
    return -3;
#endif
  
  dataTypeU *uField = (dataTypeU *) uField_;
  dataTypeV *vField = (dataTypeV *) vField_;
  
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < vertexNumber_; i++){
    // simulation of simplicity in 2 dimensions, need to use degree 2 polynoms
    
    uField[i] += i*uEpsilon;
    vField[i] += (i*vEpsilon)*(i*vEpsilon);
  }
      
  return 0;
}
