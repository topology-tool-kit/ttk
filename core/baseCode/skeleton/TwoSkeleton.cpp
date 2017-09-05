#include                  <TwoSkeleton.h>

TwoSkeleton::TwoSkeleton(){

}

TwoSkeleton::~TwoSkeleton(){
  
}

int TwoSkeleton::buildCellNeighborsFromVertices(const int &vertexNumber, 
  const int &cellNumber, const long long int *cellArray, 
  vector<vector<int> > &cellNeighbors, 
  vector<vector<int> > *vertexStars) const{

  Timer t;
  
  bool localVertexStarsAlloc = false;
  vector<vector<int> > *localVertexStars = vertexStars;
  
  if(!localVertexStars){
    localVertexStars = new vector<vector<int> >();
    localVertexStarsAlloc = true;
  }
  
  if(!localVertexStars->size()){
    
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.buildVertexStars(vertexNumber, cellNumber, cellArray, 
      *localVertexStars);
  }
  
  int vertexPerCell = cellArray[0];
  
  cellNeighbors.resize(cellNumber);
  for(int i = 0; i < (int) cellNeighbors.size(); i++)
    cellNeighbors[i].reserve(vertexPerCell);
  
  // pre-sort vertex stars
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < vertexNumber; i++)
    sort((*localVertexStars)[i].begin(), (*localVertexStars)[i].end());
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < cellNumber; i++){
    
    for(int j = 0; j < vertexPerCell; j++){
      
      int v0 = cellArray[(vertexPerCell + 1)*i + 1 + j];
      int v1 = cellArray[(vertexPerCell + 1)*i + 1 + (j + 1)%vertexPerCell];
      
      // perform an intersection of the 2 sorted star lists
      int pos0 = 0, pos1 = 0;
      int intersection = -1;
      
      while((pos0 < (int) (*localVertexStars)[v0].size())
        &&(pos1 < (int) (*localVertexStars)[v1].size())){
        
        int biggest = (*localVertexStars)[v0][pos0];
        if((*localVertexStars)[v1][pos1] > biggest){
          biggest = (*localVertexStars)[v1][pos1];
        }
        
        for(int l = pos0; l < (int) (*localVertexStars)[v0].size(); l++){
          if((*localVertexStars)[v0][l] < biggest){
            pos0++;
          }
          else{
            break;
          }
        }
        for(int l = pos1; l < (int) (*localVertexStars)[v1].size(); l++){
          if((*localVertexStars)[v1][l] < biggest){
            pos1++;
          }
          else{
            break;
          }
        }
        
        if((*localVertexStars)[v0][pos0]
          == (*localVertexStars)[v1][pos1]){
          
          if((*localVertexStars)[v0][pos0] != i){
            intersection = (*localVertexStars)[v0][pos0];
            break;
          }
        
          pos0++;
          pos1++;
        }
      }
      
      if(intersection != -1){
        cellNeighbors[i].push_back(intersection);
      }
    }
  }
  
  if(localVertexStarsAlloc)
    delete localVertexStars;
  
  {
    stringstream msg;
    msg << "[TwoSkeleton] Cell neighbors (" << cellNumber
      << " cells) computed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

int TwoSkeleton::buildEdgeTriangles(const int &vertexNumber, 
  const int &cellNumber, const long long int *cellArray, 
  vector<vector<int> > &edgeTriangleList, 
  vector<vector<int> > *vertexStarList,
  vector<pair<int, int> > *edgeList, 
  vector<vector<int> > *edgeStarList, 
  vector<vector<int> > *triangleList, 
  vector<vector<int> > *triangleStarList, 
  vector<vector<int> > *cellTriangleList) const{

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef withKamikaze
  if(vertexNumber <= 0)
    return -1;
  if(cellNumber <= 0)
    return -2;
  if(!cellArray)
    return -3;
#endif
  
  // NOTE:
  // vertexStarList: this guy we can compute if provided but we don't actually
  // need it.
  
  bool localEdgeListAlloc = false;
  vector<pair<int, int> > *localEdgeList = edgeList;
  
  if(!localEdgeList){
    localEdgeList = new vector<pair<int, int> >();
    localEdgeListAlloc = true;
  }
  
  bool localEdgeStarListAlloc = false;
  vector<vector<int> > *localEdgeStarList = edgeStarList;
  
  if(!localEdgeStarList){
    localEdgeStarList = new vector<vector<int> >();
    localEdgeStarListAlloc = true;
  }
  
  bool localTriangleListAlloc = false;
  vector<vector<int> > *localTriangleList = triangleList;
  
  if(!localTriangleList){
    localTriangleList = new vector<vector<int> >();
    localTriangleListAlloc = true;
  }
  
  // NOTE:
  // triangleStarList: this guy we can compute if provided but we don't actually
  // need it.
  
  bool localCellTriangleListAlloc = false;
  vector<vector<int> > *localCellTriangleList = cellTriangleList;
  
  if(!localCellTriangleList){
    localCellTriangleList = new vector<vector<int> >();
    localCellTriangleListAlloc = true;
  }
  
  OneSkeleton oneSkeleton;
  oneSkeleton.setDebugLevel(debugLevel_);
  oneSkeleton.setThreadNumber(threadNumber_);
  
  // now do the pre-computation
  if(localEdgeList->empty()){
    oneSkeleton.buildEdgeList(vertexNumber, cellNumber, cellArray, 
      (*localEdgeList));
  }
  
  if(localEdgeStarList->empty()){
    oneSkeleton.buildEdgeStars(vertexNumber, cellNumber, cellArray,
      (*localEdgeStarList), 
      localEdgeList,
      vertexStarList);
  }
  
  if((localTriangleList->empty())||(localCellTriangleList->empty())){
    buildTriangleList(vertexNumber, cellNumber, cellArray,
      localTriangleList, triangleStarList, localCellTriangleList);
  }
  
  edgeTriangleList.resize(localEdgeList->size());
  
  // alright, let's get things done now.
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) localEdgeList->size(); i++){
    
    int vertexId0, vertexId1, vertexId2;
    
    for(int j = 0; j < (int) (*localEdgeStarList)[i].size(); j++){
      int tetId = (*localEdgeStarList)[i][j];
      
      for(int k = 0; k < (int) (*localCellTriangleList)[tetId].size(); k++){
        int triangleId = (*localCellTriangleList)[tetId][k];
        
        bool isAttached = false;
        
        vertexId0 = (*localTriangleList)[triangleId][0];
        vertexId1 = (*localTriangleList)[triangleId][1];
        vertexId2 = (*localTriangleList)[triangleId][2];
        
        if((*localEdgeList)[i].first == vertexId0){
          if(((*localEdgeList)[i].second == vertexId1)
            ||((*localEdgeList)[i].second == vertexId2)){
            isAttached = true;
          }
        }
        
        if((*localEdgeList)[i].first == vertexId1){
          if(((*localEdgeList)[i].second == vertexId0)
            ||((*localEdgeList)[i].second == vertexId2)){
            isAttached = true;
          }
        }
        
        if((*localEdgeList)[i].first == vertexId2){
          if(((*localEdgeList)[i].second == vertexId1)
            ||((*localEdgeList)[i].second == vertexId0)){
            isAttached = true;
          }
        }
        
        if(isAttached){
        
          bool isIn = false;
          for(int l = 0; l < (int) edgeTriangleList[i].size(); l++){
            if(edgeTriangleList[i][l] == triangleId){
              isIn = true;
              break;
            }
          }
          if(!isIn){
            edgeTriangleList[i].push_back(triangleId);
          }
        }
      }
    }
  }
  
  int edgeNumber = localEdgeList->size();
  int triangleNumber = localTriangleList->size();
  
  if(localEdgeListAlloc){
    delete localEdgeList;
  }
  
  if(localEdgeStarListAlloc){
    delete localEdgeStarList;
  }
  
  if(localTriangleListAlloc){
    delete localTriangleList;
  }
  
  if(localCellTriangleListAlloc){
    delete localCellTriangleList;
  }
    
  {
    stringstream msg;
    msg << "[TwoSkeleton] Edge triangles (" << edgeNumber
      << " edge(s), " << triangleNumber 
      << " triangle(s)) computed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}


int TwoSkeleton::buildTriangleList(const int &vertexNumber, 
  const int &cellNumber, const long long int *cellArray,
  vector<vector<int> > *triangleList,
  vector<vector<int> > *triangleStars,
  vector<vector<int> > *cellTriangleList) const{

  Timer t;
  
  // NOTE: parallel implementation not efficient (see bench at the bottom)
  // let's force the usage of only 1 thread.
  int oldThreadNumber = threadNumber_;
  threadNumber_ = 1;
 
  int triangleNumber = 0;
  
  // check the consistency of the variables -- to adapt
#ifndef withKamikaze
  if(vertexNumber <= 0)
    return -1;
  if(cellNumber <= 0)
    return -2;
  if(!cellArray)
    return -3;
  if((!triangleList)&&(!triangleStars)&&(!cellTriangleList)){
    // we've got nothing to do here.
    return -4;
  }
#endif

  if(triangleList){
    // NOTE: 9 is pretty empirical here...
    triangleList->clear();
    triangleList->reserve(9*vertexNumber);
  }
  if(triangleStars){
    triangleStars->clear();
    triangleStars->reserve(9*vertexNumber);
  }
  if(cellTriangleList){
    cellTriangleList->resize(cellNumber);
    for(int i = 0; i < cellNumber; i++)
      // assuming tet-mesh here
      (*cellTriangleList)[i].resize(4, -1);
  }

  if(threadNumber_ == 1){
    
    // for each vertex, 
    //   list of triangles
    //    each triangle is a list vertex Id + a triangleId
    vector<vector<pair<vector<int>, int > > > triangleTable(vertexNumber);
    for(int i = 0; i < vertexNumber; i++)
      triangleTable[i].reserve(32);
    
    for(int i = 0; i < cellNumber; i++){
      
      if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){
        
        vector<int> triangle(3);
        
        for(int j = 0; j < 4; j++){
          // doing triangle j
          
          for(int k = 0; k < 3; k++){
            triangle[k] = cellArray[5*i + 1 + (j + k)%4];
          }
          sort(triangle.begin(), triangle.end());
          
          int triangleId = -1;
          for(int k = 0; 
            k < (int) triangleTable[triangle[0]].size(); k++){
              
            // processing a triangle stored for that vertex
            if((triangleTable[triangle[0]][k].first[1] == triangle[1])
              &&
              (triangleTable[triangle[0]][k].first[2] == triangle[2])
              ){
              triangleId = triangleTable[triangle[0]][k].second;
              break;
            }
          }
          if(triangleId == -1){
            // not found yet
            triangleId = triangleNumber;
            triangleTable[triangle[0]].push_back(
              pair<vector<int>, int>(triangle, triangleNumber));
            triangleNumber++;
            
            if(triangleList){
              triangleList->size();
              triangleList->push_back(triangle);
            }
            if(triangleStars){
              triangleStars->resize(triangleStars->size() + 1);
              triangleStars->back().resize(1);
              // store the tet i in the triangleStars list
              triangleStars->back().back() = i;
            }
          }
          else{
            if(triangleStars){
              // add tet i as a neighbor of triangleId
              (*triangleStars)[triangleId].push_back(i);
            }
          }
          
          if(cellTriangleList){
            // add the triangle to the cell
            for(int k = 0; k < 4; k++){
              if((*cellTriangleList)[i][k] == -1){
                (*cellTriangleList)[i][k] = triangleId;
                break;
              }
            }
          }
        }
        
        // update the progress bar of the wrapping code -- to adapt
        if(debugLevel_ > advancedInfoMsg){
          if((wrapper_)
            &&(!(i % ((cellNumber)/10)))){
            wrapper_->updateProgress((i + 1.0)
              /cellNumber);
          }
        }
      }
    }
  }
  else{
 
    Timer memTimer;
    //  for each thread, 
    //    for each vertex
    //      for each triangle
    //        vertex ids
    //        triangleId
    vector<vector<vector<pair<vector<int>, int > > > >
      threadedTriangleTable(threadNumber_);
      
    vector<int> triangleNumbers(threadNumber_, 0);
    
    vector<vector<vector<int> > > threadedTriangleStars(threadNumber_);
      
    for(int i = 0; i < threadNumber_; i++){
      threadedTriangleTable[i].resize(vertexNumber);
    }
    
    int count = 0;
    
    // the following open-mp processing is only relevant for embarassingly 
    // parallel algorithms (such as smoothing) -- to adapt
#ifdef withOpenMP
    omp_lock_t writeLock;
    omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
    for(int i = 0; i < (int) cellNumber; i++){
    
      // avoid any processing if the abort signal is sent
      if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

        int threadId = 1;
#ifdef withOpenMP
        threadId = omp_get_thread_num();
#endif     
        
        vector<int> triangle(3);
        
        for(int j = 0; j < 4; j++){
          // doing triangle j
          
          for(int k = 0; k < 3; k++){
            triangle[k] = cellArray[5*i + 1 + (j + k)%4];
          }
          sort(triangle.begin(), triangle.end());
          
          int triangleId = -1;
          for(int k = 0; 
            k < (int) threadedTriangleTable[threadId][triangle[0]].size(); k++){
              
            // processing a triangle stored for that vertex
            if((threadedTriangleTable[threadId][triangle[0]][k].first[1] 
              == triangle[1])
              &&
              (threadedTriangleTable[threadId][triangle[0]][k].first[2] 
              == triangle[2])
              ){
              triangleId = 
                threadedTriangleTable[threadId][triangle[0]][k].second;
              break;
            }
          }
          if(triangleId == -1){
            // not visited by this thread
            threadedTriangleTable[threadId][triangle[0]].push_back(
              pair<vector<int>, int>(triangle, triangleNumbers[threadId]));
            triangleNumbers[threadId]++;
            
            if(triangleStars){
              threadedTriangleStars[threadId].resize(
                threadedTriangleStars[threadId].size() + 1);
              threadedTriangleStars[threadId].back().resize(1);
              threadedTriangleStars[threadId].back().back() = i;
            }
          }
          else{
            if(triangleStars){
              threadedTriangleStars[threadId][triangleId].push_back(i);
            }
          }
        }
        
        // update the progress bar of the wrapping code -- to adapt
        if(debugLevel_ > advancedInfoMsg){
#ifdef withOpenMP
          omp_set_lock(&writeLock);
#endif
          if((wrapper_)
            &&(!(count % ((cellNumber)/10)))){
            wrapper_->updateProgress((count + 1.0)
              /cellNumber);
          }

          count++;
#ifdef withOpenMP
          omp_unset_lock(&writeLock);
#endif
        }
      }
    }
    
#ifdef withOpenMP
    omp_destroy_lock(&writeLock);
#endif
  
    Timer mergeTimer;
    // now merge the lists
    vector<vector<pair<vector<int>, int > > > mainTriangleTable(vertexNumber);
//     for(int i = 0; i < vertexNumber; i++){
//       mainTriangleTable[i].reserve(32);
//     }
    
    for(int i = 0; i < threadNumber_; i++){
      
      // vertex j
      for(int j = 0; j < (int) threadedTriangleTable[i].size(); j++){
        
        // triangle k
        for(int k = 0; k < (int) threadedTriangleTable[i][j].size(); k++){
          
          int triangleId = -1;
          // triangle l 
          for(int l = 0; l < (int) mainTriangleTable[j].size(); l++){
            if((threadedTriangleTable[i][j][k].first[0] == 
              mainTriangleTable[j][l].first[0])
              &&
              (threadedTriangleTable[i][j][k].first[1] == 
                mainTriangleTable[j][l].first[1])
              &&
              (threadedTriangleTable[i][j][k].first[2] == 
                mainTriangleTable[j][l].first[2])){
              
              triangleId = mainTriangleTable[j][l].second;
              break;
            }
          }
          if(triangleId == -1){
            
            triangleId = triangleNumber;
            mainTriangleTable[j].push_back(
              pair<vector<int>, int>(
                threadedTriangleTable[i][j][k].first, triangleNumber));
            triangleNumber++;
            
            if(triangleList)
              triangleList->push_back(mainTriangleTable[j].back().first);
            
            if(triangleStars){
              triangleStars->resize(triangleStars->size() + 1);
              triangleStars->back().resize(1);
              triangleStars->back() = 
                threadedTriangleStars[i][threadedTriangleTable[i][j][k].second];
            }
          }
          else{
            if(triangleStars){
              for(int l = 0; l < (int) threadedTriangleStars[
                i][threadedTriangleTable[i][j][k].second].size();
                l++){
                (*triangleStars)[triangleId].push_back(
                  threadedTriangleStars[
                    i][threadedTriangleTable[i][j][k].second][l]);
              }
            }
          }
          
          if(cellTriangleList){
            // add the triangle to the cell
            int cellId = -1;
            
            for(int l = 0; l < (int) threadedTriangleStars[
              i][threadedTriangleTable[i][j][k].second].size();
              l++){
              cellId = threadedTriangleStars[
                  i][threadedTriangleTable[i][j][k].second][l];
              for(int m = 0; m < 4; m++){
                if((*cellTriangleList)[cellId][m] == -1){
                  (*cellTriangleList)[cellId][m] = triangleId;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  
  {
    stringstream msg;
    msg << "[TwoSkeleton] Triangle list (" << triangleNumber
      << " triangles) computed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
 
  threadNumber_ = oldThreadNumber;
 
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 58.5631 s
  // 24 threads: 87.5816 s (~)
  
  // ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - no tet adj
  // 1 thread: 5.7427 s
  // 4 threads: 9.14764 s (~)
  
  // ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - tet adjacency
  // 1 thread: 8.59854 s
  // 4 threads: 15.837 s (~)
  
  // ethaneDiol.vtu, 8.7Mtets, vger (2coresHT) - tet adjacency only
  // 1 thread: 6.85897 s
  // 4 threads: 14.78 s (~)
   
  return 0;
}

int TwoSkeleton::buildTriangleEdgeList(const int &vertexNumber, 
  const int &cellNumber, const long long int *cellArray, 
  vector<vector<int> > &triangleEdgeList, 
  vector<vector<int> > *vertexEdgeList, 
  vector<pair<int, int> > *edgeList, 
  vector<vector<int> > *triangleList,
  vector<vector<int> > *triangleStarList,
  vector<vector<int> > *cellTriangleList) const{

  Timer t;
  
  bool localEdgeListAlloc = false;
  vector<pair<int, int> > *localEdgeList = edgeList;
  if(!localEdgeList){
    localEdgeList = new vector<pair<int, int> >();
    localEdgeListAlloc = true;
  }
  
  if(!localEdgeList->size()){
    OneSkeleton oneSkeleton;
    oneSkeleton.setDebugLevel(debugLevel_);
    oneSkeleton.setThreadNumber(threadNumber_);
    
    oneSkeleton.buildEdgeList(vertexNumber, cellNumber, cellArray,
      (*localEdgeList));
  }
  
  bool localVertexEdgeListAlloc = false;
  vector<vector<int> > *localVertexEdgeList = vertexEdgeList;
  if(!localVertexEdgeList){
    localVertexEdgeList = new vector<vector<int> >();
    localVertexEdgeListAlloc = true;
  }
  
  if(!localVertexEdgeList->size()){
    
    ZeroSkeleton zeroSkeleton;
    zeroSkeleton.setDebugLevel(debugLevel_);
    zeroSkeleton.setThreadNumber(threadNumber_);
    
    zeroSkeleton.buildVertexEdges(vertexNumber, (*localEdgeList), 
      (*localVertexEdgeList));
  }
  
  // NOTE:
  // triangleStarList and cellTriangleList: we do not need these guys but we 
  // can compute them for free optionally.
  
  bool localTriangleListAlloc = false;
  vector<vector<int> > *localTriangleList = triangleList;
  if(!localTriangleList){
    localTriangleList = new vector<vector<int> >();
    localTriangleListAlloc = true;
  }
  if(!localTriangleList->size()){
    
    buildTriangleList(vertexNumber, cellNumber, cellArray,
      localTriangleList,
      triangleStarList,
      cellTriangleList);
  }
  
  triangleEdgeList.resize(localTriangleList->size());
  
  // now for each triangle, grab its vertices, add the edges in the triangle
  // with no duplicate
  // let's do the real stuff
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) localTriangleList->size(); i++){
    int vertexId = -1;
    for(int j = 0; j < (int) (*localTriangleList)[i].size(); j++){
      vertexId = (*localTriangleList)[i][j];
      
      for(int k = 0; k < (int) (*localVertexEdgeList)[vertexId].size(); k++){
        int edgeId = (*localVertexEdgeList)[vertexId][k];
        
        int otherVertexId = (*localEdgeList)[edgeId].first;
        
        if(otherVertexId == vertexId){
          otherVertexId = (*localEdgeList)[edgeId].second;
        }
        
        bool isInTriangle = false;
        for(int l = 0; l < (int) (*localTriangleList)[i].size(); l++){
          if((*localTriangleList)[i][l] == otherVertexId){
            isInTriangle = true;
            break;
          }
        }
        
        if(isInTriangle){
          bool isIn = false;
          for(int l = 0; l < (int) triangleEdgeList[i].size(); l++){
            if(triangleEdgeList[i][l] == edgeId){
              isIn = true;
              break;
            }
          }
          if(!isIn){
            triangleEdgeList[i].push_back(edgeId);
          }
        }
      }
    }
  }
  
  int triangleNumber = localTriangleList->size();
  int edgeNumber = localEdgeList->size();
  
  if(localEdgeListAlloc){
    delete localEdgeList;
  }
  
  if(localVertexEdgeListAlloc){
    delete localVertexEdgeList;
  }
  
  if(localTriangleListAlloc){
    delete localTriangleList;
  }
  
  {
    stringstream msg;
    msg << "[TwoSkeleton] Triangle edge list (" 
      << triangleNumber << " triangle(s), "
      << edgeNumber << " edge(s)) computed in "
      << t.getElapsedTime() << " s. ("
      << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
    
  return 0;
}

int TwoSkeleton::buildTriangleLinks(const vector<vector<int> > &triangleList,
  const vector<vector<int> > &triangleStars,
  const long long int *cellArray,
  vector<vector<int> > &triangleLinks) const{

#ifndef withKamikaze
  if(triangleList.empty())
    return -1;
  if((triangleStars.empty())||(triangleStars.size() != triangleList.size()))
    return -2;
  if(!cellArray)
    return -3;
#endif
    
  Timer t;
    
  triangleLinks.resize(triangleList.size());
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) triangleList.size(); i++){
    
    for(int j = 0; j < (int) triangleStars[i].size(); j++){
      
      for(int k = 0; k < 4; k++){
        int vertexId = cellArray[5*triangleStars[i][j] + 1 + k];
        
        if((vertexId != triangleList[i][0])
          &&(vertexId != triangleList[i][1])
          &&(vertexId != triangleList[i][2])){
          triangleLinks[i].push_back(vertexId);
          break;
        }
      }
    }
  }
  
  {
    stringstream msg;
    msg << "[TwoSkeleton] Triangle links built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
    
  return 0;
}

int TwoSkeleton::buildVertexTriangles(
  const int &vertexNumber,
  const vector<vector<int> > &triangleList,
  vector<vector<int> > &vertexTriangleList) const{
    
  Timer t;

  int oldThreadNumber = threadNumber_;
  
  threadNumber_ = 1;
  
  if(threadNumber_ != 1){
    vector<vector<vector<int> > > threadedLists(threadNumber_);
    for(int i = 0; i < threadNumber_; i++){
      threadedLists[i].resize(vertexNumber);
    }

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) triangleList.size(); i++){
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif

      for(int j = 0; j < (int) triangleList[i].size(); j++){
        threadedLists[threadId][triangleList[i][j]].push_back(i);
      }
    }

    // merge back
    vertexTriangleList.resize(vertexNumber, vector<int>());
    for(int i = 0; i < threadNumber_; i++){
      for(int j = 0; j < vertexNumber; j++){
        for(int k = 0; k < (int) threadedLists[i][j].size(); k++){
          vertexTriangleList[j].push_back(
            threadedLists[i][j][k]);
        }
      }
    }
  }
  else{
    vertexTriangleList.resize(vertexNumber);
    for(int i = 0; i < (int) triangleList.size(); i++){
      for(int j = 0; j < (int) triangleList[i].size(); j++){
        vertexTriangleList[triangleList[i][j]].push_back(i);
      }
    }
  }

    
  {
    stringstream msg;
    msg << "[TwoSkeleton] Vertex triangle list (" << vertexNumber
      << " vertices) computed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  threadNumber_ = oldThreadNumber;
  
  return 0;
}
