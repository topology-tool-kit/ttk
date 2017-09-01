#include<MorseSmaleComplex3D.h>

MorseSmaleComplex3D::MorseSmaleComplex3D():
  AbstractMorseSmaleComplex()
{}

MorseSmaleComplex3D::~MorseSmaleComplex3D(){
}

int MorseSmaleComplex3D::getSeparatrices1(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==1 or criticalPoint.dim_==2)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const int numberOfSeparatrices=4*numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle=criticalPoints[saddleIndex];

    // add ascending vpaths
    if(ComputeAscendingSeparatrices1 and saddle.dim_==2){
      const Cell& saddle2=saddle;

      const int starNumber=inputTriangulation_->getTriangleStarNumber(saddle2.id_);
      for(int j=0; j<starNumber; ++j){
        const int shift=j;

        int tetraId;
        inputTriangulation_->getTriangleStar(saddle2.id_, j, tetraId);

        vector<Cell> vpath{saddle2};
        discreteGradient_.getAscendingPath(Cell(3,tetraId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==3 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }

    // add descending vpaths
    if(ComputeDescendingSeparatrices1 and saddle.dim_==1){
      const Cell& saddle1=saddle;

      for(int j=0; j<2; ++j){
        const int shift=j+2;

        int vertexId;
        inputTriangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

        vector<Cell> vpath{saddle1};
        discreteGradient_.getDescendingPath(Cell(0,vertexId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==0 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex3D::getSaddleConnectors(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  int descendingWallId=1;
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2){
      const Cell& saddle2=criticalPoint;

      set<int> saddles1;
      const int savedDescendingWallId=descendingWallId;
      discreteGradient_.getDescendingWall(descendingWallId, saddle2, isVisited, nullptr, &saddles1);
      ++descendingWallId;

      for(const int saddle1Id : saddles1){
        const Cell& saddle1=Cell(1,saddle1Id);

        vector<Cell> vpath;
        const bool isMultiConnected=discreteGradient_.getAscendingPathThroughWall(savedDescendingWallId, saddle1, saddle2, isVisited, &vpath);

        const Cell& lastCell=vpath.back();
        if(!isMultiConnected and lastCell.dim_==saddle2.dim_ and lastCell.id_==saddle2.id_){
          const int separatrixIndex=separatrices.size();
          separatricesGeometry.push_back(std::move(vpath));
          separatrices.push_back(std::move(Separatrix(true,saddle1,saddle2,false,separatrixIndex)));
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex3D::getAscendingSeparatrices2(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry,
    vector<set<int>>& separatricesSaddles) const{
  const Cell emptyCell;

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==1)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls = numberOfSaddles
  const int numberOfSeparatrices=numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
  vector<wallId_t> isVisited(numberOfEdges, 0);

  // apriori: by default construction, the separatrices are not valid
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle1=criticalPoints[saddleIndex];

    vector<Cell> wall;
    discreteGradient_.getAscendingWall(saddle1.id_, saddle1, isVisited, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i]=std::move(wall);
    separatrices[i]=std::move(Separatrix(true,saddle1,emptyCell,false,i));
  }

  return 0;
}

int MorseSmaleComplex3D::getDescendingSeparatrices2(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry,
    vector<set<int>>& separatricesSaddles) const{
  const Cell emptyCell;

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls = numberOfSaddles
  const int numberOfSeparatrices=numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  // apriori: by default construction, the separatrices are not valid
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle2=criticalPoints[saddleIndex];

    vector<Cell> wall;
    discreteGradient_.getDescendingWall(saddle2.id_, saddle2, isVisited, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i]=std::move(wall);
    separatrices[i]=std::move(Separatrix(true,saddle2,emptyCell,false,i));
  }

  return 0;
}

int MorseSmaleComplex3D::getDualPolygon(const int edgeId, vector<int>& polygon) const{
  // primal: star of edgeId -> dual: vertices of polygon

  // get the vertices of the polygon
  const int starNumber=inputTriangulation_->getEdgeStarNumber(edgeId);
  polygon.resize(starNumber);
  for(int i=0; i<starNumber; ++i){
    int starId;
    inputTriangulation_->getEdgeStar(edgeId, i, starId);

    polygon[i]=starId;
  }

  return 0;
}

int MorseSmaleComplex3D::sortDualPolygonVertices(vector<int>& polygon) const{
  // sort the vertices of the polygon to be clockwise
  const int vertexNumber=polygon.size();
  for(int i=1; i<vertexNumber; ++i){
    const int previousId=polygon[i-1];

    // find neighbor of previous one
    bool isFound=false;
    int index=-1;
    for(int j=i; j<vertexNumber; j++){
      const int currentId=polygon[j];

      // check if current is the neighbor
      const int neighborNumber=inputTriangulation_->getCellNeighborNumber(previousId);
      for(int k=0; k<neighborNumber; ++k){
        int neighborId;
        inputTriangulation_->getCellNeighbor(previousId, k, neighborId);

        if(neighborId==currentId){
          isFound=true;
          index=j;
          break;
        }
      }
      if(isFound) break;
    }

    // swap foundId and currentId
    if(isFound){
      const int tmpId=polygon[index];
      polygon[index]=polygon[i];
      polygon[i]=tmpId;
    }
  }

  return 0;
}

int MorseSmaleComplex3D::setDescendingSegmentation(const vector<Cell>& criticalPoints,
    int* const morseSmaleManifold,
    int& numberOfMinima){
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  // get the seeds : minima
  vector<int> minSeeds;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==0)
      minSeeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=minSeeds.size();
  numberOfMinima=numberOfSeeds;

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=minSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int vertexId=bfs.front();
      bfs.pop();

      if(morseSmaleManifold[vertexId]==-1){
        morseSmaleManifold[vertexId]=i;

        const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(vertexId);
        for(int j=0; j<edgeNumber; ++j){
          int edgeId;
          inputTriangulation_->getVertexEdge(vertexId, j, edgeId);

          for(int k=0; k<2; ++k){
            int neighborId;
            inputTriangulation_->getEdgeVertex(edgeId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(0, neighborId));

            if(pairedCellId==edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex3D::setAscendingSegmentation(const vector<Cell>& criticalPoints,
    vector<int>& maxSeeds,
    int* const morseSmaleManifold,
    int& numberOfMaxima){
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  vector<int> morseSmaleManifoldOnCells(numberOfCells, -1);

  // get the seeds : maxima
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==3)
      maxSeeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=maxSeeds.size();
  numberOfMaxima=numberOfSeeds;

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=maxSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int tetraId=bfs.front();
      bfs.pop();

      if(morseSmaleManifoldOnCells[tetraId]==-1){
        morseSmaleManifoldOnCells[tetraId]=i;

        for(int j=0; j<4; ++j){
          int triangleId;
          inputTriangulation_->getCellTriangle(tetraId, j, triangleId);

          const int starNumber=inputTriangulation_->getTriangleStarNumber(triangleId);
          for(int k=0; k<starNumber; ++k){
            int neighborId;
            inputTriangulation_->getTriangleStar(triangleId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(3, neighborId), true);

            if(pairedCellId==triangleId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  // put segmentation infos from cells to points
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfVertices; ++i){
    int starId;
    inputTriangulation_->getVertexStar(i, 0, starId);
    morseSmaleManifold[i]=morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

int MorseSmaleComplex3D::setFinalSegmentation(const int numberOfMaxima,
    const int numberOfMinima,
    const int* const ascendingManifold,
    const int* const descendingManifold,
    int* const morseSmaleManifold) const{
  vector<vector<pair<int,int>>> minTable(numberOfMinima);

  int id{};
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  for(int i=0; i<numberOfVertices; ++i){
    const int d=ascendingManifold[i];
    const int a=descendingManifold[i];

    if(a==-1 or d==-1){
      morseSmaleManifold[i]=-1;
      continue;
    }

    vector<pair<int,int>>& table=minTable[a];
    int foundId=-1;
    for(const pair<int,int>& p : table){
      if(p.first == d)
        foundId=p.second;
    }

    // add new association (a,d)
    if(foundId==-1){
      table.push_back(make_pair(d,id));
      morseSmaleManifold[i]=id;
      ++id;
    }
    // update to saved associationId
    else
      morseSmaleManifold[i]=foundId;
  }

  return 0;
}
