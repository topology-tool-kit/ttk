#include<MorseSmaleComplex3D.h>

MorseSmaleComplex3D::MorseSmaleComplex3D():
  AbstractMorseSmaleComplex()
{}

MorseSmaleComplex3D::~MorseSmaleComplex3D(){
}

int MorseSmaleComplex3D::getAscendingSeparatrices1(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const int numberOfSeparatrices=4*numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle=criticalPoints[saddleIndex];

    // add ascending vpaths
    {
      const Cell& saddle2=saddle;

      const int starNumber=inputTriangulation_->getTriangleStarNumber(saddle2.id_);
      for(int j=0; j<starNumber; ++j){
        const int shift=j;

        int tetraId;
        inputTriangulation_->getTriangleStar(saddle2.id_, j, tetraId);

        vector<Cell> vpath;
        vpath.push_back(saddle2);
        discreteGradient_.getAscendingPath(Cell(3,tetraId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==3 and discreteGradient_.isCellCritical(lastCell)){
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
#ifdef TTK_ENABLE_OPENMP
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
#ifdef TTK_ENABLE_OPENMP
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
