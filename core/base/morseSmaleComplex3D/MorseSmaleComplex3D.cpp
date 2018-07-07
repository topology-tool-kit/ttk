#include<MorseSmaleComplex3D.h>

using namespace std;
using namespace ttk;
using namespace dcg;

MorseSmaleComplex3D::MorseSmaleComplex3D():
  AbstractMorseSmaleComplex()
{}

MorseSmaleComplex3D::~MorseSmaleComplex3D(){
}

int MorseSmaleComplex3D::getAscendingSeparatrices1(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  vector<simplexId_t> saddleIndexes;
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2)
      saddleIndexes.push_back(i);
  }
  const simplexId_t numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const simplexId_t numberOfSeparatrices=4*numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(simplexId_t i=0; i<numberOfSaddles; ++i){
    const simplexId_t saddleIndex=saddleIndexes[i];
    const Cell& saddle=criticalPoints[saddleIndex];

    // add ascending vpaths
    {
      const Cell& saddle2=saddle;

      const simplexId_t starNumber=inputTriangulation_->getTriangleStarNumber(saddle2.id_);
      for(simplexId_t j=0; j<starNumber; ++j){
        const simplexId_t shift=j;

        simplexId_t tetraId;
        inputTriangulation_->getTriangleStar(saddle2.id_, j, tetraId);

        vector<Cell> vpath;
        vpath.push_back(saddle2);
        discreteGradient_.getAscendingPath(Cell(3,tetraId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==3 and discreteGradient_.isCellCritical(lastCell)){
          const simplexId_t separatrixIndex=4*i+shift;

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

  const simplexId_t numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  wallId_t descendingWallId=1;
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2){
      const Cell& saddle2=criticalPoint;

      set<simplexId_t> saddles1;
      const wallId_t savedDescendingWallId=descendingWallId;
      discreteGradient_.getDescendingWall(descendingWallId, saddle2, isVisited, nullptr, &saddles1);
      ++descendingWallId;

      for(const simplexId_t saddle1Id : saddles1){
        const Cell& saddle1=Cell(1,saddle1Id);

        vector<Cell> vpath;
        const bool isMultiConnected=discreteGradient_.getAscendingPathThroughWall(savedDescendingWallId, saddle1, saddle2, isVisited, &vpath);

        const Cell& lastCell=vpath.back();
        if(!isMultiConnected and lastCell.dim_==saddle2.dim_ and lastCell.id_==saddle2.id_){
          const simplexId_t separatrixIndex=separatrices.size();
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
    vector<set<simplexId_t>>& separatricesSaddles) const{
  const Cell emptyCell;

  vector<simplexId_t> saddleIndexes;
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==1)
      saddleIndexes.push_back(i);
  }
  const simplexId_t numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls = numberOfSaddles
  const simplexId_t numberOfSeparatrices=numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const simplexId_t numberOfEdges=inputTriangulation_->getNumberOfEdges();
  vector<wallId_t> isVisited(numberOfEdges, 0);

  // apriori: by default construction, the separatrices are not valid
  for(simplexId_t i=0; i<numberOfSaddles; ++i){
    const simplexId_t saddleIndex=saddleIndexes[i];
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
    vector<set<simplexId_t>>& separatricesSaddles) const{
  const Cell emptyCell;

  vector<simplexId_t> saddleIndexes;
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==2)
      saddleIndexes.push_back(i);
  }
  const simplexId_t numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls = numberOfSaddles
  const simplexId_t numberOfSeparatrices=numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const simplexId_t numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  // apriori: by default construction, the separatrices are not valid
  for(simplexId_t i=0; i<numberOfSaddles; ++i){
    const simplexId_t saddleIndex=saddleIndexes[i];
    const Cell& saddle2=criticalPoints[saddleIndex];

    vector<Cell> wall;
    discreteGradient_.getDescendingWall(saddle2.id_, saddle2, isVisited, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i]=std::move(wall);
    separatrices[i]=std::move(Separatrix(true,saddle2,emptyCell,false,i));
  }

  return 0;
}

int MorseSmaleComplex3D::getDualPolygon(const simplexId_t edgeId, vector<simplexId_t>& polygon) const{
  // primal: star of edgeId -> dual: vertices of polygon

  // get the vertices of the polygon
  const simplexId_t starNumber=inputTriangulation_->getEdgeStarNumber(edgeId);
  polygon.resize(starNumber);
  for(simplexId_t i=0; i<starNumber; ++i){
    simplexId_t starId;
    inputTriangulation_->getEdgeStar(edgeId, i, starId);

    polygon[i]=starId;
  }

  return 0;
}

int MorseSmaleComplex3D::sortDualPolygonVertices(vector<simplexId_t>& polygon) const{
  // sort the vertices of the polygon to be clockwise
  const simplexId_t vertexNumber=polygon.size();
  for(simplexId_t i=1; i<vertexNumber; ++i){
    const simplexId_t previousId=polygon[i-1];

    // find neighbor of previous one
    bool isFound=false;
    simplexId_t index=-1;
    for(simplexId_t j=i; j<vertexNumber; j++){
      const simplexId_t currentId=polygon[j];

      // check if current is the neighbor
      const simplexId_t neighborNumber=inputTriangulation_->getCellNeighborNumber(previousId);
      for(simplexId_t k=0; k<neighborNumber; ++k){
        simplexId_t neighborId;
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
      const simplexId_t tmpId=polygon[index];
      polygon[index]=polygon[i];
      polygon[i]=tmpId;
    }
  }

  return 0;
}
