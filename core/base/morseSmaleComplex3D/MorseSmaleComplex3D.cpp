#include <MorseSmaleComplex3D.h>

using namespace std;
using namespace ttk;
using namespace dcg;

MorseSmaleComplex3D::MorseSmaleComplex3D() : AbstractMorseSmaleComplex() {
}

MorseSmaleComplex3D::~MorseSmaleComplex3D() {
}

int MorseSmaleComplex3D::getAscendingSeparatrices1(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry) const {

  vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 2)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori :
  // numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const SimplexId numberOfSeparatrices = 4 * numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle = criticalPoints[saddleIndex];

    // add ascending vpaths
    {
      const Cell &saddle2 = saddle;

      const SimplexId starNumber
        = inputTriangulation_->getTriangleStarNumber(saddle2.id_);
      for(SimplexId j = 0; j < starNumber; ++j) {
        const SimplexId shift = j;

        SimplexId tetraId;
        inputTriangulation_->getTriangleStar(saddle2.id_, j, tetraId);

        vector<Cell> vpath;
        vpath.push_back(saddle2);
        discreteGradient_.getAscendingPath(Cell(3, tetraId), vpath);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 3 and discreteGradient_.isCellCritical(lastCell)) {
          const SimplexId separatrixIndex = 4 * i + shift;

          separatricesGeometry[separatrixIndex] = std::move(vpath);
          separatrices[separatrixIndex]
            = Separatrix(true, saddle, lastCell, false, separatrixIndex);
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex3D::getSaddleConnectors(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry) const {

  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
  wallId_t descendingWallId = 1;
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 2) {
      const Cell &saddle2 = criticalPoint;

      set<SimplexId> saddles1;
      const wallId_t savedDescendingWallId = descendingWallId;
      discreteGradient_.getDescendingWall(
        descendingWallId, saddle2, isVisited, nullptr, &saddles1);
      ++descendingWallId;

      for(const SimplexId saddle1Id : saddles1) {
        const Cell &saddle1 = Cell(1, saddle1Id);

        vector<Cell> vpath;
        const bool isMultiConnected
          = discreteGradient_.getAscendingPathThroughWall(
            savedDescendingWallId, saddle1, saddle2, isVisited, &vpath);

        const Cell &lastCell = vpath.back();
        if(!isMultiConnected and lastCell.dim_ == saddle2.dim_
           and lastCell.id_ == saddle2.id_) {
          const SimplexId separatrixIndex = separatrices.size();
          separatricesGeometry.push_back(std::move(vpath));
          separatrices.push_back(
            Separatrix(true, saddle1, saddle2, false, separatrixIndex));
        }
      }
    }
  }

  return 0;
}

int MorseSmaleComplex3D::getAscendingSeparatrices2(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry,
  vector<set<SimplexId>> &separatricesSaddles) const {
  const Cell emptyCell;

  vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 1)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  vector<wallId_t> isVisited(numberOfEdges, 0);

  // apriori: by default construction, the separatrices are not valid
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle1 = criticalPoints[saddleIndex];

    vector<Cell> wall;
    discreteGradient_.getAscendingWall(
      saddle1.id_, saddle1, isVisited, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle1, emptyCell, false, i);
  }

  return 0;
}

int MorseSmaleComplex3D::getDescendingSeparatrices2(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry,
  vector<set<SimplexId>> &separatricesSaddles) const {
  const Cell emptyCell;

  vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 2)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  // apriori: by default construction, the separatrices are not valid
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle2 = criticalPoints[saddleIndex];

    vector<Cell> wall;
    discreteGradient_.getDescendingWall(
      saddle2.id_, saddle2, isVisited, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle2, emptyCell, false, i);
  }

  return 0;
}

int MorseSmaleComplex3D::getDualPolygon(const SimplexId edgeId,
                                        vector<SimplexId> &polygon) const {
  // primal: star of edgeId -> dual: vertices of polygon

  // get the vertices of the polygon
  const SimplexId starNumber = inputTriangulation_->getEdgeStarNumber(edgeId);
  polygon.resize(starNumber);
  for(SimplexId i = 0; i < starNumber; ++i) {
    SimplexId starId;
    inputTriangulation_->getEdgeStar(edgeId, i, starId);

    polygon[i] = starId;
  }

  return 0;
}

int MorseSmaleComplex3D::sortDualPolygonVertices(
  vector<SimplexId> &polygon) const {
  // sort the vertices of the polygon to be clockwise
  const SimplexId vertexNumber = polygon.size();
  for(SimplexId i = 1; i < vertexNumber; ++i) {
    const SimplexId previousId = polygon[i - 1];

    // find neighbor of previous one
    bool isFound = false;
    SimplexId index = -1;
    for(SimplexId j = i; j < vertexNumber; j++) {
      const SimplexId currentId = polygon[j];

      // check if current is the neighbor
      const SimplexId neighborNumber
        = inputTriangulation_->getCellNeighborNumber(previousId);
      for(SimplexId k = 0; k < neighborNumber; ++k) {
        SimplexId neighborId;
        inputTriangulation_->getCellNeighbor(previousId, k, neighborId);

        if(neighborId == currentId) {
          isFound = true;
          index = j;
          break;
        }
      }
      if(isFound)
        break;
    }

    // swap foundId and currentId
    if(isFound) {
      const SimplexId tmpId = polygon[index];
      polygon[index] = polygon[i];
      polygon[i] = tmpId;
    }
  }

  return 0;
}
