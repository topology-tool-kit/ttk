#include <AbstractMorseSmaleComplex.h>

using namespace std;
using namespace ttk;
using namespace dcg;

AbstractMorseSmaleComplex::AbstractMorseSmaleComplex()
  : // all options are disabled by default
    ReverveSaddleMaximumConnection{false}, ReverveSaddleSaddleConnection{false},
    ComputeAscendingSeparatrices1{true}, ComputeDescendingSeparatrices1{true},
    ComputeSaddleConnectors{false}, ComputeAscendingSeparatrices2{false},
    ComputeDescendingSeparatrices2{false}, ReturnSaddleConnectors{false},

    // other class members are value-initialized
    inputScalarField_{}, inputTriangulation_{}, inputOffsets_{},

    outputCriticalPoints_numberOfPoints_{}, outputCriticalPoints_points_{},
    outputCriticalPoints_points_cellDimensions_{},
    outputCriticalPoints_points_cellIds_{},
    outputCriticalPoints_points_cellScalars_{},
    outputCriticalPoints_points_isOnBoundary_{},
    outputCriticalPoints_points_PLVertexIdentifiers_{},
    outputCriticalPoints_points_manifoldSize_{},

    outputSeparatrices1_numberOfPoints_{}, outputSeparatrices1_points_{},
    outputSeparatrices1_points_smoothingMask_{},
    outputSeparatrices1_points_cellDimensions_{},
    outputSeparatrices1_points_cellIds_{}, outputSeparatrices1_numberOfCells_{},
    outputSeparatrices1_cells_{}, outputSeparatrices1_cells_sourceIds_{},
    outputSeparatrices1_cells_destinationIds_{},
    outputSeparatrices1_cells_separatrixIds_{},
    outputSeparatrices1_cells_separatrixTypes_{},
    outputSeparatrices1_cells_separatrixFunctionMaxima_{},
    outputSeparatrices1_cells_separatrixFunctionMinima_{},
    outputSeparatrices1_cells_separatrixFunctionDiffs_{},
    outputSeparatrices1_cells_isOnBoundary_{},

    outputSeparatrices2_numberOfPoints_{}, outputSeparatrices2_points_{},
    outputSeparatrices2_numberOfCells_{}, outputSeparatrices2_cells_{},
    outputSeparatrices2_cells_sourceIds_{},
    outputSeparatrices2_cells_separatrixIds_{},
    outputSeparatrices2_cells_separatrixTypes_{},
    outputSeparatrices2_cells_separatrixFunctionMaxima_{},
    outputSeparatrices2_cells_separatrixFunctionMinima_{},
    outputSeparatrices2_cells_separatrixFunctionDiffs_{},
    outputSeparatrices2_cells_isOnBoundary_{},

    outputAscendingManifold_{}, outputDescendingManifold_{},
    outputMorseSmaleManifold_{} {
}

AbstractMorseSmaleComplex::~AbstractMorseSmaleComplex() {
}

int AbstractMorseSmaleComplex::getDescendingSeparatrices1(
  const vector<Cell> &criticalPoints,
  vector<Separatrix> &separatrices,
  vector<vector<Cell>> &separatricesGeometry) const {

  vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 1)
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

    // add descending vpaths
    {
      const Cell &saddle1 = saddle;

      for(int j = 0; j < 2; ++j) {
        const int shift = j + 2;

        SimplexId vertexId;
        inputTriangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

        vector<Cell> vpath;
        vpath.push_back(saddle1);
        discreteGradient_.getDescendingPath(Cell(0, vertexId), vpath);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 0 and discreteGradient_.isCellCritical(lastCell)) {
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

int AbstractMorseSmaleComplex::setAscendingSegmentation(
  const vector<Cell> &criticalPoints,
  vector<SimplexId> &maxSeeds,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMaxima) const {
  const SimplexId numberOfVertices = inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  const SimplexId numberOfCells = inputTriangulation_->getNumberOfCells();
  vector<SimplexId> morseSmaleManifoldOnCells(numberOfCells, -1);
  const int cellDim = inputTriangulation_->getDimensionality();

  // get the seeds : maxima
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == cellDim)
      maxSeeds.push_back(criticalPoint.id_);
  }
  const SimplexId numberOfSeeds = maxSeeds.size();
  numberOfMaxima = numberOfSeeds;

  // Triangulation method pointers for 3D
  auto getCellFace = &Triangulation::getCellTriangle;
  auto getFaceStarNumber = &Triangulation::getTriangleStarNumber;
  auto getFaceStar = &Triangulation::getTriangleStar;
  if(cellDim == 2) {
    // Triangulation method pointers for 2D
    getCellFace = &Triangulation::getCellEdge;
    getFaceStarNumber = &Triangulation::getEdgeStarNumber;
    getFaceStar = &Triangulation::getEdgeStar;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(SimplexId i = 0; i < numberOfSeeds; ++i) {
    queue<SimplexId> bfs;

    // push the seed
    {
      const SimplexId seedId = maxSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()) {
      const SimplexId cofacetId = bfs.front();
      bfs.pop();

      if(morseSmaleManifoldOnCells[cofacetId] == -1) {
        morseSmaleManifoldOnCells[cofacetId] = i;

        for(int j = 0; j < (cellDim + 1); ++j) {
          SimplexId facetId = -1;
          (inputTriangulation_->*getCellFace)(cofacetId, j, facetId);

          SimplexId starNumber
            = (inputTriangulation_->*getFaceStarNumber)(facetId);
          for(SimplexId k = 0; k < starNumber; ++k) {
            SimplexId neighborId = -1;
            (inputTriangulation_->*getFaceStar)(facetId, k, neighborId);

            if(neighborId == cofacetId) {
              continue;
            }

            const SimplexId pairedCellId = discreteGradient_.getPairedCell(
              Cell(cellDim, neighborId), true);

            if(pairedCellId == facetId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  // put segmentation infos from cells to points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfVertices; ++i) {
    SimplexId starId;
    inputTriangulation_->getVertexStar(i, 0, starId);
    morseSmaleManifold[i] = morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

int AbstractMorseSmaleComplex::setDescendingSegmentation(
  const vector<Cell> &criticalPoints,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMinima) const {
  const SimplexId numberOfVertices = inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  // get the seeds : minima
  vector<SimplexId> seeds;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 0)
      seeds.push_back(criticalPoint.id_);
  }
  const SimplexId numberOfSeeds = seeds.size();
  numberOfMinima = numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(SimplexId i = 0; i < numberOfSeeds; ++i) {
    queue<SimplexId> bfs;

    // push the seed
    {
      const SimplexId seedId = seeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()) {
      const SimplexId vertexId = bfs.front();
      bfs.pop();

      if(morseSmaleManifold[vertexId] == -1) {
        morseSmaleManifold[vertexId] = i;

        const SimplexId edgeNumber
          = inputTriangulation_->getVertexEdgeNumber(vertexId);
        for(SimplexId j = 0; j < edgeNumber; ++j) {
          SimplexId edgeId;
          inputTriangulation_->getVertexEdge(vertexId, j, edgeId);

          for(int k = 0; k < 2; ++k) {
            SimplexId neighborId;
            inputTriangulation_->getEdgeVertex(edgeId, k, neighborId);

            if(neighborId == vertexId) {
              continue;
            }

            const SimplexId pairedCellId
              = discreteGradient_.getPairedCell(Cell(0, neighborId));

            if(pairedCellId == edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

int AbstractMorseSmaleComplex::setFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold) const {

  const size_t nVerts = inputTriangulation_->getNumberOfVertices();

  // associate a unique "sparse region id" to each (ascending, descending) pair

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    const auto d = ascendingManifold[i];
    const auto a = descendingManifold[i];
    if(a == -1 || d == -1) {
      morseSmaleManifold[i] = -1;
    } else {
      morseSmaleManifold[i] = a * numberOfMaxima + d;
    }
  }

  // store the "sparse region ids" by copying the morseSmaleManifold output
  std::vector<SimplexId> sparseRegionIds(
    morseSmaleManifold, morseSmaleManifold + nVerts);

  // get unique "sparse region ids"
  PSORT(sparseRegionIds.begin(), sparseRegionIds.end());
  const auto last = std::unique(sparseRegionIds.begin(), sparseRegionIds.end());
  sparseRegionIds.erase(last, sparseRegionIds.end());

  // "sparse region id" -> "dense region id"
  std::map<SimplexId, size_t> sparseToDenseRegionId{};

  for(size_t i = 0; i < sparseRegionIds.size(); ++i) {
    sparseToDenseRegionId[sparseRegionIds[i]] = i;
  }

  // update region id on all vertices: "sparse id" -> "dense id"

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = sparseToDenseRegionId[morseSmaleManifold[i]];
  }

  return 0;
}
