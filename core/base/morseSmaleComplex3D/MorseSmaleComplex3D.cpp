#include <MorseSmaleComplex3D.h>
#include <iterator>

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

  const auto nTriangles = inputTriangulation_->getNumberOfTriangles();
  // visited triangles (one vector per thread)
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  std::vector<std::vector<SimplexId>> visitedTriangles(this->threadNumber_);

  for(auto &vec : isVisited) {
    // resize threads outside of main loop
    vec.resize(nTriangles, false);
  }

  // list of 2-saddles
  std::vector<Cell> saddles2{};
  // copy cells instead of taking a reference?
  std::copy_if(criticalPoints.begin(), criticalPoints.end(),
               std::back_inserter(saddles2),
               [](const Cell &c) -> bool { return c.dim_ == 2; });

  using Vpath = std::vector<Cell>;
  using SepSads = std::pair<Cell, Cell>;

  std::vector<std::vector<SepSads>> sepsByThread(saddles2.size());
  std::vector<std::vector<Vpath>> sepsGeomByThread(saddles2.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    const auto &s2{saddles2[i]};

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    std::set<SimplexId> saddles1{};
    VisitedMask mask{isVisited[tid], visitedTriangles[tid]};
    discreteGradient_.getDescendingWall(s2, mask, nullptr, &saddles1);

    for(const auto saddle1Id : saddles1) {
      const Cell s1{1, saddle1Id};

      Vpath vpath;
      const bool isMultiConnected
        = discreteGradient_.getAscendingPathThroughWall(
          s1, s2, isVisited[tid], &vpath);
      const auto &last = vpath.back();

      if(!isMultiConnected && last.dim_ == s2.dim_ && last.id_ == s2.id_) {
        sepsGeomByThread[i].emplace_back(std::move(vpath));
        sepsByThread[i].emplace_back(s1, s2);
      }
    }
  }

  // count total number of separatrices in sepsByThread
  std::vector<size_t> partialSepsId(sepsByThread.size() + 1, 0);

  for(size_t i = 0; i < sepsByThread.size(); ++i) {
    partialSepsId[i + 1] = partialSepsId[i] + sepsByThread[i].size();
  }

  // pre-allocate output vectors
  separatrices.resize(partialSepsId.back());
  separatricesGeometry.resize(partialSepsId.back());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sepsByThread.size(); ++i) {
    for(size_t j = 0; j < sepsByThread[i].size(); ++j) {
      const auto &sads = sepsByThread[i][j];
      const size_t k = partialSepsId[i] + j;
      separatrices[k] = Separatrix{
        true, sads.first, sads.second, false, static_cast<SimplexId>(k)};
      separatricesGeometry[k] = std::move(sepsGeomByThread[i][j]);
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
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  for(auto &vec : isVisited) {
    vec.resize(numberOfEdges, false);
  }
  std::vector<std::vector<SimplexId>> visitedEdges(this->threadNumber_);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle1 = criticalPoints[saddleIndex];

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    vector<Cell> wall;
    VisitedMask mask{isVisited[tid], visitedEdges[tid]};
    discreteGradient_.getAscendingWall(
      saddle1, mask, &wall, &separatricesSaddles[i]);

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
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  for(auto &vec : isVisited) {
    vec.resize(numberOfTriangles, false);
  }
  std::vector<std::vector<SimplexId>> visitedTriangles(this->threadNumber_);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle2 = criticalPoints[saddleIndex];

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    vector<Cell> wall;
    VisitedMask mask{isVisited[tid], visitedTriangles[tid]};
    discreteGradient_.getDescendingWall(
      saddle2, mask, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle2, emptyCell, false, i);
  }

  return 0;
}

int MorseSmaleComplex3D::getDualPolygon(const SimplexId edgeId,
                                        vector<SimplexId> &polygon) const {

  polygon.resize(inputTriangulation_->getEdgeStarNumber(edgeId));
  for(size_t i = 0; i < polygon.size(); ++i) {
    SimplexId starId;
    inputTriangulation_->getEdgeStar(edgeId, i, starId);
    polygon[i] = starId;
  }

  return 0;
}

int MorseSmaleComplex3D::sortDualPolygonVertices(
  vector<SimplexId> &polygon) const {

  for(size_t i = 1; i < polygon.size(); ++i) {

    // find polygon[i - 1] neighboring tetra in polygon[i..]
    bool isFound = false;
    size_t j = i;
    for(; j < polygon.size(); ++j) {
      // check if current is the neighbor
      for(SimplexId k = 0;
          k < inputTriangulation_->getCellNeighborNumber(polygon[i - 1]); ++k) {
        SimplexId neighborId{};
        inputTriangulation_->getCellNeighbor(polygon[i - 1], k, neighborId);
        if(neighborId == polygon[j]) {
          isFound = true;
          break;
        }
      }
      if(isFound)
        break;
    }

    // place polygon[j] next to polygon[i - 1]
    if(isFound) {
      std::swap(polygon[j], polygon[i]);
    }
  }

  return 0;
}

void MorseSmaleComplex3D::flattenSeparatricesVectors(
  std::vector<std::vector<ttk::Separatrix>> &separatrices,
  std::vector<std::vector<std::vector<ttk::dcg::Cell>>> &separatricesGeometry)
  const {

  std::vector<size_t> partialSizes{0};
  for(const auto &sep : separatrices) {
    partialSizes.emplace_back(partialSizes.back() + sep.size());
  }
  separatrices[0].resize(partialSizes.back());
  separatricesGeometry[0].resize(partialSizes.back());

  for(size_t i = 1; i < separatrices.size(); ++i) {
    const auto offset = partialSizes[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < separatrices[i].size(); ++j) {
      // shift separatrices geometry_
      for(size_t k = 0; k < separatrices[i][j].geometry_.size(); ++k) {
        separatrices[i][j].geometry_[k] += offset;
      }
      // flatten separatrices1 and separatricesGeometry1
      separatrices[0][offset + j] = std::move(separatrices[i][j]);
      separatricesGeometry[0][offset + j]
        = std::move(separatricesGeometry[i][j]);
    }
  }
}
