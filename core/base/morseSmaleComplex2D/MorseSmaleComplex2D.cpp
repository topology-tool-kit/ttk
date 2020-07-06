#include <MorseSmaleComplex2D.h>

using namespace std;
using namespace ttk;
using namespace dcg;

MorseSmaleComplex2D::MorseSmaleComplex2D() : AbstractMorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex2D");
}

MorseSmaleComplex2D::~MorseSmaleComplex2D() {
}

template <typename dataType, typename idType, typename triangulationType>
int ttk::MorseSmaleComplex2D::execute(const triangulationType &triangulation) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField_) {
    this->printErr("Input scalar field pointer is null.");
    return -1;
  }

  if(!inputOffsets_) {
    this->printErr("Input offset field pointer is null.");
    return -1;
  }
#endif
  Timer t;

  // nullptr_t is implicitly convertible and comparable to any pointer type
  // or pointer-to-member type.
  SimplexId *ascendingManifold
    = static_cast<SimplexId *>(outputAscendingManifold_);
  SimplexId *descendingManifold
    = static_cast<SimplexId *>(outputDescendingManifold_);
  SimplexId *morseSmaleManifold
    = static_cast<SimplexId *>(outputMorseSmaleManifold_);

  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.setDebugLevel(debugLevel_);
  {
    Timer tmp;
    discreteGradient_.buildGradient<dataType, idType>(triangulation);

    this->printMsg("Discrete gradient computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints, triangulation);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getDescendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry, triangulation);
    setSeparatrices1<dataType>(
      separatrices, separatricesGeometry, triangulation);

    this->printMsg("Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getAscendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry, triangulation);
    setSeparatrices1<dataType>(
      separatrices, separatricesGeometry, triangulation);

    this->printMsg("Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  std::vector<SimplexId> maxSeeds;
  {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ascendingManifold)
      setAscendingSegmentation(criticalPoints, maxSeeds, ascendingManifold,
                               numberOfMaxima, triangulation);

    if(descendingManifold)
      setDescendingSegmentation(
        criticalPoints, descendingManifold, numberOfMinima, triangulation);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold,
                           descendingManifold, morseSmaleManifold,
                           triangulation);

    if(ascendingManifold or descendingManifold) {
      this->printMsg("Segmentation computed", 1.0, tmp.getElapsedTime(),
                     this->threadNumber_);
    }
  }

  if(outputCriticalPoints_numberOfPoints_ and outputCriticalPoints_points_) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints<dataType>(
      criticalPoints, nCriticalPointsByDim, triangulation);

    discreteGradient_.fetchOutputCriticalPoints(
      outputCriticalPoints_numberOfPoints_, outputCriticalPoints_points_,
      outputCriticalPoints_points_cellDimensions_,
      outputCriticalPoints_points_cellIds_,
      outputCriticalPoints_points_isOnBoundary_,
      outputCriticalPoints_points_PLVertexIdentifiers_);

    if(ascendingManifold and descendingManifold) {
      discreteGradient_.setManifoldSize(criticalPoints, nCriticalPointsByDim,
                                        maxSeeds, ascendingManifold,
                                        descendingManifold);
      discreteGradient_.fetchOutputManifoldSize(
        outputCriticalPoints_points_manifoldSize_);
    }
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex2D::getAscendingSeparatrices1(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  std::vector<SimplexId> saddleIndexes;
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

    // add ascending vpaths
    {
      const SimplexId starNumber = triangulation.getEdgeStarNumber(saddle.id_);
      for(SimplexId j = 0; j < starNumber; ++j) {
        const SimplexId shift = j;

        SimplexId triangleId;
        triangulation.getEdgeStar(saddle.id_, j, triangleId);

        std::vector<Cell> vpath;
        vpath.push_back(saddle);
        discreteGradient_.getAscendingPath(
          Cell(2, triangleId), vpath, triangulation);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 2 and discreteGradient_.isCellCritical(lastCell)) {
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

// explicit specializations to reduce compile time
#define MSC2D_SPECIALIZE(DATATYPE)                                        \
  template int ttk::MorseSmaleComplex2D::execute<DATATYPE, SimplexId,     \
                                                 ExplicitTriangulation>(  \
    const ExplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex2D::execute<DATATYPE, SimplexId,     \
                                                 ImplicitTriangulation>(  \
    const ImplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex2D::execute<                         \
    DATATYPE, SimplexId, PeriodicImplicitTriangulation>(                  \
    const PeriodicImplicitTriangulation &);                               \
  template int ttk::MorseSmaleComplex2D::execute<DATATYPE, LongSimplexId, \
                                                 ExplicitTriangulation>(  \
    const ExplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex2D::execute<DATATYPE, LongSimplexId, \
                                                 ImplicitTriangulation>(  \
    const ImplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex2D::execute<                         \
    DATATYPE, LongSimplexId, PeriodicImplicitTriangulation>(              \
    const PeriodicImplicitTriangulation &)

MSC2D_SPECIALIZE(double);
MSC2D_SPECIALIZE(float);
MSC2D_SPECIALIZE(long long);
MSC2D_SPECIALIZE(unsigned long long);
MSC2D_SPECIALIZE(long);
MSC2D_SPECIALIZE(unsigned long);
MSC2D_SPECIALIZE(int);
MSC2D_SPECIALIZE(unsigned int);
MSC2D_SPECIALIZE(short);
MSC2D_SPECIALIZE(unsigned short);
MSC2D_SPECIALIZE(char);
MSC2D_SPECIALIZE(signed char);
MSC2D_SPECIALIZE(unsigned char);
