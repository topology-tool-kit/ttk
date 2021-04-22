/// \ingroup base
/// \class ttk::MorseSmaleComplex2D
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex2D processing package.
///
/// %MorseSmaleComplex2D is a TTK processing package that takes a scalar field
/// on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex2D.cpp %for a usage example.

#pragma once

// base code includes
#include <AbstractMorseSmaleComplex.h>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 2D triangulation.
   */
  class MorseSmaleComplex2D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex2D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename triangulationType>
    int execute(const triangulationType &triangulation);

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;
  };
} // namespace ttk

template <typename triangulationType>
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
    discreteGradient_.buildGradient<triangulationType>(triangulation);

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
    setSeparatrices1(separatrices, separatricesGeometry, triangulation);

    this->printMsg("Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getAscendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry, triangulation);
    setSeparatrices1(separatrices, separatricesGeometry, triangulation);

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

  if(outputCriticalPoints_points_ != nullptr) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints(
      criticalPoints, nCriticalPointsByDim, *outputCriticalPoints_points_,
      *outputCriticalPoints_points_cellDimensions_,
      *outputCriticalPoints_points_cellIds_,
      *outputCriticalPoints_points_isOnBoundary_,
      *outputCriticalPoints_points_PLVertexIdentifiers_, triangulation);

    if(ascendingManifold and descendingManifold) {
      discreteGradient_.setManifoldSize(
        criticalPoints, nCriticalPointsByDim, maxSeeds, ascendingManifold,
        descendingManifold, *outputCriticalPoints_points_manifoldSize_);
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
