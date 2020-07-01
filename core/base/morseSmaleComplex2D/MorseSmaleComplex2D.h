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
    ~MorseSmaleComplex2D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename dataType, typename idType>
    int execute();

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const;
  };
} // namespace ttk

template <typename dataType, typename idType>
int ttk::MorseSmaleComplex2D::execute() {
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
    discreteGradient_.buildGradient<dataType, idType>(*inputTriangulation_);

    this->printMsg("Discrete gradient computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints, *inputTriangulation_);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getDescendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    this->printMsg("Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getAscendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    this->printMsg("Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  std::vector<SimplexId> maxSeeds;
  {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ascendingManifold)
      setAscendingSegmentation(
        criticalPoints, maxSeeds, ascendingManifold, numberOfMaxima);

    if(descendingManifold)
      setDescendingSegmentation(
        criticalPoints, descendingManifold, numberOfMinima);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold,
                           descendingManifold, morseSmaleManifold);

    if(ascendingManifold or descendingManifold) {
      this->printMsg("Segmentation computed", 1.0, tmp.getElapsedTime(),
                     this->threadNumber_);
    }
  }

  if(outputCriticalPoints_numberOfPoints_ and outputCriticalPoints_points_) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints<dataType>(
      criticalPoints, nCriticalPointsByDim, *inputTriangulation_);

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
                   + std::to_string(inputTriangulation_->getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}
