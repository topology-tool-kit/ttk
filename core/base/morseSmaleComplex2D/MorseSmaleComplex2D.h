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

#ifndef _MORSESMALECOMPLEX2D_H
#define _MORSESMALECOMPLEX2D_H

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
    std::cerr
      << "[MorseSmaleComplex2D] Error: input scalar field pointer is null."
      << std::endl;
    return -1;
  }

  if(!inputOffsets_) {
    std::cerr
      << "[MorseSmaleComplex2D] Error: input offset field pointer is null."
      << std::endl;
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
    discreteGradient_.buildGradient<dataType, idType>();
    discreteGradient_.buildGradient2<dataType, idType>();

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex2D] Discrete gradient overall computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }
  discreteGradient_.reverseGradient<dataType, idType>();

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getDescendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex2D] Descending 1-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getAscendingSeparatrices1(
      criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex2D] Ascending 1-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
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
      std::stringstream msg;
      msg << "[MorseSmaleComplex2D] Segmentation computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(outputCriticalPoints_numberOfPoints_ and outputCriticalPoints_points_) {
    if(ascendingManifold and descendingManifold)
      discreteGradient_.setAugmentedCriticalPoints<dataType, idType>(
        criticalPoints, maxSeeds, ascendingManifold, descendingManifold);
    else
      discreteGradient_.setCriticalPoints<dataType, idType>(criticalPoints);
  }

  {
    const SimplexId numberOfVertices
      = inputTriangulation_->getNumberOfVertices();
    std::stringstream msg;
    msg << "[MorseSmaleComplex2D] Data-set (" << numberOfVertices
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // MORSESMALECOMPLEX2D_H
