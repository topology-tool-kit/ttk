/// \ingroup base
/// \class ttk::MorseSmaleComplex2D
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex2D processing package.
///
/// %MorseSmaleComplex2D is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex2D.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX2D_H
#define _MORSESMALECOMPLEX2D_H

// base code includes
#include<AbstractMorseSmaleComplex.h>

namespace ttk{

  /**
   * Class specialized in building the Morse-Smale complex
   * of 2D triangulation.
   */
  class MorseSmaleComplex2D : public AbstractMorseSmaleComplex{

    public:

      MorseSmaleComplex2D();
      ~MorseSmaleComplex2D();

      /**
       * Main function for computing the whole Morse-Smale complex.
       */
      template<typename dataType>
        int execute();

      /**
       * Compute the descending 1-separatrices by reading into the discrete
       * gradient.
       */
      int getAscendingSeparatrices1(const std::vector<Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<Cell>>& separatricesGeometry) const;

  };
}

template<typename dataType>
int ttk::MorseSmaleComplex2D::execute(){
  Timer t;

  int* ascendingManifold=static_cast<int*>(outputAscendingManifold_);
  int* descendingManifold=static_cast<int*>(outputDescendingManifold_);
  int* morseSmaleManifold=static_cast<int*>(outputMorseSmaleManifold_);

  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.setDebugLevel(debugLevel_);
  discreteGradient_.buildGradient<dataType>();
  discreteGradient_.buildGradient2<dataType>();
  discreteGradient_.reverseGradient<dataType>();

  std::vector<Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    getDescendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Descending 1-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices1){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    getAscendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Ascending 1-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  std::vector<int> maxSeeds;
  {
    Timer tmp;

    int numberOfMaxima{};
    int numberOfMinima{};

    if(ComputeAscendingSegmentation)
      setAscendingSegmentation(criticalPoints, maxSeeds, ascendingManifold, numberOfMaxima);

    if(ComputeDescendingSegmentation)
      setDescendingSegmentation(criticalPoints, descendingManifold, numberOfMinima);

    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation and ComputeFinalSegmentation)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold, descendingManifold, morseSmaleManifold);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Segmentation computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation)
    discreteGradient_.setAugmentedCriticalPoints<dataType>(criticalPoints,
        maxSeeds,
        ascendingManifold,
        descendingManifold);
  else
    discreteGradient_.setCriticalPoints<dataType>(criticalPoints);

  {
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    std::stringstream msg;
    msg << "[MorseSmaleComplex2D] Data-set (" << numberOfVertices
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // MORSESMALECOMPLEX2D_H
