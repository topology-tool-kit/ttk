/// \ingroup base
/// \class ttk::MorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK processing package for the computation of Morse-Smale complexes.
///
/// \b Related \b publications \n
/// "The Topology ToolKit" \n
/// Julien Tierny, Guillaume Favelier, Joshua Levine, Charles Gueunet, Michael
/// Michaux \n
/// IEEE Transactions on Visualization and Computer Graphics (Proc. of IEEE VIS
/// 2017) \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX_H
#define _MORSESMALECOMPLEX_H

// base code includes
#include <MorseSmaleComplex2D.h>
#include <MorseSmaleComplex3D.h>

/*
 * Morse-Smale complex developer quick guide:
 *
 * What is the architecture?
 * -------------------------
 * The DiscreteGradient class contains the basic structures to define and build
 * a discrete gradient. It also has several functions that decrease the number
 * of unpaired cells as parallel post-processing steps. Even more work is done
 * on the gradient with a sequential simplification step. Finally, it is able
 * to build critical points and vpaths.
 * Its files are:
 *   DiscreteGradient.cpp
 *   DiscreteGradient.h
 *
 * The AbstractMorseSmaleComplex class contains whatever is common
 * between the MorseSmaleComplex2D class and the MorseSmaleComplex3D class
 * e.g. parameters, configuration functions, input and output data pointers.
 * In particular, it contains a DiscreteGradient attribute and a Triangulation
 * attribute.
 * Its files are:
 *   AbstractMorseSmaleComplex.cpp
 *   AbstractMorseSmaleComplex.h
 *
 * The MorseSmaleComplex2D class is derived from the AbstractMorseSmaleComplex
 * class. It is specialized in building the Morse-Smale complex on 2D
 * triangulations. This class uses the DiscreteGradient attribute to build a
 * valid discrete gradient before building the MSC outputs itself i.e.
 * critical points, 1-separatrices, segmentation.
 * Its files are:
 *   MorseSmaleComplex2D.cpp
 *   MorseSmaleComplex2D.h
 *
 * The MorseSmaleComplex3D class is derived from the AbstractMorseSmaleComplex
 * class. It does the same job as the MorseSmaleComplex2D class but on 3D
 * triangulations. Note that this class has a function
 * computePersistencePairs() to get the saddle-saddle pairs of the data.
 * It adds the saddle-connectors to the 1-separatrices and adds another output
 * for 2-separatrices.
 * Its files are:
 *   MorseSmaleComplex3D.cpp
 *   MorseSmaleComplex3D.h
 *
 * The MorseSmaleComplex class is derived from the AbstractMorseSmaleComplex
 * class. It is a convenience class that detects the dimensionality of the data
 * and uses the correct implementation of the Morse-Smale complex computation
 * (between MorseSmaleComplex2D and MorseSmaleComplex3D).
 * Its files are:
 *   MorseSmaleComplex.cpp
 *   MorseSmaleComplex.h
 *
 * How to build the gradient?
 * --------------------------
 * Everything that concerns the gradient is in the DiscreteGradient class.
 * In order to build a valid discrete gradient you need to first set the data
 * pointers to the input scalar field, offset field and triangulation. You need
 * to set also the data pointers to the output critical points, 1-separatrices,
 * 2-separatrices and segmentation. Additional parameters can be configured
 * like an iteration threshold, options to have PL-Compliant extrema and
 * saddles, an option to enable collecting of persistence pairs or
 * post-processing of the saddle-connectors. Note that they all have default
 * values that correspond to a standard scenario. Like any other TTK
 * module, the level of debug and the number of threads can be adjusted to suit
 * your needs.
 * Once all the parameters and data pointers are set, the function
 * buildGradient() builds the discrete gradient. As a substantial number of
 * unpaired cells is expected, it is strongly recommended to use after this
 * function the function buildGradient2() and after that the buildGradient3()
 * function if the input dataset is in the 3D domain.
 * Finally, you can apply reverseGradient() to auto-detect the PL critical
 * points and impose that the gradient is PL-Compliant (except on the
 * boundary). See the related publication "The Topology ToolKit" for further
 * details.
 * Examples of such usage of the DiscreteGradient class can be found in
 * the execute() function of the MorseSmaleComplex2D class and
 * MorseSmaleComplex3D class as these classes need to compute a discrete
 * gradient to build their own outputs.
 *
 * Where is the simplification algorithm?
 * --------------------------------------
 * The main steps of the gradient simplification algorithm are stored in the
 * reverseGradient() function in the DiscreteGradient class. More informations
 * can be found in each simplify-like function as the process is slightly
 * different depending on the index of the critical points involved:
 *   simplifySaddleMaximumConnections() for reversal of (saddle,...,maximum)
 *   vpaths.
 *   simplifySaddleSaddleConnections1() for reversal of
 *   (2-saddle,...,1-saddle) vpaths.
 *   simplifySaddleSaddleConnections2() for reversal of
 *   (1-saddle,...,2-saddle) vpaths.
 *
 * How to add a scalar field on any output geometry?
 * -------------------------------------------------
 * First, go to AbstractMorseSmaleComplex.h and add the void* pointer to a STL
 * container (e.g. vector) as a class attribute. In the code, the attributes
 * of the same output are grouped together and are prefixed by its name.
 * So for example, the outputSeparatrices1_points_smoothingMask_ variable
 * represents the smoothingMask scalar field that is associated to the points
 * of the 1-separatrices of the Morse-Smale complex. As you just added a new
 * attribute, you need to update the setter corresponding to the output to
 * add this new element:
 * The setters available are:
 *   setOutputCriticalPoints()
 *   setOutputSeparatrices1()
 *   setOutputSeparatrices2()
 *   setOutputMorseComplexes()
 * As you added a new attribute in the class, you need to give it a default
 * value in the constructor in AbstractMorseSmaleComplex.cpp (typically nullptr
 * for a pointer).
 * Now, you need to overload the same setter than previously but in
 * MorseSmaleComplex.h this time in order to propagate the new field
 * information to the concrete implementations. As the MorseSmaleComplex2D and
 * MorseSmaleComplex3D classes are both derived from AbstractMorseSmaleComplex
 * they already have the updated version of the function.
 * Now that you have access to the data pointer inside the actual
 * implementation you can cast it from void* to dataType* and modify it to
 * your convenience.
 *
 * What part of the code is parallel?
 * ----------------------------------
 * From a global point a view, the part of the code that builds the gradient
 * as well as the two post-processing steps are accelerated by OpenMP (if
 * available).
 * The gradient simplification step is mostly sequential except the
 * initialization of internal structures in the initialize-like functions which
 * is done is parallel with OpenMP.
 * Then, each output of the Morse-Smale complex is computed in parallel with
 * OpenMP except the 2-separatrices of a 3D dataset that require heavy
 * synchronisation.
 * Finally, the generation of the geometry for the visualization is done
 * sequentially.
 * Complete list:
 * In the DiscreteGradient class:
 *   buildGradient()
 *   buildGradient2()
 *   buildGradient3()
 *   initializeSaddleMaximumConnections()
 *   initializeSaddleSaddleConnections1()
 *   initializeSaddleSaddleConnections2()
 * In the MorseSmaleComplex2D class:
 *   getSeparatrices1()
 *   setAscendingSegmentation()
 *   setDescendingSegmentation()
 * In the MorseSmaleComplex3D class:
 *   getSeparatrices1()
 *   getAscendingSeparatrices2()
 *   getDescendingSeparatrices2()
 *   setAscendingSegmentation()
 *   setDescendingSegmentation()
 */

namespace ttk {

  class MorseSmaleComplex {

  public:
    MorseSmaleComplex();
    ~MorseSmaleComplex();

    int setIterationThreshold(const int iterationThreshold) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setIterationThreshold(
        iterationThreshold);
    }

    int setReverseSaddleMaximumConnection(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setReverveSaddleMaximumConnection(
        state);
    }

    int setReverseSaddleSaddleConnection(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setReverveSaddleSaddleConnection(
        state);
    }

    int setComputeAscendingSeparatrices1(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setComputeAscendingSeparatrices1(
        state);
    }

    int setComputeDescendingSeparatrices1(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setComputeDescendingSeparatrices1(
        state);
    }

    int setComputeSaddleConnectors(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setComputeSaddleConnectors(state);
    }

    int setComputeAscendingSeparatrices2(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setComputeAscendingSeparatrices2(
        state);
    }

    int setComputeDescendingSeparatrices2(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setComputeDescendingSeparatrices2(
        state);
    }

    int setReturnSaddleConnectors(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setReturnSaddleConnectors(state);
    }

    int setSaddleConnectorsPersistenceThreshold(const double threshold) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_
        ->setSaddleConnectorsPersistenceThreshold(threshold);
    }

    int setPrioritizeSpeedOverMemory(const bool state) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      return abstractMorseSmaleComplex_->setPrioritizeSpeedOverMemory(state);
    }

    int setupTriangulation(Triangulation *const data) {
      dimensionality_ = data->getCellVertexNumber(0) - 1;

      switch(dimensionality_) {
        case 2:
          abstractMorseSmaleComplex_ = &morseSmaleComplex2D_;
          break;

        case 3:
          abstractMorseSmaleComplex_ = &morseSmaleComplex3D_;
          break;
      }

      abstractMorseSmaleComplex_->setupTriangulation(data);
      return 0;
    }

    inline int setDebugLevel(const int &debugLevel) {
      morseSmaleComplex2D_.setDebugLevel(debugLevel);
      morseSmaleComplex3D_.setDebugLevel(debugLevel);
      return 0;
    }

    inline int setThreadNumber(const int &threadNumber) {
      morseSmaleComplex2D_.setThreadNumber(threadNumber);
      morseSmaleComplex3D_.setThreadNumber(threadNumber);
      return 0;
    }

    inline int setWrapper(const Wrapper *const wrapper) {
      morseSmaleComplex2D_.setWrapper(wrapper);
      morseSmaleComplex3D_.setWrapper(wrapper);
      return 0;
    }

    inline int setInputScalarField(void *const data) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setInputScalarField(data);
      return 0;
    }

    inline int setInputOffsets(void *const data) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setInputOffsets(data);
      return 0;
    }

    inline int setOutputCriticalPoints(
      SimplexId *const criticalPoints_numberOfPoints,
      std::vector<float> *const criticalPoints_points,
      std::vector<char> *const criticalPoints_points_cellDimensons,
      std::vector<SimplexId> *const criticalPoints_points_cellIds,
      void *const criticalPoints_points_cellScalars,
      std::vector<char> *const criticalPoints_points_isOnBoundary,
      std::vector<SimplexId> *const criticalPoints_points_PLVertexIdentifiers,
      std::vector<SimplexId> *const criticalPoints_points_manifoldSize) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setOutputCriticalPoints(
        criticalPoints_numberOfPoints, criticalPoints_points,
        criticalPoints_points_cellDimensons, criticalPoints_points_cellIds,
        criticalPoints_points_cellScalars, criticalPoints_points_isOnBoundary,
        criticalPoints_points_PLVertexIdentifiers,
        criticalPoints_points_manifoldSize);

      return 0;
    }

    inline int setOutputSeparatrices1(
      SimplexId *const separatrices1_numberOfPoints,
      std::vector<float> *const separatrices1_points,
      std::vector<char> *const separatrices1_points_smoothingMask,
      std::vector<char> *const separatrices1_points_cellDimensions,
      std::vector<SimplexId> *const separatrices1_points_cellIds,
      SimplexId *const separatrices1_numberOfCells,
      std::vector<SimplexId> *const separatrices1_cells,
      std::vector<SimplexId> *const separatrices1_cells_sourceIds,
      std::vector<SimplexId> *const separatrices1_cells_destinationIds,
      std::vector<SimplexId> *const separatrices1_cells_separatrixIds,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      void *const separatrices1_cells_separatrixFunctionMaxima,
      void *const separatrices1_cells_separatrixFunctionMinima,
      void *const separatrices1_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices1_cells_isOnBoundary) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setOutputSeparatrices1(
        separatrices1_numberOfPoints, separatrices1_points,
        separatrices1_points_smoothingMask, separatrices1_points_cellDimensions,
        separatrices1_points_cellIds, separatrices1_numberOfCells,
        separatrices1_cells, separatrices1_cells_sourceIds,
        separatrices1_cells_destinationIds, separatrices1_cells_separatrixIds,
        separatrices1_cells_separatrixTypes,
        separatrices1_cells_separatrixFunctionMaxima,
        separatrices1_cells_separatrixFunctionMinima,
        separatrices1_cells_separatrixFunctionDiffs,
        separatrices1_cells_isOnBoundary);

      return 0;
    }

    inline int setOutputSeparatrices2(
      SimplexId *const separatrices2_numberOfPoints,
      std::vector<float> *const separatrices2_points,
      SimplexId *const separatrices2_numberOfCells,
      std::vector<SimplexId> *const separatrices2_cells,
      std::vector<SimplexId> *const separatrices2_cells_sourceIds,
      std::vector<SimplexId> *const separatrices2_cells_separatrixIds,
      std::vector<char> *const separatrices2_cells_separatrixTypes,
      void *const separatrices2_cells_separatrixFunctionMaxima,
      void *const separatrices2_cells_separatrixFunctionMinima,
      void *const separatrices2_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices2_cells_isOnBoundary) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setOutputSeparatrices2(
        separatrices2_numberOfPoints, separatrices2_points,
        separatrices2_numberOfCells, separatrices2_cells,
        separatrices2_cells_sourceIds, separatrices2_cells_separatrixIds,
        separatrices2_cells_separatrixTypes,
        separatrices2_cells_separatrixFunctionMaxima,
        separatrices2_cells_separatrixFunctionMinima,
        separatrices2_cells_separatrixFunctionDiffs,
        separatrices2_cells_isOnBoundary);

      return 0;
    }

    inline int setOutputMorseComplexes(void *const ascendingManifold,
                                       void *const descendingManifold,
                                       void *const morseSmaleManifold) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!abstractMorseSmaleComplex_) {
        return -1;
      }
#endif
      abstractMorseSmaleComplex_->setOutputMorseComplexes(
        ascendingManifold, descendingManifold, morseSmaleManifold);

      return 0;
    }

    template <typename dataType, typename idType>
    int execute() {
      switch(dimensionality_) {
        case 2:
          morseSmaleComplex2D_.execute<dataType, idType>();
          break;

        case 3:
          morseSmaleComplex3D_.execute<dataType, idType>();
          break;
      }
      return 0;
    }

  protected:
    int dimensionality_;
    AbstractMorseSmaleComplex *abstractMorseSmaleComplex_;
    MorseSmaleComplex2D morseSmaleComplex2D_;
    MorseSmaleComplex3D morseSmaleComplex3D_;
  };
} // namespace ttk

#endif // MORSESMALECOMPLEX_H
