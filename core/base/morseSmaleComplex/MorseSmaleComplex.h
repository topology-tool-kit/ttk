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

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <Triangulation.h>

#include <queue>

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
  class MorseSmaleComplex : public virtual Debug {
  public:
    MorseSmaleComplex();

    /**
     * Main function for computing the Morse-Smale complex.
     */
    template <typename dataType, typename triangulationType>
    inline int execute(const triangulationType &triangulation) {
      switch(triangulation.getDimensionality()) {
        case 1:
        case 2:
          return this->execute2D(triangulation);
        case 3:
          return this->execute3D<dataType>(triangulation);
      }
      return 0;
    }

    /**
     * Main function for computing the Morse-Smale complex on 2D triangulations.
     */
    template <typename triangulationType>
    int execute2D(const triangulationType &triangulation);

    /**
     * Main function for computing the Morse-Smale complex on 3D triangulations.
     */
    template <typename dataType, typename triangulationType>
    int execute3D(const triangulationType &triangulation);

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the ascending manifolds of the critical points
     * (disabled by default).
     */
    inline void setComputeAscendingSeparatrices1(const bool state) {
      ComputeAscendingSeparatrices1 = state;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the descending manifolds of the critical points
     * (disabled by default).
     */
    inline void setComputeDescendingSeparatrices1(const bool state) {
      ComputeDescendingSeparatrices1 = state;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the visible saddle-connectors (disabled by default).
     */
    inline void setComputeSaddleConnectors(const bool state) {
      ComputeSaddleConnectors = state;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the ascending 2-separatrices (disabled by default).
     */
    inline void setComputeAscendingSeparatrices2(const bool state) {
      ComputeAscendingSeparatrices2 = state;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the descending 2-separatrices (disabled by default).
     */
    inline void setComputeDescendingSeparatrices2(const bool state) {
      ComputeDescendingSeparatrices2 = state;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      discreteGradient_.preconditionTriangulation(data);
      data->preconditionCellEdges();
      data->preconditionCellNeighbors();
    }

    /**
     * Set the input scalar field associated on the points of the data set.
     */
    inline void setInputScalarField(const void *const data) {
      inputScalarField_ = data;
      discreteGradient_.setInputScalarField(inputScalarField_);
    }

    /**
     * Set the input offset field associated on the points of the data set
     * (if none, identifiers are used instead).
     *
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p data buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline void setInputOffsets(const SimplexId *const data) {
      inputOffsets_ = data;
      discreteGradient_.setInputOffsets(inputOffsets_);
    }

    /**
     * Set the output critical points data pointers.
     */
    inline void setOutputCriticalPoints(
      std::vector<std::array<float, 3>> *const criticalPoints_points,
      std::vector<char> *const criticalPoints_points_cellDimensons,
      std::vector<SimplexId> *const criticalPoints_points_cellIds,
      std::vector<char> *const criticalPoints_points_isOnBoundary,
      std::vector<SimplexId> *const criticalPoints_points_PLVertexIdentifiers,
      std::vector<SimplexId> *criticalPoints_points_manifoldSize) {
      outputCriticalPoints_points_ = criticalPoints_points;
      outputCriticalPoints_points_cellDimensions_
        = criticalPoints_points_cellDimensons;
      outputCriticalPoints_points_cellIds_ = criticalPoints_points_cellIds;
      outputCriticalPoints_points_isOnBoundary_
        = criticalPoints_points_isOnBoundary;
      outputCriticalPoints_points_PLVertexIdentifiers_
        = criticalPoints_points_PLVertexIdentifiers;
      outputCriticalPoints_points_manifoldSize_
        = criticalPoints_points_manifoldSize;
    }

    /**
     * Set the data pointers to the output 1-separatrices.
     */
    inline void setOutputSeparatrices1(
      SimplexId *const separatrices1_numberOfPoints,
      std::vector<float> *const separatrices1_points,
      std::vector<char> *const separatrices1_points_smoothingMask,
      std::vector<char> *const separatrices1_points_cellDimensions,
      std::vector<SimplexId> *const separatrices1_points_cellIds,
      SimplexId *const separatrices1_numberOfCells,
      std::vector<SimplexId> *const separatrices1_cells_connectivity,
      std::vector<SimplexId> *const separatrices1_cells_sourceIds,
      std::vector<SimplexId> *const separatrices1_cells_destinationIds,
      std::vector<SimplexId> *const separatrices1_cells_separatrixIds,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      std::vector<SimplexId> *const s1_cells_separatrixFunctionMaximaId,
      std::vector<SimplexId> *const s1_cells_separatrixFunctionMinimaId,
      std::vector<char> *const separatrices1_cells_isOnBoundary) {
      outputSeparatrices1_numberOfPoints_ = separatrices1_numberOfPoints;
      outputSeparatrices1_points_ = separatrices1_points;
      outputSeparatrices1_points_smoothingMask_
        = separatrices1_points_smoothingMask;
      outputSeparatrices1_points_cellDimensions_
        = separatrices1_points_cellDimensions;
      outputSeparatrices1_points_cellIds_ = separatrices1_points_cellIds;
      outputSeparatrices1_numberOfCells_ = separatrices1_numberOfCells;
      outputSeparatrices1_cells_connectivity_
        = separatrices1_cells_connectivity;
      outputSeparatrices1_cells_sourceIds_ = separatrices1_cells_sourceIds;
      outputSeparatrices1_cells_destinationIds_
        = separatrices1_cells_destinationIds;
      outputSeparatrices1_cells_separatrixIds_
        = separatrices1_cells_separatrixIds;
      outputSeparatrices1_cells_separatrixTypes_
        = separatrices1_cells_separatrixTypes;
      outputS1_cells_separatrixFunctionMaximaId_
        = s1_cells_separatrixFunctionMaximaId;
      outputS1_cells_separatrixFunctionMinimaId_
        = s1_cells_separatrixFunctionMinimaId;
      outputSeparatrices1_cells_isOnBoundary_
        = separatrices1_cells_isOnBoundary;
    }

    /**
     * Set the data pointers to the output 2-separatrices.
     */
    inline void setOutputSeparatrices2(
      SimplexId *const separatrices2_numberOfPoints,
      std::vector<float> *const separatrices2_points,
      SimplexId *const separatrices2_numberOfCells,
      std::vector<SimplexId> *const separatrices2_cells_offsets,
      std::vector<SimplexId> *const separatrices2_cells_connectivity,
      std::vector<SimplexId> *const separatrices2_cells_sourceIds,
      std::vector<SimplexId> *const separatrices2_cells_separatrixIds,
      std::vector<char> *const separatrices2_cells_separatrixTypes,
      std::vector<SimplexId> *const s2_cells_separatrixFunctionMaximaId,
      std::vector<SimplexId> *const s2_cells_separatrixFunctionMinimaId,
      std::vector<char> *const separatrices2_cells_isOnBoundary) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_offsets_ = separatrices2_cells_offsets;
      outputSeparatrices2_cells_connectivity_
        = separatrices2_cells_connectivity;
      outputSeparatrices2_cells_sourceIds_ = separatrices2_cells_sourceIds;
      outputSeparatrices2_cells_separatrixIds_
        = separatrices2_cells_separatrixIds;
      outputSeparatrices2_cells_separatrixTypes_
        = separatrices2_cells_separatrixTypes;
      outputS2_cells_separatrixFunctionMaximaId_
        = s2_cells_separatrixFunctionMaximaId;
      outputS2_cells_separatrixFunctionMinimaId_
        = s2_cells_separatrixFunctionMinimaId;
      outputSeparatrices2_cells_isOnBoundary_
        = separatrices2_cells_isOnBoundary;
    }

    /**
     * Set the data pointers to the output segmentation scalar fields.
     */
    inline void setOutputMorseComplexes(void *const ascendingManifold,
                                        void *const descendingManifold,
                                        void *const morseSmaleManifold) {
      outputAscendingManifold_ = ascendingManifold;
      outputDescendingManifold_ = descendingManifold;
      outputMorseSmaleManifold_ = morseSmaleManifold;
    }

  protected:
    /**
     * Utility class representing Ridge lines, Valley lines
     * and Saddle connectors.
     */
    struct Separatrix {
      // default :
      explicit Separatrix()
        : isValid_{}, source_{}, destination_{}, isReversed_{}, geometry_{} {
      }

      // initialization with one segment :
      explicit Separatrix(const bool isValid,
                          const dcg::Cell &saddle,
                          const dcg::Cell &extremum,
                          const bool isSegmentReversed,
                          const SimplexId segmentGeometry)
        : isValid_{isValid}, source_{saddle}, destination_{extremum} {
        isReversed_.push_back(isSegmentReversed);
        geometry_.push_back(segmentGeometry);
      }

      // initialization with multiple segments :
      explicit Separatrix(const bool isValid,
                          const dcg::Cell &saddle,
                          const dcg::Cell &extremum,
                          const std::vector<char> &isReversed,
                          const std::vector<SimplexId> &geometry)
        : isValid_{isValid}, source_{saddle}, destination_{extremum},
          isReversed_{isReversed}, geometry_{geometry} {
      }

      explicit Separatrix(const Separatrix &separatrix)
        : isValid_{separatrix.isValid_}, source_{separatrix.source_},
          destination_{separatrix.destination_},
          isReversed_{separatrix.isReversed_}, geometry_{separatrix.geometry_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : isValid_{separatrix.isValid_}, source_{separatrix.source_},
          destination_{separatrix.destination_}, isReversed_{std::move(
                                                   separatrix.isReversed_)},
          geometry_{std::move(separatrix.geometry_)} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        isValid_ = separatrix.isValid_;
        source_ = separatrix.source_;
        destination_ = separatrix.destination_;
        isReversed_ = std::move(separatrix.isReversed_);
        geometry_ = std::move(separatrix.geometry_);

        return *this;
      }

      /**
       * Flag indicating if this separatrix can be processed.
       */
      bool isValid_;

      /**
       * Source cell of the separatrix.
       */
      dcg::Cell source_;

      /**
       * Destination cell of the separatrix.
       */
      dcg::Cell destination_;

      /**
       * Container of flags, isReversed[i] indicates if the
       * element stored at id=geometry_[i] can be reversed.
       */
      std::vector<char> isReversed_;

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;
    };

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the ascending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the 1-separatrices. This
     * function needs the following internal pointers to be set:
     * outputSeparatrices1_numberOfPoints_
     * outputSeparatrices1_points_
     * outputSeparatrices1_numberOfCells_
     * outputSeparatrices1_cells_
     * inputScalarField_
     */
    template <typename triangulationType>
    int setSeparatrices1(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the ascending manifold of the maxima.
     */
    template <typename triangulationType>
    int setAscendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                 std::vector<SimplexId> &maxSeeds,
                                 SimplexId *const morseSmaleManifold,
                                 SimplexId &numberOfMaxima,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the descending manifold of the minima.
     */
    template <typename triangulationType>
    int setDescendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                  SimplexId *const morseSmaleManifold,
                                  SimplexId &numberOfMinima,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the final combinatorial Morse-Smale complex
     * segmentation.
     */
    template <typename triangulationType>
    int setFinalSegmentation(const SimplexId numberOfMaxima,
                             const SimplexId numberOfMinima,
                             const SimplexId *const ascendingManifold,
                             const SimplexId *const descendingManifold,
                             SimplexId *const morseSmaleManifold,
                             const triangulationType &triangulation) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getSaddleConnectors(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the maxima.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the descending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename triangulationType>
    int setDescendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Find all tetras in the star of edgeId
     *
     * (primal: star of edgeId -> dual: vertices of polygon)
     */
    template <typename triangulationType>
    int getDualPolygon(const SimplexId edgeId,
                       SimplexId *const polygon,
                       const size_t polSize,
                       const triangulationType &triangulation) const;

    /**
     * Sort the polygon vertices to be clockwise
     */
    template <typename triangulationType>
    int sortDualPolygonVertices(SimplexId *const polygon,
                                const size_t polSize,
                                const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the minima.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the ascending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename triangulationType>
    int setAscendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<Separatrix>> &separatrices,
      std::vector<std::vector<std::vector<ttk::dcg::Cell>>>
        &separatricesGeometry) const;

  protected:
    bool ComputeAscendingSeparatrices1{true};
    bool ComputeDescendingSeparatrices1{true};
    bool ComputeSaddleConnectors{false};
    bool ComputeAscendingSeparatrices2{false};
    bool ComputeDescendingSeparatrices2{false};

    dcg::DiscreteGradient discreteGradient_{};

    const void *inputScalarField_{};
    const SimplexId *inputOffsets_{};

    std::vector<std::array<float, 3>> *outputCriticalPoints_points_{};
    std::vector<char> *outputCriticalPoints_points_cellDimensions_{};
    std::vector<SimplexId> *outputCriticalPoints_points_cellIds_{};
    std::vector<char> *outputCriticalPoints_points_isOnBoundary_{};
    std::vector<SimplexId> *outputCriticalPoints_points_PLVertexIdentifiers_{};
    std::vector<SimplexId> *outputCriticalPoints_points_manifoldSize_{};

    SimplexId *outputSeparatrices1_numberOfPoints_{};
    std::vector<float> *outputSeparatrices1_points_{};
    std::vector<char> *outputSeparatrices1_points_smoothingMask_{};
    std::vector<char> *outputSeparatrices1_points_cellDimensions_{};
    std::vector<SimplexId> *outputSeparatrices1_points_cellIds_{};
    SimplexId *outputSeparatrices1_numberOfCells_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_connectivity_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_sourceIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_destinationIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_separatrixIds_{};
    std::vector<char> *outputSeparatrices1_cells_separatrixTypes_{};
    std::vector<SimplexId> *outputS1_cells_separatrixFunctionMaximaId_{};
    std::vector<SimplexId> *outputS1_cells_separatrixFunctionMinimaId_{};
    std::vector<char> *outputSeparatrices1_cells_isOnBoundary_{};

    SimplexId *outputSeparatrices2_numberOfPoints_{};
    std::vector<float> *outputSeparatrices2_points_{};
    SimplexId *outputSeparatrices2_numberOfCells_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_offsets_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_connectivity_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_sourceIds_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_separatrixIds_{};
    std::vector<char> *outputSeparatrices2_cells_separatrixTypes_{};
    std::vector<SimplexId> *outputS2_cells_separatrixFunctionMaximaId_{};
    std::vector<SimplexId> *outputS2_cells_separatrixFunctionMinimaId_{};
    std::vector<char> *outputSeparatrices2_cells_isOnBoundary_{};

    void *outputAscendingManifold_{};
    void *outputDescendingManifold_{};
    void *outputMorseSmaleManifold_{};
  };
} // namespace ttk

template <typename triangulationType>
int ttk::MorseSmaleComplex::getDescendingSeparatrices1(
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

    // add descending vpaths
    {
      const Cell &saddle1 = saddle;

      for(int j = 0; j < 2; ++j) {
        const int shift = j + 2;

        SimplexId vertexId;
        triangulation.getEdgeVertex(saddle1.id_, j, vertexId);

        std::vector<Cell> vpath;
        vpath.push_back(saddle1);
        discreteGradient_.getDescendingPath(
          Cell(0, vertexId), vpath, triangulation);

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

template <typename triangulationType>
int ttk::MorseSmaleComplex::setSeparatrices1(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices1_numberOfPoints_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices1_points_ == nullptr) {
    this->printErr("1-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices1_numberOfCells_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices1_cells_connectivity_ == nullptr) {
    this->printErr("1-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr(
      " 1-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  auto separatrixFunctionMaxima = outputS1_cells_separatrixFunctionMaximaId_;
  auto separatrixFunctionMinima = outputS1_cells_separatrixFunctionMinimaId_;

  // max existing separatrix id + 1 or 0
  const SimplexId separatrixId
    = (outputSeparatrices1_cells_separatrixIds_ != nullptr
       && !outputSeparatrices1_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
                            outputSeparatrices1_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices1_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices1_numberOfCells_)};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // points beginning id for each separatrix geometry
  std::vector<size_t> geomPointsBegId{npoints};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      const auto sepSize = separatricesGeometry[geomId].size();
      npoints += sepSize;
      ncells += sepSize - 1;
      geomPointsBegId.emplace_back(npoints);
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;

  // resize arrays
  outputSeparatrices1_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize(2 * ncells);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  if(outputSeparatrices1_points_smoothingMask_ != nullptr)
    outputSeparatrices1_points_smoothingMask_->resize(npoints);
  if(outputSeparatrices1_points_cellDimensions_ != nullptr)
    outputSeparatrices1_points_cellDimensions_->resize(npoints);
  if(outputSeparatrices1_points_cellIds_ != nullptr)
    outputSeparatrices1_points_cellIds_->resize(npoints);
  if(outputSeparatrices1_cells_sourceIds_ != nullptr)
    outputSeparatrices1_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices1_cells_destinationIds_ != nullptr)
    outputSeparatrices1_cells_destinationIds_->resize(ncells);
  if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
    outputSeparatrices1_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
    outputSeparatrices1_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(separatrixId + validGeomIds.size());
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(separatrixId + validGeomIds.size());
  if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
    outputSeparatrices1_cells_isOnBoundary_->resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    // saddle (asc/desc sep) or saddle1 (saddle connector)
    const dcg::Cell &src = sep.source_;
    // extremum (asc/desc sep) or saddle2 (saddle connector)
    const dcg::Cell &dst = sep.destination_;

    // get separatrix type
    const auto saddleConnector
      = dimensionality == 3 && src.dim_ == 1 && dst.dim_ == 2;
    const char sepType
      = saddleConnector ? 1 : std::min(dst.dim_, dimensionality - 1);

    // compute separatrix function diff
    const auto vertsOrder = [this](const SimplexId a, const SimplexId b) {
      return this->inputOffsets_[a] < this->inputOffsets_[b];
    };
    const std::array<SimplexId, 2> gVerts{
      discreteGradient_.getCellGreaterVertex(src, triangulation),
      discreteGradient_.getCellGreaterVertex(dst, triangulation)};
    const auto sepFuncMax
      = *std::max_element(gVerts.begin(), gVerts.end(), vertsOrder);
    const std::array<SimplexId, 2> lVerts{
      discreteGradient_.getCellLowerVertex(src, triangulation),
      discreteGradient_.getCellLowerVertex(dst, triangulation)};
    const auto sepFuncMin
      = *std::min_element(lVerts.begin(), lVerts.end(), vertsOrder);
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[sepId] = sepFuncMax;
    if(separatrixFunctionMinima != nullptr)
      (*separatrixFunctionMinima)[sepId] = sepFuncMin;

    // get boundary condition
    const auto onBoundary
      = static_cast<char>(discreteGradient_.isBoundary(src, triangulation))
        + static_cast<char>(discreteGradient_.isBoundary(dst, triangulation));

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      std::array<float, 3> pt{};
      triangulation.getCellIncenter(cell.id_, cell.dim_, pt.data());

      // index of current point in point data arrays
      const auto k = geomPointsBegId[i] + j;

      points[3 * k + 0] = pt[0];
      points[3 * k + 1] = pt[1];
      points[3 * k + 2] = pt[2];

      if(outputSeparatrices1_points_smoothingMask_ != nullptr)
        (*outputSeparatrices1_points_smoothingMask_)[k]
          = (j == 0 || j == sepGeom.size() - 1) ? 0 : 1;
      if(outputSeparatrices1_points_cellDimensions_ != nullptr)
        (*outputSeparatrices1_points_cellDimensions_)[k] = cell.dim_;
      if(outputSeparatrices1_points_cellIds_ != nullptr)
        (*outputSeparatrices1_points_cellIds_)[k] = cell.id_;

      // skip filling cell data for first geometry point
      if(j == 0)
        continue;

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j - 1;

      cellsConn[2 * l + 0] = k - 1;
      cellsConn[2 * l + 1] = k;

      if(outputSeparatrices1_cells_sourceIds_ != nullptr)
        (*outputSeparatrices1_cells_sourceIds_)[l] = src.id_;
      if(outputSeparatrices1_cells_destinationIds_ != nullptr)
        (*outputSeparatrices1_cells_destinationIds_)[l] = dst.id_;
      if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
        (*outputSeparatrices1_cells_separatrixIds_)[l] = sepId;
      if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices1_cells_separatrixTypes_)[l] = sepType;
      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[l] = onBoundary;
    }
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = ncells;

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getAscendingSeparatrices1(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  const auto dim{triangulation.getDimensionality()};

  // Triangulation method pointers for 3D
  auto getFaceStarNumber = &triangulationType::getTriangleStarNumber;
  auto getFaceStar = &triangulationType::getTriangleStar;
  if(dim == 2) {
    // Triangulation method pointers for 2D
    getFaceStarNumber = &triangulationType::getEdgeStarNumber;
    getFaceStar = &triangulationType::getEdgeStar;
  }

  std::vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == dim - 1)
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
    const auto starNumber{(triangulation.*getFaceStarNumber)(saddle.id_)};
    for(SimplexId j = 0; j < starNumber; ++j) {

      SimplexId sId{};
      (triangulation.*getFaceStar)(saddle.id_, j, sId);

      std::vector<Cell> vpath{saddle};
      discreteGradient_.getAscendingPath(Cell(dim, sId), vpath, triangulation);

      const Cell &lastCell = vpath.back();
      if(lastCell.dim_ == dim and discreteGradient_.isCellCritical(lastCell)) {
        const SimplexId separatrixIndex = 4 * i + j;
        separatricesGeometry[separatrixIndex] = std::move(vpath);
        separatrices[separatrixIndex]
          = Separatrix(true, saddle, lastCell, false, separatrixIndex);
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setAscendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  std::vector<SimplexId> &maxSeeds,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation) const {

  const SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  const SimplexId numberOfCells = triangulation.getNumberOfCells();
  std::vector<SimplexId> morseSmaleManifoldOnCells(numberOfCells, -1);
  const int cellDim = triangulation.getDimensionality();

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
  auto getCellFace = &triangulationType::getCellTriangle;
  auto getFaceStarNumber = &triangulationType::getTriangleStarNumber;
  auto getFaceStar = &triangulationType::getTriangleStar;
  if(cellDim == 2) {
    // Triangulation method pointers for 2D
    getCellFace = &triangulationType::getCellEdge;
    getFaceStarNumber = &triangulationType::getEdgeStarNumber;
    getFaceStar = &triangulationType::getEdgeStar;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(SimplexId i = 0; i < numberOfSeeds; ++i) {
    std::queue<SimplexId> bfs;

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
          (triangulation.*getCellFace)(cofacetId, j, facetId);

          SimplexId starNumber = (triangulation.*getFaceStarNumber)(facetId);
          for(SimplexId k = 0; k < starNumber; ++k) {
            SimplexId neighborId = -1;
            (triangulation.*getFaceStar)(facetId, k, neighborId);

            if(neighborId == cofacetId) {
              continue;
            }

            const SimplexId pairedCellId = discreteGradient_.getPairedCell(
              Cell(cellDim, neighborId), triangulation, true);

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
    triangulation.getVertexStar(i, 0, starId);
    morseSmaleManifold[i] = morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setDescendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMinima,
  const triangulationType &triangulation) const {

  const SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  // get the seeds : minima
  std::vector<SimplexId> seeds;
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
    std::queue<SimplexId> bfs;

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
          = triangulation.getVertexEdgeNumber(vertexId);
        for(SimplexId j = 0; j < edgeNumber; ++j) {
          SimplexId edgeId;
          triangulation.getVertexEdge(vertexId, j, edgeId);

          for(int k = 0; k < 2; ++k) {
            SimplexId neighborId;
            triangulation.getEdgeVertex(edgeId, k, neighborId);

            if(neighborId == vertexId) {
              continue;
            }

            const SimplexId pairedCellId = discreteGradient_.getPairedCell(
              Cell(0, neighborId), triangulation);

            if(pairedCellId == edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId ttkNotUsed(numberOfMinima),
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  const size_t nVerts = triangulation.getNumberOfVertices();

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
  PSORT(this->threadNumber_)(sparseRegionIds.begin(), sparseRegionIds.end());
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

// ----------------------- //
//       2D Methods        //
// ----------------------- //

template <typename triangulationType>
int ttk::MorseSmaleComplex::execute2D(const triangulationType &triangulation) {
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

// ----------------------- //
//       3D Methods        //
// ----------------------- //

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplex::execute3D(const triangulationType &triangulation) {
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

  std::vector<std::vector<Separatrix>> separatrices1{};
  std::vector<std::vector<std::vector<dcg::Cell>>> separatricesGeometry1;

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getDescendingSeparatrices1(criticalPoints, separatrices1.back(),
                               separatricesGeometry1.back(), triangulation);

    this->printMsg("Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getAscendingSeparatrices1(criticalPoints, separatrices1.back(),
                              separatricesGeometry1.back(), triangulation);

    this->printMsg("Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  // saddle-connectors
  if(ComputeSaddleConnectors) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getSaddleConnectors(criticalPoints, separatrices1.back(),
                        separatricesGeometry1.back(), triangulation);

    this->printMsg("Saddle connectors computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  if(ComputeDescendingSeparatrices1 || ComputeAscendingSeparatrices1
     || ComputeSaddleConnectors) {
    Timer tmp{};

    flattenSeparatricesVectors(separatrices1, separatricesGeometry1);
    setSeparatrices1(separatrices1[0], separatricesGeometry1[0], triangulation);

    this->printMsg(
      "1-separatrices set", 1.0, tmp.getElapsedTime(), this->threadNumber_);
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices,
                               separatricesGeometry, separatricesSaddles,
                               triangulation);
    setDescendingSeparatrices2(
      separatrices, separatricesGeometry, separatricesSaddles, triangulation);

    this->printMsg("Descending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices,
                              separatricesGeometry, separatricesSaddles,
                              triangulation);
    setAscendingSeparatrices2(
      separatrices, separatricesGeometry, separatricesSaddles, triangulation);

    this->printMsg("Ascending 2-separatrices computed", 1.0,
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
    std::vector<size_t> nCriticalPointsByDim;
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
int ttk::MorseSmaleComplex::getSaddleConnectors(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  const auto nTriangles = triangulation.getNumberOfTriangles();
  // visited triangles (one vector per thread)
  std::vector<bool> isVisited(nTriangles, false);
  std::vector<SimplexId> visitedTriangles{};

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
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(isVisited, visitedTriangles)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    const auto &s2{saddles2[i]};

    std::set<SimplexId> saddles1{};
    VisitedMask mask{isVisited, visitedTriangles};
    discreteGradient_.getDescendingWall(
      s2, mask, triangulation, nullptr, &saddles1);

    for(const auto saddle1Id : saddles1) {
      const Cell s1{1, saddle1Id};

      Vpath vpath;
      const bool isMultiConnected
        = discreteGradient_.getAscendingPathThroughWall(
          s1, s2, isVisited, &vpath, triangulation);
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

template <typename triangulationType>
int ttk::MorseSmaleComplex::getAscendingSeparatrices2(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  std::vector<SimplexId> saddleIndexes;
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

  const auto nEdges = triangulation.getNumberOfEdges();
  std::vector<bool> isVisited(nEdges, false);
  std::vector<SimplexId> visitedEdges{};

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(isVisited, visitedEdges)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle1 = criticalPoints[saddleIndex];

    std::vector<Cell> wall;
    VisitedMask mask{isVisited, visitedEdges};
    discreteGradient_.getAscendingWall(
      saddle1, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle1, emptyCell, false, i);
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getDescendingSeparatrices2(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  std::vector<SimplexId> saddleIndexes;
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

  const auto nTriangles = triangulation.getNumberOfTriangles();
  std::vector<bool> isVisited(nTriangles, false);
  std::vector<SimplexId> visitedTriangles{};

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(isVisited, visitedTriangles)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle2 = criticalPoints[saddleIndex];

    std::vector<Cell> wall;
    VisitedMask mask{isVisited, visitedTriangles};
    discreteGradient_.getDescendingWall(
      saddle2, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle2, emptyCell, false, i);
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getDualPolygon(
  const SimplexId edgeId,
  SimplexId *const polygon,
  const size_t polSize,
  const triangulationType &triangulation) const {

  for(size_t i = 0; i < polSize; ++i) {
    SimplexId starId;
    triangulation.getEdgeStar(edgeId, i, starId);
    polygon[i] = starId;
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::sortDualPolygonVertices(
  SimplexId *const polygon,
  const size_t polSize,
  const triangulationType &triangulation) const {

  for(size_t i = 1; i < polSize; ++i) {

    // find polygon[i - 1] neighboring tetra in polygon[i..]
    bool isFound = false;
    size_t j = i;
    for(; j < polSize; ++j) {
      // check if current is the neighbor
      for(SimplexId k = 0;
          k < triangulation.getCellNeighborNumber(polygon[i - 1]); ++k) {
        SimplexId neighborId{};
        triangulation.getCellNeighbor(polygon[i - 1], k, neighborId);
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

template <typename triangulationType>
int ttk::MorseSmaleComplex::setAscendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    this->printErr("2-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices2_cells_connectivity_ == nullptr
     || outputSeparatrices2_cells_offsets_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  const auto offsets = inputOffsets_;
  auto separatrixFunctionMaxima = outputS2_cells_separatrixFunctionMaximaId_;
  auto separatrixFunctionMinima = outputS2_cells_separatrixFunctionMinimaId_;

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = (outputSeparatrices2_cells_separatrixIds_ != nullptr
       && !outputSeparatrices2_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
                            outputSeparatrices2_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices2_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices2_numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outputSeparatrices2_cells_connectivity_->size()};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      ncells += separatricesGeometry[geomId].size();
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  // store the separatrices info (one per separatrix)
  std::vector<SimplexId> sepSourceIds(validGeomIds.size());
  std::vector<SimplexId> sepIds(validGeomIds.size());
  std::vector<char> sepOnBoundary(validGeomIds.size());
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(separatrixId + validGeomIds.size());
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(separatrixId + validGeomIds.size());
  // store the polygonal cells tetras SimplexId
  std::vector<SimplexId> polygonNTetras(ncells - noldcells);
  std::vector<SimplexId> polygonEdgeIds(ncells - noldcells);
  std::vector<SimplexId> polygonSepInfosIds(ncells - noldcells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto &sepSaddles = separatricesSaddles[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle1

    // compute separatrix function diff
    const auto sepFuncMin
      = discreteGradient_.getCellLowerVertex(src, triangulation);
    const auto maxId = *std::max_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
        return offsets[discreteGradient_.getCellGreaterVertex(
                 Cell{2, a}, triangulation)]
               < offsets[discreteGradient_.getCellGreaterVertex(
                 Cell{2, b}, triangulation)];
      });
    const auto sepFuncMax
      = discreteGradient_.getCellGreaterVertex(Cell{2, maxId}, triangulation);

    // get boundary condition
    const char onBoundary
      = std::count_if(sepSaddles.begin(), sepSaddles.end(),
                      [&triangulation](const SimplexId a) {
                        return triangulation.isTriangleOnBoundary(a);
                      })
        + triangulation.isEdgeOnBoundary(src.id_);

    sepIds[i] = sepId;
    sepSourceIds[i] = src.id_;
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[sepId] = sepFuncMax;
    if(separatrixFunctionMinima != nullptr)
      (*separatrixFunctionMinima)[sepId] = sepFuncMin;
    sepOnBoundary[i] = onBoundary;

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      // index of current cell in cell data arrays
      const auto k = geomCellsBegId[i] + j - noldcells;

      polygonNTetras[k] = triangulation.getEdgeStarNumber(cell.id_);

      if(polygonNTetras[k] > 2) {
        polygonEdgeIds[k] = cell.id_;
        polygonSepInfosIds[k] = i;
      }
    }
  }

  // indices of valid polygon tetras
  std::vector<SimplexId> validTetraIds{};
  validTetraIds.reserve(polygonNTetras.size());

  for(size_t i = 0; i < polygonNTetras.size(); ++i) {
    if(polygonNTetras[i] > 2) {
      validTetraIds.emplace_back(i);
    }
  }

  // count number of valid new cells and new points
  size_t nnewpoints{};
  std::vector<SimplexId> pointsPerCell(validTetraIds.size() + 1);
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    nnewpoints += polygonNTetras[validTetraIds[i]];
    pointsPerCell[i + 1] = nnewpoints;
  }

  // resize connectivity array
  outputSeparatrices2_cells_connectivity_->resize(firstCellId + nnewpoints);
  auto cellsConn = &outputSeparatrices2_cells_connectivity_->at(firstCellId);
  // copy of cell connectivity array (for removing duplicates vertices)
  std::vector<SimplexId> cellVertsIds(nnewpoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    const auto k = validTetraIds[i];

    // get tetras in edge star
    getDualPolygon(polygonEdgeIds[k], &cellVertsIds[pointsPerCell[i]],
                   polygonNTetras[k], triangulation);
    // sort tetras (in-place)
    sortDualPolygonVertices(
      &cellVertsIds[pointsPerCell[i]], polygonNTetras[k], triangulation);

    for(SimplexId j = 0; j < polygonNTetras[k]; ++j) {
      cellsConn[pointsPerCell[i] + j] = cellVertsIds[pointsPerCell[i] + j];
    }
  }

  PSORT(this->threadNumber_)(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<SimplexId> vertId2PointsId(triangulation.getNumberOfCells());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  ncells = noldcells + validTetraIds.size();

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);
  outputSeparatrices2_cells_offsets_->resize(ncells + 1);
  outputSeparatrices2_cells_offsets_->at(0) = 0;
  auto cellsOff = &outputSeparatrices2_cells_offsets_->at(noldcells);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
    outputSeparatrices2_cells_isOnBoundary_->resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    triangulation.getTetraIncenter(cellVertsIds[i], &points[3 * i]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    const auto m = validTetraIds[i];
    const auto k = pointsPerCell[i];
    for(SimplexId j = 0; j < polygonNTetras[m]; ++j) {
      cellsConn[k + j] = vertId2PointsId[cellsConn[k + j]];
    }
    const auto l = i + noldcells;
    const auto n = polygonSepInfosIds[m];

    if(outputSeparatrices2_cells_sourceIds_ != nullptr)
      (*outputSeparatrices2_cells_sourceIds_)[l] = sepSourceIds[n];
    if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
      (*outputSeparatrices2_cells_separatrixIds_)[l] = sepIds[n];
    if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
      (*outputSeparatrices2_cells_separatrixTypes_)[l] = 1;
    if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices2_cells_isOnBoundary_)[l] = sepOnBoundary[n];
  }

  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    // fill offsets sequentially (due to iteration dependencies)
    cellsOff[i + 1] = cellsOff[i] + polygonNTetras[validTetraIds[i]];
  }

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setDescendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    this->printErr("2-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices2_cells_connectivity_ == nullptr
     || outputSeparatrices2_cells_offsets_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  const auto offsets = inputOffsets_;
  auto separatrixFunctionMaxima = outputS2_cells_separatrixFunctionMaximaId_;
  auto separatrixFunctionMinima = outputS2_cells_separatrixFunctionMinimaId_;

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = (outputSeparatrices2_cells_separatrixIds_ != nullptr
       && !outputSeparatrices2_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
                            outputSeparatrices2_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices2_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices2_numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outputSeparatrices2_cells_connectivity_->size()};

  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      ncells += separatricesGeometry[geomId].size();
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  // resize arrays
  outputSeparatrices2_cells_offsets_->resize(ncells + 1);
  outputSeparatrices2_cells_offsets_->at(0) = 0;
  outputSeparatrices2_cells_connectivity_->resize(
    firstCellId + 3 * (ncells - noldcells)); // triangles cells
  auto cellsOff = &outputSeparatrices2_cells_offsets_->at(noldcells);
  auto cellsConn = &outputSeparatrices2_cells_connectivity_->at(firstCellId);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(separatrixId + validGeomIds.size());
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(separatrixId + validGeomIds.size());
  if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
    outputSeparatrices2_cells_isOnBoundary_->resize(ncells);

  // store the cells/triangles vertices vertexId
  std::vector<SimplexId> cellVertsIds(3 * (ncells - noldcells));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto &sepSaddles = separatricesSaddles[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle2
    const char sepType = 2;

    // compute separatrix function diff
    const auto sepFuncMax
      = discreteGradient_.getCellGreaterVertex(src, triangulation);
    const auto minId = *std::min_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
        return offsets[discreteGradient_.getCellLowerVertex(
                 Cell{1, a}, triangulation)]
               < offsets[discreteGradient_.getCellLowerVertex(
                 Cell{1, b}, triangulation)];
      });
    const auto sepFuncMin
      = discreteGradient_.getCellLowerVertex(Cell{1, minId}, triangulation);
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[sepId] = sepFuncMax;
    if(separatrixFunctionMinima != nullptr)
      (*separatrixFunctionMinima)[sepId] = sepFuncMin;

    // get boundary condition
    const char onBoundary
      = std::count_if(sepSaddles.begin(), sepSaddles.end(),
                      [&triangulation](const SimplexId a) {
                        return triangulation.isEdgeOnBoundary(a);
                      })
        + triangulation.isTriangleOnBoundary(src.id_);

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];

      // first store the SimplexId of the cell/triangle vertices
      SimplexId v0{}, v1{}, v2{};
      triangulation.getTriangleVertex(cell.id_, 0, v0);
      triangulation.getTriangleVertex(cell.id_, 1, v1);
      triangulation.getTriangleVertex(cell.id_, 2, v2);

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j;
      // index of current cell among all new cells
      const auto m = l - noldcells;

      cellsConn[3 * m + 0] = v0;
      cellsConn[3 * m + 1] = v1;
      cellsConn[3 * m + 2] = v2;
      cellVertsIds[3 * m + 0] = v0;
      cellVertsIds[3 * m + 1] = v1;
      cellVertsIds[3 * m + 2] = v2;

      if(outputSeparatrices2_cells_sourceIds_ != nullptr)
        (*outputSeparatrices2_cells_sourceIds_)[l] = src.id_;
      if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
        (*outputSeparatrices2_cells_separatrixIds_)[l] = sepId;
      if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices2_cells_separatrixTypes_)[l] = sepType;
      if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices2_cells_isOnBoundary_)[l] = onBoundary;
    }
  }

  // reduce the cell vertices ids
  // (cells are triangles sharing two vertices)
  PSORT(this->threadNumber_)(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<size_t> vertId2PointsId(triangulation.getNumberOfVertices());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    triangulation.getVertexPoint(
      cellVertsIds[i], points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

  const auto lastOffset = noldcells == 0 ? 0 : cellsOff[-1];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < ncells - noldcells; ++i) {
    cellsOff[i] = 3 * i + lastOffset;
    cellsConn[3 * i + 0] = vertId2PointsId[cellsConn[3 * i + 0]];
    cellsConn[3 * i + 1] = vertId2PointsId[cellsConn[3 * i + 1]];
    cellsConn[3 * i + 2] = vertId2PointsId[cellsConn[3 * i + 2]];
  }

  cellsOff[ncells - noldcells] = cellsOff[ncells - noldcells - 1] + 3;

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}
