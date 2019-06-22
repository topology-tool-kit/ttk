/// \ingroup baseCode
/// \class ttk::DiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2016.
///
/// \brief TTK %discreteGradient processing package.
///
/// %DiscreteGradient is a TTK processing package that handles discrete gradient
/// (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation

#ifndef _DISCRETEGRADIENT_H
#define _DISCRETEGRADIENT_H

// base code includes
#include <FTMTree.h>
#include <Geometry.h>
#include <ScalarFieldCriticalPoints.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include <algorithm>
#include <set>

namespace ttk {
  namespace dcg {
    /**
     * Type of the identifiers of the 2-separatrices.
     * Must be the biggest integer type because it will provide more identifiers
     * for 2-separatrices.
     */
    using wallId_t = unsigned long long int;

    /**
     * Basic concept of cell, so it must be able to identify any cell of any
     * dimension.
     */
    struct Cell {
      explicit Cell() : dim_{-1}, id_{-1} {
      }

      explicit Cell(const int dim, const SimplexId id) : dim_{dim}, id_{id} {
      }

      explicit Cell(const Cell &cell) : dim_{cell.dim_}, id_{cell.id_} {
      }

      Cell &operator=(const Cell &cell) {
        dim_ = cell.dim_;
        id_ = cell.id_;
        return *this;
      }

      int dim_;
      SimplexId id_;
    };

    /**
     * Low-level structure storing a succession of cells. The orientation tells
     * whether the segment has been reversed or not.
     */
    struct Segment {
      explicit Segment() : orientation_{}, isValid_{} {
      }

      explicit Segment(const bool orientation,
                       const std::vector<Cell> &cells,
                       const bool isValid)
        : orientation_{orientation}, cells_{cells}, isValid_{isValid} {
      }

      explicit Segment(const bool orientation,
                       std::vector<Cell> &&cells,
                       const bool isValid)
        : orientation_{orientation}, cells_{cells}, isValid_{isValid} {
      }

      explicit Segment(const Segment &segment)
        : orientation_{segment.orientation_}, cells_{segment.cells_},
          isValid_{segment.isValid_} {
      }

      explicit Segment(Segment &&segment)
        : orientation_{segment.orientation_}, cells_{segment.cells_},
          isValid_{segment.isValid_} {
      }

      /**
       * Invalidate this segment so that it is ignored in the simplification
       * process, free memory for economy reasons.
       */
      int invalidate() {
        isValid_ = false;
        clear();

        return 0;
      }

      /**
       * Free internal cells' memory.
       */
      int clear() {
        cells_.clear();

        return 0;
      }

      bool orientation_;
      std::vector<Cell> cells_;
      bool isValid_;
    };

    /**
     * Sequence of cells such that two consecutive cells differ in dimension by
     * one.
     */
    struct VPath {
      explicit VPath()
        : isValid_{}, source_{-1}, destination_{-1}, sourceSlot_{-1},
          destinationSlot_{-1}, persistence_{} {
      }

      explicit VPath(const bool isValid,
                     const SimplexId segmentId,
                     const SimplexId source,
                     const SimplexId destination,
                     const SimplexId sourceSlot,
                     const SimplexId destinationSlot,
                     const double persistence)
        : isValid_{isValid}, states_{1}, segments_{segmentId}, source_{source},
          destination_{destination}, sourceSlot_{sourceSlot},
          destinationSlot_{destinationSlot}, persistence_{persistence} {
      }

      explicit VPath(const bool isValid,
                     const std::vector<char> &states,
                     const std::vector<SimplexId> &segments,
                     const SimplexId source,
                     const SimplexId destination,
                     const SimplexId sourceSlot,
                     const SimplexId destinationSlot,
                     const double persistence)
        : isValid_{isValid}, states_{states}, segments_{segments},
          source_{source}, destination_{destination}, sourceSlot_{sourceSlot},
          destinationSlot_{destinationSlot}, persistence_{persistence} {
      }

      explicit VPath(const bool isValid,
                     std::vector<char> &&states,
                     std::vector<SimplexId> &&segments,
                     const SimplexId source,
                     const SimplexId destination,
                     const SimplexId sourceSlot,
                     const SimplexId destinationSlot,
                     const double persistence)
        : isValid_{isValid}, states_{states}, segments_{segments},
          source_{source}, destination_{destination}, sourceSlot_{sourceSlot},
          destinationSlot_{destinationSlot}, persistence_{persistence} {
      }

      explicit VPath(const VPath &vpath)
        : isValid_{vpath.isValid_}, states_{vpath.states_},
          segments_{vpath.segments_}, source_{vpath.source_},
          destination_{vpath.destination_}, sourceSlot_{vpath.sourceSlot_},
          destinationSlot_{vpath.destinationSlot_}, persistence_{
                                                      vpath.persistence_} {
      }

      explicit VPath(VPath &&vpath)
        : isValid_{vpath.isValid_}, states_{vpath.states_},
          segments_{vpath.segments_}, source_{vpath.source_},
          destination_{vpath.destination_}, sourceSlot_{vpath.sourceSlot_},
          destinationSlot_{vpath.destinationSlot_}, persistence_{
                                                      vpath.persistence_} {
      }

      /**
       * Invalidate this vpath so that it is ignored in the simplification
       * process, free memory for economy reasons.
       */
      int invalidate() {
        isValid_ = false;
        source_ = -1;
        destination_ = -1;
        persistence_ = -1;
        clear();

        return 0;
      }

      /**
       * Free internal segments' memory.
       */
      int clear() {
        states_.clear();
        segments_.clear();

        return 0;
      }

      bool isValid_;
      std::vector<char> states_;
      std::vector<SimplexId> segments_;
      SimplexId source_;
      SimplexId destination_;
      SimplexId sourceSlot_;
      SimplexId destinationSlot_;
      double persistence_;
    };

    /**
     * Limit point of integral lines in the gradient.
     */
    struct CriticalPoint {
      explicit CriticalPoint() : numberOfSlots_{} {
      }

      explicit CriticalPoint(const Cell &cell) : cell_{cell}, numberOfSlots_{} {
      }

      explicit CriticalPoint(const Cell &cell,
                             const std::vector<SimplexId> &vpaths)
        : cell_{cell}, vpaths_{vpaths}, numberOfSlots_{} {
      }

      explicit CriticalPoint(const Cell &cell, std::vector<SimplexId> &&vpaths)
        : cell_{cell}, vpaths_{vpaths}, numberOfSlots_{} {
      }

      explicit CriticalPoint(const CriticalPoint &criticalPoint)
        : cell_{criticalPoint.cell_}, vpaths_{criticalPoint.vpaths_},
          numberOfSlots_{criticalPoint.numberOfSlots_} {
      }

      explicit CriticalPoint(CriticalPoint &&criticalPoint)
        : cell_{criticalPoint.cell_}, vpaths_{criticalPoint.vpaths_},
          numberOfSlots_{criticalPoint.numberOfSlots_} {
      }

      /**
       * Increase the connectivity of the critical point with openmp
       * acceleration if enabled.
       */
      SimplexId omp_addSlot() {
        SimplexId numberOfSlots = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
        numberOfSlots = (numberOfSlots_++);

        return numberOfSlots;
      }

      /**
       * Increase the connectivity of the critical point.
       */
      SimplexId addSlot() {
        return (numberOfSlots_++);
      }

      /**
       * Free the connectivity of the critical point from the memory.
       */
      int clear() {
        vpaths_.clear();

        return 0;
      }

      Cell cell_;
      std::vector<SimplexId> vpaths_;
      SimplexId numberOfSlots_;
    };

    /**
     * Comparator of VPaths, first compare persistence values then vpaths'
     * identifiers.
     */
    template <typename dataType>
    struct SaddleMaximumVPathComparator {
      bool operator()(const std::pair<dataType, SimplexId> &v1,
                      const std::pair<dataType, SimplexId> &v2) const {
        const dataType persistence1 = v1.first;
        const dataType persistence2 = v2.first;

        const SimplexId vpathId1 = v1.second;
        const SimplexId vpathId2 = v2.second;

        if(persistence1 != persistence2)
          return (persistence1 < persistence2);

        return (vpathId1 < vpathId2);
      };
    };

    /**
     * Comparator of saddle-connectors, first compare persistence values then
     * saddle identifiers and finally vpaths' identifiers.
     */
    template <typename dataType>
    struct SaddleSaddleVPathComparator {
      bool
        operator()(const std::tuple<dataType, SimplexId, SimplexId> &v1,
                   const std::tuple<dataType, SimplexId, SimplexId> &v2) const {
        const dataType persistence1 = std::get<0>(v1);
        const dataType persistence2 = std::get<0>(v2);

        const SimplexId vpathId1 = std::get<1>(v1);
        const SimplexId vpathId2 = std::get<1>(v2);

        const SimplexId saddleId1 = std::get<2>(v1);
        const SimplexId saddleId2 = std::get<2>(v2);

        if(persistence1 != persistence2)
          return (persistence1 < persistence2);

        if(saddleId1 != saddleId2)
          return (saddleId1 < saddleId2);

        return (vpathId1 < vpathId2);
      };
    };

    /**
     * Compute and manage a discrete gradient of a function on a triangulation.
     * TTK assumes that the input dataset is made of only one connected
     * component.
     */
    class DiscreteGradient : public Debug {

    public:
      explicit DiscreteGradient();

      ~DiscreteGradient();

      /**
       * Impose a threshold on the number of simplification passes.
       */
      int setIterationThreshold(const int iterationThreshold) {
        IterationThreshold = iterationThreshold;
        return 0;
      }

      /**
       * Enable/Disable gradient reversal of saddle to maximum VPaths.
       */
      int setReverseSaddleMaximumConnection(const bool state) {
        ReverseSaddleMaximumConnection = state;
        return 0;
      }

      /**
       * Enable/Disable gradient reversal of saddle/saddle VPaths.
       */
      int setReverseSaddleSaddleConnection(const bool state) {
        ReverseSaddleSaddleConnection = state;
        return 0;
      }

      /**
       * Enable/Disable collecting of persistence pairs during the
simplification.
       */
      int setCollectPersistencePairs(const bool state) {
        CollectPersistencePairs = state;
        return 0;
      }

      /**
       * Enable/Disable returning saddle-connectors as post-process.
       */
      int setReturnSaddleConnectors(const bool state) {
        ReturnSaddleConnectors = state;
        return 0;
      }

      /**
       * Set the persistence threshold value of the post-processing on the
saddle-connectors.
       */
      int setSaddleConnectorsPersistenceThreshold(const double threshold) {
        SaddleConnectorsPersistenceThreshold = threshold;
        return 0;
      }

      /**
       * Set the output data pointer to the container of the persistence pairs
       * (collecting of persistence pairs must be enabled).
       */
      int setOutputPersistencePairs(
        std::vector<std::tuple<Cell, Cell>> *const data) {
        outputPersistencePairs_ = data;
        return 0;
      }

      /**
       * Return the scalar value of the point in the cell which has the highest
function value.
       */
      template <typename dataType>
      dataType scalarMax(const Cell &cell, const dataType *const scalars) const;

      /**
       * Return the scalar value of the point in the cell which has the lowest
function value.
       */
      template <typename dataType>
      dataType scalarMin(const Cell &cell, const dataType *const scalars) const;

      /**
       * Compute the difference of function values of a pair of cells.
       */
      template <typename dataType>
      dataType getPersistence(const Cell &up,
                              const Cell &down,
                              const dataType *const scalars) const;

      /**
       * Body of AssignGradient algorithm from "Parallel Computation of 3D
Morse-Smale Complexes",
       * N. Shivashankar and V. Natarajan.
       * Compute the initial gradient field of the input scalar function for a
given dimension.
       */
      template <typename dataType, typename idType>
      int assignGradient(const int alphaDim,
                         const dataType *const scalars,
                         const idType *const offsets,
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
                         std::vector<std::vector<char>> &gradient) const;
#else
                         std::vector<std::vector<SimplexId>> &gradient) const;
#endif

      /**
       * Body of AssignGradient2 algorithm from "Parallel Computation of 3D
Morse-Smale Complexes",
       * N. Shivashankar and V. Natarajan.
       * Second pass of AssignGradient algorithm, minimize the number of
unpaired cells.
       */
      template <typename dataType, typename idType>
      int assignGradient2(const int alphaDim,
                          const dataType *const scalars,
                          const idType *const offsets,
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
                          std::vector<std::vector<char>> &gradient) const;
#else
                          std::vector<std::vector<SimplexId>> &gradient) const;
#endif

      /**
       * Brand new pass on the discrete gradient designed specifically for this
project,
       * the goal is to minimize the number of unpaired cells further (3D
triangulation only).
       */
      template <typename dataType, typename idType>
      int assignGradient3(const int alphaDim,
                          const dataType *const scalars,
                          const idType *const offsets,
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
                          std::vector<std::vector<char>> &gradient) const;
#else
                          std::vector<std::vector<SimplexId>> &gradient) const;
#endif

      /**
       * Compute the initial gradient field of the input scalar function on the
triangulation.
       */
      template <typename dataType, typename idType>
      int buildGradient();

      /**
       * Minimize the number of unpaired cells of any dimensions.
       * Assume that buildGradient() has been called before.
       */
      template <typename dataType, typename idType>
      int buildGradient2();

      /**
       * Minimize further the number of unpaired cells of any dimensions (3D
triangulation only).
       * Assume that buildGradient2() has been called before.
       */
      template <typename dataType, typename idType>
      int buildGradient3();

      /**
       * Get the list of maxima candidates for simplification.
       */
      template <typename dataType>
      int getRemovableMaxima(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableMaximum,
        std::vector<SimplexId> &pl2dmt_maximum);

      /**
       * Get the list of 1-saddles candidates for simplification.
       */
      template <typename dataType>
      int getRemovableSaddles1(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableSaddle,
        std::vector<SimplexId> &pl2dmt_saddle);

      /**
       * Get the list of 2-saddles candidates for simplification.
       */
      template <typename dataType>
      int getRemovableSaddles2(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableSaddle,
        std::vector<SimplexId> &pl2dmt_saddle);

      /**
       * Create initial Morse-Smale Complex structure and initialize the
(saddle,...,maximum)
       * vpaths to the simplification process.
       */
      template <typename dataType>
      int initializeSaddleMaximumConnections(
        std::vector<char> &isRemovableMaximum,
        std::vector<char> &isRemovableSaddle,
        const bool allowBruteForce,
        std::vector<Segment> &segments,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints) const;

      /**
       * Order the (saddle,...,maximum) vpaths by persistence value.
       */
      template <typename dataType>
      int orderSaddleMaximumConnections(
        const std::vector<VPath> &vpaths,
        std::set<std::pair<dataType, SimplexId>,
                 SaddleMaximumVPathComparator<dataType>> &S);

      /**
       * Compute simple algebra on the vpaths to minimize the number of gradient
paths reversal.
       * Two representations are available for the accumulation vector, a dense
and a sparse one (default).
       */
      template <typename dataType>
      int computeCoefficients(const bool isDense,
                              std::vector<char> &denseCoefficients,
                              std::vector<Segment> &segments,
                              const CriticalPoint &source,
                              VPath &newVPath,
                              const std::vector<VPath> &vpaths) const;

      /**
       * Core of the simplification process, tag the (saddle,...,maximum) vpaths
to be reversed.
       */
      template <typename dataType>
      int processSaddleMaximumConnections(
        const int iterationThreshold,
        const std::vector<char> &isPL,
        const bool allowBoundary,
        const bool allowBruteForce,
        std::set<std::pair<dataType, SimplexId>,
                 SaddleMaximumVPathComparator<dataType>> &S,
        std::vector<SimplexId> &pl2dmt_saddle,
        std::vector<SimplexId> &pl2dmt_maximum,
        std::vector<Segment> &segments,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints);

      /**
       * Actually reverse the so-tagged (saddle,...,maximum) vpaths to simplify
the discrete gradient.
       * The gradient is modified during this step.
       */
      template <typename dataType>
      int reverseSaddleMaximumConnections(const std::vector<Segment> &segments);

      /**
       * High-level function that manages the global simplification of
(saddle,...,maximum) vpaths.
       */
      template <typename dataType>
      int simplifySaddleMaximumConnections(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const std::vector<char> &isPL,
        const int iterationThreshold,
        const bool allowBoundary,
        const bool allowBruteForce);

      /**
       * Create initial Morse-Smale Complex structure and initialize the
(2-saddle,...,1-saddle) vpaths
       * to the simplification process.
       */
      template <typename dataType>
      int initializeSaddleSaddleConnections1(
        const std::vector<char> &isRemovableSaddle1,
        const std::vector<char> &isRemovableSaddle2,
        const bool allowBruteForce,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index) const;

      /**
       * Order the (2-saddle,...,1-saddle) vpaths by persistence value.
       */
      template <typename dataType>
      int orderSaddleSaddleConnections1(
        const std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::set<std::tuple<dataType, SimplexId, SimplexId>,
                 SaddleSaddleVPathComparator<dataType>> &S);

      /**
       * Core of the simplification process, modify the gradient and
       * reverse the selected (2-saddle,...,1-saddle) vpaths to simplify.
       */
      template <typename dataType>
      int processSaddleSaddleConnections1(
        const int iterationThreshold,
        const std::vector<char> &isPL,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors,
        std::set<std::tuple<dataType, SimplexId, SimplexId>,
                 SaddleSaddleVPathComparator<dataType>> &S,
        std::vector<SimplexId> &pl2dmt_saddle1,
        std::vector<SimplexId> &pl2dmt_saddle2,
        std::vector<char> &isRemovableSaddle1,
        std::vector<char> &isRemovableSaddle2,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index);

      /**
       * High-level function that manages the global simplification of
(2-saddle,...,1-saddle) vpaths.
       */
      template <typename dataType>
      int simplifySaddleSaddleConnections1(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const std::vector<char> &isPL,
        const int iterationThreshold,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors);

      /**
       * Create initial Morse-Smale Complex structure and initialize the
(1-saddle,...,2-saddle) vpaths
       * to the simplification process.
       */
      template <typename dataType>
      int initializeSaddleSaddleConnections2(
        const std::vector<char> &isRemovableSaddle1,
        const std::vector<char> &isRemovableSaddle2,
        const bool allowBruteForce,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index) const;

      /**
       * Order the (1-saddle,...,2-saddle) vpaths by persistence value.
       */
      template <typename dataType>
      int orderSaddleSaddleConnections2(
        const std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::set<std::tuple<dataType, SimplexId, SimplexId>,
                 SaddleSaddleVPathComparator<dataType>> &S);

      /**
       * Core of the simplification process, modify the gradient and
       * reverse the selected (1-saddle,...,2-saddle) vpaths to simplify.
       */
      template <typename dataType>
      int processSaddleSaddleConnections2(
        const int iterationThreshold,
        const std::vector<char> &isPL,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors,
        std::set<std::tuple<dataType, SimplexId, SimplexId>,
                 SaddleSaddleVPathComparator<dataType>> &S,
        std::vector<SimplexId> &pl2dmt_saddle1,
        std::vector<SimplexId> &pl2dmt_saddle2,
        std::vector<char> &isRemovableSaddle1,
        std::vector<char> &isRemovableSaddle2,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index);

      /**
       * High-level function that manages the global simplification of
(1-saddle,...,2-saddle) vpaths.
       */
      template <typename dataType>
      int simplifySaddleSaddleConnections2(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const std::vector<char> &isPL,
        const int iterationThreshold,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors);

      /**
       * Build the dense representation of the PL critical point list.
       */
      int getCriticalPointMap(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        std::vector<char> &isPL);

      /**
       * Process the saddle connectors by increasing value of persistence until
a given threshold is met.
       */
      template <typename dataType, typename idType>
      int filterSaddleConnectors(const bool allowBoundary);

      /**
       * Highest-level simplification function, manage all the simplification
steps
       * compliant to the critical points given by the user.
       */
      template <typename dataType, typename idType>
      int reverseGradient(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints);

      /**
       * Automatic detection of the PL critical points and simplification
according to them.
       */
      template <typename dataType, typename idType>
      int reverseGradient();

      /**
       * Set the input scalar function.
       */
      inline int setInputScalarField(const void *const data) {
        inputScalarField_ = data;
        return 0;
      }

      /**
       * Preprocess all the required connectivity requests on the triangulation.
       */
      inline int setupTriangulation(Triangulation *const data) {
        inputTriangulation_ = data;
        if(inputTriangulation_) {
          dimensionality_ = inputTriangulation_->getCellVertexNumber(0) - 1;
          numberOfVertices_ = inputTriangulation_->getNumberOfVertices();

          inputTriangulation_->preprocessBoundaryVertices();
          inputTriangulation_->preprocessBoundaryEdges();
          inputTriangulation_->preprocessVertexNeighbors();
          inputTriangulation_->preprocessVertexEdges();
          inputTriangulation_->preprocessVertexStars();
          inputTriangulation_->preprocessEdges();
          inputTriangulation_->preprocessEdgeStars();
          if(dimensionality_ == 2) {
            inputTriangulation_->preprocessCellEdges();
          } else if(dimensionality_ == 3) {
            inputTriangulation_->preprocessBoundaryTriangles();
            inputTriangulation_->preprocessVertexTriangles();
            inputTriangulation_->preprocessEdgeTriangles();
            inputTriangulation_->preprocessTriangles();
            inputTriangulation_->preprocessTriangleEdges();
            inputTriangulation_->preprocessTriangleStars();
            inputTriangulation_->preprocessCellTriangles();
          }
        }

        return 0;
      }

      /**
       * Set the input offset function.
       */
      inline int setInputOffsets(const void *const data) {
        inputOffsets_ = data;
        return 0;
      }

      /**
       * Set the output data pointer to the critical points.
       */
      inline int setOutputCriticalPoints(
        SimplexId *const criticalPoints_numberOfPoints,
        std::vector<float> *const criticalPoints_points,
        std::vector<char> *const criticalPoints_points_cellDimensons,
        std::vector<SimplexId> *const criticalPoints_points_cellIds,
        void *const criticalPoints_points_cellScalars,
        std::vector<char> *const criticalPoints_points_isOnBoundary,
        std::vector<SimplexId> *const criticalPoints_points_PLVertexIdentifiers,
        std::vector<SimplexId> *const criticalPoints_points_manifoldSize) {
        outputCriticalPoints_numberOfPoints_ = criticalPoints_numberOfPoints;
        outputCriticalPoints_points_ = criticalPoints_points;

        outputCriticalPoints_points_cellDimensions_
          = criticalPoints_points_cellDimensons;
        outputCriticalPoints_points_cellIds_ = criticalPoints_points_cellIds;

        outputCriticalPoints_points_cellScalars_
          = criticalPoints_points_cellScalars;

        outputCriticalPoints_points_isOnBoundary_
          = criticalPoints_points_isOnBoundary;

        outputCriticalPoints_points_PLVertexIdentifiers_
          = criticalPoints_points_PLVertexIdentifiers;

        outputCriticalPoints_points_manifoldSize_
          = criticalPoints_points_manifoldSize;
        return 0;
      }

      /**
       * Set the output data pointer to the gradient glyphs.
       */
      inline int setOutputGradientGlyphs(
        SimplexId *const gradientGlyphs_numberOfPoints,
        std::vector<float> *const gradientGlyphs_points,
        std::vector<char> *const gradientGlyphs_points_pairOrigins,
        SimplexId *const gradientGlyphs_numberOfCells,
        std::vector<SimplexId> *const gradientGlyphs_cells,
        std::vector<char> *const gradientGlyphs_cells_pairTypes) {
        outputGradientGlyphs_numberOfPoints_ = gradientGlyphs_numberOfPoints;
        outputGradientGlyphs_points_ = gradientGlyphs_points;

        outputGradientGlyphs_points_pairOrigins_
          = gradientGlyphs_points_pairOrigins;
        outputGradientGlyphs_numberOfCells_ = gradientGlyphs_numberOfCells;
        outputGradientGlyphs_cells_ = gradientGlyphs_cells;
        outputGradientGlyphs_cells_pairTypes_ = gradientGlyphs_cells_pairTypes;
        return 0;
      }

      /**
       * Get the dimensionality of the triangulation.
       */
      int getDimensionality() const;

      /**
       * Get the number of dimensions available for the cells in the
triangulation (equal to dimensionality+1).
       */
      int getNumberOfDimensions() const;

      /**
       * Get the number of cells of the given dimension.
       */
      SimplexId getNumberOfCells(const int dimension) const;

      /**
       * Return true if the given cell is a minimum regarding the discrete
gradient, false otherwise.
       */
      bool isMinimum(const Cell &cell) const;

      /**
       * Return true if the given cell is a 1-saddle regarding the discrete
gradient, false otherwise.
       */
      bool isSaddle1(const Cell &cell) const;

      /**
       * Return true if the given cell is a 2-saddle regarding the discrete
gradient, false otherwise.
       */
      bool isSaddle2(const Cell &cell) const;

      /**
       * Return true if the given cell is a maximum regarding the discrete
gradient, false otherwise.
       */
      bool isMaximum(const Cell &cell) const;

      /**
       * Return true if the given cell is a critical point regarding the
discrete gradient, false otherwise.
       */
      bool isCellCritical(const int cellDim, const SimplexId cellId) const;
      bool isCellCritical(const Cell &cell) const;

      /**
       * Return true if the given cell is at boundary, false otherwise.
       */
      bool isBoundary(const Cell &cell) const;

      /**
       * Return the identifier of the cell paired to the cell given by the user
in the gradient.
       */
      SimplexId getPairedCell(const Cell &cell, bool isReverse = false) const;

      /**
       * Get the output critical points as a STL vector of cells.
       */
      int getCriticalPoints(std::vector<Cell> &criticalPoints) const;

      /**
       * Return the VPath coming from the given cell.
       */
      int getAscendingPath(const Cell &cell,
                           std::vector<Cell> &vpath,
                           const bool enableCycleDetector = false) const;

      /**
       * Return the VPath terminating at the given cell.
       */
      int getDescendingPath(const Cell &cell, std::vector<Cell> &vpath) const;

      /**
       * Return the VPath terminating at the given 2-saddle restricted to the
2-separatrice of the 1-saddle.
       */
      int getDescendingPathThroughWall(const wallId_t wallId,
                                       const Cell &saddle2,
                                       const Cell &saddle1,
                                       const std::vector<wallId_t> &isVisited,
                                       std::vector<Cell> *const vpath,
                                       const bool enableCycleDetector
                                       = false) const;

      /**
       * Return the VPath coming from the given 1-saddle restricted to the
2-separatrice of the 2-saddle.
       */
      bool getAscendingPathThroughWall(const wallId_t wallId,
                                       const Cell &saddle1,
                                       const Cell &saddle2,
                                       const std::vector<wallId_t> &isVisited,
                                       std::vector<Cell> *const vpath,
                                       const bool enableCycleDetector
                                       = false) const;

      /**
       * Return the 2-separatrice terminating at the given 2-saddle.
       */
      int getDescendingWall(const wallId_t wallId,
                            const Cell &cell,
                            std::vector<wallId_t> &isVisited,
                            std::vector<Cell> *const wall = nullptr,
                            std::set<SimplexId> *const saddles = nullptr) const;

      /**
       * Return the 2-separatrice coming from the given 1-saddle.
       */
      int getAscendingWall(const wallId_t wallId,
                           const Cell &cell,
                           std::vector<wallId_t> &isVisited,
                           std::vector<Cell> *const wall = nullptr,
                           std::set<SimplexId> *const saddles = nullptr) const;

      /**
       * Reverse the given ascending VPath.
       */
      int reverseAscendingPath(const std::vector<Cell> &vpath);

      /**
       * Reverse the given ascending VPath restricted on a 2-separatrice.
       */
      int reverseAscendingPathOnWall(const std::vector<Cell> &vpath);

      /**
       * Reverse the given descending VPath restricted on a 2-separatrice.
       */
      int reverseDescendingPathOnWall(const std::vector<Cell> &vpath);

      /**
       * Compute the barycenter of the points of the given edge identifier.
       */
      int getEdgeIncenter(SimplexId edgeId, float incenter[3]) const;

      /**
       * Compute the incenter of the points of the given triangle identifier.
       */
      int getTriangleIncenter(SimplexId triangleId, float incenter[3]) const;

      /**
       * Compute the barycenter of the incenters of the triangles of the given
tetra identifier.
       */
      int getTetraIncenter(SimplexId tetraId, float incenter[3]) const;

      /**
       * Compute the geometric barycenter of a given cell.
       */
      int getCellIncenter(const Cell &cell, float incenter[3]) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output data pointers are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename dataType, typename idType>
      int setCriticalPoints(const std::vector<Cell> &criticalPoints) const;

      /**
       * Detect the critical points and build their geometric embedding.
       * The output data pointers are modified accordingly.
       */
      template <typename dataType, typename idType>
      int setCriticalPoints() const;

      /**
       * Build the geometric embedding of the given STL vector of cells and add
       * global information as scalar fields.
       * The output data pointers are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename dataType, typename idType>
      int setAugmentedCriticalPoints(const std::vector<Cell> &criticalPoints,
                                     std::vector<SimplexId> &maxSeeds,
                                     SimplexId *ascendingManifold,
                                     SimplexId *descendingManifold) const;

      /**
       * Build the glyphs representing the discrete gradient vector field.
       */
      int setGradientGlyphs() const;

    protected:
      int IterationThreshold;
      bool ReverseSaddleMaximumConnection;
      bool ReverseSaddleSaddleConnection;
      bool CollectPersistencePairs;
      bool ReturnSaddleConnectors;
      double SaddleConnectorsPersistenceThreshold;

      int dimensionality_;
      SimplexId numberOfVertices_;
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      std::vector<std::vector<std::vector<char>>> gradient_;
#else
      std::vector<std::vector<std::vector<SimplexId>>> gradient_;
#endif
      std::vector<SimplexId> dmtMax2PL_;
      std::vector<SimplexId> dmt1Saddle2PL_;
      std::vector<SimplexId> dmt2Saddle2PL_;

      const void *inputScalarField_;
      const void *inputOffsets_;
      Triangulation *inputTriangulation_;

      SimplexId *outputCriticalPoints_numberOfPoints_;
      std::vector<float> *outputCriticalPoints_points_;
      std::vector<char> *outputCriticalPoints_points_cellDimensions_;
      std::vector<SimplexId> *outputCriticalPoints_points_cellIds_;
      void *outputCriticalPoints_points_cellScalars_;
      std::vector<char> *outputCriticalPoints_points_isOnBoundary_;
      std::vector<SimplexId> *outputCriticalPoints_points_PLVertexIdentifiers_;
      std::vector<SimplexId> *outputCriticalPoints_points_manifoldSize_;

      SimplexId *outputGradientGlyphs_numberOfPoints_;
      std::vector<float> *outputGradientGlyphs_points_;
      std::vector<char> *outputGradientGlyphs_points_pairOrigins_;
      SimplexId *outputGradientGlyphs_numberOfCells_;
      std::vector<SimplexId> *outputGradientGlyphs_cells_;
      std::vector<char> *outputGradientGlyphs_cells_pairTypes_;

      std::vector<std::tuple<Cell, Cell>> *outputPersistencePairs_;
    };

#include <DiscreteGradient_Template.h>
  } // namespace dcg
} // namespace ttk

#endif // DISCRETEGRADIENT_H
