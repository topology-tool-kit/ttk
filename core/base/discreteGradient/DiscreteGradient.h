/// \ingroup baseCode
/// \class ttk::dcg::DiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Yizhe Wang <wangyizhe3518@gmail.com>
/// \date November 2016.
///
/// \brief TTK %discreteGradient processing package.
///
/// %DiscreteGradient is a TTK processing package that handles discrete gradient
/// (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <FTMTree.h>
#include <Geometry.h>
#include <Triangulation.h>
#include <VisitedMask.h>

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <set>
#include <utility>

namespace ttk {

  namespace dcg {
    /**
     * Basic concept of cell, so it must be able to identify any cell of any
     * dimension.
     */
    struct Cell {
      explicit Cell() = default;

      explicit Cell(const int dim, const SimplexId id) : dim_{dim}, id_{id} {
      }

      int dim_{-1};
      SimplexId id_{-1};
    };

    /**
     * @brief Extended Cell structure for processLowerStars
     */
    struct CellExt : Cell {
      explicit CellExt(const int dim, const SimplexId id) : Cell{dim, id} {
      }
      explicit CellExt(const int dim,
                       const SimplexId id,
                       const std::array<SimplexId, 3> &lowVerts,
                       const std::array<uint8_t, 3> &faces)
        : Cell{dim, id}, lowVerts_{lowVerts}, faces_{faces} {
      }

      // (order field value on) lower vertices in current lower star
      // (1 for edges, 2 for triangles, 3 for tetras)
      const std::array<SimplexId, 3> lowVerts_{};
      // indices of faces (cells of dimensions dim_ - 1) in lower star
      // structure, only applicable for triangles (2 edge faces)  and tetras (3
      // triangle faces)
      const std::array<uint8_t, 3> faces_{};
      // if cell has been paired with another in current lower star
      bool paired_{false};
    };

    /**
     * Low-level structure storing a succession of cells. The orientation tells
     * whether the segment has been reversed or not.
     */
    struct Segment {
      explicit Segment() = default;

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

      bool orientation_{};
      std::vector<Cell> cells_{};
      bool isValid_{};
    };

    /**
     * Sequence of cells such that two consecutive cells differ in dimension by
     * one.
     */
    struct VPath {
      explicit VPath() = default;

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

      bool isValid_{};
      std::vector<char> states_{};
      std::vector<SimplexId> segments_{};
      SimplexId source_{-1};
      SimplexId destination_{-1};
      SimplexId sourceSlot_{-1};
      SimplexId destinationSlot_{-1};
      double persistence_{};
    };

    /**
     * Limit point of integral lines in the gradient.
     */
    struct CriticalPoint {
      explicit CriticalPoint() = default;

      explicit CriticalPoint(const Cell &cell) : cell_{cell}, numberOfSlots_{} {
      }

      explicit CriticalPoint(const Cell &cell,
                             const std::vector<SimplexId> &vpaths)
        : cell_{cell}, vpaths_{vpaths}, numberOfSlots_{} {
      }

      explicit CriticalPoint(const Cell &cell, std::vector<SimplexId> &&vpaths)
        : cell_{cell}, vpaths_{vpaths}, numberOfSlots_{} {
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

      Cell cell_{};
      std::vector<SimplexId> vpaths_{};
      SimplexId numberOfSlots_{};
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

        if(persistence1 != persistence2) {
          return (persistence1 < persistence2);
        }

        if(saddleId1 != saddleId2) {
          return (saddleId1 < saddleId2);
        }

        return (vpathId1 < vpathId2);
      }
    };

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
    using gradIdType = char;
#else
    using gradIdType = SimplexId;
#endif

    /**
     * @brief Discrete gradient struct
     *
     * 0: paired edge id per vertex
     * 1: paired vertex id per edge
     * 2: paired triangle id per edge
     * 3: paired edge id per triangle
     * 4: paired tetra id per triangle
     * 5: paired triangle id per tetra
     * -1 if critical or paired to a cell of another dimension
     */
    using gradientType = std::array<std::vector<gradIdType>, 6>;

    /**
     * Compute and manage a discrete gradient of a function on a triangulation.
     * TTK assumes that the input dataset is made of only one connected
     * component.
     */
    class DiscreteGradient : virtual public Debug {

    public:
      DiscreteGradient() {
        this->setDebugMsgPrefix("DiscreteGradient");
      }

      /**
       * Impose a threshold on the number of simplification passes.
       */
      int setIterationThreshold(const int iterationThreshold) {
        IterationThreshold = iterationThreshold;
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
        std::vector<std::array<Cell, 2>> *const data) {
        outputPersistencePairs_ = data;
        return 0;
      }

      /**
       * Compute the initial gradient field of the input scalar function on the
triangulation.
       */
      template <typename triangulationType>
      int buildGradient(const triangulationType &triangulation);

      /**
       * Automatic detection of the PL critical points and simplification
according to them.
       */
      template <typename dataType, typename triangulationType>
      int reverseGradient(const triangulationType &triangulation,
                          const bool detectCriticalPoints = true);

      /**
       * Compute the (saddle1, saddle2) pairs not detected by the
       * contour tree.
       */
      template <typename dataType, typename triangulationType>
      void computeSaddleSaddlePersistencePairs(
        std::vector<std::tuple<SimplexId, SimplexId, dataType>>
          &pl_saddleSaddlePairs,
        const triangulationType &triangulation);

      /**
       * Set the input scalar function.
       */
      inline void setInputScalarField(const void *const data) {
        inputScalarField_ = data;
      }

      /**
       * Preprocess all the required connectivity requests on the triangulation.
       */
      inline void preconditionTriangulation(AbstractTriangulation *const data) {
        if(data != nullptr) {
          dimensionality_ = data->getCellVertexNumber(0) - 1;
          numberOfVertices_ = data->getNumberOfVertices();

          data->preconditionBoundaryVertices();
          data->preconditionVertexNeighbors();
          data->preconditionVertexEdges();
          data->preconditionVertexStars();
          data->preconditionEdges();
          data->preconditionEdgeStars();
          if(dimensionality_ >= 2) {
            data->preconditionBoundaryEdges();
          }
          if(dimensionality_ == 2) {
            data->preconditionCellEdges();
          } else if(dimensionality_ == 3) {
            data->preconditionBoundaryTriangles();
            data->preconditionVertexTriangles();
            data->preconditionEdgeTriangles();
            data->preconditionTriangles();
            data->preconditionTriangleEdges();
            data->preconditionTriangleStars();
            data->preconditionCellTriangles();
            // for filterSaddleConnectors
            contourTree_.preconditionTriangulation(data);
          }
          this->initMemory(*data);
        }
      }

      /**
       * Set the input offset function.
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
      template <typename triangulationType>
      SimplexId getNumberOfCells(const int dimension,
                                 const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is at boundary, false otherwise.
       */
      template <typename triangulationType>
      bool isBoundary(const Cell &cell,
                      const triangulationType &triangulation) const;

      /**
       * Return true if the given cell is a critical point regarding the
discrete gradient, false otherwise.
       */
      bool isCellCritical(const int cellDim, const SimplexId cellId) const;
      bool isCellCritical(const Cell &cell) const;

      /**
       * Return the identifier of the cell paired to the cell given by the user
in the gradient.
       */
      template <typename triangulationType>
      SimplexId getPairedCell(const Cell &cell,
                              const triangulationType &triangulation,
                              bool isReverse = false) const;

      /**
       * Return the VPath coming from the given cell.
       */
      template <typename triangulationType>
      int getAscendingPath(const Cell &cell,
                           std::vector<Cell> &vpath,
                           const triangulationType &triangulation,
                           const bool enableCycleDetector = false) const;

      /**
       * Return the VPath terminating at the given cell.
       */
      template <typename triangulationType>
      int getDescendingPath(const Cell &cell,
                            std::vector<Cell> &vpath,
                            const triangulationType &triangulation) const;

      /**
       * Return the VPath terminating at the given 2-saddle restricted to the
2-separatrice of the 1-saddle.
       */
      template <typename triangulationType>
      bool getDescendingPathThroughWall(const Cell &saddle2,
                                        const Cell &saddle1,
                                        const std::vector<bool> &isVisited,
                                        std::vector<Cell> *const vpath,
                                        const triangulationType &triangulation,
                                        const bool stopIfMultiConnected = false,
                                        const bool enableCycleDetector
                                        = false) const;

      /**
       * Return the VPath coming from the given 1-saddle restricted to the
2-separatrice of the 2-saddle.
       */
      template <typename triangulationType>
      bool getAscendingPathThroughWall(const Cell &saddle1,
                                       const Cell &saddle2,
                                       const std::vector<bool> &isVisited,
                                       std::vector<Cell> *const vpath,
                                       const triangulationType &triangulation,
                                       const bool stopIfMultiConnected = false,
                                       const bool enableCycleDetector
                                       = false) const;

      /**
       * Return the 2-separatrice terminating at the given 2-saddle.
       */
      template <typename triangulationType>
      int getDescendingWall(const Cell &cell,
                            VisitedMask &mask,
                            const triangulationType &triangulation,
                            std::vector<Cell> *const wall = nullptr,
                            std::set<SimplexId> *const saddles = nullptr) const;

      /**
       * Return the 2-separatrice coming from the given 1-saddle.
       */
      template <typename triangulationType>
      int getAscendingWall(const Cell &cell,
                           VisitedMask &mask,
                           const triangulationType &triangulation,
                           std::vector<Cell> *const wall = nullptr,
                           std::set<SimplexId> *const saddles = nullptr) const;

      /**
       * Get the vertex id of with the maximum scalar field value on
       * the given cell.
       */
      template <typename triangulationType>
      SimplexId
        getCellGreaterVertex(const Cell c,
                             const triangulationType &triangulation) const;

      /**
       * Get the vertex id of with the minimum scalar field value on
       * the given cell.
       */
      template <typename triangulationType>
      SimplexId
        getCellLowerVertex(const Cell c,
                           const triangulationType &triangulation) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output vectors are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename triangulationType>
      int setCriticalPoints(const std::vector<Cell> &criticalPoints,
                            std::vector<size_t> &nCriticalPointsByDim,
                            std::vector<std::array<float, 3>> &points,
                            std::vector<char> &cellDimensions,
                            std::vector<SimplexId> &cellIds,
                            std::vector<char> &isOnBoundary,
                            std::vector<SimplexId> &PLVertexIdentifiers,
                            const triangulationType &triangulation) const;

      /**
       * Detect the critical points and build their geometric embedding.
       * The output vectors are modified accordingly.
       */
      template <typename triangulationType>
      int setCriticalPoints(std::vector<std::array<float, 3>> &points,
                            std::vector<char> &cellDimensions,
                            std::vector<SimplexId> &cellIds,
                            std::vector<char> &isOnBoundary,
                            std::vector<SimplexId> &PLVertexIdentifiers,
                            const triangulationType &triangulation) const;

      /**
       * Get the output critical points as a STL vector of cells.
       */
      template <typename triangulationType>
      int getCriticalPoints(std::vector<Cell> &criticalPoints,
                            const triangulationType &triangulation) const;

      /**
       * Compute manifold size for critical extrema
       */
      int setManifoldSize(const size_t nCritPoints,
                          const std::vector<size_t> &nCriticalPointsByDim,
                          const SimplexId *const ascendingManifold,
                          const SimplexId *const descendingManifold,
                          std::vector<SimplexId> &manifoldSize) const;

      /**
       * Build the glyphs representing the discrete gradient vector field.
       */
      template <typename triangulationType>
      int setGradientGlyphs(std::vector<std::array<float, 3>> &points,
                            std::vector<char> &points_pairOrigins,
                            std::vector<char> &cells_pairTypes,
                            const triangulationType &triangulation) const;

    private:
      /**
       * Type alias for lower stars of a given cell
       */
      using lowerStarType = std::array<std::vector<CellExt>, 4>;

      /**
       * @brief Store the subcomplexes around vertex for which offset
       * at vertex is maximum
       *
       * @param[in] a Vertex Id
       *
       * @return Lower star as 4 sets of cells (0-cells, 1-cells, 2-cells and
       * 3-cells)
       */
      template <typename triangulationType>
      inline void lowerStar(lowerStarType &ls,
                            const SimplexId a,
                            const SimplexId *const offsets,
                            const triangulationType &triangulation) const;

      /**
       * @brief Return the number of unpaired faces of a given cell in
       * a lower star
       *
       * @param[in] c Input cell
       * @param[in] ls Input lower star
       *
       * @return Number of unpaired faces and a face id
       */
      std::pair<size_t, SimplexId>
        numUnpairedFaces(const CellExt &c, const lowerStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTriangle(const CellExt &c,
                                 const lowerStarType &ls) const;
      std::pair<size_t, SimplexId>
        numUnpairedFacesTetra(const CellExt &c, const lowerStarType &ls) const;

      /**
       * @brief Return the critical type corresponding to given
       * dimension
       */
      CriticalType criticalTypeFromCellDimension(const int dim) const;

      /**
       * @brief Pair cells into discrete gradient field
       *
       * @param[in] alpha Cell of lower dimension
       * @param[in] beta Cell of higher dimension
       */
      template <typename triangulationType>
      inline void pairCells(CellExt &alpha,
                            CellExt &beta,
                            const triangulationType &triangulation);

      /**
       * Implements the ProcessLowerStars algorithm from "Theory and
       * Algorithms for Constructing Discrete Morse Complexes from
       * Grayscale Digital Images", V. Robins, P. J. Wood,
       * A. P. Sheppard
       */
      template <typename triangulationType>
      int processLowerStars(const SimplexId *const offsets,
                            const triangulationType &triangulation);

      /**
       * Compute the difference of function values of a pair of cells.
       */
      template <typename dataType, typename triangulationType>
      dataType getPersistence(const Cell &up,
                              const Cell &down,
                              const dataType *const scalars,
                              const triangulationType &triangulation) const;

      /**
       * Get the list of maxima candidates for simplification.
       */
      template <typename dataType, typename triangulationType>
      int getRemovableMaxima(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableMaximum,
        std::vector<SimplexId> &pl2dmt_maximum,
        const triangulationType &triangulation);

      /**
       * Get the list of 1-saddles candidates for simplification.
       */
      template <typename dataType, typename triangulationType>
      int getRemovableSaddles1(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableSaddle,
        std::vector<SimplexId> &pl2dmt_saddle,
        const triangulationType &triangulation);

      /**
       * Get the list of 2-saddles candidates for simplification.
       */
      template <typename dataType, typename triangulationType>
      int getRemovableSaddles2(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const bool allowBoundary,
        std::vector<char> &isRemovableSaddle,
        std::vector<SimplexId> &pl2dmt_saddle,
        const triangulationType &triangulation);

      /**
       * Create initial Morse-Smale Complex structure and initialize the
(2-saddle,...,1-saddle) vpaths
       * to the simplification process.
       */
      template <typename dataType, typename triangulationType>
      int initializeSaddleSaddleConnections1(
        const std::vector<char> &isRemovableSaddle1,
        const std::vector<char> &isRemovableSaddle2,
        const bool allowBruteForce,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index,
        const triangulationType &triangulation) const;

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
      template <typename dataType, typename triangulationType>
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
        std::vector<SimplexId> &saddle2Index,
        const triangulationType &triangulation);

      /**
       * High-level function that manages the global simplification of
(2-saddle,...,1-saddle) vpaths.
       */
      template <typename dataType, typename triangulationType>
      int simplifySaddleSaddleConnections1(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const std::vector<char> &isPL,
        const int iterationThreshold,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors,
        const triangulationType &triangulation);

      /**
       * Create initial Morse-Smale Complex structure and initialize the
(1-saddle,...,2-saddle) vpaths
       * to the simplification process.
       */
      template <typename dataType, typename triangulationType>
      int initializeSaddleSaddleConnections2(
        const std::vector<char> &isRemovableSaddle1,
        const std::vector<char> &isRemovableSaddle2,
        const bool allowBruteForce,
        std::vector<VPath> &vpaths,
        std::vector<CriticalPoint> &criticalPoints,
        std::vector<SimplexId> &saddle1Index,
        std::vector<SimplexId> &saddle2Index,
        const triangulationType &triangulation) const;

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
      template <typename dataType, typename triangulationType>
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
        std::vector<SimplexId> &saddle2Index,
        const triangulationType &triangulation);

      /**
       * High-level function that manages the global simplification of
(1-saddle,...,2-saddle) vpaths.
       */
      template <typename dataType, typename triangulationType>
      int simplifySaddleSaddleConnections2(
        const std::vector<std::pair<SimplexId, char>> &criticalPoints,
        const std::vector<char> &isPL,
        const int iterationThreshold,
        const bool allowBoundary,
        const bool allowBruteForce,
        const bool returnSaddleConnectors,
        const triangulationType &triangulation);

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
      template <typename dataType, typename triangulationType>
      int filterSaddleConnectors(const bool allowBoundary,
                                 const triangulationType &triangulation);

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
       * Reverse the given ascending VPath.
       */
      template <typename triangulationType>
      int reverseAscendingPath(const std::vector<Cell> &vpath,
                               const triangulationType &triangulation);

      /**
       * Reverse the given descending VPath.
       */
      template <typename triangulationType>
      int reverseDescendingPath(const std::vector<Cell> &vpath,
                                const triangulationType &triangulation);

      /**
       * Reverse the given ascending VPath restricted on a 2-separatrice.
       */
      template <typename triangulationType>
      int reverseAscendingPathOnWall(const std::vector<Cell> &vpath,
                                     const triangulationType &triangulation);

      /**
       * Reverse the given descending VPath restricted on a 2-separatrice.
       */
      template <typename triangulationType>
      int reverseDescendingPathOnWall(const std::vector<Cell> &vpath,
                                      const triangulationType &triangulation);

      /**
       * @brief Initialize/Allocate discrete gradient memory
       */
      void initMemory(const AbstractTriangulation &triangulation);

    protected:
      ftm::FTMTree contourTree_{};

      int IterationThreshold{-1};
      bool CollectPersistencePairs{false};
      bool ReturnSaddleConnectors{false};
      double SaddleConnectorsPersistenceThreshold{0.0};

      int dimensionality_{-1};
      SimplexId numberOfVertices_{};
      gradientType gradient_{};
      std::vector<SimplexId> dmtMax2PL_{};
      std::vector<SimplexId> dmt1Saddle2PL_{};
      std::vector<SimplexId> dmt2Saddle2PL_{};
      std::vector<std::array<Cell, 2>> *outputPersistencePairs_{};

      const void *inputScalarField_{};
      const SimplexId *inputOffsets_{};
    };

  } // namespace dcg
} // namespace ttk

#include <DiscreteGradient_Template.h>
