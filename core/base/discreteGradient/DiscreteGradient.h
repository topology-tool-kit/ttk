/// \ingroup baseCode
/// \class ttk::DiscreteGradient
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
#include <Wrapper.h>

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <set>
#include <utility>

#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#endif

namespace ttk {
#if defined(_GLIBCXX_PARALLEL_FEATURES_H) && defined(TTK_ENABLE_OPENMP)
#define PSORT                               \
  omp_set_num_threads(this->threadNumber_); \
  __gnu_parallel::sort
#else
#define PSORT std::sort
#endif // _GLIBCXX_PARALLEL_FEATURES_H && TTK_ENABLE_OPENMP

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

      // if cell has been paired with another in current lower star
      bool paired_{false};
      // lower vertices in current lower star (1 for edges, 2 for triangles, 3
      // for tetras)
      const std::array<SimplexId, 3> lowVerts_{};
      // indices of faces (cells of dimensions dim_ - 1) in lower star
      // structure, only applicable for triangles (2 edge faces)  and tetras (3
      // triangle faces)
      const std::array<uint8_t, 3> faces_{};
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
      };
    };

    struct VisitedMask {
      std::vector<bool> &isVisited_;
      std::vector<SimplexId> &visitedIds_;

      VisitedMask(std::vector<bool> &isVisited,
                  std::vector<SimplexId> &visitedIds)
        : isVisited_{isVisited}, visitedIds_{visitedIds} {
      }
      ~VisitedMask() {
        // use RAII to clean & reset referenced vectors
        for(const auto id : this->visitedIds_) {
          this->isVisited_[id] = false;
        }
        // set size to 0 but keep allocated memory
        this->visitedIds_.clear();
      }
    };

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
    using gradientType = std::vector<std::vector<std::vector<char>>>;
#else
    using gradientType = std::vector<std::vector<std::vector<SimplexId>>>;
#endif

    /**
     * Compute and manage a discrete gradient of a function on a triangulation.
     * TTK assumes that the input dataset is made of only one connected
     * component.
     */
    class DiscreteGradient : public Debug {

    public:
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
                              const dataType *scalars) const;

    private:
      template <typename scalarType, typename offsetType>
      void sortVertices(const SimplexId vertexNumber,
                        std::vector<size_t> &vertsOrder,
                        const scalarType *const scalarField,
                        const offsetType *const offsetField) const {

        std::vector<SimplexId> sortedVertices(vertexNumber);
        vertsOrder.resize(vertexNumber);

        // fill with numbers from 0 to vertexNumber - 1
        std::iota(sortedVertices.begin(), sortedVertices.end(), 0);

        // sort vertices in ascending order following scalarfield / offsets
        PSORT(sortedVertices.begin(), sortedVertices.end(),
              [&](const SimplexId a, const SimplexId b) {
                return (scalarField[a] < scalarField[b])
                       || (scalarField[a] == scalarField[b]
                           && offsetField[a] < offsetField[b]);
              });

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t i = 0; i < vertsOrder.size(); ++i) {
          vertsOrder[sortedVertices[i]] = i;
        }
      }

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
      inline lowerStarType lowerStar(const SimplexId a) const;

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
      inline void pairCells(CellExt &alpha, CellExt &beta);

      /**
       * Implements the ProcessLowerStars algorithm from "Theory and
       * Algorithms for Constructing Discrete Morse Complexes from
       * Grayscale Digital Images", V. Robins, P. J. Wood,
       * A. P. Sheppard
       */
      int processLowerStars();

    public:
      /**
       * Compute the initial gradient field of the input scalar function on the
triangulation.
       */
      template <typename dataType, typename idType>
      int buildGradient();

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
       * Automatic detection of the PL critical points and simplification
according to them.
       */
      template <typename dataType, typename idType>
      int reverseGradient(const bool detectCriticalPoints = true);

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
        if(inputTriangulation_ != nullptr) {
          dimensionality_ = inputTriangulation_->getCellVertexNumber(0) - 1;
          numberOfVertices_ = inputTriangulation_->getNumberOfVertices();

          inputTriangulation_->preconditionBoundaryVertices();
          inputTriangulation_->preconditionBoundaryEdges();
          inputTriangulation_->preconditionVertexNeighbors();
          inputTriangulation_->preconditionVertexEdges();
          inputTriangulation_->preconditionVertexStars();
          inputTriangulation_->preconditionEdges();
          inputTriangulation_->preconditionEdgeStars();
          if(dimensionality_ == 2) {
            inputTriangulation_->preconditionCellEdges();
          } else if(dimensionality_ == 3) {
            inputTriangulation_->preconditionBoundaryTriangles();
            inputTriangulation_->preconditionVertexTriangles();
            inputTriangulation_->preconditionEdgeTriangles();
            inputTriangulation_->preconditionTriangles();
            inputTriangulation_->preconditionTriangleEdges();
            inputTriangulation_->preconditionTriangleStars();
            inputTriangulation_->preconditionCellTriangles();
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
      bool getDescendingPathThroughWall(const Cell &saddle2,
                                        const Cell &saddle1,
                                        const std::vector<bool> &isVisited,
                                        std::vector<Cell> *const vpath,
                                        const bool stopIfMultiConnected = false,
                                        const bool enableCycleDetector
                                        = false) const;

      /**
       * Return the VPath coming from the given 1-saddle restricted to the
2-separatrice of the 2-saddle.
       */
      bool getAscendingPathThroughWall(const Cell &saddle1,
                                       const Cell &saddle2,
                                       const std::vector<bool> &isVisited,
                                       std::vector<Cell> *const vpath,
                                       const bool stopIfMultiConnected = false,
                                       const bool enableCycleDetector
                                       = false) const;

      /**
       * Return the 2-separatrice terminating at the given 2-saddle.
       */
      int getDescendingWall(const Cell &cell,
                            VisitedMask &mask,
                            std::vector<Cell> *const wall = nullptr,
                            std::set<SimplexId> *const saddles = nullptr) const;

      /**
       * Return the 2-separatrice coming from the given 1-saddle.
       */
      int getAscendingWall(const Cell &cell,
                           VisitedMask &mask,
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
       * Get the vertex id of with the maximum scalar field value on
       * the given cell. Compare offsets if scalar field is constant.
       */
      SimplexId getCellGreaterVertex(const Cell c) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output data pointers are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename dataType>
      int setCriticalPoints(const std::vector<Cell> &criticalPoints,
                            std::vector<size_t> &nCriticalPointsByDim) const;

      /**
       * Detect the critical points and build their geometric embedding.
       * The output data pointers are modified accordingly.
       */
      template <typename dataType>
      int setCriticalPoints() const;

      /**
       * Compute manifold size for critical extrema
       */
      int setManifoldSize(const std::vector<Cell> &criticalPoints,
                          const std::vector<size_t> &nCriticalPointsByDim,
                          const std::vector<SimplexId> &maxSeeds,
                          const SimplexId *const ascendingManifold,
                          const SimplexId *const descendingManifold) const;

      /**
       * Build the glyphs representing the discrete gradient vector field.
       */
      int setGradientGlyphs(SimplexId &numberOfPoints,
                            std::vector<float> &points,
                            std::vector<char> &points_pairOrigins,
                            SimplexId &numberOfCells,
                            std::vector<SimplexId> &cells,
                            std::vector<char> &cells_pairTypes) const;

    protected:
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

      const void *inputScalarField_{};
      const void *inputOffsets_{};
      Triangulation *inputTriangulation_{};

      SimplexId *outputCriticalPoints_numberOfPoints_{};
      std::vector<float> *outputCriticalPoints_points_{};
      std::vector<char> *outputCriticalPoints_points_cellDimensions_{};
      std::vector<SimplexId> *outputCriticalPoints_points_cellIds_{};
      void *outputCriticalPoints_points_cellScalars_{};
      std::vector<char> *outputCriticalPoints_points_isOnBoundary_{};
      std::vector<SimplexId>
        *outputCriticalPoints_points_PLVertexIdentifiers_{};
      std::vector<SimplexId> *outputCriticalPoints_points_manifoldSize_{};

      std::vector<std::array<Cell, 2>> *outputPersistencePairs_{};

      // index of vertices sorted by ascending order
      std::vector<size_t> vertsOrder_{};
    };

  } // namespace dcg
} // namespace ttk

#include <DiscreteGradient_Template.h>
