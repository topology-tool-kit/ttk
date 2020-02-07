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
                       const std::array<SimplexId, 3> &&lowVerts,
                       const std::array<uint8_t, 3> &&faces)
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
                              const dataType *scalars) const;

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
       * @param[in] scalars Scalar field
       * @param[in] offset Offset field (for comparing vertices when on a scalar
       * field plateau)
       *
       * @return Lower star as 4 sets of cells (0-cells, 1-cells, 2-cells and
       * 3-cells)
       */
      template <typename dataType, typename idType>
      inline lowerStarType lowerStar(const SimplexId a,
                                     const dataType *const scalars,
                                     const idType *const offsets) const {
        lowerStarType res{};

        // a belongs to its lower star
        res[0].emplace_back(CellExt{0, a});

        const auto sosGreaterThan
          = [&scalars, &offsets](const SimplexId m, const SimplexId n) {
              if(scalars[m] != scalars[n]) {
                return scalars[m] > scalars[n];
              } else {
                return offsets[m] > offsets[n];
              }
            };

        // store lower edges
        const auto nedges = inputTriangulation_->getVertexEdgeNumber(a);
        res[1].reserve(nedges);
        for(SimplexId i = 0; i < nedges; i++) {
          SimplexId edgeId;
          inputTriangulation_->getVertexEdge(a, i, edgeId);
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(edgeId, 0, vertexId);
          if(vertexId == a) {
            inputTriangulation_->getEdgeVertex(edgeId, 1, vertexId);
          }
          if(!sosGreaterThan(vertexId, a)) {
            res[1].emplace_back(CellExt{1, edgeId, {vertexId}, {}});
          }
        }

        if(res[1].size() < 2) {
          // at least two edges in the lower star for one triangle
          return res;
        }

        const auto processTriangle = [&](const SimplexId triangleId,
                                         const SimplexId v0, const SimplexId v1,
                                         const SimplexId v2) {
          std::array<SimplexId, 3> lowVerts{};
          if(v0 == a) {
            lowVerts[0] = v1;
            lowVerts[1] = v2;
          } else if(v1 == a) {
            lowVerts[0] = v0;
            lowVerts[1] = v2;
          } else if(v2 == a) {
            lowVerts[0] = v0;
            lowVerts[1] = v1;
          }
          if(sosGreaterThan(a, lowVerts[0]) && sosGreaterThan(a, lowVerts[1])) {
            uint8_t j{}, k{};
            // store edges indices of current triangle
            std::array<uint8_t, 3> faces{};
            for(const auto &e : res[1]) {
              if(e.lowVerts_[0] == lowVerts[0]
                 || e.lowVerts_[0] == lowVerts[1]) {
                faces[k++] = j;
              }
              j++;
            }
            res[2].emplace_back(
              CellExt{2, triangleId, std::move(lowVerts), std::move(faces)});
          }
        };

        if(dimensionality_ == 2) {
          // store lower triangles

          // use optimised triangulation methods:
          // getVertexStar instead of getVertexTriangle
          // getCellVertex instead of getTriangleVertex
          const auto ncells = inputTriangulation_->getVertexStarNumber(a);
          res[2].reserve(ncells);
          for(SimplexId i = 0; i < ncells; ++i) {
            SimplexId cellId;
            inputTriangulation_->getVertexStar(a, i, cellId);
            SimplexId v0{}, v1{}, v2{};
            inputTriangulation_->getCellVertex(cellId, 0, v0);
            inputTriangulation_->getCellVertex(cellId, 1, v1);
            inputTriangulation_->getCellVertex(cellId, 2, v2);
            processTriangle(cellId, v0, v1, v2);
          }
        } else if(dimensionality_ == 3) {
          // store lower triangles
          const auto ntri = inputTriangulation_->getVertexTriangleNumber(a);
          res[2].reserve(ntri);
          for(SimplexId i = 0; i < ntri; i++) {
            SimplexId triangleId;
            inputTriangulation_->getVertexTriangle(a, i, triangleId);
            SimplexId v0{}, v1{}, v2{};
            inputTriangulation_->getTriangleVertex(triangleId, 0, v0);
            inputTriangulation_->getTriangleVertex(triangleId, 1, v1);
            inputTriangulation_->getTriangleVertex(triangleId, 2, v2);
            processTriangle(triangleId, v0, v1, v2);
          }

          // at least three triangles in the lower star for one tetra
          if(res[2].size() >= 3) {
            // store lower tetra
            const auto ncells = inputTriangulation_->getVertexStarNumber(a);
            res[3].reserve(ncells);
            for(SimplexId i = 0; i < ncells; ++i) {
              SimplexId cellId;
              inputTriangulation_->getVertexStar(a, i, cellId);
              std::array<SimplexId, 3> lowVerts{};
              SimplexId v0{}, v1{}, v2{}, v3{};
              inputTriangulation_->getCellVertex(cellId, 0, v0);
              inputTriangulation_->getCellVertex(cellId, 1, v1);
              inputTriangulation_->getCellVertex(cellId, 2, v2);
              inputTriangulation_->getCellVertex(cellId, 3, v3);
              if(v0 == a) {
                lowVerts[0] = v1;
                lowVerts[1] = v2;
                lowVerts[2] = v3;
              } else if(v1 == a) {
                lowVerts[0] = v0;
                lowVerts[1] = v2;
                lowVerts[2] = v3;
              } else if(v2 == a) {
                lowVerts[0] = v0;
                lowVerts[1] = v1;
                lowVerts[2] = v3;
              } else if(v3 == a) {
                lowVerts[0] = v0;
                lowVerts[1] = v1;
                lowVerts[2] = v2;
              }
              if(sosGreaterThan(a, lowVerts[0])
                 && sosGreaterThan(a, lowVerts[1])
                 && sosGreaterThan(a, lowVerts[2])) {
                uint8_t j{}, k{};
                // store triangles indices of current tetra
                std::array<uint8_t, 3> faces{};
                for(const auto &t : res[2]) {
                  if((t.lowVerts_[0] == lowVerts[0]
                      || t.lowVerts_[0] == lowVerts[1]
                      || t.lowVerts_[0] == lowVerts[2])
                     && (t.lowVerts_[1] == lowVerts[0]
                         || t.lowVerts_[1] == lowVerts[1]
                         || t.lowVerts_[1] == lowVerts[2])) {
                    faces[k++] = j;
                  }
                  j++;
                }

                res[3].emplace_back(
                  CellExt{3, cellId, std::move(lowVerts), std::move(faces)});
              }
            }
          }
        }

        return res;
      }

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
      inline void pairCells(CellExt &alpha, CellExt &beta) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        char localBId{0}, localAId{0};
        SimplexId a{}, b{};

        if(beta.dim_ == 1) {

          for(SimplexId i = 0; i < 2; ++i) {
            inputTriangulation_->getEdgeVertex(beta.id_, i, a);
            if(a == alpha.id_) {
              localAId = i;
              break;
            }
          }
          const auto nedges
            = inputTriangulation_->getVertexEdgeNumber(alpha.id_);
          for(SimplexId i = 0; i < nedges; ++i) {
            inputTriangulation_->getVertexEdge(alpha.id_, i, b);
            if(b == beta.id_) {
              localBId = i;
            }
          }
        } else if(beta.dim_ == 2) {
          for(SimplexId i = 0; i < 3; ++i) {
            inputTriangulation_->getTriangleEdge(beta.id_, i, a);
            if(a == alpha.id_) {
              localAId = i;
              break;
            }
          }
          const auto ntri
            = inputTriangulation_->getEdgeTriangleNumber(alpha.id_);
          for(SimplexId i = 0; i < ntri; ++i) {
            inputTriangulation_->getEdgeTriangle(alpha.id_, i, b);
            if(b == beta.id_) {
              localBId = i;
            }
          }
        } else {
          for(SimplexId i = 0; i < 4; ++i) {
            inputTriangulation_->getCellTriangle(beta.id_, i, a);
            if(a == alpha.id_) {
              localAId = i;
              break;
            }
          }
          const auto ntetra
            = inputTriangulation_->getTriangleStarNumber(alpha.id_);
          for(SimplexId i = 0; i < ntetra; ++i) {
            inputTriangulation_->getTriangleStar(alpha.id_, i, b);
            if(b == beta.id_) {
              localBId = i;
            }
          }
        }
        gradient_[alpha.dim_][alpha.dim_][alpha.id_] = localBId;
        gradient_[alpha.dim_][alpha.dim_ + 1][beta.id_] = localAId;
#else
        gradient_[alpha.dim_][alpha.dim_][alpha.id_] = beta.id_;
        gradient_[alpha.dim_][alpha.dim_ + 1][beta.id_] = alpha.id_;
#endif // TTK_ENABLE_DCG_OPTIMIZE_MEMORY
        alpha.paired_ = true;
        beta.paired_ = true;
      }

      /**
       * Implements the ProcessLowerStars algorithm from "Theory and
       * Algorithms for Constructing Discrete Morse Complexes from
       * Grayscale Digital Images", V. Robins, P. J. Wood,
       * A. P. Sheppard
       */
      template <typename dataType, typename idType>
      int processLowerStars(const dataType *scalars, const idType *offsets);

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
      bool getDescendingPathThroughWall(const wallId_t wallId,
                                        const Cell &saddle2,
                                        const Cell &saddle1,
                                        const std::vector<wallId_t> &isVisited,
                                        std::vector<Cell> *const vpath,
                                        const bool stopIfMultiConnected = false,
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
                                       const bool stopIfMultiConnected = false,
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
       * Get the vertex id of with the maximum scalar field value on
       * the given cell. Compare offsets if scalar field is constant.
       */
      template <typename dataType, typename idType>
      inline SimplexId getCellGreaterVertex(const Cell c) const;

      /**
       * Build the geometric embedding of the given STL vector of cells.
       * The output data pointers are modified accordingly. This
       * function needs the following internal pointers to be set:
       * outputCriticalPoints_numberOfPoints_
       * outputCriticalPoints_points_
       * inputScalarField_
       */
      template <typename dataType, typename idType>
      int setCriticalPoints(const std::vector<Cell> &criticalPoints,
                            std::vector<size_t> &nCriticalPointsByDim) const;

      /**
       * Detect the critical points and build their geometric embedding.
       * The output data pointers are modified accordingly.
       */
      template <typename dataType, typename idType>
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

      std::vector<std::tuple<Cell, Cell>> *outputPersistencePairs_{};
    };

  } // namespace dcg
} // namespace ttk

#include <DiscreteGradient_Template.h>
