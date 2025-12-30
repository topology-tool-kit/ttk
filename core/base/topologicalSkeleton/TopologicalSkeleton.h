/// \ingroup base
/// \class ttk::TopologicalSkeleton
/// \author Tanner Finken <finkent@arizona.edu>
/// \author Joshua A. Levine <josh@cs.arizona.edu>
/// \date Summer 2024.
///
/// \brief TTK processing package for the computation of Topological Skeleton of
/// Vector Fields. The skeleton consists of four components: critical points,
/// 1-separatrices(including orbits), 2-separatrices, and segmentation of where
/// cells flow to and from.
///
/// The code implementation is largely based on the MorseSmaleComplex filter.
/// Additionally, the user is allowed to specify a threshold to simplify down
/// to a certain number of critical points in the DiscreteVectorField.
///
/// \sa ttk::MorseSmaleComplex
/// \sa ttk::dcvf::DiscreteVectorField
/// \sa ttk::VectorSimplification
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE VIS 2024.
///
/// \sa ttk::TopologicalSkeleton
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

// base code includes
#include <DiscreteVectorField.h>
#include <Triangulation.h>
#include <VectorSimplification.h>

#include <queue>

namespace ttk {
  class TopologicalSkeleton : public virtual Debug {
  public:
    TopologicalSkeleton();

    /** @brief Critical points data arrays */
    struct OutputCriticalPoints {
      std::vector<std::array<float, 3>> points_{};
      std::vector<char> cellDimensions_{};
      std::vector<SimplexId> cellIds_{};
      std::vector<char> isOnBoundary_{};
      std::vector<SimplexId> PLVertexIdentifiers_{};
      std::vector<SimplexId> manifoldSize_{};
      void clear() {
        *this = {};
      };
    };

    /** @brief 1-Separatrices point and cell data arrays */
    struct Output1Separatrices {
      struct {
        SimplexId numberOfPoints_{};
        std::vector<float> points_{};
        std::vector<char> smoothingMask_{};
        std::vector<char> cellDimensions_{};
        std::vector<ttk::SimplexId> cellIds_{};
      } pt{}; // point data arrays
      struct {
        SimplexId numberOfCells_{};
        std::vector<ttk::SimplexId> connectivity_{};
        std::vector<ttk::SimplexId> sourceIds_{};
        std::vector<ttk::SimplexId> destinationIds_{};
        std::vector<ttk::SimplexId> separatrixIds_{};
        std::vector<char> separatrixTypes_{};
        std::vector<char> isOnBoundary_{};
        std::vector<char> isOrbit_{};
      } cl{}; // cell data arrays
      void clear() {
        *this = {};
      };
    };

    /** @brief 2-Separatrices point and cell data arrays */
    struct Output2Separatrices {
      struct {
        SimplexId numberOfPoints_{};
        std::vector<float> points_{};
      } pt{}; // point data arrays
      struct {
        SimplexId numberOfCells_{};
        std::vector<ttk::SimplexId> offsets_{};
        std::vector<ttk::SimplexId> connectivity_{};
        std::vector<ttk::SimplexId> sourceIds_{};
        std::vector<ttk::SimplexId> separatrixIds_{};
        std::vector<char> separatrixTypes_{};
        std::vector<char> isOnBoundary_{};
      } cl{}; // cell data arrays
      void clear() {
        *this = {};
      };
    };

    /** @brief Pointers to pre-allocated segmentation point data arrays */
    struct OutputManifold {
      SimplexId *ascending_;
      SimplexId *descending_;
      SimplexId *morseSmale_;
    };

    /**
     * Main function for computing the Topological Skeleton.
     *
     * @pre TopologicalSkeleton::preconditionTriangulation must be
     * called prior to this.
     */
    template <typename dataType, typename triangulationType>
    inline int execute(OutputCriticalPoints &outCP,
                       Output1Separatrices &outSeps1,
                       Output2Separatrices &outSeps2,
                       OutputManifold &outManifold,
                       const dataType *const vectors,
                       const size_t vectorsMTime,
                       const triangulationType &triangulation);

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the critical points.
     */
    inline void setComputeCriticalPoints(const bool state) {
      this->ComputeCriticalPoints = state;
    }
    /**
     * Enable/Disable computation of the geometrical embedding of
     * the 1-separatrices.
     */
    inline void setComputeSeparatrices1(const bool doAscending,
                                        const bool doDescending,
                                        const bool doSaddleConnectors) {
      this->ComputeAscendingSeparatrices1 = doAscending;
      this->ComputeDescendingSeparatrices1 = doDescending;
      this->ComputeSaddleConnectors = doSaddleConnectors;
    }

    /**
     * Enable/Disable computation of geometrical embedding of 1-cycles
     */
    inline void setComputeCycles1(const bool doAttracting,
                                  const bool doRepelling) {
      this->ComputeAttractingCycles1 = doAttracting;
      this->ComputeRepellingCycles1 = doRepelling;
    }
    /**
     * Enable/Disable computation of the geometrical embedding of
     * the 2-separatrices (disabled by default).
     */
    inline void setComputeSeparatrices2(const bool doAscending,
                                        const bool doDescending) {
      this->ComputeAscendingSeparatrices2 = doAscending;
      this->ComputeDescendingSeparatrices2 = doDescending;
    }
    /**
     * Enable/Disable computation of the geometrical embedding of the
     * manifolds of the critical points.
     */
    inline void setComputeSegmentation(const bool doAscending,
                                       const bool doDescending,
                                       const bool doMorseSmale) {
      this->ComputeAscendingSegmentation = doAscending;
      this->ComputeDescendingSegmentation = doDescending;
      this->ComputeFinalSegmentation = doMorseSmale;
    }

    /**
     * Enable/Disable post-processing vector field reversal of
     * the (saddle,...,saddle) vpaths under a given persistence
     * threshold (disabled by default).
     */
    inline void setRunSimplification(const bool state) {
      RunSimplification = state;
    }

    inline void setFullOrbits(const bool fullOrbits) {
      this->simplifierField_.setFullOrbitSimplification(fullOrbits);
    }

    /**
     * Set the threshold value for post-processing of
     * (saddle,...,saddle) vpaths vector field reversal
     * (default value is 0.0).
     */
    inline void setSimplificationThreshold(const double threshold) {
      SimplificationThreshold = threshold;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      this->simplifierField_.preconditionTriangulation(data);
      data->preconditionCellEdges();
      data->preconditionCellNeighbors();
    }

  protected:
    /**
     * Utility class representing Ridge lines, Valley lines
     * Oribits, and Saddle connectors.
     */
    struct Separatrix {
      /** Source cell of the separatrix. */
      dcg::Cell source_;
      /** Destination cell of the separatrix. */
      dcg::Cell destination_;
      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<dcg::Cell> geometry_;
    };

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * vector field.
     */
    template <typename dataType, typename triangulationType>
    int
      getDescendingSeparatrices1(const std::vector<SimplexId> &saddles,
                                 std::vector<Separatrix> &separatrices,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the ascending 1-separatrices by reading into the discrete
     * vector field.
     */
    template <typename dataType, typename triangulationType>
    int getAscendingSeparatrices1(const std::vector<SimplexId> &saddles,
                                  std::vector<Separatrix> &separatrices,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * vector field.
     */
    template <typename triangulationType>
    int getSaddleConnectors(const std::vector<SimplexId> &saddles2,
                            std::vector<Separatrix> &separatrices,
                            const triangulationType &triangulation) const;

    /**
     * Compute the attracting 1-cycles by reading into the discrete
     * vector field.
     */
    template <typename triangulationType>
    int getAttractingCycles1(std::vector<Separatrix> &separatrices,
                             const triangulationType &triangulation) const;

    /**
     * Compute the repelling 1-cycles by reading into the discrete
     * vector field.
     */
    template <typename triangulationType>
    int getRepellingCycles1(std::vector<Separatrix> &separatrices,
                            const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the 1-separatrices.
     */
    template <typename dataType, typename triangulationType>
    int setSeparatrices1(Output1Separatrices &outSeps1,
                         const std::vector<Separatrix> &separatrices,
                         const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * vector field from the maxima.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices2(
      const std::vector<SimplexId> &saddles2,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the descending
     * 2-separatrices.
     */
    template <typename triangulationType>
    int setDescendingSeparatrices2(
      Output2Separatrices &outSeps2,
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<SimplexId>> &separatricesSaddles,
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
     * vector field from the minima.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices2(
      const std::vector<SimplexId> &saddles1,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the ascending
     * 2-separatrices.
     */
    template <typename triangulationType>
    int setAscendingSeparatrices2(
      Output2Separatrices &outSeps2,
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<Separatrix>> &separatrices) const;

    /**
     * Compute the ascending manifold of the 'sources'
     * (maxima for gradient fields).
     */
    template <typename triangulationType>
    int setAscendingSegmentation(const std::vector<SimplexId> &maxima,
                                 std::vector<Separatrix> &repellingOrbits,
                                 SimplexId *const morseSmaleManifold,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the descending manifold of the 'sinks'
     * (minima for gradient fields).
     */
    template <typename triangulationType>
    int setDescendingSegmentation(const std::vector<SimplexId> &minima,
                                  std::vector<Separatrix> &attractingOrbits,
                                  SimplexId *const morseSmaleManifold,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the final combinatorial Morse-Smale complex
     * segmentation(intersection of outward and inward flowing segmentation).
     */
    template <typename triangulationType>
    int setFinalSegmentation(const SimplexId numberOfMaxima,
                             const SimplexId *const ascendingManifold,
                             const SimplexId *const descendingManifold,
                             SimplexId *const morseSmaleManifold,
                             const triangulationType &triangulation) const;

    VectorSimplification simplifierField_{};

    bool ComputeCriticalPoints{true};
    bool ComputeAscendingSeparatrices1{true};
    bool ComputeDescendingSeparatrices1{true};
    bool ComputeSaddleConnectors{false};
    bool ComputeAttractingCycles1{true};
    bool ComputeRepellingCycles1{true};
    bool ComputeAscendingSeparatrices2{false};
    bool ComputeDescendingSeparatrices2{false};
    bool ComputeAscendingSegmentation{true};
    bool ComputeDescendingSegmentation{true};
    bool ComputeFinalSegmentation{true};

    bool RunSimplification{false};
    double SimplificationThreshold{};
    bool ReverseFullOrbit{true};
  };
} // namespace ttk

// ---------------- //
//  Execute method  //
// ---------------- //

template <typename dataType, typename triangulationType>
int ttk::TopologicalSkeleton::execute(OutputCriticalPoints &outCP,
                                      Output1Separatrices &outSeps1,
                                      Output2Separatrices &outSeps2,
                                      OutputManifold &outManifold,
                                      const dataType *const vectors,
                                      const size_t vectorsMTime,
                                      const triangulationType &triangulation) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vectors == nullptr) {
    this->printErr("Input vector field pointer is null.");
    return -1;
  }
#endif
  TTK_FORCE_USE(outManifold); // Remove unused param error.
  Timer t;

  outCP.clear();
  outSeps1.clear();
  outSeps2.clear();
  const auto dim = triangulation.getDimensionality();

  this->simplifierField_.setThreadNumber(threadNumber_);
  this->simplifierField_.setDebugLevel(debugLevel_);
  this->simplifierField_.buildField<dataType, triangulationType>(
    vectors, vectorsMTime, triangulation);
  // First simplify the discrete field as desired
  if(this->RunSimplification) {
    int persistenceThreshold = static_cast<int>(this->SimplificationThreshold);
    std::vector<ttk::VectorSimplification::PlotPoint> emptyPlot;
    this->simplifierField_.performSimplification<dataType, triangulationType>(
      persistenceThreshold, false, emptyPlot, triangulation);
  }

  std::array<std::vector<SimplexId>, 4> criticalPoints{};
  if(ComputeCriticalPoints) {
    Timer tm{};
    this->simplifierField_.dcvf_.getCriticalPoints(
      criticalPoints, triangulation);
    this->printMsg("  Critical points extracted", 1.0, tm.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);
  }

  std::vector<std::vector<Separatrix>> separatrices1{};

  // 1-separatrices
  Timer tm1sep{};

  if(dim > 1 && ComputeDescendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    getDescendingSeparatrices1<dataType, triangulationType>(
      criticalPoints[1], separatrices1.back(), triangulation);

    this->printMsg("  Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(dim > 1 && ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();

    getAscendingSeparatrices1<dataType, triangulationType>(
      criticalPoints[dim - 1], separatrices1.back(), triangulation);

    this->printMsg("  Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  std::vector<Separatrix> attractingCycles{};

  if(dim > 1 && (ComputeAttractingCycles1 || ComputeDescendingSegmentation)) {
    Timer tmp;
    // separatrices1.emplace_back();

    getAttractingCycles1(attractingCycles, triangulation);
    if(ComputeAttractingCycles1)
      separatrices1.emplace_back(attractingCycles);
    this->printMsg("  Attracting 1-cycles computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);
  }

  std::vector<Separatrix> repellingCycles{};

  if(dim > 1 && (ComputeRepellingCycles1 || ComputeAscendingSegmentation)) {
    Timer tmp;
    // separatrices1.emplace_back();

    getRepellingCycles1(repellingCycles, triangulation);
    if(ComputeRepellingCycles1)
      separatrices1.emplace_back(repellingCycles);
    this->printMsg("  Repelling 1-cycles computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);
  }

  // saddle-connectors
  if(dim == 3 && ComputeSaddleConnectors) {
    Timer tmp;
    separatrices1.emplace_back();

    getSaddleConnectors(criticalPoints[2], separatrices1.back(), triangulation);

    this->printMsg("  Saddle connectors computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);
  }

  if(dim > 1
     && (ComputeDescendingSeparatrices1 || ComputeAscendingSeparatrices1
         || ComputeSaddleConnectors || ComputeAttractingCycles1
         || ComputeRepellingCycles1)) {
    Timer tmp{};
    if(separatrices1.size() > 0) {
      this->flattenSeparatricesVectors(separatrices1);
      setSeparatrices1<dataType, triangulationType>(
        outSeps1, separatrices1[0], triangulation);
    }

    this->printMsg("  1-separatrices set", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);

    this->printMsg("1-separatrices computed", 1.0, tm1sep.getElapsedTime(),
                   this->threadNumber_);
  }

  // 2-separatrices
  Timer tm2sep{};

  if(dim == 3 && ComputeDescendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<SimplexId>> separatricesSaddles;
    getDescendingSeparatrices2(
      criticalPoints[2], separatrices, separatricesSaddles, triangulation);
    setDescendingSeparatrices2(
      outSeps2, separatrices, separatricesSaddles, triangulation);

    this->printMsg("  Descending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(dim == 3 && ComputeAscendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<SimplexId>> separatricesSaddles;
    getAscendingSeparatrices2(
      criticalPoints[1], separatrices, separatricesSaddles, triangulation);
    setAscendingSeparatrices2(
      outSeps2, separatrices, separatricesSaddles, triangulation);

    this->printMsg("  Ascending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(this->ComputeAscendingSeparatrices2
     || this->ComputeDescendingSeparatrices2) {
    this->printMsg("2-separatrices computed", 1.0, tm2sep.getElapsedTime(),
                   this->threadNumber_);
  }

  if(ComputeAscendingSegmentation || ComputeDescendingSegmentation) {
    Timer tmp;

    if(ComputeAscendingSegmentation && criticalPoints[dim].size() > 0) {
      setAscendingSegmentation(criticalPoints[dim], repellingCycles,
                               outManifold.ascending_, triangulation);
    }
    if(ComputeDescendingSegmentation && criticalPoints[0].size() > 0) {
      setDescendingSegmentation(criticalPoints[0], attractingCycles,
                                outManifold.descending_, triangulation);
    }
    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation
       && ComputeFinalSegmentation && criticalPoints[dim].size() > 0
       && criticalPoints[0].size() > 0) {
      SimplexId numSources = static_cast<SimplexId>(criticalPoints[dim].size()
                                                    + repellingCycles.size());
      setFinalSegmentation(numSources, outManifold.ascending_,
                           outManifold.descending_, outManifold.morseSmale_,
                           triangulation);
    }

    this->printMsg(
      "Segmentation computed", 1.0, tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeCriticalPoints) {
    this->simplifierField_.dcvf_.setCriticalPoints<dataType, triangulationType>(
      criticalPoints, outCP.points_, outCP.cellDimensions_, outCP.cellIds_,
      outCP.isOnBoundary_, outCP.PLVertexIdentifiers_, triangulation);

    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation) {
      this->simplifierField_.dcvf_.setManifoldSize(
        criticalPoints, outManifold.ascending_, outManifold.descending_,
        outCP.manifoldSize_);
    }
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

// ---------------- //
//  1-Separatrices  //
// ---------------- //

template <typename dataType, typename triangulationType>
int ttk::TopologicalSkeleton::getDescendingSeparatrices1(
  const std::vector<SimplexId> &saddles,
  std::vector<Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  const SimplexId numberOfSaddles = saddles.size();

  // only 2 descending separatrices per 1-saddle
  const SimplexId numberOfSeparatrices = 2 * numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const Cell saddle{1, saddles[i]};

    // add descending vpaths
    {
      const Cell &saddle1 = saddle;

      for(int j = 0; j < 2; ++j) {
        SimplexId vertexId;
        triangulation.getEdgeVertex(saddle1.id_, j, vertexId);

        std::vector<Cell> vpath;
        vpath.emplace_back(saddle1);
        simplifierField_.dcvf_.getDescendingPath<dataType, triangulationType>(
          Cell(0, vertexId), vpath, triangulation, true);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 0
           and simplifierField_.dcvf_.isCellCritical(lastCell)) {
          separatrices[2 * i + j].source_ = saddle;
          separatrices[2 * i + j].destination_ = lastCell;
          separatrices[2 * i + j].geometry_ = std::move(vpath);
        } else { // A cycle
          separatrices[2 * i + j].source_ = saddle;
          separatrices[2 * i + j].destination_ = Cell(0, -1); // Denote a cycle
          separatrices[2 * i + j].geometry_ = std::move(vpath);
        }
      }
    }
  }

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::TopologicalSkeleton::getAscendingSeparatrices1(
  const std::vector<SimplexId> &saddles,
  std::vector<Separatrix> &separatrices,
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

  const SimplexId numberOfSaddles = saddles.size();

  std::vector<std::vector<Separatrix>> sepsPerSaddle(numberOfSaddles);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const Cell saddle{dim - 1, saddles[i]};

    // add ascending vpaths
    const auto starNumber{(triangulation.*getFaceStarNumber)(saddle.id_)};
    for(SimplexId j = 0; j < starNumber; ++j) {

      SimplexId sId{};
      (triangulation.*getFaceStar)(saddle.id_, j, sId);

      std::vector<Cell> vpath{saddle};
      simplifierField_.dcvf_.getAscendingPath<dataType, triangulationType>(
        Cell(dim, sId), vpath, triangulation, true);

      const Cell &lastCell = vpath.back();
      if(lastCell.dim_ == dim
         and simplifierField_.dcvf_.isCellCritical(lastCell)) {
        sepsPerSaddle[i].emplace_back();
        sepsPerSaddle[i].back().source_ = saddle;
        sepsPerSaddle[i].back().destination_ = lastCell;
        sepsPerSaddle[i].back().geometry_ = std::move(vpath);
      } else { // A Cycle or border
        sepsPerSaddle[i].emplace_back();
        sepsPerSaddle[i].back().source_ = saddle;
        sepsPerSaddle[i].back().destination_ = Cell(dim, -1);
        sepsPerSaddle[i].back().geometry_ = std::move(vpath);
      }
    }
  }

  if(numberOfSaddles != 0) {
    this->flattenSeparatricesVectors(sepsPerSaddle);

    separatrices = std::move(sepsPerSaddle[0]);
  }
  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::getSaddleConnectors(
  const std::vector<SimplexId> &saddles2,
  std::vector<Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  const auto nTriangles = triangulation.getNumberOfTriangles();
  // visited triangles (one vector per thread)
  std::vector<bool> isVisited(nTriangles, false);
  std::vector<SimplexId> visitedTriangles{};

  using Vpath = std::vector<Cell>;

  const auto dim{triangulation.getDimensionality()};

  std::vector<std::vector<Separatrix>> sepsByThread(saddles2.size());
  std::vector<SimplexId> saddles1{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(isVisited, visitedTriangles, saddles1)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    const Cell s2{dim - 1, saddles2[i]};

    VisitedMask mask{isVisited, visitedTriangles};
    simplifierField_.dcvf_.getDescendingWall(
      s2, mask, triangulation, nullptr, &saddles1);

    for(const auto saddle1Id : saddles1) {
      const Cell s1{1, saddle1Id};

      Vpath vpath;
      simplifierField_.dcvf_.getAscendingPathThroughWall(
        s1, s2, isVisited, &vpath, triangulation);

      if(vpath.empty()) {
        // safety, should be unreachable
        continue;
      }
      const auto &last = vpath.back();

      if(last.dim_ == s2.dim_ && last.id_ == s2.id_) {
        sepsByThread[i].emplace_back();
        sepsByThread[i].back().source_ = s1;
        sepsByThread[i].back().destination_ = s2;
        sepsByThread[i].back().geometry_ = std::move(vpath);
      }
    }
  }
  if(saddles2.size() != 0) {
    this->flattenSeparatricesVectors(sepsByThread);

    separatrices = std::move(sepsByThread[0]);
  }

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::getAttractingCycles1(
  std::vector<Separatrix> &separatrices,
  const triangulationType &triangulation) const {
  // Store the vpaths of cycles found by minimum vertex id
  std::map<SimplexId, std::vector<Cell>> minToCycles;

  const auto nVerts{triangulation.getNumberOfVertices()};

  std::vector<SimplexId> visited{};
  std::vector<char> isCycle;
  isCycle.resize(nVerts, 0);
  std::vector<char> hasChecked;
  hasChecked.resize(nVerts, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(visited, isCycle)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    if(hasChecked[i] != 0) {
      continue;
    }
    visited.clear();
    auto curr{i};
    while(hasChecked[curr] == 0) {
      // Check if a cycle has occured
      if(isCycle[curr] == 1) {
        std::vector<Cell> cyclePath{Cell{0, curr}};
        SimplexId currentMin{curr};
        while(visited.back() != curr) {
          cyclePath.emplace_back(Cell{0, visited.back()});

          hasChecked[visited.back()] = 1;
          if(visited.back() < currentMin) {
            currentMin = visited.back();
          }
          visited.pop_back();
        }
        cyclePath.emplace_back(Cell{0, curr});
// Critical section to safely update minToCycles
#pragma omp critical
        {
          if(minToCycles[currentMin].size() == 0) { // Not assigned key
            minToCycles[currentMin] = cyclePath;
          }
        }
        break;
      }
      if(this->simplifierField_.dcvf_.isCellCritical(Cell{0, curr})) {
        break;
      }
      // follow a V-path till an already checked vertex is reached
      const auto pairedEdge{this->simplifierField_.dcvf_.getPairedCell(
        Cell{0, curr}, triangulation)};
      SimplexId next{};
      triangulation.getEdgeVertex(pairedEdge, 0, next);
      if(next == curr) {
        triangulation.getEdgeVertex(pairedEdge, 1, next);
      }
      visited.emplace_back(curr);
      isCycle[curr] = 1;
      curr = next;
    }
    for(const auto el : visited) {
      hasChecked[el] = 1;
      isCycle[el] = 0;
    }
  } // End of Parallel for loop

  // Extract VPaths(iterate over cycles)
  for(auto it = minToCycles.begin(); it != minToCycles.end(); ++it) {
    auto vpath = it->second;
    auto startCell = vpath.back();
    separatrices.emplace_back();
    separatrices.back().source_ = startCell;
    separatrices.back().destination_ = startCell;
    separatrices.back().geometry_ = std::move(vpath);
  }

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::getRepellingCycles1(
  std::vector<Separatrix> &separatrices,
  const triangulationType &triangulation) const {
  // Store the vpaths of cycles found by minimum cell id
  std::map<SimplexId, std::vector<Cell>> minToCycles;

  const auto nCells{triangulation.getNumberOfCells()};

  const auto dim{triangulation.getDimensionality()};

  // Triangulation method pointers for 3D
  auto getFaceStarNumber = &triangulationType::getTriangleStarNumber;
  auto getFaceStar = &triangulationType::getTriangleStar;
  if(dim == 2) {
    // Triangulation method pointers for 2D
    getFaceStarNumber = &triangulationType::getEdgeStarNumber;
    getFaceStar = &triangulationType::getEdgeStar;
  } else if(dim == 1) {
    // Triangulation method pointers for 1D
    getFaceStarNumber = &triangulationType::getVertexStarNumber;
    getFaceStar = &triangulationType::getVertexStar;
  }

  // cells visited during the propagation alongside one integral line
  std::vector<SimplexId> visited{};
  // all marked cells
  std::vector<char> hasChecked;
  hasChecked.resize(nCells, 0);
  // cycle detection array
  std::vector<char> isCycle;
  isCycle.resize(nCells, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(visited, isCycle)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nCells; ++i) {
    if(hasChecked[i] != 0) {
      continue;
    }
    visited.clear();
    auto curr{i};
    while(hasChecked[curr] == 0) {
      if(isCycle[curr] == 1) {
        std::vector<Cell> cyclePath{Cell{dim, curr}};
        SimplexId currentMin = curr;
        while(visited.back() != curr) {
          cyclePath.emplace_back(Cell{dim, visited.back()});
          hasChecked[visited.back()] = 1;
          if(visited.back() < currentMin) {
            currentMin = visited.back();
          }
          visited.pop_back();
        }
        cyclePath.emplace_back(Cell{dim, curr});
// Critical section to safely update minToCycles
#pragma omp critical
        {
          if(minToCycles[currentMin].size() == 0) { // Not assigned key
            minToCycles[currentMin] = cyclePath;
          }
        }
        break;
      }
      // Handle critical cell
      if(this->simplifierField_.dcvf_.isCellCritical(Cell{dim, curr})) {
        break;
      }
      // follow a V-path till an already marked cell is reached
      // Break when reaching a cycle
      const auto paired{this->simplifierField_.dcvf_.getPairedCell(
        Cell{dim, curr}, triangulation, true)};
      SimplexId next{curr};
      const auto nStars{(triangulation.*getFaceStarNumber)(paired)};
      for(SimplexId j = 0; j < nStars; ++j) {
        (triangulation.*getFaceStar)(paired, j, next);
        // get the first cell != curr (what of non-manifold datasets?)
        if(next != curr) {
          break;
        }
      }
      visited.emplace_back(curr);
      isCycle[curr] = 1;
      if(next == curr) {
        // on the boundary?
        break;
      }
      curr = next;
    }
    for(const auto el : visited) {
      hasChecked[el] = 1;
      isCycle[el] = 0;
    }
  } // End of Parallel For Loop

  // Extract VPaths(iterate over cycles)
  for(auto it = minToCycles.begin(); it != minToCycles.end(); ++it) {
    auto vpath = it->second;
    auto startCell = vpath.back();
    separatrices.emplace_back();
    separatrices.back().source_ = startCell;
    separatrices.back().destination_ = startCell;
    separatrices.back().geometry_ = std::move(vpath);
  }

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::TopologicalSkeleton::setSeparatrices1(
  Output1Separatrices &outSeps1,
  const std::vector<Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  // max existing separatrix id + 1 or 0
  const SimplexId separatrixId
    = !outSeps1.cl.separatrixIds_.empty()
        ? *std::max_element(outSeps1.cl.separatrixIds_.begin(),
                            outSeps1.cl.separatrixIds_.end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(outSeps1.pt.numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(outSeps1.cl.numberOfCells_)};
  // points beginning id for each separatrix geometry
  std::vector<size_t> geomPointsBegId{npoints};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    const auto sepSize = sep.geometry_.size();
    npoints += sepSize;
    ncells += sepSize - 1;
    geomPointsBegId.emplace_back(npoints);
    geomCellsBegId.emplace_back(ncells);
  }

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;

  // resize arrays
  outSeps1.pt.points_.resize(3 * npoints);
  auto &points = outSeps1.pt.points_;
  outSeps1.cl.connectivity_.resize(2 * ncells);
  auto &cellsConn = outSeps1.cl.connectivity_;
  outSeps1.pt.smoothingMask_.resize(npoints);
  outSeps1.pt.cellDimensions_.resize(npoints);
  outSeps1.pt.cellIds_.resize(npoints);
  outSeps1.cl.sourceIds_.resize(ncells);
  outSeps1.cl.destinationIds_.resize(ncells);
  outSeps1.cl.separatrixIds_.resize(ncells);
  outSeps1.cl.separatrixTypes_.resize(ncells);
  outSeps1.cl.isOnBoundary_.resize(ncells);
  outSeps1.cl.isOrbit_.resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    const auto &sepGeom = sep.geometry_;
    const auto sepId = separatrixId + i;
    // saddle (asc/desc sep) or saddle1 (saddle connector)
    const dcg::Cell &src = sep.source_;
    // extremum/ cycle (asc/desc sep) or saddle2 (saddle connector)
    const dcg::Cell &dst = sep.destination_;

    // get separatrix type
    const auto saddleConnector
      = dimensionality == 3 && src.dim_ == 1 && dst.dim_ == 2;
    const char sepType
      = saddleConnector ? 1 : std::min(dst.dim_, dimensionality - 1);

    // get boundary condition
    const auto onBoundary
      = static_cast<char>(
          simplifierField_.dcvf_.isBoundary<dataType, triangulationType>(
            src, triangulation))
        + static_cast<char>(
          simplifierField_.dcvf_.isBoundary<dataType, triangulationType>(
            dst, triangulation));

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      std::array<float, 3> pt{};
      triangulation.getCellIncenter(cell.id_, cell.dim_, pt.data());

      // index of current point in point data arrays
      const auto k = geomPointsBegId[i] + j;

      points[3 * k + 0] = pt[0];
      points[3 * k + 1] = pt[1];
      points[3 * k + 2] = pt[2];

      outSeps1.pt.smoothingMask_[k]
        = (j == 0 || j == sepGeom.size() - 1) ? 0 : 1;
      outSeps1.pt.cellDimensions_[k] = cell.dim_;
      outSeps1.pt.cellIds_[k] = cell.id_;

      // skip filling cell data for first geometry point
      if(j == 0)
        continue;

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j - 1;

      cellsConn[2 * l + 0] = k - 1;
      cellsConn[2 * l + 1] = k;

      outSeps1.cl.sourceIds_[l] = src.id_;
      outSeps1.cl.destinationIds_[l] = dst.id_;
      outSeps1.cl.separatrixIds_[l] = sepId;
      outSeps1.cl.separatrixTypes_[l] = sepType;
      outSeps1.cl.isOnBoundary_[l] = onBoundary;
      outSeps1.cl.isOrbit_[l] = static_cast<char>(src.id_ == dst.id_);
    }
  }

  // update pointers
  outSeps1.pt.numberOfPoints_ = npoints;
  outSeps1.cl.numberOfCells_ = ncells;

  return 0;
}

// ---------------- //
//  2-Separatrices  //
// ---------------- //

template <typename triangulationType>
int ttk::TopologicalSkeleton::getAscendingSeparatrices2(
  const std::vector<SimplexId> &saddles1,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  const SimplexId numberOfSaddles = saddles1.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
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
    const Cell saddle1{1, saddles1[i]};

    std::vector<Cell> wall;
    VisitedMask mask{isVisited, visitedEdges};
    simplifierField_.dcvf_.getAscendingWall(
      saddle1, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatrices[i].source_ = saddle1;
    separatrices[i].destination_ = emptyCell;
    separatrices[i].geometry_ = std::move(wall);
  }

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::getDescendingSeparatrices2(
  const std::vector<SimplexId> &saddles2,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  const SimplexId numberOfSaddles = saddles2.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const auto nTriangles = triangulation.getNumberOfTriangles();
  std::vector<bool> isVisited(nTriangles, false);
  std::vector<SimplexId> visitedTriangles{};

  const auto dim{triangulation.getDimensionality()};

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(isVisited, visitedTriangles)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const Cell saddle2{dim - 1, saddles2[i]};

    std::vector<Cell> wall;
    VisitedMask mask{isVisited, visitedTriangles};
    simplifierField_.dcvf_.getDescendingWall(
      saddle2, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatrices[i].source_ = saddle2;
    separatrices[i].destination_ = emptyCell;
    separatrices[i].geometry_ = std::move(wall);
  }

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::getDualPolygon(
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
int ttk::TopologicalSkeleton::sortDualPolygonVertices(
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
int ttk::TopologicalSkeleton::setAscendingSeparatrices2(
  Output2Separatrices &outSeps2,
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = !outSeps2.cl.separatrixIds_.empty()
        ? *std::max_element(outSeps2.cl.separatrixIds_.begin(),
                            outSeps2.cl.separatrixIds_.end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(outSeps2.pt.numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(outSeps2.cl.numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outSeps2.cl.connectivity_.size()};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    ncells += sep.geometry_.size();
    geomCellsBegId.emplace_back(ncells);
  }

  // store the separatrices info (one per separatrix)
  std::vector<SimplexId> sepSourceIds(separatrices.size());
  std::vector<SimplexId> sepIds(separatrices.size());
  std::vector<char> sepOnBoundary(separatrices.size());

  // store the polygonal cells tetras SimplexId
  std::vector<SimplexId> polygonNTetras(ncells - noldcells);
  std::vector<SimplexId> polygonEdgeIds(ncells - noldcells);
  std::vector<SimplexId> polygonSepInfosIds(ncells - noldcells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    const auto &sepGeom = sep.geometry_;
    const auto &sepSaddles = separatricesSaddles[i];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle1

    // get boundary condition
    const char onBoundary
      = (sepSaddles.empty()
           ? 0
           : std::count_if(sepSaddles.begin(), sepSaddles.end(),
                           [&triangulation](const SimplexId a) {
                             return triangulation.isTriangleOnBoundary(a);
                           }))
        + triangulation.isEdgeOnBoundary(src.id_);

    sepIds[i] = sepId;
    sepSourceIds[i] = src.id_;
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
  outSeps2.cl.connectivity_.resize(firstCellId + nnewpoints);
  auto cellsConn = &outSeps2.cl.connectivity_[firstCellId];
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

  TTK_PSORT(this->threadNumber_, cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<SimplexId> vertId2PointsId(triangulation.getNumberOfCells());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  ncells = noldcells + validTetraIds.size();

  // resize arrays
  outSeps2.pt.points_.resize(3 * npoints);
  auto points = &outSeps2.pt.points_[3 * noldpoints];
  outSeps2.cl.offsets_.resize(ncells + 1);
  outSeps2.cl.offsets_[0] = 0;
  auto cellsOff = &outSeps2.cl.offsets_[noldcells];
  outSeps2.cl.sourceIds_.resize(ncells);
  outSeps2.cl.separatrixIds_.resize(ncells);
  outSeps2.cl.separatrixTypes_.resize(ncells);
  outSeps2.cl.isOnBoundary_.resize(ncells);

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
    outSeps2.cl.sourceIds_[l] = sepSourceIds[n];
    outSeps2.cl.separatrixIds_[l] = sepIds[n];
    outSeps2.cl.separatrixTypes_[l] = 1;
    outSeps2.cl.isOnBoundary_[l] = sepOnBoundary[n];
  }

  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    // fill offsets sequentially (due to iteration dependencies)
    cellsOff[i + 1] = cellsOff[i] + polygonNTetras[validTetraIds[i]];
  }

  outSeps2.pt.numberOfPoints_ = npoints;
  outSeps2.cl.numberOfCells_ = ncells;

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::setDescendingSeparatrices2(
  Output2Separatrices &outSeps2,
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = !outSeps2.cl.separatrixIds_.empty()
        ? *std::max_element(outSeps2.cl.separatrixIds_.begin(),
                            outSeps2.cl.separatrixIds_.end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(outSeps2.pt.numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(outSeps2.cl.numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outSeps2.cl.connectivity_.size()};

  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    ncells += sep.geometry_.size();
    geomCellsBegId.emplace_back(ncells);
  }

  // resize arrays
  outSeps2.cl.offsets_.resize(ncells + 1);
  outSeps2.cl.offsets_[0] = 0;
  outSeps2.cl.connectivity_.resize(
    firstCellId + 3 * (ncells - noldcells)); // triangles cells
  auto cellsOff = &outSeps2.cl.offsets_[noldcells];
  auto cellsConn = &outSeps2.cl.connectivity_[firstCellId];
  outSeps2.cl.sourceIds_.resize(ncells);
  outSeps2.cl.separatrixIds_.resize(ncells);
  outSeps2.cl.separatrixTypes_.resize(ncells);
  outSeps2.cl.isOnBoundary_.resize(ncells);

  // store the cells/triangles vertices vertexId
  std::vector<SimplexId> cellVertsIds(3 * (ncells - noldcells));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    const auto &sepGeom = sep.geometry_;
    const auto &sepSaddles = separatricesSaddles[i];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle2
    const char sepType = 2;

    // get boundary condition
    const char onBoundary
      = (sepSaddles.empty()
           ? 0
           : std::count_if(sepSaddles.begin(), sepSaddles.end(),
                           [&triangulation](const SimplexId a) {
                             return triangulation.isEdgeOnBoundary(a);
                           }))
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

      outSeps2.cl.sourceIds_[l] = src.id_;
      outSeps2.cl.separatrixIds_[l] = sepId;
      outSeps2.cl.separatrixTypes_[l] = sepType;
      outSeps2.cl.isOnBoundary_[l] = onBoundary;
    }
  }

  // reduce the cell vertices ids
  // (cells are triangles sharing two vertices)
  TTK_PSORT(this->threadNumber_, cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<size_t> vertId2PointsId(triangulation.getNumberOfVertices());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  outSeps2.pt.points_.resize(3 * npoints);
  auto points = &outSeps2.pt.points_[3 * noldpoints];

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

  outSeps2.pt.numberOfPoints_ = npoints;
  outSeps2.cl.numberOfCells_ = ncells;

  return 0;
}

// ---------------- //
//   Segmentation   //
// ---------------- //

template <typename triangulationType>
int ttk::TopologicalSkeleton::setAscendingSegmentation(
  const std::vector<SimplexId> &maxima,
  std::vector<Separatrix> &repellingOrbits,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  const auto thisDim{triangulation.getDimensionality()};

  if(morseSmaleManifold == nullptr) {
    this->printErr("Could not compute ascending segmentation");
    return 1;
  }

  Timer tm{};

  const auto nVerts{triangulation.getNumberOfVertices()};
  std::fill(morseSmaleManifold, morseSmaleManifold + nVerts, -1);
  if(maxima.empty()) {
    // shortcut for elevation
    return 0;
  }

  const auto nCells{triangulation.getNumberOfCells()};
  std::vector<SimplexId> morseSmaleManifoldOnCells(nCells, -1);

  size_t nMax{};
  for(const auto &id : maxima) {
    // mark the maxima
    morseSmaleManifoldOnCells[id] = nMax++;
  }
  for(const auto &orbit : repellingOrbits) {
    // mark the orbits
    auto cycleCells = orbit.geometry_;
    auto cycleNumber = nMax++;
    for(auto &cell : cycleCells) {
      if(cell.dim_ == thisDim) {
        morseSmaleManifoldOnCells[cell.id_] = cycleNumber;
      }
    }
  }

  const auto dim{triangulation.getDimensionality()};

  // Triangulation method pointers for 3D
  auto getFaceStarNumber = &triangulationType::getTriangleStarNumber;
  auto getFaceStar = &triangulationType::getTriangleStar;
  if(dim == 2) {
    // Triangulation method pointers for 2D
    getFaceStarNumber = &triangulationType::getEdgeStarNumber;
    getFaceStar = &triangulationType::getEdgeStar;
  } else if(dim == 1) {
    // Triangulation method pointers for 1D
    getFaceStarNumber = &triangulationType::getVertexStarNumber;
    getFaceStar = &triangulationType::getVertexStar;
  }

  // cells visited during the propagation alongside one integral line
  std::vector<SimplexId> visited{};
  // all marked cells
  std::vector<uint8_t> isMarked(nCells, 0);
  // cycle detection array
  std::vector<char> isCycle;
  isCycle.resize(nCells, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(visited, isCycle)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nCells; ++i) {
    if(isMarked[i] == 1) {
      continue;
    }
    visited.clear();
    auto curr{i};
    while(morseSmaleManifoldOnCells[curr] == -1) {
      if(isMarked[curr] == 1) {
        break;
      } else if(isCycle[curr] == 1) { // Unmarked cycle(shouldn't occur)
        morseSmaleManifoldOnCells[curr]
          = -2; // Backup behavior(handled like boundary)
        break;
      }
      // follow a V-path till an already marked cell is reached
      // Break when reaching a cycle
      const auto paired{this->simplifierField_.dcvf_.getPairedCell(
        Cell{dim, curr}, triangulation, true)};
      SimplexId next{curr};
      const auto nStars{(triangulation.*getFaceStarNumber)(paired)};
      for(SimplexId j = 0; j < nStars; ++j) {
        (triangulation.*getFaceStar)(paired, j, next);
        // get the first cell != curr (what of non-manifold datasets?)
        if(next != curr) {
          break;
        }
      }
      visited.emplace_back(curr);
      isCycle[curr] = 1;
      if(next == curr) {
        // on the boundary?
        break;
      }
      curr = next;
    }
    for(const auto el : visited) {
      morseSmaleManifoldOnCells[el] = morseSmaleManifoldOnCells[curr];
      isMarked[el] = 1;
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    if(triangulation.getVertexStarNumber(i) < 1) {
      // handle non-manifold datasets?
      continue;
    }
    SimplexId starId{};
    triangulation.getVertexStar(i, 0, starId);
    // put segmentation infos from cells to points
    morseSmaleManifold[i] = morseSmaleManifoldOnCells[starId];
  }

  this->printMsg("  Ascending segmentation computed", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::setDescendingSegmentation(
  const std::vector<SimplexId> &minima,
  std::vector<Separatrix> &attractingOrbits,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  if(morseSmaleManifold == nullptr) {
    this->printErr("Could not compute descending segmentation");
    return 1;
  }

  Timer tm{};

  const auto nVerts{triangulation.getNumberOfVertices()};

  if(minima.size() == 1) {
    // shortcut for elevation
    std::fill(morseSmaleManifold, morseSmaleManifold + nVerts, 0);
    return 0;
  }

  std::fill(morseSmaleManifold, morseSmaleManifold + nVerts, -1);

  size_t nMin{};
  for(const auto &cp : minima) {
    // mark the minima
    morseSmaleManifold[cp] = nMin++;
  }

  for(const auto &orbit : attractingOrbits) {
    auto cycleCells = orbit.geometry_;
    auto cycleNumber = nMin++;
    for(auto &cell : cycleCells) {
      if(cell.dim_ == 0) {
        morseSmaleManifold[cell.id_] = cycleNumber;
      }
    }
  }

  std::vector<SimplexId> visited{};
  std::vector<char> isCycle;
  isCycle.resize(nVerts, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(visited, isCycle)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    if(morseSmaleManifold[i] != -1) {
      continue;
    }
    visited.clear();
    auto curr{i};
    while(morseSmaleManifold[curr] == -1) {
      // Check if a cycle has occured
      if(isCycle[curr] == 1) { // Unmarked cycle (shouldn't occur)
        morseSmaleManifold[curr]
          = -2; // Backup behavior (handled like boundary)
        break;
      }
      // follow a V-path till an already marked vertex is reached
      const auto pairedEdge{this->simplifierField_.dcvf_.getPairedCell(
        Cell{0, curr}, triangulation)};
      SimplexId next{};
      triangulation.getEdgeVertex(pairedEdge, 0, next);
      if(next == curr) {
        triangulation.getEdgeVertex(pairedEdge, 1, next);
      }
      visited.emplace_back(curr);
      isCycle[curr] = 1;
      curr = next;
    }
    for(const auto el : visited) {
      morseSmaleManifold[el] = morseSmaleManifold[curr];
    }
  }

  this->printMsg("  Descending segmentation computed", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

template <typename triangulationType>
int ttk::TopologicalSkeleton::setFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  if(ascendingManifold == nullptr || descendingManifold == nullptr
     || morseSmaleManifold == nullptr) {
    this->printErr("Could not compute final segmentation");
    return 1;
  }

  Timer tm{};

  const size_t nVerts = triangulation.getNumberOfVertices();

  // associate a unique "sparse region id" to each (ascending, descending) pair

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    const auto d = ascendingManifold[i];
    const auto a = descendingManifold[i];
    if(a == -1 || d == -1 || a == -2 || d == -2) {
      morseSmaleManifold[i] = -1;
    } else {
      morseSmaleManifold[i] = a * numberOfMaxima + d;
    }
  }

  // store the "sparse region ids" by copying the morseSmaleManifold output
  std::vector<SimplexId> sparseRegionIds(
    morseSmaleManifold, morseSmaleManifold + nVerts);

  // get unique "sparse region ids"
  TTK_PSORT(
    this->threadNumber_, sparseRegionIds.begin(), sparseRegionIds.end());
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

  this->printMsg("  Final segmentation computed", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}
