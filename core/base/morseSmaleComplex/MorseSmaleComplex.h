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

namespace ttk {
  class MorseSmaleComplex : public virtual Debug {
  public:
    MorseSmaleComplex();

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
        std::vector<SimplexId> sepFuncMaxId_{};
        std::vector<SimplexId> sepFuncMinId_{};
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
        std::vector<SimplexId> sepFuncMaxId_{};
        std::vector<SimplexId> sepFuncMinId_{};
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
     * Main function for computing the Morse-Smale complex.
     *
     * @pre MorseSmaleComplex::preconditionTriangulation must be
     * called prior to this.
     */
    template <typename dataType, typename triangulationType>
    inline int execute(OutputCriticalPoints &outCP,
                       Output1Separatrices &outSeps1,
                       Output2Separatrices &outSeps2,
                       OutputManifold &outManifold,
                       const dataType *const scalars,
                       const SimplexId *const offsets,
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
     * Set the threshold for the iterative gradient reversal process.
     * Disable thresholding with -1 (default).
     */
    int setIterationThreshold(const int iterationThreshold) {
      discreteGradient_.setIterationThreshold(iterationThreshold);
      return 0;
    }

    /**
     * Enable/Disable post-processing gradient reversal of
     * the (saddle,...,saddle) vpaths under a given persistence
     * threshold (disabled by default).
     */
    int setReturnSaddleConnectors(const bool state) {
      ReturnSaddleConnectors = state;
      discreteGradient_.setReturnSaddleConnectors(state);
      return 0;
    }

    /**
     * Set the threshold value for post-processing of
     * (saddle,...,saddle) vpaths gradient reversal
     * (default value is 0.0).
     */
    int setSaddleConnectorsPersistenceThreshold(const double threshold) {
      SaddleConnectorsPersistenceThreshold = threshold;
      discreteGradient_.setSaddleConnectorsPersistenceThreshold(threshold);
      return 0;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      this->discreteGradient_.preconditionTriangulation(data);
      data->preconditionCellEdges();
      data->preconditionCellNeighbors();
    }

  protected:
    /**
     * Utility class representing Ridge lines, Valley lines
     * and Saddle connectors.
     */
    struct Separatrix {
      explicit Separatrix() = default;

      // initialization with one segment
      explicit Separatrix(const bool isValid,
                          const dcg::Cell &saddle,
                          const dcg::Cell &extremum,
                          const bool isSegmentReversed,
                          const SimplexId segmentGeometry)
        : isValid_{isValid}, source_{saddle}, destination_{extremum} {
        isReversed_.push_back(isSegmentReversed);
        geometry_.push_back(segmentGeometry);
      }

      /** Flag indicating if this separatrix can be processed. */
      bool isValid_;
      /** Source cell of the separatrix. */
      dcg::Cell source_;
      /** Destination cell of the separatrix. */
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
     * Compute the geometrical embedding of the 1-separatrices.
     */
    template <typename triangulationType>
    int setSeparatrices1(
      Output1Separatrices &outSeps1,
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const SimplexId *const offsets,
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
     * 2-separatrices.
     */
    template <typename triangulationType>
    int setDescendingSeparatrices2(
      Output2Separatrices &outSeps2,
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const SimplexId *const offsets,
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
     * 2-separatrices.
     */
    template <typename triangulationType>
    int setAscendingSeparatrices2(
      Output2Separatrices &outSeps2,
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const SimplexId *const offsets,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<Separatrix>> &separatrices,
      std::vector<std::vector<std::vector<ttk::dcg::Cell>>>
        &separatricesGeometry) const;

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

  protected:
    dcg::DiscreteGradient discreteGradient_{};

    bool ComputeCriticalPoints{true};
    bool ComputeAscendingSeparatrices1{true};
    bool ComputeDescendingSeparatrices1{true};
    bool ComputeSaddleConnectors{true};
    bool ComputeAscendingSeparatrices2{false};
    bool ComputeDescendingSeparatrices2{false};
    bool ComputeAscendingSegmentation{true};
    bool ComputeDescendingSegmentation{true};
    bool ComputeFinalSegmentation{true};

    bool ReturnSaddleConnectors{false};
    double SaddleConnectorsPersistenceThreshold{};
  };
} // namespace ttk

// ---------------- //
//  Execute method  //
// ---------------- //

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplex::execute(OutputCriticalPoints &outCP,
                                    Output1Separatrices &outSeps1,
                                    Output2Separatrices &outSeps2,
                                    OutputManifold &outManifold,
                                    const dataType *const scalars,
                                    const SimplexId *const offsets,
                                    const triangulationType &triangulation) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(scalars == nullptr) {
    this->printErr("Input scalar field pointer is null.");
    return -1;
  }

  if(offsets == nullptr) {
    this->printErr("Input offset field pointer is null.");
    return -1;
  }
#endif
  Timer t;

  outCP.clear();
  outSeps1.clear();
  outSeps2.clear();
  const auto dim = triangulation.getDimensionality();

  this->discreteGradient_.setThreadNumber(threadNumber_);
  this->discreteGradient_.setDebugLevel(debugLevel_);
  this->discreteGradient_.setInputScalarField(scalars);
  this->discreteGradient_.setInputOffsets(offsets);
  this->discreteGradient_.buildGradient(triangulation);

  if(dim == 3 && ReturnSaddleConnectors) {
    discreteGradient_.reverseGradient<dataType>(triangulation);
  }

  std::vector<dcg::Cell> criticalPoints{};
  discreteGradient_.getCriticalPoints(criticalPoints, triangulation);

  std::vector<std::vector<Separatrix>> separatrices1{};
  std::vector<std::vector<std::vector<dcg::Cell>>> separatricesGeometry1;

  // 1-separatrices
  Timer tm1sep{};

  if(dim > 1 && ComputeDescendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getDescendingSeparatrices1(criticalPoints, separatrices1.back(),
                               separatricesGeometry1.back(), triangulation);

    this->printMsg("  Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(dim > 1 && ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getAscendingSeparatrices1(criticalPoints, separatrices1.back(),
                              separatricesGeometry1.back(), triangulation);

    this->printMsg("  Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  // saddle-connectors
  if(dim == 3 && ComputeSaddleConnectors) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getSaddleConnectors(criticalPoints, separatrices1.back(),
                        separatricesGeometry1.back(), triangulation);

    this->printMsg("  Saddle connectors computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_, debug::LineMode::NEW,
                   debug::Priority::DETAIL);
  }

  if(dim > 1
     && (ComputeDescendingSeparatrices1 || ComputeAscendingSeparatrices1
         || ComputeSaddleConnectors)) {
    Timer tmp{};

    flattenSeparatricesVectors(separatrices1, separatricesGeometry1);
    setSeparatrices1(outSeps1, separatrices1[0], separatricesGeometry1[0],
                     offsets, triangulation);

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
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices,
                               separatricesGeometry, separatricesSaddles,
                               triangulation);
    setDescendingSeparatrices2(outSeps2, separatrices, separatricesGeometry,
                               separatricesSaddles, offsets, triangulation);

    this->printMsg("  Descending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(dim == 3 && ComputeAscendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices,
                              separatricesGeometry, separatricesSaddles,
                              triangulation);
    setAscendingSeparatrices2(outSeps2, separatrices, separatricesGeometry,
                              separatricesSaddles, offsets, triangulation);

    this->printMsg("  Ascending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  this->printMsg("2-separatrices computed", 1.0, tm2sep.getElapsedTime(),
                 this->threadNumber_);

  if(ComputeAscendingSegmentation || ComputeDescendingSegmentation) {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ComputeAscendingSegmentation) {
      std::vector<SimplexId> maxSeeds{};
      setAscendingSegmentation(criticalPoints, maxSeeds, outManifold.ascending_,
                               numberOfMaxima, triangulation);
    }
    if(ComputeDescendingSegmentation) {
      setDescendingSegmentation(
        criticalPoints, outManifold.descending_, numberOfMinima, triangulation);
    }
    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation
       && ComputeFinalSegmentation) {
      setFinalSegmentation(numberOfMaxima, numberOfMinima,
                           outManifold.ascending_, outManifold.descending_,
                           outManifold.morseSmale_, triangulation);
    }

    this->printMsg(
      "Segmentation computed", 1.0, tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeCriticalPoints) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints(
      criticalPoints, nCriticalPointsByDim, outCP.points_,
      outCP.cellDimensions_, outCP.cellIds_, outCP.isOnBoundary_,
      outCP.PLVertexIdentifiers_, triangulation);

    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation) {
      discreteGradient_.setManifoldSize(
        criticalPoints.size(), nCriticalPointsByDim, outManifold.ascending_,
        outManifold.descending_, outCP.manifoldSize_);
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
int ttk::MorseSmaleComplex::setSeparatrices1(
  Output1Separatrices &outSeps1,
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  auto &separatrixFunctionMaxima = outSeps1.cl.sepFuncMaxId_;
  auto &separatrixFunctionMinima = outSeps1.cl.sepFuncMinId_;

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
  separatrixFunctionMaxima.resize(separatrixId + validGeomIds.size());
  separatrixFunctionMinima.resize(separatrixId + validGeomIds.size());
  outSeps1.cl.isOnBoundary_.resize(ncells);

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
    const auto vertsOrder = [offsets](const SimplexId a, const SimplexId b) {
      return offsets[a] < offsets[b];
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
    separatrixFunctionMaxima[sepId] = sepFuncMax;
    separatrixFunctionMinima[sepId] = sepFuncMin;

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
  Output2Separatrices &outSeps2,
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  auto &separatrixFunctionMaxima = outSeps2.cl.sepFuncMaxId_;
  auto &separatrixFunctionMinima = outSeps2.cl.sepFuncMinId_;

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
  separatrixFunctionMaxima.resize(separatrixId + validGeomIds.size());
  separatrixFunctionMinima.resize(separatrixId + validGeomIds.size());
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
    separatrixFunctionMaxima[sepId] = sepFuncMax;
    separatrixFunctionMinima[sepId] = sepFuncMin;
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
int ttk::MorseSmaleComplex::setDescendingSeparatrices2(
  Output2Separatrices &outSeps2,
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  auto &separatrixFunctionMaxima = outSeps2.cl.sepFuncMaxId_;
  auto &separatrixFunctionMinima = outSeps2.cl.sepFuncMinId_;

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
  outSeps2.cl.offsets_.resize(ncells + 1);
  outSeps2.cl.offsets_[0] = 0;
  outSeps2.cl.connectivity_.resize(
    firstCellId + 3 * (ncells - noldcells)); // triangles cells
  auto cellsOff = &outSeps2.cl.offsets_[noldcells];
  auto cellsConn = &outSeps2.cl.connectivity_[firstCellId];
  outSeps2.cl.sourceIds_.resize(ncells);
  outSeps2.cl.separatrixIds_.resize(ncells);
  outSeps2.cl.separatrixTypes_.resize(ncells);
  separatrixFunctionMaxima.resize(separatrixId + validGeomIds.size());
  separatrixFunctionMinima.resize(separatrixId + validGeomIds.size());
  outSeps2.cl.isOnBoundary_.resize(ncells);

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
    separatrixFunctionMaxima[sepId] = sepFuncMax;
    separatrixFunctionMinima[sepId] = sepFuncMin;

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
int ttk::MorseSmaleComplex::setAscendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  std::vector<SimplexId> &maxSeeds,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation) const {

  if(morseSmaleManifold == nullptr) {
    this->printErr("Could not compute ascending segmentation");
    return 1;
  }

  Timer tm{};

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
  } else if(cellDim == 1) {
    // Triangulation method pointers for 1D
    getCellFace = &triangulationType::getEdgeVertex;
    getFaceStarNumber = &triangulationType::getVertexStarNumber;
    getFaceStar = &triangulationType::getVertexStar;
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

  this->printMsg("  Ascending segmentation computed", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setDescendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMinima,
  const triangulationType &triangulation) const {

  if(morseSmaleManifold == nullptr) {
    this->printErr("Could not compute descending segmentation");
    return 1;
  }

  Timer tm{};

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

  this->printMsg("  Descending segmentation computed", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

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
