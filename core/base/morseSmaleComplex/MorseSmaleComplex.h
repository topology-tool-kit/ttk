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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions/">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseMolecule/">Morse
///   molecule example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 0 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering1/">Persistence
///   clustering 1 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering2/">Persistence
///   clustering 2 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering3/">Persistence
///   clustering 3 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering4/">Persistence
///   clustering 4 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_at/">Persistent
///   Generators AT example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_darkSky/">Persistent
///   Generators DarkSky example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tribute/">Tribute
///   example</a> \n
///

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <DiscreteMorseSandwich.h>
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
                       const size_t scalarsMTime,
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
     * Enable/Disable post-processing gradient reversal of
     * the (saddle,...,saddle) vpaths under a given persistence
     * threshold (disabled by default).
     */
    inline void setReturnSaddleConnectors(const bool state) {
      ReturnSaddleConnectors = state;
    }

    /**
     * Set the threshold value for post-processing of
     * (saddle,...,saddle) vpaths gradient reversal
     * (default value is 0.0).
     */
    inline void
      setSaddleConnectorsPersistenceThreshold(const double threshold) {
      SaddleConnectorsPersistenceThreshold = threshold;
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
     * gradient.
     */
    template <typename triangulationType>
    int
      getDescendingSeparatrices1(const std::vector<SimplexId> &saddles,
                                 std::vector<Separatrix> &separatrices,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the ascending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(const std::vector<SimplexId> &saddles,
                                  std::vector<Separatrix> &separatrices,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getSaddleConnectors(const std::vector<SimplexId> &saddles2,
                            std::vector<Separatrix> &separatrices,
                            const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the 1-separatrices.
     */
    template <typename triangulationType>
    int setSeparatrices1(Output1Separatrices &outSeps1,
                         const std::vector<Separatrix> &separatrices,
                         const SimplexId *const offsets,
                         const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the maxima.
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
      const SimplexId *const offsets,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<Separatrix>> &separatrices) const;

    /**
     * Compute the ascending manifold of the maxima.
     */
    template <typename triangulationType>
    int setAscendingSegmentation(const std::vector<SimplexId> &maxima,
                                 SimplexId *const morseSmaleManifold,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the descending manifold of the minima.
     */
    template <typename triangulationType>
    int setDescendingSegmentation(const std::vector<SimplexId> &minima,
                                  SimplexId *const morseSmaleManifold,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the final combinatorial Morse-Smale complex
     * segmentation.
     */
    template <typename triangulationType>
    int setFinalSegmentation(const SimplexId numberOfMaxima,
                             const SimplexId *const ascendingManifold,
                             const SimplexId *const descendingManifold,
                             SimplexId *const morseSmaleManifold,
                             const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int returnSaddleConnectors(const double persistenceThreshold,
                               const dataType *const scalars,
                               const SimplexId *const offsets,
                               const triangulationType &triangulation);

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
    bool ThresholdIsAbsolute{false};
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
                                    const size_t scalarsMTime,
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
  this->discreteGradient_.setInputScalarField(scalars, scalarsMTime);
  this->discreteGradient_.setInputOffsets(offsets);
  this->discreteGradient_.buildGradient(
    triangulation, this->ReturnSaddleConnectors);

  if(this->ReturnSaddleConnectors) {
    auto persistenceThreshold{this->SaddleConnectorsPersistenceThreshold};
    if(!this->ThresholdIsAbsolute) {
      const auto nVerts{triangulation.getNumberOfVertices()};
      // global extrema are (generally) faster computed on offsets
      // than on scalar field
      const auto pair{std::minmax_element(offsets, offsets + nVerts)};
      // global extrema vertex ids
      const auto globmin = std::distance(offsets, pair.first);
      const auto globmax = std::distance(offsets, pair.second);
      persistenceThreshold *= (scalars[globmax] - scalars[globmin]);
      this->printMsg("Absolute saddle connectors persistence threshold is "
                       + std::to_string(persistenceThreshold),
                     debug::Priority::DETAIL);
    }

    this->returnSaddleConnectors(
      persistenceThreshold, scalars, offsets, triangulation);
  }

  std::array<std::vector<SimplexId>, 4> criticalPoints{};
  {
    Timer tm{};
    discreteGradient_.getCriticalPoints(criticalPoints, triangulation);
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

    getDescendingSeparatrices1(
      criticalPoints[1], separatrices1.back(), triangulation);

    this->printMsg("  Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  if(dim > 1 && ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();

    getAscendingSeparatrices1(
      criticalPoints[dim - 1], separatrices1.back(), triangulation);

    this->printMsg("  Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
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
         || ComputeSaddleConnectors)) {
    Timer tmp{};

    this->flattenSeparatricesVectors(separatrices1);
    setSeparatrices1(outSeps1, separatrices1[0], offsets, triangulation);

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
      outSeps2, separatrices, separatricesSaddles, offsets, triangulation);

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
      outSeps2, separatrices, separatricesSaddles, offsets, triangulation);

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

    if(ComputeAscendingSegmentation) {
      setAscendingSegmentation(
        criticalPoints[dim], outManifold.ascending_, triangulation);
    }
    if(ComputeDescendingSegmentation) {
      setDescendingSegmentation(
        criticalPoints[0], outManifold.descending_, triangulation);
    }
    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation
       && ComputeFinalSegmentation) {
      setFinalSegmentation(criticalPoints[dim].size(), outManifold.ascending_,
                           outManifold.descending_, outManifold.morseSmale_,
                           triangulation);
    }

    this->printMsg(
      "Segmentation computed", 1.0, tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeCriticalPoints) {
    discreteGradient_.setCriticalPoints(
      criticalPoints, outCP.points_, outCP.cellDimensions_, outCP.cellIds_,
      outCP.isOnBoundary_, outCP.PLVertexIdentifiers_, triangulation);

    if(ComputeAscendingSegmentation && ComputeDescendingSegmentation) {
      discreteGradient_.setManifoldSize(criticalPoints, outManifold.ascending_,
                                        outManifold.descending_,
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

template <typename triangulationType>
int ttk::MorseSmaleComplex::getDescendingSeparatrices1(
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
        vpath.push_back(saddle1);
        discreteGradient_.getDescendingPath(
          Cell(0, vertexId), vpath, triangulation);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 0 and discreteGradient_.isCellCritical(lastCell)) {
          separatrices[2 * i + j].source_ = saddle;
          separatrices[2 * i + j].destination_ = lastCell;
          separatrices[2 * i + j].geometry_ = std::move(vpath);
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getAscendingSeparatrices1(
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
      discreteGradient_.getAscendingPath(Cell(dim, sId), vpath, triangulation);

      const Cell &lastCell = vpath.back();
      if(lastCell.dim_ == dim and discreteGradient_.isCellCritical(lastCell)) {
        sepsPerSaddle[i].emplace_back();
        sepsPerSaddle[i].back().source_ = saddle;
        sepsPerSaddle[i].back().destination_ = lastCell;
        sepsPerSaddle[i].back().geometry_ = std::move(vpath);
      }
    }
  }

  this->flattenSeparatricesVectors(sepsPerSaddle);

  separatrices = std::move(sepsPerSaddle[0]);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getSaddleConnectors(
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
    discreteGradient_.getDescendingWall(
      s2, mask, triangulation, nullptr, &saddles1);

    for(const auto saddle1Id : saddles1) {
      const Cell s1{1, saddle1Id};

      Vpath vpath;
      const bool isMultiConnected
        = discreteGradient_.getAscendingPathThroughWall(
          s1, s2, isVisited, &vpath, triangulation);

      if(vpath.empty()) {
        // safety, should be unreachable
        continue;
      }
      const auto &last = vpath.back();

      if(!isMultiConnected && last.dim_ == s2.dim_ && last.id_ == s2.id_) {
        sepsByThread[i].emplace_back();
        sepsByThread[i].back().source_ = s1;
        sepsByThread[i].back().destination_ = s2;
        sepsByThread[i].back().geometry_ = std::move(vpath);
      }
    }
  }

  this->flattenSeparatricesVectors(sepsByThread);

  separatrices = std::move(sepsByThread[0]);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::setSeparatrices1(
  Output1Separatrices &outSeps1,
  const std::vector<Separatrix> &separatrices,
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
  separatrixFunctionMaxima.resize(separatrixId + separatrices.size());
  separatrixFunctionMinima.resize(separatrixId + separatrices.size());
  outSeps1.cl.isOnBoundary_.resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    const auto &sepGeom = sep.geometry_;
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
    discreteGradient_.getAscendingWall(
      saddle1, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatrices[i].source_ = saddle1;
    separatrices[i].destination_ = emptyCell;
    separatrices[i].geometry_ = std::move(wall);
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex::getDescendingSeparatrices2(
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
    discreteGradient_.getDescendingWall(
      saddle2, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatrices[i].source_ = saddle2;
    separatrices[i].destination_ = emptyCell;
    separatrices[i].geometry_ = std::move(wall);
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
  const std::vector<std::vector<SimplexId>> &separatricesSaddles,
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
  separatrixFunctionMaxima.resize(separatrixId + separatrices.size());
  separatrixFunctionMinima.resize(separatrixId + separatrices.size());
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

    // compute separatrix function diff
    const auto sepFuncMin
      = discreteGradient_.getCellLowerVertex(src, triangulation);
    SimplexId sepFuncMax{};
    if(!sepSaddles.empty()) {
      // find minimum vertex on the critical triangles of the 2-separatrix
      const auto maxId = *std::max_element(
        sepSaddles.begin(), sepSaddles.end(),
        [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
          return offsets[discreteGradient_.getCellGreaterVertex(
                   Cell{2, a}, triangulation)]
                 < offsets[discreteGradient_.getCellGreaterVertex(
                   Cell{2, b}, triangulation)];
        });
      sepFuncMax
        = discreteGradient_.getCellGreaterVertex(Cell{2, maxId}, triangulation);
    } else {
      // find maximum vertex by iterating over all the edges in the
      // 2-separatrix
      const auto maxId = *std::max_element(
        sepGeom.begin(), sepGeom.end(),
        [&triangulation, offsets, this](const Cell &a, const Cell &b) {
          return offsets[discreteGradient_.getCellGreaterVertex(
                   a, triangulation)]
                 < offsets[discreteGradient_.getCellGreaterVertex(
                   b, triangulation)];
        });
      sepFuncMax = discreteGradient_.getCellGreaterVertex(maxId, triangulation);
    }

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
  const std::vector<std::vector<SimplexId>> &separatricesSaddles,
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
  separatrixFunctionMaxima.resize(separatrixId + separatrices.size());
  separatrixFunctionMinima.resize(separatrixId + separatrices.size());
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

    // compute separatrix function diff
    const auto sepFuncMax
      = discreteGradient_.getCellGreaterVertex(src, triangulation);
    SimplexId sepFuncMin{};
    if(!sepSaddles.empty()) {
      // find minimum vertex on the critical edges of the 2-separatrix
      const auto minId = *std::min_element(
        sepSaddles.begin(), sepSaddles.end(),
        [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
          return offsets[discreteGradient_.getCellLowerVertex(
                   Cell{1, a}, triangulation)]
                 < offsets[discreteGradient_.getCellLowerVertex(
                   Cell{1, b}, triangulation)];
        });
      sepFuncMin
        = discreteGradient_.getCellLowerVertex(Cell{1, minId}, triangulation);
    } else {
      // find minimum vertex by iterating over all the triangles in the
      // 2-separatrix
      const auto minId = *std::min_element(
        sepGeom.begin(), sepGeom.end(),
        [&triangulation, offsets, this](const Cell &a, const Cell &b) {
          return offsets[discreteGradient_.getCellLowerVertex(a, triangulation)]
                 < offsets[discreteGradient_.getCellLowerVertex(
                   b, triangulation)];
        });
      sepFuncMin = discreteGradient_.getCellLowerVertex(minId, triangulation);
    }
    separatrixFunctionMaxima[sepId] = sepFuncMax;
    separatrixFunctionMinima[sepId] = sepFuncMin;

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
int ttk::MorseSmaleComplex::setAscendingSegmentation(
  const std::vector<SimplexId> &maxima,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

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

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(visited)
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
      }
      // follow a V-path till an already marked cell is reached
      const auto paired{this->discreteGradient_.getPairedCell(
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
int ttk::MorseSmaleComplex::setDescendingSegmentation(
  const std::vector<SimplexId> &minima,
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

  std::vector<SimplexId> visited{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) firstprivate(visited)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    if(morseSmaleManifold[i] != -1) {
      continue;
    }
    visited.clear();
    auto curr{i};
    while(morseSmaleManifold[curr] == -1) {
      // follow a V-path till an already marked vertex is reached
      const auto pairedEdge{
        this->discreteGradient_.getPairedCell(Cell{0, curr}, triangulation)};
      SimplexId next{};
      triangulation.getEdgeVertex(pairedEdge, 0, next);
      if(next == curr) {
        triangulation.getEdgeVertex(pairedEdge, 1, next);
      }
      visited.emplace_back(curr);
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
int ttk::MorseSmaleComplex::setFinalSegmentation(
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

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplex::returnSaddleConnectors(
  const double persistenceThreshold,
  const dataType *const scalars,
  const SimplexId *const offsets,
  const triangulationType &triangulation) {

  Timer tm{};

  const auto dim{triangulation.getDimensionality()};
  if(dim != 3) {
    this->printWrn("Can't return saddle connectors without a 3D dataset");
    return 0;
  }

  // compute saddle-saddle Persistence Pairs from DiscreteMorseSandwich
  ttk::DiscreteMorseSandwich dms{};
  dms.setThreadNumber(this->threadNumber_);
  dms.setDebugLevel(this->debugLevel_);
  dms.setGradient(std::move(this->discreteGradient_));

  using PersPairType = DiscreteMorseSandwich::PersistencePair;

  std::vector<PersPairType> dms_pairs{};
  dms.computePersistencePairs(dms_pairs, offsets, triangulation, false, true);
  this->discreteGradient_ = dms.getGradient();
  // reset gradient pointer to local storage
  this->discreteGradient_.setLocalGradient();

  const auto getPersistence
    = [this, &triangulation, scalars](const PersPairType &p) {
        return this->discreteGradient_.getPersistence(
          Cell{2, p.death}, Cell{1, p.birth}, scalars, triangulation);
      };

  // saddle-saddle pairs should be in one continuous block inside dms_pairs
  const auto firstSadSadPair{std::distance(
    dms_pairs.begin(),
    std::find_if(dms_pairs.begin(), dms_pairs.end(),
                 [](const auto &pair) { return pair.type == 1; }))};

  // DMS pairs output are already sorted by 2-saddles in ascending order

  std::vector<bool> isVisited(triangulation.getNumberOfTriangles(), false);
  std::vector<SimplexId> visitedTriangles{};

  size_t nReturned{};

  // offset from firstSadSadPair should index s2Children and simplifyS2
  const auto &s2Children{dms.get2SaddlesChildren()};
  std::vector<bool> simplifyS2(s2Children.size(), true);

  for(size_t i = firstSadSadPair; i < dms_pairs.size(); ++i) {
    const auto &pair{dms_pairs[i]};

    if(pair.type != 1 || getPersistence(pair) > persistenceThreshold) {
      continue;
    }

    const auto o{i - firstSadSadPair};
    const Cell birth{1, pair.birth};
    const Cell death{2, pair.death};
    bool skip{false};
    for(const auto child : s2Children[o]) {
      if(!simplifyS2[child]) {
        skip = true;
        break;
      }
    }
    if(skip) {
      this->printMsg("Skipping saddle connector " + birth.to_string() + " -> "
                       + death.to_string(),
                     debug::Priority::DETAIL);
      simplifyS2[o] = false;
      continue;
    }
    // 1. get the 2-saddle wall
    VisitedMask mask{isVisited, visitedTriangles};
    this->discreteGradient_.getDescendingWall(death, mask, triangulation);
    // 2. get the saddle connector
    std::vector<Cell> vpath{};
    this->discreteGradient_.getAscendingPathThroughWall(
      birth, death, isVisited, &vpath, triangulation, true);
    // 3. reverse the gradient on the saddle connector path
    if(vpath.back() == death) {
      this->discreteGradient_.reverseAscendingPathOnWall(vpath, triangulation);
      nReturned++;
    } else {
      this->printMsg("Could not return saddle connector " + birth.to_string()
                       + " -> " + death.to_string(),
                     debug::Priority::DETAIL);
      simplifyS2[o] = false;
    }
  }

  this->printMsg("Returned " + std::to_string(nReturned) + " saddle connectors",
                 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}
