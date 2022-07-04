/// \ingroup baseCode
/// \class ttk::DiscreteMorseSandwich
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date January 2021.
///
/// \brief TTK %DiscreteMorseSandwich processing package.
///
/// %DiscreteMorseSandwich computes a Persistence Diagram by using the
/// %Discrete Morse-Theory %DiscreteGradient algorithms.
///
/// \b Related \b publication \n
/// "Discrete Morse Sandwich: Fast Computation of Persistence Diagrams for
/// Scalar Data -- An Algorithm and A Benchmark" \n
/// Pierre Guillou, Jules Vidal, Julien Tierny \n
/// Technical Report, arXiv:2206.13932, 2022
///
///
/// \sa ttk::dcg::DiscreteGradient

#pragma once

#include <DiscreteGradient.h>

#include <numeric>

namespace ttk {
  class DiscreteMorseSandwich : virtual public Debug {
  public:
    DiscreteMorseSandwich();

    /**
     * @brief Persistence pair struct as exported by DiscreteGradient
     */
    struct PersistencePair {
      /** first (lower/birth) simplex cell id */
      SimplexId birth;
      /** second (higher/death) simplex cell id */
      SimplexId death;
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      int type;

      PersistencePair(SimplexId b, SimplexId d, int t)
        : birth{b}, death{d}, type{t} {
      }
    };

    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      this->dg_.preconditionTriangulation(data);
    }

    inline void setInputOffsets(const SimplexId *const offsets) {
      this->dg_.setInputOffsets(offsets);
    }

    template <typename triangulationType>
    inline int buildGradient(const void *const scalars,
                             const size_t scalarsMTime,
                             const SimplexId *const offsets,
                             const triangulationType &triangulation) {
      this->dg_.setDebugLevel(this->debugLevel_);
      this->dg_.setThreadNumber(this->threadNumber_);
      this->dg_.setInputOffsets(offsets);
      this->dg_.setInputScalarField(scalars, scalarsMTime);
      return this->dg_.buildGradient(triangulation);
    }

    template <typename triangulationType>
    inline SimplexId
      getCellGreaterVertex(const dcg::Cell &c,
                           const triangulationType &triangulation) {
      return this->dg_.getCellGreaterVertex(c, triangulation);
    }

    /**
     * @brief Compute the persistence pairs from the discrete gradient
     *
     * @pre @ref buildGradient and @ref preconditionTriangulation
     * should be called prior to this function
     *
     * @param[out] pairs Output persistence pairs
     * @param[in] offsets Order field
     * @param[in] triangulation Preconditionned triangulation
     * @param[in] ignoreBoundary Ignore the boundary component
     *
     * @return 0 when success
     */
    template <typename triangulationType>
    int computePersistencePairs(std::vector<PersistencePair> &pairs,
                                const SimplexId *const offsets,
                                const triangulationType &triangulation,
                                const bool ignoreBoundary);

    /**
     * @brief Type for exporting persistent generators
     *
     * A generator = a 2-saddle index + vector of edges with 1-saddle
     * at index 0.
     */
    struct GeneratorType {
      /** Generator edges beginning with the 1-saddle */
      std::vector<SimplexId> boundary;
      /** Critical triangle index (-1 if infinite) */
      SimplexId critTriangleId;
      /** Vertex indices for the critical triangle (or global max) and
          the critical edge */
      std::array<SimplexId, 2> critVertsIds;
    };

  protected:
    /**
     * @brief Follow the descending 1-separatrices to compute the saddles ->
     * minima association
     *
     * @param[in] criticalEdges Critical edges identifiers
     * @param[in] triangulation Triangulation
     *
     * @return a vector of minima per 1-saddle
     */
    template <typename triangulationType>
    std::vector<std::vector<SimplexId>>
      getSaddle1ToMinima(const std::vector<SimplexId> &criticalEdges,
                         const triangulationType &triangulation) const;

    /**
     * @brief Follow the ascending 1-separatrices to compute the saddles ->
     * maxima association
     *
     * @param[in] criticalCells Critical cells identifiers
     * @param[in] getFaceStar Either getEdgeStar (in 2D) or getTriangleStar
     * (in 3D)
     * @param[in] getFaceStarNumber Either getEdgeStarNumber (in 2D) or
     * getTriangleStarNumber (in 3D)
     * @param[in] isOnBoundary Either isEdgeOnBoundary (in 2D) or
     * isTriangleOnBoundary (in 3D)
     * @param[in] triangulation Triangulation
     * @param[in] ignoreBoundary FTM compatibility (ignore saddles on
     * boundary)
     *
     * @return a vector of maxima per 2-saddle
     */
    template <typename triangulationType,
              typename GFS,
              typename GFSN,
              typename OB>
    std::vector<std::vector<SimplexId>>
      getSaddle2ToMaxima(const std::vector<SimplexId> &criticalCells,
                         const GFS &getFaceStar,
                         const GFSN &getFaceStarNumber,
                         const OB &isOnBoundary,
                         const triangulationType &triangulation,
                         const bool ignoreBoundary) const;

    /**
     * @brief Compute the pairs of dimension 0
     *
     * @param[out] pairs Output persistence pairs
     * @param[in] pairedMinima If minima are paired
     * @param[in] paired1Saddles If 1-saddles (or maxima in 1D) are paired
     * @param[in] criticalEdges List of 1-saddles (or maxima in 1D)
     * @param[in] critEdgesOrder Filtration order on critical edges
     * @param[in] offsets Vertex offset field
     * @param[in] triangulation Triangulation
     */
    template <typename triangulationType>
    void getMinSaddlePairs(std::vector<PersistencePair> &pairs,
                           std::vector<bool> &pairedMinima,
                           std::vector<bool> &paired1Saddles,
                           const std::vector<SimplexId> &criticalEdges,
                           const std::vector<SimplexId> &critEdgesOrder,
                           const SimplexId *const offsets,
                           const triangulationType &triangulation) const;

    /**
     * @brief Compute the pairs of dimension dim - 1
     *
     * @param[out] pairs Output persistence pairs
     * @param[in] pairedMaxima If maxima are paired
     * @param[in] pairedSaddles If 2-saddles (or 1-saddles in 2D) are paired
     * @param[in] criticalSaddles List of 2-saddles (or 1-saddles in 2D)
     * @param[in] critSaddlesOrder Filtration order on critical saddles
     * @param[in] critMaxsOrder Filtration order on maxima
     * @param[in] triangulation Triangulation
     * @param[in] ignoreBoundary Ignore the boundary component
     */
    template <typename triangulationType>
    void getMaxSaddlePairs(std::vector<PersistencePair> &pairs,
                           std::vector<bool> &pairedMaxima,
                           std::vector<bool> &pairedSaddles,
                           const std::vector<SimplexId> &criticalSaddles,
                           const std::vector<SimplexId> &critSaddlesOrder,
                           const std::vector<SimplexId> &critMaxsOrder,
                           const triangulationType &triangulation,
                           const bool ignoreBoundary) const;

    /**
     * @brief Compute the saddle-saddle pairs (in 3D)
     *
     * @param[out] pairs Output persistence pairs
     * @param[in] paired1Saddles If 1-saddles are paired
     * @param[in] paired2Saddles If 2-saddles are paired
     * @param[in] exportBoundaries If 2-saddles boundaries must be exported
     * @param[out] boundaries Vector of 2-saddles boundaries
     * @param[in] critical1Saddles Full list of 1-saddles
     * @param[in] critical2Saddles Full list of 2-saddles
     * @param[in] crit1SaddlesOrder Filtration order on 1-saddles
     * @param[in] triangulation Triangulation
     */
    template <typename triangulationType>
    void getSaddleSaddlePairs(std::vector<PersistencePair> &pairs,
                              std::vector<bool> &paired1Saddles,
                              std::vector<bool> &paired2Saddles,
                              const bool exportBoundaries,
                              std::vector<GeneratorType> &boundaries,
                              const std::vector<SimplexId> &critical1Saddles,
                              const std::vector<SimplexId> &critical2Saddles,
                              const std::vector<SimplexId> &crit1SaddlesOrder,
                              const triangulationType &triangulation) const;

    /**
     * @brief Extract & sort critical cell from the DiscreteGradient
     *
     * @param[out] criticalCellsByDim Store critical cells ids per dimension
     * @param[out] critCellsOrder Filtration order on critical cells
     * @param[in] offsets Vertex offset field
     * @param[in] triangulation Triangulation
     * @param[in] sortEdges Sort all edges vs. only 1-saddles
     */
    template <typename triangulationType>
    void extractCriticalCells(
      std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
      std::array<std::vector<SimplexId>, 4> &critCellsOrder,
      const SimplexId *const offsets,
      const triangulationType &triangulation,
      const bool sortEdges) const;

    /**
     * @brief Print number of pairs, critical cells per dimension & unpaired
     * cells
     *
     * @param[in] pairs Computed persistence pairs
     * @param[in] criticalCellsByDim Store critical cells ids per dimension
     * @param[in] pairedMinima If minima are paired
     * @param[in] paired1Saddles If 1-saddles are paired
     * @param[in] paired2Saddles If 2-saddles are paired
     * @param[in] pairedMaxima If maxima are paired
     */
    void displayStats(
      const std::vector<PersistencePair> &pairs,
      const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
      const std::vector<bool> &pairedMinima,
      const std::vector<bool> &paired1Saddles,
      const std::vector<bool> &paired2Saddles,
      const std::vector<bool> &pairedMaxima) const;

    /**
     * @brief Triplet type for persistence pairs
     *
     * [0]: saddle cell id
     * [1]: extremum 1 cell id
     * [2]: extremum 2 cell id
     */
    using tripletType = std::array<SimplexId, 3>;

    /**
     * @brief Compute persistence pairs from triplets
     *
     * @param[out] pairs Store generated persistence pairs
     * @param[in,out] pairedExtrema If critical extrema are paired
     * @param[in,out] pairedSaddles If critical saddles are paired
     * @param[in,out] reps Extrema representatives
     * @param[in] triplets Input triplets (saddle, extremum, extremum)
     * @param[in] saddlesOrder Order on saddles
     * @param[in] extremaOrder Order on extrema
     * @param[in] pairDim Pair birth simplex dimension
     */
    void tripletsToPersistencePairs(std::vector<PersistencePair> &pairs,
                                    std::vector<bool> &pairedExtrema,
                                    std::vector<bool> &pairedSaddles,
                                    std::vector<SimplexId> &reps,
                                    std::vector<tripletType> &triplets,
                                    const SimplexId *const saddlesOrder,
                                    const SimplexId *const extremaOrder,
                                    const SimplexId pairDim) const;

    /**
     * @brief Detect 1-saddles paired to a given 2-saddle
     *
     * Adapted version of ttk::PersistentSimplexPairs::eliminateBoundaries()
     *
     * @param[in] s2 Input 2-saddle (critical triangle)
     * @param[in,out] onBoundary Propagation mask
     * @param[in,out] s2Boundaries Boundaries storage (compact)
     * @param[in] s2Mapping From (critical) triangle id to compact id in
     * s2Boundaries
     * @param[in] partners Get 2-saddles paired to 1-saddles on boundary
     * @param[in] triangulation Simplicial complex
     *
     * @return Identifier of paired 1-saddle or -1
     */
    template <typename triangulationType, typename Container>
    SimplexId
      eliminateBoundariesSandwich(const SimplexId s2,
                                  std::vector<bool> &onBoundary,
                                  std::vector<Container> &s2Boundaries,
                                  const std::vector<SimplexId> &s2Mapping,
                                  const std::vector<SimplexId> &partners,
                                  const triangulationType &triangulation) const;

    /**
     * @brief Ad-hoc struct for sorting simplices
     *
     * Adapted version of ttk::PersistentSimplexPairs::Simplex
     */
    template <size_t n>
    struct Simplex {
      /** Index in the triangulation */
      SimplexId id_{};
      /** Order field value of the simplex vertices, sorted in
          decreasing order */
      std::array<SimplexId, n> vertsOrder_{};
      /** To compare two vertices according to the filtration (lexicographic
       * order) */
      friend bool operator<(const Simplex<n> &lhs, const Simplex<n> &rhs) {
        return lhs.vertsOrder_ < rhs.vertsOrder_;
      }
    };

    /**
     * @brief \ref Simplex adaptation for edges
     */
    struct EdgeSimplex : Simplex<2> {
      template <typename triangulationType>
      void fillEdge(const SimplexId id,
                    const SimplexId *const offsets,
                    const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getEdgeVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getEdgeVertex(id, 1, this->vertsOrder_[1]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

    /**
     * @brief \ref Simplex adaptation for triangles
     */
    struct TriangleSimplex : Simplex<3> {
      template <typename triangulationType>
      void fillTriangle(const SimplexId id,
                        const SimplexId *const offsets,
                        const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getTriangleVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getTriangleVertex(id, 1, this->vertsOrder_[1]);
        triangulation.getTriangleVertex(id, 2, this->vertsOrder_[2]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offsets[this->vertsOrder_[2]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

    /**
     * @brief \ref Simplex adaptation for tetrahedra
     */
    struct TetraSimplex : Simplex<4> {
      template <typename triangulationType>
      void fillTetra(const SimplexId id,
                     const SimplexId *const offsets,
                     const triangulationType &triangulation) {
        this->id_ = id;
        triangulation.getCellVertex(id, 0, this->vertsOrder_[0]);
        triangulation.getCellVertex(id, 1, this->vertsOrder_[1]);
        triangulation.getCellVertex(id, 2, this->vertsOrder_[2]);
        triangulation.getCellVertex(id, 3, this->vertsOrder_[3]);
        this->vertsOrder_[0] = offsets[this->vertsOrder_[0]];
        this->vertsOrder_[1] = offsets[this->vertsOrder_[1]];
        this->vertsOrder_[2] = offsets[this->vertsOrder_[2]];
        this->vertsOrder_[3] = offsets[this->vertsOrder_[3]];
        // sort vertices in decreasing order
        std::sort(this->vertsOrder_.rbegin(), this->vertsOrder_.rend());
      }
    };

    template <typename triangulationType>
    void alloc(const triangulationType &triangulation) {
      Timer tm{};
      const auto dim{this->dg_.getDimensionality()};
      this->firstRepMin_.resize(triangulation.getNumberOfVertices());
      if(dim > 1) {
        this->firstRepMax_.resize(triangulation.getNumberOfCells());
      }
      if(dim > 2) {
        this->critEdges_.resize(triangulation.getNumberOfEdges());
        this->edgeTrianglePartner_.resize(triangulation.getNumberOfEdges(), -1);
        this->onBoundary_.resize(triangulation.getNumberOfEdges(), false);
        this->s2Mapping_.resize(triangulation.getNumberOfTriangles(), -1);
      }
      for(int i = 0; i < dim + 1; ++i) {
        this->pairedCritCells_[i].resize(
          this->dg_.getNumberOfCells(i, triangulation), false);
      }
      for(int i = 1; i < dim + 1; ++i) {
        this->critCellsOrder_[i].resize(
          this->dg_.getNumberOfCells(i, triangulation), -1);
      }
      this->printMsg("Memory allocations", 1.0, tm.getElapsedTime(), 1,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    void clear() {
      Timer tm{};
      this->firstRepMin_ = {};
      this->firstRepMax_ = {};
      this->edgeTrianglePartner_ = {};
      this->s2Mapping_ = {};
      this->critEdges_ = {};
      this->pairedCritCells_ = {};
      this->onBoundary_ = {};
      this->critCellsOrder_ = {};
      this->printMsg("Memory cleanup", 1.0, tm.getElapsedTime(), 1,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    dcg::DiscreteGradient dg_{};

    // factor memory allocations outside computation loops
    mutable std::vector<SimplexId> firstRepMin_{}, firstRepMax_{},
      edgeTrianglePartner_{}, s2Mapping_{};
    mutable std::vector<EdgeSimplex> critEdges_{};
    mutable std::array<std::vector<bool>, 4> pairedCritCells_{};
    mutable std::vector<bool> onBoundary_{};
    mutable std::array<std::vector<SimplexId>, 4> critCellsOrder_{};
  };
} // namespace ttk

template <typename triangulationType>
std::vector<std::vector<SimplexId>>
  ttk::DiscreteMorseSandwich::getSaddle1ToMinima(
    const std::vector<SimplexId> &criticalEdges,
    const triangulationType &triangulation) const {

  Timer tm{};

  std::vector<std::vector<SimplexId>> res(criticalEdges.size());

  // follow vpaths from 1-saddles to minima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < criticalEdges.size(); ++i) {
    auto &mins = res[i];

    const auto followVPath = [this, &mins, &triangulation](const SimplexId v) {
      std::vector<Cell> vpath{};
      this->dg_.getDescendingPath(Cell{0, v}, vpath, triangulation);
      Cell &lastCell = vpath.back();
      if(lastCell.dim_ == 0 && this->dg_.isCellCritical(lastCell)) {
        mins.emplace_back(lastCell.id_);
      }
    };

    // critical edge vertices
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(criticalEdges[i], 0, v0);
    triangulation.getEdgeVertex(criticalEdges[i], 1, v1);

    // follow vpath from each vertex of the critical edge
    followVPath(v0);
    followVPath(v1);
  }

  this->printMsg("Computed the descending 1-separatrices", 1.0,
                 tm.getElapsedTime(), this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return res;
}

template <typename triangulationType, typename GFS, typename GFSN, typename OB>
std::vector<std::vector<SimplexId>>
  ttk::DiscreteMorseSandwich::getSaddle2ToMaxima(
    const std::vector<SimplexId> &criticalCells,
    const GFS &getFaceStar,
    const GFSN &getFaceStarNumber,
    const OB &isOnBoundary,
    const triangulationType &triangulation,
    const bool ignoreBoundary) const {

  Timer tm{};

  const auto dim = this->dg_.getDimensionality();
  std::vector<std::vector<SimplexId>> res(criticalCells.size());

  // follow vpaths from 2-saddles to maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < criticalCells.size(); ++i) {
    const auto sid = criticalCells[i];
    auto &maxs = res[i];

    const auto followVPath
      = [this, dim, &maxs, &triangulation](const SimplexId v) {
          std::vector<Cell> vpath{};
          this->dg_.getAscendingPath(Cell{dim, v}, vpath, triangulation);
          Cell &lastCell = vpath.back();
          if(lastCell.dim_ == dim && this->dg_.isCellCritical(lastCell)) {
            maxs.emplace_back(lastCell.id_);
          } else if(lastCell.dim_ == dim - 1) {
            maxs.emplace_back(-1);
          }
        };

    const auto starNumber = getFaceStarNumber(sid);

    for(SimplexId j = 0; j < starNumber; ++j) {
      SimplexId cellId{};
      getFaceStar(sid, j, cellId);
      followVPath(cellId);
    }

    // ignoreBoundary: skip triplets with saddles on boundary
    if(!ignoreBoundary && isOnBoundary(sid)) {
      // critical saddle is on boundary
      maxs.emplace_back(-1);
    }
  }

  this->printMsg("Computed the ascending 1-separatrices", 1.0,
                 tm.getElapsedTime(), this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return res;
}

template <typename triangulationType>
void ttk::DiscreteMorseSandwich::getMinSaddlePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<bool> &pairedMinima,
  std::vector<bool> &paired1Saddles,
  const std::vector<SimplexId> &criticalEdges,
  const std::vector<SimplexId> &critEdgesOrder,
  const SimplexId *const offsets,
  const triangulationType &triangulation) const {

  Timer tm{};

  auto saddle1ToMinima = getSaddle1ToMinima(criticalEdges, triangulation);

  Timer tmseq{};

  auto &firstRep{this->firstRepMin_};
  std::iota(firstRep.begin(), firstRep.end(), 0);
  std::vector<tripletType> sadMinTriplets{};

  for(size_t i = 0; i < saddle1ToMinima.size(); ++i) {
    auto &mins = saddle1ToMinima[i];
    const auto s1 = criticalEdges[i];
    // remove duplicates
    TTK_PSORT(this->threadNumber_, mins.begin(), mins.end());
    const auto last = std::unique(mins.begin(), mins.end());
    mins.erase(last, mins.end());
    if(mins.size() != 2) {
      continue;
    }
    sadMinTriplets.emplace_back(tripletType{s1, mins[0], mins[1]});
  }

  tripletsToPersistencePairs(pairs, pairedMinima, paired1Saddles, firstRep,
                             sadMinTriplets, critEdgesOrder.data(), offsets, 0);

  const auto nMinSadPairs = pairs.size();

  this->printMsg(
    "Computed " + std::to_string(nMinSadPairs) + " min-saddle pairs", 1.0,
    tm.getElapsedTime(), this->threadNumber_);

  this->printMsg("min-saddle pairs sequential part", 1.0,
                 tmseq.getElapsedTime(), 1, debug::LineMode::NEW,
                 debug::Priority::VERBOSE);
}

template <typename triangulationType>
void ttk::DiscreteMorseSandwich::getMaxSaddlePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<bool> &pairedMaxima,
  std::vector<bool> &pairedSaddles,
  const std::vector<SimplexId> &criticalSaddles,
  const std::vector<SimplexId> &critSaddlesOrder,
  const std::vector<SimplexId> &critMaxsOrder,
  const triangulationType &triangulation,
  const bool ignoreBoundary) const {

  Timer tm{};

  const auto dim = this->dg_.getDimensionality();

  auto saddle2ToMaxima
    = dim == 3
        ? getSaddle2ToMaxima(
          criticalSaddles,
          [&triangulation](const SimplexId a, const SimplexId i, SimplexId &r) {
            return triangulation.getTriangleStar(a, i, r);
          },
          [&triangulation](const SimplexId a) {
            return triangulation.getTriangleStarNumber(a);
          },
          [&triangulation](const SimplexId a) {
            return triangulation.isTriangleOnBoundary(a);
          },
          triangulation, ignoreBoundary)
        : getSaddle2ToMaxima(
          criticalSaddles,
          [&triangulation](const SimplexId a, const SimplexId i, SimplexId &r) {
            return triangulation.getEdgeStar(a, i, r);
          },
          [&triangulation](const SimplexId a) {
            return triangulation.getEdgeStarNumber(a);
          },
          [&triangulation](const SimplexId a) {
            return triangulation.isEdgeOnBoundary(a);
          },
          triangulation, ignoreBoundary);

  Timer tmseq{};

  auto &firstRep{this->firstRepMax_};
  std::iota(firstRep.begin(), firstRep.end(), 0);
  std::vector<tripletType> sadMaxTriplets{};

  for(size_t i = 0; i < saddle2ToMaxima.size(); ++i) {
    auto &maxs = saddle2ToMaxima[i];
    // remove duplicates
    TTK_PSORT(this->threadNumber_, maxs.begin(), maxs.end(),
              [](const SimplexId a, const SimplexId b) {
                // positive values (actual maxima) before negative ones
                // (boundary component id)
                if(a * b >= 0) {
                  return a < b;
                } else {
                  return a > b;
                }
              });
    const auto last = std::unique(maxs.begin(), maxs.end());
    maxs.erase(last, maxs.end());

    // remove "doughnut" configurations: two ascending separatrices
    // leading to the same maximum/boundary component
    if(maxs.size() != 2) {
      continue;
    }

    const auto s2 = criticalSaddles[i];
    if(!pairedSaddles[s2]) {
      sadMaxTriplets.emplace_back(tripletType{s2, maxs[0], maxs[1]});
    }
  }

  const auto nMinSadPairs = pairs.size();

  tripletsToPersistencePairs(pairs, pairedMaxima, pairedSaddles, firstRep,
                             sadMaxTriplets, critSaddlesOrder.data(),
                             critMaxsOrder.data(), dim - 1);

  const auto nSadMaxPairs = pairs.size() - nMinSadPairs;

  this->printMsg(
    "Computed " + std::to_string(nSadMaxPairs) + " saddle-max pairs", 1.0,
    tm.getElapsedTime(), this->threadNumber_);

  this->printMsg("saddle-max pairs sequential part", 1.0,
                 tmseq.getElapsedTime(), 1, debug::LineMode::NEW,
                 debug::Priority::VERBOSE);
}

template <typename triangulationType, typename Container>
SimplexId ttk::DiscreteMorseSandwich::eliminateBoundariesSandwich(
  const SimplexId s2,
  std::vector<bool> &onBoundary,
  std::vector<Container> &s2Boundaries,
  const std::vector<SimplexId> &s2Mapping,
  const std::vector<SimplexId> &partners,
  const triangulationType &triangulation) const {

  auto &boundaryIds{s2Boundaries[s2Mapping[s2]]};
  const auto addBoundary = [&boundaryIds, &onBoundary](const SimplexId e) {
    // add edge e to boundaryIds/onBoundary modulo 2
    if(!onBoundary[e]) {
      boundaryIds.emplace(e);
      onBoundary[e] = true;
    } else {
      const auto it = boundaryIds.find(e);
      boundaryIds.erase(it);
      onBoundary[e] = false;
    }
  };

  if(!boundaryIds.empty()) {
    // restore previously computed s2 boundary
    for(const auto e : boundaryIds) {
      onBoundary[e] = true;
    }
  } else {
    // init cascade with s2 triangle boundary (3 edges)
    for(SimplexId i = 0; i < 3; ++i) {
      SimplexId e{};
      triangulation.getTriangleEdge(s2, i, e);
      addBoundary(e);
    }
  }

  while(!boundaryIds.empty()) {
    // tau: youngest edge on boundary
    const auto tau{*boundaryIds.begin()};
    // use the Discrete Gradient to find a triangle paired to tau
    auto pTau{this->dg_.getPairedCell(Cell{1, tau}, triangulation)};
    bool critical{false};
    if(pTau == -1) {
      // maybe tau is critical and paired to a critical triangle
      pTau = partners[tau];
      critical = true;
    }
    if(pTau == -1) {
      // tau is critical and not paired
      return tau;
    }
    // expand boundary
    if(critical && s2Mapping[pTau] != -1) { // pTau is a non-paired 2-saddle
      // merge pTau boundary into s2 boundary
      for(const auto e : s2Boundaries[s2Mapping[pTau]]) {
        addBoundary(e);
      }
    } else { // pTau is a regular triangle
      // add pTau triangle boundary (3 edges)
      for(SimplexId i = 0; i < 3; ++i) {
        SimplexId e{};
        triangulation.getTriangleEdge(pTau, i, e);
        addBoundary(e);
      }
    }
  }

  return -1;
}

template <typename triangulationType>
void ttk::DiscreteMorseSandwich::getSaddleSaddlePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<bool> &paired1Saddles,
  std::vector<bool> &paired2Saddles,
  const bool exportBoundaries,
  std::vector<GeneratorType> &boundaries,
  const std::vector<SimplexId> &critical1Saddles,
  const std::vector<SimplexId> &critical2Saddles,
  const std::vector<SimplexId> &crit1SaddlesOrder,
  const triangulationType &triangulation) const {

  Timer tm2{};
  const auto nSadExtrPairs = pairs.size();

  // 1- and 2-saddles yet to be paired
  std::vector<SimplexId> saddles1{}, saddles2{};
  // filter out already paired 1-saddles (edge id)
  for(const auto s1 : critical1Saddles) {
    if(!paired1Saddles[s1]) {
      saddles1.emplace_back(s1);
    }
  }
  // filter out already paired 2-saddles (triangle id)
  for(const auto s2 : critical2Saddles) {
    if(!paired2Saddles[s2]) {
      saddles2.emplace_back(s2);
    }
  }

  Timer tmpar{};

  // sort every triangulation edges by filtration order
  const auto &edgesFiltrOrder{crit1SaddlesOrder};

  auto &onBoundary{this->onBoundary_};
  auto &edgeTrianglePartner{this->edgeTrianglePartner_};

  const auto cmpEdges
    = [&edgesFiltrOrder](const SimplexId a, const SimplexId b) {
        return edgesFiltrOrder[a] > edgesFiltrOrder[b];
      };
  using Container = std::set<SimplexId, decltype(cmpEdges)>;
  std::vector<Container> s2Boundaries(saddles2.size(), Container(cmpEdges));
  // unpaired critical triangle id -> index in saddle2 vector
  auto &s2Mapping{this->s2Mapping_};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    s2Mapping[saddles2[i]] = i;
  }

  // first parallel pass to pre-determine the boundary of every
  // unpaired 2-saddle to the first unpaired 1-saddle on its wall
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
  firstprivate(onBoundary)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    // 2-saddles sorted in increasing order
    const auto s2 = saddles2[i];
    this->eliminateBoundariesSandwich(s2, onBoundary, s2Boundaries, s2Mapping,
                                      edgeTrianglePartner, triangulation);
    // reset onBoundary to false
    for(const auto e : s2Boundaries[i]) {
      onBoundary[e] = false;
    }
  }

  this->printMsg("Precomputed 2-saddle boundaries", 1.0, tmpar.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  Timer tmseq{};

  for(size_t i = 0; i < saddles2.size(); ++i) {
    // 2-saddles sorted in increasing order
    const auto s2 = saddles2[i];
    SimplexId s1{-1};

    if(!s2Boundaries[i].empty()
       && edgeTrianglePartner[*s2Boundaries[i].begin()] == -1) {
      // use shortcut if first 1-saddle on wall is non-paired
      s1 = *s2Boundaries[i].begin();
    } else {
      s1 = this->eliminateBoundariesSandwich(s2, onBoundary, s2Boundaries,
                                             s2Mapping, edgeTrianglePartner,
                                             triangulation);
    }

    if(s1 != -1) {
      // we found a pair
      pairs.emplace_back(s1, s2, 1);
      paired1Saddles[s1] = true;
      paired2Saddles[s2] = true;
      edgeTrianglePartner[s1] = s2;
      // reset onBoundary to false
      for(const auto e : s2Boundaries[i]) {
        onBoundary[e] = false;
      }
    }
  }

  if(exportBoundaries) {
    boundaries.resize(s2Boundaries.size());
    for(size_t i = 0; i < boundaries.size(); ++i) {
      const auto &boundSet{s2Boundaries[i]};
      if(boundSet.empty()) {
        continue;
      }
      boundaries[i] = {
        {boundSet.begin(), boundSet.end()},
        saddles2[i],
        std::array<SimplexId, 2>{
          this->dg_.getCellGreaterVertex(Cell{2, saddles2[i]}, triangulation),
          this->dg_.getCellGreaterVertex(
            Cell{1, *boundSet.begin()}, triangulation),
        }};
    }
  }

  const auto nSadSadPairs = pairs.size() - nSadExtrPairs;

  this->printMsg(
    "Computed " + std::to_string(nSadSadPairs) + " saddle-saddle pairs", 1.0,
    tm2.getElapsedTime(), this->threadNumber_);

  this->printMsg("saddle-saddle pairs sequential part", 1.0,
                 tmseq.getElapsedTime(), 1, debug::LineMode::NEW,
                 debug::Priority::VERBOSE);
}

template <typename triangulationType>
void ttk::DiscreteMorseSandwich::extractCriticalCells(
  std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
  std::array<std::vector<SimplexId>, 4> &critCellsOrder,
  const SimplexId *const offsets,
  const triangulationType &triangulation,
  const bool sortEdges) const {

  Timer tm{};
  const auto dim = this->dg_.getDimensionality();

  for(int i = 0; i < dim + 1; ++i) {

    // map: store critical cell per dimension per thread
    std::vector<std::vector<SimplexId>> critCellsPerThread(this->threadNumber_);

    const SimplexId numberOfCells
      = this->dg_.getNumberOfCells(i, triangulation);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId j = 0; j < numberOfCells; ++j) {
#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
      if(this->dg_.isCellCritical(i, j)) {
        critCellsPerThread[tid].emplace_back(j);
      }
    }

    // reduce: aggregate critical cells per thread
    criticalCellsByDim[i] = std::move(critCellsPerThread[0]);
    for(size_t j = 1; j < critCellsPerThread.size(); ++j) {
      const auto &vec{critCellsPerThread[j]};
      criticalCellsByDim[i].insert(
        criticalCellsByDim[i].end(), vec.begin(), vec.end());
    }
  }

  this->printMsg("Extracted critical cells", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::VERBOSE);

  // memory allocations
  auto &critEdges{this->critEdges_};
  if(!sortEdges) {
    critEdges.resize(criticalCellsByDim[1].size());
  }
  std::vector<TriangleSimplex> critTriangles(criticalCellsByDim[2].size());
  std::vector<TetraSimplex> critTetras(criticalCellsByDim[3].size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
    if(sortEdges) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < critEdges.size(); ++i) {
        critEdges[i].fillEdge(i, offsets, triangulation);
      }
    } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < critEdges.size(); ++i) {
        critEdges[i].fillEdge(criticalCellsByDim[1][i], offsets, triangulation);
      }
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critTriangles.size(); ++i) {
      critTriangles[i].fillTriangle(
        criticalCellsByDim[2][i], offsets, triangulation);
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critTetras.size(); ++i) {
      critTetras[i].fillTetra(criticalCellsByDim[3][i], offsets, triangulation);
    }
  }

  TTK_PSORT(this->threadNumber_, critEdges.begin(), critEdges.end());
  TTK_PSORT(this->threadNumber_, critTriangles.begin(), critTriangles.end());
  TTK_PSORT(this->threadNumber_, critTetras.begin(), critTetras.end());

  if(sortEdges) {
    TTK_PSORT(this->threadNumber_, criticalCellsByDim[1].begin(),
              criticalCellsByDim[1].end(),
              [&critCellsOrder](const SimplexId a, const SimplexId b) {
                return critCellsOrder[1][a] < critCellsOrder[1][b];
              });
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critEdges.size(); ++i) {
      criticalCellsByDim[1][i] = critEdges[i].id_;
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critEdges.size(); ++i) {
      critCellsOrder[1][critEdges[i].id_] = i;
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critTriangles.size(); ++i) {
      criticalCellsByDim[2][i] = critTriangles[i].id_;
      critCellsOrder[2][critTriangles[i].id_] = i;
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < critTetras.size(); ++i) {
      criticalCellsByDim[3][i] = critTetras[i].id_;
      critCellsOrder[3][critTetras[i].id_] = i;
    }
  }

  this->printMsg("Extracted & sorted critical cells", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}

template <typename triangulationType>
int ttk::DiscreteMorseSandwich::computePersistencePairs(
  std::vector<PersistencePair> &pairs,
  const SimplexId *const offsets,
  const triangulationType &triangulation,
  const bool ignoreBoundary) {

  // allocate memory
  this->alloc(triangulation);

  Timer tm{};
  pairs.clear();
  const auto dim = this->dg_.getDimensionality();

  // get every critical cell sorted them by dimension
  std::array<std::vector<SimplexId>, 4> criticalCellsByDim{};
  // holds the critical cells order
  auto &critCellsOrder{this->critCellsOrder_};

  this->extractCriticalCells(
    criticalCellsByDim, critCellsOrder, offsets, triangulation, dim == 3);

  // if minima are paired
  auto &pairedMinima{this->pairedCritCells_[0]};
  // if 1-saddles are paired
  auto &paired1Saddles{this->pairedCritCells_[1]};
  // if 2-saddles are paired
  auto &paired2Saddles{this->pairedCritCells_[dim - 1]};
  // if maxima are paired
  auto &pairedMaxima{this->pairedCritCells_[dim]};

  // minima - saddle pairs
  this->getMinSaddlePairs(pairs, pairedMinima, paired1Saddles,
                          criticalCellsByDim[1], critCellsOrder[1], offsets,
                          triangulation);

  // connected components (global min/max pair)
  size_t nConnComp{};
  for(const auto min : criticalCellsByDim[0]) {
    if(!pairedMinima[min]) {
      pairs.emplace_back(min, -1, 0);
      pairedMinima[min] = true;
      nConnComp++;
    }
  }

  if(dim == 1) {
    // early return in 1D
    this->printMsg(
      "Computed " + std::to_string(pairs.size()) + " persistence pairs", 1.0,
      tm.getElapsedTime(), this->threadNumber_);
    return 0;
  }

  // saddle - maxima pairs
  this->getMaxSaddlePairs(pairs, pairedMaxima, paired2Saddles,
                          criticalCellsByDim[dim - 1], critCellsOrder[dim - 1],
                          critCellsOrder[dim], triangulation, ignoreBoundary);

  // saddle - saddle pairs
  if(dim == 3 && !criticalCellsByDim[1].empty()
     && !criticalCellsByDim[2].empty()) {
    std::vector<GeneratorType> tmp{};
    this->getSaddleSaddlePairs(
      pairs, paired1Saddles, paired2Saddles, false, tmp, criticalCellsByDim[1],
      criticalCellsByDim[2], critCellsOrder[1], triangulation);
  }

  if(std::is_same<triangulationType, ttk::ExplicitTriangulation>::value) {
    // create infinite pairs from non-paired 1-saddles
    size_t nHandles{}, nCavities{}, nNonPairedMax{};
    if((dim == 2 && !ignoreBoundary) || dim == 3) {
      for(const auto s1 : criticalCellsByDim[1]) {
        if(!paired1Saddles[s1]) {
          paired1Saddles[s1] = true;
          // topological handles
          pairs.emplace_back(s1, -1, 1);
          nHandles++;
        }
      }
    }
    if(dim == 3 && !ignoreBoundary) {
      for(const auto s2 : criticalCellsByDim[2]) {
        if(!paired2Saddles[s2]) {
          paired2Saddles[s2] = true;
          // cavities
          pairs.emplace_back(s2, -1, 2);
          nCavities++;
        }
      }
    }
    if(dim == 2 && !ignoreBoundary) {
      for(const auto max : criticalCellsByDim[dim]) {
        if(!pairedMaxima[max]) {
          pairs.emplace_back(max, -1, 2);
          nNonPairedMax++;
        }
      }
    }

    // print Betti numbers
    std::vector<std::vector<std::string>> rows{
      {" #Connected components", std::to_string(nConnComp)},
      {" #Topological handles", std::to_string(nHandles)},
      {" #Cavities", std::to_string(nCavities)},
      {" #Boundary components", std::to_string((dim == 3 ? nCavities : nHandles)
                                               + nConnComp - nNonPairedMax)},
    };

    this->printMsg(rows, debug::Priority::DETAIL);
  }

  this->printMsg(
    "Computed " + std::to_string(pairs.size()) + " persistence pairs", 1.0,
    tm.getElapsedTime(), this->threadNumber_);

  this->displayStats(pairs, criticalCellsByDim, pairedMinima, paired1Saddles,
                     paired2Saddles, pairedMaxima);

  // free memory
  this->clear();

  return 0;
}
