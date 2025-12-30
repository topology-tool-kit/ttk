/// \ingroup baseCode
/// \class ttk::VectorSimplification
/// \author Tanner Finken <finkent@arizona.edu>
/// \author Josh A. Levine <josh@cs.arizona.edu>
/// \date November 2024.
///
/// \brief TTK %VectorSimplification processing package.
///
/// %VectorSimplification computes a Weight Curve by using the
/// %Discrete Morse-Theory %DiscreteVectorField simplification
///  with flipping paths based on the lowest weight value.
///  The implementation is slightly different than described in
///  the paper as this will 'shortcut' the saddle-orbit cancellation pairs
///  by tracing through the orbit by flipping dimensionality (e.g. 0-1 to 1-2
///  V-Path) and either including the orbit or not depending on the given flag.
///
/// \b Related \b publication \n
/// "Localized Evaluation for Constructing Discrete Vector Fields" \n
/// Tanner Finken, Julien Tierny, Joshua A. Levine \n
/// IEEE Vis 2024.
///
/// Additionally, a large portion of the tracing code has been adapted from
/// ttk::DiscreteMorseSandwich.
///
/// \sa ttk::DiscreteMorseSandwich
/// \sa ttk::dcvf::DiscreteVectorField
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/">Discrete
///   Vector Field Topology example</a> \n
///

#pragma once

#include <DiscreteVectorField.h>

using ttk::dcg::Cell;

#include <algorithm>
#include <numeric>

namespace ttk {
  class VectorSimplification : virtual public Debug {
  public:
    VectorSimplification();

    /***
     * @brief The discrete vector field this class will be building and
     * simplifing Can be used to access information about the simplified field.
     */
    dcvf::DiscreteVectorField dcvf_{};

    /**
     * @brief Candidate pair struct for information of pairs connected by
     * V-Paths We need to know the start and end points, type of connection,
     * weight of connection and which cell to follow for ascending or descending
     * path.
     */
    struct CandidatePair {
      /** first (lower dim/birth) simplex cell id */
      SimplexId birth;
      /** Next cell to follow for asc/des path call */
      Cell nextCell;
      /** second (higher dim/death) simplex cell id */
      SimplexId death;
      /** Number of alternations performed while tracing V-Path  */
      int alternations{0};
      /** Will simplifying this pair generate an orbit
       * (i.e., does the saddle connect to the same spot)
       */
      bool generateOrbit{false};

      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      int type;
      /** Weight of simplification pair in discrete vector field
       *  (computed using the alternating sum along pairs and anti-pairs)
       */
      float weight;

      CandidatePair(SimplexId b, Cell next, SimplexId d, int t, float w)
        : birth{b}, nextCell{next}, death{d}, type{t}, weight{w} {
      }
    };

    struct PlotPoint {
      /* Number of Critical points remaining*/
      SimplexId numCP;
      /* Number of Sink(dim=0) Critical points remaining*/
      SimplexId numSinks;
      /* Number of Source(dim=dim)Critical points remaining*/
      SimplexId numSources;
      /* Number of Orbits generated so far */
      SimplexId orbitsAdded{0};
      /* Number of Orbits removed so far */
      SimplexId orbitsRemoved{0};
      /* Weight Threshold */
      double weight;

      PlotPoint(SimplexId num, SimplexId sinks, SimplexId sources, double w)
        : numCP{num}, numSinks{sinks}, numSources{sources}, weight{w} {
      }
    };

    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      this->dcvf_.preconditionTriangulation(data);
    }

    template <typename dataType, typename triangulationType>
    inline int buildField(const void *const vectors,
                          const size_t vectorsMTime,
                          const triangulationType &triangulation) {
      this->dcvf_.setDebugLevel(this->debugLevel_);
      this->dcvf_.setThreadNumber(this->threadNumber_);
      this->dcvf_.setInputVectorField(vectors, vectorsMTime);
      return this->dcvf_.buildField<dataType, triangulationType>(triangulation);
    }

    template <typename dataType, typename triangulationType>
    inline int performSimplification(const int criticalThreshold,
                                     const bool storePlotPoints,
                                     std::vector<PlotPoint> &listOfPlotPoints,
                                     const triangulationType &triangulation) {
      Timer tm{}; // Time the entire simplification
      // Extract Critical Simplices
      std::array<std::vector<SimplexId>, 4> criticalCellsByDim;
      this->dcvf_.getCriticalPoints(criticalCellsByDim, triangulation);
      int dim = this->dcvf_.getDimensionality();
      int numCriticalPoints = criticalCellsByDim[0].size()
                              + criticalCellsByDim[1].size()
                              + criticalCellsByDim[2].size();
      int numSinks = criticalCellsByDim[0].size();
      int numSources = criticalCellsByDim[2].size();
      int orbitsAdded{0}, orbitsRemoved{0};
      double weightThres{-1000000};
      if(numCriticalPoints < criticalThreshold)
        return 0;
      // Trace the saddles and track the candidate pairs
      std::vector<CandidatePair> pairs;
      if(dim == 2) {
        // saddle downward pairs
        this->getDescSaddlePairs<dataType, triangulationType>(
          pairs, criticalCellsByDim[1], triangulation);
        // saddle upward pairs
        this->getAscSaddlePairs<dataType, triangulationType>(
          pairs, criticalCellsByDim[dim - 1], triangulation);
      } else if(dim == 3) {
        this->printErr("This filter can not simplify for 3D yet.");
      }
      // Sort the candidate pairs to remove the smallest ones first and update
      // the paths as needed
      if(dim == 2) {
        const auto orderPairs
          = [&](const CandidatePair &a, const CandidatePair &b) -> bool {
          return a.weight > b.weight;
        };

        // Type alias for the priority queue
        using pqType
          = std::priority_queue<CandidatePair, std::vector<CandidatePair>,
                                decltype(orderPairs)>;
        pqType options{orderPairs};
        for(auto pair : pairs) {
          options.push(pair);
        }
        pairs.clear();
        // Simplify based on smallest options
        while(!options.empty()
              and numCriticalPoints - 2
                    >= static_cast<int>(criticalThreshold)) {
          CandidatePair bestOption = options.top();
          // Ensure the critical pair is still valid  (need the end(is it a
          // cycle that still exists?))
          if(bestOption.type == 0) {
            if(!this->dcvf_.isCellCritical(1, bestOption.death)) {
              options.pop();
              continue;
            }
            if(!this->dcvf_.isCellCritical(0, bestOption.birth)) {
              options.pop();
              // Trace the saddle (based on ascending(upward) or
              // descending(downward))
              std::vector<SimplexId> saddle{bestOption.death};
              if(bestOption.nextCell.dim_ == 0) {
                this->getDescSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              } else {
                this->getAscSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              }
              for(auto pair : pairs) {
                options.push(pair);
              }
              pairs.clear();
              continue;
            }

          } else if(bestOption.type == 2) {
            if(!this->dcvf_.isCellCritical(1, bestOption.birth)) {
              options.pop();
              continue;
            }
            if(!this->dcvf_.isCellCritical(2, bestOption.death)) {
              options.pop();
              // Trace the saddle (based on ascending(upward) or
              // descending(downward))
              std::vector<SimplexId> saddle{bestOption.birth};
              if(bestOption.nextCell.dim_ == 2) {
                this->getAscSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              } else {
                this->getDescSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              }
              for(auto pair : pairs) {
                options.push(pair);
              }
              pairs.clear();
              continue;
            }
          }
          // Trace the path to confirm it still ends at the same CP
          std::vector<Cell> vpath;
          // Add the saddle
          if(bestOption.type == 0) {
            vpath.emplace_back(Cell(1, bestOption.death));
          } else if(bestOption.type == 2) {
            vpath.emplace_back(Cell(1, bestOption.birth));
          }
          // Then trace desc/asc path
          if(bestOption.nextCell.dim_ == 0) {
            this->dcvf_.getDescendingPath<dataType, triangulationType>(
              bestOption.nextCell, vpath, triangulation, false);
          } else {
            this->dcvf_.getAscendingPath<dataType, triangulationType>(
              bestOption.nextCell, vpath, triangulation, false);
          }
          if(bestOption.type == 0) {
            if(bestOption.birth
               != vpath.back().id_) { // Ends elsewhere (needs retraced)
              options.pop();
              // Trace the new saddle destinations
              std::vector<SimplexId> saddle{bestOption.death};
              if(bestOption.nextCell.dim_ == 0) {
                this->getDescSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              } else {
                this->getAscSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              }
              for(auto pair : pairs) {
                options.push(pair);
              }
              pairs.clear();
              continue;
            }
          } else if(bestOption.type == 2) {
            if(bestOption.death
               != vpath.back().id_) { // Ends elsewhere (needs retraced)
              options.pop();
              // Trace the new saddle destinations
              std::vector<SimplexId> saddle{bestOption.birth};
              if(bestOption.nextCell.dim_ == 2) {
                this->getAscSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              } else {
                this->getDescSaddlePairs<dataType, triangulationType>(
                  pairs, saddle, triangulation);
              }
              for(auto pair : pairs) {
                options.push(pair);
              }
              pairs.clear();
              continue;
            }
          }

          // Path has been Validated, so we can Trace and simplify

          // Update weight threshold to be max of newest weight and previous
          // weight threshold Since lower weights are possible
          weightThres
            = std::max(weightThres, static_cast<double>(bestOption.weight));

          // First add plot point if desired (before)
          if(storePlotPoints) {
            listOfPlotPoints.emplace_back(
              numCriticalPoints, numSinks, numSources, weightThres);
            listOfPlotPoints.back().orbitsAdded = orbitsAdded;
            listOfPlotPoints.back().orbitsRemoved = orbitsRemoved;
          }

          // Then update tracking values
          if(bestOption.type == 0) {
            numSinks--;
          } else if(bestOption.type == 2) {
            numSources--;
          }
          if(bestOption.generateOrbit) {
            orbitsAdded++;
          }
          if(bestOption.alternations > 0) {
            orbitsRemoved = orbitsRemoved + bestOption.alternations;
          }
          // Then add plot point if desired (after)
          if(storePlotPoints) {
            listOfPlotPoints.emplace_back(
              numCriticalPoints - 2, numSinks, numSources, weightThres);
            listOfPlotPoints.back().orbitsAdded = orbitsAdded;
            listOfPlotPoints.back().orbitsRemoved = orbitsRemoved;
          }

          if(bestOption.nextCell.dim_ == 0) {
            if(isAlternatingVpath(
                 vpath)) { // Check if the Vpath alternating type 1-2 to 0-1
              this->dcvf_.reverseAlternatingPath(vpath, triangulation);
            } else {
              this->dcvf_.reverseDescendingPath(vpath, triangulation);
            }
            numCriticalPoints -= 2;
          } else {
            if(isAlternatingVpath(
                 vpath)) { // Check if the Vpath alternating type 0-1 to 1-2
              this->dcvf_.reverseAlternatingPath(vpath, triangulation);
            } else {
              this->dcvf_.reverseAscendingPath(vpath, triangulation);
            }
            numCriticalPoints -= 2;
          }
          options.pop();
        }
      }
      this->printMsg("Simplified to " + std::to_string(numCriticalPoints)
                       + " critical point(s)",
                     1.0, tm.getElapsedTime(), this->threadNumber_);

      return 0;
    }

    /**
     * @brief Ugly hack to avoid a call to buildField()
     *
     * An externally computed field can be retrofitted into this
     * class using move semantics with setField().
     * The internal field can be fetched back with getField()
     * once the simplification is performed.
     *
     * @param[in] dcvf External discrete vector field instance
     */
    inline void setField(ttk::dcvf::DiscreteVectorField &&dcvf) {
      this->dcvf_ = std::move(dcvf);
    }

    void setFullOrbitSimplification(bool doFullOrbit) {
      this->dcvf_.setReverseFullOrbit(doFullOrbit);
    }

    inline ttk::dcvf::DiscreteVectorField &&getField() {
      // WARNING, this function will change the location of the discrete field
      return std::move(this->dcvf_);
    }

  protected:
    /**
     * @brief Follow the descending 1-separatrices to compute the saddles ->
     * extrema association
     *
     * @param[in] criticalEdges Critical edges identifiers
     * @param[in] triangulation Triangulation
     *
     * @return a vector of Candidate Pairs per 1-saddle
     */
    template <typename dataType, typename triangulationType>
    std::vector<std::vector<CandidatePair>>
      getSaddle1ToDescPair(const std::vector<SimplexId> &criticalEdges,
                           const triangulationType &triangulation) const;

    /**
     * @brief Follow the ascending 1-separatrices to compute the saddles ->
     * extrema association
     *
     * @param[in] criticalCells Critical cells identifiers
     * @param[in] getFaceStar Either getEdgeStar (in 2D) or getTriangleStar
     * (in 3D)
     * @param[in] getFaceStarNumber Either getEdgeStarNumber (in 2D) or
     * getTriangleStarNumber (in 3D)
     * @param[in] isOnBoundary Either isEdgeOnBoundary (in 2D) or
     * isTriangleOnBoundary (in 3D)
     * @param[in] triangulation Triangulation
     * @param[in] dummyVariable provided for compiler to establish dataType
     *
     * @return a vector of Candidate Pairs per 2-saddle
     */
    template <typename dataType,
              typename triangulationType,
              typename GFS,
              typename GFSN,
              typename OB>
    std::vector<std::vector<CandidatePair>>
      getSaddle2ToAscPair(const std::vector<SimplexId> &criticalCells,
                          const GFS &getFaceStar,
                          const GFSN &getFaceStarNumber,
                          const OB &isOnBoundary,
                          const triangulationType &triangulation,
                          const dataType dummyVariable) const;

    /**
     * @brief Compute the candidate pairs from descending paths from saddles
     *
     * @param[out] pairs Output candidate pairs
     * @param[in] criticalEdges List of 1-saddles (or maxima in 1D)
     * @param[in] triangulation Triangulation
     */
    template <typename dataType, typename triangulationType>
    void getDescSaddlePairs(std::vector<CandidatePair> &pairs,
                            const std::vector<SimplexId> &criticalEdges,
                            const triangulationType &triangulation) const;

    /**
     * @brief Compute the candidate pairs from ascending paths from saddles
     *
     * @param[out] pairs Output candidate pairs
     * @param[in] criticalSaddles List of 2-saddles (or 1-saddles in 2D)
     * @param[in] triangulation Triangulation
     */
    template <typename dataType, typename triangulationType>
    void getAscSaddlePairs(std::vector<CandidatePair> &pairs,
                           const std::vector<SimplexId> &criticalSaddles,
                           const triangulationType &triangulation) const;

    /**
     * @brief Determine if the VPath has 'alternating' behavior of dimension
     *  ex: 0-1 to 1-2.
     *
     * @param[in] vpath Collection of cells along the VPath
     *
     * @return true if rotates (i.e., 0-1 to 1-2 VPaths)
     */
    inline bool isAlternatingVpath(std::vector<Cell> &vpath);

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
      const std::vector<CandidatePair> &pairs,
      const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
      const std::vector<bool> &pairedMinima,
      const std::vector<bool> &paired1Saddles,
      const std::vector<bool> &paired2Saddles,
      const std::vector<bool> &pairedMaxima) const;
  };
} // namespace ttk

bool ttk::VectorSimplification::isAlternatingVpath(std::vector<Cell> &vpath) {
  int pathSize = vpath.size();
  // Assume it starts with saddle
  if(pathSize > 2) {
    int previousDim = vpath[1].dim_;
    int previousSecondDim = vpath[2].dim_;
    for(int i = 1; i < pathSize; i++) {
      if((i % 2) == 0) {
        if(vpath[i].dim_ != previousSecondDim) {
          return true;
        }
      } else {
        if(vpath[i].dim_ != previousDim) {
          return true;
        }
      }
    }
  }

  return false;
}

template <typename dataType, typename triangulationType>
std::vector<std::vector<ttk::VectorSimplification::CandidatePair>>
  ttk::VectorSimplification::getSaddle1ToDescPair(
    const std::vector<SimplexId> &criticalEdges,
    const triangulationType &triangulation) const {

  const auto dim = this->dcvf_.getDimensionality();
  std::vector<std::vector<CandidatePair>> res(criticalEdges.size());

  // follow vpaths from 1-saddles to descending paths (typically minima)
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < criticalEdges.size(); ++i) {
    auto &mins = res[i];
    const auto sid = criticalEdges[i];

    const auto followVPath = [this, dim, sid, &mins,
                              &triangulation](const SimplexId v) {
      std::vector<Cell> vpath{};
      vpath.emplace_back(Cell{1, sid});
      int alternates
        = this->dcvf_.getDescendingPath<dataType, triangulationType>(
          Cell{0, v}, vpath, triangulation, false);
      Cell &lastCell = vpath.back();
      if(lastCell.dim_ == 0 && this->dcvf_.isCellCritical(lastCell)) {
        auto weight = this->dcvf_.getPersistence<dataType, triangulationType>(
          vpath, triangulation);
        mins.emplace_back(lastCell.id_, Cell(0, v), sid, 0, weight);
        mins.back().alternations = alternates;
      } else if(lastCell.dim_ == dim && this->dcvf_.isCellCritical(lastCell)) {
        // Might change dimensions of end CP
        auto weight = this->dcvf_.getPersistence<dataType, triangulationType>(
          vpath, triangulation);
        mins.emplace_back(sid, Cell{0, v}, lastCell.id_, 2, weight);
        mins.back().alternations = alternates;
      }
    };
#ifndef TTK_ENABLE_KAMIKAZE
    // Test for valid sid ()
    if(sid < 0 || sid > triangulation.getNumberOfEdges()) {
      std::cout << "[WARNING] Not valid sid " << sid << std::endl;
      continue;
    }
#endif
    // critical edge vertices
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(sid, 0, v0);
    triangulation.getEdgeVertex(sid, 1, v1);

    // follow vpath from each vertex of the critical edge
    followVPath(v0);
    followVPath(v1);
    // Check for same end point (which would make an orbit if simplified)
    //  End point depends on the type (birth(type=0)/death(type=2))
    if(mins.size() >= 2 && mins[0].type == mins[1].type) {
      if(mins[0].type == 0 && mins[0].birth == mins[1].birth) {
        mins[0].generateOrbit = true;
        mins[1].generateOrbit = true;
      } else if(mins[0].type == 2 && mins[0].death == mins[1].death) {
        mins[0].generateOrbit = true;
        mins[1].generateOrbit = true;
      }
    }
  }

  return res;
}

template <typename dataType,
          typename triangulationType,
          typename GFS,
          typename GFSN,
          typename OB>
std::vector<std::vector<ttk::VectorSimplification::CandidatePair>>
  ttk::VectorSimplification::getSaddle2ToAscPair(
    const std::vector<SimplexId> &criticalCells,
    const GFS &getFaceStar,
    const GFSN &getFaceStarNumber,
    const OB &isOnBoundary,
    const triangulationType &triangulation,
    const dataType dummyVariable) const {
  // To eliminate unused variable warning
  TTK_FORCE_USE(dummyVariable);
  TTK_FORCE_USE(isOnBoundary);

  const auto dim = this->dcvf_.getDimensionality();
  std::vector<std::vector<CandidatePair>> res(criticalCells.size());

  // follow vpaths from 2-saddles to maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < criticalCells.size(); ++i) {
    const auto sid = criticalCells[i];
    auto &maxs = res[i];

    const auto followVPath = [this, dim, sid, &maxs,
                              &triangulation](const SimplexId v) {
      std::vector<Cell> vpath{};
      vpath.emplace_back(Cell{dim - 1, sid});
      int alternates
        = this->dcvf_.getAscendingPath<dataType, triangulationType>(
          Cell{dim, v}, vpath, triangulation, false);
      Cell &lastCell = vpath.back();
      if(lastCell.dim_ == dim && this->dcvf_.isCellCritical(lastCell)) {
        auto weight = this->dcvf_.getPersistence<dataType, triangulationType>(
          vpath, triangulation);
        maxs.emplace_back(sid, Cell{dim, v}, lastCell.id_, 2, weight);
        maxs.back().alternations = alternates;
      } else if(lastCell.dim_ == 0 && this->dcvf_.isCellCritical(lastCell)) {
        // Might change dimensions of end CP
        auto weight = this->dcvf_.getPersistence<dataType, triangulationType>(
          vpath, triangulation);
        maxs.emplace_back(lastCell.id_, Cell{dim, v}, sid, 0, weight);
        maxs.back().alternations = alternates;
      }
    };

    const auto starNumber = getFaceStarNumber(sid);

    for(SimplexId j = 0; j < starNumber; ++j) {
      SimplexId cellId{};
      getFaceStar(sid, j, cellId);
      followVPath(cellId);
    }

    // Not sure what this does other than add an extra saddle-max pair:
    if(!isOnBoundary(sid) && maxs.size() >= 2) {
      // Check for same end point (which would make an orbit if simplified)
      //  End point depends on the type (birth(type=0)/death(type=2))
      if(maxs.size() >= 2 && maxs[0].type == maxs[1].type) {
        if(maxs[0].type == 0 && maxs[0].birth == maxs[1].birth) {
          maxs[0].generateOrbit = true;
          maxs[1].generateOrbit = true;
        } else if(maxs[0].type == 2 && maxs[0].death == maxs[1].death) {
          maxs[0].generateOrbit = true;
          maxs[1].generateOrbit = true;
        }
      }
    }
  }

  return res;
}

template <typename dataType, typename triangulationType>
void ttk::VectorSimplification::getDescSaddlePairs(
  std::vector<CandidatePair> &pairs,
  const std::vector<SimplexId> &criticalEdges,
  const triangulationType &triangulation) const {

  auto saddle1ToMinima = getSaddle1ToDescPair<dataType, triangulationType>(
    criticalEdges, triangulation);

  for(size_t i = 0; i < saddle1ToMinima.size(); ++i) {
    auto &mins = saddle1ToMinima[i];
    // Add to pairs
    for(size_t j = 0; j < mins.size(); ++j) {
      pairs.emplace_back(mins[j]);
    }
  }
}

template <typename dataType, typename triangulationType>
void ttk::VectorSimplification::getAscSaddlePairs(
  std::vector<CandidatePair> &pairs,
  const std::vector<SimplexId> &criticalSaddles,
  const triangulationType &triangulation) const {

  const auto dim = this->dcvf_.getDimensionality();

  auto saddle2ToMaxima
    = dim == 3
        ? getSaddle2ToAscPair(
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
          triangulation, static_cast<dataType>(0.0))
        : getSaddle2ToAscPair(
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
          triangulation, static_cast<dataType>(0.0));

  for(size_t i = 0; i < saddle2ToMaxima.size(); ++i) {
    auto &maxs = saddle2ToMaxima[i];
    // Add to pairs
    for(size_t j = 0; j < maxs.size(); ++j) {
      pairs.emplace_back(maxs[j]);
    }
  }
}
