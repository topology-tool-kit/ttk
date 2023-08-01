/// \ingroup base
/// \class ttk::MergeTreePrincipalGeodesics
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021
///
/// This module defines the %MergeTreePrincipalGeodesics class that computes
/// Principal Geodesic Analysis on the space of merge trees or persistence
/// diagrams, that is, a set of orthogonal geodesic axes defining a basis with
/// the barycenter as origin.
///
/// \b Related \b publication: \n
/// "Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Jules Vidal, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2022
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <MergeTreePrincipalGeodesicsBase.h>

namespace ttk {

  /**
   * The MergeTreePrincipalGeodesics class provides methods to compute
   * Principal Geodesic Analysis on the space of merge trees or persistence
   * diagrams, that is, a set of orthogonal geodesic axes defining a basis with
   * the barycenter as origin.
   */
  class MergeTreePrincipalGeodesics : virtual public Debug,
                                      public MergeTreePrincipalGeodesicsBase {

  protected:
    bool doComputeReconstructionError_ = false;
    // TODO keepState works only when enabled before first computation
    bool keepState_ = false;
    unsigned int noProjectionStep_ = 2;

    // Advanced parameters
    bool projectInitializedVectors_ = true;

    // Old/Testing
    double t_vectorCopy_time_ = 0.0, t_allVectorCopy_time_ = 0.0;

    // Filled by the algorithm
    std::vector<double> inputToBaryDistances_;
    std::vector<std::vector<double>> inputToGeodesicsDistances_;
    ftm::MergeTree<double> barycenter_, barycenterInput2_, barycenterBDT_,
      barycenterInput2BDT_;
    bool barycenterWasComputed_ = false;

    int newVectorOffset_ = 0;
    double cumulVariance_ = 0.0, cumulTVariance_ = 0.0;

  public:
    MergeTreePrincipalGeodesics() {
      // inherited from Debug: prefix will be printed at the beginning of every
      // msg
      this->setDebugMsgPrefix("MergeTreePrincipalGeodesics");
#ifdef TTK_ENABLE_OPENMP
      omp_set_nested(1);
#endif
    }

    unsigned int getGeodesicNumber() {
      return vS_.size();
    }

    //----------------------------------------------------------------------------
    // OpenMP Reduction
    //----------------------------------------------------------------------------
    struct Compare {
      double bestDistance = std::numeric_limits<double>::max();
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> bestMatching,
        bestMatching2;
      int bestIndex = 0;
    };

    //----------------------------------------------------------------------------
    // Init
    //----------------------------------------------------------------------------
    template <class dataType>
    void initVectors(int geodesicNumber,
                     ftm::MergeTree<dataType> &barycenter,
                     std::vector<ftm::MergeTree<dataType>> &trees,
                     ftm::MergeTree<dataType> &barycenter2,
                     std::vector<ftm::MergeTree<dataType>> &trees2,
                     std::vector<std::vector<double>> &v1,
                     std::vector<std::vector<double>> &v2,
                     std::vector<std::vector<double>> &trees2V1,
                     std::vector<std::vector<double>> &trees2V2) {
      auto initializedVectorsProjection
        = [=](int _geodesicNumber, ftm::MergeTree<dataType> &_barycenter,
              std::vector<std::vector<double>> &_v,
              std::vector<std::vector<double>> &_v2,
              std::vector<std::vector<std::vector<double>>> &_vS,
              std::vector<std::vector<std::vector<double>>> &_v2s,
              ftm::MergeTree<dataType> &_barycenter2,
              std::vector<std::vector<double>> &_trees2V,
              std::vector<std::vector<double>> &_trees2V2,
              std::vector<std::vector<std::vector<double>>> &_trees2Vs,
              std::vector<std::vector<std::vector<double>>> &_trees2V2s,
              bool _useSecondInput, unsigned int _noProjectionStep) {
            return this->projectionStep(_geodesicNumber, _barycenter, _v, _v2,
                                        _vS, _v2s, _barycenter2, _trees2V,
                                        _trees2V2, _trees2Vs, _trees2V2s,
                                        _useSecondInput, _noProjectionStep);
          };

      MergeTreeAxesAlgorithmBase::initVectors<dataType>(
        geodesicNumber, barycenter, trees, barycenter2, trees2, v1, v2,
        trees2V1, trees2V2, newVectorOffset_, inputToBaryDistances_,
        baryMatchings_, baryMatchings2_, inputToGeodesicsDistances_, vS_, v2s_,
        trees2Vs_, trees2V2s_, projectInitializedVectors_,
        initializedVectorsProjection);
    }

    //----------------------------------------------------------------------------
    // Costs
    //----------------------------------------------------------------------------
    double orthogonalCost(std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s,
                          std::vector<std::vector<double>> &v,
                          std::vector<std::vector<double>> &v2) {
      return verifyOrthogonality(vS, v2s, v, v2, false);
    }

    double regularizerCost(std::vector<std::vector<double>> &v,
                           std::vector<std::vector<double>> &v2) {
      auto cost = ttk::Geometry::dotProductFlatten(v, v2)
                  - ttk::Geometry::magnitudeFlatten(v)
                      * ttk::Geometry::magnitudeFlatten(v2);
      return cost * cost;
    }

    double projectionCost(std::vector<std::vector<double>> &v,
                          std::vector<std::vector<double>> &v2,
                          std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s,
                          double optMapCost) {
      return regularizerCost(v, v2) + orthogonalCost(vS, v2s, v, v2)
             + optMapCost;
    }

    //----------------------------------------------------------------------------
    // Projection
    //----------------------------------------------------------------------------
    template <class dataType>
    double barycentricProjection(ftm::MergeTree<dataType> &barycenter,
                                 ftm::MergeTree<dataType> &extremity,
                                 std::vector<std::vector<double>> &v,
                                 bool isV1,
                                 bool useDoubleInput = false,
                                 bool isFirstInput = true) {
      ftm::FTMTree_MT *barycenterTree = &(barycenter.tree);
      ftm::FTMTree_MT *extremityTree = &(extremity.tree);
      double const t = (isV1 ? -1.0 : 1.0);

      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      dataType distance;
      std::vector<ftm::idNode> matchingVector;

      if(extremityTree->getRealNumberOfNodes() != 0) {
        computeOneDistance(barycenter, extremity, matching, distance, true,
                           useDoubleInput, isFirstInput);
        getMatchingVector(barycenter, extremity, matching, matchingVector);
      } else
        matchingVector.resize(barycenterTree->getNumberOfNodes(),
                              std::numeric_limits<ftm::idNode>::max());

      std::vector<std::vector<double>> const oriV = v;
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        auto matched = matchingVector[i];
        auto birthDeathBary
          = getParametrizedBirthDeath<dataType>(barycenterTree, i);
        dataType birthBary = std::get<0>(birthDeathBary);
        dataType deathBary = std::get<1>(birthDeathBary);
        std::vector<double> newV{0.0, 0.0};
        if(matched != std::numeric_limits<ftm::idNode>::max()) {
          auto birthDeathMatched
            = getParametrizedBirthDeath<dataType>(extremityTree, matched);
          newV[0] = std::get<0>(birthDeathMatched);
          newV[1] = std::get<1>(birthDeathMatched);
        } else {
          dataType projec = (birthBary + deathBary) / 2.0;
          newV[0] = projec;
          newV[1] = projec;
        }
        newV[0] = (newV[0] - birthBary) * t;
        newV[1] = (newV[1] - deathBary) * t;
        v[i] = newV;
      }

      // Compute distance between old and new extremity
      double const cost = ttk::Geometry::distanceFlatten(v, oriV);
      return cost;
    }

    template <class dataType>
    double
      optimalMappingSetProjection(ftm::MergeTree<dataType> &barycenter,
                                  std::vector<std::vector<double>> &v,
                                  std::vector<std::vector<double>> &v2,
                                  ftm::MergeTree<dataType> &barycenter2,
                                  std::vector<std::vector<double>> &trees2V,
                                  std::vector<std::vector<double>> &trees2V2,
                                  bool useSecondInput = false) {
      std::vector<ftm::MergeTree<dataType>> extremities(
        (useSecondInput ? 4 : 2));
      getInterpolation<dataType>(barycenter, v, v2, 0.0, extremities[0]);
      getInterpolation<dataType>(barycenter, v, v2, 1.0, extremities[1]);
      if(useSecondInput) {
        getInterpolation<dataType>(
          barycenter2, trees2V, trees2V2, 0.0, extremities[2]);
        getInterpolation<dataType>(
          barycenter2, trees2V, trees2V2, 1.0, extremities[3]);
      }
      double cost = 0.0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < extremities.size(); ++i) {
        bool isFirstInput = (i < 2);
        ftm::MergeTree<dataType> &baryToUse
          = (i < 2 ? barycenter : barycenter2);
        std::vector<std::vector<double>> &vToUse
          = (i == 0 ? v : (i == 1 ? v2 : (i == 2 ? trees2V : trees2V2)));
        cost
          += barycentricProjection(baryToUse, extremities[i], vToUse,
                                   (i % 2 == 0), useSecondInput, isFirstInput);
      }
      return cost;
    }

    // Collinearity constraint
    void
      trueGeneralizedGeodesicProjection(std::vector<std::vector<double>> &v1,
                                        std::vector<std::vector<double>> &v2) {
      std::vector<double> v1_flatten, v2_flatten;
      ttk::Geometry::flattenMultiDimensionalVector(v1, v1_flatten);
      ttk::Geometry::flattenMultiDimensionalVector(v2, v2_flatten);
      double const v1_norm = ttk::Geometry::magnitude(v1_flatten);
      double const v2_norm = ttk::Geometry::magnitude(v2_flatten);
      double const beta = v2_norm / (v1_norm + v2_norm);
      std::vector<double> v;
      ttk::Geometry::addVectors(v1_flatten, v2_flatten, v);
      ttk::Geometry::scaleVector(v, (1 - beta), v1_flatten);
      ttk::Geometry::scaleVector(v, beta, v2_flatten);
      ttk::Geometry::unflattenMultiDimensionalVector(v1_flatten, v1, 2);
      ttk::Geometry::unflattenMultiDimensionalVector(v2_flatten, v2, 2);
    }

    void
      orthogonalProjection(std::vector<std::vector<double>> &v1,
                           std::vector<std::vector<double>> &v2,
                           std::vector<std::vector<std::vector<double>>> &vS,
                           std::vector<std::vector<std::vector<double>>> &v2s) {
      // Multi flatten and sum vS and v2s
      std::vector<std::vector<double>> sumVs;
      ttk::Geometry::multiAddVectorsFlatten(vS, v2s, sumVs);

      // Flatten v1 and v2
      std::vector<double> v1_flatten, v2_flatten, v1_proj, v2_proj;
      ttk::Geometry::flattenMultiDimensionalVector(v1, v1_flatten);
      ttk::Geometry::flattenMultiDimensionalVector(v2, v2_flatten);

      // Call Gram Schmidt
      callGramSchmidt(sumVs, v1_flatten, v1_proj);
      callGramSchmidt(sumVs, v2_flatten, v2_proj);

      // Unflatten the resulting vectors
      ttk::Geometry::unflattenMultiDimensionalVector(v1_proj, v1, 2);
      ttk::Geometry::unflattenMultiDimensionalVector(v2_proj, v2, 2);
    }

    // TODO avoid copying vectors
    template <class dataType>
    double
      projectionStep(int geodesicNumber,
                     ftm::MergeTree<dataType> &barycenter,
                     std::vector<std::vector<double>> &v,
                     std::vector<std::vector<double>> &v2,
                     std::vector<std::vector<std::vector<double>>> &vS,
                     std::vector<std::vector<std::vector<double>>> &v2s,
                     ftm::MergeTree<dataType> &barycenter2,
                     std::vector<std::vector<double>> &trees2V,
                     std::vector<std::vector<double>> &trees2V2,
                     std::vector<std::vector<std::vector<double>>> &trees2Vs,
                     std::vector<std::vector<std::vector<double>>> &trees2V2s,
                     bool useSecondInput,
                     unsigned int noProjectionStep) {
      std::vector<std::vector<std::vector<double>>> vSConcat, v2sConcat;
      if(useSecondInput) {
        Timer t_vectorCopy;
        vSConcat = vS;
        v2sConcat = v2s;
        for(unsigned int j = 0; j < vS.size(); ++j) {
          vSConcat[j].insert(
            vSConcat[j].end(), trees2Vs[j].begin(), trees2Vs[j].end());
          v2sConcat[j].insert(
            v2sConcat[j].end(), trees2V2s[j].begin(), trees2V2s[j].end());
        }
        t_vectorCopy_time_ += t_vectorCopy.getElapsedTime();
      }

      double optMapCost = 0.0;
      for(unsigned i = 0; i < noProjectionStep; ++i) {
        std::vector<std::vector<double>> vOld = v, v2Old = v2;

        // --- Optimal mapping set projecton
        printMsg("OMS Proj.", 0, 0, threadNumber_, debug::LineMode::REPLACE,
                 debug::Priority::DETAIL);
        Timer t_optMap;
        optMapCost = optimalMappingSetProjection(
          barycenter, v, v2, barycenter2, trees2V, trees2V2, useSecondInput);
        printMsg("OMS Proj.", 1, t_optMap.getElapsedTime(), threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

        // ---
        std::vector<std::vector<double>> vConcat, v2Concat;
        if(useSecondInput) {
          Timer t_vectorCopy;
          vConcat = v;
          vConcat.insert(vConcat.end(), trees2V.begin(), trees2V.end());
          v2Concat = v2;
          v2Concat.insert(v2Concat.end(), trees2V2.begin(), trees2V2.end());
          t_vectorCopy_time_ += t_vectorCopy.getElapsedTime();
        }

        // --- True generalized geodesic projection
        printMsg("TGG Proj.", 0, 0, threadNumber_, debug::LineMode::REPLACE,
                 debug::Priority::DETAIL);
        Timer t_trueGeod;
        if(useSecondInput)
          trueGeneralizedGeodesicProjection(vConcat, v2Concat);
        else
          trueGeneralizedGeodesicProjection(v, v2);
        printMsg("TGG Proj.", 1, t_trueGeod.getElapsedTime(), threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

        // --- Orthogonal projection
        if(geodesicNumber != 0) {
          printMsg("Orth. Proj.", 0, 0, threadNumber_, debug::LineMode::REPLACE,
                   debug::Priority::DETAIL);
          Timer t_ortho;
          if(useSecondInput) {
            orthogonalProjection(vConcat, v2Concat, vSConcat, v2sConcat);
          } else
            orthogonalProjection(v, v2, vS, v2s);
          printMsg("Orth. Proj.", 1, t_ortho.getElapsedTime(), threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
        }

        // ---
        if(useSecondInput) {
          Timer t_vectorCopy;
          for(unsigned int j = 0; j < v.size(); ++j) {
            v[j] = vConcat[j];
            v2[j] = v2Concat[j];
          }
          for(unsigned int j = 0; j < trees2V.size(); ++j) {
            trees2V[j] = vConcat[v.size() + j];
            trees2V2[j] = v2Concat[v.size() + j];
          }
          t_vectorCopy_time_ += t_vectorCopy.getElapsedTime();
        }
      }
      return optMapCost;
    }

    //----------------------------------------------------------------------------
    // Assignment
    //----------------------------------------------------------------------------
    template <class dataType>
    void assignmentImpl(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      ftm::MergeTree<dataType> &barycenter2,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<std::vector<double>> &trees2V,
      std::vector<std::vector<double>> &trees2V2,
      std::vector<std::vector<double>> &allTreesTs,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<std::vector<double>>> &trees2Vs,
      std::vector<std::vector<std::vector<double>>> &trees2V2s,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      std::vector<double> &ts,
      std::vector<double> &distances) {
      std::vector<std::vector<Compare>> best(
        trees.size(), std::vector<Compare>(k_));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) if(parallelize_) \
  shared(best)
      {
#pragma omp single nowait
        {
#endif
          for(unsigned int k = 0; k < k_; ++k) {
            for(unsigned int i = 0; i < trees.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task shared(best) firstprivate(i, k)
              {
#endif
                double const kT = (k % 2 == 0 ? k / 2 : k_ - 1 - (int)(k / 2));
                double t = 1.0 / (k_ - 1) * kT;

                dataType distance, distance2;
                ftm::MergeTree<dataType> interpolated;
                auto tsToUse = (allTreesTs.size() == 0 ? std::vector<double>()
                                                       : allTreesTs[i]);
                getMultiInterpolation(
                  barycenter, vS, v2s, v, v2, tsToUse, t, interpolated);
                if(interpolated.tree.getRealNumberOfNodes() != 0) {
                  computeOneDistance<dataType>(interpolated, trees[i],
                                               best[i][kT].bestMatching,
                                               distance, true, useDoubleInput_);
                  if(trees2.size() != 0) {
                    ftm::MergeTree<dataType> interpolated2;
                    getMultiInterpolation(barycenter2, trees2Vs, trees2V2s,
                                          trees2V, trees2V2, tsToUse, t,
                                          interpolated2);
                    computeOneDistance<dataType>(
                      interpolated2, trees2[i], best[i][kT].bestMatching2,
                      distance2, true, useDoubleInput_, false);
                    distance = mixDistances(distance, distance2);
                  }
                  best[i][kT].bestDistance = distance;
                  best[i][kT].bestIndex = kT;
                }
#ifdef TTK_ENABLE_OPENMP
              } // pragma omp task
#endif
            }
          }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

          // Reduction
          for(unsigned int i = 0; i < trees.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i)
            {
#endif
              double bestDistance = std::numeric_limits<double>::max();
              int bestIndex = 0;
              for(unsigned int k = 0; k < k_; ++k) {
                if(best[i][k].bestDistance < bestDistance) {
                  bestIndex = k;
                  bestDistance = best[i][k].bestDistance;
                }
              }
              matchings[i] = best[i][bestIndex].bestMatching;
              if(trees2.size() != 0)
                matchings2[i] = best[i][bestIndex].bestMatching2;
              ts[i] = best[i][bestIndex].bestIndex * 1.0 / (k_ - 1);
              distances[i] = best[i][bestIndex].bestDistance;
#ifdef TTK_ENABLE_OPENMP
            } // pragma omp task
#endif
          }
#ifdef TTK_ENABLE_OPENMP
        } // pragma omp single nowait
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void assignmentStep(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      ftm::MergeTree<dataType> &barycenter2,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<std::vector<double>> &trees2V,
      std::vector<std::vector<double>> &trees2V2,
      std::vector<std::vector<double>> &allTreesTs,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<std::vector<double>>> &trees2Vs,
      std::vector<std::vector<std::vector<double>>> &trees2V2s,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      std::vector<double> &ts,
      std::vector<double> &distances) {

      matchings.resize(trees.size());
      matchings2.resize(trees2.size());
      ts.resize(trees.size());
      distances.resize(trees.size());

      // Assignment
      assignmentImpl<dataType>(barycenter, trees, v, v2, barycenter2, trees2,
                               trees2V, trees2V2, allTreesTs, vS, v2s, trees2Vs,
                               trees2V2s, matchings, matchings2, ts, distances);
    }

    //----------------------------------------------------------------------------
    // Update
    //----------------------------------------------------------------------------
    template <class dataType>
    void updateClosedForm(
      int geodesicNumber,
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<ftm::MergeTree<dataType>> &allInterpolated,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<double>> &tss,
      std::vector<std::vector<double>> &vR,
      std::vector<std::vector<double>> &vR2,
      std::vector<bool> &isUniform) {

      // Init
      ftm::FTMTree_MT *barycenterTree = &(barycenter.tree);
      std::vector<ftm::FTMTree_MT *> ftmTrees, allInterpolatedTrees;
      ttk::ftm::mergeTreeToFTMTree<dataType>(trees, ftmTrees);
      if(geodesicNumber != 0) {
        ttk::ftm::mergeTreeToFTMTree<dataType>(
          allInterpolated, allInterpolatedTrees);
      }

      // Get matching matrix
      std::vector<std::vector<ftm::idNode>> matchingMatrix;
      getMatchingMatrix(barycenter, trees, matchings, matchingMatrix);

      // Update
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;

        // Verify that ts is not uniform
        if(isUniform[i]) {
          v[i] = vR[i];
          v2[i] = vR2[i];
          continue;
        }

        // Compute projection
        auto birthDeathBary
          = getParametrizedBirthDeath<dataType>(barycenterTree, i);
        dataType birthBary = std::get<0>(birthDeathBary);
        dataType deathBary = std::get<1>(birthDeathBary);
        dataType projec = (birthBary + deathBary) / 2.0;
        std::vector<dataType> allBirthBary(trees.size(), birthBary);
        std::vector<dataType> allDeathBary(trees.size(), deathBary);
        std::vector<dataType> allProjec(trees.size(), projec);

        if(geodesicNumber != 0)
          for(unsigned int j = 0; j < trees.size(); ++j) {
            auto birthDeathInterpol
              = getParametrizedBirthDeath<dataType>(allInterpolatedTrees[j], i);
            allBirthBary[j] = std::get<0>(birthDeathInterpol);
            allDeathBary[j] = std::get<1>(birthDeathInterpol);
            allProjec[j] = (allBirthBary[j] + allDeathBary[j]) / 2.0;
          }

        // Compute all matched values
        std::vector<std::vector<dataType>> allMatched(trees.size());
        for(unsigned int j = 0; j < trees.size(); ++j) {
          dataType birth = allProjec[j];
          dataType death = allProjec[j];
          if(matchingMatrix[i][j] != std::numeric_limits<ftm::idNode>::max()) {
            auto birthDeath = getParametrizedBirthDeath<dataType>(
              ftmTrees[j], matchingMatrix[i][j]);
            birth = std::get<0>(birthDeath);
            death = std::get<1>(birthDeath);
          }
          allMatched[j].resize(2);
          allMatched[j][0] = birth;
          allMatched[j][1] = death;
        }

        // Compute general terms
        double ti_squared = 0.0, one_min_ti_squared = 0.0, ti_one_min_ti = 0.0;
        for(auto t : tss[i]) {
          ti_squared += t * t;
          one_min_ti_squared += (1 - t) * (1 - t);
          ti_one_min_ti += t * (1 - t);
        }

        // Compute multiplier
        double multBirthV1 = 0.0, multDeathV1 = 0.0, multBirthV2 = 0.0,
               multDeathV2 = 0.0;
        for(unsigned int j = 0; j < trees.size(); ++j) {
          multBirthV1
            += tss[i][j] * (allMatched[j][0] - allBirthBary[j]) / ti_squared;
          multDeathV1
            += tss[i][j] * (allMatched[j][1] - allDeathBary[j]) / ti_squared;
          multBirthV2 += (1 - tss[i][j]) * (-allMatched[j][0] + allBirthBary[j])
                         / one_min_ti_squared;
          multDeathV2 += (1 - tss[i][j]) * (-allMatched[j][1] + allDeathBary[j])
                         / one_min_ti_squared;
        }

        // Compute new birth death
        double newBirthV1 = 0.0, newDeathV1 = 0.0, newBirthV2 = 0.0,
               newDeathV2 = 0.0;
        for(unsigned int j = 0; j < trees.size(); ++j) {
          newBirthV1 += (1 - tss[i][j])
                        * (-allMatched[j][0] + allBirthBary[j]
                           + tss[i][j] * multBirthV1);
          newDeathV1 += (1 - tss[i][j])
                        * (-allMatched[j][1] + allDeathBary[j]
                           + tss[i][j] * multDeathV1);
          newBirthV2 += tss[i][j]
                        * (allMatched[j][0] - allBirthBary[j]
                           + (1 - tss[i][j]) * multBirthV2);
          newDeathV2 += tss[i][j]
                        * (allMatched[j][1] - allDeathBary[j]
                           + (1 - tss[i][j]) * multDeathV2);
        }
        double const divisorV1
          = one_min_ti_squared - ti_one_min_ti * ti_one_min_ti / ti_squared;
        double const divisorV2
          = ti_squared - ti_one_min_ti * ti_one_min_ti / one_min_ti_squared;
        newBirthV1 /= divisorV1;
        newDeathV1 /= divisorV1;
        newBirthV2 /= divisorV2;
        newDeathV2 /= divisorV2;

        // Update vectors
        v[i][0] = newBirthV1;
        v[i][1] = newDeathV1;
        v2[i][0] = newBirthV2;
        v2[i][1] = newDeathV2;
      }
    }

    template <class dataType>
    void
      manageIndividualTs(int geodesicNumber,
                         ftm::MergeTree<dataType> &barycenter,
                         std::vector<ftm::MergeTree<dataType>> &trees,
                         std::vector<std::vector<double>> &v,
                         std::vector<std::vector<double>> &v2,
                         std::vector<std::vector<std::vector<double>>> &vS,
                         std::vector<std::vector<std::vector<double>>> &v2s,
                         std::vector<double> &ts,
                         std::vector<std::vector<double>> &allTreesTs,
                         std::vector<ftm::MergeTree<dataType>> &allInterpolated,
                         std::vector<bool> &isUniform,
                         std::vector<std::vector<double>> &tss,
                         unsigned int &noUniform,
                         bool &foundAllUniform) {
      // Get multi interpolation
      allInterpolated.resize(trees.size());
      if(geodesicNumber != 0) {
        for(unsigned int i = 0; i < trees.size(); ++i)
          getMultiInterpolation(
            barycenter, vS, v2s, allTreesTs[i], allInterpolated[i]);
      }

      // Manage individuals t
      noUniform = 0;
      foundAllUniform = true;
      isUniform.resize(barycenter.tree.getNumberOfNodes(), false);
      tss.resize(barycenter.tree.getNumberOfNodes());
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        tss[i] = ts;
        for(unsigned int j = 0; j < tss[i].size(); ++j) {
          auto &treeToUse
            = (geodesicNumber != 0 ? allInterpolated[j] : barycenter);
          tss[i][j] = getTNew<dataType>(treeToUse, v, v2, i, ts[j]);
        }
        isUniform[i] = ttk::Geometry::isVectorUniform(tss[i]);
        noUniform += isUniform[i];
        foundAllUniform &= isUniform[i];
      }
    }

    template <class dataType>
    bool updateClosedFormStep(
      int geodesicNumber,
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      ftm::MergeTree<dataType> &barycenter2,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<std::vector<double>> &trees2V,
      std::vector<std::vector<double>> &trees2V2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      std::vector<std::vector<std::vector<double>>> &trees2Vs,
      std::vector<std::vector<std::vector<double>>> &trees2V2s,
      std::vector<double> &ts,
      std::vector<std::vector<double>> &allTreesTs) {

      std::vector<ftm::MergeTree<dataType>> allInterpolated, allInterpolated2;
      std::vector<bool> isUniform, isUniform2;
      std::vector<std::vector<double>> tss, tss2;
      unsigned int noUniform;
      bool foundAllUniform;
      manageIndividualTs(geodesicNumber, barycenter, trees, v, v2, vS, v2s, ts,
                         allTreesTs, allInterpolated, isUniform, tss, noUniform,
                         foundAllUniform);
      if(trees2.size() != 0) {
        unsigned int noUniform2;
        bool foundAllUniform2;
        manageIndividualTs(geodesicNumber, barycenter2, trees2, trees2V,
                           trees2V2, trees2Vs, trees2V2s, ts, allTreesTs,
                           allInterpolated2, isUniform2, tss2, noUniform2,
                           foundAllUniform2);
        noUniform += noUniform2;
        foundAllUniform &= foundAllUniform2;
      }

      if(foundAllUniform) {
        printMsg("All projection coefficients are the same.");
        printMsg("New vectors will be initialized.");
        newVectorOffset_ += 1;
        initVectors(geodesicNumber, barycenter, trees, barycenter2, trees2, v,
                    v2, trees2V, trees2V2);
        return true;
      }
      std::vector<std::vector<double>> vR, vR2, trees2VR, trees2VR2;
      if(noUniform != 0) {
        printMsg("Found " + std::to_string(noUniform)
                 + " uniform coefficients.");
        initVectors(geodesicNumber, barycenter, trees, barycenter2, trees2, vR,
                    vR2, trees2VR, trees2VR2);
      }

      updateClosedForm(geodesicNumber, barycenter, trees, allInterpolated, v,
                       v2, matchings, tss, vR, vR2, isUniform);
      if(trees2.size() != 0) {
        updateClosedForm(geodesicNumber, barycenter2, trees2, allInterpolated2,
                         trees2V, trees2V2, matchings2, tss2, trees2VR,
                         trees2VR2, isUniform2);
        copyMinMaxPairVector(v, v2, trees2V, trees2V2);
      }
      return false;
    }

    template <class dataType>
    bool updateStep(
      int geodesicNumber,
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      ftm::MergeTree<dataType> &barycenter2,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<std::vector<double>> &trees2V,
      std::vector<std::vector<double>> &trees2V2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      std::vector<std::vector<std::vector<double>>> &trees2Vs,
      std::vector<std::vector<std::vector<double>>> &trees2V2s,
      std::vector<double> &ts,
      std::vector<std::vector<double>> &allTreesTs) {
      return updateClosedFormStep<dataType>(
        geodesicNumber, barycenter, trees, v, v2, matchings, vS, v2s,
        barycenter2, trees2, trees2V, trees2V2, matchings2, trees2Vs, trees2V2s,
        ts, allTreesTs);
    }

    //----------------------------------------------------------------------------
    // Main functions
    //----------------------------------------------------------------------------
    template <class dataType>
    bool convergenceStep(std::vector<double> &distances,
                         std::vector<std::vector<double>> &v,
                         std::vector<std::vector<double>> &v2,
                         dataType &oldFrechetEnergy,
                         dataType &minFrechetEnergy,
                         int &cptBlocked,
                         bool &converged,
                         double optMapCost) {
      bool isBestEnergy = false;

      // Reconstruction cost (main energy)
      double frechetEnergy = 0;
      for(unsigned int i = 0; i < distances.size(); ++i)
        frechetEnergy += distances[i] * distances[i] / distances.size();
      std::stringstream ssEnergy;
      ssEnergy << "Energy      = " << frechetEnergy;
      printMsg(ssEnergy.str());

      // Prop. cost
      std::stringstream ssReg;
      auto reg = regularizerCost(v, v2);
      ssReg << "Prop. cost  = " << reg;
      printMsg(ssReg.str());

      // Ortho. cost
      std::stringstream ssOrthoCost;
      auto orthoCost = orthogonalCost(vS_, v2s_, v, v2);
      ssOrthoCost << "Ortho. cost = " << orthoCost;
      printMsg(ssOrthoCost.str());

      // Map. cost
      std::stringstream ssOptMapCost;
      ssOptMapCost << "Map. cost   = " << optMapCost;
      printMsg(ssOptMapCost.str());

      // Detect convergence
      double tol = 0.01;
      tol = oldFrechetEnergy / 125.0;
      converged = std::abs(frechetEnergy - oldFrechetEnergy) < tol;
      oldFrechetEnergy = frechetEnergy;

      if(frechetEnergy + ENERGY_COMPARISON_TOLERANCE < minFrechetEnergy) {
        minFrechetEnergy = frechetEnergy;
        cptBlocked = 0;
        isBestEnergy = true;
      }
      if(not converged) {
        cptBlocked += (minFrechetEnergy < frechetEnergy) ? 1 : 0;
        converged = (cptBlocked >= 10);
      }

      return isBestEnergy;
    }

    template <class dataType>
    void
      computePrincipalGeodesic(unsigned int geodesicNumber,
                               ftm::MergeTree<dataType> &barycenter,
                               std::vector<ftm::MergeTree<dataType>> &trees,
                               ftm::MergeTree<dataType> &barycenter2,
                               std::vector<ftm::MergeTree<dataType>> &trees2) {
      // ----- Init Parameters
      printMsg("Init", 0, 0, threadNumber_, debug::LineMode::REPLACE);
      Timer t_init;
      std::vector<std::vector<double>> v, v2, trees2V, trees2V2;
      initVectors<dataType>(geodesicNumber, barycenter, trees, barycenter2,
                            trees2, v, v2, trees2V, trees2V2);
      newVectorOffset_ = 0;
      printMsg("Init", 1, t_init.getElapsedTime(), threadNumber_);

      std::vector<std::vector<double>> bestV, bestV2, bestTrees2V, bestTrees2V2;
      std::vector<double> bestTs, bestDistances;
      int bestIteration = 0;

      // ----- Init Loop
      dataType oldFrechetEnergy, minFrechetEnergy;
      int cptBlocked, iteration = 0;
      auto initLoop = [&]() {
        oldFrechetEnergy = -1;
        minFrechetEnergy = std::numeric_limits<dataType>::max();
        cptBlocked = 0;
        iteration = 0;
      };
      initLoop();

      // ----- Algorithm
      double optMapCost = 0.0;
      bool converged = false;
      while(not converged) {
        std::stringstream ss;
        ss << "Iteration " << iteration;
        printMsg(debug::Separator::L2);
        printMsg(ss.str());

        // --- Assignment
        printMsg("Assignment", 0, 0, threadNumber_, debug::LineMode::REPLACE);
        Timer t_assignment;
        std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
          matchings, matchings2;
        std::vector<double> ts, distances;
        assignmentStep(barycenter, trees, v, v2, barycenter2, trees2, trees2V,
                       trees2V2, allTreesTs_, vS_, v2s_, trees2Vs_, trees2V2s_,
                       matchings, matchings2, ts, distances);
        inputToGeodesicsDistances_[geodesicNumber] = distances;
        allTs_[geodesicNumber] = ts;
        printMsg("Assignment", 1, t_assignment.getElapsedTime(), threadNumber_);

        // --- Convergence
        bool isBest = convergenceStep(distances, v, v2, oldFrechetEnergy,
                                      minFrechetEnergy, cptBlocked, converged,
                                      optMapCost);
        if(isBest) {
          Timer t_copy;
          bestV = v;
          bestV2 = v2;
          bestTrees2V = trees2V;
          bestTrees2V2 = trees2V2;
          bestTs = ts;
          bestDistances = distances;
          bestIteration = iteration;
          t_allVectorCopy_time_ += t_copy.getElapsedTime();
        }
        if(converged)
          break;

        // --- Update
        printMsg("Update", 0, 0, threadNumber_, debug::LineMode::REPLACE);
        Timer t_update;
        bool reset
          = updateStep(geodesicNumber, barycenter, trees, v, v2, matchings, vS_,
                       v2s_, barycenter2, trees2, trees2V, trees2V2, matchings2,
                       trees2Vs_, trees2V2s_, ts, allTreesTs_);
        printMsg("Update", 1, t_update.getElapsedTime(), threadNumber_);
        if(reset) {
          initLoop();
          continue;
        }

        // --- Projection
        printMsg("Projection", 0, 0, threadNumber_, debug::LineMode::REPLACE);
        Timer t_projection;
        optMapCost
          = projectionStep(geodesicNumber, barycenter, v, v2, vS_, v2s_,
                           barycenter2, trees2V, trees2V2, trees2Vs_,
                           trees2V2s_, (trees2.size() != 0), noProjectionStep_);
        auto projectionTime
          = t_projection.getElapsedTime() - t_vectorCopy_time_;
        t_allVectorCopy_time_ += t_vectorCopy_time_;
        t_vectorCopy_time_ = 0.0;
        printMsg("Projection", 1, projectionTime, threadNumber_);

        ++iteration;
      }
      printMsg(debug::Separator::L2);
      printMsg("Best energy is " + std::to_string(minFrechetEnergy)
               + " (iteration " + std::to_string(bestIteration) + " / "
               + std::to_string(iteration) + ")");
      printMsg(debug::Separator::L2);

      Timer t_copy;
      v = bestV;
      v2 = bestV2;
      trees2V = bestTrees2V;
      trees2V2 = bestTrees2V2;
      inputToGeodesicsDistances_[geodesicNumber] = bestDistances;
      allTs_[geodesicNumber] = bestTs;

      vS_.push_back(v);
      v2s_.push_back(v2);
      trees2Vs_.push_back(trees2V);
      trees2V2s_.push_back(trees2V2);
      ttk::Geometry::transposeMatrix(allTs_, allTreesTs_);
      t_allVectorCopy_time_ += t_copy.getElapsedTime();
    }

    template <class dataType>
    void
      computePrincipalGeodesics(std::vector<ftm::MergeTree<dataType>> &trees,
                                std::vector<ftm::MergeTree<dataType>> &trees2) {
      // --- Compute barycenter
      ftm::MergeTree<dataType> barycenter, barycenter2;
      if(not keepState_ or not barycenterWasComputed_) {
        Timer t_barycenter;
        printMsg("Barycenter", 0, t_barycenter.getElapsedTime(), threadNumber_,
                 debug::LineMode::REPLACE);
        computeOneBarycenter<dataType>(trees, barycenter, baryMatchings_,
                                       inputToBaryDistances_, useDoubleInput_);
        mergeTreeTemplateToDouble(barycenter, barycenterBDT_);
        if(trees2.size() != 0) {
          std::vector<double> baryDistances2;
          computeOneBarycenter<dataType>(trees2, barycenter2, baryMatchings2_,
                                         baryDistances2, useDoubleInput_,
                                         false);
          mergeTreeTemplateToDouble(barycenter2, barycenterInput2BDT_);
          for(unsigned int i = 0; i < inputToBaryDistances_.size(); ++i)
            inputToBaryDistances_[i]
              = mixDistances(inputToBaryDistances_[i], baryDistances2[i]);

          verifyMinMaxPair(barycenter, barycenter2);
        }
        printMsg("Barycenter", 1, t_barycenter.getElapsedTime(), threadNumber_);
        barycenterWasComputed_ = true;
      } else {
        printMsg("KeepState is enabled and barycenter was already computed");
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
          barycenterBDT_, barycenter);
        if(trees2.size() != 0) {
          ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
            barycenterInput2BDT_, barycenter2);
        }
      }
      printMsg(barycenter.tree.printTreeStats().str());
      mergeTreeTemplateToDouble(barycenter, barycenter_);
      if(trees2.size() != 0) {
        printMsg(barycenter2.tree.printTreeStats().str());
        mergeTreeTemplateToDouble(barycenter2, barycenterInput2_);
      }

      // --- Compute global variance
      double const globalVariance
        = computeVarianceFromDistances(inputToBaryDistances_);

      // --- Manage maximum number of geodesics
      unsigned int maxNoGeodesics = barycenter.tree.getRealNumberOfNodes() * 2;
      if(trees2.size() != 0)
        maxNoGeodesics += barycenter2.tree.getRealNumberOfNodes() * 2;
      if(maxNoGeodesics < numberOfGeodesics_) {
        std::stringstream ss;
        ss << numberOfGeodesics_ << " principal geodesics are asked but only "
           << maxNoGeodesics << " can be computed.";
        printMsg(ss.str());
        printMsg("(the maximum is twice the number of persistence pairs in the "
                 "barycenter)");
        numberOfGeodesics_ = maxNoGeodesics;
      }

      // --- Init
      unsigned int const oldNoGeod = allTs_.size();
      if(not keepState_) {
        allTs_.resize(numberOfGeodesics_, std::vector<double>(trees.size()));
        inputToGeodesicsDistances_.resize(
          numberOfGeodesics_, std::vector<double>(trees.size()));
        vS_.clear();
        v2s_.clear();
        trees2Vs_.clear();
        trees2V2s_.clear();
        allTreesTs_.clear();
        srand(deterministic_ ? 7 : time(nullptr));
        t_vectorCopy_time_ = 0.0;
        t_allVectorCopy_time_ = 0.0;
        cumulVariance_ = 0.0;
        cumulTVariance_ = 0.0;
      } else {
        allTs_.resize(numberOfGeodesics_, std::vector<double>(trees.size()));
        ttk::Geometry::transposeMatrix(allTs_, allTreesTs_);
        inputToGeodesicsDistances_.resize(
          numberOfGeodesics_, std::vector<double>(trees.size()));
        if(oldNoGeod != 0)
          printMsg(
            "KeepState is enabled, restart the computation at geodesic number "
            + std::to_string(oldNoGeod));
      }

      // --- Compute each geodesic
      for(unsigned int geodNum = oldNoGeod; geodNum < numberOfGeodesics_;
          ++geodNum) {
        printMsg(debug::Separator::L1);
        std::stringstream ss;
        ss << "Compute geodesic " << geodNum;
        printMsg(ss.str());

        // - Compute geodesic
        computePrincipalGeodesic<dataType>(
          geodNum, barycenter, trees, barycenter2, trees2);

        // - Compute explained variance
        printIterationVariances(
          barycenter, trees, barycenter2, trees2, geodNum, globalVariance);
      }
    }

    template <class dataType>
    void execute(std::vector<ftm::MergeTree<dataType>> &trees,
                 std::vector<ftm::MergeTree<dataType>> &trees2) {
      // --- Preprocessing
      Timer t_preprocess;
      preprocessingTrees<dataType>(trees, treesNodeCorr_);
      if(trees2.size() != 0)
        preprocessingTrees<dataType>(trees2, trees2NodeCorr_);
      printMsg(
        "Preprocessing", 1, t_preprocess.getElapsedTime(), threadNumber_);
      useDoubleInput_ = (trees2.size() != 0);

      // --- Compute principal geodesics
      Timer t_total;
      computePrincipalGeodesics<dataType>(trees, trees2);
      auto totalTime = t_total.getElapsedTime() - t_allVectorCopy_time_;
      printMsg(debug::Separator::L1);
      printMsg("Total time", 1, totalTime, threadNumber_);
      ftm::MergeTree<dataType> barycenter;
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(barycenter_, barycenter);

      // - Compute merge tree geodesic extremities
      computeGeodesicExtremities<dataType>();

      // - Compute branches correlation matrix
      computeBranchesCorrelationMatrix<dataType>(barycenter, trees);

      // --- Reconstruction
      if(doComputeReconstructionError_) {
        auto reconstructionError = computeReconstructionError(
          barycenter, trees, vS_, v2s_, allTreesTs_);
        std::stringstream ss;
        ss << "Reconstruction Error = " << reconstructionError;
        printMsg(ss.str());
      }

      // --- Postprocessing
      if(normalizedWasserstein_) { // keep BDT if input is a PD
        postprocessingPipeline<double>(&(barycenter_.tree));
        if(trees2.size() != 0)
          postprocessingPipeline<double>(&(barycenterInput2_.tree));
      }

      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(barycenter_, barycenter);
      for(unsigned int i = 0; i < trees.size(); ++i) {
        postprocessingPipeline<dataType>(&(trees[i].tree));
        convertBranchDecompositionMatching<dataType>(
          &(barycenter.tree), &(trees[i].tree), baryMatchings_[i]);
      }
      for(unsigned int i = 0; i < trees2.size(); ++i)
        postprocessingPipeline<dataType>(&(trees2[i].tree));
    }

    // ----------------------------------------
    // End functions
    // ----------------------------------------
    void copyMinMaxPairVector(std::vector<std::vector<double>> &v,
                              std::vector<std::vector<double>> &v2,
                              std::vector<std::vector<double>> &trees2V,
                              std::vector<std::vector<double>> &trees2V2) {
      auto root = barycenter_.tree.getRoot();
      auto root2 = barycenterInput2_.tree.getRoot();
      trees2V[root2] = v[root];
      trees2V2[root2] = v2[root];
    }

    template <class dataType>
    void computeGeodesicExtremities() {
      allScaledTs_.resize(allTs_.size(), std::vector<double>(allTs_[0].size()));
      ftm::MergeTree<dataType> barycenter, barycenter2;
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(barycenter_, barycenter);
      if(trees2NodeCorr_.size() != 0)
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
          barycenterInput2_, barycenter2);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < numberOfGeodesics_; ++i) {
        ftm::MergeTree<dataType> extremityV1, extremityV2;
        getInterpolation<dataType>(
          barycenter, vS_[i], v2s_[i], 0.0, extremityV1);
        getInterpolation<dataType>(
          barycenter, vS_[i], v2s_[i], 1.0, extremityV2);
        // Get distance
        dataType distance;
        computeOneDistance(
          extremityV1, extremityV2, distance, true, useDoubleInput_);
        if(trees2NodeCorr_.size() != 0) {
          ftm::MergeTree<dataType> extremity2V1, extremity2V2;
          getInterpolation<dataType>(
            barycenter2, trees2Vs_[i], trees2V2s_[i], 0.0, extremity2V1);
          getInterpolation<dataType>(
            barycenter2, trees2Vs_[i], trees2V2s_[i], 1.0, extremity2V2);
          // Get distance
          dataType distance2;
          computeOneDistance(extremity2V1, extremity2V2, distance2, true,
                             useDoubleInput_, false);
          distance = mixDistances(distance, distance2);
        }
        for(unsigned int j = 0; j < allTs_[i].size(); ++j)
          allScaledTs_[i][j] = allTs_[i][j] * distance;
      }
    }

    template <class dataType>
    void computeBranchesCorrelationMatrix(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees) {
      ttk::MergeTreeAxesAlgorithmBase::computeBranchesCorrelationMatrix(
        barycenter, trees, baryMatchings_, allTs_, branchesCorrelationMatrix_,
        persCorrelationMatrix_);
    }

    // ----------------------------------------
    // Message
    // ----------------------------------------
    template <class dataType>
    void printIterationVariances(ftm::MergeTree<dataType> &barycenter,
                                 std::vector<ftm::MergeTree<dataType>> &trees,
                                 ftm::MergeTree<dataType> &barycenter2,
                                 std::vector<ftm::MergeTree<dataType>> &trees2,
                                 int geodesicNumber,
                                 double globalVariance) {
      bool const printOriginalVariances = false;
      bool const printSurfaceVariance = false;
      bool const printTVariances = true;

      if(printOriginalVariances) {
        // Variance
        double variance = computeExplainedVariance<dataType>(
          barycenter, trees, vS_[geodesicNumber], v2s_[geodesicNumber],
          allTs_[geodesicNumber]);
        double const variancePercent = variance / globalVariance * 100.0;
        std::stringstream ssVariance, ssCumul;
        ssVariance << "Variance explained            : "
                   << round(variancePercent * 100.0) / 100.0 << " %";
        printMsg(ssVariance.str());

        // Cumul Variance
        cumulVariance_ += variance;
        double const cumulVariancePercent
          = cumulVariance_ / globalVariance * 100.0;
        ssCumul << "Cumulative explained variance : "
                << round(cumulVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssCumul.str());
      }

      if(printSurfaceVariance) {
        // Surface Variance
        double surfaceVariance = computeSurfaceExplainedVariance<dataType>(
          barycenter, trees, vS_, v2s_, allTs_);
        double const surfaceVariancePercent
          = surfaceVariance / globalVariance * 100.0;
        std::stringstream ssSurface;
        ssSurface << "Surface Variance explained    : "
                  << round(surfaceVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssSurface.str());
      }

      if(printTVariances) {
        // T-Variance
        double tVariance;
        if(trees2.size() != 0) {
          tVariance = computeExplainedVarianceT(
            barycenter, vS_[geodesicNumber], v2s_[geodesicNumber], barycenter2,
            trees2Vs_[geodesicNumber], trees2V2s_[geodesicNumber],
            allTs_[geodesicNumber]);
        } else
          tVariance = computeExplainedVarianceT(barycenter, vS_[geodesicNumber],
                                                v2s_[geodesicNumber],
                                                allTs_[geodesicNumber]);
        double const tVariancePercent = tVariance / globalVariance * 100.0;
        std::stringstream ssTVariance, ssCumulT;
        ssTVariance << "Explained T-Variance            : "
                    << round(tVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssTVariance.str());

        // Cumul T-Variance
        cumulTVariance_ += tVariance;
        double const cumulTVariancePercent
          = cumulTVariance_ / globalVariance * 100.0;
        ssCumulT << "Cumulative explained T-Variance : "
                 << round(cumulTVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssCumulT.str());
      }
    }

    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    double
      verifyOrthogonality(std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s,
                          bool doPrint = true);
    double
      verifyOrthogonality(std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s,
                          std::vector<std::vector<double>> &v,
                          std::vector<std::vector<double>> &v2,
                          bool doPrint = true);

    template <class dataType>
    dataType computeVarianceFromDistances(std::vector<dataType> &distances);

    template <class dataType>
    double
      computeExplainedVariance(ftm::MergeTree<dataType> &barycenter,
                               std::vector<ftm::MergeTree<dataType>> &trees,
                               std::vector<std::vector<double>> &v,
                               std::vector<std::vector<double>> &v2,
                               std::vector<double> &ts,
                               bool computeGlobalVariance = false);

    template <class dataType>
    double computeGlobalVariance(ftm::MergeTree<dataType> &barycenter,
                                 std::vector<ftm::MergeTree<dataType>> &trees);

    template <class dataType>
    double computeSurfaceExplainedVariance(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<double>> &ts);

    template <class dataType>
    void computeProjectionDistances(ftm::MergeTree<dataType> &barycenter,
                                    std::vector<std::vector<double>> &v,
                                    std::vector<std::vector<double>> &v2,
                                    std::vector<double> &ts,
                                    std::vector<double> &distances,
                                    bool useDoubleInput = false,
                                    bool isFirstInput = true);

    template <class dataType>
    double computeExplainedVarianceT(ftm::MergeTree<dataType> &barycenter,
                                     std::vector<std::vector<double>> &v,
                                     std::vector<std::vector<double>> &v2,
                                     std::vector<double> &ts);

    template <class dataType>
    double computeExplainedVarianceT(ftm::MergeTree<dataType> &barycenter,
                                     std::vector<std::vector<double>> &v,
                                     std::vector<std::vector<double>> &v2,
                                     ftm::MergeTree<dataType> &barycenter2,
                                     std::vector<std::vector<double>> &trees2V,
                                     std::vector<std::vector<double>> &trees2V2,
                                     std::vector<double> &ts);

    // ----------------------------------------
    // Testing
    // ----------------------------------------
    template <class dataType>
    void verifyMinMaxPair(ftm::MergeTree<dataType> &mTree1,
                          ftm::MergeTree<dataType> &mTree) {
      if(not normalizedWasserstein_) // isPersistenceDiagram
        return;
      ftm::FTMTree_MT *tree = &(mTree.tree);
      std::stringstream ss;
      auto birthDeathRoot = tree->getMergedRootBirthDeath<dataType>();
      auto birthRoot = std::get<0>(birthDeathRoot);
      auto deathRoot = std::get<1>(birthDeathRoot);
      bool found = false;
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(tree->isNodeAlone(i))
          continue;
        auto birthDeath = tree->getBirthDeath<dataType>(i);
        auto birth = std::get<0>(birthDeath);
        auto death = std::get<1>(birthDeath);
        if(birth < birthRoot or birth > deathRoot or death < birthRoot
           or death > deathRoot) {
          ss << tree->printNode2<dataType>(i).str() << std::endl;
          found = true;
        }
      }
      if(found) {
        ftm::FTMTree_MT *tree1 = &(mTree1.tree);
        printMsg("Tree1 root:");
        printMsg(tree1->printNode2<dataType>(tree1->getRoot()).str());
        printMsg("Tree root:");
        printMsg(tree->printNode2<dataType>(tree->getRoot()).str());
        printMsg("Tree merged root:");
        printMsg(tree->printMergedRoot<dataType>().str());
        printMsg("Bad pairs:");
        printMsg(ss.str());
        printErr("[computePrincipalGeodesics] tree root is not min max.");
      }
    }
  }; // MergeTreePrincipalGeodesics class
} // namespace ttk

#include <MergeTreePrincipalGeodesicsUtils.h>
