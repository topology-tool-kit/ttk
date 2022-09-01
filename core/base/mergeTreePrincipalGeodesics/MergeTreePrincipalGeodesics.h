/// \ingroup base
/// \class ttk::MergeTreePrincipalGeodesics
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021
///
/// This module defines the %MergeTreePrincipalGeodesics class that computes
/// Principal Geodesic Analysis on the space of merge trees or persistence
/// diagrams, that is, a set of orthognal geodesic axes defining a basis with
/// the barycenter as origin.
///
/// \b Related \b publication: \n
/// "Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Jules Vidal, Julien Tierny.\n

#pragma once

// ttk common includes
#include <Debug.h>
#include <MergeTreeBarycenter.h>
#include <MergeTreePrincipalGeodesicsBase.h>

namespace ttk {

  /**
   * The MergeTreePrincipalGeodesics class provides methods to compute
   * Principal Geodesic Analysis on the space of merge trees or persistence
   * diagrams, that is, a set of orthognal geodesic axes defining a basis with
   * the barycenter as origin.
   */
  class MergeTreePrincipalGeodesics : virtual public Debug,
                                      public MergeTreePrincipalGeodesicsBase {

  protected:
    bool deterministic_ = true;
    unsigned int numberOfGeodesics_ = 1;
    unsigned int noProjectionStep_ = 2;
    // bool rebootDegeneratedGeod_ = true;
    bool doComputeReconstructionError_ = false;
    double barycenterSizeLimitPercent_ = 0.0;
    // TODO keepState works only when enabled before first computation
    bool keepState_ = false;

    // Gradient descent [Seguy and Cuturi, 2015]
    // gradient descent does not work anymore
    bool useGradientDescent_ = false;
    double beta_ = 1e-3; // gradient step
    double lambda_ = 1; // regularizer

    // Proximal Gradient Descent [Cazelles et al., 2018]
    // not implemented
    bool useProximalGradientDescent_ = false;
    double t0_ = 0.0;

    // Advanced parameters
    bool projectInitializedVectors_ = true;
    int parallelMode_ = 0; // 0 : Double _ 1 : Simple _ 3 : Double triangle
                           // inequality _ 2 : Simple triangle inequality

    // Old/Testing
    bool moveVectorsAlongGeodesics_ = true;
    double weightPropCost_ = 1.0;
    double weightOrthoCost_ = 1.0;
    double weightMapCost_ = 1.0;
    double t_vectorCopy_time_ = 0.0, t_allVectorCopy_time_ = 0.0;

    // Filled by the algorithm
    std::vector<double> baryDistances_;
    std::vector<std::vector<double>> allDistances_;
    ftm::MergeTree<double> barycenter_, barycenter2_, barycenterBD_,
      barycenterBD2_;
    bool barycenterWasComputed_ = false;

    int newVectorOffset_ = 0;
    double cumulVariance_ = 0.0, cumulTVariance_ = 0.0;

    // Clean correspondence
    std::vector<std::vector<int>> trees2NodeCorr_;

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
#ifdef TTK_ENABLE_OPENMP
#pragma omp declare reduction( \
  minimum                      \
  : struct Compare             \
  : omp_out = omp_in.bestDistance < omp_out.bestDistance ? omp_in : omp_out)
#pragma omp declare reduction( \
  maximum                      \
  : struct Compare             \
  : omp_out = omp_in.bestDistance > omp_out.bestDistance ? omp_in : omp_out)
#endif

    //----------------------------------------------------------------------------
    // Init
    //----------------------------------------------------------------------------
    template <class dataType>
    void initVectorFromMatching(
      ftm::MergeTree<dataType> &barycenter,
      ftm::MergeTree<dataType> &tree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<std::vector<double>> &v) {
      ftm::FTMTree_MT *barycenterTree = &(barycenter.tree);
      ftm::FTMTree_MT *treeTree = &(tree.tree);

      std::vector<ftm::idNode> matchingVector;
      getMatchingVector<dataType>(barycenter, tree, matching, matchingVector);

      v = std::vector<std::vector<double>>(
        barycenter.tree.getNumberOfNodes(), std::vector<double>(2, 0));
      for(unsigned int j = 0; j < barycenter.tree.getNumberOfNodes(); ++j) {
        if(barycenter.tree.isNodeAlone(j))
          continue;
        auto birthDeathBary
          = getParametrizedBirthDeath<dataType>(barycenterTree, j);
        std::tuple<dataType, dataType> birthDeath;
        if((int)matchingVector[j] != -1) {
          birthDeath
            = getParametrizedBirthDeath<dataType>(treeTree, matchingVector[j]);
        } else {
          dataType projec
            = (std::get<0>(birthDeathBary) + std::get<1>(birthDeathBary)) / 2.0;
          birthDeath = std::make_tuple(projec, projec);
        }
        v[j][0] = std::get<0>(birthDeath) - std::get<0>(birthDeathBary);
        v[j][1] = std::get<1>(birthDeath) - std::get<1>(birthDeathBary);
      }
    }

    template <class dataType>
    void initRandomVector(ftm::MergeTree<dataType> &barycenter,
                          std::vector<std::vector<double>> &v,
                          std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s) {
      // Get average norm of the previous vectors
      std::vector<std::vector<double>> sumVs;
      multiSumVectorFlatten(vS, v2s, sumVs);
      double newNorm = 0;
      for(auto sumVi : sumVs)
        newNorm += norm(sumVi) / sumVs.size();

      // Generate random vector
      v = std::vector<std::vector<double>>(
        barycenter.tree.getNumberOfNodes(), std::vector<double>(2, 0));
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        v[i][0] = (double)rand() / RAND_MAX * newNorm * 2 - newNorm;
        v[i][1] = (double)rand() / RAND_MAX * newNorm * 2 - newNorm;
      }

      // Change the norm of the random vector to be the average norm
      double normV = normFlatten(v);
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        v[i][0] = v[i][0] / normV * newNorm;
        v[i][1] = v[i][1] / normV * newNorm;
      }
    }

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
      bool doOffset = (newVectorOffset_ != 0);
      // Get best distance, best matching and best index
      dataType bestDistance = std::numeric_limits<dataType>::lowest();
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> bestMatching,
        bestMatching2;
      std::vector<std::tuple<double, unsigned int>> distancesAndIndexes(
        trees.size());
      int bestIndex = -1;
      for(unsigned int i = 0; i < trees.size(); ++i) {
        dataType distance = 0.0, distance2 = 0.0;
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching,
          matching2;
        if(geodesicNumber == 0) {
          if(baryDistances_.size() == 0) {
            computeOneDistance<dataType>(
              barycenter, trees[i], matching, distance, false, useDoubleInput_);
            if(trees2.size() != 0) {
              computeOneDistance<dataType>(barycenter2, trees2[i], matching2,
                                           distance2, false, useDoubleInput_,
                                           false);
              distance = mixDistances(distance, distance2);
            }
          } else {
            distance = baryDistances_[i];
            matching = baryMatchings_[i];
            if(trees2.size() != 0)
              matching2 = baryMatchings2_[i];
          }
        } else {
          for(unsigned j = 0; j < allDistances_.size(); ++j)
            distance += allDistances_[j][i];
          distancesAndIndexes[i] = std::make_tuple(distance, i);
        }
        if(distance > bestDistance) {
          bestDistance = distance;
          bestMatching = matching;
          bestMatching2 = matching2;
          bestIndex = i;
        }
      }

      // Sort all distances and their respective indexes
      if(geodesicNumber != 0)
        std::sort(distancesAndIndexes.begin(), distancesAndIndexes.end(),
                  [](const std::tuple<double, unsigned int> &a,
                     const std::tuple<double, unsigned int> &b) -> bool {
                    return (std::get<0>(a) > std::get<0>(b));
                  });

      // Init vectors according farest input
      // (repeat with the ith farest until projection gives non null vector)
      unsigned int i = 0;
      bool foundGoodIndex = false;
      while(not foundGoodIndex) {
        // Get matching of the ith farest input
        if(bestIndex >= 0 and bestIndex < (int)trees.size()) {
          if(geodesicNumber != 0) {
            dataType distance;
            computeOneDistance<dataType>(barycenter, trees[bestIndex],
                                         bestMatching, distance, false,
                                         useDoubleInput_);
            if(trees2.size() != 0)
              computeOneDistance<dataType>(barycenter2, trees2[bestIndex],
                                           bestMatching2, distance, false,
                                           useDoubleInput_, false);
          }

          // Init vectors from matching
          initVectorFromMatching<dataType>(
            barycenter, trees[bestIndex], bestMatching, v1);
          v2 = v1;
          if(trees2.size() != 0) {
            initVectorFromMatching<dataType>(
              barycenter2, trees2[bestIndex], bestMatching2, trees2V1);
            trees2V2 = trees2V1;
          }
        } else {
          initRandomVector(barycenter, v1, vS_, v2s_);
          v2 = v1;
          if(trees2.size() != 0) {
            initRandomVector(barycenter2, trees2V1, trees2Vs_, trees2V2s_);
            trees2V2 = trees2V1;
          }
        }

        // Verify that the two highest persistence pairs are not destroyed
        /*verifyExtremitiesValidity<dataType>(barycenter, trees, v1, v2);
        if(trees2.size() != 0)
          verifyExtremitiesValidity<dataType>(
            barycenter, trees2, trees2V1, trees2V2);*/

        // Project initialized vectors to satisfy constraints
        if(projectInitializedVectors_) {
          projectionStep<dataType>(
            geodesicNumber, barycenter, v1, v2, vS_, v2s_, barycenter2,
            trees2V1, trees2V2, trees2Vs_, trees2V2s_, (trees2.size() != 0), 1);
        }

        // Check if the initialized vectors are good
        foundGoodIndex = (geodesicNumber == 0 or not isVectorNullFlatten(v1));

        // Init next bestIndex
        if(not foundGoodIndex) {
          i += 1;
          if(i < distancesAndIndexes.size())
            bestIndex = std::get<1>(distancesAndIndexes[i]);
          else
            bestIndex = -1;
        }

        // If newVector jump to the next valid bestIndex
        if(foundGoodIndex and doOffset and bestIndex >= 0) {
          bestIndex += newVectorOffset_;
          if(bestIndex >= (int)trees.size())
            bestIndex = -1;
          foundGoodIndex = false;
          doOffset = false;
        }
      }
    }

    //----------------------------------------------------------------------------
    // Barycenter / Interpolation
    //----------------------------------------------------------------------------
    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<double> &finalDistances,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      MergeTreeBarycenter mergeTreeBary;
      mergeTreeBary.setDebugLevel(2);
      // mergeTreeBary.setDebugLevel(4);
      mergeTreeBary.setPreprocess(false);
      mergeTreeBary.setPostprocess(false);
      mergeTreeBary.setBranchDecomposition(true);
      mergeTreeBary.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeBary.setKeepSubtree(false);
      mergeTreeBary.setAssignmentSolver(assignmentSolverID_);
      // mergeTreeBary.setIsCalled(true);
      mergeTreeBary.setThreadNumber(this->threadNumber_);
      mergeTreeBary.setDeterministic(deterministic_);
      mergeTreeBary.setBarycenterSizeLimitPercent(barycenterSizeLimitPercent_);

      matchings = std::vector<
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>(
        trees.size());
      mergeTreeBary.execute<dataType>(
        trees, matchings, baryMergeTree, useDoubleInput, isFirstInput);
      finalDistances = mergeTreeBary.getFinalDistances();
    }

    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings) {
      std::vector<double> finalDistances;
      computeOneBarycenter<dataType>(
        trees, baryMergeTree, matchings, finalDistances);
    }

    template <class dataType>
    void computeOneBarycenter(std::vector<ftm::MergeTree<dataType>> &trees,
                              ftm::MergeTree<dataType> &baryMergeTree) {
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        matchings;
      computeOneBarycenter<dataType>(trees, baryMergeTree, matchings);
    }

    template <class dataType>
    void getParametrizedInterpolation(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<double> &ts,
      double t,
      ftm::MergeTree<dataType> &interpolated) {
      if(not moveVectorsAlongGeodesics_)
        getInterpolation<dataType>(barycenter, v, v2, t, interpolated);
      else
        getMultiInterpolation(barycenter, vS, v2s, v, v2, ts, t, interpolated);
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
      return std::pow(
        (scalarProductFlatten(v, v2) - normFlatten(v) * normFlatten(v2)), 2);
    }

    double projectionCost(std::vector<std::vector<double>> &v,
                          std::vector<std::vector<double>> &v2,
                          std::vector<std::vector<std::vector<double>>> &vS,
                          std::vector<std::vector<std::vector<double>>> &v2s,
                          double optMapCost) {
      return weightPropCost_ * regularizerCost(v, v2)
             + weightOrthoCost_ * orthogonalCost(vS, v2s, v, v2)
             + weightMapCost_ * optMapCost;
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
      double t = (isV1 ? -1.0 : 1.0);

      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      dataType distance;
      std::vector<ftm::idNode> matchingVector;

      if(extremityTree->getRealNumberOfNodes() != 0) {
        computeOneDistance(barycenter, extremity, matching, distance, true,
                           useDoubleInput, isFirstInput);
        getMatchingVector(barycenter, extremity, matching, matchingVector);
      } else
        matchingVector
          = std::vector<ftm::idNode>(barycenterTree->getNumberOfNodes(), -1);

      std::vector<std::vector<double>> oriV = v;
      bool hasChanged = false;
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        auto matched = matchingVector[i];
        auto birthDeathBary
          = getParametrizedBirthDeath<dataType>(barycenterTree, i);
        dataType birthBary = std::get<0>(birthDeathBary);
        dataType deathBary = std::get<1>(birthDeathBary);
        std::vector<double> newV{0.0, 0.0};
        if((int)matched != -1) {
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
        hasChanged |= (i != matched);
      }

      // Compute distance between old and new extremity
      double cost = distanceL2Flatten(v, oriV);
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

    void positivelyProportionnalProjection(std::vector<double> &v1,
                                           std::vector<double> &v2,
                                           std::vector<double> &v1_V1OnV2,
                                           std::vector<double> &v2_V2OnV1,
                                           std::vector<double> &v1_V1AndV2,
                                           std::vector<double> &v2_V1AndV2) {
      // --- Compute all projections
      // Project V1 on V2
      vectorProjection(v2, v1, v1_V1OnV2);

      // Project V2 on V1
      vectorProjection(v1, v2, v2_V2OnV1);

      // Project V1 and V2
      sumProjection(v1, v2, v1_V1AndV2, v2_V1AndV2);
    }

    // Colinearity constraint
    void
      trueGeneralizedGeodesicProjection(std::vector<std::vector<double>> &v1,
                                        std::vector<std::vector<double>> &v2) {
      std::vector<double> v1_flatten, v2_flatten;
      flatten(v1, v1_flatten);
      flatten(v2, v2_flatten);
      double v1_norm = norm(v1_flatten);
      double v2_norm = norm(v2_flatten);
      double beta = v2_norm / (v1_norm + v2_norm);
      std::vector<double> v;
      sumVector(v1_flatten, v2_flatten, v);
      multVectorByScalar(v, (1 - beta), v1_flatten);
      multVectorByScalar(v, beta, v2_flatten);
      unflatten(v1_flatten, v1);
      unflatten(v2_flatten, v2);
    }

    void
      orthogonalProjection(std::vector<std::vector<double>> &v1,
                           std::vector<std::vector<double>> &v2,
                           std::vector<std::vector<std::vector<double>>> &vS,
                           std::vector<std::vector<std::vector<double>>> &v2s) {
      // Multi flatten and sum vS and v2s
      std::vector<std::vector<double>> sumVs;
      multiSumVectorFlatten(vS, v2s, sumVs);

      // Flatten v1 and v2
      std::vector<double> v1_flatten, v2_flatten, v1_proj, v2_proj;
      flatten(v1, v1_flatten);
      flatten(v2, v2_flatten);

      // Call Gram Schmidt
      gramSchmidt(sumVs, v1_flatten, v1_proj);
      gramSchmidt(sumVs, v2_flatten, v2_proj);

      // Unflatten the resulting vectors
      unflatten(v1_proj, v1);
      unflatten(v2_proj, v2);
    }

    bool projectionStepConvergence(
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<double>> &vOld,
      std::vector<std::vector<double>> &v2Old,
      double optMapCost,
      double &oldCost,
      double &oldCost2) {
      bool converged = false;

      // Cost 1
      std::vector<double> v_flat, v2_flat, vOld_flat, v2Old_flat;
      flatten(v, v_flat);
      flatten(v2, v2_flat);
      flatten(vOld, vOld_flat);
      flatten(v2Old, v2Old_flat);
      std::vector<double> subV, subV2;
      subVector(v_flat, vOld_flat, subV);
      subVector(v2_flat, v2Old_flat, subV2);
      double cost
        = norm(subV) / norm(vOld_flat) + norm(subV2) / norm(v2Old_flat);
      // std::cout << "cost = " << cost << std::endl;
      double tol = oldCost / 125.0;
      if(std::abs(cost - oldCost) < tol) {
        // std::cout << "converged cost" << std::endl;
        converged = true;
      }
      oldCost = cost;

      // Cost 2
      double cost2 = projectionCost(v, v2, vS, v2s, optMapCost);
      // std::cout << "cost2 = " << cost2 << std::endl;
      double tol2 = oldCost2 / 125.0;
      if(std::abs(cost2 - oldCost2) < tol2) {
        // std::cout << "converged cost2" << std::endl;
        converged = true;
      }
      oldCost2 = cost2;

      return converged;
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

      // double oldCost = -1, oldCost2 = -1;
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
        if(not useGradientDescent_) {
          printMsg("TGG Proj.", 0, 0, threadNumber_, debug::LineMode::REPLACE,
                   debug::Priority::DETAIL);
          Timer t_trueGeod;
          if(useSecondInput)
            trueGeneralizedGeodesicProjection(vConcat, v2Concat);
          else
            trueGeneralizedGeodesicProjection(v, v2);
          printMsg("TGG Proj.", 1, t_trueGeod.getElapsedTime(), threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
        }

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

        // --- Projection convergence
        /*if(projectionStepConvergence(geodesicNumber, v, v2, vS, v2s, vOld,
                                     v2Old, optMapCost, oldCost, oldCost2))
          break;*/
      }
      return optMapCost;
    }

    //----------------------------------------------------------------------------
    // Assignment
    //----------------------------------------------------------------------------
    template <class dataType>
    void assignmentDoublePara(
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

      std::vector<double> tMin(trees.size(), 0.0), tMax(trees.size(), 1.0);
      dataType geodDistance;
      if(parallelMode_ == 2) {
        ttk::ftm::MergeTree<dataType> extremity1, extremity2;
        getInterpolation<dataType>(barycenter, v, v2, 0.0, extremity1);
        getInterpolation<dataType>(barycenter, v, v2, 1.0, extremity2);
        computeOneDistance<dataType>(
          extremity1, extremity2, geodDistance, false, useDoubleInput_);
        if(trees2.size() != 0) {
          dataType geodDistance2;
          ttk::ftm::MergeTree<dataType> extremity1_2, extremity2_2;
          getInterpolation<dataType>(
            barycenter2, trees2V, trees2V2, 0.0, extremity1_2);
          getInterpolation<dataType>(
            barycenter2, trees2V, trees2V2, 1.0, extremity2_2);
          computeOneDistance<dataType>(extremity1_2, extremity2_2,
                                       geodDistance2, false, useDoubleInput_,
                                       false);
          geodDistance = mixDistances(geodDistance, geodDistance2);
        }
      }

      // Timer t_para;
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
#pragma omp task shared(best, tMin, tMax) firstprivate(i, k)
              {
#endif
                double kT = (k % 2 == 0 ? k / 2 : k_ - 1 - (int)(k / 2));
                double t = 1.0 / (k_ - 1) * kT;

                if(parallelMode_ == 2 and (t < tMin[i] or t > tMax[i])) {
                  std::stringstream ss;
                  ss << i << " _ " << kT << " skipped" << std::endl;
                  printMsg(ss.str());
                } else {
                  dataType distance, distance2;
                  ftm::MergeTree<dataType> interpolated;
                  auto tsToUse = (allTreesTs.size() == 0 ? std::vector<double>()
                                                         : allTreesTs[i]);
                  getParametrizedInterpolation(
                    barycenter, vS, v2s, v, v2, tsToUse, t, interpolated);
                  if(interpolated.tree.getRealNumberOfNodes() != 0) {
                    computeOneDistance<dataType>(
                      interpolated, trees[i], best[i][kT].bestMatching,
                      distance, true, useDoubleInput_);
                    if(trees2.size() != 0) {
                      ftm::MergeTree<dataType> interpolated2;
                      getParametrizedInterpolation(barycenter2, trees2Vs,
                                                   trees2V2s, trees2V, trees2V2,
                                                   tsToUse, t, interpolated2);
                      computeOneDistance<dataType>(
                        interpolated2, trees2[i], best[i][kT].bestMatching2,
                        distance2, true, useDoubleInput_, false);
                      distance = mixDistances(distance, distance2);
                    }
                    best[i][kT].bestDistance = distance;
                    best[i][kT].bestIndex = kT;

                    if(parallelMode_ == 2) {
                      auto bound
                        = (2 * distance + t * geodDistance) / geodDistance;
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
                      {
#endif
                        tMin[i] = std::max(tMin[i], -bound);
                        tMax[i] = std::min(tMax[i], bound);
#ifdef TTK_ENABLE_OPENMP
                      } // pragma omp critical
#endif
                    }
                  }
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

      // std::cout << "t_para : " << t_para.getElapsedTime() << std::endl;
    }

    template <class dataType>
    void assignmentSimplePara(
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

      dataType geodDistance;
      if(parallelMode_ == 3) {
        ttk::ftm::MergeTree<dataType> extremity1, extremity2;
        getInterpolation<dataType>(barycenter, v, v2, 0.0, extremity1);
        getInterpolation<dataType>(barycenter, v, v2, 1.0, extremity2);
        computeOneDistance<dataType>(
          extremity1, extremity2, geodDistance, false, useDoubleInput_);
        if(trees2.size() != 0) {
          dataType geodDistance2;
          ttk::ftm::MergeTree<dataType> extremity1_2, extremity2_2;
          getInterpolation<dataType>(
            barycenter2, trees2V, trees2V2, 0.0, extremity1_2);
          getInterpolation<dataType>(
            barycenter2, trees2V, trees2V2, 1.0, extremity2_2);
          computeOneDistance<dataType>(extremity1_2, extremity2_2,
                                       geodDistance2, false, useDoubleInput_,
                                       false);
          geodDistance = mixDistances(geodDistance, geodDistance2);
        }
      }

      // Timer t_para;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) if(parallelize_)
      {
#pragma omp single nowait
        {
#endif
          for(unsigned int i = 0; i < trees.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i)
            {
#endif
              int noSkipped = 0;
              dataType bestDistance = std::numeric_limits<dataType>::max();
              std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                bestMatching, bestMatching2;
              double bestT = 0.0;

              double tMin = 0.0, tMax = 1.0;

              for(unsigned int k = 0; k < k_; ++k) {
                double kT = (k % 2 == 0 ? k / 2 : k_ - 1 - (int)(k / 2));
                double t = 1.0 / (k_ - 1) * kT;

                if(parallelMode_ == 3 and (t < tMin or t > tMax)) {
                  ++noSkipped;
                  continue;
                }

                std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                  matching, matching2;
                dataType distance, distance2;
                ftm::MergeTree<dataType> interpolated;
                auto tsToUse = (allTreesTs.size() == 0 ? std::vector<double>()
                                                       : allTreesTs[i]);
                getParametrizedInterpolation(
                  barycenter, vS, v2s, v, v2, tsToUse, t, interpolated);
                if(interpolated.tree.getRealNumberOfNodes() != 0) {
                  computeOneDistance<dataType>(interpolated, trees[i], matching,
                                               distance, true, useDoubleInput_);
                  if(trees2.size() != 0) {
                    ftm::MergeTree<dataType> interpolated2;
                    getParametrizedInterpolation(barycenter2, trees2Vs,
                                                 trees2V2s, trees2V, trees2V2,
                                                 tsToUse, t, interpolated2);
                    computeOneDistance<dataType>(interpolated2, trees2[i],
                                                 matching2, distance2, true,
                                                 useDoubleInput_, false);
                    distance = mixDistances(distance, distance2);
                  }
                  if(distance < bestDistance) {
                    bestDistance = distance;
                    bestMatching = matching;
                    bestMatching2 = matching2;
                    bestT = t;
                  }

                  if(parallelMode_ == 3) {
                    auto bound
                      = (2 * distance + t * geodDistance) / geodDistance;
                    tMin = std::max(tMin, -bound);
                    tMax = std::min(tMax, bound);
                  }
                }
              }
              matchings[i] = bestMatching;
              if(trees2.size() != 0)
                matchings2[i] = bestMatching2;
              ts[i] = bestT;
              distances[i] = bestDistance;

              std::stringstream ss;
              ss << i << " _ noSkipped = " << noSkipped << " / " << k_
                 << std::endl;
              printMsg(ss.str());
#ifdef TTK_ENABLE_OPENMP
            } // pragma omp task
#endif
          } // end first for loop
#ifdef TTK_ENABLE_OPENMP
        } // pragma omp single nowait
#pragma omp taskwait
      } // pragma omp parallel
#endif

      // std::cout << "t_para : " << t_para.getElapsedTime() << std::endl;
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

      // Init output
      matchings = std::vector<
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>(
        trees.size());
      matchings2 = std::vector<
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>(
        trees2.size());
      ts = std::vector<double>(trees.size());
      distances = std::vector<double>(trees.size());

      // Assignment
      if(parallelMode_ == 0 or parallelMode_ == 2)
        assignmentDoublePara<dataType>(barycenter, trees, v, v2, barycenter2,
                                       trees2, trees2V, trees2V2, allTreesTs,
                                       vS, v2s, trees2Vs, trees2V2s, matchings,
                                       matchings2, ts, distances);
      else
        assignmentSimplePara<dataType>(barycenter, trees, v, v2, barycenter2,
                                       trees2, trees2V, trees2V2, allTreesTs,
                                       vS, v2s, trees2Vs, trees2V2s, matchings,
                                       matchings2, ts, distances);
    }

    //----------------------------------------------------------------------------
    // Update
    //----------------------------------------------------------------------------
    void diagramAugmentation(std::vector<std::vector<double>> &tree1Matrix,
                             std::vector<std::vector<double>> &tree2Matrix) {
      unsigned int oriTree1Size = tree1Matrix.size(),
                   oriTree2Size = tree2Matrix.size();
      tree1Matrix.resize(oriTree1Size + oriTree2Size);
      for(unsigned int j = 0; j < oriTree2Size; ++j) {
        double projec = (tree2Matrix[j][0] + tree2Matrix[j][1]) / 2.0;
        tree1Matrix[oriTree1Size + j] = std::vector<double>{projec, projec};
      }
      tree2Matrix.resize(oriTree1Size + oriTree2Size);
      for(unsigned int j = 0; j < oriTree1Size; ++j) {
        double projec = (tree1Matrix[j][0] + tree1Matrix[j][1]) / 2.0;
        tree2Matrix[oriTree2Size + j] = std::vector<double>{projec, projec};
      }
    }

    template <class dataType>
    void updateGradientDescent(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<double>> &v,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<double> &ts) {
      std::vector<double> v_flat, v2_flat;
      flatten(v, v_flat);
      flatten(v2, v2_flat);

      // Gradient
      // std::cout << "// Gradient" << std::endl;
      std::vector<std::vector<double>> gradientV1(
        v.size(), std::vector<double>(v[0].size(), 0.0)),
        gradientV2(v2.size(), std::vector<double>(v2[0].size(), 0.0));
      for(unsigned int i = 0; i < trees.size(); ++i) {
        // Matrices
        // std::cout << "// Matrices" << std::endl;
        ftm::MergeTree<dataType> interpolated;
        getInterpolation<dataType>(barycenter, v, v2, ts[i], interpolated);

        std::vector<std::vector<double>> interpolatedMatrix, treeMatrix;
        std::vector<int> interpolatedCorr, treeCorr;
        getTreeMatrix(interpolated, interpolatedMatrix, interpolatedCorr);
        getTreeMatrix(trees[i], treeMatrix, treeCorr);

        std::vector<std::vector<double>> transportationMatrix;
        getTransportationMatrix(interpolated, trees[i], interpolatedCorr,
                                treeCorr, matchings[i], transportationMatrix);

        std::vector<std::vector<double>> weightMatrix(
          transportationMatrix.size(),
          std::vector<double>(transportationMatrix.size(), 0));
        for(unsigned int j = 0; j < transportationMatrix.size(); ++j)
          weightMatrix[j][j] = transportationMatrix.size();

        diagramAugmentation(interpolatedMatrix, treeMatrix);

        // Gradient update
        // std::cout << "// Gradient update" << std::endl;
        std::vector<std::vector<double>> interpolatedMatrixT, treeMatrixT,
          transportationMatrixT, multiplier, temp, temp2;
        transpose(treeMatrix, treeMatrixT);
        transpose(transportationMatrix, transportationMatrixT);
        transpose(interpolatedMatrix, interpolatedMatrixT);
        matrixDot(treeMatrixT, transportationMatrixT, temp);
        matrixDot(temp, weightMatrix, temp2);
        subMatrix(interpolatedMatrixT, temp2, multiplier);

        std::vector<std::vector<double>> updateGradV1, updateGradV2,
          updateGradV1T, updateGradV2T;
        multMatrix(multiplier, 2.0 * (ts[i] - 1), updateGradV1);
        multMatrix(multiplier, 2.0 * ts[i], updateGradV2);
        transpose(updateGradV1, updateGradV1T);
        transpose(updateGradV2, updateGradV2T);

        std::vector<std::vector<double>> realUpdateGradV1(
          v.size(), std::vector<double>(v[0].size(), 0.0)),
          realUpdateGradV2(v2.size(), std::vector<double>(v2[0].size(), 0.0));
        int cpt = 0;
        for(unsigned int j = 0; j < barycenter.tree.getNumberOfNodes(); ++j) {
          if(barycenter.tree.isNodeAlone(j))
            continue;
          /*std::cout << realUpdateGradV1.size() << " _ " << j << std::endl;
          std::cout << realUpdateGradV2.size() << " _ " << j << std::endl;
          std::cout << updateGradV1T.size() << " _ " << cpt << std::endl;
          std::cout << updateGradV2T.size() << " _ " << cpt << std::endl;*/
          realUpdateGradV1[j] = updateGradV1T[cpt];
          realUpdateGradV2[j] = updateGradV2T[cpt];
          ++cpt;
        }

        std::vector<std::vector<double>> tempGradV1 = gradientV1,
                                         tempGradV2 = gradientV2;
        sumMatrix(tempGradV1, realUpdateGradV1, gradientV1);
        sumMatrix(tempGradV2, realUpdateGradV2, gradientV2);
      }

      // Regularizer
      // std::cout << "// Regularizer" << std::endl;
      std::vector<std::vector<double>> regV1(
        v.size(), std::vector<double>(v[0].size(), 0.0)),
        regV2(v2.size(), std::vector<double>(v2[0].size(), 0.0));
      double v1_norm = norm(v_flat);
      double v2_norm = norm(v2_flat);
      double regularizerMult
        = scalarProduct(v_flat, v2_flat) - v1_norm * v2_norm;
      double normDivV1 = v2_norm / v1_norm;
      double normDivV2 = v1_norm / v2_norm;

      std::vector<std::vector<double>> v1_mult, v2_mult, v2_min_v1, v1_min_v2;
      multMatrix(v, normDivV1, v1_mult);
      multMatrix(v2, normDivV2, v2_mult);
      subMatrix(v2, v1_mult, v2_min_v1);
      subMatrix(v, v2_mult, v1_min_v2);
      multMatrix(v2_min_v1, lambda_ * 2 * regularizerMult, regV1);
      multMatrix(v1_min_v2, lambda_ * 2 * regularizerMult, regV2);

      std::vector<std::vector<double>> tempGradV1 = gradientV1,
                                       tempGradV2 = gradientV2;
      sumMatrix(tempGradV1, regV1, gradientV1);
      sumMatrix(tempGradV2, regV2, gradientV2);

      // Vector update
      // std::cout << "// Vector update" << std::endl;
      multMatrix(gradientV1, beta_, gradientV1);
      multMatrix(gradientV2, beta_, gradientV2);
      std::vector<std::vector<double>> tempV = v, tempV2 = v2;
      subMatrix(tempV, gradientV1, v);
      subMatrix(tempV2, gradientV2, v2);
    }

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
      if(moveVectorsAlongGeodesics_ and geodesicNumber != 0) {
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

        if(moveVectorsAlongGeodesics_ and geodesicNumber != 0)
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
          if((int)matchingMatrix[i][j] != -1) {
            auto birthDeath = getParametrizedBirthDeath<dataType>(
              ftmTrees[j], matchingMatrix[i][j]);
            birth = std::get<0>(birthDeath);
            death = std::get<1>(birthDeath);
          }
          allMatched[j] = std::vector<dataType>{birth, death};
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
        double divisorV1
          = one_min_ti_squared - ti_one_min_ti * ti_one_min_ti / ti_squared;
        double divisorV2
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
      allInterpolated = std::vector<ftm::MergeTree<dataType>>(trees.size());
      if(moveVectorsAlongGeodesics_ and geodesicNumber != 0) {
        for(unsigned int i = 0; i < trees.size(); ++i)
          getMultiInterpolation(
            barycenter, vS, v2s, allTreesTs[i], allInterpolated[i]);
      }

      // Manage individuals t
      noUniform = 0;
      foundAllUniform = true;
      isUniform = std::vector<bool>(barycenter.tree.getNumberOfNodes(), false);
      tss
        = std::vector<std::vector<double>>(barycenter.tree.getNumberOfNodes());
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        tss[i] = ts;
        for(unsigned int j = 0; j < tss[i].size(); ++j) {
          auto &treeToUse = (moveVectorsAlongGeodesics_ and geodesicNumber != 0
                               ? allInterpolated[j]
                               : barycenter);
          tss[i][j] = getTNew<dataType>(treeToUse, v, v2, i, ts[j]);
        }
        isUniform[i] = isVectorUniform(tss[i]);
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
      if(useGradientDescent_) {
        updateGradientDescent<dataType>(
          barycenter, trees, v, v2, matchings, ts);
        return false;
      } else
        return updateClosedFormStep<dataType>(
          geodesicNumber, barycenter, trees, v, v2, matchings, vS, v2s,
          barycenter2, trees2, trees2V, trees2V2, matchings2, trees2Vs,
          trees2V2s, ts, allTreesTs);
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
      double lambda = (useGradientDescent_ ? lambda_ : weightPropCost_);
      auto reg = lambda * regularizerCost(v, v2);
      ssReg << "Prop. cost  = " << reg;
      printMsg(ssReg.str());
      if(useGradientDescent_)
        frechetEnergy += reg;

      // Ortho. cost
      std::stringstream ssOrthoCost;
      auto orthoCost = weightOrthoCost_ * orthogonalCost(vS_, v2s_, v, v2);
      ssOrthoCost << "Ortho. cost = " << orthoCost;
      printMsg(ssOrthoCost.str());

      // Map. cost
      std::stringstream ssOptMapCost;
      optMapCost *= weightMapCost_;
      ssOptMapCost << "Map. cost   = " << optMapCost;
      printMsg(ssOptMapCost.str());

      // Total cost
      if(useGradientDescent_) {
        std::stringstream ssTotal;
        ssTotal << "Total       = " << frechetEnergy;
        printMsg(ssTotal.str());
      }

      // Detect convergence
      double tol = 0.01;
      tol = oldFrechetEnergy / 125.0;
      // tol = oldFrechetEnergy / 200.0;
      converged = std::abs(frechetEnergy - oldFrechetEnergy) < tol;
      oldFrechetEnergy = frechetEnergy;

      if(frechetEnergy + 1e-6 < minFrechetEnergy) {
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
        allDistances_[geodesicNumber] = distances;
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
      //--iteration;
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
      allDistances_[geodesicNumber] = bestDistances;
      allTs_[geodesicNumber] = bestTs;

      vS_.push_back(v);
      v2s_.push_back(v2);
      trees2Vs_.push_back(trees2V);
      trees2V2s_.push_back(trees2V2);
      transpose(allTs_, allTreesTs_);
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
        computeOneBarycenter<dataType>(
          trees, barycenter, baryMatchings_, baryDistances_, useDoubleInput_);
        mergeTreeTemplateToDouble(barycenter, barycenterBD_);
        if(trees2.size() != 0) {
          std::vector<double> baryDistances2;
          computeOneBarycenter<dataType>(trees2, barycenter2, baryMatchings2_,
                                         baryDistances2, useDoubleInput_,
                                         false);
          // Copy min max pair if not useMinMaxPair
          // copyMinMaxPair(barycenter, barycenter2);
          mergeTreeTemplateToDouble(barycenter2, barycenterBD2_);
          for(unsigned int i = 0; i < baryDistances_.size(); ++i)
            baryDistances_[i]
              = mixDistances(baryDistances_[i], baryDistances2[i]);

          verifyMinMaxPair(barycenter, barycenter2);
        }
        printMsg("Barycenter", 1, t_barycenter.getElapsedTime(), threadNumber_);
        barycenterWasComputed_ = true;
      } else {
        printMsg("KeepState is enabled and barycenter was already computed");
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
          barycenterBD_, barycenter);
        if(trees2.size() != 0) {
          ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
            barycenterBD2_, barycenter2);
        }
      }
      printMsg(barycenter.tree.printTreeStats().str());
      mergeTreeTemplateToDouble(barycenter, barycenter_);
      if(trees2.size() != 0) {
        printMsg(barycenter2.tree.printTreeStats().str());
        mergeTreeTemplateToDouble(barycenter2, barycenter2_);
      }

      // --- Compute global variance
      double globalVariance = computeVarianceFromDistances(baryDistances_);

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
      unsigned int oldNoGeod = allTs_.size();
      if(not keepState_) {
        allTs_ = std::vector<std::vector<double>>(
          numberOfGeodesics_, std::vector<double>(trees.size()));
        allDistances_ = std::vector<std::vector<double>>(
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
        transpose(allTs_, allTreesTs_);
        allDistances_.resize(
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
      // --- Testing
      // makeMyTest<dataType>(trees, trees2);
      // normalizedWasserstein_ = false;
      // normalizedWasserstein_ = true;
      // this->threadNumber_ = 1;

      // --- Preprocessing
      Timer t_preprocess;
      treesNodeCorr_ = std::vector<std::vector<int>>(trees.size());
      preprocessingTrees<dataType>(trees, treesNodeCorr_);
      trees2NodeCorr_ = std::vector<std::vector<int>>(trees2.size());
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
          postprocessingPipeline<double>(&(barycenter2_.tree));
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
      auto root2 = barycenter2_.tree.getRoot();
      trees2V[root2] = v[root];
      trees2V2[root2] = v2[root];
    }

    template <class dataType>
    void computeGeodesicExtremities() {
      allScaledTs_ = std::vector<std::vector<double>>(
        allTs_.size(), std::vector<double>(allTs_[0].size()));
      ftm::MergeTree<dataType> barycenter, barycenter2;
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(barycenter_, barycenter);
      if(trees2NodeCorr_.size() != 0)
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
          barycenter2_, barycenter2);
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
        // std::cout << distance << std::endl;
        for(unsigned int j = 0; j < allTs_[i].size(); ++j)
          allScaledTs_[i][j] = allTs_[i][j] * distance;
      }
    }

    template <class dataType>
    void computeBranchesCorrelationMatrix(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees) {
      branchesCorrelationMatrix_ = std::vector<std::vector<double>>(
        barycenter.tree.getNumberOfNodes(),
        std::vector<double>(numberOfGeodesics_, 0.0));
      persCorrelationMatrix_ = branchesCorrelationMatrix_;

      // m[i][j] contains the node in trees[j] matched to the node i in the
      // barycenter
      std::vector<std::vector<ftm::idNode>> matchingMatrix;
      getMatchingMatrix(barycenter, trees, baryMatchings_, matchingMatrix);

      std::queue<ftm::idNode> queue;
      queue.emplace(barycenter.tree.getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        // Get births and deaths array
        std::vector<double> births(trees.size(), 0.0),
          deaths(trees.size(), 0.0), pers(trees.size(), 0.0);
        for(unsigned int i = 0; i < trees.size(); ++i) {
          auto matched = matchingMatrix[node][i];
          std::tuple<dataType, dataType> birthDeath;
          if((int)matched == -1) {
            birthDeath = barycenter.tree.template getBirthDeath<dataType>(node);
            auto projec
              = (std::get<0>(birthDeath) + std::get<1>(birthDeath)) / 2.0;
            birthDeath = std::make_tuple(projec, projec);
          } else
            birthDeath
              = trees[i].tree.template getBirthDeath<dataType>(matched);
          births[i] = std::get<0>(birthDeath);
          deaths[i] = std::get<1>(birthDeath);
          pers[i] = deaths[i] - births[i];
        }

        // Compute correlations
        for(unsigned int g = 0; g < numberOfGeodesics_; ++g) {
          double birthCorr = corr(births, allTs_[g]);
          double deathCorr = corr(deaths, allTs_[g]);
          double persCorr = corr(pers, allTs_[g]);

          if(std::isnan(birthCorr))
            birthCorr = 0.0;
          if(std::isnan(deathCorr))
            deathCorr = 0.0;
          if(std::isnan(persCorr))
            persCorr = 0.0;

          auto birthDeathNode
            = barycenter.tree.template getBirthDeathNode<dataType>(node);
          auto birthNode = std::get<0>(birthDeathNode);
          auto deathNode = std::get<1>(birthDeathNode);
          branchesCorrelationMatrix_[birthNode][g] = birthCorr;
          branchesCorrelationMatrix_[deathNode][g] = deathCorr;
          persCorrelationMatrix_[birthNode][g] = persCorr;
          persCorrelationMatrix_[deathNode][g] = persCorr;
        }

        // Push children to the queue
        std::vector<ftm::idNode> children;
        barycenter.tree.getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
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
      bool printOriginalVariances = false;
      bool printSurfaceVariance = false;
      bool printTVariances = true;

      if(printOriginalVariances) {
        // Variance
        double variance = computeExplainedVariance<dataType>(
          barycenter, trees, vS_[geodesicNumber], v2s_[geodesicNumber],
          allTs_[geodesicNumber]);
        double variancePercent = variance / globalVariance * 100.0;
        std::stringstream ssVariance, ssCumul;
        ssVariance << "Variance explained            : "
                   << round(variancePercent * 100.0) / 100.0 << " %";
        printMsg(ssVariance.str());

        // Cumul Variance
        cumulVariance_ += variance;
        double cumulVariancePercent = cumulVariance_ / globalVariance * 100.0;
        ssCumul << "Cumulative explained variance : "
                << round(cumulVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssCumul.str());
      }

      if(printSurfaceVariance) {
        // Surface Variance
        double surfaceVariance = computeSurfaceExplainedVariance<dataType>(
          barycenter, trees, vS_, v2s_, allTs_);
        double surfaceVariancePercent
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
        double tVariancePercent = tVariance / globalVariance * 100.0;
        std::stringstream ssTVariance, ssCumulT;
        ssTVariance << "Explained T-Variance            : "
                    << round(tVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssTVariance.str());

        // Cumul T-Variance
        cumulTVariance_ += tVariance;
        double cumulTVariancePercent = cumulTVariance_ / globalVariance * 100.0;
        ssCumulT << "Cumulative explained T-Variance : "
                 << round(cumulTVariancePercent * 100.0) / 100.0 << " %";
        printMsg(ssCumulT.str());
      }
    }

    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    template <class dataType>
    void getTransportationMatrix(
      ftm::MergeTree<dataType> &tree1,
      ftm::MergeTree<dataType> &tree2,
      std::vector<int> &tree1Corr,
      std::vector<int> &tree2Corr,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<std::vector<double>> &transportationMatrix);

    template <class dataType>
    void getTreeMatrix(ftm::MergeTree<dataType> &tree,
                       std::vector<std::vector<double>> &treeMatrix,
                       std::vector<int> &treeCorr);

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

    template <class dataType>
    std::tuple<dataType, dataType>
      getParametrizedBirthDeath(ftm::FTMTree_MT *tree, ftm::idNode node);

    // ----------------------------------------
    // Testing
    // ----------------------------------------
    template <class dataType>
    void verifyExtremitiesValidity(ftm::MergeTree<dataType> &barycenter,
                                   std::vector<ftm::MergeTree<dataType>> &trees,
                                   std::vector<std::vector<double>> &v1,
                                   std::vector<std::vector<double>> &v2) {
      auto verifyPair = [&](ttk::ftm::FTMTree_MT *tree, ftm::idNode node,
                            double &alphaV1, double &alphaV2) {
        std::tuple<dataType, dataType> bd;
        if(tree->isRoot(node))
          bd = tree->getBirthDeath<dataType>(node);
        else
          bd = getParametrizedBirthDeath<dataType>(tree, node);
        auto newBirthV1 = std::get<0>(bd) - v1[node][0];
        auto newDeathV1 = std::get<1>(bd) - v1[node][1];
        auto newBirthV2 = std::get<0>(bd) + v2[node][0];
        auto newDeathV2 = std::get<1>(bd) + v2[node][1];

        std::cout << newBirthV1 << " _ " << newDeathV1 << std::endl;
        std::cout << newBirthV2 << " _ " << newDeathV2 << std::endl;

        if(newBirthV1 >= newDeathV1) {
          std::cout << "newBirthV1 >= newDeathV1" << std::endl;
          double alpha
            = (std::get<1>(bd) - std::get<0>(bd)) / (v1[node][1] - v1[node][0]);
          alpha /= 2;
          alphaV1 = std::min(alphaV1, alpha);

          std::cout << "    alpha : " << alpha << std::endl;
          std::cout << "    " << std::get<0>(bd) - alpha * v1[node][0] << " _ "
                    << std::get<1>(bd) - alpha * v1[node][1] << std::endl;
        }
        if(newBirthV2 >= newDeathV2) {
          std::cout << "newBirthV2 >= newDeathV2" << std::endl;
          double alpha
            = (std::get<1>(bd) - std::get<0>(bd)) / (v2[node][0] - v2[node][1]);
          alpha /= 2;
          alphaV2 = std::min(alphaV2, alpha);
        }
      };

      double alphaV1 = 1.0, alphaV2 = 1.0;
      for(unsigned int i = 0; i < trees.size(); ++i) {
        ftm::MergeTree<dataType> interpolated;
        ttk::ftm::FTMTree_MT *tree = &(barycenter.tree);
        if(allTreesTs_.size() != 0) {
          getMultiInterpolation(
            barycenter, vS_, v2s_, allTreesTs_[i], interpolated);
          tree = &(interpolated.tree);

          std::cout << "i : " << i << std::endl;
          printMsg(interpolated.tree.printTree().str());
          printMsg(interpolated.tree.template printPairsFromTree<dataType>(true)
                     .str());
        }
        auto nodeSecMax = tree->getSecondMaximumPersistenceNode<dataType>();
        verifyPair(tree, tree->getRoot(), alphaV1, alphaV2);
        verifyPair(tree, nodeSecMax, alphaV1, alphaV2);

        if(allTreesTs_.size() != 0) {
          persistenceThresholding<dataType>(tree, 0.001);
          printMsg(tree->printTree().str());
        }

        if(allTreesTs_.size() == 0)
          break;
      }

      std::cout << "alphaV1 : " << alphaV1 << " _ alphaV2 : " << alphaV2
                << std::endl;
      multVectorByScalarFlatten(v1, alphaV1, v1);
      multVectorByScalarFlatten(v2, alphaV2, v2);

      printVectorOfVector(v1);
      printVectorOfVector(v2);
    }

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

    template <class dataType>
    void isSomePointsBelowDiagonal(ftm::MergeTree<dataType> &barycenter,
                                   std::vector<std::vector<double>> &v,
                                   std::vector<std::vector<double>> &v2,
                                   double t) {
      ftm::FTMTree_MT *barycenterTree = &(barycenter.tree);
      std::vector<dataType> scalarsVector;
      ttk::ftm::getTreeScalars<dataType>(barycenter, scalarsVector);
      std::vector<dataType> interpolationVector(scalarsVector.size());
      int cptBelowDiagonal = 0;

      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        auto iOrigin = barycenterTree->getNode(i)->getOrigin();
        auto iValue = barycenterTree->getValue<dataType>(i);
        auto iOriginValue = barycenterTree->getValue<dataType>(iOrigin);
        auto nodeBirth = (iValue < iOriginValue ? i : iOrigin);
        auto nodeDeath = (iValue < iOriginValue ? iOrigin : i);
        interpolationVector[nodeBirth]
          = scalarsVector[nodeBirth] - v[i][0] + t * (v[i][0] + v2[i][0]);
        interpolationVector[nodeDeath]
          = scalarsVector[nodeDeath] - v[i][1] + t * (v[i][1] + v2[i][1]);
        if(interpolationVector[nodeBirth] > interpolationVector[nodeDeath]) {
          printErr(std::to_string(interpolationVector[nodeBirth]) + " > "
                   + std::to_string(interpolationVector[nodeDeath]));
          ++cptBelowDiagonal;
        }
      }

      if(cptBelowDiagonal != 0)
        printErr(std::to_string(cptBelowDiagonal)
                 + " points in interpolation below diagonal.");
    }

    template <class dataType>
    void makeMyTest(std::vector<ftm::MergeTree<dataType>> &trees,
                    std::vector<ftm::MergeTree<dataType>> &trees2) {
      trees.clear();
      trees2.clear();
      std::vector<std::vector<double>> allScalarsVector;
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>> allArcs;

      int example = 0;

      // Example 1.2
      if(example == 0) {
        int noTrees = 5;
        allScalarsVector = std::vector<std::vector<double>>(noTrees);
        allArcs
          = std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(
            noTrees);
        for(int i = 0; i < noTrees; ++i)
          allArcs[i] = std::vector<std::tuple<ftm::idNode, ftm::idNode>>{
            std::make_tuple(1, 0), std::make_tuple(2, 1), std::make_tuple(3, 1),
            std::make_tuple(4, 3), std::make_tuple(5, 3)};
        allScalarsVector[0] = std::vector<double>{
          1.8800202512849888, 2.9551038837928356, 8.642432224416988,
          2.001355144655213,  9.744227844325652,  10.48293106347959};
        allScalarsVector[1] = std::vector<double>{
          0.910762771193335, 2.6194811162215252, 9.405209753906986,
          2.818094257015912, 9.952873224743906,  10.38959595264684};
        allScalarsVector[2] = std::vector<double>{
          1.689474497866048, 3.114915658665831, 9.273077611933017,
          2.578954627292262, 9.580422748247427, 9.871009661747955};
        allScalarsVector[3] = std::vector<double>{
          1.9153721495469198, 3.9100959288339396, 9.278534823150661,
          4.5328469743884945, 9.253373026738721,  9.648281581657754};
        allScalarsVector[4] = std::vector<double>{
          -0.2205962234257033, 3.914344539876729, 8.917818621169264,
          2.281341122587819,   9.045743067198401, 10.450424079708139};
      }

      // Create trees
      trees = std::vector<ftm::MergeTree<dataType>>(allScalarsVector.size());
      for(unsigned int i = 0; i < allScalarsVector.size(); ++i) {
        std::vector<dataType> scalarsVector(allScalarsVector[i].size());
        for(unsigned int j = 0; j < allScalarsVector[i].size(); ++j)
          scalarsVector[j] = allScalarsVector[i][j];
        trees[i] = makeFakeMergeTree(scalarsVector, allArcs[i]);
      }
    }

  }; // MergeTreePrincipalGeodesics class
} // namespace ttk

#include <MergeTreePrincipalGeodesicsUtils.h>
