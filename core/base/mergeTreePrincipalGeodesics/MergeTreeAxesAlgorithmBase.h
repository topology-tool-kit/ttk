/// \ingroup base
/// \class ttk::MergeTreeAxesAlgorithmBase
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023
///
/// \b Related \b publication: \n
/// "Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Jules Vidal, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2022

#pragma once

#include <Geometry.h>
#include <MergeTreeBarycenter.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>
#include <MergeTreeUtils.h>
#include <Statistics.h>

#define ENERGY_COMPARISON_TOLERANCE 1e-6

namespace ttk {

  class MergeTreeAxesAlgorithmBase : virtual public Debug,
                                     public MergeTreeBase {

  protected:
    bool deterministic_ = true;
    unsigned int numberOfGeodesics_ = 1;
    unsigned int k_ = 10;
    double barycenterSizeLimitPercent_ = 0.0;

    // Clean correspondence
    std::vector<std::vector<int>> trees2NodeCorr_;

  public:
    MergeTreeAxesAlgorithmBase() {
      this->setDebugMsgPrefix(
        "MergeTreeAxesAlgorithmBase"); // inherited from Debug: prefix will
                                       // be printed at the beginning of
                                       // every msg
    }

    //----------------------------------------------------------------------------
    // Matching / Distance
    //----------------------------------------------------------------------------
    template <class dataType>
    void computeOneDistance(
      ftm::MergeTree<dataType> &tree1,
      ftm::MergeTree<dataType> &tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      dataType &distance,
      bool isCalled = false,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      MergeTreeDistance mergeTreeDistance;
      mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
      mergeTreeDistance.setPreprocess(false);
      mergeTreeDistance.setPostprocess(false);
      mergeTreeDistance.setBranchDecomposition(true);
      mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeDistance.setKeepSubtree(false);
      mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
      mergeTreeDistance.setIsCalled(isCalled);
      mergeTreeDistance.setThreadNumber(this->threadNumber_);
      mergeTreeDistance.setDistanceSquaredRoot(true); // squared root
      mergeTreeDistance.setNodePerTask(nodePerTask_);
      if(useDoubleInput) {
        double const weight = mixDistancesMinMaxPairWeight(isFirstInput);
        mergeTreeDistance.setMinMaxPairWeight(weight);
      }
      distance = mergeTreeDistance.computeDistance<dataType>(
        &(tree1.tree), &(tree2.tree), matching);
    }

    template <class dataType>
    void computeOneDistance(ftm::MergeTree<dataType> &tree1,
                            ftm::MergeTree<dataType> &tree2,
                            dataType &distance,
                            bool isCalled = false,
                            bool useDoubleInput = false,
                            bool isFirstInput = true) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
      computeOneDistance<dataType>(tree1, tree2, matching, distance, isCalled,
                                   useDoubleInput, isFirstInput);
    }

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

      v.resize(barycenter.tree.getNumberOfNodes(), std::vector<double>(2, 0));
      for(unsigned int j = 0; j < barycenter.tree.getNumberOfNodes(); ++j) {
        if(barycenter.tree.isNodeAlone(j))
          continue;
        auto birthDeathBary
          = getParametrizedBirthDeath<dataType>(barycenterTree, j);
        std::tuple<dataType, dataType> birthDeath;
        if(matchingVector[j] != std::numeric_limits<ftm::idNode>::max()) {
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
      ttk::Geometry::multiAddVectorsFlatten(vS, v2s, sumVs);
      double newNorm = 0;
      for(auto &sumVi : sumVs)
        newNorm += ttk::Geometry::magnitude(sumVi) / sumVs.size();

      // Generate random vector
      v.resize(barycenter.tree.getNumberOfNodes(), std::vector<double>(2, 0));
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        v[i][0] = (double)rand() / RAND_MAX * newNorm * 2 - newNorm;
        v[i][1] = (double)rand() / RAND_MAX * newNorm * 2 - newNorm;
      }

      // Change the norm of the random vector to be the average norm
      double const normV = ttk::Geometry::magnitudeFlatten(v);
      for(unsigned int i = 0; i < barycenter.tree.getNumberOfNodes(); ++i) {
        if(barycenter.tree.isNodeAlone(i))
          continue;
        v[i][0] = v[i][0] / normV * newNorm;
        v[i][1] = v[i][1] / normV * newNorm;
      }
    }

    template <class dataType, typename F>
    int initVectors(
      int geodesicNumber,
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &barycenter2,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<std::vector<double>> &v1,
      std::vector<std::vector<double>> &v2,
      std::vector<std::vector<double>> &trees2V1,
      std::vector<std::vector<double>> &trees2V2,
      int newVectorOffset,
      std::vector<double> &inputToOriginDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2,
      std::vector<std::vector<double>> &inputToGeodesicsDistances,
      std::vector<std::vector<std::vector<double>>> &vS,
      std::vector<std::vector<std::vector<double>>> &v2s,
      std::vector<std::vector<std::vector<double>>> &trees2Vs,
      std::vector<std::vector<std::vector<double>>> &trees2V2s,
      bool projectInitializedVectors,
      F initializedVectorsProjection) {

      bool doOffset = (newVectorOffset != 0);
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
          if(inputToOriginDistances.size() == 0) {
            computeOneDistance<dataType>(
              barycenter, trees[i], matching, distance, false, useDoubleInput_);
            if(trees2.size() != 0) {
              computeOneDistance<dataType>(barycenter2, trees2[i], matching2,
                                           distance2, false, useDoubleInput_,
                                           false);
              distance = mixDistances(distance, distance2);
            }
          } else {
            distance = inputToOriginDistances[i];
            matching = baryMatchings[i];
            if(trees2.size() != 0)
              matching2 = baryMatchings2[i];
          }
        } else {
          for(unsigned j = 0; j < inputToGeodesicsDistances.size(); ++j)
            distance += inputToGeodesicsDistances[j][i];
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

      // Init vectors according farthest input
      // (repeat with the ith farthest until projection gives non null vector)
      unsigned int i = 0;
      bool foundGoodIndex = false;
      while(not foundGoodIndex) {
        // Get matching of the ith farthest input
        if(bestIndex >= 0 and bestIndex < (int)trees.size()) {
          if(geodesicNumber != 0) {
            if(baryMatchings.size() == 0
               and (baryMatchings.size() == 0 or trees2.size() == 0)) {
              dataType distance;
              computeOneDistance<dataType>(barycenter, trees[bestIndex],
                                           bestMatching, distance, false,
                                           useDoubleInput_);
              if(trees2.size() != 0)
                computeOneDistance<dataType>(barycenter2, trees2[bestIndex],
                                             bestMatching2, distance, false,
                                             useDoubleInput_, false);
            } else {
              bestMatching = baryMatchings[bestIndex];
              if(trees2.size() != 0)
                bestMatching2 = baryMatchings2[bestIndex];
            }
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
          initRandomVector(barycenter, v1, vS, v2s);
          v2 = v1;
          if(trees2.size() != 0) {
            initRandomVector(barycenter2, trees2V1, trees2Vs, trees2V2s);
            trees2V2 = trees2V1;
          }
        }

        // Project initialized vectors to satisfy constraints
        if(projectInitializedVectors) {
          initializedVectorsProjection(
            geodesicNumber, barycenter, v1, v2, vS, v2s, barycenter2, trees2V1,
            trees2V2, trees2Vs, trees2V2s, (trees2.size() != 0), 1);
        }

        // Check if the initialized vectors are good
        foundGoodIndex
          = (geodesicNumber == 0 or not ttk::Geometry::isVectorNullFlatten(v1));

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
          bestIndex += newVectorOffset;
          if(bestIndex >= (int)trees.size())
            bestIndex = -1;
          foundGoodIndex = false;
          doOffset = false;
        }
      }
      return bestIndex;
    }

    //----------------------------------------------------------------------------
    // Barycenter
    //----------------------------------------------------------------------------
    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<double> &finalDistances,
      double barycenterSizeLimitPercent,
      unsigned int barycenterMaximumNumberOfPairs,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      MergeTreeBarycenter mergeTreeBary;
      mergeTreeBary.setDebugLevel(std::min(debugLevel_, 2));
      mergeTreeBary.setPreprocess(false);
      mergeTreeBary.setPostprocess(false);
      mergeTreeBary.setBranchDecomposition(true);
      mergeTreeBary.setNormalizedWasserstein(normalizedWasserstein_);
      mergeTreeBary.setKeepSubtree(false);
      mergeTreeBary.setAssignmentSolver(assignmentSolverID_);
      mergeTreeBary.setThreadNumber(this->threadNumber_);
      mergeTreeBary.setDeterministic(deterministic_);
      mergeTreeBary.setBarycenterSizeLimitPercent(barycenterSizeLimitPercent);
      mergeTreeBary.setBarycenterMaximumNumberOfPairs(
        barycenterMaximumNumberOfPairs);

      matchings.resize(trees.size());
      mergeTreeBary.execute<dataType>(
        trees, matchings, baryMergeTree, useDoubleInput, isFirstInput);
      finalDistances = mergeTreeBary.getFinalDistances();
    }

    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<double> &finalDistances,
      double barycenterSizeLimitPercent,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      unsigned int const barycenterMaximumNumberOfPairs = 0;
      computeOneBarycenter(trees, baryMergeTree, matchings, finalDistances,
                           barycenterSizeLimitPercent,
                           barycenterMaximumNumberOfPairs, useDoubleInput,
                           isFirstInput);
    }

    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::MergeTree<dataType>> &trees,
      ftm::MergeTree<dataType> &baryMergeTree,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<double> &finalDistances,
      bool useDoubleInput = false,
      bool isFirstInput = true) {
      computeOneBarycenter(trees, baryMergeTree, matchings, finalDistances,
                           barycenterSizeLimitPercent_, useDoubleInput,
                           isFirstInput);
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

    //----------------------------------------------------------------------------
    // Preprocessing
    //----------------------------------------------------------------------------
    template <class dataType>
    void preprocessingTrees(std::vector<ftm::MergeTree<dataType>> &trees,
                            std::vector<std::vector<int>> &nodeCorr,
                            bool useMinMaxPairT = true) {
      nodeCorr.resize(trees.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i) {
        preprocessingPipeline<dataType>(
          trees[i], epsilonTree1_, epsilon2Tree1_, epsilon3Tree1_,
          branchDecomposition_, useMinMaxPairT, cleanTree_, nodeCorr[i]);
        if(trees.size() < 40)
          printTreeStats(trees[i]);
      }
      if(trees.size() != 0)
        printTreesStats(trees);
    }

    template <class dataType>
    void preprocessingTrees(std::vector<ftm::MergeTree<dataType>> &trees,
                            bool useMinMaxPairT = true) {
      std::vector<std::vector<int>> nodeCorr(trees.size());
      preprocessingTrees(trees, nodeCorr, useMinMaxPairT);
    }

    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    // v[i] contains the node in tree matched to the node i in barycenter
    template <class dataType>
    void getMatchingVector(
      ftm::MergeTree<dataType> &barycenter,
      ftm::MergeTree<dataType> &tree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<ftm::idNode> &matchingVector) {
      matchingVector.clear();
      matchingVector.resize(barycenter.tree.getNumberOfNodes(),
                            std::numeric_limits<ftm::idNode>::max());
      for(unsigned int j = 0; j < matchings.size(); ++j) {
        auto match0 = std::get<0>(matchings[j]);
        auto match1 = std::get<1>(matchings[j]);
        if(match0 < barycenter.tree.getNumberOfNodes()
           and match1 < tree.tree.getNumberOfNodes())
          matchingVector[match0] = match1;
      }
    }

    // v[i] contains the node in barycenter matched to the node i in tree
    template <class dataType>
    void getInverseMatchingVector(
      ftm::MergeTree<dataType> &barycenter,
      ftm::MergeTree<dataType> &tree,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<ftm::idNode> &matchingVector) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> invMatchings(
        matchings.size());
      for(unsigned int i = 0; i < matchings.size(); ++i)
        invMatchings[i] = std::make_tuple(std::get<1>(matchings[i]),
                                          std::get<0>(matchings[i]),
                                          std::get<2>(matchings[i]));
      getMatchingVector(tree, barycenter, invMatchings, matchingVector);
    }

    // m[i][j] contains the node in trees[j] matched to the node i in the
    // barycenter
    template <class dataType>
    void getMatchingMatrix(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<ftm::idNode>> &matchingMatrix) {
      matchingMatrix.clear();
      matchingMatrix.resize(
        barycenter.tree.getNumberOfNodes(),
        std::vector<ftm::idNode>(
          trees.size(), std::numeric_limits<ftm::idNode>::max()));
      for(unsigned int i = 0; i < trees.size(); ++i) {
        std::vector<ftm::idNode> matchingVector;
        getMatchingVector<dataType>(
          barycenter, trees[i], matchings[i], matchingVector);
        for(unsigned int j = 0; j < matchingVector.size(); ++j)
          matchingMatrix[j][i] = matchingVector[j];
      }
    }

    template <class dataType>
    std::tuple<dataType, dataType>
      getParametrizedBirthDeath(ftm::FTMTree_MT *tree, ftm::idNode node) {
      return ttk::getParametrizedBirthDeath<dataType>(
        tree, node, normalizedWasserstein_);
    }

    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    template <class dataType>
    void computeBranchesCorrelationMatrix(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<double>> &allTs,
      std::vector<std::vector<double>> &branchesCorrelationMatrix,
      std::vector<std::vector<double>> &persCorrelationMatrix) {
      branchesCorrelationMatrix.resize(barycenter.tree.getNumberOfNodes(),
                                       std::vector<double>(allTs.size(), 0.0));
      persCorrelationMatrix = branchesCorrelationMatrix;

      // m[i][j] contains the node in trees[j] matched to the node i in the
      // barycenter
      std::vector<std::vector<ftm::idNode>> matchingMatrix;
      getMatchingMatrix(barycenter, trees, baryMatchings, matchingMatrix);

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
          if(matched == std::numeric_limits<ftm::idNode>::max()) {
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
        for(unsigned int g = 0; g < allTs.size(); ++g) {
          double birthCorr = ttk::Statistics::corr(births, allTs[g]);
          double deathCorr = ttk::Statistics::corr(deaths, allTs[g]);
          double persCorr = ttk::Statistics::corr(pers, allTs[g]);

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
          branchesCorrelationMatrix[birthNode][g] = birthCorr;
          branchesCorrelationMatrix[deathNode][g] = deathCorr;
          persCorrelationMatrix[birthNode][g] = persCorr;
          persCorrelationMatrix[deathNode][g] = persCorr;
        }

        // Push children to the queue
        std::vector<ftm::idNode> children;
        barycenter.tree.getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }
    }

    //----------------------------------------------------------------------------
    // Output Utils
    //----------------------------------------------------------------------------
    void zeroPadding(std::string &colName,
                     const size_t numberCols,
                     const size_t colIdx) {
      std::string const max{std::to_string(numberCols - 1)};
      std::string const cur{std::to_string(colIdx)};
      std::string const zer(max.size() - cur.size(), '0');
      colName.append(zer).append(cur);
    }

    std::string getTableCoefficientName(int noGeodesics, int geodesicNum) {
      std::string name{"T"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableCoefficientNormName(int noGeodesics, int geodesicNum) {
      std::string name{"TNorm"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableVectorName(int noGeodesics,
                                   int geodesicNum,
                                   int vId,
                                   int vComp,
                                   bool isSecondInput = false) {
      std::string indexString{};
      zeroPadding(indexString, noGeodesics, geodesicNum);
      std::string const prefix{(isSecondInput ? "T2_" : "")};
      std::string name{prefix + "V" + indexString + "_" + std::to_string(vId)
                       + "_" + std::to_string(vComp)};
      return name;
    }

    std::string getTableCorrelationName(int noGeodesics, int geodesicNum) {
      std::string name{"Corr"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableCorrelationPersName(int noGeodesics, int geodesicNum) {
      std::string name{"CorrPers"};
      zeroPadding(name, noGeodesics, geodesicNum);
      return name;
    }

    std::string getTableCorrelationTreeName(int noTrees, int treeNum) {
      std::string name{"Tree"};
      zeroPadding(name, noTrees, treeNum);
      return name;
    }
  }; // MergeTreeAxesAlgorithmBase class

} // namespace ttk
