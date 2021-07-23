/// \ingroup base
/// \class MergeTreeClustering
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// This module defines the %MergeTreeClustering class that computes
/// a clustering of an ensemble of merge trees in k clusters
/// It is actually a version of KMeans (dynamic clustering algorithm)
/// where centroids are merge trees
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021

#ifndef _MERGETREECLUSTERING_H
#define _MERGETREECLUSTERING_H

#define treesMatchingVector \
  std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
#define matchingVector std::vector<treesMatchingVector>

#pragma once

#include <random>

// ttk common includes
#include <Debug.h>

#include "MergeTreeBarycenter.h"

namespace ttk {

  /**
   * The %MergeTreeClustering class that computes
   * a clustering of an ensemble of merge trees in k clusters
   */
  // TODO rename dataType2 to dataType and remove template everywhere else in
  // this class
  template <class dataType2>
  class MergeTreeClustering : virtual public Debug, public MergeTreeBarycenter {

  private:
    bool parallelizeUpdate_ = true;

    unsigned int noCentroids_ = 2;
    double mixtureCoefficient_ = 0.5;

    // Progressive parameters
    int noIterationC_ = 0;
    double addDeletedNodesTime_ = 0;

    // Accelerated KMeans
    bool acceleratedInitialized_ = false;
    std::vector<std::vector<double>> lowerBound_;
    std::vector<double> upperBound_;
    std::vector<int> bestCentroid_, oldBestCentroid_;
    std::vector<double> bestDistance_;
    std::vector<bool> recompute_;
    std::vector<MergeTree<dataType2>> oldCentroids_, oldCentroids2_;

    // Clean correspondence
    std::vector<std::vector<int>> trees2NodeCorr_;

  public:
    MergeTreeClustering() {
      this->setDebugMsgPrefix(
        "MergeTreeClustering"); // inherited from Debug: prefix will be printed
                                // at the beginning of every msg
    };
    ~MergeTreeClustering(){};

    void setNoCentroids(unsigned int noCentroidsT) {
      noCentroids_ = noCentroidsT;
    }

    void setMixtureCoefficient(double coef) {
      mixtureCoefficient_ = coef;
    }

    std::vector<std::vector<int>> getTrees2NodeCorr() {
      return trees2NodeCorr_;
    }

    /**
     * Implementation of the algorithm.
     */

    template <class dataType>
    double mixDistances(dataType distanceJoin, dataType distanceSplit) {
      return mixtureCoefficient_ * distanceJoin
             + (1 - mixtureCoefficient_) * distanceSplit;
    }

    // KMeans++ init
    template <class dataType>
    void initCentroids(
      std::vector<ftm::FTMTree_MT *> &trees,
      std::vector<ftm::FTMTree_MT *> &trees2,
      std::vector<std::vector<MergeTree<dataType>>> &allCentroids) {
      allCentroids = std::vector<std::vector<MergeTree<dataType>>>(
        2, std::vector<MergeTree<dataType>>(noCentroids_));
      std::vector<dataType> distances(
        trees.size(), std::numeric_limits<dataType>::max());
      for(unsigned int i = 0; i < noCentroids_; ++i) {
        int bestIndex = -1;
        if(i == 0) {
          bestIndex = getBestInitTreeIndex<dataType>(trees, trees2, false);
          allCentroids[0][i] = copyMergeTree<dataType>(trees[bestIndex], true);
          cleanMergeTree<dataType>(allCentroids[0][i]);
          if(trees2.size() != 0) {
            allCentroids[1][i]
              = copyMergeTree<dataType>(trees2[bestIndex], true);
            cleanMergeTree<dataType>(allCentroids[1][i]);
          }
        } else {
          // Create vector of probabilities
          dataType sum = 0;
          for(auto val : distances)
            sum += val;
          dataType bestValue = std::numeric_limits<dataType>::lowest();
          std::vector<dataType> probabilities(trees.size());
          for(unsigned int j = 0; j < distances.size(); ++j) {
            probabilities[j] = distances[j] / sum;
            if(probabilities[j] > bestValue) {
              bestValue = probabilities[j];
              bestIndex = j;
            }
          }
          if(not deterministic_) {
            std::random_device rd;
            std::default_random_engine generator(rd());
            std::discrete_distribution<int> distribution(
              probabilities.begin(), probabilities.end());
            bestIndex = distribution(generator);
          }
          // Create new centroid
          allCentroids[0][i] = copyMergeTree<dataType>(trees[bestIndex], true);
          cleanMergeTree<dataType>(allCentroids[0][i]);
          if(trees2.size() != 0) {
            allCentroids[1][i]
              = copyMergeTree<dataType>(trees2[bestIndex], true);
            cleanMergeTree<dataType>(allCentroids[1][i]);
          }
        }

        if(i == noCentroids_ - 1)
          continue;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) shared(allCentroids) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
        for(unsigned int j = 0; j < trees.size(); ++j) {
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching,
            matching2;
          dataType distanceT, distanceT2;
          computeOneDistance<dataType>(
            trees[j], allCentroids[0][i], matching, distanceT);
          if(trees2.size() != 0) {
            computeOneDistance<dataType>(
              trees2[j], allCentroids[1][i], matching2, distanceT2);
            distanceT = mixDistances<dataType>(distanceT, distanceT2);
          }
          distances[j] = std::min(distances[j], distanceT);
        }
      }
    }

    template <class dataType>
    void initNewCentroid(std::vector<ftm::FTMTree_MT *> &trees,
                         MergeTree<dataType> &centroid) {
      int bestIndex = -1;
      dataType bestValue = std::numeric_limits<dataType>::lowest();
      for(unsigned int i = 0; i < bestDistance_.size(); ++i) {
        if(bestValue < bestDistance_[i]) {
          bestValue = bestDistance_[i];
          bestIndex = bestCentroid_[i];
        }
      }
      centroid = copyMergeTree<dataType>(trees[bestIndex], true);
      cleanMergeTree<dataType>(centroid);
    }

    template <class dataType>
    void
      assignmentCentroidsNaive(std::vector<ftm::FTMTree_MT *> &trees,
                               std::vector<MergeTree<dataType>> &centroids,
                               std::vector<std::tuple<int, int>> &assignmentC,
                               std::vector<dataType> &bestDistanceT,
                               std::vector<ftm::FTMTree_MT *> &trees2,
                               std::vector<MergeTree<dataType>> &centroids2) {
      std::vector<int> bestCentroidT(trees.size(), -1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) shared(centroids, centroids2) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i) {
        for(unsigned int j = 0; j < centroids.size(); ++j) {
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching;
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching2;
          dataType distance, distance2;
          computeOneDistance<dataType>(
            trees[i], centroids[j], matching, distance);
          if(trees2.size() != 0) {
            computeOneDistance<dataType>(
              trees2[i], centroids2[j], matching2, distance2);
            distance = mixDistances<dataType>(distance, distance2);
          }
          if(distance < bestDistanceT[i]) {
            bestDistanceT[i] = distance;
            bestDistance_[i] = distance;
            bestCentroidT[i] = j;
            bestCentroid_[i] = j;
          }
        }
      }

      for(unsigned int i = 0; i < bestCentroidT.size(); ++i)
        assignmentC.push_back(std::make_tuple(bestCentroidT[i], i));
    }

    template <class dataType>
    void getCentroidsDistanceMatrix(
      std::vector<MergeTree<dataType>> &centroids,
      std::vector<std::vector<double>> &distanceMatrix) {
      std::vector<ftm::FTMTree_MT *> trees;
      for(MergeTree<dataType> &centroid : centroids)
        trees.push_back(&(centroid.tree));
      getDistanceMatrix<dataType>(trees, distanceMatrix);
    }

    void hadamardProduct(std::vector<std::vector<double>> &distanceMatrix,
                         std::vector<std::vector<double>> &distanceMatrix2) {
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i)
        for(unsigned int j = 0; j < distanceMatrix[i].size(); ++j)
          distanceMatrix[i][j]
            = mixDistances<double>(distanceMatrix[i][j], distanceMatrix2[i][j]);
    }

    template <class dataType>
    void
      initAcceleratedKMeansVectors(std::vector<ftm::FTMTree_MT *> &trees,
                                   std::vector<MergeTree<dataType>> &centroids,
                                   std::vector<ftm::FTMTree_MT *> &trees2) {
      lowerBound_ = std::vector<std::vector<double>>(
        trees.size(), std::vector<double>(centroids.size(), 0));
      upperBound_
        = std::vector<double>(trees.size(), std::numeric_limits<double>::max());
      bestCentroid_ = std::vector<int>(trees.size(), -1);
      oldBestCentroid_ = std::vector<int>(trees.size(), -1);
      bestDistance_
        = std::vector<double>(trees.size(), std::numeric_limits<double>::max());
      recompute_ = std::vector<bool>(trees.size(), true);
    }

    template <class dataType>
    void initAcceleratedKMeans(std::vector<ftm::FTMTree_MT *> &trees,
                               std::vector<MergeTree<dataType>> &centroids,
                               std::vector<ftm::FTMTree_MT *> &trees2,
                               std::vector<MergeTree<dataType>> &centroids2) {
      acceleratedInitialized_ = true;
      std::vector<std::tuple<int, int>> assignmentC;
      std::vector<dataType> bestDistanceT(
        trees.size(), std::numeric_limits<dataType>::max());
      assignmentCentroidsNaive<dataType>(
        trees, centroids, assignmentC, bestDistanceT, trees2, centroids2);
      for(unsigned int i = 0; i < bestDistanceT.size(); ++i)
        bestDistance_[i] = bestDistanceT[i];
      for(auto asgn : assignmentC)
        bestCentroid_[std::get<1>(asgn)] = std::get<0>(asgn);
      for(unsigned int i = 0; i < bestDistance_.size(); ++i)
        upperBound_[i] = bestDistance_[i];
    }

    template <class dataType>
    void copyCentroids(std::vector<MergeTree<dataType>> &centroids,
                       std::vector<MergeTree<dataType>> &oldCentroids) {
      oldCentroids.clear();
      for(unsigned int i = 0; i < centroids.size(); ++i)
        oldCentroids.push_back(copyMergeTree<dataType>(centroids[i]));
    }

    template <class dataType>
    void assignmentCentroidsAccelerated(
      std::vector<ftm::FTMTree_MT *> &trees,
      std::vector<MergeTree<dataType>> &centroids,
      std::vector<std::tuple<int, int>> &assignmentC,
      std::vector<dataType> &bestDistanceT,
      std::vector<ftm::FTMTree_MT *> &trees2,
      std::vector<MergeTree<dataType>> &centroids2) {
      if(not acceleratedInitialized_) {
        initAcceleratedKMeans<dataType>(trees, centroids, trees2, centroids2);
      } else {
        // Compute distance between old and new corresponding centroids
        std::vector<dataType> distanceShift(centroids.size()),
          distanceShift2(centroids2.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)                     \
  shared(centroids, centroids2, oldCentroids_, oldCentroids2_) \
    num_threads(this->threadNumber_) if(parallelize_)
#endif
        for(unsigned int i = 0; i < centroids.size(); ++i) {
          std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> matching,
            matching2;
          computeOneDistance<dataType>(
            centroids[i], oldCentroids_[i], matching, distanceShift[i]);
          if(trees2.size() != 0) {
            computeOneDistance<dataType>(
              centroids2[i], oldCentroids2_[i], matching2, distanceShift2[i]);
            distanceShift[i]
              = mixDistances<dataType>(distanceShift[i], distanceShift2[i]);
          }
        }

        // Step 5
        for(unsigned int i = 0; i < trees.size(); ++i)
          for(unsigned int c = 0; c < centroids.size(); ++c)
            lowerBound_[i][c]
              = std::max(lowerBound_[i][c] - distanceShift[c], 0.0);

        // Step 6
        for(unsigned int i = 0; i < trees.size(); ++i) {
          upperBound_[i] = upperBound_[i] + distanceShift[bestCentroid_[i]];
          recompute_[i] = true;
        }
      }

      // Step 1
      std::vector<std::vector<double>> centroidsDistance, centroidsDistance2;
      getCentroidsDistanceMatrix<dataType>(centroids, centroidsDistance);
      if(trees2.size() != 0) {
        getCentroidsDistanceMatrix<dataType>(centroids2, centroidsDistance2);
        hadamardProduct(centroidsDistance, centroidsDistance2);
      }
      std::vector<double> centroidScore(
        centroids.size(), std::numeric_limits<double>::max());
      for(unsigned int i = 0; i < centroids.size(); ++i)
        for(unsigned int j = i + 1; j < centroids.size(); ++j) {
          if(0.5 * centroidsDistance[i][j] < centroidScore[i])
            centroidScore[i] = 0.5 * centroidsDistance[i][j];
          if(0.5 * centroidsDistance[i][j] < centroidScore[j])
            centroidScore[j] = 0.5 * centroidsDistance[i][j];
        }

      // Step 2
      std::vector<bool> identified(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i)
        identified[i] = (upperBound_[i] <= centroidScore[bestCentroid_[i]]);

        // Step 3
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) shared(centroids, centroids2) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i)
        for(unsigned int c = 0; c < centroids.size(); ++c) {
          if(not identified[i] and (int) c != bestCentroid_[i]
             and upperBound_[i] > lowerBound_[i][c]
             and upperBound_[i]
                   > 0.5 * centroidsDistance[bestCentroid_[i]][c]) {
            // Step 3a
            if(recompute_[i]) {
              std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                matching, matching2;
              dataType distance, distance2;
              computeOneDistance<dataType>(
                trees[i], centroids[bestCentroid_[i]], matching, distance);
              if(trees2.size() != 0) {
                computeOneDistance<dataType>(trees2[i],
                                             centroids2[bestCentroid_[i]],
                                             matching2, distance2);
                distance = mixDistances<dataType>(distance, distance2);
              }
              recompute_[i] = false;
              lowerBound_[i][bestCentroid_[i]] = distance;
              upperBound_[i] = distance;
              bestDistance_[i] = distance;
            } else {
              bestDistance_[i] = upperBound_[i];
            }
            // Step 3b
            if(bestDistance_[i] > lowerBound_[i][c]
               and bestDistance_[i]
                     > 0.5 * centroidsDistance[bestCentroid_[i]][c]) {
              std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                matching, matching2;
              dataType distance, distance2;
              computeOneDistance<dataType>(
                trees[i], centroids[c], matching, distance);
              if(trees2.size() != 0) {
                computeOneDistance<dataType>(
                  trees2[i], centroids2[c], matching2, distance2);
                distance = mixDistances<dataType>(distance, distance2);
              }
              lowerBound_[i][c] = distance;
              if(distance < bestDistance_[i]) {
                bestCentroid_[i] = c;
                upperBound_[i] = distance;
                bestDistance_[i] = distance;
              }
            }
          }
        }

      // Copy centroids for next step
      copyCentroids<dataType>(centroids, oldCentroids_);
      if(trees2.size() != 0)
        copyCentroids<dataType>(centroids2, oldCentroids2_);

      // Manage output
      for(unsigned int i = 0; i < bestDistance_.size(); ++i)
        bestDistanceT[i] = bestDistance_[i];
      for(unsigned int i = 0; i < bestCentroid_.size(); ++i)
        assignmentC.push_back(std::make_tuple(bestCentroid_[i], i));
    }

    template <class dataType>
    void assignmentCentroids(std::vector<ftm::FTMTree_MT *> &trees,
                             std::vector<MergeTree<dataType>> &centroids,
                             std::vector<std::tuple<int, int>> &assignmentC,
                             std::vector<dataType> &bestDistanceT,
                             std::vector<ftm::FTMTree_MT *> &trees2,
                             std::vector<MergeTree<dataType>> &centroids2) {
      oldBestCentroid_ = bestCentroid_;
      if(normalizedWasserstein_ and rescaledWasserstein_)
        assignmentCentroidsNaive<dataType>(
          trees, centroids, assignmentC, bestDistanceT, trees2, centroids2);
      else
        assignmentCentroidsAccelerated<dataType>(
          trees, centroids, assignmentC, bestDistanceT, trees2, centroids2);
    }

    template <class dataType>
    void
      finalAssignmentCentroids(std::vector<ftm::FTMTree_MT *> &trees,
                               std::vector<MergeTree<dataType>> &centroids,
                               matchingVector &matchingsC,
                               std::vector<std::tuple<int, int>> &assignmentC,
                               std::vector<dataType> &bestDistanceT,
                               std::vector<ftm::FTMTree_MT *> &trees2,
                               std::vector<MergeTree<dataType>> &centroids2,
                               matchingVector &matchingsC2) {
      int noC = centroids.size();
      std::vector<std::vector<ftm::FTMTree_MT *>> assignedTrees(noC),
        assignedTrees2(noC);
      std::vector<std::vector<int>> assignedTreesIndex(noC);

      for(auto asgn : assignmentC) {
        assignedTreesIndex[std::get<0>(asgn)].push_back(std::get<1>(asgn));
        assignedTrees[std::get<0>(asgn)].push_back(trees[std::get<1>(asgn)]);
        if(trees2.size() != 0)
          assignedTrees2[std::get<0>(asgn)].push_back(
            trees2[std::get<1>(asgn)]);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) shared(centroids, centroids2) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
      for(unsigned int i = 0; i < centroids.size(); ++i) {
        std::vector<dataType> distances(assignedTrees[i].size(), 0);
        std::vector<dataType> distances2(assignedTrees[i].size(), 0);
        treesMatchingVector matching(trees.size()), matching2(trees2.size());
        assignment<dataType>(
          assignedTrees[i], centroids[i], matching, distances);
        matchingsC[i] = matching;
        if(trees2.size() != 0) {
          assignment<dataType>(
            assignedTrees2[i], centroids2[i], matching2, distances2);
          matchingsC2[i] = matching2;
          for(unsigned int j = 0; j < assignedTreesIndex[i].size(); ++j)
            distances[j] = mixDistances<dataType>(distances[j], distances2[j]);
        }
        for(unsigned int j = 0; j < assignedTreesIndex[i].size(); ++j) {
          int index = assignedTreesIndex[i][j];
          bestDistanceT[index] = distances[j];
        }
      }
    }

    void matchingCorrespondence(treesMatchingVector &matchingT,
                                std::vector<int> &nodeCorr,
                                std::vector<int> assignedTreesIndex) {
      for(int i : assignedTreesIndex) {
        std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> newMatching;
        for(auto tup : matchingT[i])
          newMatching.push_back(std::make_tuple(
            nodeCorr[std::get<0>(tup)], std::get<1>(tup), std::get<2>(tup)));
        matchingT[i] = newMatching;
      }
    }

    bool samePreviousAssignment(int clusterId) {
      for(unsigned int i = 0; i < bestCentroid_.size(); ++i)
        if(bestCentroid_[i] == clusterId
           and bestCentroid_[i] != oldBestCentroid_[i])
          return false;
      return true;
    }

    template <class dataType>
    bool updateCentroids(std::vector<ftm::FTMTree_MT *> &trees,
                         std::vector<MergeTree<dataType>> &centroids,
                         std::vector<double> &alphas,
                         std::vector<std::tuple<int, int>> &assignmentC) {
      bool oneCentroidUpdated = false;
      int noC = centroids.size();
      std::vector<std::vector<ftm::FTMTree_MT *>> assignedTrees(noC);
      std::vector<std::vector<int>> assignedTreesIndex(noC);
      std::vector<std::vector<double>> assignedAlphas(noC);

      for(auto asgn : assignmentC) {
        assignedTrees[std::get<0>(asgn)].push_back(trees[std::get<1>(asgn)]);
        assignedTreesIndex[std::get<0>(asgn)].push_back(std::get<1>(asgn));
        assignedAlphas[std::get<0>(asgn)].push_back(alphas[std::get<1>(asgn)]);
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) \
  shared(centroids) if(parallelize_ and parallelizeUpdate_)
      {
#pragma omp single nowait
        {
#endif
          for(unsigned int i = 0; i < centroids.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i) shared(centroids)
            {
#endif
              // Init new centroid if no trees are assigned to it
              if(assignedTrees[i].size() == 0) {
                initNewCentroid<dataType>(trees, centroids[i]);
                for(unsigned int t = 0; t < trees.size(); ++t)
                  lowerBound_[t][i] = 0;
              } else if(assignedTrees[i].size() == 1) {
                centroids[i]
                  = copyMergeTree<dataType>(assignedTrees[i][0], true);
                cleanMergeTree<dataType>(centroids[i]);
              } else if(not samePreviousAssignment(
                          i)) { // Do not update if same previous assignment
                oneCentroidUpdated = true;
                // Otherwise compute barycenter of the assigned trees
                double alphasSum = 0;
                for(unsigned int j = 0; j < assignedAlphas[i].size(); ++j)
                  alphasSum += assignedAlphas[i][j];
                for(unsigned int j = 0; j < assignedAlphas[i].size(); ++j)
                  assignedAlphas[i][j] /= alphasSum;
                treesMatchingVector matching(assignedTrees[i].size());
                computeOneBarycenter<dataType>(
                  assignedTrees[i], centroids[i], assignedAlphas[i], matching);
                std::vector<ftm::idNode> deletedNodesT;
                persistenceThresholding<dataType>(
                  &(centroids[i].tree), 0, deletedNodesT);
                cleanMergeTree<dataType>(centroids[i]);
              }
#ifdef TTK_ENABLE_OPENMP
            } // pragma omp task
#endif
          }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
        } // pragma omp single nowait
      } // pragma omp parallel
#endif
      return oneCentroidUpdated;
    }

    template <class dataType>
    void computeOneBarycenter(
      std::vector<ftm::FTMTree_MT *> &trees,
      MergeTree<dataType> &baryMergeTree,
      std::vector<double> &alphas,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &finalMatchings) {
      MergeTreeBarycenter ftmTreeEditDistanceBary;
      ftmTreeEditDistanceBary.setDebugLevel(2);
      ftmTreeEditDistanceBary.setProgressiveComputation(false);
      ftmTreeEditDistanceBary.setBranchDecomposition(true);
      ftmTreeEditDistanceBary.setNormalizedWasserstein(normalizedWasserstein_);
      ftmTreeEditDistanceBary.setNormalizedWassersteinReg(
        normalizedWassersteinReg_);
      ftmTreeEditDistanceBary.setRescaledWasserstein(rescaledWasserstein_);
      ftmTreeEditDistanceBary.setKeepSubtree(keepSubtree_);
      ftmTreeEditDistanceBary.setAssignmentSolver(assignmentSolverID_);
      ftmTreeEditDistanceBary.setIsCalled(true);
      ftmTreeEditDistanceBary.setThreadNumber(this->threadNumber_);
      ftmTreeEditDistanceBary.setDistanceSquared(true); // squared root
      ftmTreeEditDistanceBary.setProgressiveBarycenter(progressiveBarycenter_);
      ftmTreeEditDistanceBary.setDeterministic(deterministic_);
      ftmTreeEditDistanceBary.setTol(tol_);

      ftmTreeEditDistanceBary.computeBarycenter<dataType>(
        trees, baryMergeTree, alphas, finalMatchings);

      addDeletedNodesTime_ += ftmTreeEditDistanceBary.getAddDeletedNodesTime();
    }

    // ----------------------------------------
    // Main Functions
    // ----------------------------------------
    template <class dataType>
    void computeCentroids(std::vector<ftm::FTMTree_MT *> &trees,
                          std::vector<MergeTree<dataType>> &centroids,
                          matchingVector &outputMatching,
                          std::vector<double> &alphas,
                          std::vector<int> &clusteringAssignment,
                          std::vector<ftm::FTMTree_MT *> &trees2,
                          std::vector<MergeTree<dataType>> &centroids2,
                          matchingVector &outputMatching2) {
      Timer t_clust;

      // Run
      int noCentroidsT = centroids.size();
      bool converged = false;
      dataType inertia = -1;
      dataType minInertia = std::numeric_limits<dataType>::max();
      int cptBlocked = 0;
      noIterationC_ = 0;
      std::vector<std::tuple<int, int>> assignmentC;
      std::vector<dataType> bestDistanceT(
        trees.size(), std::numeric_limits<dataType>::max());
      while(not converged) {
        ++noIterationC_;

        printMsg(debug::Separator::L1);
        std::stringstream ssIter;
        ssIter << "Itration " << noIterationC_;
        printMsg(ssIter.str());

        // --- Assignment
        Timer t_assignment;
        assignmentCentroids<dataType>(
          trees, centroids, assignmentC, bestDistanceT, trees2, centroids2);
        auto t_assignment_time = t_assignment.getElapsedTime();
        printMsg("Assignment", 1, t_assignment_time, this->threadNumber_);

        // --- Update
        Timer t_update;
        bool trees1Updated = true, trees2Updated = true;
        trees1Updated
          = updateCentroids<dataType>(trees, centroids, alphas, assignmentC);
        if(trees2.size() != 0)
          trees2Updated = updateCentroids<dataType>(
            trees2, centroids2, alphas, assignmentC);
        auto t_update_time = t_update.getElapsedTime();
        printMsg("Update", 1, t_update_time, this->threadNumber_);

        // --- Check convergence
        dataType currentInertia = 0;
        for(auto distance : bestDistanceT)
          currentInertia += distance * distance;
        converged = std::abs((double)(inertia - currentInertia)) < 0.01;
        inertia = currentInertia;
        std::stringstream ss3;
        ss3 << "Inertia : " << inertia;
        printMsg(ss3.str());

        minInertia = std::min(minInertia, inertia);
        if(not converged) {
          cptBlocked += (minInertia < inertia) ? 1 : 0;
          converged = (cptBlocked >= 10);
        }

        // Converged if barycenters were not updated (same assignment than last
        // iteration)
        converged = converged or (not trees1Updated and not trees2Updated);

        // --- Reset vectors
        if(not converged) {
          assignmentC.clear();
          bestDistanceT = std::vector<dataType>(
            trees.size(), std::numeric_limits<dataType>::max());
        }
      }

      // Final processing
      printMsg(debug::Separator::L1);
      printMsg("Final assignment");
      matchingVector matchingsC(noCentroidsT);
      matchingVector matchingsC2(noCentroidsT);
      finalAssignmentCentroids<dataType>(trees, centroids, matchingsC,
                                         assignmentC, bestDistanceT, trees2,
                                         centroids2, matchingsC2);
      for(auto dist : bestDistanceT)
        finalDistances_.push_back(dist);
      dataType currentInertia = 0;
      for(auto distance : bestDistanceT)
        currentInertia += distance * distance;
      std::stringstream ss;
      ss << "Inertia : " << currentInertia;
      printMsg(ss.str());

      // Manage output
      std::vector<int> cptCentroid(centroids.size(), 0);
      for(auto asgn : assignmentC) {
        int centroid = std::get<0>(asgn);
        int tree = std::get<1>(asgn);
        // std::cout << centroid << " " << tree << std::endl;
        clusteringAssignment[tree] = centroid;
        outputMatching[centroid][tree]
          = matchingsC[centroid][cptCentroid[centroid]];
        if(trees2.size() != 0)
          outputMatching2[centroid][tree]
            = matchingsC2[centroid][cptCentroid[centroid]];
        ++cptCentroid[centroid];
      }

      auto clusteringTime = t_clust.getElapsedTime() - addDeletedNodesTime_;
      printMsg("Total", 1, clusteringTime, this->threadNumber_);
    }

    template <class dataType>
    void execute(std::vector<MergeTree<dataType>> &trees,
                 matchingVector &outputMatching,
                 std::vector<double> &alphas,
                 std::vector<int> &clusteringAssignment,
                 std::vector<MergeTree<dataType>> &trees2,
                 matchingVector &outputMatching2,
                 std::vector<MergeTree<dataType>> &centroids) {
      // --- Preprocessing
      // std::vector<ftm::FTMTree_MT*> oldTrees, oldTrees2;
      treesNodeCorr_ = std::vector<std::vector<int>>(trees.size());
      preprocessingClustering<dataType>(trees, treesNodeCorr_);
      if(trees2.size() != 0) {
        trees2NodeCorr_ = std::vector<std::vector<int>>(trees2.size());
        preprocessingClustering<dataType>(trees2, trees2NodeCorr_, false);
      }
      std::vector<ftm::FTMTree_MT *> treesT;
      mergeTreeToFTMTree<dataType>(trees, treesT);
      std::vector<ftm::FTMTree_MT *> treesT2;
      mergeTreeToFTMTree<dataType>(trees2, treesT2);

      // --- Init centroids
      std::vector<std::vector<MergeTree<dataType>>> allCentroids;
      initCentroids<dataType>(treesT, treesT2, allCentroids);
      centroids = allCentroids[0];
      std::vector<MergeTree<dataType>> centroids2;
      if(trees2.size() != 0)
        centroids2 = allCentroids[1];
      /*for(unsigned int i = 0; i < centroids.size(); ++i){
        verifyBranchDecompositionInconsistency<dataType>(centroids[i]->tree);
        if(trees2.size() != 0)
          verifyBranchDecompositionInconsistency<dataType>(centroids2[i]->tree);
      }*/

      // --- Init accelerated kmeans
      initAcceleratedKMeansVectors<dataType>(treesT, centroids, treesT2);

      // --- Execute
      computeCentroids<dataType>(treesT, centroids, outputMatching, alphas,
                                 clusteringAssignment, treesT2, centroids2,
                                 outputMatching2);

      // --- Postprocessing
      if(postprocess_) {
        fixMergedRootOriginClustering<dataType>(centroids);
        postprocessingClustering<dataType>(
          trees, centroids, outputMatching, clusteringAssignment);
        // TODO fix min max pair in trees2
        /*if(trees2.size() != 0){
          putBackMinMaxPair<dataType>(centroids, centroids2);
          postprocessingClustering<dataType>(trees2, centroids2,
        outputMatching2, clusteringAssignment);
        }*/
      }
    }

    template <class dataType>
    void execute(std::vector<MergeTree<dataType>> &trees,
                 matchingVector &outputMatching,
                 std::vector<int> &clusteringAssignment,
                 std::vector<MergeTree<dataType>> &trees2,
                 matchingVector &outputMatching2,
                 std::vector<MergeTree<dataType>> &centroids) {
      if(trees2.size() != 0)
        printMsg("Use join and split trees");

      std::vector<double> alphas;
      for(unsigned int i = 0; i < trees.size(); ++i)
        alphas.push_back(1.0 / trees.size());

      execute<dataType>(trees, outputMatching, alphas, clusteringAssignment,
                        trees2, outputMatching2, centroids);
    }

    template <class dataType>
    void execute(std::vector<MergeTree<dataType>> &trees,
                 matchingVector &outputMatching,
                 std::vector<int> &clusteringAssignment,
                 std::vector<MergeTree<dataType>> &centroids) {
      std::vector<MergeTree<dataType>> trees2
        = std::vector<MergeTree<dataType>>();
      matchingVector outputMatching2 = matchingVector();
      execute<dataType>(trees, outputMatching, clusteringAssignment, trees2,
                        outputMatching2, centroids);
    }

    // ----------------------------------------
    // Preprocessing
    // ----------------------------------------
    template <class dataType>
    void preprocessingClustering(std::vector<MergeTree<dataType>> &trees,
                                 std::vector<std::vector<int>> &nodeCorr,
                                 bool useMinMaxPairT = true) {
      for(unsigned int i = 0; i < trees.size(); ++i) {
        preprocessingPipeline<dataType>(
          trees[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          branchDecomposition_, useMinMaxPairT, cleanTree_, nodeCorr[i]);
        if(trees.size() < 40)
          printTreeStats(trees[i]);
      }
      printTreesStats(trees);
    }

    // ----------------------------------------
    // Postprocessing
    // ----------------------------------------
    template <class dataType>
    void fixMergedRootOriginClustering(
      std::vector<MergeTree<dataType>> &centroids) {
      for(unsigned int i = 0; i < centroids.size(); ++i)
        fixMergedRootOriginBarycenter<dataType>(centroids[i]);
    }

    template <class dataType>
    void updateNodesAndScalars(
      MergeTree<dataType> &mTree1,
      std::vector<dataType> &newScalarsVector,
      std::vector<std::tuple<ftm::idNode, bool>> &origins) {
      ftm::FTMTree_MT *tree1 = &(mTree1.tree);

      // Create new tree
      MergeTree<dataType> mTreeNew
        = createEmptyMergeTree<dataType>(newScalarsVector.size());
      setTreeScalars<dataType>(mTreeNew, newScalarsVector);
      ftm::FTMTree_MT *treeNew = &(mTreeNew.tree);

      // Copy the old tree structure
      copyMergeTreeStructure(&(mTreeNew.tree), tree1);

      // Manage new origins
      int noNodes = treeNew->getNumberOfNodes();
      for(int i = noNodes; i < newScalarsVector.size(); ++i) {
        treeNew->makeNode(i);
        auto tup = origins[i - noNodes];
        ftm::idNode origin = std::get<0>(tup);
        bool reverse = std::get<1>(tup);
        treeNew->getNode(i)->setOrigin(origin);
        if(reverse)
          treeNew->getNode(origin)->setOrigin(i);
      }

      // Copy new tree
      mTree1 = mTreeNew;
    }

    // TODO manage if join tree?
    template <class dataType>
    void putBackMinMaxPair(std::vector<MergeTree<dataType>> &centroids,
                           std::vector<MergeTree<dataType>> &centroids2) {
      for(unsigned int i = 0; i < centroids2.size(); ++i) {
        // Get min max pair
        ftm::FTMTree_MT *centroidITree = &(centroids[i].tree);
        ftm::idNode root = centroidITree->getRoot();
        dataType newMax = centroidITree->getValue<dataType>(root);
        ftm::idNode rootOrigin = centroidITree->getNode(root)->getOrigin();
        dataType newMin = centroidITree->getValue<dataType>(rootOrigin);

        // Update tree
        ftm::idNode root2 = centroids2[i].tree.getRoot();
        std::vector<dataType> newScalarsVector
          = getTreeScalars<dataType>(centroids2[i]);
        newScalarsVector[root2] = newMax;
        newScalarsVector.push_back(newMin);
        std::vector<std::tuple<ftm::idNode, bool>> origins;
        origins.push_back(std::make_tuple(root2, true));
        updateNodesAndScalars<dataType>(
          centroids2[i], newScalarsVector, origins);
      }
    }

    template <class dataType>
    void postprocessingClustering(std::vector<MergeTree<dataType>> &trees,
                                  std::vector<MergeTree<dataType>> &centroids,
                                  matchingVector &outputMatching,
                                  std::vector<int> &clusteringAssignment) {
      for(unsigned int i = 0; i < trees.size(); ++i)
        postprocessingPipeline<dataType>(&(trees[i].tree));
      for(unsigned int i = 0; i < centroids.size(); ++i)
        postprocessingPipeline<dataType>(&(centroids[i].tree));
      for(unsigned int c = 0; c < centroids.size(); ++c)
        for(unsigned int i = 0; i < trees.size(); ++i)
          if(clusteringAssignment[i] == (int)c)
            convertBranchDecompositionMatching<dataType>(
              &(centroids[c].tree), &(trees[i].tree), outputMatching[c][i]);
    }

  }; // MergeTreeClustering class

} // namespace ttk

#endif
