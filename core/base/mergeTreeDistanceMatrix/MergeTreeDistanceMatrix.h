/// \ingroup base
/// \class ttk::MergeTreeDistanceMatrix
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021.
///
/// This VTK filter uses the ttk::MergeTreeDistanceMatrix module to compute the
/// distance matrix of a group of merge trees.

#pragma once

// ttk common includes
#include <Debug.h>

#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>

namespace ttk {

  /**
   * The MergeTreeDistanceMatrix class provides methods to compute the
   * distance between multiple merge trees and output a distance matrix.
   */
  class MergeTreeDistanceMatrix : virtual public Debug,
                                  virtual public MergeTreeBase {
  public:
    MergeTreeDistanceMatrix() {
      this->setDebugMsgPrefix(
        "MergeTreeDistanceMatrix"); // inherited from Debug: prefix will be
                                    // printed at the
      // beginning of every msg
    }
    ~MergeTreeDistanceMatrix() override = default;

    /**
     * Implementation of the algorithm.
     */
    template <class dataType>
    void execute(std::vector<ftm::MergeTree<dataType>> &trees,
                 std::vector<ftm::MergeTree<dataType>> &trees2,
                 std::vector<std::vector<double>> &distanceMatrix) {
      executePara<dataType>(trees, distanceMatrix);
      if(trees2.size() != 0) {
        useDoubleInput_ = true;
        std::vector<std::vector<double>> distanceMatrix2(
          trees2.size(), std::vector<double>(trees2.size()));
        executePara<dataType>(trees2, distanceMatrix2, false);
        mixDistancesMatrix(distanceMatrix, distanceMatrix2);
      }
    }

    template <class dataType>
    void executePara(std::vector<ftm::MergeTree<dataType>> &trees,
                     std::vector<std::vector<double>> &distanceMatrix,
                     bool isFirstInput = true) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp single nowait
#endif
        executeParaImpl<dataType>(trees, distanceMatrix, isFirstInput);
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void executeParaImpl(std::vector<ftm::MergeTree<dataType>> &trees,
                         std::vector<std::vector<double>> &distanceMatrix,
                         bool isFirstInput = true) {
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(i) UNTIED() shared(distanceMatrix, trees)
        {
#endif
          if(i % std::max(int(distanceMatrix.size() / 10), 1) == 0) {
            std::stringstream stream;
            stream << i << " / " << distanceMatrix.size();
            printMsg(stream.str());
          }

          distanceMatrix[i][i] = 0.0;
          for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
            // Execute
            MergeTreeDistance mergeTreeDistance;
            mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
            mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
            mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
            mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
            mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
            mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
            mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
            mergeTreeDistance.setProgressiveComputation(
              progressiveComputation_);
            mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
            mergeTreeDistance.setParallelize(parallelize_);
            mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
            mergeTreeDistance.setDebugLevel(2);
            mergeTreeDistance.setThreadNumber(this->threadNumber_);
            mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
            mergeTreeDistance.setRescaledWasserstein(rescaledWasserstein_);
            mergeTreeDistance.setNormalizedWassersteinReg(
              normalizedWassersteinReg_);
            mergeTreeDistance.setKeepSubtree(keepSubtree_);
            mergeTreeDistance.setDistanceSquared(distanceSquared_);
            mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
            mergeTreeDistance.setSaveTree(true);
            mergeTreeDistance.setCleanTree(true);
            mergeTreeDistance.setIsCalled(true);
            mergeTreeDistance.setPostprocess(false);
            if(useDoubleInput_) {
              double weight = mixDistancesMinMaxPairWeight(isFirstInput);
              mergeTreeDistance.setMinMaxPairWeight(weight);
              mergeTreeDistance.setDistanceSquared(true);
            }
            std::vector<std::tuple<ftm::idNode, ftm::idNode>> outputMatching;
            distanceMatrix[i][j] = mergeTreeDistance.execute<dataType>(
              trees[i], trees[j], outputMatching);
            // distance matrix is symmetric
            distanceMatrix[j][i] = distanceMatrix[i][j];
          } // end for j
#ifdef TTK_ENABLE_OPENMP
        } // end task
#endif
      } // end for i
    }

  }; // MergeTreeDistanceMatrix class

} // namespace ttk
