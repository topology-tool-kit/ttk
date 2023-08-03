/// \ingroup base
/// \class ttk::MergeTreeDistanceMatrix
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021.
///
/// This VTK filter uses the ttk::MergeTreeDistanceMatrix module to compute the
/// distance matrix of a group of merge trees.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>

#include <BranchMappingDistance.h>
#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>
#include <PathMappingDistance.h>

namespace ttk {

  /**
   * The MergeTreeDistanceMatrix class provides methods to compute the
   * distance between multiple merge trees and output a distance matrix.
   */
  class MergeTreeDistanceMatrix : virtual public Debug,
                                  virtual public MergeTreeBase {
  protected:
    int baseModule_ = 0;
    int branchMetric_ = 0;
    int pathMetric_ = 0;

  public:
    MergeTreeDistanceMatrix() {
      this->setDebugMsgPrefix(
        "MergeTreeDistanceMatrix"); // inherited from Debug: prefix will be
                                    // printed at the
      // beginning of every msg
    }
    ~MergeTreeDistanceMatrix() override = default;

    void setBaseModule(int m) {
      baseModule_ = m;
    }

    void setBranchMetric(int m) {
      branchMetric_ = m;
    }

    void setPathMetric(int m) {
      pathMetric_ = m;
    }

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
    void execute(std::vector<ftm::MergeTree<dataType>> &ftmtrees,
                 std::vector<std::vector<double>> &distanceMatrix) {
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
        if(i % std::max(int(distanceMatrix.size() / 10), 1) == 0) {
          std::stringstream stream;
          stream << i << " / " << distanceMatrix.size();
          printMsg(stream.str());
        }

        BranchMappingDistance branchDist;
        branchDist.setBaseMetric(branchMetric_);
        branchDist.setAssignmentSolver(assignmentSolverID_);
        branchDist.setSquared(not distanceSquaredRoot_);
        PathMappingDistance pathDist;
        pathDist.setBaseMetric(pathMetric_);
        pathDist.setAssignmentSolver(assignmentSolverID_);
        pathDist.setSquared(not distanceSquaredRoot_);
        pathDist.setComputeMapping(true);

        distanceMatrix[i][i] = 0.0;
        // compareTrees(trees[i],&(ftmtrees[i].tree));
        for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
          // Execute
          if(baseModule_ == 0) {
            distanceMatrix[i][j] = 0;
          } else if(baseModule_ == 1) {
            dataType dist = branchDist.editDistance_branch<dataType>(
              &(ftmtrees[i].tree), &(ftmtrees[j].tree));
            distanceMatrix[i][j] = static_cast<double>(dist);
          } else if(baseModule_ == 2) {
            dataType dist = pathDist.editDistance_path<dataType>(
              &(ftmtrees[i].tree), &(ftmtrees[j].tree));
            distanceMatrix[i][j] = static_cast<double>(dist);
          }
          // distance matrix is symmetric
          distanceMatrix[j][i] = distanceMatrix[i][j];
        } // end for j
      } // end for i
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
            if(baseModule_ == 0) {
              MergeTreeDistance mergeTreeDistance;
              mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
              mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
              mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
              mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
              mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
              mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
              mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
              mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
              mergeTreeDistance.setParallelize(parallelize_);
              mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
              mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
              mergeTreeDistance.setThreadNumber(this->threadNumber_);
              mergeTreeDistance.setNormalizedWasserstein(
                normalizedWasserstein_);
              mergeTreeDistance.setKeepSubtree(keepSubtree_);
              mergeTreeDistance.setDistanceSquaredRoot(distanceSquaredRoot_);
              mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
              mergeTreeDistance.setSaveTree(true);
              mergeTreeDistance.setCleanTree(true);
              mergeTreeDistance.setIsCalled(true);
              mergeTreeDistance.setPostprocess(false);
              mergeTreeDistance.setIsPersistenceDiagram(isPersistenceDiagram_);
              if(useDoubleInput_) {
                double const weight
                  = mixDistancesMinMaxPairWeight(isFirstInput);
                mergeTreeDistance.setMinMaxPairWeight(weight);
                mergeTreeDistance.setDistanceSquaredRoot(true);
              }
              std::vector<std::tuple<ftm::idNode, ftm::idNode>> outputMatching;
              distanceMatrix[i][j] = mergeTreeDistance.execute<dataType>(
                trees[i], trees[j], outputMatching);
            } else {
              distanceMatrix[i][j] = 0;
            }
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
