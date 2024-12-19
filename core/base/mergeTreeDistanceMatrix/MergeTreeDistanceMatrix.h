/// \ingroup base
/// \class ttk::MergeTreeDistanceMatrix
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2021.
///
/// This VTK filter uses the ttk::MergeTreeDistanceMatrix module to compute the
/// distance matrix of a group of merge trees.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \b Related \b publication \n
/// "Edit Distance between Merge Trees" \n
/// R. Sridharamurthy, T. B. Masood, A. Kamakshidasan and V. Natarajan. \n
/// IEEE Transactions on Visualization and Computer Graphics, 2020.
///
/// \b Related \b publication \n
/// "Branch Decomposition-Independent Edit Distances for Merge Trees." \n
/// Florian Wetzels, Heike Leitte, and Christoph Garth. \n
/// Computer Graphics Forum, 2022.
///
/// \b Related \b publication \n
/// "A Deformation-based Edit Distance for Merge Trees" \n
/// Florian Wetzels, Christoph Garth. \n
/// TopoInVis 2022.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramPGA/">Persistence
///   Diagram Principal Geodesic Analysis example</a> \n

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
      treesNodeCorr_.resize(trees.size());
      for(unsigned int i = 0; i < trees.size(); ++i) {
        preprocessingPipeline<dataType>(
          trees[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_, true,
          treesNodeCorr_[i], true, baseModule_ == 2);
      }
      executePara<dataType>(trees, distanceMatrix);
      if(trees2.size() != 0) {
        std::vector<std::vector<int>> trees2NodeCorr(trees2.size());
        for(unsigned int i = 0; i < trees.size(); ++i) {
          preprocessingPipeline<dataType>(
            trees2[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
            baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_,
            true, treesNodeCorr_[i], true, baseModule_ == 2);
        }
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
              mergeTreeDistance.setPreprocess(false);
              // mergeTreeDistance.setSaveTree(true);
              mergeTreeDistance.setSaveTree(false);
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
            } else if(baseModule_ == 1) {
              BranchMappingDistance branchDist;
              branchDist.setBaseMetric(branchMetric_);
              branchDist.setAssignmentSolver(assignmentSolverID_);
              branchDist.setSquared(distanceSquaredRoot_);
              branchDist.setEpsilonTree1(epsilonTree1_);
              branchDist.setEpsilonTree2(epsilonTree2_);
              branchDist.setEpsilon2Tree1(epsilon2Tree1_);
              branchDist.setEpsilon2Tree2(epsilon2Tree2_);
              branchDist.setEpsilon3Tree1(epsilon3Tree1_);
              branchDist.setEpsilon3Tree2(epsilon3Tree2_);
              branchDist.setPersistenceThreshold(persistenceThreshold_);
              branchDist.setPreprocess(false);
              // branchDist.setSaveTree(true);
              branchDist.setSaveTree(false);
              dataType dist = branchDist.execute<dataType>(trees[i], trees[j]);
              distanceMatrix[i][j] = static_cast<double>(dist);
            } else if(baseModule_ == 2) {
              PathMappingDistance pathDist;
              pathDist.setBaseMetric(pathMetric_);
              pathDist.setAssignmentSolver(assignmentSolverID_);
              pathDist.setSquared(distanceSquaredRoot_);
              pathDist.setComputeMapping(true);
              pathDist.setEpsilonTree1(epsilonTree1_);
              pathDist.setEpsilonTree2(epsilonTree2_);
              pathDist.setEpsilon2Tree1(epsilon2Tree1_);
              pathDist.setEpsilon2Tree2(epsilon2Tree2_);
              pathDist.setEpsilon3Tree1(epsilon3Tree1_);
              pathDist.setEpsilon3Tree2(epsilon3Tree2_);
              pathDist.setPersistenceThreshold(persistenceThreshold_);
              pathDist.setPreprocess(false);
              // pathDist.setSaveTree(true);
              pathDist.setSaveTree(false);
              dataType dist = pathDist.execute<dataType>(trees[i], trees[j]);
              distanceMatrix[i][j] = static_cast<double>(dist);
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
