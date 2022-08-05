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
    int baseModule = 0;
    int branchMetric = 0;
    int pathMetric = 0;

  public:
    MergeTreeDistanceMatrix() {
      this->setDebugMsgPrefix(
        "MergeTreeDistanceMatrix"); // inherited from Debug: prefix will be
                                    // printed at the
      // beginning of every msg
    }
    ~MergeTreeDistanceMatrix() override = default;

    void setBaseModule(int m) {
      baseModule = m;
    }

    void setBranchMetric(int m) {
      branchMetric = m;
    }

    void setPathMetric(int m) {
      pathMetric = m;
    }

    template <class dataType>
    void compareTrees(
      std::tuple<std::vector<dataType>, std::vector<std::vector<int>>, int>
        &adjList,
      ftm::FTMTree_MT* ftmtree) {
      auto nodeScalars = std::get<0>(adjList);
      auto childrenLists = std::get<1>(adjList);
      auto rootID = std::get<2>(adjList);

      std::cout << "\n-------------------\n";
      std::cout << "Adjacency list:\n";
      for(int nIdx=0; nIdx<nodeScalars.size(); nIdx++){
        std::cout << nIdx << ": " << nodeScalars[nIdx] << ", ";
        for(int cIdx : childrenLists[nIdx]){
          std::cout << cIdx << " ";
        }
        std::cout << "\n";
      }
      std::cout << "\n-------------------\n" << std::endl;
      std::cout << "FTM:\n";
      for(int nIdx=0; nIdx<ftmtree->getNumberOfNodes(); nIdx++){
        double v = ftmtree->getValue<dataType>(nIdx);
        std::cout << nIdx << ": " << v << ", ";
        std::vector<ftm::idNode> children;
        ftmtree->getChildren(nIdx,children);
        for(int cIdx : children){
          std::cout << cIdx << " ";
        }
        std::cout << "\n";
      }
      std::cout << "\n-------------------\n" << std::endl;
      std::stack<std::pair<int,int>> stack;
      stack.push(std::make_pair(rootID,ftmtree->getRoot()));
      while(!stack.empty()) {
        int n1 = stack.top().first;
        int n2 = stack.top().second;
        stack.pop();
        std::cout << nodeScalars[n1] << " " << ftmtree->getValue<dataType>(n2) << "\n";
        auto children1 = childrenLists[n1];
        if(children1.size()==2){
          if(nodeScalars[children1[0]]>nodeScalars[children1[1]]){
            int temp = children1[0];
            children1[0] = children1[1];
            children1[0] = temp;
          }
        }
        std::vector<ftm::idNode> children2;
        ftmtree->getChildren(n2,children2);
        if(children2.size()==2){
          auto v1 = ftmtree->getValue<dataType>(children2[0]);
          auto v2 = ftmtree->getValue<dataType>(children2[1]);
          if(v1>v2){
            int temp = children2[0];
            children2[0] = children2[1];
            children2[0] = temp;
          }
          stack.push(std::make_pair(children1[0],children2[0]));
          stack.push(std::make_pair(children1[1],children2[1]));
        }
        if(children2.size()==1){
          
          stack.push(std::make_pair(children1[0],children2[0]));
        }
      }
      std::cout << "\n-------------------\n" << std::endl;
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
    void execute(
      std::vector<
        std::tuple<std::vector<dataType>, std::vector<std::vector<int>>, int>>
        &trees,
      std::vector<ftm::MergeTree<dataType>> &ftmtrees,
      std::vector<std::vector<double>> &distanceMatrix) {
      for(unsigned int i = 0; i < distanceMatrix.size(); ++i) {
        if(i % std::max(int(distanceMatrix.size() / 10), 1) == 0) {
          std::stringstream stream;
          stream << i << " / " << distanceMatrix.size();
          printMsg(stream.str());
        }

        BranchMappingDistance branchDist;
        branchDist.setBaseMetric(branchMetric);
        branchDist.setAssignmentSolver(assignmentSolverID_);
        branchDist.setSquared(distanceSquared_);
        PathMappingDistance pathDist;
        pathDist.setBaseMetric(pathMetric);
        pathDist.setAssignmentSolver(assignmentSolverID_);
        pathDist.setSquared(distanceSquared_);

        distanceMatrix[i][i] = 0.0;
        //compareTrees(trees[i],&(ftmtrees[i].tree));
        for(unsigned int j = i + 1; j < distanceMatrix[0].size(); ++j) {
          // Execute
          if(baseModule == 0) {
            distanceMatrix[i][j] = 0;
          } else if(baseModule == 1) {
            auto nodes1 = std::get<0>(trees[i]);
            auto topo1 = std::get<1>(trees[i]);
            auto rootID1 = std::get<2>(trees[i]);
            auto nodes2 = std::get<0>(trees[j]);
            auto topo2 = std::get<1>(trees[j]);
            auto rootID2 = std::get<2>(trees[j]);
            auto start = std::chrono::high_resolution_clock::now();
            dataType dist = branchDist.editDistance_branch<dataType>(
              nodes1, topo1, rootID1, nodes2, topo2, rootID2);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            start = std::chrono::high_resolution_clock::now();
            dataType dist_ftm = branchDist.editDistance_branch<dataType>(&(ftmtrees[i].tree),&(ftmtrees[j].tree));
            end = std::chrono::high_resolution_clock::now();
            auto duration_ftm = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            distanceMatrix[i][j] = static_cast<double>(dist);
            std::cout << dist << " ; " << dist_ftm  << std::endl;
            std::cout << duration.count()*0.000001 << " seconds" << " ; " << duration_ftm.count()*0.000001 << " seconds"  << std::endl;
          } else if(baseModule == 2) {
            auto nodes1 = std::get<0>(trees[i]);
            auto topo1 = std::get<1>(trees[i]);
            auto rootID1 = std::get<2>(trees[i]);
            auto nodes2 = std::get<0>(trees[j]);
            auto topo2 = std::get<1>(trees[j]);
            auto rootID2 = std::get<2>(trees[j]);
            auto start = std::chrono::high_resolution_clock::now();
            dataType dist = pathDist.editDistance_path<dataType>(
              nodes1, topo1, rootID1, nodes2, topo2, rootID2);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            start = std::chrono::high_resolution_clock::now();
            dataType dist_ftm = pathDist.editDistance_path<dataType>(&(ftmtrees[i].tree),&(ftmtrees[j].tree));
            end = std::chrono::high_resolution_clock::now();
            auto duration_ftm = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            distanceMatrix[i][j] = static_cast<double>(dist);
            std::cout << dist << " ; " << dist_ftm  << std::endl;
            std::cout << duration.count()*0.000001 << " seconds" << " ; " << duration_ftm.count()*0.000001 << " seconds"  << std::endl;
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
            if(baseModule == 0) {
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
              mergeTreeDistance.setNormalizedWasserstein(
                normalizedWasserstein_);
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
