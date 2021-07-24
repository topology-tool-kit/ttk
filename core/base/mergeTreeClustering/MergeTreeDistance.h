/// \ingroup base
/// \class MergeTreeDistance
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// This module defines the %MergeTreeDistance class that computes distance
/// between two merge trees.
///
/// \b Related \b publication \n
/// "Wasserstein Distances, Geodesics and Barycenters of Merge Trees" \n
/// Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.\n
/// Proc. of IEEE VIS 2021.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021

#ifndef _MERGETREEDISTANCE_H
#define _MERGETREEDISTANCE_H

#pragma once

#include <stack>
#include <thread>

// ttk common includes
#include <Debug.h>

#include "MergeTreeBase.h"
#include <AssignmentAuction.h>
#include <AssignmentExhaustive.h>
#include <AssignmentMunkres.h>
#include <AssignmentSolver.h>

namespace ttk {

  /**
   * The MergeTreeDistance class provides methods to compute distance
   * between two merge trees.
   */
  class MergeTreeDistance : virtual public Debug, public MergeTreeBase {

  private:
    double t_assignment_time_ = 0;

    bool preprocess_ = true;
    bool postprocess_ = true;
    bool saveTree_ = false;
    bool onlyEmptyTreeDistance_ = false;

    bool isCalled_ = false;

    double auctionEpsilon_ = -1;
    double auctionEpsilonDiviser_ = 0;
    int auctionRound_ = -1;

    // Just to get some stats about run
    // std::map<int, int> assignmentProblemSize, assignmentProblemIter;

    bool testing_ = true;

    // Parallel version data
    std::vector<std::vector<ftm::idNode>> tree2LevelToNode_;
    std::vector<int> tree1Level_, tree2Level_;

  public:
    MergeTreeDistance() {
      this->setDebugMsgPrefix(
        "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                              // the beginning of every msg
#ifdef TTK_ENABLE_OPENMP
      omp_set_nested(1);
#endif
    };
    ~MergeTreeDistance(){};

    void setIsCalled(bool ic) {
      isCalled_ = ic;
    }

    void setPreprocess(bool preproc) {
      if(not progressiveComputation_)
        preprocess_ = preproc;
    }

    void setPostprocess(bool postproc) {
      postprocess_ = postproc;
    }

    void setTesting(bool test) {
      testing_ = test;
    }

    void setSaveTree(bool save) {
      saveTree_ = save;
    }

    void setAuctionEpsilon(double aucEpsilon) {
      auctionEpsilon_ = aucEpsilon;
    }

    void setAuctionEpsilonDiviser(double aucEpsilonDiviser) {
      auctionEpsilonDiviser_ = aucEpsilonDiviser;
    }

    void setAuctionNoRounds(double aucNoRounds) {
      auctionRound_ = aucNoRounds;
    }

    void setOnlyEmptyTreeDistance(double only) {
      onlyEmptyTreeDistance_ = only;
    }

    /**
     * Implementation of the algorithm.
     */

    // ----------------------------------------
    // Assignment Problem
    // ----------------------------------------
    template <class dataType>
    void
      runAssignmentProblemSolver(std::vector<std::vector<dataType>> &costMatrix,
                                 std::vector<asgnMatchingTuple> &matchings) {
      AssignmentSolver<dataType> *assignmentSolver;
      AssignmentExhaustive<dataType> solverExhaustive;
      AssignmentMunkres<dataType> solverMunkres;
      AssignmentAuction<dataType> solverAuction;
      switch(assignmentSolverID_) {
        case 1:
          solverExhaustive = AssignmentExhaustive<dataType>();
          assignmentSolver = &solverExhaustive;
          break;
        case 2:
          solverMunkres = AssignmentMunkres<dataType>();
          assignmentSolver = &solverMunkres;
          break;
        case 0:
        default:
          solverAuction = AssignmentAuction<dataType>();
          solverAuction.setEpsilon(auctionEpsilon_);
          solverAuction.setEpsilonDiviserMultiplier(auctionEpsilonDiviser_);
          solverAuction.setNumberOfRounds(auctionRound_);
          assignmentSolver = &solverAuction;
      }
      assignmentSolver->setInput(costMatrix);
      assignmentSolver->setBalanced(false);
      assignmentSolver->run(matchings);
    }

    template <class dataType>
    void createCostMatrix(std::vector<std::vector<dataType>> &treeTable,
                          std::vector<ftm::idNode> &children1,
                          std::vector<ftm::idNode> &children2,
                          std::vector<std::vector<dataType>> &costMatrix) {
      unsigned int nRows = children1.size(), nCols = children2.size();
      for(unsigned int i = 0; i < nRows; ++i) {
        int forestTableI = children1[i] + 1;
        for(unsigned int j = 0; j < nCols; ++j) {
          int forestTableJ = children2[j] + 1;
          // Cost of assigning i and j
          costMatrix[i][j] = treeTable[forestTableI][forestTableJ];
          if(tree1Level_[children1[i]] != tree2Level_[children2[j]]
             and not keepSubtree_)
            printErr("different levels!"); // should be impossible
        }
        // Cost of not assigning i
        costMatrix[i][nCols] = treeTable[forestTableI][0];
      }
      for(unsigned int j = 0; j < nCols; ++j) {
        int forestTableJ = children2[j] + 1;
        // Cost of not assigning j
        costMatrix[nRows][j] = treeTable[0][forestTableJ];
      }
      costMatrix[nRows][nCols] = 0;
    }

    template <class dataType>
    dataType forestAssignmentProblem(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<ftm::idNode> &children1,
      std::vector<ftm::idNode> &children2,
      std::vector<std::tuple<int, int>> &forestAssignment) {
      // --- Create cost matrix
      int nRows = children1.size(), nCols = children2.size();
      std::vector<std::vector<dataType>> costMatrix(
        nRows + 1, std::vector<dataType>(nCols + 1));
      createCostMatrix(treeTable, children1, children2, costMatrix);

      // assignmentProblemSize[costMatrix.size()*costMatrix[0].size()]++;

      // --- Solve assignment problem
      std::vector<asgnMatchingTuple> matchings;
      runAssignmentProblemSolver(costMatrix, matchings);

      // --- Postprocess matching to create output assignment
      dataType cost = postprocessAssignment<dataType>(
        matchings, children1, children2, forestAssignment);

      return cost;
    }

    template <class dataType>
    void computeEquation13(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      int i,
      int j,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      std::vector<ftm::idNode> &children1,
      std::vector<ftm::idNode> &children2) {
      if(children1.size() != 0 && children2.size() != 0) {
        dataType forestTerm1, forestTerm2, forestTerm3;
        std::tuple<dataType, ftm::idNode> forestCoTerm1, forestCoTerm2;
        // Term 1
        forestCoTerm1
          = computeTerm1_2<dataType>(children2, i, forestTable, true);
        forestTerm1 = forestTable[0][j] + std::get<0>(forestCoTerm1);

        // Term2
        forestCoTerm2
          = computeTerm1_2<dataType>(children1, j, forestTable, false);
        forestTerm2 = forestTable[i][0] + std::get<0>(forestCoTerm2);

        // Term 3
        Timer t_assignment;
        std::vector<std::tuple<int, int>> forestAssignment;
        forestTerm3 = forestAssignmentProblem<dataType>(
          tree1, tree2, treeTable, children1, children2, forestAssignment);
        if(not parallelize_)
          t_assignment_time_ += t_assignment.getElapsedTime();

        // Compute table value
        forestTable[i][j] = keepSubtree_ ? std::min(
                              std::min(forestTerm1, forestTerm2), forestTerm3)
                                         : forestTerm3;

        // Add backtracking information
        if(forestTable[i][j] == forestTerm3) {
          forestBackTable[i][j] = forestAssignment;
        } else if(forestTable[i][j] == forestTerm2) {
          forestBackTable[i][j].push_back(
            std::make_tuple(std::get<1>(forestCoTerm2), j));
        } else {
          forestBackTable[i][j].push_back(
            std::make_tuple(i, std::get<1>(forestCoTerm1)));
        }
      } else {
        // If one of the forest is empty we get back to equation 8 or 10
        forestTable[i][j]
          = (children1.size() == 0) ? forestTable[0][j] : forestTable[i][0];
      }
    }

    // ----------------------------------------
    // Main Functions
    // ----------------------------------------
    template <class dataType>
    dataType
      computeDistance(FTMTree_MT *tree1,
                      FTMTree_MT *tree2,
                      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                        &outputMatching) {
      // ---------------------
      // ----- Init dynamic progamming tables
      // --------------------
      size_t nRows = tree1->getNumberOfNodes() + 1;
      size_t nCols = tree2->getNumberOfNodes() + 1;
      std::vector<std::vector<dataType>> treeTable(
        nRows, std::vector<dataType>(nCols));
      std::vector<std::vector<dataType>> forestTable(
        nRows, std::vector<dataType>(nCols));

      // Backtracking tables (output matching)
      std::vector<std::vector<std::tuple<int, int>>> treeBackTable(
        nRows, std::vector<std::tuple<int, int>>(nCols));
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        forestBackTable(
          nRows, std::vector<std::vector<std::tuple<int, int>>>(nCols));

      int indR = tree1->getRoot() + 1;
      int indC = tree2->getRoot() + 1;

      tree1->getAllNodeLevel(tree1Level_);
      tree2->getAllNodeLevel(tree2Level_);
      tree2->getLevelToNode(tree2LevelToNode_);

      // ---------------------
      // ----- Compute edit distance
      // --------------------
      if(progressiveComputation_) {
        /*MergeTreeDistanceProgressive<dataType> editDistanceProgressive;
        editDistanceProgressive.computeEditDistanceProgressive(tree1, tree2,
                            treeTable, forestTable, treeBackTable,
        forestBackTable);*/
      } else {
        computeEditDistance(tree1, tree2, treeTable, forestTable, treeBackTable,
                            forestBackTable, nRows, nCols);
      }
      dataType distance = treeTable[indR][indC];
      if(onlyEmptyTreeDistance_)
        distance = treeTable[indR][0];
      if(distanceSquared_)
        distance = std::sqrt(distance);

      // ---------------------
      // ----- Compute matching
      // --------------------
      computeMatching<dataType>(tree1, tree2, treeBackTable, forestBackTable,
                                outputMatching, indR, indC);

      return distance;
    }

    template <class dataType>
    dataType computeDistance(
      FTMTree_MT *tree1,
      FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        realOutputMatching;
      dataType res
        = computeDistance<dataType>(tree1, tree2, realOutputMatching);
      for(auto tup : realOutputMatching)
        outputMatching.push_back(
          std::make_tuple(std::get<0>(tup), std::get<1>(tup)));
      return res;
    }

    template <class dataType>
    dataType execute(MergeTree<dataType> &mTree1,
                     MergeTree<dataType> &mTree2,
                     std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                       &outputMatching) {
      Memory m;
      Timer t_total;

      ftm::FTMTree_MT *tree1 = &(mTree1.tree);
      ftm::FTMTree_MT *tree2 = &(mTree2.tree);

      // ---------------------
      // ----- Testing
      // --------------------
      testing_ = false;

      // ---------------------
      // ----- Preprocessing
      // --------------------
      MergeTree<dataType> tree1Ori = mTree1;
      MergeTree<dataType> tree2Ori = mTree2;
      if(saveTree_) {
        mTree1 = copyMergeTree<dataType>(mTree1);
        mTree2 = copyMergeTree<dataType>(mTree2);
        tree1 = &(mTree1.tree);
        tree2 = &(mTree2.tree);
      }
      if(not isCalled_) {
        verifyMergeTreeStructure<dataType>(tree1);
        verifyMergeTreeStructure<dataType>(tree2);
      }
      if(preprocess_) {
        treesNodeCorr_ = std::vector<std::vector<int>>(2);
        preprocessingPipeline<dataType>(
          mTree1, epsilonTree1_, epsilon2Tree1_, epsilon3Tree1_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[0]);
        preprocessingPipeline<dataType>(
          mTree2, epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[1]);
      }
      tree1 = &(mTree1.tree);
      tree2 = &(mTree2.tree);

      // ---------------------
      // ----- Compute Distance
      // --------------------
      dataType distance
        = computeDistance<dataType>(tree1, tree2, outputMatching);

      // ---------------------
      // ----- Postprocessing
      // --------------------
      if(postprocess_) {
        postprocessingPipeline<dataType>(tree1);
        postprocessingPipeline<dataType>(tree2);
        if(branchDecomposition_)
          convertBranchDecompositionMatching<dataType>(
            tree1, tree2, outputMatching);
      }

      // std::cout << "TIME COMP.MATC. = " << t_match_time << std::endl;
      printMsg("Total", 1, t_total.getElapsedTime(), this->threadNumber_,
               debug::LineMode::NEW, debug::Priority::INFO);
      printMsg(debug::Separator::L2);
      std::stringstream ss2;
      ss2 << "DISTANCE²       = " << distance;
      printMsg(ss2.str());
      std::stringstream ss3;
      ss3 << "DISTANCE        = " << std::sqrt(distance);
      printMsg(ss3.str());
      printMsg(debug::Separator::L2);
      std::stringstream ss4;
      ss4 << "MEMORY          = " << m.getElapsedUsage();
      printMsg(ss4.str());
      printMsg(debug::Separator::L2);

      if(saveTree_) {
        mTree1 = tree1Ori;
        mTree2 = tree2Ori;
      }

      return distance;
    }

    template <class dataType>
    dataType execute(
      MergeTree<dataType> &tree1,
      MergeTree<dataType> &tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        realOutputMatching;
      dataType res = execute<dataType>(tree1, tree2, realOutputMatching);
      for(auto tup : realOutputMatching)
        outputMatching.push_back(
          std::make_tuple(std::get<0>(tup), std::get<1>(tup)));
      return res;
    }

    template <class dataType>
    void computeEditDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      int nRows,
      int nCols) {
      Timer t_dyn;
      t_assignment_time_ = 0;

      if(parallelize_) {
        parallelEditDistance(tree1, tree2, treeTable, forestTable,
                             treeBackTable, forestBackTable, nRows, nCols);
      } else {
        // Distance T1 to empty tree
        classicEditDistance(tree1, tree2, true, true, tree1->getRoot(),
                            tree2->getRoot(), treeTable, forestTable,
                            treeBackTable, forestBackTable, nRows, nCols);
        if(onlyEmptyTreeDistance_)
          return;
        // Distance T2 to empty tree
        classicEditDistance(tree1, tree2, false, true, tree1->getRoot(),
                            tree2->getRoot(), treeTable, forestTable,
                            treeBackTable, forestBackTable, nRows, nCols);
        // Distance T1 to T2
        classicEditDistance(tree1, tree2, true, false, tree1->getRoot(),
                            tree2->getRoot(), treeTable, forestTable,
                            treeBackTable, forestBackTable, nRows, nCols);
      }

      printMsg("Dynamic programing", 1, t_dyn.getElapsedTime(),
               this->threadNumber_, debug::LineMode::NEW,
               debug::Priority::INFO);
      if(not parallelize_)
        printMsg("Assignment problems", 1, t_assignment_time_,
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::INFO);
    }

    template <class dataType>
    void classicEditDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool processTree1,
      bool computeEmptyTree,
      ftm::idNode nodeI,
      ftm::idNode nodeJ,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      int nRows,
      int nCols) {
      if(processTree1) {
        std::vector<ftm::idNode> childrens;
        tree1->getChildren(nodeI, childrens);
        for(auto children : childrens)
          classicEditDistance(tree1, tree2, processTree1, computeEmptyTree,
                              children, nodeJ, treeTable, forestTable,
                              treeBackTable, forestBackTable, nRows, nCols);
      } else {
        std::vector<ftm::idNode> childrens;
        tree2->getChildren(nodeJ, childrens);
        for(auto children : childrens)
          classicEditDistance(tree1, tree2, processTree1, computeEmptyTree,
                              nodeI, children, treeTable, forestTable,
                              treeBackTable, forestBackTable, nRows, nCols);
      }

      if(processTree1) {
        if(computeEmptyTree) {
          int i = nodeI + 1;
          // --- Equation 8
          computeEquation8(tree1, nodeI, i, treeTable, forestTable);

          // --- Equation 9
          computeEquation9(tree1, nodeI, i, treeTable, forestTable);
        } else
          classicEditDistance(tree1, tree2, false, false, nodeI,
                              tree2->getRoot(), treeTable, forestTable,
                              treeBackTable, forestBackTable, nRows, nCols);
      } else {
        int j = nodeJ + 1;
        if(computeEmptyTree) {
          // --- Equation 10
          computeEquation10(tree2, nodeJ, j, treeTable, forestTable);

          // --- Equation 11
          computeEquation11(tree2, nodeJ, j, treeTable, forestTable);
          //}else{
        } else if(keepSubtree_ or tree1Level_[nodeI] == tree2Level_[nodeJ]) {
          int i = nodeI + 1;
          std::vector<ftm::idNode> children1;
          tree1->getChildren(nodeI, children1);
          std::vector<ftm::idNode> children2;
          tree2->getChildren(nodeJ, children2);
          // --- Equation 13
          computeEquation13(tree1, tree2, i, j, treeTable, forestTable,
                            forestBackTable, children1, children2);

          // --- Equation 12
          computeEquation12(tree1, tree2, i, j, nodeI, nodeJ, treeTable,
                            forestTable, treeBackTable, children1, children2);
        }
      }
    }

    // ----------------------------------------
    // Parallel version
    // ----------------------------------------
    template <class dataType>
    void parallelEditDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      int nRows,
      int nCols) {
      std::vector<int> tree1NodeChildSize, tree2NodeChildSize;
      for(unsigned int i = 0; i < tree1->getNumberOfNodes(); ++i) {
        std::vector<idNode> children;
        tree1->getChildren(i, children);
        tree1NodeChildSize.push_back(children.size());
      }
      for(unsigned int j = 0; j < tree2->getNumberOfNodes(); ++j) {
        std::vector<idNode> children;
        tree2->getChildren(j, children);
        tree2NodeChildSize.push_back(children.size());
      }

      // Get trees data
      std::vector<ftm::idNode> tree1Leaves;
      tree1->getLeavesFromTree(tree1Leaves);
      std::vector<ftm::idNode> tree2Leaves;
      tree2->getLeavesFromTree(tree2Leaves);

      // Distance T1 to empty tree
      parallelEmptyTreeDistance_v2(tree1, true, tree1Leaves, tree1NodeChildSize,
                                   treeTable, forestTable, treeBackTable,
                                   forestBackTable);
      if(onlyEmptyTreeDistance_)
        return;
      // Distance T2 to empty tree
      parallelEmptyTreeDistance_v2(tree2, false, tree2Leaves,
                                   tree2NodeChildSize, treeTable, forestTable,
                                   treeBackTable, forestBackTable);
      // Distance T1 to T2
      parallelTreeDistance_v2(tree1, tree2, true, 0, tree1Leaves,
                              tree1NodeChildSize, tree2Leaves,
                              tree2NodeChildSize, treeTable, forestTable,
                              treeBackTable, forestBackTable, true);
    }

    // Equation 12, 13
    template <class dataType>
    void parallelTreeDistance_v2(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool isTree1,
      int i,
      std::vector<ftm::idNode> &tree1Leaves,
      std::vector<int> &tree1NodeChildSize,
      std::vector<ftm::idNode> &tree2Leaves,
      std::vector<int> &tree2NodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      bool firstCall = false) {
      ftm::idNode nodeT = -1;
      ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
      std::vector<int> treeChildDone(treeT->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(treeT->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;

      if(isTree1)
        for(ftm::idNode leaf : tree1Leaves)
          treeQueue.emplace(leaf);
      else if(keepSubtree_)
        for(ftm::idNode leaf : tree2Leaves)
          treeQueue.emplace(leaf);
      else if(tree1Level_[i - 1] < (int)tree2LevelToNode_.size())
        for(ftm::idNode node : tree2LevelToNode_[tree1Level_[i - 1]])
          treeQueue.emplace(node);

      if(not isCalled_) // and firstCall)
        parallelTreeDistancePara(tree1, tree2, isTree1, i, tree1Leaves,
                                 tree1NodeChildSize, tree2Leaves,
                                 tree2NodeChildSize, treeTable, forestTable,
                                 treeBackTable, forestBackTable, firstCall,
                                 nodeT, treeChildDone, treeNodeDone, treeQueue);
      else
        parallelTreeDistanceTask(tree1, tree2, isTree1, i, tree1Leaves,
                                 tree1NodeChildSize, tree2Leaves,
                                 tree2NodeChildSize, treeTable, forestTable,
                                 treeBackTable, forestBackTable, nodeT,
                                 treeChildDone, treeNodeDone, treeQueue);
    }

    // TODO verify use of first call for the distance computation only
    // (isCalled_=false)
    template <class dataType>
    void parallelTreeDistancePara(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool isTree1,
      int i,
      std::vector<ftm::idNode> &tree1Leaves,
      std::vector<int> &tree1NodeChildSize,
      std::vector<ftm::idNode> &tree2Leaves,
      std::vector<int> &tree2NodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      bool firstCall,
      ftm::idNode nodeT,
      std::vector<int> &treeChildDone,
      std::vector<bool> &treeNodeDone,
      std::queue<ftm::idNode> &treeQueue) {
#ifdef TTK_ENABLE_OPENMP
      unsigned int nthreads = std::thread::hardware_concurrency();
#pragma omp parallel num_threads( \
  this->threadNumber_) if(firstCall or (int) nthreads == this->threadNumber_)
      {
#pragma omp single nowait
#endif
        parallelTreeDistanceTask(tree1, tree2, isTree1, i, tree1Leaves,
                                 tree1NodeChildSize, tree2Leaves,
                                 tree2NodeChildSize, treeTable, forestTable,
                                 treeBackTable, forestBackTable, nodeT,
                                 treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void parallelTreeDistanceTask(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool isTree1,
      int i,
      std::vector<ftm::idNode> &tree1Leaves,
      std::vector<int> &tree1NodeChildSize,
      std::vector<ftm::idNode> &tree2Leaves,
      std::vector<int> &tree2NodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      ftm::idNode nodeT,
      std::vector<int> &treeChildDone,
      std::vector<bool> &treeNodeDone,
      std::queue<ftm::idNode> &treeQueue) {
      int nodePerTask = nodePerTask_;
      while(!treeQueue.empty()) {
        std::queue<ftm::idNode> taskQueue;
        nodePerTask = nodePerTask > (int)treeQueue.size() ? treeQueue.size()
                                                          : nodePerTask;
        for(int j = 0; j < nodePerTask; ++j) {
          nodeT = treeQueue.front();
          treeQueue.pop();
          taskQueue.emplace(nodeT);
        }
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(taskQueue, nodeT)                         \
  untied shared(treeTable, forestTable, treeBackTable, forestBackTable, \
                treeChildDone, treeNodeDone) if(isTree1)
        {
#endif
          ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
          // while(nodeT != -1){
          while(!taskQueue.empty()) {
            nodeT = taskQueue.front();
            taskQueue.pop();
            int t = nodeT + 1;
            ftm::idNode nodeI = i - 1;

            if(isTree1) {
              parallelTreeDistance_v2(
                tree1, tree2, false, t, tree1Leaves, tree1NodeChildSize,
                tree2Leaves, tree2NodeChildSize, treeTable, forestTable,
                treeBackTable, forestBackTable, false);
              //}else{
            } else if(keepSubtree_
                      or tree1Level_[nodeI] == tree2Level_[nodeT]) {
              int j = nodeT + 1;
              std::vector<ftm::idNode> children1;
              tree1->getChildren(nodeI, children1);
              std::vector<ftm::idNode> children2;
              tree2->getChildren(nodeT, children2);
              // --- Equation 13
              computeEquation13(tree1, tree2, i, j, treeTable, forestTable,
                                forestBackTable, children1, children2);

              // --- Equation 12
              computeEquation12(tree1, tree2, i, j, nodeI, nodeT, treeTable,
                                forestTable, treeBackTable, children1,
                                children2);
            }

            if(not isTree1 and not keepSubtree_
               and tree1Level_[nodeI] == tree2Level_[nodeT])
              continue;

            // Manage parent
            ftm::idNode nodeTParent = treeT->getParentSafe(nodeT);
            int childSize = (isTree1) ? tree1NodeChildSize[nodeTParent]
                                      : tree2NodeChildSize[nodeTParent];
            int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
            {
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
            } // pragma omp atomic capture
#endif
            if(not treeNodeDone[nodeTParent]
               and oldTreeChildDone + 1 == childSize) {
              // nodeT = nodeTParent;
              taskQueue.emplace(nodeTParent);
              treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
            } else
              nodeT = -1;

          } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
        } // pragma omp task
#endif
      } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
    }

    // Equation 8, 9, 10, 11
    template <class dataType>
    void parallelEmptyTreeDistance_v2(
      ftm::FTMTree_MT *tree,
      bool isTree1,
      std::vector<ftm::idNode> &treeLeaves,
      std::vector<int> &treeNodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable) {
      ftm::idNode nodeT = -1;
      std::vector<int> treeChildDone(tree->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(tree->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;
      for(ftm::idNode leaf : treeLeaves)
        treeQueue.emplace(leaf);
      if(not isCalled_)
        parallelEmptyTreeDistancePara(tree, isTree1, treeLeaves,
                                      treeNodeChildSize, treeTable, forestTable,
                                      treeBackTable, forestBackTable, nodeT,
                                      treeChildDone, treeNodeDone, treeQueue);
      else
        parallelEmptyTreeDistanceTask(tree, isTree1, treeLeaves,
                                      treeNodeChildSize, treeTable, forestTable,
                                      treeBackTable, forestBackTable, nodeT,
                                      treeChildDone, treeNodeDone, treeQueue);
    }

    template <class dataType>
    void parallelEmptyTreeDistancePara(
      ftm::FTMTree_MT *tree,
      bool isTree1,
      std::vector<ftm::idNode> &treeLeaves,
      std::vector<int> &treeNodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      ftm::idNode nodeT,
      std::vector<int> &treeChildDone,
      std::vector<bool> &treeNodeDone,
      std::queue<ftm::idNode> &treeQueue) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp single nowait
#endif
        parallelEmptyTreeDistanceTask(tree, isTree1, treeLeaves,
                                      treeNodeChildSize, treeTable, forestTable,
                                      treeBackTable, forestBackTable, nodeT,
                                      treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void parallelEmptyTreeDistanceTask(
      ftm::FTMTree_MT *tree,
      bool isTree1,
      std::vector<ftm::idNode> &treeLeaves,
      std::vector<int> &treeNodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      ftm::idNode nodeT,
      std::vector<int> &treeChildDone,
      std::vector<bool> &treeNodeDone,
      std::queue<ftm::idNode> &treeQueue) {
      while(!treeQueue.empty()) {
        nodeT = treeQueue.front();
        treeQueue.pop();

#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT)                                    \
  untied shared(treeTable, forestTable, treeBackTable, forestBackTable, \
                treeChildDone, treeNodeDone)
        {
#endif
          while((int)nodeT != -1) {
            if(isTree1) {
              int i = nodeT + 1;
              // --- Equation 8
              computeEquation8(tree, nodeT, i, treeTable, forestTable);

              // --- Equation 9
              computeEquation9(tree, nodeT, i, treeTable, forestTable);
            } else {
              int j = nodeT + 1;
              // --- Equation 10
              computeEquation10(tree, nodeT, j, treeTable, forestTable);

              // --- Equation 11
              computeEquation11(tree, nodeT, j, treeTable, forestTable);
            }

            // Manage parent
            ftm::idNode nodeTParent = tree->getParentSafe(nodeT);
            int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
            {
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
            } // pragma omp atomic capture
#endif
            if(not treeNodeDone[nodeTParent]
               and oldTreeChildDone + 1 == treeNodeChildSize[nodeTParent]) {
              nodeT = nodeTParent;
              treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
            } else
              nodeT = -1;

          } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
        } // pragma omp task
#endif
      } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
    }

    // ----------------------------------------
    // OLD v1
    // ----------------------------------------
    // Equation 12, 13
    template <class dataType>
    void parallelTreeDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      bool isTree1,
      int i,
      std::vector<ftm::idNode> &tree1Leaves,
      std::vector<int> &tree1NodeChildSize,
      std::vector<ftm::idNode> &tree2Leaves,
      std::vector<int> &tree2NodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      bool firstCall = false) {
      ftm::idNode nodeT = -1;
      ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
      std::vector<int> treeChildDone(treeT->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(treeT->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;
      if(isTree1)
        for(ftm::idNode leaf : tree1Leaves)
          treeQueue.emplace(leaf);
      else
        for(ftm::idNode leaf : tree2Leaves)
          treeQueue.emplace(leaf);

#ifdef TTK_ENABLE_OPENMP
      const int nthreads = std::thread::hardware_concurrency();
#pragma omp parallel num_threads(                                        \
  this->threadNumber_) if((firstCall or nthreads == this->threadNumber_) \
                          && !isCalled_)
      {
#pragma omp single nowait
#endif
        while(!treeQueue.empty()) {
          nodeT = treeQueue.front();
          treeQueue.pop();
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) \
  untied // shared(treeTable, forestTable)//, treeBackTable, forestBackTable)
          {
#endif
            std::stringstream ss;
            ss << &treeTable[0] << " _ " << tree1 << " _ " << tree2
               << std::endl;
            // std::cout << ss.str();
            while(static_cast<int>(nodeT) != -1) {
              int t = nodeT + 1;

              if(isTree1) {
                parallelTreeDistance(tree1, tree2, false, t, tree1Leaves,
                                     tree1NodeChildSize, tree2Leaves,
                                     tree2NodeChildSize, treeTable, forestTable,
                                     treeBackTable, forestBackTable);
              } else {
                int j = nodeT + 1;
                ftm::idNode nodeI = i - 1;
                std::vector<ftm::idNode> children1;
                tree1->getChildren(nodeI, children1);
                std::vector<ftm::idNode> children2;
                tree2->getChildren(nodeT, children2);
                // --- Equation 13
                computeEquation13(tree1, tree2, i, j, treeTable, forestTable,
                                  forestBackTable, children1, children2);

                // --- Equation 12
                computeEquation12(tree1, tree2, i, j, nodeI, nodeT, treeTable,
                                  forestTable, treeBackTable, children1,
                                  children2);
              }

              // Manage parent
              ftm::idNode nodeTParent = treeT->getParentSafe(nodeT);
              int childSize = (isTree1) ? tree1NodeChildSize[nodeTParent]
                                        : tree2NodeChildSize[nodeTParent];
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
              {
#endif
                oldTreeChildDone = treeChildDone[nodeTParent];
                treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
              } // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent]
                 and oldTreeChildDone + 1 == childSize) {
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
              } else
                nodeT = -1;

            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
          } // pragma omp task
#endif
        } // while loop
#ifdef TTK_ENABLE_OPENMP
      } // pragma omp parallel
#endif
    }

    // Equation 8, 9, 10, 11
    template <class dataType>
    void parallelEmptyTreeDistance(
      ftm::FTMTree_MT *tree,
      bool isTree1,
      std::vector<ftm::idNode> &treeLeaves,
      std::vector<int> &treeNodeChildSize,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable) {
      ftm::idNode nodeT = -1;
      std::vector<int> treeChildDone(tree->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(tree->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;
      for(ftm::idNode leaf : treeLeaves) {
        treeQueue.emplace(leaf);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_) if(not isCalled_)
      {
#pragma omp single nowait
#endif
        while(!treeQueue.empty()) {
          nodeT = treeQueue.front();
          treeQueue.pop();

#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(nodeT) \
  untied // shared(treeTable, forestTable, treeBackTable, forestBackTable)
          {
#endif
            while(static_cast<int>(nodeT) != -1) {
              if(isTree1) {
                int i = nodeT + 1;
                // --- Equation 8
                computeEquation8(tree, nodeT, i, treeTable, forestTable);

                // --- Equation 9
                computeEquation9(tree, nodeT, i, treeTable, forestTable);
              } else {
                int j = nodeT + 1;
                // --- Equation 10
                computeEquation10(tree, nodeT, j, treeTable, forestTable);

                // --- Equation 11
                computeEquation11(tree, nodeT, j, treeTable, forestTable);
              }

              // Manage parent
              ftm::idNode nodeTParent = tree->getParentSafe(nodeT);
              int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
              {
#endif
                oldTreeChildDone = treeChildDone[nodeTParent];
                treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP
              } // pragma omp atomic capture
#endif
              if(not treeNodeDone[nodeTParent]
                 and oldTreeChildDone + 1 == treeNodeChildSize[nodeTParent]) {
                nodeT = nodeTParent;
                treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
              } else
                nodeT = -1;

            } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP
          } // pragma omp task
#endif
        } // while loop
#ifdef TTK_ENABLE_OPENMP
      } // pragma omp parallel
#endif
    }

    // ----------------------------------------
    // Utils
    // ----------------------------------------
    void printMapIntInt(std::map<int, int> theMap) {
      for(auto itr = theMap.begin(); itr != theMap.end(); ++itr) {
        std::stringstream ss;
        ss << '\t' << itr->first << '\t' << itr->second;
        printMsg(ss.str());
      }
      printMsg("");
    }

    template <class dataType>
    void verifyMergeTreeStructure(FTMTree_MT *tree) {
      bool problem = false;

      bool isJT = tree->isJoinTree<dataType>();
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> problemNodes;
      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        if(!tree->isRoot(node)) {
          bool thisProblem;
          if(isJT)
            thisProblem = tree->getValue<dataType>(node)
                          > tree->getValue<dataType>(tree->getParentSafe(node));
          else
            thisProblem = tree->getValue<dataType>(node)
                          < tree->getValue<dataType>(tree->getParentSafe(node));

          if(thisProblem)
            problemNodes.push_back(
              std::make_tuple(node, tree->getParentSafe(node)));

          problem |= thisProblem;
        }

        std::vector<idNode> children;
        tree->getChildren(node, children);
        for(auto c : children)
          queue.emplace(c);
      }

      if(problem) {
        printMsg("merge tree in input is not valid");
        for(auto tup : problemNodes) {
          std::stringstream ss;
          ss << std::get<0>(tup) << " _ " << std::get<1>(tup);
          printMsg(ss.str());
        }
        tree->printTree();
        tree->printTreeScalars<dataType>();
      }
    }

    // ----------------------------------------
    // Testing
    // ----------------------------------------
    template <class dataType>
    void classicalPersistenceAssignmentProblem(ftm::FTMTree_MT *tree1,
                                               ftm::FTMTree_MT *tree2) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs1,
        pairs2;
      tree1->getPersistencePairsFromTree(pairs1);
      tree2->getPersistencePairsFromTree(pairs2);
      std::vector<std::vector<dataType>> costMatrix(
        pairs1.size() + 1, std::vector<dataType>(pairs2.size() + 1));
      std::stringstream ss;
      ss << costMatrix.size() << " _ " << costMatrix[0].size();
      printMsg(ss.str());
      for(unsigned int i = 0; i < costMatrix.size() - 1; ++i) {
        dataType nodeIValue = tree1->getValue<dataType>(std::get<0>(pairs1[i]));
        dataType nodeIOriginValue
          = tree1->getValue<dataType>(std::get<1>(pairs1[i]));
        for(unsigned int j = 0; j < costMatrix[0].size() - 1; ++j) {
          dataType nodeJValue
            = tree2->getValue<dataType>(std::get<0>(pairs2[j]));
          dataType nodeJOriginValue
            = tree2->getValue<dataType>(std::get<1>(pairs2[j]));
          costMatrix[i][j] = std::pow(nodeIValue - nodeJValue, 2)
                             + std::pow(nodeIOriginValue - nodeJOriginValue, 2);
        }
        costMatrix[i][costMatrix[0].size() - 1]
          = 2 * std::pow(std::get<2>(pairs1[i]), 2) / (std::pow(2, 2));
      }
      for(unsigned int j = 0; j < costMatrix[0].size() - 1; ++j)
        costMatrix[costMatrix.size() - 1][j]
          = 2 * std::pow(std::get<2>(pairs2[j]), 2) / (std::pow(2, 2));
      std::vector<asgnMatchingTuple> matchings;
      forestAssignmentProblemMunkres(costMatrix, matchings);
      dataType cost = 0;
      for(auto tuple : matchings)
        cost += std::get<2>(tuple);
      std::stringstream ss2;
      ss2 << "cost      = " << cost;
      printMsg(ss2.str());
      std::stringstream ss3;
      ss3 << "cost sqrt = " << std::sqrt(cost);
      printMsg(ss3.str());
    }

  }; // MergeTreeDistance class

} // namespace ttk

#endif
