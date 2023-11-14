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

    double minMaxPairWeight_ = 1.0;

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
#ifdef TTK_ENABLE_OPENMP4
      omp_set_nested(1);
#endif
    }
    ~MergeTreeDistance() override = default;

    void setIsCalled(bool ic) {
      isCalled_ = ic;
    }

    void setPreprocess(bool preproc) {
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

    void setMinMaxPairWeight(double weight) {
      minMaxPairWeight_ = weight;
    }

    /**
     * Implementation of the algorithm.
     */

    // ------------------------------------------------------------------------
    // Assignment Problem
    // ------------------------------------------------------------------------
    template <class dataType>
    void
      runAssignmentProblemSolver(std::vector<std::vector<dataType>> &costMatrix,
                                 std::vector<MatchingType> &matchings) {
      AssignmentSolver<dataType> *assignmentSolver;
      AssignmentExhaustive<dataType> solverExhaustive;
      AssignmentMunkres<dataType> solverMunkres;
      AssignmentAuction<dataType> solverAuction;

      int const nRows = costMatrix.size() - 1;
      int const nCols = costMatrix[0].size() - 1;
      int const max_dim = std::max(nRows, nCols);
      int const min_dim = std::min(nRows, nCols);

      int assignmentSolverID = assignmentSolverID_;
      if((min_dim <= 2 and max_dim <= 2) or (min_dim <= 1 and max_dim <= 6))
        assignmentSolverID = 1;

      switch(assignmentSolverID) {
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
        int const forestTableI = children1[i] + 1;
        for(unsigned int j = 0; j < nCols; ++j) {
          int const forestTableJ = children2[j] + 1;
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
        int const forestTableJ = children2[j] + 1;
        // Cost of not assigning j
        costMatrix[nRows][j] = treeTable[0][forestTableJ];
      }
      costMatrix[nRows][nCols] = 0;
    }

    template <class dataType>
    dataType postprocessAssignment(
      std::vector<MatchingType> &matchings,
      std::vector<ftm::idNode> &children1,
      std::vector<ftm::idNode> &children2,
      std::vector<std::tuple<int, int>> &forestAssignment) {
      dataType cost = 0;
      for(const auto &mTuple : matchings) {
        cost += std::get<2>(mTuple);
        if(std::get<0>(mTuple) >= (int)children1.size()
           || std::get<1>(mTuple) >= (int)children2.size())
          continue;
        int const tableId1 = children1[std::get<0>(mTuple)] + 1;
        int const tableId2 = children2[std::get<1>(mTuple)] + 1;
        forestAssignment.emplace_back(tableId1, tableId2);
      }
      return cost;
    }

    template <class dataType>
    dataType forestAssignmentProblem(
      ftm::FTMTree_MT *ttkNotUsed(tree1),
      ftm::FTMTree_MT *ttkNotUsed(tree2),
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
      std::vector<MatchingType> matchings;
      runAssignmentProblemSolver(costMatrix, matchings);

      // --- Postprocess matching to create output assignment
      dataType cost = postprocessAssignment<dataType>(
        matchings, children1, children2, forestAssignment);

      return cost;
    }

    template <class dataType>
    void computeForestsDistance(
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
        dataType forestTerm3;

        // Term 3
        Timer t_assignment;
        std::vector<std::tuple<int, int>> forestAssignment;
        forestTerm3 = forestAssignmentProblem<dataType>(
          tree1, tree2, treeTable, children1, children2, forestAssignment);
        if(not parallelize_)
          t_assignment_time_ += t_assignment.getElapsedTime();

        if(not keepSubtree_) {
          // Compute table value
          forestTable[i][j] = forestTerm3;
          // Add backtracking information
          forestBackTable[i][j] = forestAssignment;
        } else {
          dataType forestTerm1, forestTerm2;
          std::tuple<dataType, ftm::idNode> forestCoTerm1, forestCoTerm2;

          // Term 1
          forestCoTerm1
            = computeTerm1_2<dataType>(children2, i, forestTable, true);
          forestTerm1 = forestTable[0][j] + std::get<0>(forestCoTerm1);

          // Term2
          forestCoTerm2
            = computeTerm1_2<dataType>(children1, j, forestTable, false);
          forestTerm2 = forestTable[i][0] + std::get<0>(forestCoTerm2);

          // Compute table value
          forestTable[i][j]
            = std::min(std::min(forestTerm1, forestTerm2), forestTerm3);

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
        }
      } else {
        // If one of the forest is empty we get back to equation 8 or 10
        forestTable[i][j]
          = (children1.size() == 0) ? forestTable[0][j] : forestTable[i][0];
      }
    }

    // ------------------------------------------------------------------------
    // Edit Distance Dynamic Programming Equations
    // ------------------------------------------------------------------------
    template <class dataType>
    void computeForestToEmptyDistance(
      ftm::FTMTree_MT *tree1,
      ftm::idNode nodeI,
      int i,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable) {
      std::vector<ftm::idNode> children;
      tree1->getChildren(nodeI, children);
      forestTable[i][0] = 0;
      for(ftm::idNode const child : children)
        forestTable[i][0] += treeTable[child + 1][0];
    }

    template <class dataType>
    void computeSubtreeToEmptyDistance(
      ftm::FTMTree_MT *tree1,
      ftm::idNode nodeI,
      int i,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable) {
      treeTable[i][0] = forestTable[i][0] + deleteCost<dataType>(tree1, nodeI);
    }

    template <class dataType>
    void computeEmptyToForestDistance(
      ftm::FTMTree_MT *tree2,
      ftm::idNode nodeJ,
      int j,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable) {
      std::vector<ftm::idNode> children;
      tree2->getChildren(nodeJ, children);
      forestTable[0][j] = 0;
      for(ftm::idNode const child : children)
        forestTable[0][j] += treeTable[0][child + 1];
    }

    template <class dataType>
    void computeEmptyToSubtreeDistance(
      ftm::FTMTree_MT *tree2,
      ftm::idNode nodeJ,
      int j,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable) {
      treeTable[0][j] = forestTable[0][j] + insertCost<dataType>(tree2, nodeJ);
    }

    // Compute first or second term of forests and subtrees distance
    template <class dataType>
    std::tuple<dataType, ftm::idNode>
      computeTerm1_2(std::vector<ftm::idNode> &childrens,
                     int ind,
                     std::vector<std::vector<dataType>> &table,
                     bool computeTerm1) {
      dataType tempMin = (childrens.size() == 0)
                           ? ((computeTerm1) ? table[ind][0] : table[0][ind])
                           : std::numeric_limits<dataType>::max();
      ftm::idNode bestIdNode = 0;
      for(ftm::idNode children : childrens) {
        children += 1;
        dataType temp;
        if(computeTerm1) {
          temp = table[ind][children] - table[0][children];
        } else {
          temp = table[children][ind] - table[children][0];
        }
        if(temp < tempMin) {
          tempMin = temp;
          bestIdNode = children;
        }
      }
      return std::make_tuple(tempMin, bestIdNode);
    }

    template <class dataType>
    void computeSubtreesDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      int i,
      int j,
      ftm::idNode nodeI,
      ftm::idNode nodeJ,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<ftm::idNode> &children1,
      std::vector<ftm::idNode> &children2) {
      dataType treeTerm3;

      // Term 3
      treeTerm3
        = forestTable[i][j] + relabelCost<dataType>(tree1, nodeI, tree2, nodeJ);

      if(not keepSubtree_) {
        // Compute table value
        treeTable[i][j] = treeTerm3;
        // Add backtracking information
        treeBackTable[i][j] = std::make_tuple(i, j);
      } else {
        dataType treeTerm1, treeTerm2;
        std::tuple<dataType, ftm::idNode> treeCoTerm1, treeCoTerm2;

        // Term 1
        treeCoTerm1 = computeTerm1_2<dataType>(children2, i, treeTable, true);
        treeTerm1 = treeTable[0][j] + std::get<0>(treeCoTerm1);

        // Term 2
        treeCoTerm2 = computeTerm1_2<dataType>(children1, j, treeTable, false);
        treeTerm2 = treeTable[i][0] + std::get<0>(treeCoTerm2);

        // Compute table value
        treeTable[i][j] = std::min(std::min(treeTerm1, treeTerm2), treeTerm3);

        // Add backtracking information
        if(treeTable[i][j] == treeTerm3) {
          treeBackTable[i][j] = std::make_tuple(i, j);
        } else if(treeTable[i][j] == treeTerm2) {
          treeBackTable[i][j] = std::make_tuple(std::get<1>(treeCoTerm2), j);
        } else {
          treeBackTable[i][j] = std::make_tuple(i, std::get<1>(treeCoTerm1));
        }
      }
    }

    // --------------------------------------------------------------------------------
    // Output Matching
    // --------------------------------------------------------------------------------
    template <class dataType>
    void computeMatching(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &outputMatching,
      int startR,
      int startC) {
      outputMatching.clear();
      std::queue<std::tuple<int, int, bool>> backQueue;
      backQueue.emplace(startR, startC, true);
      while(!backQueue.empty()) {
        std::tuple<int, int, bool> elem = backQueue.front();
        backQueue.pop();
        bool useTreeTable = std::get<2>(elem);
        int const i = std::get<0>(elem);
        int const j = std::get<1>(elem);

        if(useTreeTable) {
          int const tupleI = std::get<0>(treeBackTable[i][j]);
          int const tupleJ = std::get<1>(treeBackTable[i][j]);
          if(tupleI != 0 && tupleJ != 0) {
            useTreeTable = (tupleI != i || tupleJ != j);
            backQueue.emplace(tupleI, tupleJ, useTreeTable);
            if(not useTreeTable) { // We have matched i and j
              ftm::idNode const tree1Node = tupleI - 1;
              ftm::idNode const tree2Node = tupleJ - 1;
              double cost = 0;
              dataType costT
                = relabelCost<dataType>(tree1, tree1Node, tree2, tree2Node);
              cost = static_cast<double>(costT);
              outputMatching.emplace_back(tree1Node, tree2Node, cost);
            }
          }
        } else {
          for(std::tuple<int, int> forestBackElem : forestBackTable[i][j]) {
            int const tupleI = std::get<0>(forestBackElem);
            int const tupleJ = std::get<1>(forestBackElem);
            if(tupleI != 0 && tupleJ != 0) {
              useTreeTable = (tupleI != i && tupleJ != j);
              backQueue.emplace(tupleI, tupleJ, useTreeTable);
            }
          }
        }
      }
    }

    // ------------------------------------------------------------------------
    // Main Functions
    // ------------------------------------------------------------------------
    template <class dataType>
    dataType
      computeDistance(ftm::FTMTree_MT *tree1,
                      ftm::FTMTree_MT *tree2,
                      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                        &outputMatching) {
      // ---------------------
      // ----- Init dynamic programming tables
      // --------------------
      size_t const nRows = tree1->getNumberOfNodes() + 1;
      size_t const nCols = tree2->getNumberOfNodes() + 1;
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

      int const indR = tree1->getRoot() + 1;
      int const indC = tree2->getRoot() + 1;

      tree1->getAllNodeLevel(tree1Level_);
      tree2->getAllNodeLevel(tree2Level_);
      tree2->getLevelToNode(tree2LevelToNode_);

      // ---------------------
      // ----- Compute edit distance
      // --------------------
      computeEditDistance(tree1, tree2, treeTable, forestTable, treeBackTable,
                          forestBackTable, nRows, nCols);
      dataType distance = treeTable[indR][indC];
      if(onlyEmptyTreeDistance_)
        distance = treeTable[indR][0];
      if(branchDecomposition_) {
        if(not useMinMaxPair_) {
          if(onlyEmptyTreeDistance_)
            distance -= deleteCost<dataType>(tree1, tree1->getRoot());
          else
            distance -= relabelCost<dataType>(
              tree1, tree1->getRoot(), tree2, tree2->getRoot());
        } else {
          if(minMaxPairWeight_ != 1.0) {
            auto cost = relabelCost<dataType>(
              tree1, tree1->getRoot(), tree2, tree2->getRoot());
            distance = distance - cost + minMaxPairWeight_ * cost;
          }
        }
      }
      if(distanceSquaredRoot_)
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
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        realOutputMatching;
      dataType res
        = computeDistance<dataType>(tree1, tree2, realOutputMatching);
      for(auto tup : realOutputMatching)
        outputMatching.emplace_back(std::get<0>(tup), std::get<1>(tup));
      return res;
    }

    template <class dataType>
    dataType execute(ftm::MergeTree<dataType> &mTree1,
                     ftm::MergeTree<dataType> &mTree2,
                     std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
                       &outputMatching) {
      Memory m;
      Timer t_total;

      // ---------------------
      // ----- Testing
      // --------------------
      testing_ = false;

      // ---------------------
      // ----- Preprocessing
      // --------------------
      ftm::MergeTree<dataType> mTree1Copy;
      ftm::MergeTree<dataType> mTree2Copy;
      if(saveTree_) {
        mTree1Copy = ftm::copyMergeTree<dataType>(mTree1);
        mTree2Copy = ftm::copyMergeTree<dataType>(mTree2);
      }
      ftm::MergeTree<dataType> &mTree1Int = (saveTree_ ? mTree1Copy : mTree1);
      ftm::MergeTree<dataType> &mTree2Int = (saveTree_ ? mTree2Copy : mTree2);
      ftm::FTMTree_MT *tree1 = &(mTree1Int.tree);
      ftm::FTMTree_MT *tree2 = &(mTree2Int.tree);
      if(not isCalled_) {
        verifyMergeTreeStructure<dataType>(tree1);
        verifyMergeTreeStructure<dataType>(tree2);
      }
      if(preprocess_) {
        treesNodeCorr_.resize(2);
        preprocessingPipeline<dataType>(
          mTree1Int, epsilonTree1_, epsilon2Tree1_, epsilon3Tree1_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[0]);
        preprocessingPipeline<dataType>(
          mTree2Int, epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          branchDecomposition_, useMinMaxPair_, cleanTree_, treesNodeCorr_[1]);
      }
      tree1 = &(mTree1Int.tree);
      tree2 = &(mTree2Int.tree);

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

      // std::cout << "TIME COMP.MATCH. = " << t_match_time << std::endl;
      printMsg("Total", 1, t_total.getElapsedTime(), this->threadNumber_,
               debug::LineMode::NEW, debug::Priority::INFO);
      printMsg(debug::Separator::L2);
      std::stringstream ss2;
      ss2 << "DISTANCE²       = "
          << (distanceSquaredRoot_ ? distance * distance : distance);
      printMsg(ss2.str());
      std::stringstream ss3;
      ss3 << "DISTANCE        = "
          << (distanceSquaredRoot_ ? distance : std::sqrt(distance));
      printMsg(ss3.str());
      printMsg(debug::Separator::L2);
      std::stringstream ss4;
      ss4 << "MEMORY          = " << m.getElapsedUsage();
      printMsg(ss4.str());
      printMsg(debug::Separator::L2);

      return distance;
    }

    template <class dataType>
    dataType execute(
      ftm::MergeTree<dataType> &tree1,
      ftm::MergeTree<dataType> &tree2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> &outputMatching) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>
        realOutputMatching;
      dataType res = execute<dataType>(tree1, tree2, realOutputMatching);
      for(auto tup : realOutputMatching)
        outputMatching.emplace_back(std::get<0>(tup), std::get<1>(tup));
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
          int const i = nodeI + 1;
          // --- Forest to empty tree distance
          computeForestToEmptyDistance(tree1, nodeI, i, treeTable, forestTable);

          // --- Subtree to empty tree distance
          computeSubtreeToEmptyDistance(
            tree1, nodeI, i, treeTable, forestTable);
        } else
          classicEditDistance(tree1, tree2, false, false, nodeI,
                              tree2->getRoot(), treeTable, forestTable,
                              treeBackTable, forestBackTable, nRows, nCols);
      } else {
        int const j = nodeJ + 1;
        if(computeEmptyTree) {
          // --- Empty tree to forest distance
          computeEmptyToForestDistance(tree2, nodeJ, j, treeTable, forestTable);

          // --- Empty tree to subtree distance
          computeEmptyToSubtreeDistance(
            tree2, nodeJ, j, treeTable, forestTable);
          //}else{
        } else if(keepSubtree_ or tree1Level_[nodeI] == tree2Level_[nodeJ]) {
          int const i = nodeI + 1;
          std::vector<ftm::idNode> children1;
          tree1->getChildren(nodeI, children1);
          std::vector<ftm::idNode> children2;
          tree2->getChildren(nodeJ, children2);
          // --- Forests distance
          computeForestsDistance(tree1, tree2, i, j, treeTable, forestTable,
                                 forestBackTable, children1, children2);

          // --- Subtrees distance
          computeSubtreesDistance(tree1, tree2, i, j, nodeI, nodeJ, treeTable,
                                  forestTable, treeBackTable, children1,
                                  children2);
        }
      }
    }

    // ------------------------------------------------------------------------
    // Parallel version
    // ------------------------------------------------------------------------
    template <class dataType>
    void parallelEditDistance(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      std::vector<std::vector<dataType>> &treeTable,
      std::vector<std::vector<dataType>> &forestTable,
      std::vector<std::vector<std::tuple<int, int>>> &treeBackTable,
      std::vector<std::vector<std::vector<std::tuple<int, int>>>>
        &forestBackTable,
      int ttkNotUsed(nRows),
      int ttkNotUsed(nCols)) {
      std::vector<int> tree1NodeChildSize, tree2NodeChildSize;
      for(unsigned int i = 0; i < tree1->getNumberOfNodes(); ++i) {
        std::vector<ftm::idNode> children;
        tree1->getChildren(i, children);
        tree1NodeChildSize.push_back(children.size());
      }
      for(unsigned int j = 0; j < tree2->getNumberOfNodes(); ++j) {
        std::vector<ftm::idNode> children;
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

    // Forests and subtrees distances
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
      ftm::idNode const nodeT = -1;
      ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
      std::vector<int> treeChildDone(treeT->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(treeT->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;

      if(isTree1)
        for(ftm::idNode const leaf : tree1Leaves)
          treeQueue.emplace(leaf);
      else if(keepSubtree_)
        for(ftm::idNode const leaf : tree2Leaves)
          treeQueue.emplace(leaf);
      else if(tree1Level_[i - 1] < (int)tree2LevelToNode_.size())
        for(ftm::idNode const node : tree2LevelToNode_[tree1Level_[i - 1]])
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
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel num_threads(this->threadNumber_) if(firstCall)
      {
#pragma omp single nowait
#endif
        parallelTreeDistanceTask(tree1, tree2, isTree1, i, tree1Leaves,
                                 tree1NodeChildSize, tree2Leaves,
                                 tree2NodeChildSize, treeTable, forestTable,
                                 treeBackTable, forestBackTable, nodeT,
                                 treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP4
      } // pragma omp parallel
#endif

      TTK_FORCE_USE(firstCall);
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
#ifdef TTK_ENABLE_OPENMP4
#pragma omp task firstprivate(taskQueue, nodeT) UNTIED()         \
  shared(treeTable, forestTable, treeBackTable, forestBackTable, \
         treeChildDone, treeNodeDone) if(isTree1)
        {
#endif
          ftm::FTMTree_MT *treeT = (isTree1) ? tree1 : tree2;
          // while(nodeT != -1){
          while(!taskQueue.empty()) {
            nodeT = taskQueue.front();
            taskQueue.pop();
            int const t = nodeT + 1;
            ftm::idNode const nodeI = i - 1;

            if(isTree1) {
              parallelTreeDistance_v2(
                tree1, tree2, false, t, tree1Leaves, tree1NodeChildSize,
                tree2Leaves, tree2NodeChildSize, treeTable, forestTable,
                treeBackTable, forestBackTable, false);
              //}else{
            } else if(keepSubtree_
                      or tree1Level_[nodeI] == tree2Level_[nodeT]) {
              int const j = nodeT + 1;
              std::vector<ftm::idNode> children1;
              tree1->getChildren(nodeI, children1);
              std::vector<ftm::idNode> children2;
              tree2->getChildren(nodeT, children2);
              // --- Forests distance
              computeForestsDistance(tree1, tree2, i, j, treeTable, forestTable,
                                     forestBackTable, children1, children2);

              // --- Subtrees distance
              computeSubtreesDistance(tree1, tree2, i, j, nodeI, nodeT,
                                      treeTable, forestTable, treeBackTable,
                                      children1, children2);
            }

            if(not isTree1 and not keepSubtree_
               and tree1Level_[nodeI] == tree2Level_[nodeT])
              continue;

            // Manage parent
            ftm::idNode const nodeTParent = treeT->getParentSafe(nodeT);
            int const childSize = (isTree1) ? tree1NodeChildSize[nodeTParent]
                                            : tree2NodeChildSize[nodeTParent];
            int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP4
#pragma omp atomic capture
            {
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP4
            } // pragma omp atomic capture
#endif
            if(not treeNodeDone[nodeTParent]
               and oldTreeChildDone + 1 == childSize) {
              // nodeT = nodeTParent;
              taskQueue.emplace(nodeTParent);
              treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP4
#pragma omp taskyield
#endif
            } else
              nodeT = -1;

          } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP4
        } // pragma omp task
#endif
      } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP4
#pragma omp taskwait
#endif
    }

    // Subtree/Forest with empty tree distances
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
      ftm::idNode const nodeT = -1;
      std::vector<int> treeChildDone(tree->getNumberOfNodes(), 0);
      std::vector<bool> treeNodeDone(tree->getNumberOfNodes(), false);
      std::queue<ftm::idNode> treeQueue;
      for(ftm::idNode const leaf : treeLeaves)
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
#ifdef TTK_ENABLE_OPENMP4
#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp single nowait
#endif
        parallelEmptyTreeDistanceTask(tree, isTree1, treeLeaves,
                                      treeNodeChildSize, treeTable, forestTable,
                                      treeBackTable, forestBackTable, nodeT,
                                      treeChildDone, treeNodeDone, treeQueue);
#ifdef TTK_ENABLE_OPENMP4
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void parallelEmptyTreeDistanceTask(
      ftm::FTMTree_MT *tree,
      bool isTree1,
      std::vector<ftm::idNode> &ttkNotUsed(treeLeaves),
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

#ifdef TTK_ENABLE_OPENMP4
#pragma omp task firstprivate(nodeT) UNTIED()                    \
  shared(treeTable, forestTable, treeBackTable, forestBackTable, \
         treeChildDone, treeNodeDone)
        {
#endif
          while((int)nodeT != -1) {
            if(isTree1) {
              int const i = nodeT + 1;
              // --- Forest to empty tree distance
              computeForestToEmptyDistance(
                tree, nodeT, i, treeTable, forestTable);

              // --- Subtree to empty tree distance
              computeSubtreeToEmptyDistance(
                tree, nodeT, i, treeTable, forestTable);
            } else {
              int const j = nodeT + 1;
              // --- Empty tree to forest distance
              computeEmptyToForestDistance(
                tree, nodeT, j, treeTable, forestTable);

              // --- Empty tree to subtree distance
              computeEmptyToSubtreeDistance(
                tree, nodeT, j, treeTable, forestTable);
            }

            // Manage parent
            ftm::idNode const nodeTParent = tree->getParentSafe(nodeT);
            int oldTreeChildDone;
#ifdef TTK_ENABLE_OPENMP4
#pragma omp atomic capture
            {
#endif
              oldTreeChildDone = treeChildDone[nodeTParent];
              treeChildDone[nodeTParent]++;
#ifdef TTK_ENABLE_OPENMP4
            } // pragma omp atomic capture
#endif
            if(not treeNodeDone[nodeTParent]
               and oldTreeChildDone + 1 == treeNodeChildSize[nodeTParent]) {
              nodeT = nodeTParent;
              treeNodeDone[nodeTParent] = true;
#ifdef TTK_ENABLE_OPENMP4
#pragma omp taskyield
#endif
            } else
              nodeT = -1;

          } // while nodeI loop
#ifdef TTK_ENABLE_OPENMP4
        } // pragma omp task
#endif
      } // while treeQueue loop
#ifdef TTK_ENABLE_OPENMP4
#pragma omp taskwait
#endif

      TTK_FORCE_USE(treeBackTable);
      TTK_FORCE_USE(forestBackTable);
    }

    // ------------------------------------------------------------------------
    // Utils
    // ------------------------------------------------------------------------
    void printMapIntInt(std::map<int, int> theMap) {
      for(auto itr = theMap.begin(); itr != theMap.end(); ++itr) {
        std::stringstream ss;
        ss << '\t' << itr->first << '\t' << itr->second;
        printMsg(ss.str());
      }
      printMsg("");
    }

    template <class dataType>
    void verifyMergeTreeStructure(ftm::FTMTree_MT *tree) {
      bool problem = false;

      bool const isJT = tree->isJoinTree<dataType>();
      std::vector<std::tuple<ftm::idNode, ftm::idNode>> problemNodes;
      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode const node = queue.front();
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
            problemNodes.emplace_back(node, tree->getParentSafe(node));

          problem |= thisProblem;
        }

        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        for(auto c : children)
          queue.emplace(c);
      }

      if(problem) {
        printErr("merge tree in input is not valid");
        for(auto tup : problemNodes) {
          std::stringstream ss;
          ss << std::get<0>(tup) << " _ " << std::get<1>(tup);
          printMsg(ss.str());
        }
        printMsg(tree->printTree().str());
        printMsg(tree->printTreeScalars<dataType>().str());
      }
    }

    // ------------------------------------------------------------------------
    // Testing
    // ------------------------------------------------------------------------
    template <class dataType>
    void classicalPersistenceAssignmentProblem(ftm::FTMTree_MT *tree1,
                                               ftm::FTMTree_MT *tree2) {
      std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> pairs1,
        pairs2;
      tree1->getPersistencePairsFromTree(pairs1);
      tree2->getPersistencePairsFromTree(pairs2);
      std::vector<std::vector<dataType>> costMatrix(
        pairs1.size() + 1, std::vector<dataType>(pairs2.size() + 1));
      std::stringstream const ss;
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
      std::vector<MatchingType> matchings;
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
