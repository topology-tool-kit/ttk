/// \ingroup base
/// \class BranchMappingDistance
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2022.
///
/// This module defines the %BranchMappingDistance class that computes distances
/// between two merge trees.
///
/// \b Related \b publication \n
/// "Branch Decomposition-Independent Edit Distances for Merge Trees." \n
/// Florian Wetzels, Heike Leitte, and Christoph Garth. \n
/// Computer Graphics Forum, 2022.

#pragma once

#include <set>
#include <vector>

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <vector>

// ttk common includes
#include <AssignmentAuction.h>
#include <AssignmentExhaustive.h>
#include <AssignmentMunkres.h>
#include <Debug.h>
#include <FTMTree_MT.h>

namespace ttk {

  class BranchMappingDistance : virtual public Debug {

  private:
    int baseMetric_ = 0;
    int assignmentSolverID_ = 0;
    bool squared_ = false;

    template <class dataType>
    inline dataType editCost_Wasserstein1(int n1,
                                          int p1,
                                          int n2,
                                          int p2,
                                          ftm::FTMTree_MT *tree1,
                                          ftm::FTMTree_MT *tree2) {
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(p2);
        dataType b2 = (b1 + d1) * 0.5;
        dataType d2 = (b1 + d1) * 0.5;
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = db + dd;
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = (b1 + d1) * 0.5;
        dataType d2 = (b1 + d1) * 0.5;
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = db + dd;
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(p2);
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = db + dd;
      }
      return squared_ ? d * d : d;
    }

    template <class dataType>
    inline dataType editCost_Wasserstein2(int n1,
                                          int p1,
                                          int n2,
                                          int p2,
                                          ftm::FTMTree_MT *tree1,
                                          ftm::FTMTree_MT *tree2) {
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(p2);
        dataType b2 = (b1 + d1) * 0.5;
        dataType d2 = (b1 + d1) * 0.5;
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = std::sqrt(db * db + dd * dd);
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = (b1 + d1) * 0.5;
        dataType d2 = (b1 + d1) * 0.5;
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = std::sqrt(db * db + dd * dd);
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(p2);
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dd = d1 > d2 ? d1 - d2 : d2 - d1;
        d = std::sqrt(db * db + dd * dd);
      }
      return squared_ ? d * d : d;
    }

    template <class dataType>
    inline dataType editCost_Persistence(int n1,
                                         int p1,
                                         int n2,
                                         int p2,
                                         ftm::FTMTree_MT *tree1,
                                         ftm::FTMTree_MT *tree2) {
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(p2);
        d = d1 > b1 ? d1 - b1 : b1 - d1;
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        d = d1 > b1 ? d1 - b1 : b1 - d1;
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(p2);
        dataType dist1 = d1 > b1 ? d1 - b1 : b1 - d1;
        dataType dist2 = d2 > b2 ? d2 - b2 : b2 - d2;
        d = dist1 > dist2 ? dist1 - dist2 : dist2 - dist1;
      }
      return squared_ ? d * d : d;
    }

    template <class dataType>
    inline dataType editCost_Shifting(int n1,
                                      int p1,
                                      int n2,
                                      int p2,
                                      ftm::FTMTree_MT *tree1,
                                      ftm::FTMTree_MT *tree2) {
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(p2);
        d = d1 > b1 ? d1 - b1 : b1 - d1;
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        d = d1 > b1 ? d1 - b1 : b1 - d1;
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(p2);
        dataType pers1 = d1 > b1 ? d1 - b1 : b1 - d1;
        dataType pers2 = d2 > b2 ? d2 - b2 : b2 - d2;
        dataType db = b1 > b2 ? b1 - b2 : b2 - b1;
        dataType dp = pers1 > pers2 ? pers1 - pers2 : pers2 - pers1;
        d = db + dp;
      }
      return squared_ ? d * d : d;
    }

  public:
    BranchMappingDistance() {
      this->setDebugMsgPrefix(
        "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                              // the beginning of every msg
    }
    ~BranchMappingDistance() override = default;

    void setBaseMetric(int m) {
      baseMetric_ = m;
    }

    void setAssignmentSolver(int assignmentSolver) {
      assignmentSolverID_ = assignmentSolver;
    }

    void setSquared(bool s) {
      squared_ = s;
    }

    template <class dataType>
    dataType editDistance_branch(ftm::FTMTree_MT *tree1,
                                 ftm::FTMTree_MT *tree2) {

      // initialize memoization tables

      std::vector<std::vector<int>> predecessors1(tree1->getNumberOfNodes());
      std::vector<std::vector<int>> predecessors2(tree2->getNumberOfNodes());
      int const rootID1 = tree1->getRoot();
      int const rootID2 = tree2->getRoot();
      std::vector<int> preorder1(tree1->getNumberOfNodes());
      std::vector<int> preorder2(tree2->getNumberOfNodes());

      int depth1 = 0;
      int depth2 = 0;
      std::stack<int> stack;
      stack.push(rootID1);
      int count = tree1->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int const nIdx = stack.top();
        stack.pop();
        preorder1[count] = nIdx;
        count--;
        depth1 = std::max((int)predecessors1[nIdx].size(), depth1);
        std::vector<ftm::idNode> children;
        tree1->getChildren(nIdx, children);
        for(int const cIdx : children) {
          stack.push(cIdx);
          predecessors1[cIdx].reserve(predecessors1[nIdx].size() + 1);
          predecessors1[cIdx].insert(predecessors1[cIdx].end(),
                                     predecessors1[nIdx].begin(),
                                     predecessors1[nIdx].end());
          predecessors1[cIdx].push_back(nIdx);
        }
      }
      stack.push(rootID2);
      count = tree2->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int const nIdx = stack.top();
        stack.pop();
        preorder2[count] = nIdx;
        count--;
        depth2 = std::max((int)predecessors2[nIdx].size(), depth2);
        std::vector<ftm::idNode> children;
        tree2->getChildren(nIdx, children);
        for(int const cIdx : children) {
          stack.push(cIdx);
          predecessors2[cIdx].reserve(predecessors2[nIdx].size() + 1);
          predecessors2[cIdx].insert(predecessors2[cIdx].end(),
                                     predecessors2[nIdx].begin(),
                                     predecessors2[nIdx].end());
          predecessors2[cIdx].push_back(nIdx);
        }
      }

      size_t nn1 = tree1->getNumberOfNodes();
      size_t nn2 = tree2->getNumberOfNodes();
      size_t const dim1 = 1;
      size_t const dim2 = (nn1 + 1) * dim1;
      size_t const dim3 = (depth1 + 1) * dim2;
      size_t const dim4 = (nn2 + 1) * dim3;

      std::vector<dataType> memT((nn1 + 1) * (depth1 + 1) * (nn2 + 1)
                                 * (depth2 + 1));

      memT[nn1 + 0 * dim2 + nn2 * dim3 + 0 * dim4] = 0;
      for(size_t i = 0; i < nn1; i++) {
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);
        for(size_t l = 1; l <= predecessors1[preorder1[i]].size(); l++) {
          int parent1 = predecessors1[preorder1[i]]
                                     [predecessors1[preorder1[i]].size() - l];

          //-----------------------------------------------------------------------
          // If first subtree has only one branch, return deletion cost of this
          // branch
          if(tree1->getNumberOfChildren(curr1) == 0) {
            memT[curr1 + l * dim2 + nn2 * dim3 + 0 * dim4]
              = this->baseMetric_ == 0 ? editCost_Wasserstein1<dataType>(
                  curr1, parent1, -1, -1, tree1, tree2)
                : this->baseMetric_ == 1 ? editCost_Wasserstein2<dataType>(
                    curr1, parent1, -1, -1, tree1, tree2)
                : this->baseMetric_ == 2
                  ? editCost_Persistence<dataType>(
                    curr1, parent1, -1, -1, tree1, tree2)
                  : editCost_Shifting<dataType>(
                    curr1, parent1, -1, -1, tree1, tree2);
          }
          //-----------------------------------------------------------------------
          // If first subtree has more than one branch, try all decompositions
          else {
            dataType c = std::numeric_limits<dataType>::max();
            for(auto child1_mb : children1) {
              dataType c_
                = memT[child1_mb + (l + 1) * dim2 + nn2 * dim3 + 0 * dim4];
              for(auto child1 : children1) {
                if(child1 == child1_mb) {
                  continue;
                }
                c_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
              }
              c = std::min(c, c_);
            }
            memT[curr1 + l * dim2 + nn2 * dim3 + 0 * dim4] = c;
          }
        }
      }
      for(size_t j = 0; j < nn2; j++) {
        int curr2 = preorder2[j];
        std::vector<ftm::idNode> children2;
        tree2->getChildren(curr2, children2);
        for(size_t l = 1; l <= predecessors2[preorder2[j]].size(); l++) {
          int parent2 = predecessors2[preorder2[j]]
                                     [predecessors2[preorder2[j]].size() - l];

          //-----------------------------------------------------------------------
          // If first subtree has only one branch, return deletion cost of this
          // branch
          if(tree2->getNumberOfChildren(curr2) == 0) {
            memT[nn1 + 0 * dim2 + curr2 * dim3 + l * dim4]
              = this->baseMetric_ == 0 ? editCost_Wasserstein1<dataType>(
                  -1, -1, curr2, parent2, tree1, tree2)
                : this->baseMetric_ == 1 ? editCost_Wasserstein2<dataType>(
                    -1, -1, curr2, parent2, tree1, tree2)
                : this->baseMetric_ == 2
                  ? editCost_Persistence<dataType>(
                    -1, -1, curr2, parent2, tree1, tree2)
                  : editCost_Shifting<dataType>(
                    -1, -1, curr2, parent2, tree1, tree2);
          }
          //-----------------------------------------------------------------------
          // If first subtree has more than one branch, try all decompositions
          else {
            dataType c = std::numeric_limits<dataType>::max();
            for(auto child2_mb : children2) {
              dataType c_
                = memT[nn1 + 0 * dim2 + child2_mb * dim3 + (l + 1) * dim4];
              for(auto child2 : children2) {
                if(child2 == child2_mb) {
                  continue;
                }
                c_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
              }
              c = std::min(c, c_);
            }
            memT[nn1 + 0 * dim2 + curr2 * dim3 + l * dim4] = c;
          }
        }
      }

      for(size_t i = 0; i < nn1; i++) {
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);
        for(size_t j = 0; j < nn2; j++) {
          int curr2 = preorder2[j];
          std::vector<ftm::idNode> children2;
          tree2->getChildren(curr2, children2);
          for(size_t l1 = 1; l1 <= predecessors1[preorder1[i]].size(); l1++) {
            int parent1
              = predecessors1[preorder1[i]]
                             [predecessors1[preorder1[i]].size() - l1];
            for(size_t l2 = 1; l2 <= predecessors2[preorder2[j]].size(); l2++) {
              int parent2
                = predecessors2[preorder2[j]]
                               [predecessors2[preorder2[j]].size() - l2];

              //===============================================================================
              // If both trees not empty, find optimal edit operation

              //---------------------------------------------------------------------------
              // If both trees only have one branch, return edit cost between
              // the two branches
              if(tree1->getNumberOfChildren(curr1) == 0
                 and tree2->getNumberOfChildren(curr2) == 0) {
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4]
                  = this->baseMetric_ == 0 ? editCost_Wasserstein1<dataType>(
                      curr1, parent1, curr2, parent2, tree1, tree2)
                    : this->baseMetric_ == 1 ? editCost_Wasserstein2<dataType>(
                        curr1, parent1, curr2, parent2, tree1, tree2)
                    : this->baseMetric_ == 2
                      ? editCost_Persistence<dataType>(
                        curr1, parent1, curr2, parent2, tree1, tree2)
                      : editCost_Shifting<dataType>(
                        curr1, parent1, curr2, parent2, tree1, tree2);
              }
              //---------------------------------------------------------------------------
              // If first tree only has one branch, try all decompositions of
              // second tree
              else if(children1.size() == 0) {
                dataType d = std::numeric_limits<dataType>::max();
                for(auto child2_mb : children2) {
                  dataType d_ = memT[curr1 + l1 * dim2 + child2_mb * dim3
                                     + (l2 + 1) * dim4];
                  for(auto child2 : children2) {
                    if(child2 == child2_mb) {
                      continue;
                    }
                    d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
              //---------------------------------------------------------------------------
              // If second tree only has one branch, try all decompositions of
              // first tree
              else if(children2.size() == 0) {
                dataType d = std::numeric_limits<dataType>::max();
                for(auto child1_mb : children1) {
                  dataType d_ = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3
                                     + l2 * dim4];
                  for(auto child1 : children1) {
                    if(child1 == child1_mb) {
                      continue;
                    }
                    d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
              //---------------------------------------------------------------------------
              // If both trees have more than one branch, try all decompositions
              // of both trees
              else {
                dataType d = std::numeric_limits<dataType>::max();
                //-----------------------------------------------------------------------
                // Try all possible main branches of first tree (child1_mb) and
                // all possible main branches of second tree (child2_mb) Then
                // try all possible matchings of subtrees
                if(children1.size() == 2 && children2.size() == 2) {
                  int const child11 = children1[0];
                  int const child12 = children1[1];
                  int const child21 = children2[0];
                  int const child22 = children2[1];
                  d = std::min<dataType>(
                    d,
                    memT[child11 + (l1 + 1) * dim2 + child21 * dim3
                         + (l2 + 1) * dim4]
                      + memT[child12 + 1 * dim2 + child22 * dim3 + 1 * dim4]);
                  d = std::min<dataType>(
                    d,
                    memT[child12 + (l1 + 1) * dim2 + child22 * dim3
                         + (l2 + 1) * dim4]
                      + memT[child11 + 1 * dim2 + child21 * dim3 + 1 * dim4]);
                  d = std::min<dataType>(
                    d,
                    memT[child11 + (l1 + 1) * dim2 + child22 * dim3
                         + (l2 + 1) * dim4]
                      + memT[child12 + 1 * dim2 + child21 * dim3 + 1 * dim4]);
                  d = std::min<dataType>(
                    d,
                    memT[child12 + (l1 + 1) * dim2 + child21 * dim3
                         + (l2 + 1) * dim4]
                      + memT[child11 + 1 * dim2 + child22 * dim3 + 1 * dim4]);
                } else {
                  for(auto child1_mb : children1) {
                    auto topo1_ = children1;
                    topo1_.erase(
                      std::remove(topo1_.begin(), topo1_.end(), child1_mb),
                      topo1_.end());
                    for(auto child2_mb : children2) {
                      auto topo2_ = children2;
                      topo2_.erase(
                        std::remove(topo2_.begin(), topo2_.end(), child2_mb),
                        topo2_.end());

                      auto f = [&](unsigned r, unsigned c) {
                        int const c1 = r < topo1_.size() ? topo1_[r] : -1;
                        int const c2 = c < topo2_.size() ? topo2_[c] : -1;
                        return memT[c1 + 1 * dim2 + c2 * dim3 + 1 * dim4];
                      };
                      int size = std::max(topo1_.size(), topo2_.size()) + 1;
                      auto costMatrix = std::vector<std::vector<dataType>>(
                        size, std::vector<dataType>(size, 0));
                      std::vector<MatchingType> matching;
                      for(int r = 0; r < size; r++) {
                        for(int c = 0; c < size; c++) {
                          costMatrix[r][c] = f(r, c);
                        }
                      }

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
                          assignmentSolver = &solverAuction;
                      }
                      assignmentSolver->setInput(costMatrix);
                      assignmentSolver->setBalanced(true);
                      assignmentSolver->run(matching);
                      dataType d_ = memT[child1_mb + (l1 + 1) * dim2
                                         + child2_mb * dim3 + (l2 + 1) * dim4];
                      for(auto m : matching)
                        d_ += std::get<2>(m);
                      d = std::min(d, d_);
                    }
                  }
                }
                //-----------------------------------------------------------------------
                // Try to continue main branch on one child of first tree and
                // delete all other subtrees Then match continued branch to
                // current branch in second tree
                for(auto child1_mb : children1) {
                  dataType d_ = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3
                                     + l2 * dim4];
                  for(auto child1 : children1) {
                    if(child1 == child1_mb) {
                      continue;
                    }
                    d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
                  }
                  d = std::min(d, d_);
                }
                //-----------------------------------------------------------------------
                // Try to continue main branch on one child of second tree and
                // delete all other subtrees Then match continued branch to
                // current branch in first tree
                for(auto child2_mb : children2) {
                  dataType d_ = memT[curr1 + l1 * dim2 + child2_mb * dim3
                                     + (l2 + 1) * dim4];
                  for(auto child2 : children2) {
                    if(child2 == child2_mb) {
                      continue;
                    }
                    d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
            }
          }
        }
      }

      std::vector<ftm::idNode> children1;
      tree1->getChildren(rootID1, children1);
      std::vector<ftm::idNode> children2;
      tree2->getChildren(rootID2, children2);

      dataType res
        = memT[children1[0] + 1 * dim2 + children2[0] * dim3 + 1 * dim4];

      return squared_ ? std::sqrt(res) : res;
    }
  };

} // namespace ttk