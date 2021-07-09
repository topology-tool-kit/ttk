/// \ingroup base
/// \class ttk::FTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#pragma once

#include <FTMTree.h>

//#include <ttkUtils.h>

using namespace ttk;
using namespace ftm;

// --------------------
// Is
// --------------------

bool isNodeOriginDefined(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isRoot(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isLeaf(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isNodeAlone(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

bool isFullMerge(ftm::FTMTree_MT *tree);

bool isBranchOrigin(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

// --------------------
// Get
// --------------------

ftm::idNode getRoot(ftm::FTMTree_MT *tree);

ftm::idNode getParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

std::vector<ftm::idNode> getChildren(ftm::FTMTree_MT *tree, ftm::idNode nodeId);

std::vector<ttk::ftm::idNode> getLeaves(ttk::ftm::FTMTree_MT *tree);

int getNumberOfLeaves(ftm::FTMTree_MT *tree);

int getNumberOfNodeAlone(ftm::FTMTree_MT *tree);

int getRealNumberOfNodes(ftm::FTMTree_MT *tree);

std::tuple<std::vector<ftm::idNode>, std::vector<ftm::idNode>>
  getBranchOriginsFromThisBranch(FTMTree_MT *tree, ftm::idNode node);

void getTreeBranching(FTMTree_MT *tree,
                      std::vector<idNode> &branching,
                      std::vector<int> &branchingID,
                      std::vector<std::vector<ftm::idNode>> &nodeBranching);

void getTreeBranching(FTMTree_MT *tree,
                      std::vector<idNode> &branching,
                      std::vector<int> &branchingID);

// ----------------------------------------
// Template functions
// ----------------------------------------

// --------------------
// Get
// --------------------

template <class dataType>
std::tuple<dataType, dataType> getBirthDeath(ftm::FTMTree_MT *tree,
                                             ftm::idNode nodeId) {
  ftm::idNode originId = tree->getNode(nodeId)->getOrigin();
  if(isNodeOriginDefined(
       tree, nodeId)) { // Avoid error if origin is not defined
    dataType pers1 = tree->getValue<dataType>(nodeId);
    dataType pers2 = tree->getValue<dataType>(originId);
    dataType birth = std::min(pers1, pers2);
    dataType death = std::max(pers1, pers2);
    return std::make_tuple(birth, death);
  }
  return std::make_tuple(0.0, 0.0);
}

template <class dataType>
dataType getBirth(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return std::get<0>(getBirthDeath<dataType>(tree, nodeId));
}

template <class dataType>
dataType getNodePersistence(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  std::tuple<dataType, dataType> birthDeath
    = getBirthDeath<dataType>(tree, nodeId);
  return std::get<1>(birthDeath) - std::get<0>(birthDeath);
}

// --------------------
// Is
// --------------------

template <class dataType>
bool isJoinTree(ftm::FTMTree_MT *tree) {
  auto root = getRoot(tree);
  ftm::idNode child = getChildren(tree, root)[0];
  if(isFullMerge(tree)) {
    dataType min = std::numeric_limits<dataType>::max();
    for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
      dataType value = tree->getValue<dataType>(i);
      if(not isNodeAlone(tree, i) and value < min) {
        min = value;
        child = i;
      }
    }
  }
  return tree->getValue<dataType>(root) > tree->getValue<dataType>(child);
}

template <class dataType>
bool isImportantPair(ftm::FTMTree_MT *tree,
                     ftm::idNode nodeId,
                     double threshold) {
  dataType rootPers = getNodePersistence<dataType>(tree, getRoot(tree));
  if(threshold > 1)
    threshold /= 100.0;
  threshold = rootPers * threshold;
  return getNodePersistence<dataType>(tree, nodeId) > threshold;
}

// --------------------
// Utils
// --------------------

template <typename type>
static type myAbs2(const type var) {
  return (var >= 0) ? var : -var;
}

template <class dataType>
bool isEqual(dataType first, dataType two, double eps = 1e-6) {
  return myAbs2<dataType>(first - two)
         < eps * std::max(myAbs2<dataType>(first), myAbs2<dataType>(two));
}

#endif
