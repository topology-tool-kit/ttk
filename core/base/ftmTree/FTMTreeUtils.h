/// \ingroup base
/// \class FTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#pragma once

// --------------------
// Is
// --------------------

bool isNodeOriginDefined(idNode nodeId);

bool isRoot(idNode nodeId);

bool isLeaf(idNode nodeId);

bool isNodeAlone(idNode nodeId);

bool isFullMerge();

bool isBranchOrigin(idNode nodeId);

// --------------------
// Get
// --------------------

idNode getRoot();

idNode getParent(idNode nodeId);

std::vector<idNode> getChildren(idNode nodeId);

std::vector<idNode> getLeavesFromTree();

int getNumberOfLeavesFromTree();

int getNumberOfNodeAlone();

int getRealNumberOfNodes();

std::tuple<std::vector<idNode>, std::vector<idNode>>
  getBranchOriginsFromThisBranch(idNode node);

void getTreeBranching(std::vector<idNode> &branching,
                      std::vector<int> &branchingID,
                      std::vector<std::vector<idNode>> &nodeBranching);

void getTreeBranching(std::vector<idNode> &branching,
                      std::vector<int> &branchingID);

// ----------------------------------------
// Template functions
// ----------------------------------------

// --------------------
// Get
// --------------------

template <class dataType>
std::tuple<dataType, dataType> getBirthDeath(idNode nodeId);

template <class dataType>
dataType getBirth(idNode nodeId);

template <class dataType>
dataType getNodePersistence(idNode nodeId);

// --------------------
// Is
// --------------------

template <class dataType>
bool isJoinTree();

template <class dataType>
bool isImportantPair(idNode nodeId, double threshold);

#endif
