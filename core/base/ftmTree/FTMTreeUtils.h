/// \ingroup base
/// \class FTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#pragma once

#include <FTMTree.h>

namespace ttk {
  namespace ftm {

    struct MergeTree {
      ftm::Scalars scalars;
      ftm::Params params;
      ftm::FTMTree_MT tree;
      MergeTree(ftm::Scalars scalarsT, ftm::Params paramsT)
        : scalars(scalarsT), params(paramsT),
          tree(&params, &scalars, params.treeType) {
        tree.makeAlloc();
      }
    };

    // --------------------
    // Is
    // --------------------

    bool isNodeOriginDefined(FTMTree_MT *tree, idNode nodeId);

    bool isRoot(FTMTree_MT *tree, idNode nodeId);

    bool isLeaf(FTMTree_MT *tree, idNode nodeId);

    bool isNodeAlone(FTMTree_MT *tree, idNode nodeId);

    bool isFullMerge(FTMTree_MT *tree);

    bool isBranchOrigin(FTMTree_MT *tree, idNode nodeId);

    // --------------------
    // Get
    // --------------------

    idNode getRoot(FTMTree_MT *tree);

    idNode getParent(FTMTree_MT *tree, idNode nodeId);

    std::vector<idNode> getChildren(FTMTree_MT *tree, idNode nodeId);

    std::vector<idNode> getLeaves(FTMTree_MT *tree);

    int getNumberOfLeaves(FTMTree_MT *tree);

    int getNumberOfNodeAlone(FTMTree_MT *tree);

    int getRealNumberOfNodes(FTMTree_MT *tree);

    std::tuple<std::vector<idNode>, std::vector<idNode>>
      getBranchOriginsFromThisBranch(FTMTree_MT *tree, idNode node);

    void getTreeBranching(FTMTree_MT *tree,
                          std::vector<idNode> &branching,
                          std::vector<int> &branchingID,
                          std::vector<std::vector<idNode>> &nodeBranching);

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
    std::tuple<dataType, dataType> getBirthDeath(FTMTree_MT *tree,
                                                 idNode nodeId) {
      idNode originId = tree->getNode(nodeId)->getOrigin();
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
    dataType getBirth(FTMTree_MT *tree, idNode nodeId) {
      return std::get<0>(getBirthDeath<dataType>(tree, nodeId));
    }

    template <class dataType>
    dataType getNodePersistence(FTMTree_MT *tree, idNode nodeId) {
      std::tuple<dataType, dataType> birthDeath
        = getBirthDeath<dataType>(tree, nodeId);
      return std::get<1>(birthDeath) - std::get<0>(birthDeath);
    }

    // --------------------
    // Is
    // --------------------

    template <class dataType>
    bool isJoinTree(FTMTree_MT *tree) {
      auto root = getRoot(tree);
      idNode child = getChildren(tree, root)[0];
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
    bool isImportantPair(FTMTree_MT *tree, idNode nodeId, double threshold) {
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
    static type myAbs(const type var) {
      return (var >= 0) ? var : -var;
    }

    template <class dataType>
    bool isEqual(dataType first, dataType two, double eps = 1e-6) {
      return myAbs<dataType>(first - two)
             < eps * std::max(myAbs<dataType>(first), myAbs<dataType>(two));
    }

  } // namespace ftm
} // namespace ttk

#endif
