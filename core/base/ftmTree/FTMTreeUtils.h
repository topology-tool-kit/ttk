/// \ingroup base
/// \class FTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#pragma once

#include <FTMTree_MT.h>
#include <stack>

namespace ttk {
  namespace ftm {

    // --------------------
    // Utils
    // --------------------
    void printTreesStats(std::vector<ftm::FTMTree_MT *> &trees);

    template <class dataType>
    void printTree(MergeTree<dataType> &tree, bool doPrint = true) {
      tree.tree.printTree(doPrint);
    }

    template <class dataType>
    void printTreeStats(MergeTree<dataType> &tree) {
      tree.tree.printTreeStats();
    }

    template <class dataType>
    void printTreeScalars(MergeTree<dataType> &tree,
                          bool printNodeAlone = true) {
      printTreeScalars<dataType>(&(tree.tree), printNodeAlone);
    }

    template <class dataType>
    void mergeTreeToFTMTree(std::vector<MergeTree<dataType>> &trees,
                            std::vector<ftm::FTMTree_MT *> &treesT) {
      treesT.clear();
      for(MergeTree<dataType> &t : trees)
        treesT.push_back(&(t.tree));
    }

    template <class dataType>
    void mergeTreeTemplateToDouble(MergeTree<dataType> &mt,
                                   MergeTree<double> &newMt) {
      std::vector<double> newScalarsValues;
      for(auto val : mt.scalarsValues)
        newScalarsValues.push_back(static_cast<double>(val));
      newMt = MergeTree<double>(mt.scalars, newScalarsValues, mt.params);
      newMt.tree.copyMergeTreeStructure(&(mt.tree));
    }

    template <class dataType>
    void mergeTreesTemplateToDouble(std::vector<MergeTree<dataType>> &mts,
                                    std::vector<MergeTree<double>> &newMts) {
      newMts.clear();
      for(auto &mt : mts) {
        MergeTree<double> newMt;
        mergeTreeTemplateToDouble<dataType>(mt, newMt);
        newMts.push_back(newMt);
      }
    }

    template <class dataType>
    void mergeTreesTemplateToDouble(
      std::vector<std::vector<MergeTree<dataType>>> &mts,
      std::vector<std::vector<MergeTree<double>>> &newMts) {
      newMts.clear();
      for(auto &mt : mts) {
        std::vector<MergeTree<double>> newMt;
        mergeTreesTemplateToDouble<dataType>(mt, newMt);
        newMts.push_back(newMt);
      }
    }

    template <class dataType>
    void mergeTreeDoubleToTemplate(MergeTree<double> &mt,
                                   MergeTree<dataType> &newMt) {
      std::vector<dataType> newScalarsValues;
      for(auto val : mt.scalarsValues)
        newScalarsValues.push_back(static_cast<dataType>(val));
      newMt = MergeTree<dataType>(mt.scalars, newScalarsValues, mt.params);
      newMt.tree.copyMergeTreeStructure(&(mt.tree));
    }

    template <class dataType>
    void mergeTreesDoubleToTemplate(std::vector<MergeTree<double>> &mts,
                                    std::vector<MergeTree<dataType>> &newMts) {
      newMts.clear();
      for(auto &mt : mts) {
        MergeTree<dataType> newMt;
        mergeTreeDoubleToTemplate<dataType>(mt, newMt);
        newMts.push_back(newMt);
      }
    }

    template <class dataType>
    void mergeTreesDoubleToTemplate(
      std::vector<std::vector<MergeTree<double>>> &mts,
      std::vector<std::vector<MergeTree<dataType>>> &newMts) {
      newMts.clear();
      for(auto &mt : mts) {
        std::vector<MergeTree<dataType>> newMt;
        mergeTreesDoubleToTemplate<dataType>(mt, newMt);
        newMts.push_back(newMt);
      }
    }

    // --------------------
    // Make tree utils
    // --------------------
    void manageInconsistentArcsMultiParent(FTMTree_MT *tree);

    void removeSelfLink(FTMTree_MT *tree);

    // --------------------
    // MergeTree
    // --------------------
    template <class dataType>
    MergeTree<dataType> createEmptyMergeTree(int scalarSize) {
      // Init Scalars
      ftm::Scalars scalars;
      scalars.size = scalarSize;
      dataType *scalarsValues = nullptr;
      scalars.values = (void *)scalarsValues;

      // Init Params
      ftm::Params params;
      params.treeType = ftm::Join_Split;

      // Init tree
      MergeTree<dataType> mergeTree(scalars, params);

      return mergeTree;
    }

    template <class dataType>
    void setTreeScalars(MergeTree<dataType> &mergeTree,
                        std::vector<dataType> &scalarsVector) {
      mergeTree.scalarsValues = scalarsVector;
      mergeTree.scalars.values = (void *)(mergeTree.scalarsValues.data());
      mergeTree.scalars.size = mergeTree.scalarsValues.size();
    }

    template <class dataType>
    void getTreeScalars(ftm::FTMTree_MT *tree,
                        std::vector<dataType> &scalarsVector) {
      scalarsVector.clear();
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        scalarsVector.push_back(tree->getValue<dataType>(i));
    }

    template <class dataType>
    void getTreeScalars(MergeTree<dataType> &mergeTree,
                        std::vector<dataType> &scalarsVector) {
      getTreeScalars<dataType>(&(mergeTree.tree), scalarsVector);
    }

    template <class dataType>
    MergeTree<dataType> copyMergeTree(ftm::FTMTree_MT *tree,
                                      bool doSplitMultiPersPairs = false) {
      std::vector<dataType> scalarsVector;
      getTreeScalars<dataType>(tree, scalarsVector);

      // Get multi persistence pairs
      std::vector<ftm::idNode> multiPersOrigins;
      if(doSplitMultiPersPairs) {
        multiPersOrigins = tree->getMultiPersOrigins<dataType>(true);
        for(ftm::idNode nodeOrigin : multiPersOrigins) {
          scalarsVector[nodeOrigin]
            = tree->getValue<dataType>(tree->getNode(nodeOrigin)->getOrigin());
          scalarsVector.push_back(tree->getValue<dataType>(nodeOrigin));
        }
      }

      // Create new tree
      MergeTree<dataType> mergeTree
        = createEmptyMergeTree<dataType>(scalarsVector.size());
      ftm::FTMTree_MT *treeNew = &(mergeTree.tree);
      setTreeScalars<dataType>(mergeTree, scalarsVector);

      // Copy tree structure
      treeNew->copyMergeTreeStructure(tree);

      // Add multi persistence nodes origins
      if(doSplitMultiPersPairs) {
        for(ftm::idNode nodeOrigin : multiPersOrigins) {
          int nodeCpt = treeNew->getNumberOfNodes();
          treeNew->makeNode(nodeCpt);
          treeNew->getNode(nodeCpt)->setOrigin(nodeOrigin);
          treeNew->getNode(nodeOrigin)->setOrigin(nodeCpt);
        }
      }

      return mergeTree;
    }

    template <class dataType>
    MergeTree<dataType> copyMergeTree(MergeTree<dataType> &mergeTree,
                                      bool doSplitMultiPersPairs = false) {
      return copyMergeTree<dataType>(&(mergeTree.tree), doSplitMultiPersPairs);
    }

    // Remove unused nodes
    template <class dataType>
    MergeTree<dataType> cleanMergeTree(ftm::FTMTree_MT *tree,
                                       std::vector<int> &nodeCorr,
                                       bool useBD = true) {
      // Create new tree
      int newNoNodes = tree->getRealNumberOfNodes() * 2;
      MergeTree<dataType> mTreeNew = createEmptyMergeTree<dataType>(newNoNodes);
      ftm::FTMTree_MT *treeNew = &(mTreeNew.tree);
      std::vector<dataType> newScalarsVector(newNoNodes, 0);

      // Copy the old tree structure
      std::vector<int> nodeDone(tree->getNumberOfNodes(), 0);
      nodeCorr = std::vector<int>(tree->getNumberOfNodes(), -1);
      std::vector<std::vector<ftm::idNode>> treeMultiPers;
      if(not useBD)
        tree->getMultiPersOriginsVectorFromTree(treeMultiPers);
      std::queue<ftm::idNode> queue;
      std::vector<idNode> leaves;
      tree->getLeavesFromTree(leaves);
      for(auto leaf : leaves)
        queue.push(leaf);
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();
        ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
        if(tree->isRoot(node) and tree->isFullMerge())
          nodeOrigin = tree->getMergedRootOrigin<dataType>();

        if(useBD) {
          int nodeOriginIndex = treeNew->getNumberOfNodes();
          if(nodeCorr[nodeOrigin] == -1)
            treeNew->makeNode(nodeOriginIndex);
          else
            nodeOriginIndex = nodeCorr[nodeOrigin];
          int nodeIndex = treeNew->getNumberOfNodes();
          if(nodeCorr[node] == -1)
            treeNew->makeNode(nodeIndex);
          else
            nodeIndex = nodeCorr[node];

          if(nodeCorr[nodeOrigin] == -1)
            treeNew->getNode(nodeOriginIndex)->setOrigin(nodeIndex);
          treeNew->getNode(nodeIndex)->setOrigin(nodeOriginIndex);

          newScalarsVector[nodeOriginIndex]
            = tree->getValue<dataType>(nodeOrigin);
          newScalarsVector[nodeIndex] = tree->getValue<dataType>(node);
          nodeCorr[nodeOrigin] = nodeOriginIndex;
          nodeCorr[node] = nodeIndex;
        } else {
          int nodeCpt = treeNew->getNumberOfNodes();
          treeNew->makeNode(nodeCpt);
          if(!tree->isLeaf(node)) {
            treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
            if(not(tree->isRoot(node) and node == nodeOrigin))
              treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
            for(auto nodeMultiPers : treeMultiPers[node])
              treeNew->getNode(nodeCorr[nodeMultiPers])->setOrigin(nodeCpt);
          } else if(tree->isNodeAlone(nodeOrigin)) { // saddle merged
            treeNew->makeNode(nodeCpt + 1);
            newScalarsVector[nodeCpt + 1]
              = tree->getValue<dataType>(nodeOrigin);
            nodeCorr[nodeOrigin] = nodeCpt + 1;
            treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
            treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
          }
          newScalarsVector[nodeCpt] = tree->getValue<dataType>(node);
          nodeCorr[node] = nodeCpt;
        }

        std::vector<idNode> children;
        tree->getChildren(node, children);
        for(idNode child : children)
          treeNew->makeSuperArc(nodeCorr[child], nodeCorr[node]);

        if(!tree->isRoot(node)) {
          ftm::idNode parent = tree->getParentSafe(node);
          nodeDone[parent] += 1;
          if(nodeDone[parent] == tree->getNumberOfChildren(parent))
            queue.push(parent);
        }
      }

      // Manage full merge
      auto treeRoot = tree->getRoot();
      if(tree->getNode(treeRoot)->getOrigin() == (int)treeRoot) {
        auto treeNewRoot = treeNew->getRoot();
        auto treeNewRootOrigin = treeNew->getNode(treeNewRoot)->getOrigin();
        newScalarsVector[treeNewRootOrigin]
          = tree->getValue<dataType>(tree->getMergedRootOrigin<dataType>());
        treeNew->getNode(treeNewRoot)->setOrigin(treeNewRoot);
      }

      // Set new scalars
      setTreeScalars<dataType>(mTreeNew, newScalarsVector);

      // Return new tree
      return mTreeNew;
    }

    template <class dataType>
    void cleanMergeTree(MergeTree<dataType> &mTree,
                        std::vector<int> &nodeCorr,
                        bool useBD = true) {
      mTree = cleanMergeTree<dataType>(&(mTree.tree), nodeCorr, useBD);
    }

    template <class dataType>
    void cleanMergeTree(MergeTree<dataType> &mTree, bool useBD = true) {
      std::vector<int> nodeCorr;
      cleanMergeTree<dataType>(mTree, nodeCorr, useBD);
    }

  } // namespace ftm
} // namespace ttk

#endif
