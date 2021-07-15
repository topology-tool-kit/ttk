#include <FTMTreeUtils.h>

namespace ttk {
  namespace ftm {

    // --------------------
    // Is
    // --------------------

    bool isNodeOriginDefined(FTMTree_MT *tree, idNode nodeId) {
      unsigned int origin = (unsigned int)tree->getNode(nodeId)->getOrigin();
      return origin != nullNodes and origin < tree->getNumberOfNodes()
             and origin >= 0;
      ;
    }

    bool isRoot(FTMTree_MT *tree, idNode nodeId) {
      /*if(nodeId >= tree->getNumberOfNodes())
        std::cout << nodeId << " _____ " << tree->getNumberOfNodes() <<
        std::endl;*/
      return tree->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
    }

    bool isLeaf(FTMTree_MT *tree, idNode nodeId) {
      return tree->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
    }

    bool isNodeAlone(FTMTree_MT *tree, idNode nodeId) {
      return isRoot(tree, nodeId) and isLeaf(tree, nodeId);
    }

    bool isFullMerge(FTMTree_MT *tree) {
      idNode treeRoot = getRoot(tree);
      return (unsigned int)tree->getNode(treeRoot)->getOrigin() == treeRoot;
    }

    bool isBranchOrigin(FTMTree_MT *tree, idNode nodeId) {
      return getParent(tree, tree->getNode(nodeId)->getOrigin()) != nodeId;
    }

    // --------------------
    // Get
    // --------------------

    idNode getRoot(FTMTree_MT *tree) {
      for(idNode node = 0; node < tree->getNumberOfNodes(); ++node)
        if(isRoot(tree, node) and !isLeaf(tree, node))
          return node;
      return nullNodes;
    }

    idNode getParent(FTMTree_MT *tree, idNode nodeId) {
      if(!isRoot(tree, nodeId)) {
        // _ Nodes in merge trees should have only one parent
        idSuperArc arcId = tree->getNode(nodeId)->getUpSuperArcId(0);
        idNode parentNodeId = tree->getSuperArc(arcId)->getUpNodeId();
        return parentNodeId;
      }
      return nodeId;
    }

    std::vector<idNode> getChildren(FTMTree_MT *tree, idNode nodeId) {
      std::vector<idNode> childrens;
      for(idSuperArc i = 0;
          i < tree->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i) {
        idSuperArc arcId = tree->getNode(nodeId)->getDownSuperArcId(i);
        childrens.push_back(tree->getSuperArc(arcId)->getDownNodeId());
      }
      return childrens;
    }

    std::vector<idNode> getLeaves(FTMTree_MT *tree) {
      std::vector<idNode> treeLeaves;
      for(idNode i = 0; i < tree->getNumberOfNodes(); ++i) {
        if(isLeaf(tree, i) and !isRoot(tree, i))
          treeLeaves.push_back(i);
      }
      return treeLeaves;
    }

    int getNumberOfLeaves(FTMTree_MT *tree) {
      auto leaves = getLeaves(tree);
      return leaves.size();
    }

    int getNumberOfNodeAlone(FTMTree_MT *tree) {
      int cpt = 0;
      for(idNode i = 0; i < tree->getNumberOfNodes(); ++i)
        cpt += isNodeAlone(tree, i) ? 1 : 0;
      return cpt;
    }

    int getRealNumberOfNodes(FTMTree_MT *tree) {
      return tree->getNumberOfNodes() - getNumberOfNodeAlone(tree);
    }

    std::tuple<std::vector<idNode>, std::vector<idNode>>
      getBranchOriginsFromThisBranch(FTMTree_MT *tree, idNode node) {
      std::vector<idNode> branchOrigins, nonBranchOrigins;

      idNode nodeOrigin = tree->getNode(node)->getOrigin();
      idNode nodeParent = getParent(tree, nodeOrigin);
      while(nodeParent != node) {
        if(isBranchOrigin(tree, nodeParent))
          branchOrigins.push_back(nodeParent);
        else
          nonBranchOrigins.push_back(nodeParent);
        nodeParent = getParent(tree, nodeParent);
      }

      return std::make_tuple(branchOrigins, nonBranchOrigins);
    }

    void getTreeBranching(FTMTree_MT *tree,
                          std::vector<idNode> &branching,
                          std::vector<int> &branchingID,
                          std::vector<std::vector<idNode>> &nodeBranching) {
      branching = std::vector<idNode>(tree->getNumberOfNodes());
      branchingID = std::vector<int>(tree->getNumberOfNodes(), -1);
      nodeBranching
        = std::vector<std::vector<idNode>>(tree->getNumberOfNodes());
      int branchID = 0;
      std::queue<idNode> queue;
      queue.push(getRoot(tree));
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        if(isLeaf(tree, node))
          continue;
        auto nodeOrigin = tree->getNode(node)->getOrigin();
        idNode parentNodeOrigin = nodeOrigin;
        // idNode oldParentNodeOrigin;
        while(parentNodeOrigin != node) {
          branching[parentNodeOrigin] = node;
          branchingID[parentNodeOrigin] = branchID;
          nodeBranching[node].push_back(parentNodeOrigin);
          // oldParentNodeOrigin = parentNodeOrigin;
          parentNodeOrigin = getParent(tree, parentNodeOrigin);
          /*if(oldParentNodeOrigin == parentNodeOrigin) {
            printTree(tree);
            std::cout << getRoot(tree) << " _ " << parentNodeOrigin <<
          std::endl; std::cout << node << " _ " << nodeOrigin << std::endl;
            std::cout << "getTreeBranching" << std::endl;
            myPause();
          }*/
        }
        if(isRoot(tree, node)) {
          branching[node] = node;
          branchingID[node] = branchID;
        }
        ++branchID;
        for(auto child : getChildren(tree, node))
          queue.push(child);
      }
    }

    void getTreeBranching(FTMTree_MT *tree,
                          std::vector<idNode> &branching,
                          std::vector<int> &branchingID) {
      std::vector<std::vector<idNode>> nodeBranching;
      getTreeBranching(tree, branching, branchingID, nodeBranching);
    }

  } // namespace ftm
} // namespace ttk
