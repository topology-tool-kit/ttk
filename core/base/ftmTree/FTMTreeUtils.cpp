#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <iostream>

namespace ttk {
  namespace ftm {

    void printTreesStats(std::vector<ftm::FTMTree_MT *> &trees) {
      for(auto tree : trees)
        tree->printTreeStats();
    }

    // --------------------
    // Is
    // --------------------
    bool FTMTree_MT::isNodeOriginDefined(idNode nodeId) {
      unsigned int origin = (unsigned int)this->getNode(nodeId)->getOrigin();
      return origin != nullNodes and origin < this->getNumberOfNodes()
             and origin >= 0;
      ;
    }

    bool FTMTree_MT::isRoot(idNode nodeId) {
      return this->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
    }

    bool FTMTree_MT::isLeaf(idNode nodeId) {
      return this->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
    }

    bool FTMTree_MT::isNodeAlone(idNode nodeId) {
      return this->isRoot(nodeId) and this->isLeaf(nodeId);
    }

    bool FTMTree_MT::isFullMerge() {
      idNode treeRoot = this->getRoot();
      return (unsigned int)this->getNode(treeRoot)->getOrigin() == treeRoot;
    }

    bool FTMTree_MT::isBranchOrigin(idNode nodeId) {
      return this->getParentSafe(this->getNode(nodeId)->getOrigin()) != nodeId;
    }

    bool FTMTree_MT::isNodeMerged(idNode nodeId) {
      bool merged = this->isNodeAlone(nodeId)
                    or this->isNodeAlone(this->getNode(nodeId)->getOrigin());
      auto nodeIdOrigin = this->getNode(nodeId)->getOrigin();
      merged
        = merged or nodeIdOrigin == this->getNode(nodeIdOrigin)->getOrigin();
      return merged;
    }

    bool FTMTree_MT::isNodeIdInconsistent(idNode nodeId) {
      return nodeId >= this->getNumberOfNodes() or nodeId < 0;
    }

    bool FTMTree_MT::isThereOnlyOnePersistencePair() {
      idNode treeRoot = this->getRoot();
      unsigned int cptNodeAlone = 0;
      idNode otherNode = treeRoot;
      for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i)
        if(this->isNodeAlone(i))
          cptNodeAlone++;
        else if(i != treeRoot)
          otherNode = i;
      // unsigned int origin = (unsigned
      // int)tree->getNode(otherNode)->getOrigin();
      idNode treeRootOrigin = this->getNode(treeRoot)->getOrigin();
      return (otherNode != treeRoot
              and this->getNumberOfNodes() - cptNodeAlone == 2
              and (treeRootOrigin == otherNode or treeRootOrigin == treeRoot));
      /*return (otherNode != treeRoot and cptNodeAlone ==
        tree->getNumberOfNodes()-2 and (origin == treeRoot or otherNode ==
        origin));*/
    }

    // Do not normalize node is if root or son of a merged root
    bool FTMTree_MT::notNeedToNormalize(idNode nodeId) {
      auto nodeIdParent = this->getParentSafe(nodeId);
      return this->isRoot(nodeId)
             or (this->isRoot(nodeIdParent)
                 and nodeIdParent
                       == (unsigned int)this->getNode(nodeIdParent)
                            ->getOrigin());
      // and nodeIdOrigin == nodeIdParent) )
    }

    bool FTMTree_MT::isMultiPersPair(idNode nodeId) {
      auto nodeOriginOrigin
        = (unsigned int)this->getNode(this->getNode(nodeId)->getOrigin())
            ->getOrigin();
      return nodeOriginOrigin != nodeId;
    }

    // --------------------
    // Get
    // --------------------
    idNode FTMTree_MT::getRoot() {
      for(idNode node = 0; node < this->getNumberOfNodes(); ++node)
        if(this->isRoot(node) and !this->isLeaf(node))
          return node;
      return nullNodes;
    }

    idNode FTMTree_MT::getParentSafe(idNode nodeId) {
      if(!this->isRoot(nodeId)) {
        // _ Nodes in merge trees should have only one parent
        idSuperArc arcId = this->getNode(nodeId)->getUpSuperArcId(0);
        idNode parentNodeId = this->getSuperArc(arcId)->getUpNodeId();
        return parentNodeId;
      }
      return nodeId;
    }

    std::vector<idNode> FTMTree_MT::getChildren(idNode nodeId) {
      std::vector<idNode> childrens;
      for(idSuperArc i = 0;
          i < this->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i) {
        idSuperArc arcId = this->getNode(nodeId)->getDownSuperArcId(i);
        childrens.push_back(this->getSuperArc(arcId)->getDownNodeId());
      }
      return childrens;
    }

    std::vector<idNode> FTMTree_MT::getLeavesFromTree() {
      std::vector<idNode> treeLeaves;
      for(idNode i = 0; i < this->getNumberOfNodes(); ++i) {
        if(this->isLeaf(i) and !this->isRoot(i))
          treeLeaves.push_back(i);
      }
      return treeLeaves;
    }

    int FTMTree_MT::getNumberOfLeavesFromTree() {
      auto leaves = this->getLeaves();
      return leaves.size();
    }

    int FTMTree_MT::getNumberOfNodeAlone() {
      int cpt = 0;
      for(idNode i = 0; i < this->getNumberOfNodes(); ++i)
        cpt += this->isNodeAlone(i) ? 1 : 0;
      return cpt;
    }

    int FTMTree_MT::getRealNumberOfNodes() {
      return this->getNumberOfNodes() - this->getNumberOfNodeAlone();
    }

    std::tuple<std::vector<idNode>, std::vector<idNode>>
      FTMTree_MT::getBranchOriginsFromThisBranch(idNode node) {
      std::vector<idNode> branchOrigins, nonBranchOrigins;

      idNode nodeOrigin = this->getNode(node)->getOrigin();
      idNode nodeParent = this->getParentSafe(nodeOrigin);
      while(nodeParent != node) {
        if(this->isBranchOrigin(nodeParent))
          branchOrigins.push_back(nodeParent);
        else
          nonBranchOrigins.push_back(nodeParent);
        nodeParent = this->getParentSafe(nodeParent);
      }

      return std::make_tuple(branchOrigins, nonBranchOrigins);
    }

    void FTMTree_MT::getTreeBranching(
      std::vector<idNode> &branching,
      std::vector<int> &branchingID,
      std::vector<std::vector<idNode>> &nodeBranching) {
      branching = std::vector<idNode>(this->getNumberOfNodes());
      branchingID = std::vector<int>(this->getNumberOfNodes(), -1);
      nodeBranching
        = std::vector<std::vector<idNode>>(this->getNumberOfNodes());
      int branchID = 0;
      std::queue<idNode> queue;
      queue.emplace(this->getRoot());
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        if(this->isLeaf(node))
          continue;
        auto nodeOrigin = this->getNode(node)->getOrigin();
        idNode parentNodeOrigin = nodeOrigin;
        while(parentNodeOrigin != node) {
          branching[parentNodeOrigin] = node;
          branchingID[parentNodeOrigin] = branchID;
          nodeBranching[node].push_back(parentNodeOrigin);
          parentNodeOrigin = this->getParentSafe(parentNodeOrigin);
        }
        if(this->isRoot(node)) {
          branching[node] = node;
          branchingID[node] = branchID;
        }
        ++branchID;
        for(auto child : this->getChildren(node))
          queue.emplace(child);
      }
    }

    void FTMTree_MT::getTreeBranching(std::vector<idNode> &branching,
                                      std::vector<int> &branchingID) {
      std::vector<std::vector<idNode>> nodeBranching;
      this->getTreeBranching(branching, branchingID, nodeBranching);
    }

    std::vector<idNode> FTMTree_MT::getAllRoots() {
      std::vector<idNode> roots;
      for(idNode node = 0; node < this->getNumberOfNodes(); ++node)
        if(this->isRoot(node) and !this->isLeaf(node))
          roots.push_back(node);
      return roots;
    }

    int FTMTree_MT::getNumberOfRoot() {
      int noRoot = 0;
      for(idNode node = 0; node < this->getNumberOfNodes(); ++node)
        if(this->isRoot(node) and !this->isLeaf(node))
          ++noRoot;
      return noRoot;
    }

    int FTMTree_MT::getNumberOfChildren(idNode nodeId) {
      return this->getNode(nodeId)->getNumberOfDownSuperArcs();
    }

    int FTMTree_MT::getTreeDepth() {
      int maxDepth = 0;
      std::queue<std::tuple<idNode, int>> queue;
      queue.push(std::make_tuple(this->getRoot(), 0));
      while(!queue.empty()) {
        auto tup = queue.front();
        queue.pop();
        idNode node = std::get<0>(tup);
        int depth = std::get<1>(tup);
        maxDepth = std::max(maxDepth, depth);
        for(idNode child : this->getChildren(node))
          queue.push(std::make_tuple(child, depth + 1));
      }
      return maxDepth;
    }

    int FTMTree_MT::getNodeLevel(idNode nodeId) {
      int level = 0;
      auto root = this->getRoot();
      int noRoot = this->getNumberOfRoot();
      if(noRoot != 1) {
        std::stringstream ss;
        ss << "problem, there is " << noRoot << " root(s)";
        printErr(ss.str());
        this->printTree2();
        this->printTree();
      }
      if(this->isNodeAlone(nodeId))
        return 0;
      while(nodeId != root) {
        nodeId = this->getParentSafe(nodeId);
        ++level;
      }
      return level;
    }

    std::vector<int> FTMTree_MT::getAllNodeLevel() {
      std::vector<int> allNodeLevel(this->getNumberOfNodes());
      std::queue<std::tuple<idNode, int>> queue;
      queue.push(std::make_tuple(this->getRoot(), 0));
      while(!queue.empty()) {
        auto tup = queue.front();
        queue.pop();
        idNode node = std::get<0>(tup);
        int level = std::get<1>(tup);
        allNodeLevel[node] = level;
        for(idNode child : this->getChildren(node))
          queue.push(std::make_tuple(child, level + 1));
      }
      return allNodeLevel;
    }

    std::vector<std::vector<idNode>> FTMTree_MT::getLevelToNode() {
      std::vector<int> allNodeLevel = this->getAllNodeLevel();
      int maxLevel = *max_element(allNodeLevel.begin(), allNodeLevel.end());
      std::vector<std::vector<idNode>> levelToNode(maxLevel + 1);
      for(unsigned int i = 0; i < allNodeLevel.size(); ++i) {
        levelToNode[allNodeLevel[i]].push_back(i);
      }
      return levelToNode;
    }

    std::vector<idNode>
      FTMTree_MT::getBranchSubtree(std::vector<idNode> &branching,
                                   idNode branchRoot) {
      std::vector<idNode> branchSubtree;
      std::queue<idNode> queue;
      queue.push(branchRoot);
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();

        if(branching[node] != branchRoot
           and this->getParentSafe(node) == branchRoot and node != branchRoot)
          continue;

        branchSubtree.push_back(node);
        for(auto child : this->getChildren(node))
          queue.push(child);
      }
      return branchSubtree;
    }

    // --------------------
    // Persistence
    // --------------------
    std::vector<std::vector<idNode>>
      FTMTree_MT::getMultiPersOriginsVectorFromTree() {
      std::vector<std::vector<idNode>> treeMultiPers(this->getNumberOfNodes());
      for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i)
        if(this->isLeaf(i) and this->isNodeOriginDefined(i)
           and not this->isNodeAlone(i)
           and this->getNode(this->getNode(i)->getOrigin())->getOrigin()
                 != (int)i)
          treeMultiPers[this->getNode(i)->getOrigin()].push_back(i);
      return treeMultiPers;
    }

    // --------------------
    // Set
    // --------------------
    void FTMTree_MT::setParent(idNode nodeId, idNode newParentNodeId) {
      this->deleteParent(nodeId);
      this->makeSuperArc(nodeId, newParentNodeId);
    }

    // --------------------
    // Delete
    // --------------------
    // Delete node by keeping subtree
    void FTMTree_MT::deleteNode(idNode nodeId) {
      if(this->isRoot(nodeId) and !this->isLeaf(nodeId))
        printErr("deletion of root!");

      idNode parentNodeId = this->getParentSafe(nodeId);
      if(!this->isRoot(nodeId)) {
        // Delete down arc from parent node
        // _ Nodes in trees should have only one parent
        idSuperArc nodeArcId = this->getNode(nodeId)->getUpSuperArcId(0);
        this->getNode(parentNodeId)->removeDownSuperArc(nodeArcId);
      }
      // Delete up arc from child nodes
      for(idSuperArc i = 0;
          i < this->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i) {
        idSuperArc arcId = this->getNode(nodeId)->getDownSuperArcId(i);
        idNode childNodeId = this->getSuperArc(arcId)->getDownNodeId();
        this->getNode(childNodeId)->removeUpSuperArc(arcId);
        if(!this->isRoot(nodeId))
          this->makeSuperArc(childNodeId, parentNodeId);
      }
      // Reset deleted node
      this->getNode(nodeId)->clearDownSuperArcs();
      this->getNode(nodeId)->clearUpSuperArcs();
    }

    void FTMTree_MT::deleteIthUpArc(idNode nodeId, int arcIth) {
      idSuperArc nodeArcId = this->getNode(nodeId)->getUpSuperArcId(arcIth);
      // Delete down arc from old parent
      idNode parentNodeId = this->getSuperArc(nodeArcId)->getUpNodeId();
      this->getNode(parentNodeId)->removeDownSuperArc(nodeArcId);
      // Delete up arc from node
      this->getNode(nodeId)->removeUpSuperArc(nodeArcId);
    }

    void FTMTree_MT::deleteParent(idNode nodeId) {
      if(!this->isRoot(nodeId)) {
        // _ Nodes in trees should have only one parent
        this->deleteIthUpArc(nodeId, 0);
      }
    }

    void FTMTree_MT::deleteSubtree(idNode nodeId) {
      std::queue<idNode> queue;
      queue.push(nodeId);
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        auto children = this->getChildren(node);
        for(auto child : children)
          queue.push(child);
        this->deleteNode(node);
      }
    }

    // --------------------
    // Create/Delete/Modify Tree
    // --------------------
    void FTMTree_MT::copyMergeTreeStructure(FTMTree_MT *tree) {
      // Add Nodes
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        this->makeNode(i);
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
        this->getNode(i)->setOrigin(tree->getNode(i)->getOrigin());

      // Add Arcs
      for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i) {
        for(unsigned int j = 0;
            j < tree->getNode(i)->getNumberOfDownSuperArcs(); ++j) {
          auto superArcId = tree->getNode(i)->getDownSuperArcId(j);
          this->makeSuperArc(tree->getSuperArc(superArcId)->getDownNodeId(), i);
        }
      }
    }

    // --------------------
    // Utils
    // --------------------
    void FTMTree_MT::printNodeSS(idNode node, std::stringstream &ss) {
      ss << "(" << node << ") \\ ";
      for(auto child : this->getChildren(node))
        ss << "+" << child << " ";

      if(!this->isRoot(node))
        ss << " / +" << this->getParentSafe(node);
      ss << std::endl;
    }

    std::stringstream FTMTree_MT::printTree(bool doPrint) {
      std::stringstream ss;
      std::vector<idNode> allRoots = this->getAllRoots();
      if(allRoots.size() != 1)
        ss << allRoots.size() << " roots" << std::endl;
      for(unsigned int i = 0; i < allRoots.size(); ++i) {
        if(allRoots.size() != 1)
          ss << i << " _ ";
        ss << "Nodes----------" << std::endl;
        std::queue<idNode> queue;
        queue.push(allRoots[i]);
        while(!queue.empty()) {
          idNode node = queue.front();
          queue.pop();

          printNodeSS(node, ss);

          for(auto child : this->getChildren(node))
            queue.push(child);
        }
      }
      if(allRoots.size() == 0)
        for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i)
          // if(not tree->isNodeAlone(i))
          printNodeSS(i, ss);

      if(doPrint)
        printMsg(ss.str());
      return ss;
    }

    void FTMTree_MT::printTreeStats() {
      auto noNodesT = this->getNumberOfNodes();
      auto noNodes = this->getRealNumberOfNodes();
      std::stringstream ss;
      ss << "tree [node: " << noNodes << " / " << noNodesT;
      ss << ", depth: " << this->getTreeDepth() << "]";
      printMsg(ss.str());
    }

    std::stringstream
      FTMTree_MT::printMultiPersOriginsVectorFromTree(bool doPrint) {
      std::stringstream ss;
      auto vec = this->getMultiPersOriginsVectorFromTree();
      for(unsigned int i = 0; i < vec.size(); ++i)
        if(vec[i].size() != 0) {
          ss << i << " : ";
          for(auto t : vec[i])
            ss << t << " ";
          ss << std::endl;
        }

      if(doPrint)
        printMsg(ss.str());
      return ss;
    }

  } // namespace ftm
} // namespace ttk
