#include <FTMTreeUtils.h>

using namespace ttk;
using namespace ftm;

// --------------------
// Is
// --------------------

bool isNodeOriginDefined(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  unsigned int origin = (unsigned int)tree->getNode(nodeId)->getOrigin();
  return origin != ftm::nullNodes and origin < tree->getNumberOfNodes()
         and origin >= 0;
  ;
}

bool isRoot(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  /*if(nodeId >= tree->getNumberOfNodes())
    std::cout << nodeId << " _____ " << tree->getNumberOfNodes() << std::endl;*/
  return tree->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
}

bool isLeaf(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return tree->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
}

bool isNodeAlone(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return isRoot(tree, nodeId) and isLeaf(tree, nodeId);
}

bool isFullMerge(ftm::FTMTree_MT *tree) {
  ftm::idNode treeRoot = getRoot(tree);
  return (unsigned int)tree->getNode(treeRoot)->getOrigin() == treeRoot;
}

bool isBranchOrigin(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return getParent(tree, tree->getNode(nodeId)->getOrigin()) != nodeId;
}

// --------------------
// Get
// --------------------

ftm::idNode getRoot(ftm::FTMTree_MT *tree) {
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node)
    if(isRoot(tree, node) and !isLeaf(tree, node))
      return node;
  return ftm::nullNodes;
}

ftm::idNode getParent(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  if(!isRoot(tree, nodeId)) {
    // _ Nodes in merge trees should have only one parent
    ftm::idSuperArc arcId = tree->getNode(nodeId)->getUpSuperArcId(0);
    ftm::idNode parentNodeId = tree->getSuperArc(arcId)->getUpNodeId();
    return parentNodeId;
  }
  return nodeId;
}

std::vector<ftm::idNode> getChildren(ftm::FTMTree_MT *tree,
                                     ftm::idNode nodeId) {
  std::vector<ftm::idNode> childrens;
  for(ftm::idSuperArc i = 0;
      i < tree->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i) {
    ftm::idSuperArc arcId = tree->getNode(nodeId)->getDownSuperArcId(i);
    childrens.push_back(tree->getSuperArc(arcId)->getDownNodeId());
  }
  return childrens;
}

std::vector<ttk::ftm::idNode> getLeaves(ttk::ftm::FTMTree_MT *tree) {
  std::vector<ttk::ftm::idNode> treeLeaves;
  for(ttk::ftm::idNode i = 0; i < tree->getNumberOfNodes(); ++i) {
    if(isLeaf(tree, i) and !isRoot(tree, i))
      treeLeaves.push_back(i);
  }
  return treeLeaves;
}

int getNumberOfLeaves(ftm::FTMTree_MT *tree) {
  auto leaves = getLeaves(tree);
  return leaves.size();
}

int getNumberOfNodeAlone(ftm::FTMTree_MT *tree) {
  int cpt = 0;
  for(ftm::idNode i = 0; i < tree->getNumberOfNodes(); ++i)
    cpt += isNodeAlone(tree, i) ? 1 : 0;
  return cpt;
}

int getRealNumberOfNodes(ftm::FTMTree_MT *tree) {
  return tree->getNumberOfNodes() - getNumberOfNodeAlone(tree);
}

std::tuple<std::vector<ftm::idNode>, std::vector<ftm::idNode>>
  getBranchOriginsFromThisBranch(FTMTree_MT *tree, ftm::idNode node) {
  std::vector<ftm::idNode> branchOrigins, nonBranchOrigins;

  ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
  ftm::idNode nodeParent = getParent(tree, nodeOrigin);
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
                      std::vector<std::vector<ftm::idNode>> &nodeBranching) {
  branching = std::vector<idNode>(tree->getNumberOfNodes());
  branchingID = std::vector<int>(tree->getNumberOfNodes(), -1);
  nodeBranching
    = std::vector<std::vector<ftm::idNode>>(tree->getNumberOfNodes());
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
        std::cout << getRoot(tree) << " _ " << parentNodeOrigin << std::endl;
        std::cout << node << " _ " << nodeOrigin << std::endl;
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
  std::vector<std::vector<ftm::idNode>> nodeBranching;
  getTreeBranching(tree, branching, branchingID, nodeBranching);
}
