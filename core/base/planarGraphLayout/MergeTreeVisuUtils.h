/// \ingroup base
/// \class ttk::MergeTreeVisuUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
///
/// Temporary file until mergeTreeClustering is in TTK

#ifndef _MERGETREEVISUUTILS_H
#define _MERGETREEVISUUTILS_H

#pragma once

#include <FTMTree.h>
#include <FTMTreePP.h>

#include <ttkUtils.h>

#include <vtkCellData.h>

using namespace ttk;
using namespace ftm;

bool isNodeOriginDefined(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  unsigned int origin = (unsigned int)tree->getNode(nodeId)->getOrigin();
  return origin != ftm::nullNodes and origin < tree->getNumberOfNodes()
         and origin >= 0;
  ;
}

bool isRoot(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  /*if(nodeId >= tree->getNumberOfNodes()){
    std::stringstream ss;
    ss << nodeId << " _____ " << tree->getNumberOfNodes();
    //printMsg(ss.str());
  }*/
  return tree->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
}

bool isLeaf(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return tree->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
}

bool isNodeAlone(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return isRoot(tree, nodeId) and isLeaf(tree, nodeId);
}

ftm::idNode getRoot(ftm::FTMTree_MT *tree) {
  for(ftm::idNode node = 0; node < tree->getNumberOfNodes(); ++node)
    if(isRoot(tree, node) and !isLeaf(tree, node))
      return node;
  return ftm::nullNodes;
}

bool isFullMerge(ftm::FTMTree_MT *tree) {
  ftm::idNode treeRoot = getRoot(tree);
  return (unsigned int)tree->getNode(treeRoot)->getOrigin() == treeRoot;
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

bool isBranchOrigin(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return getParent(tree, tree->getNode(nodeId)->getOrigin()) != nodeId;
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

template <class dataType>
bool isJoinTree2(ftm::FTMTree_MT *tree) {
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

template <class dataType>
std::tuple<dataType, dataType> getBirthDeath2(ftm::FTMTree_MT *tree,
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
dataType getBirth2(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  return std::get<0>(getBirthDeath2<dataType>(tree, nodeId));
}

template <class dataType>
dataType getNodePersistence2(ftm::FTMTree_MT *tree, ftm::idNode nodeId) {
  std::tuple<dataType, dataType> birthDeath
    = getBirthDeath2<dataType>(tree, nodeId);
  return std::get<1>(birthDeath) - std::get<0>(birthDeath);
}

template <class dataType>
void getPersistencePairs2(
  ftm::FTMTree_MT *tree,
  // std::vector<std::tuple<SimplexId, SimplexId>> &pairs){
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> &pairs) {
  /*ttk::ftm::FTMTreePPCompute pairsCompute;
  pairsCompute.computePersistencePairs<dataType>(tree, pairs);*/
  ttk::ftm::FTMTreePP pairsCompute;
  pairsCompute.setCustomTree(tree);
  pairsCompute.computePersistencePairs<dataType>(
    pairs, isJoinTree2<dataType>(tree));
}

template <class dataType>
std::vector<std::tuple<SimplexId, SimplexId, dataType>>
  computePersistencePairs2(ftm::FTMTree_MT *tree) {
  // std::vector<std::tuple<SimplexId, SimplexId>> pairs;
  std::vector<std::tuple<SimplexId, SimplexId, dataType>> pairs;
  getPersistencePairs2<dataType>(tree, pairs);
  for(auto pair : pairs) {
    if(tree->getNode(std::get<0>(pair))->getOrigin() < std::get<0>(pair)
       and tree->getNode(std::get<0>(pair))->getOrigin() >= 0)
      tree->getNode(tree->getNode(std::get<0>(pair))->getOrigin())
        ->setOrigin(std::get<1>(pair));

    tree->getNode(std::get<0>(pair))->setOrigin(std::get<1>(pair));
    tree->getNode(std::get<1>(pair))->setOrigin(std::get<0>(pair));
  }
  return pairs;
}

template <class dataType>
bool isImportantPair(ftm::FTMTree_MT *tree,
                     ftm::idNode nodeId,
                     double threshold) {
  dataType rootPers = getNodePersistence2<dataType>(tree, getRoot(tree));
  if(threshold > 1)
    threshold /= 100.0;
  threshold = rootPers * threshold;
  return getNodePersistence2<dataType>(tree, nodeId) > threshold;
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
        //printMsg( getRoot(tree) << " _ " << parentNodeOrigin );
        //printMsg( node << " _ " << nodeOrigin );
        //printMsg( "getTreeBranching" );
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

template <typename type>
static type myAbs2(const type var) {
  return (var >= 0) ? var : -var;
}

template <class dataType>
bool isEqual2(dataType first, dataType two, double eps = 1e-6) {
  return myAbs2<dataType>(first - two)
         < eps * std::max(myAbs2<dataType>(first), myAbs2<dataType>(two));
}

struct MergeTree2 {
  ftm::Scalars scalars;
  ftm::Params params;
  ftm::FTMTree_MT tree;
  MergeTree2(ftm::Scalars scalarsT, ftm::Params paramsT)
    : scalars(scalarsT), params(paramsT), tree(&params, &scalars, Join_Split) {
  }
};

MergeTree2 makeTree2(vtkUnstructuredGrid *treeNodes,
                     vtkUnstructuredGrid *treeArcs) {
  // Init Scalars
  Scalars scalars;
  vtkSmartPointer<vtkDataArray> nodesScalar
    = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
  scalars.size = nodesScalar->GetNumberOfTuples();
  scalars.values = ttkUtils::GetVoidPointer(nodesScalar);

  // Init Tree
  Params params;
  params.treeType = Join_Split;
  MergeTree2 mergeTree(scalars, params);
  mergeTree.tree.makeAlloc();

  // Add Nodes
  vtkSmartPointer<vtkDataArray> nodesId
    = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
  vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
  for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
    mergeTree.tree.makeNode(i);
  }

  // Add Arcs
  vtkSmartPointer<vtkDataArray> arcsUp
    = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
  vtkSmartPointer<vtkDataArray> arcsDown
    = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
  vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
  std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
  for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
    double downId = arcsDown->GetTuple1(i);
    double upId = arcsUp->GetTuple1(i);
    auto it = added_arcs.find(std::make_tuple(downId, upId));
    if(it == added_arcs.end()) { // arc not added yet
      mergeTree.tree.makeSuperArc(downId, upId); // (down, Up)
      added_arcs.insert(std::make_tuple(downId, upId));
    }
  }

  // Manage inconsistent arcs
  // manageInconsistentArcsMultiParent(tree);;

  // Remove self link
  // removeSelfLink(tree);

  // tree->printTree2();

  return mergeTree;
}

#endif
