#include <CTA_contourtree.h>
#include <algorithm>
#include <cmath>
#include <float.h>

using ttk::cta::ContourTree;

ContourTree::ContourTree(float *scalars,
                         int *regionSizes,
                         int *segmentationIds,
                         long long *topology,
                         size_t nVertices,
                         size_t nEdges,
                         std::vector<std::vector<int>> regions) {

  nodes = std::vector<std::shared_ptr<ttk::cta::CTNode>>();
  arcs = std::vector<std::shared_ptr<ttk::cta::CTEdge>>();

  float minVal = FLT_MAX;
  float maxVal = -FLT_MAX;

  for(size_t i = 0; i < nVertices; i++) {
    std::shared_ptr<ttk::cta::CTNode> node(new ttk::cta::CTNode);
    node->scalarValue = scalars[i];
    node->edgeList = std::vector<int>();
    node->branchID = -1;
    nodes.push_back(node);
    if(node->scalarValue > maxVal)
      maxVal = node->scalarValue;
    if(node->scalarValue < minVal)
      minVal = node->scalarValue;
  }
  int j = 0;
  for(size_t i = 0; i < nEdges; i++) {
    std::shared_ptr<ttk::cta::CTEdge> edge(new ttk::cta::CTEdge);
    edge->area = regionSizes[i];
    edge->segId = segmentationIds[i];
    if(!regions.empty())
      edge->region = regions[segmentationIds[i]];
    std::sort(edge->region.begin(), edge->region.end());
    // if(!regions.empty()) std::cout << regionSizes[i] << " " <<
    // regions[segmentationIds[i]].size() << std::endl;
    edge->node1Idx = topology[j + 0];
    nodes[topology[j + 0]]->edgeList.push_back(i);
    edge->node2Idx = topology[j + 1];
    nodes[topology[j + 1]]->edgeList.push_back(i);
    edge->scalardistance = std::abs(nodes[edge->node1Idx]->scalarValue
                                    - nodes[edge->node2Idx]->scalarValue);
    edge->volume = edge->area * edge->scalardistance;
    arcs.push_back(edge);
    j += 2;
  }

  for(const std::shared_ptr<ttk::cta::CTNode> &node : nodes) {
    if(node->edgeList.size() == 0)
      std::cout << "wtf?\n" << std::flush;
    if(node->edgeList.size() > 1)
      node->type = saddleNode;
    else {
      std::shared_ptr<ttk::cta::CTEdge> edge = arcs[node->edgeList[0]];
      std::shared_ptr<ttk::cta::CTNode> neighbor = nodes[edge->node1Idx] == node
                                                     ? nodes[edge->node2Idx]
                                                     : nodes[edge->node1Idx];
      if(neighbor->scalarValue > node->scalarValue)
        node->type = minNode;
      else
        node->type = maxNode;
    }
  }

  binary = true;
  for(const std::shared_ptr<ttk::cta::CTNode> &node : nodes) {
    if(node->edgeList.size() > 3)
      binary = false;
  }
}

ContourTree::~ContourTree() {
}

std::shared_ptr<ttk::cta::Tree> ContourTree::computeRootedTree(
  const std::shared_ptr<ttk::cta::CTNode> &node,
  const std::shared_ptr<ttk::cta::CTEdge> &parent,
  int &id) {

  // initialize tree
  std::shared_ptr<Tree> t(new Tree);

  // set id and increment for later calls
  t->id = id;
  id++;

  // set type/label
  t->type = node->type;

  // set height and size to 0/1. For inner nodes this will be updated while
  // traversing the children
  t->size = 1;
  t->height = 0;

  // compute number of children (depends on whether current node is the root or
  // not)
  if(parent == nullptr)
    t->children = std::vector<std::shared_ptr<Tree>>(node->edgeList.size());
  else
    t->children = std::vector<std::shared_ptr<Tree>>(node->edgeList.size() - 1);

  bool parentVisited = false;

  // add neighbors to children
  for(size_t i = 0; i < node->edgeList.size(); i++) {

    std::shared_ptr<ttk::cta::CTEdge> edge = arcs[node->edgeList[i]];
    if(edge == parent) {
      parentVisited = true;
      continue;
    }

    std::shared_ptr<ttk::cta::CTNode> child = nodes[edge->node1Idx] == node
                                                ? nodes[edge->node2Idx]
                                                : nodes[edge->node1Idx];

    std::shared_ptr<Tree> childTree = computeRootedTree(child, edge, id);

    t->children[i - (parentVisited ? 1 : 0)] = childTree;

    t->size += childTree->size;

    if(childTree->height + 1 > t->height)
      t->height = childTree->height + 1;
  }

  // get Persistence of parent edge and compute volume
  if(parent == nullptr) {
    t->scalardistanceParent = 0.0001;
    t->volume = 0.0001;
  } else {
    t->scalardistanceParent = parent->scalardistance;
    t->volume = t->scalardistanceParent * parent->area;
  }

  return t;
}

std::shared_ptr<ttk::cta::BinaryTree> ContourTree::computeRootedTree_binary(
  const std::shared_ptr<ttk::cta::CTNode> &node,
  const std::shared_ptr<ttk::cta::CTEdge> &parent,
  int &id) {

  // initialize tree
  std::shared_ptr<BinaryTree> t(new BinaryTree);

  // set id and increment for later calls
  t->id = id;
  id++;

  // set type/label
  t->type = node->type;
  std::shared_ptr<ttk::cta::CTEdge> edge = arcs[node->edgeList[0]];
  int nodeIdx = nodes[edge->node1Idx] == node ? edge->node1Idx : edge->node2Idx;
  t->nodeRefs = std::vector<std::pair<int, int>>();
  t->nodeRefs.push_back(std::make_pair(-1, nodeIdx));
  t->arcRefs = std::vector<std::pair<int, int>>();
  // if(parent != nullptr)
  // t->arcRefs.push_back(std::make_pair(-1,parent->segId));
  if(parent != nullptr) {
    int arcRef = arcs[node->edgeList[0]] == parent   ? node->edgeList[0]
                 : arcs[node->edgeList[1]] == parent ? node->edgeList[1]
                                                     : node->edgeList[2];
    t->arcRefs.push_back(std::make_pair(-1, arcRef));
  }

  // set height and size to 0/1. For inner nodes this will be updated while
  // traversing the children
  t->size = 1;
  t->height = 0;

  // children at first into vector
  std::vector<std::shared_ptr<BinaryTree>> children;

  // add neighbors to children
  for(size_t i = 0; i < node->edgeList.size(); i++) {

    edge = arcs[node->edgeList[i]];

    if(edge != parent) {

      std::shared_ptr<ttk::cta::CTNode> child = nodes[edge->node1Idx] == node
                                                  ? nodes[edge->node2Idx]
                                                  : nodes[edge->node1Idx];

      std::shared_ptr<BinaryTree> childTree
        = computeRootedTree_binary(child, edge, id);

      children.push_back(childTree);

      t->size += childTree->size;

      if(childTree->height + 1 > t->height)
        t->height = childTree->height + 1;
    }
  }

  // children from vector to binary tree
  t->child1 = children.size() > 0 ? children[0] : nullptr;
  t->child2 = children.size() > 1 ? children[1] : nullptr;

  t->freq = 1;

  t->scalarValue = node->scalarValue;

  // get Persistence of parent edge and compute volume
  if(parent == nullptr) {
    t->scalardistanceParent = 10000;
    t->area = 10000;
    t->volume = 10000;
    t->region = std::vector<int>(1, -1);
  } else {
    t->scalardistanceParent = parent->scalardistance;
    t->area = parent->area;
    t->volume = t->scalardistanceParent * parent->area;
    t->region = parent->region;
  }

  return t;
}

std::shared_ptr<ttk::cta::BinaryTree> ContourTree::rootAtMax() {

  // get global maximum node to build rooted tree from there
  float maxVal = -FLT_MAX;
  std::shared_ptr<ttk::cta::CTNode> globalMax = nullptr;

  for(const std::shared_ptr<ttk::cta::CTNode> &node : nodes) {
    if(node->scalarValue > maxVal) {
      globalMax = node;
      maxVal = node->scalarValue;
    }
  }

  int id = 1;

  return computeRootedTree_binary(globalMax, nullptr, id);
}

std::shared_ptr<ttk::cta::BinaryTree>
  ContourTree::rootAtNode(const std::shared_ptr<ttk::cta::CTNode> &root) {

  int id = 1;

  // rootedTree = computeRootedTree(root,-1, id);
  return computeRootedTree_binary(root, nullptr, id);
}

bool ContourTree::isBinary() {
  return binary;
}

void ContourTree::computeBranches() {

  // find global minimum
  int minIdx = 0;
  for(size_t i = 1; i < nodes.size(); i++) {
    if(nodes[minIdx]->scalarValue > nodes[i]->scalarValue)
      minIdx = i;
  }

  // find path to global max
  int nextIdx = arcs[nodes[minIdx]->edgeList[0]]->node1Idx == minIdx
                  ? arcs[nodes[minIdx]->edgeList[0]]->node2Idx
                  : arcs[nodes[minIdx]->edgeList[0]]->node1Idx;
  std::vector<int> maxPath_ = pathToMax(nextIdx, minIdx).second;
  std::vector<int> maxPath;
  maxPath.push_back(minIdx);
  maxPath.insert(maxPath.end(), maxPath_.begin(), maxPath_.end());

  int currID = 0;
  nodes[minIdx]->branchID = 0;
  std::stack<std::vector<int>> q;
  q.push(maxPath);

  while(!q.empty()) {

    std::vector<int> path = q.top();
    q.pop();

    for(size_t i = 1; i < path.size() - 1; i++) {

      int idx = path[i];

      for(int cE : nodes[idx]->edgeList) {

        int cIdx
          = arcs[cE]->node1Idx == idx ? arcs[cE]->node2Idx : arcs[cE]->node1Idx;
        if(cIdx == path[i - 1])
          continue;
        if(cIdx == path[i + 1])
          continue;

        if(nodes[cIdx]->scalarValue > nodes[idx]->scalarValue) {
          std::vector<int> newPath_ = pathToMax(cIdx, idx).second;
          std::vector<int> newPath;
          newPath.push_back(idx);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        } else {
          std::vector<int> newPath_ = pathToMin(cIdx, idx).second;
          std::vector<int> newPath;
          newPath.push_back(idx);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        }
      }

      nodes[idx]->branchID = currID;
    }

    // nodes[path.front()]->branchID = currID;
    nodes[path.back()]->branchID = currID;
    currID++;
  }
}

std::pair<float, std::vector<int>> ContourTree::pathToMax(int root,
                                                          int parent) {

  std::vector<int> path;
  path.push_back(root);

  if(nodes[root]->edgeList.size() == 1) {
    return std::make_pair(nodes[root]->scalarValue, path);
  }

  std::vector<int> bestPath;
  float bestVal = -FLT_MAX;

  for(int cE : nodes[root]->edgeList) {

    int nextIdx
      = arcs[cE]->node1Idx == root ? arcs[cE]->node2Idx : arcs[cE]->node1Idx;
    if(parent == nextIdx)
      continue;
    if(nodes[nextIdx]->scalarValue < nodes[root]->scalarValue)
      continue;

    auto p = pathToMax(nextIdx, root);
    if(p.first > bestVal) {
      bestVal = p.first;
      bestPath = p.second;
    }
  }

  path.insert(path.end(), bestPath.begin(), bestPath.end());

  return std::make_pair(bestVal, path);
}

std::pair<float, std::vector<int>> ContourTree::pathToMin(int root,
                                                          int parent) {

  std::vector<int> path;
  path.push_back(root);

  if(nodes[root]->edgeList.size() == 1) {
    return std::make_pair(nodes[root]->scalarValue, path);
  }

  std::vector<int> bestPath;
  float bestVal = FLT_MAX;

  for(int cE : nodes[root]->edgeList) {

    int nextIdx
      = arcs[cE]->node1Idx == root ? arcs[cE]->node2Idx : arcs[cE]->node1Idx;
    if(parent == nextIdx)
      continue;
    if(nodes[nextIdx]->scalarValue > nodes[root]->scalarValue)
      continue;

    auto p = pathToMin(nextIdx, root);
    if(p.first < bestVal) {
      bestVal = p.first;
      bestPath = p.second;
    }
  }

  path.insert(path.end(), bestPath.begin(), bestPath.end());

  return std::make_pair(bestVal, path);
}

std::pair<std::vector<std::shared_ptr<ttk::cta::CTNode>>,
          std::vector<std::shared_ptr<ttk::cta::CTEdge>>>
  ContourTree::getGraph() {
  return std::make_pair(nodes, arcs);
}
