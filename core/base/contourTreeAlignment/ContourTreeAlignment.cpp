#include "ContourTreeAlignment.h"

/*
=====================================================================================================================
 branch decomposition
=====================================================================================================================
*/

void ttk::ContourTreeAlignment::computeBranches() {

  // find global minimum
  std::shared_ptr<ttk::cta::AlignmentNode> minNode = nodes[0];
  for(size_t i = 1; i < nodes.size(); i++) {
    if(minNode->scalarValue > nodes[i]->scalarValue)
      minNode = nodes[i];
  }

  // find path to global max
  std::shared_ptr<ttk::cta::AlignmentNode> nextNode
    = minNode->edgeList[0]->node1.lock() == minNode
        ? minNode->edgeList[0]->node2.lock()
        : minNode->edgeList[0]->node1.lock();
  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> maxPath_
    = pathToMax(nextNode, minNode).second;
  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> maxPath;
  maxPath.push_back(minNode);
  maxPath.insert(maxPath.end(), maxPath_.begin(), maxPath_.end());

  int currID = 0;
  minNode->branchID = 0;
  std::stack<std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>> q;
  q.push(maxPath);

  while(!q.empty()) {

    std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> path = q.top();
    q.pop();

    for(size_t i = 1; i < path.size() - 1; i++) {

      std::shared_ptr<ttk::cta::AlignmentNode> currNode = path[i];

      for(const auto &cE : currNode->edgeList) {

        std::shared_ptr<ttk::cta::AlignmentNode> cN
          = (cE->node1.lock()) == currNode ? cE->node2.lock()
                                           : cE->node1.lock();
        if(cN == path[i - 1])
          continue;
        if(cN == path[i + 1])
          continue;

        if(cN->scalarValue > currNode->scalarValue) {
          std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath_
            = pathToMax(cN, currNode).second;
          std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath;
          newPath.push_back(currNode);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        } else {
          std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath_
            = pathToMin(cN, currNode).second;
          std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath;
          newPath.push_back(currNode);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        }
      }

      currNode->branchID = currID;
    }

    for(const auto &cE : path.back()->edgeList) {

      std::shared_ptr<ttk::cta::AlignmentNode> cN
        = cE->node1.lock() == path.back() ? cE->node2.lock() : cE->node1.lock();
      if(cN == path[path.size() - 2])
        continue;

      if(cN->scalarValue > path.back()->scalarValue) {
        std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath_
          = pathToMax(cN, path.back()).second;
        std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath;
        newPath.push_back(path.back());
        newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
        q.push(newPath);
      } else {
        std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath_
          = pathToMin(cN, path.back()).second;
        std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> newPath;
        newPath.push_back(path.back());
        newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
        q.push(newPath);
      }
    }

    path.back()->branchID = currID;
    currID++;
  }
}

std::pair<float, std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>>
  ttk::ContourTreeAlignment::pathToMax(
    const std::shared_ptr<ttk::cta::AlignmentNode> &root,
    const std::shared_ptr<ttk::cta::AlignmentNode> &parent) {

  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> path;
  path.push_back(root);

  if(root->edgeList.size() == 1) {
    return std::make_pair(root->scalarValue, path);
  }

  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> bestPath;
  float bestVal = -FLT_MAX;

  for(const auto &cE : root->edgeList) {

    std::shared_ptr<ttk::cta::AlignmentNode> nextNode
      = cE->node1.lock() == root ? cE->node2.lock() : cE->node1.lock();
    if(parent == nextNode)
      continue;
    if(nextNode->scalarValue < root->scalarValue)
      continue;

    auto p = pathToMax(nextNode, root);
    if(p.first > bestVal) {
      bestVal = p.first;
      bestPath = p.second;
    }
  }

  path.insert(path.end(), bestPath.begin(), bestPath.end());

  return std::make_pair(bestVal, path);
}

std::pair<float, std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>>
  ttk::ContourTreeAlignment::pathToMin(
    const std::shared_ptr<ttk::cta::AlignmentNode> &root,
    const std::shared_ptr<ttk::cta::AlignmentNode> &parent) {

  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> path;
  path.push_back(root);

  if(root->edgeList.size() == 1) {
    return std::make_pair(root->scalarValue, path);
  }

  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> bestPath;
  float bestVal = FLT_MAX;

  for(const auto &cE : root->edgeList) {

    std::shared_ptr<ttk::cta::AlignmentNode> nextNode
      = cE->node1.lock() == root ? cE->node2.lock() : cE->node1.lock();
    if(parent == nextNode)
      continue;
    if(nextNode->scalarValue > root->scalarValue)
      continue;

    auto p = pathToMin(nextNode, root);
    if(p.first < bestVal) {
      bestVal = p.first;
      bestPath = p.second;
    }
  }

  path.insert(path.end(), bestPath.begin(), bestPath.end());

  return std::make_pair(bestVal, path);
}

/*
=====================================================================================================================
 iterated aligning
=====================================================================================================================
*/

using ttk::cta::ContourTree;

bool ttk::ContourTreeAlignment::initialize(
  const std::shared_ptr<ContourTree> &ct) {

  contourtrees.push_back(ct);

  // cancel computation if tree is not binary
  if(!ct->isBinary())
    return false;

  // compute initial Alignment from initial contourtree

  std::shared_ptr<ttk::cta::BinaryTree> t = ct->rootAtMax();

  std::queue<std::pair<std::shared_ptr<ttk::cta::BinaryTree>,
                       std::shared_ptr<ttk::cta::AlignmentNode>>>
    q;

  std::shared_ptr<ttk::cta::AlignmentNode> currNode(
    new ttk::cta::AlignmentNode());
  std::shared_ptr<ttk::cta::BinaryTree> currTree;

  currNode->freq = 1;
  currNode->type = t->type;
  currNode->branchID = -1;
  currNode->scalarValue = t->scalarValue;

  currNode->nodeRefs = std::vector<std::pair<int, int>>();
  currNode->nodeRefs.push_back(std::make_pair(0, t->nodeRefs[0].second));

  nodes.push_back(currNode);

  q.push(std::make_pair(t, currNode));

  while(!q.empty()) {

    currTree = q.front().first;
    currNode = q.front().second;
    q.pop();

    if(currTree->child1 != nullptr) {

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());
      childNode->freq = 1;
      childNode->type = currTree->child1->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child1->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child1->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());
      childEdge->area = currTree->child1->area;
      childEdge->scalardistance = currTree->child1->scalardistanceParent;
      childEdge->volume = currTree->child1->volume;
      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      q.push(std::make_pair(currTree->child1, childNode));
    }

    if(currTree->child2 != nullptr) {

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());
      childNode->freq = 1;
      childNode->type = currTree->child2->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child2->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child2->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());
      childEdge->area = currTree->child2->area;
      childEdge->scalardistance = currTree->child2->scalardistanceParent;
      childEdge->volume = currTree->child2->volume;
      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      q.push(std::make_pair(currTree->child2, childNode));
    }
  }

  return true;
}

bool ttk::ContourTreeAlignment::initialize_consistentRoot(
  const std::shared_ptr<ContourTree> &ct, int rootIdx) {

  // cancel computation if tree is not binary
  if(!ct->isBinary())
    return false;

  contourtrees.push_back(ct);

  // compute initial Alignment from initial contourtree

  std::shared_ptr<ttk::cta::BinaryTree> t = ct->rootAtMax();

  std::queue<std::pair<std::shared_ptr<ttk::cta::BinaryTree>,
                       std::shared_ptr<ttk::cta::AlignmentNode>>>
    q;

  std::shared_ptr<ttk::cta::AlignmentNode> currNode(
    new ttk::cta::AlignmentNode());
  std::shared_ptr<ttk::cta::BinaryTree> currTree;

  currNode->freq = 1;
  currNode->type = t->type;
  currNode->branchID = -1;
  currNode->scalarValue = t->scalarValue;

  currNode->nodeRefs = std::vector<std::pair<int, int>>();
  currNode->nodeRefs.push_back(std::make_pair(0, t->nodeRefs[0].second));

  nodes.push_back(currNode);

  q.push(std::make_pair(t, currNode));

  while(!q.empty()) {

    currTree = q.front().first;
    currNode = q.front().second;
    q.pop();

    if(currTree->child1 != nullptr) {

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());
      childNode->freq = 1;
      childNode->type = currTree->child1->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child1->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child1->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());
      childEdge->area = currTree->child1->area;
      childEdge->scalardistance = currTree->child1->scalardistanceParent;
      childEdge->volume = currTree->child1->volume;
      childEdge->region = currTree->child1->region;
      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      childEdge->arcRefs.push_back(
        std::make_pair(0, currTree->child1->arcRefs[0].second));

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      q.push(std::make_pair(currTree->child1, childNode));
    }

    if(currTree->child2 != nullptr) {

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());
      childNode->freq = 1;
      childNode->type = currTree->child2->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child2->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child2->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());
      childEdge->area = currTree->child2->area;
      childEdge->scalardistance = currTree->child2->scalardistanceParent;
      childEdge->volume = currTree->child2->volume;
      childEdge->region = currTree->child2->region;
      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      childEdge->arcRefs.push_back(
        std::make_pair(0, currTree->child2->arcRefs[0].second));

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      q.push(std::make_pair(currTree->child2, childNode));
    }
  }

  if(nodes.size() != ct->getGraph().first.size()) {
    printErr("wtf?");
  }

  alignmentRoot = nodes[rootIdx];
  alignmentRootIdx = rootIdx;

  return true;
}

bool ttk::ContourTreeAlignment::alignTree(
  const std::shared_ptr<ContourTree> &ct) {

  contourtrees.push_back(ct);

  // cancel computation if tree not binary
  if(!ct->isBinary())
    return false;

  // compute optimal alignment between current alignment and new tree

  std::shared_ptr<ttk::cta::BinaryTree> t1, t2;

  std::shared_ptr<ttk::cta::AlignmentTree> res = nullptr;
  float resVal = FLT_MAX;

  std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> nodes1 = nodes;
  std::vector<std::shared_ptr<ttk::cta::CTNode>> nodes2 = ct->getGraph().first;

  int i = 0;
  int j = 0;

  for(const auto &node1 : nodes1) {

    t1 = this->rootAtNode(node1);

    j = 0;

    for(const auto &node2 : nodes2) {

      t2 = ct->rootAtNode(node2);

      // compute matching

      if((node1->type == ttk::cta::maxNode && node2->type == ttk::cta::maxNode)
         || (node1->type == ttk::cta::minNode
             && node2->type == ttk::cta::minNode)) {

        std::pair<float, std::shared_ptr<ttk::cta::AlignmentTree>> match
          = getAlignmentBinary(t1, t2);

        if(match.first < resVal) {
          resVal = match.first;
          res = match.second;
        }
      }

      j++;
    }

    i++;
  }

  if(res)
    computeNewAlignmenttree(res);
  else {
    printErr("Alignment computation failed.");
    return false;
  }

  return true;
}

bool ttk::ContourTreeAlignment::alignTree_consistentRoot(
  const std::shared_ptr<ContourTree> &ct) {

  // cancel computation if tree not binary
  if(!ct->isBinary())
    return false;

  contourtrees.push_back(ct);

  // compute optimal alignment between current alignment and new tree

  std::shared_ptr<ttk::cta::BinaryTree> t1, t2;

  std::shared_ptr<ttk::cta::AlignmentTree> res = nullptr;
  float resVal = FLT_MAX;

  std::vector<std::shared_ptr<ttk::cta::CTNode>> nodes2 = ct->getGraph().first;

  int i = 0;

  t1 = this->rootAtNode(alignmentRoot);

  for(const auto &node2 : nodes2) {

    // compute matching

    if((alignmentRoot->type == ttk::cta::maxNode
        && node2->type == ttk::cta::maxNode)
       || (alignmentRoot->type == ttk::cta::minNode
           && node2->type == ttk::cta::minNode)) {

      t2 = ct->rootAtNode(node2);

      std::pair<float, std::shared_ptr<ttk::cta::AlignmentTree>> match
        = getAlignmentBinary(t1, t2);

      if(match.first < resVal) {
        resVal = match.first;
        res = match.second;
      }
    }

    i++;
  }

  if(res)
    computeNewAlignmenttree(res);
  else {
    printErr("Alignment computation failed.");
    return false;
  }

  alignmentVal += resVal;

  alignmentRoot = nodes[0];

  return true;
}

void ttk::ContourTreeAlignment::computeNewAlignmenttree(
  const std::shared_ptr<ttk::cta::AlignmentTree> &res) {

  nodes.clear();
  arcs.clear();

  std::queue<std::tuple<std::shared_ptr<ttk::cta::AlignmentTree>,
                        std::shared_ptr<ttk::cta::AlignmentNode>,
                        std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>,
                        std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>>>
    q;

  std::shared_ptr<ttk::cta::AlignmentNode> currNode;
  std::shared_ptr<ttk::cta::AlignmentTree> currTree;

  currNode
    = std::shared_ptr<ttk::cta::AlignmentNode>(new ttk::cta::AlignmentNode());

  if(res->node1 == nullptr && res->node2 == nullptr) {
    return;
  }

  currNode->freq = (res->node1 == nullptr ? 0 : res->node1->freq)
                   + (res->node2 == nullptr ? 0 : res->node2->freq);
  currNode->type = res->node1 == nullptr ? res->node2->type : res->node1->type;
  currNode->branchID = -1;

  if(alignmenttreeType == ttk::cta::lastMatchedValue) {

    currNode->scalarValue = res->node2 == nullptr ? res->node1->scalarValue
                                                  : res->node2->scalarValue;

    currNode->nodeRefs = std::vector<std::pair<int, int>>();
    if(res->node1 != nullptr)
      currNode->nodeRefs.insert(currNode->nodeRefs.end(),
                                res->node1->nodeRefs.begin(),
                                res->node1->nodeRefs.end());
    if(res->node2 != nullptr)
      currNode->nodeRefs.push_back(std::make_pair(
        (int)contourtrees.size() - 1, res->node2->nodeRefs[0].second));

  } else if(alignmenttreeType == ttk::cta::averageValues) {

    currNode->scalarValue = res->node1 == nullptr ? res->node2->scalarValue
                            : res->node2 == nullptr
                              ? res->node1->scalarValue
                              : (res->node1->scalarValue * res->node1->freq
                                 + res->node2->scalarValue)
                                  / currNode->freq;

  } else {

    if(res->node1 == nullptr)
      currNode->scalarValue = res->node2->scalarValue;
    else if(res->node2 == nullptr)
      currNode->scalarValue = res->node1->scalarValue;
    else {
      std::vector<float> values;
      for(auto p : res->node1->nodeRefs) {
        values.push_back(
          contourtrees[p.first]->getGraph().first[p.second]->scalarValue);
      }
      values.push_back(res->node2->scalarValue);
      std::sort(values.begin(), values.end());
      float newMedian = values[values.size() / 2];
      currNode->scalarValue = newMedian;
    }
  }

  currNode->nodeRefs = std::vector<std::pair<int, int>>();
  if(res->node1 != nullptr)
    currNode->nodeRefs.insert(currNode->nodeRefs.end(),
                              res->node1->nodeRefs.begin(),
                              res->node1->nodeRefs.end());
  if(res->node2 != nullptr)
    currNode->nodeRefs.push_back(std::make_pair(
      (int)contourtrees.size() - 1, res->node2->nodeRefs[0].second));

  nodes.push_back(currNode);

  q.push(std::make_tuple(
    res, currNode, std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>(),
    std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>()));

  while(!q.empty()) {

    currTree = std::get<0>(q.front());
    currNode = std::get<1>(q.front());
    auto openEdgesOld1 = std::get<2>(q.front());
    auto openEdgesNew1 = std::get<3>(q.front());
    auto openEdgesOld2 = std::get<2>(q.front());
    auto openEdgesNew2 = std::get<3>(q.front());
    q.pop();

    // type 'rootNode' not possible, otherwise has to be catched

    if(currTree->child1 != nullptr) {

      if(currTree->child1->node1 == nullptr
         && currTree->child1->node2 == nullptr) {
        continue;
      }

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());

      childNode->freq
        = (currTree->child1->node1 == nullptr ? 0
                                              : currTree->child1->node1->freq)
          + (currTree->child1->node2 == nullptr
               ? 0
               : currTree->child1->node2->freq);
      childNode->type = currTree->child1->node1 == nullptr
                          ? currTree->child1->node2->type
                          : currTree->child1->node1->type;

      childNode->branchID = -1;

      if(alignmenttreeType == ttk::cta::lastMatchedValue) {

        childNode->scalarValue = currTree->child1->node2 == nullptr
                                   ? currTree->child1->node1->scalarValue
                                   : currTree->child1->node2->scalarValue;

      } else if(alignmenttreeType == ttk::cta::averageValues) {

        childNode->scalarValue = currTree->child1->node1 == nullptr
                                   ? currTree->child1->node2->scalarValue
                                 : currTree->child1->node2 == nullptr
                                   ? currTree->child1->node1->scalarValue
                                   : (currTree->child1->node1->scalarValue
                                        * currTree->child1->node1->freq
                                      + currTree->child1->node2->scalarValue)
                                       / childNode->freq;

      } else {

        if(currTree->child1->node1 == nullptr)
          childNode->scalarValue = currTree->child1->node2->scalarValue;
        else if(currTree->child1->node2 == nullptr)
          childNode->scalarValue = currTree->child1->node1->scalarValue;
        else {
          std::vector<float> values;
          for(auto p : currTree->child1->node1->nodeRefs) {
            values.push_back(
              contourtrees[p.first]->getGraph().first[p.second]->scalarValue);
          }
          values.push_back(res->node2->scalarValue);
          std::sort(values.begin(), values.end());
          float newMedian = values[values.size() / 2];
          childNode->scalarValue = newMedian;
        }
      }

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      if(currTree->child1->node1 != nullptr)
        childNode->nodeRefs.insert(childNode->nodeRefs.end(),
                                   currTree->child1->node1->nodeRefs.begin(),
                                   currTree->child1->node1->nodeRefs.end());
      if(currTree->child1->node2 != nullptr)
        childNode->nodeRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child1->node2->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());

      if(alignmenttreeType == ttk::cta::lastMatchedValue) {

        childEdge->area = currTree->child1->node2 == nullptr
                            ? currTree->child1->node1->area
                            : currTree->child1->node2->area;

        childEdge->scalardistance
          = currTree->child1->node2 == nullptr
              ? currTree->child1->node1->scalardistanceParent
              : currTree->child1->node2->scalardistanceParent;

      } else if(alignmenttreeType == ttk::cta::averageValues) {

        childEdge->area
          = currTree->child1->node1 == nullptr ? currTree->child1->node2->area
            : currTree->child1->node2 == nullptr
              ? currTree->child1->node1->area
              : (currTree->child1->node1->area * currTree->child1->node1->freq
                 + currTree->child1->node2->area)
                  / childNode->freq;

        childEdge->scalardistance
          = currTree->child1->node1 == nullptr
              ? currTree->child1->node2->scalardistanceParent
            : currTree->child1->node2 == nullptr
              ? currTree->child1->node1->scalardistanceParent
              : (currTree->child1->node1->scalardistanceParent
                   * currTree->child1->node1->freq
                 + currTree->child1->node2->scalardistanceParent)
                  / childNode->freq;

      } else {

        if(currTree->child1->node1 == nullptr) {
          childEdge->area = currTree->child1->node2->area;
          childEdge->scalardistance
            = currTree->child1->node2->scalardistanceParent;
        } else if(currTree->child1->node2 == nullptr) {
          childEdge->area = currTree->child1->node1->scalarValue;
          childEdge->scalardistance
            = currTree->child1->node1->scalardistanceParent;
        } else {
          std::vector<float> values;
          for(auto p : currTree->child1->node1->arcRefs) {
            values.push_back(
              contourtrees[p.first]->getGraph().second[p.second]->area);
          }
          values.push_back(currTree->child1->node2->area);
          std::sort(values.begin(), values.end());
          float newMedian = values[values.size() / 2];
          childEdge->area = newMedian;
          values.clear();
          for(auto p : currTree->child1->node1->arcRefs) {
            values.push_back(contourtrees[p.first]
                               ->getGraph()
                               .second[p.second]
                               ->scalardistance);
          }
          values.push_back(currTree->child1->node2->scalardistanceParent);
          std::sort(values.begin(), values.end());
          newMedian = values[values.size() / 2];
          childEdge->scalardistance = newMedian;
        }
      }

      childEdge->region = currTree->child1->node2 == nullptr
                            ? currTree->child1->node1->region
                            : currTree->child1->node2->region;

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      if(currTree->child1->node1 != nullptr) {
        for(const auto &e : openEdgesOld1) {
          e->arcRefs.insert(e->arcRefs.end(),
                            currTree->child1->node1->arcRefs.begin(),
                            currTree->child1->node1->arcRefs.end());
        }
        openEdgesOld1.clear();
        childEdge->arcRefs.insert(childEdge->arcRefs.end(),
                                  currTree->child1->node1->arcRefs.begin(),
                                  currTree->child1->node1->arcRefs.end());
      } else
        openEdgesOld1.push_back(childEdge);
      if(currTree->child1->node2 != nullptr) {
        for(const auto &e : openEdgesNew1) {
          e->arcRefs.push_back(
            std::make_pair((int)contourtrees.size() - 1,
                           currTree->child1->node2->arcRefs[0].second));
        }
        openEdgesNew1.clear();
        childEdge->arcRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child1->node2->arcRefs[0].second));
      } else
        openEdgesNew1.push_back(childEdge);

      childEdge->volume = childEdge->area * childEdge->scalardistance;

      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      q.push(std::make_tuple(
        currTree->child1, childNode, openEdgesOld1, openEdgesNew1));
    }

    if(currTree->child2 != nullptr) {

      if(currTree->child2->node1 == nullptr
         && currTree->child2->node2 == nullptr) {
        continue;
      }

      std::shared_ptr<ttk::cta::AlignmentNode> childNode(
        new ttk::cta::AlignmentNode());

      childNode->freq
        = (currTree->child2->node1 == nullptr ? 0
                                              : currTree->child2->node1->freq)
          + (currTree->child2->node2 == nullptr
               ? 0
               : currTree->child2->node2->freq);
      childNode->type = currTree->child2->node1 == nullptr
                          ? currTree->child2->node2->type
                          : currTree->child2->node1->type;

      childNode->branchID = -1;

      if(alignmenttreeType == ttk::cta::lastMatchedValue) {

        childNode->scalarValue = currTree->child2->node2 == nullptr
                                   ? currTree->child2->node1->scalarValue
                                   : currTree->child2->node2->scalarValue;

      } else if(alignmenttreeType == ttk::cta::averageValues) {

        childNode->scalarValue = currTree->child2->node1 == nullptr
                                   ? currTree->child2->node2->scalarValue
                                 : currTree->child2->node2 == nullptr
                                   ? currTree->child2->node1->scalarValue
                                   : (currTree->child2->node1->scalarValue
                                        * currTree->child2->node1->freq
                                      + currTree->child2->node2->scalarValue)
                                       / childNode->freq;

      } else {

        if(currTree->child2->node1 == nullptr)
          childNode->scalarValue = currTree->child2->node2->scalarValue;
        else if(currTree->child2->node2 == nullptr)
          childNode->scalarValue = currTree->child2->node1->scalarValue;
        else {
          std::vector<float> values;
          for(auto p : currTree->child2->node1->nodeRefs) {
            values.push_back(
              contourtrees[p.first]->getGraph().first[p.second]->scalarValue);
          }
          values.push_back(res->node2->scalarValue);
          std::sort(values.begin(), values.end());
          float newMedian = values[values.size() / 2];
          childNode->scalarValue = newMedian;
        }
      }

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      if(currTree->child2->node1 != nullptr)
        childNode->nodeRefs.insert(childNode->nodeRefs.end(),
                                   currTree->child2->node1->nodeRefs.begin(),
                                   currTree->child2->node1->nodeRefs.end());
      if(currTree->child2->node2 != nullptr)
        childNode->nodeRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child2->node2->nodeRefs[0].second));

      std::shared_ptr<ttk::cta::AlignmentEdge> childEdge(
        new ttk::cta::AlignmentEdge());

      if(alignmenttreeType == ttk::cta::lastMatchedValue) {

        childEdge->area = currTree->child2->node2 == nullptr
                            ? currTree->child2->node1->area
                            : currTree->child2->node2->area;

        childEdge->scalardistance
          = currTree->child2->node2 == nullptr
              ? currTree->child2->node1->scalardistanceParent
              : currTree->child2->node2->scalardistanceParent;

      } else if(alignmenttreeType == ttk::cta::averageValues) {

        childEdge->area
          = currTree->child2->node1 == nullptr ? currTree->child2->node2->area
            : currTree->child2->node2 == nullptr
              ? currTree->child2->node1->area
              : (currTree->child2->node1->area * currTree->child2->node1->freq
                 + currTree->child2->node2->area)
                  / childNode->freq;

        childEdge->scalardistance
          = currTree->child2->node1 == nullptr
              ? currTree->child2->node2->scalardistanceParent
            : currTree->child2->node2 == nullptr
              ? currTree->child2->node1->scalardistanceParent
              : (currTree->child2->node1->scalardistanceParent
                   * currTree->child2->node1->freq
                 + currTree->child2->node2->scalardistanceParent)
                  / childNode->freq;

      } else {

        if(currTree->child2->node1 == nullptr) {
          childEdge->area = currTree->child2->node2->area;
          childEdge->scalardistance
            = currTree->child2->node2->scalardistanceParent;
        } else if(currTree->child2->node2 == nullptr) {
          childEdge->area = currTree->child2->node1->scalarValue;
          childEdge->scalardistance
            = currTree->child2->node1->scalardistanceParent;
        } else {
          std::vector<float> values;
          for(auto p : currTree->child2->node1->arcRefs) {
            values.push_back(
              contourtrees[p.first]->getGraph().second[p.second]->area);
          }
          values.push_back(currTree->child2->node2->area);
          std::sort(values.begin(), values.end());
          float newMedian = values[values.size() / 2];
          childEdge->area = newMedian;
          values.clear();
          for(auto p : currTree->child2->node1->arcRefs) {
            values.push_back(contourtrees[p.first]
                               ->getGraph()
                               .second[p.second]
                               ->scalardistance);
          }
          values.push_back(currTree->child2->node2->scalardistanceParent);
          std::sort(values.begin(), values.end());
          newMedian = values[values.size() / 2];
          childEdge->scalardistance = newMedian;
        }
      }

      childEdge->region = currTree->child2->node2 == nullptr
                            ? currTree->child2->node1->region
                            : currTree->child2->node2->region;

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      if(currTree->child2->node1 != nullptr) {
        for(const auto &e : openEdgesOld2) {
          e->arcRefs.insert(e->arcRefs.end(),
                            currTree->child2->node1->arcRefs.begin(),
                            currTree->child2->node1->arcRefs.end());
        }
        openEdgesOld2.clear();
        childEdge->arcRefs.insert(childEdge->arcRefs.end(),
                                  currTree->child2->node1->arcRefs.begin(),
                                  currTree->child2->node1->arcRefs.end());
      } else
        openEdgesOld2.push_back(childEdge);
      if(currTree->child2->node2 != nullptr) {
        for(const auto &e : openEdgesNew1) {
          e->arcRefs.push_back(
            std::make_pair((int)contourtrees.size() - 1,
                           currTree->child2->node2->arcRefs[0].second));
        }
        openEdgesNew2.clear();
        childEdge->arcRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child2->node2->arcRefs[0].second));
      } else
        openEdgesNew2.push_back(childEdge);

      childEdge->volume = childEdge->area * childEdge->scalardistance;

      childEdge->node1 = currNode;
      childEdge->node2 = childNode;

      childNode->edgeList.push_back(childEdge);
      currNode->edgeList.push_back(childEdge);

      nodes.push_back(childNode);
      arcs.push_back(childEdge);

      // q.push(std::make_pair(currTree->child2,childNode));
      q.push(std::make_tuple(
        currTree->child2, childNode, openEdgesOld2, openEdgesNew2));
    }
  }
}

/*
=====================================================================================================================
 getters and setters
=====================================================================================================================
*/

int ttk::ContourTreeAlignment::getAlignmentRootIdx() {
  int idx = 0;
  for(const auto &node : nodes) {
    if(node == alignmentRoot)
      break;
    idx++;
  }
  return idx;
}

std::vector<std::shared_ptr<ContourTree>>
  ttk::ContourTreeAlignment::getContourTrees() {

  return contourtrees;
}

std::vector<std::pair<std::vector<std::shared_ptr<ttk::cta::CTNode>>,
                      std::vector<std::shared_ptr<ttk::cta::CTEdge>>>>
  ttk::ContourTreeAlignment::getGraphs() {

  std::vector<std::pair<std::vector<std::shared_ptr<ttk::cta::CTNode>>,
                        std::vector<std::shared_ptr<ttk::cta::CTEdge>>>>
    trees_simplified;

  for(const auto &ct : contourtrees) {
    trees_simplified.push_back(ct->getGraph());
  }

  return trees_simplified;
}

std::pair<std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>,
          std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>>
  ttk::ContourTreeAlignment::getAlignmentGraph() {

  return std::make_pair(nodes, arcs);
}

std::shared_ptr<ttk::cta::BinaryTree>
  ttk::ContourTreeAlignment::getAlignmentGraphRooted() {

  std::shared_ptr<ttk::cta::AlignmentNode> root = nodes[0];
  float maxScalar = FLT_MIN;
  for(const auto &node : nodes) {

    // if(node->type==maxNode && node->edgeList[0]->scalardistance > maxRange){
    if(node->scalarValue > maxScalar) {
      root = node;
      maxScalar = node->scalarValue;
    }
  }

  return rootAtNode(root);
}

/*
=====================================================================================================================
 aligning two trees
=====================================================================================================================
*/

std::pair<float, std::shared_ptr<ttk::cta::AlignmentTree>>
  ttk::ContourTreeAlignment::getAlignmentBinary(
    const std::shared_ptr<ttk::cta::BinaryTree> &t1,
    const std::shared_ptr<ttk::cta::BinaryTree> &t2) {

  // initialize memoization tables
  std::vector<std::vector<float>> memT(
    t1->size + 1, std::vector<float>(t2->size + 1, -1));
  std::vector<std::vector<float>> memF(
    t1->size + 1, std::vector<float>(t2->size + 1, -1));

  // compute table of distances
  float dist = alignTreeBinary(t1, t2, memT, memF);

  // backtrace through the table to get the alignment
  std::shared_ptr<ttk::cta::AlignmentTree> res
    = traceAlignmentTree(t1, t2, memT, memF);

  return std::make_pair(dist, res);
}

float ttk::ContourTreeAlignment::alignTreeBinary(
  const std::shared_ptr<ttk::cta::BinaryTree> &t1,
  const std::shared_ptr<ttk::cta::BinaryTree> &t2,
  std::vector<std::vector<float>> &memT,
  std::vector<std::vector<float>> &memF) {

  // base cases for matching to empty tree

  if(t1 == nullptr && t2 == nullptr) {
    if(memT[0][0] < 0) {
      memT[0][0] = 0;
    }
    return memT[0][0];
  }

  else if(t1 == nullptr) {
    if(memT[0][t2->id] < 0) {
      memT[0][t2->id]
        = editCost(nullptr, t2) + alignForestBinary(nullptr, t2, memT, memF);
    }
    return memT[0][t2->id];
  }

  else if(t2 == nullptr) {
    if(memT[t1->id][0] < 0) {
      memT[t1->id][0]
        = editCost(t1, nullptr) + alignForestBinary(t1, nullptr, memT, memF);
    }
    return memT[t1->id][0];
  }

  // find optimal possible matching in other cases

  else {
    if(memT[t1->id][t2->id] < 0) {

      // match t1 to t2 and then try to match their children
      memT[t1->id][t2->id]
        = editCost(t1, t2) + alignForestBinary(t1, t2, memT, memF);

      // match t1 to blank, one of its children to t2, the other to blank (try
      // both children)
      if(t1->size > 1)
        memT[t1->id][t2->id]
          = std::min(editCost(t1, nullptr)
                       + alignTreeBinary(t1->child2, nullptr, memT, memF)
                       + alignTreeBinary(t1->child1, t2, memT, memF),
                     memT[t1->id][t2->id]);
      if(t1->size > 1)
        memT[t1->id][t2->id]
          = std::min(editCost(t1, nullptr)
                       + alignTreeBinary(t1->child1, nullptr, memT, memF)
                       + alignTreeBinary(t1->child2, t2, memT, memF),
                     memT[t1->id][t2->id]);

      // match t2 to blank, one of its children to t1, the other to blank (try
      // both children)
      if(t2->size > 1)
        memT[t1->id][t2->id]
          = std::min(editCost(nullptr, t2)
                       + alignTreeBinary(nullptr, t2->child2, memT, memF)
                       + alignTreeBinary(t1, t2->child1, memT, memF),
                     memT[t1->id][t2->id]);
      if(t2->size > 1)
        memT[t1->id][t2->id]
          = std::min(editCost(nullptr, t2)
                       + alignTreeBinary(nullptr, t2->child1, memT, memF)
                       + alignTreeBinary(t1, t2->child2, memT, memF),
                     memT[t1->id][t2->id]);
    }
    return memT[t1->id][t2->id];
  }
}

float ttk::ContourTreeAlignment::alignForestBinary(
  const std::shared_ptr<ttk::cta::BinaryTree> &t1,
  const std::shared_ptr<ttk::cta::BinaryTree> &t2,
  std::vector<std::vector<float>> &memT,
  std::vector<std::vector<float>> &memF) {

  // base cases for matching to empty tree

  if(t1 == nullptr && t2 == nullptr) {
    if(memF[0][0] < 0) {
      memF[0][0] = 0;
    }
    return memF[0][0];
  }

  else if(t1 == nullptr) {
    if(memF[0][t2->id] < 0) {
      memF[0][t2->id] = 0;
      memF[0][t2->id] += alignTreeBinary(nullptr, t2->child1, memT, memF);
      memF[0][t2->id] += alignTreeBinary(nullptr, t2->child2, memT, memF);
    }
    return memF[0][t2->id];
  }

  else if(t2 == nullptr) {
    if(memF[t1->id][0] < 0) {
      memF[t1->id][0] = 0;
      memF[t1->id][0] += alignTreeBinary(t1->child1, nullptr, memT, memF);
      memF[t1->id][0] += alignTreeBinary(t1->child2, nullptr, memT, memF);
    }
    return memF[t1->id][0];
  }

  // find optimal possible matching in other cases

  else {
    if(memF[t1->id][t2->id] < 0) {

      memF[t1->id][t2->id] = FLT_MAX;

      if(t1->child2 != nullptr && t1->child2->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(t1->child2, nullptr)
                       + alignForestBinary(t1->child2, t2, memT, memF)
                       + alignTreeBinary(t1->child1, nullptr, memT, memF));
      if(t1->child1 != nullptr && t1->child1->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(t1->child1, nullptr)
                       + alignForestBinary(t1->child1, t2, memT, memF)
                       + alignTreeBinary(t1->child2, nullptr, memT, memF));

      if(t2->child2 != nullptr && t2->child2->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(nullptr, t2->child2)
                       + alignForestBinary(t1, t2->child2, memT, memF)
                       + alignTreeBinary(nullptr, t2->child1, memT, memF));
      if(t2->child1 != nullptr && t2->child1->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(nullptr, t2->child1)
                       + alignForestBinary(t1, t2->child1, memT, memF)
                       + alignTreeBinary(nullptr, t2->child2, memT, memF));

      memF[t1->id][t2->id]
        = std::min(memF[t1->id][t2->id],
                   alignTreeBinary(t1->child1, t2->child1, memT, memF)
                     + alignTreeBinary(t1->child2, t2->child2, memT, memF));
      memF[t1->id][t2->id]
        = std::min(memF[t1->id][t2->id],
                   alignTreeBinary(t1->child1, t2->child2, memT, memF)
                     + alignTreeBinary(t1->child2, t2->child1, memT, memF));
    }
    return memF[t1->id][t2->id];
  }
}

float ttk::ContourTreeAlignment::editCost(
  const std::shared_ptr<ttk::cta::BinaryTree> &t1,
  const std::shared_ptr<ttk::cta::BinaryTree> &t2) {

  float v1 = 0, v2 = 0;
  if(t1.get() != nullptr)
    v1 = arcMatchMode == ttk::cta::persistence ? t1->scalardistanceParent
         : arcMatchMode == ttk::cta::area      ? t1->area
                                               : t1->volume;
  if(t2.get() != nullptr)
    v2 = arcMatchMode == ttk::cta::persistence ? t2->scalardistanceParent
         : arcMatchMode == ttk::cta::area      ? t2->area
                                               : t2->volume;

  if(arcMatchMode == ttk::cta::overlap) {
    // this->printMsg("overlap");
    int unionsize = 0;
    int intersectionsize = 0;
    size_t i = 0;
    size_t j = 0;
    if(t1 && t2) {
      while(i < t1->region.size() && j < t2->region.size()) {
        int vi = t1->region[i];
        int vj = t2->region[j];
        if(vi == vj) {
          intersectionsize++;
          i++;
          j++;
        } else if(vi < vj)
          i++;
        else
          j++;
        unionsize++;
      }
      unionsize += i == (t1->region.size()) ? t2->region.size() - j
                                            : t1->region.size() - i;
      if(false) {
        std::cout << "========================\nRegion 1: ";
        for(size_t k = 0; k < t1->region.size(); k++)
          std::cout << t1->region[k] << " ";
        std::cout << std::endl;
        std::cout << "Region 2: ";
        for(size_t k = 0; k < t2->region.size(); k++)
          std::cout << t2->region[k] << " ";
        std::cout << std::endl;
        std::cout << "Intersection size: " << intersectionsize
                  << "\nUnion size: " << unionsize << std::endl;
      }
    }

    if(t1 == nullptr && t2 == nullptr)
      return 0;

    else if(t1 == nullptr || t2 == nullptr)
      return weightCombinatorialMatch + weightArcMatch * 1
             + weightScalarValueMatch;

    else if(t1->type == t2->type)
      return weightArcMatch * (1 - ((float)intersectionsize / (float)unionsize))
             + weightScalarValueMatch
                 * std::abs(t1->scalarValue - t2->scalarValue);

    else
      return FLT_MAX;
  }

  if(t1 == nullptr && t2 == nullptr)
    return 0;

  else if(t1 == nullptr)
    return weightCombinatorialMatch + weightArcMatch * v2
           + weightScalarValueMatch;

  else if(t2 == nullptr)
    return weightCombinatorialMatch + weightArcMatch * v1
           + weightScalarValueMatch;

  else if(t1->type == t2->type)
    return weightArcMatch * std::abs(v1 - v2)
           + weightScalarValueMatch
               * std::abs(t1->scalarValue - t2->scalarValue);

  else
    return FLT_MAX;
}

std::shared_ptr<ttk::cta::AlignmentTree>
  ttk::ContourTreeAlignment::traceAlignmentTree(
    const std::shared_ptr<ttk::cta::BinaryTree> &t1,
    const std::shared_ptr<ttk::cta::BinaryTree> &t2,
    std::vector<std::vector<float>> &memT,
    std::vector<std::vector<float>> &memF) {

  if(t1 == nullptr)
    return traceNullAlignment(t2, false);
  if(t2 == nullptr)
    return traceNullAlignment(t1, true);

  auto id = [](std::shared_ptr<ttk::cta::BinaryTree> &t) {
    return t == nullptr ? 0 : t->id;
  };

  if(memT[t1->id][t2->id] == editCost(t1, t2) + memF[t1->id][t2->id]) {

    std::shared_ptr<ttk::cta::AlignmentTree> resNode(
      new ttk::cta::AlignmentTree);

    resNode->node1 = t1;
    resNode->node2 = t2;
    resNode->child1 = nullptr;
    resNode->child2 = nullptr;
    resNode->height = 0;
    resNode->size = 1;

    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> resChildren
      = traceAlignmentForest(t1, t2, memT, memF);

    if(resChildren.size() > 0)
      resNode->child1 = resChildren[0];
    if(resChildren.size() > 1)
      resNode->child2 = resChildren[1];

    for(const auto &child : resChildren) {

      if(child != nullptr) {
        resNode->height = std::max(resNode->height, child->height + 1);
        resNode->size += child->size;
      }
    }

    return resNode;
  }

  if(memT[t1->id][t2->id]
     == editCost(t1, nullptr) + memT[id(t1->child2)][0]
          + memT[id(t1->child1)]
                [t2->id] /* && t1->type != maxNode && t1->type != minNode */) {

    std::shared_ptr<ttk::cta::AlignmentTree> resChild1
      = traceAlignmentTree(t1->child1, t2, memT, memF);
    std::shared_ptr<ttk::cta::AlignmentTree> resChild2
      = traceNullAlignment(t1->child2, true);
    std::shared_ptr<ttk::cta::AlignmentTree> res(new ttk::cta::AlignmentTree());
    res->node1 = t1;
    res->node2 = nullptr;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == nullptr ? 0 : resChild1->size;
    res->size += resChild2 == nullptr ? 1 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == nullptr ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == nullptr ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(t1, nullptr) + memT[id(t1->child1)][0]
          + memT[id(t1->child2)]
                [t2->id] /* && t1->type != maxNode && t1->type != minNode */) {

    std::shared_ptr<ttk::cta::AlignmentTree> resChild1
      = traceAlignmentTree(t1->child2, t2, memT, memF);
    std::shared_ptr<ttk::cta::AlignmentTree> resChild2
      = traceNullAlignment(t1->child1, true);
    std::shared_ptr<ttk::cta::AlignmentTree> res(new ttk::cta::AlignmentTree());
    res->node1 = t1;
    res->node2 = nullptr;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == nullptr ? 0 : resChild1->size;
    res->size += resChild2 == nullptr ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == nullptr ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == nullptr ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(nullptr, t2) + memT[0][id(t2->child2)]
          + memT[t1->id][id(
            t2->child1)] /* && t2->type != maxNode && t2->type != minNode */) {

    std::shared_ptr<ttk::cta::AlignmentTree> resChild1
      = traceAlignmentTree(t1, t2->child1, memT, memF);
    std::shared_ptr<ttk::cta::AlignmentTree> resChild2
      = traceNullAlignment(t2->child2, false);
    std::shared_ptr<ttk::cta::AlignmentTree> res(new ttk::cta::AlignmentTree());
    res->node1 = nullptr;
    res->node2 = t2;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == nullptr ? 0 : resChild1->size;
    res->size += resChild2 == nullptr ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == nullptr ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == nullptr ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(nullptr, t2) + memT[0][id(t2->child1)]
          + memT[t1->id][id(
            t2->child2)] /* && t2->type != maxNode && t2->type != minNode */) {

    std::shared_ptr<ttk::cta::AlignmentTree> resChild1
      = traceAlignmentTree(t1, t2->child2, memT, memF);
    std::shared_ptr<ttk::cta::AlignmentTree> resChild2
      = traceNullAlignment(t2->child1, false);
    std::shared_ptr<ttk::cta::AlignmentTree> res(new ttk::cta::AlignmentTree());
    res->node1 = nullptr;
    res->node2 = t2;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == nullptr ? 0 : resChild1->size;
    res->size += resChild2 == nullptr ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == nullptr ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == nullptr ? 0 : resChild2->height + 1);
    return res;
  }

  printErr("Alignment computation failed. Traceback of memoization table not "
           "possible.");
  return std::shared_ptr<ttk::cta::AlignmentTree>(
    new ttk::cta::AlignmentTree());
}

std::vector<std::shared_ptr<ttk::cta::AlignmentTree>>
  ttk::ContourTreeAlignment::traceAlignmentForest(
    const std::shared_ptr<ttk::cta::BinaryTree> &t1,
    const std::shared_ptr<ttk::cta::BinaryTree> &t2,
    std::vector<std::vector<float>> &memT,
    std::vector<std::vector<float>> &memF) {

  if(t1 == nullptr && t2 == nullptr)
    return std::vector<std::shared_ptr<ttk::cta::AlignmentTree>>();
  if(t1 == nullptr) {
    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;
    if(t2->child1 != nullptr)
      res.push_back(traceNullAlignment(t2->child1, false));
    if(t2->child2 != nullptr)
      res.push_back(traceNullAlignment(t2->child2, false));
  }
  if(t2 == nullptr) {
    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;
    if(t1->child1 != nullptr)
      res.push_back(traceNullAlignment(t1->child1, true));
    if(t1->child2 != nullptr)
      res.push_back(traceNullAlignment(t1->child2, true));
  }

  auto id = [](const std::shared_ptr<ttk::cta::BinaryTree> &t) {
    return t == nullptr ? 0 : t->id;
  };

  if(memF[t1->id][t2->id]
     == memT[id(t1->child1)][id(t2->child1)]
          + memT[id(t1->child2)][id(t2->child2)]) {

    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;
    std::shared_ptr<ttk::cta::AlignmentTree> res1
      = traceAlignmentTree(t1->child1, t2->child1, memT, memF);
    if(res1 != nullptr)
      res.push_back(res1);
    std::shared_ptr<ttk::cta::AlignmentTree> res2
      = traceAlignmentTree(t1->child2, t2->child2, memT, memF);
    if(res2 != nullptr)
      res.push_back(res2);

    return res;
  }

  if(memF[t1->id][t2->id]
     == memT[id(t1->child1)][id(t2->child2)]
          + memT[id(t1->child2)][id(t2->child1)]) {

    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;
    std::shared_ptr<ttk::cta::AlignmentTree> res1
      = traceAlignmentTree(t1->child1, t2->child2, memT, memF);
    if(res1 != nullptr)
      res.push_back(res1);
    std::shared_ptr<ttk::cta::AlignmentTree> res2
      = traceAlignmentTree(t1->child2, t2->child1, memT, memF);
    if(res2 != nullptr)
      res.push_back(res2);

    return res;
  }

  if(memF[t1->id][t2->id]
     == editCost(t1->child1, nullptr) + memF[id(t1->child1)][t2->id]
          + memT[id(t1->child2)][0]) {

    if(t1->child1 != nullptr) {

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;

      std::shared_ptr<ttk::cta::AlignmentTree> t(new ttk::cta::AlignmentTree);
      t->node1 = t1->child1;
      t->node2 = nullptr;

      t->child1 = nullptr;
      t->child2 = nullptr;

      t->size = 1;
      t->height = 0;

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> resChildren
        = traceAlignmentForest(t1->child1, t2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(const auto &c : resChildren) {
        if(c != nullptr) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t1->child2 != nullptr)
        res.push_back(traceNullAlignment(t1->child2, true));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(t1->child2, nullptr) + memF[id(t1->child2)][t2->id]
          + memT[id(t1->child1)][0]) {

    if(t1->child2 != nullptr) {

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;

      std::shared_ptr<ttk::cta::AlignmentTree> t(new ttk::cta::AlignmentTree);
      t->node1 = t1->child2;
      t->node2 = nullptr;

      t->child1 = nullptr;
      t->child2 = nullptr;

      t->size = 1;
      t->height = 0;

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> resChildren
        = traceAlignmentForest(t1->child2, t2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(const auto &c : resChildren) {
        if(c != nullptr) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t1->child1 != nullptr)
        res.push_back(traceNullAlignment(t1->child1, true));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(nullptr, t2->child1) + memF[t1->id][id(t2->child1)]
          + memT[0][id(t2->child2)]) {

    if(t2->child1 != nullptr) {

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;

      std::shared_ptr<ttk::cta::AlignmentTree> t(new ttk::cta::AlignmentTree);
      t->node1 = nullptr;
      t->node2 = t2->child1;

      t->child1 = nullptr;
      t->child2 = nullptr;

      t->size = 1;
      t->height = 0;

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> resChildren
        = traceAlignmentForest(t1, t2->child1, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(const auto &c : resChildren) {
        if(c != nullptr) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t2->child2 != nullptr)
        res.push_back(traceNullAlignment(t2->child2, false));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(nullptr, t2->child2) + memF[t1->id][id(t2->child2)]
          + memT[0][id(t2->child1)]) {

    if(t2->child2 != nullptr) {

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> res;

      std::shared_ptr<ttk::cta::AlignmentTree> t(new ttk::cta::AlignmentTree);
      t->node1 = nullptr;
      t->node2 = t2->child2;

      t->child1 = nullptr;
      t->child2 = nullptr;

      t->size = 1;
      t->height = 0;

      std::vector<std::shared_ptr<ttk::cta::AlignmentTree>> resChildren
        = traceAlignmentForest(t1, t2->child2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(const auto &c : resChildren) {
        if(c != nullptr) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t2->child1 != nullptr)
        res.push_back(traceNullAlignment(t2->child1, false));

      return res;
    }
  }

  printErr("Alignment computation failed. Traceback of memoization table not "
           "possible.");
  return std::vector<std::shared_ptr<ttk::cta::AlignmentTree>>();
}

std::shared_ptr<ttk::cta::AlignmentTree>
  ttk::ContourTreeAlignment::traceNullAlignment(
    const std::shared_ptr<ttk::cta::BinaryTree> &t, bool first) {

  if(t == nullptr)
    return nullptr;
  std::shared_ptr<ttk::cta::AlignmentTree> at(new ttk::cta::AlignmentTree());
  at->node1 = first ? t : nullptr;
  at->node2 = first ? nullptr : t;
  at->height = t->height;
  at->size = t->size;
  at->child1 = traceNullAlignment(t->child1, first);
  at->child2 = traceNullAlignment(t->child2, first);
  return at;
}

/*
=====================================================================================================================
 helper functions
=====================================================================================================================
*/

bool ttk::ContourTreeAlignment::isBinary(
  const std::shared_ptr<ttk::cta::Tree> &t) {

  if(t->children.size() > 2)
    return false;
  else {
    bool res = true;
    for(const auto &c : t->children) {
      res = res & isBinary(c);
    }
    return res;
  }
}

std::shared_ptr<ttk::cta::BinaryTree> ttk::ContourTreeAlignment::rootAtNode(
  const std::shared_ptr<ttk::cta::AlignmentNode> &root) {

  int id = 1;
  std::shared_ptr<ttk::cta::BinaryTree> t
    = computeRootedTree(root, nullptr, id);
  return t;
}

std::shared_ptr<ttk::cta::BinaryTree>
  ttk::ContourTreeAlignment::computeRootedTree(
    const std::shared_ptr<ttk::cta::AlignmentNode> &node,
    const std::shared_ptr<ttk::cta::AlignmentEdge> &parent,
    int &id) {

  std::shared_ptr<ttk::cta::BinaryTree> t(new ttk::cta::BinaryTree);

  if(parent == nullptr) {
    t->scalardistanceParent = 10000;
    t->area = 10000;
    t->volume = 10000;
    t->region = std::vector<int>(1, -1);
  } else {
    t->scalardistanceParent = parent->scalardistance;
    t->area = parent->area;
    t->volume = t->area * t->scalardistanceParent;
    t->region = parent->region;
  }
  t->freq = node->freq;
  t->type = node->type;
  t->child1 = nullptr;
  t->child2 = nullptr;
  t->id = id;
  t->scalarValue = node->scalarValue;
  t->nodeRefs = node->nodeRefs;
  if(parent != nullptr)
    t->arcRefs = parent->arcRefs;
  else
    t->arcRefs = std::vector<std::pair<int, int>>();
  id++;

  t->size = 1;
  t->height = 0;

  std::vector<std::shared_ptr<ttk::cta::BinaryTree>> children;

  for(const auto &edge : node->edgeList) {

    if(edge != parent) {

      std::shared_ptr<ttk::cta::BinaryTree> child = computeRootedTree(
        edge->node1.lock() == node ? edge->node2.lock() : edge->node1.lock(),
        edge, id);
      children.push_back(child);
      t->size += child->size;
      if(t->height < child->height + 1)
        t->height = child->height + 1;
    }
  }

  if(children.size() > 0)
    t->child1 = children[0];
  if(children.size() > 1)
    t->child2 = children[1];

  return t;
}

std::shared_ptr<ttk::cta::BinaryTree>
  ttk::ContourTreeAlignment::computeRootedDualTree(
    const std::shared_ptr<ttk::cta::AlignmentEdge> &arc,
    bool parent1,
    int &id) {

  std::shared_ptr<ttk::cta::BinaryTree> t(new ttk::cta::BinaryTree);

  t->scalardistanceParent = arc->scalardistance;
  t->area = arc->area;
  t->volume = t->area * t->scalardistanceParent;
  t->freq = arc->freq;
  t->type = arc->node1.lock()->type == ttk::cta::maxNode
                || arc->node2.lock()->type == ttk::cta::maxNode
              ? ttk::cta::maxNode
            : arc->node1.lock()->type == ttk::cta::minNode
                || arc->node2.lock()->type == ttk::cta::minNode
              ? ttk::cta::minNode
              : ttk::cta::saddleNode;
  t->child1 = nullptr;
  t->child2 = nullptr;
  t->id = id;
  t->scalarValue = -1;
  t->arcRefs = arc->arcRefs;
  id++;

  t->size = 1;
  t->height = 0;

  std::vector<std::shared_ptr<ttk::cta::BinaryTree>> children;

  std::shared_ptr<ttk::cta::AlignmentNode> node
    = parent1 ? arc->node2.lock() : arc->node1.lock();

  for(const auto &edge : node->edgeList) {

    if(edge != arc) {

      std::shared_ptr<ttk::cta::BinaryTree> child = computeRootedDualTree(
        edge, edge->node1.lock() == node ? true : false, id);
      children.push_back(child);
      t->size += child->size;
      if(t->height < child->height + 1)
        t->height = child->height + 1;
    }
  }

  if(children.size() > 0)
    t->child1 = children[0];
  if(children.size() > 1)
    t->child2 = children[1];

  return t;
}
