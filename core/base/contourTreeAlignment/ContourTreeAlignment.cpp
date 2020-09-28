#include "ContourTreeAlignment.h"

#include "float.h"
#include <algorithm>
#include <limits.h>
#include <list>

///=====================================================================================================================
/// branch decomposition
///=====================================================================================================================

void ttk::ContourTreeAlignment::computeBranches() {

  // find global minimum
  AlignmentNode *minNode = nodes[0];
  for(size_t i = 1; i < nodes.size(); i++) {
    if(minNode->scalarValue > nodes[i]->scalarValue)
      minNode = nodes[i];
  }

  // find path to global max
  AlignmentNode *nextNode = minNode->edgeList[0]->node1 == minNode
                              ? minNode->edgeList[0]->node2
                              : minNode->edgeList[0]->node1;
  std::vector<AlignmentNode *> maxPath_ = pathToMax(nextNode, minNode).second;
  std::vector<AlignmentNode *> maxPath;
  maxPath.push_back(minNode);
  maxPath.insert(maxPath.end(), maxPath_.begin(), maxPath_.end());

  int currID = 0;
  minNode->branchID = 0;
  std::stack<std::vector<AlignmentNode *>> q;
  q.push(maxPath);

  while(!q.empty()) {

    std::vector<AlignmentNode *> path = q.top();
    q.pop();

    for(size_t i = 1; i < path.size() - 1; i++) {

      AlignmentNode *currNode = path[i];

      for(AlignmentEdge *cE : currNode->edgeList) {

        AlignmentNode *cN = cE->node1 == currNode ? cE->node2 : cE->node1;
        if(cN == path[i - 1])
          continue;
        if(cN == path[i + 1])
          continue;

        if(cN->scalarValue > currNode->scalarValue) {
          std::vector<AlignmentNode *> newPath_
            = pathToMax(cN, currNode).second;
          std::vector<AlignmentNode *> newPath;
          newPath.push_back(currNode);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        } else {
          std::vector<AlignmentNode *> newPath_
            = pathToMin(cN, currNode).second;
          std::vector<AlignmentNode *> newPath;
          newPath.push_back(currNode);
          newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
          q.push(newPath);
        }
      }

      currNode->branchID = currID;
    }

    for(AlignmentEdge *cE : path.back()->edgeList) {

      AlignmentNode *cN = cE->node1 == path.back() ? cE->node2 : cE->node1;
      if(cN == path[path.size() - 2])
        continue;

      if(cN->scalarValue > path.back()->scalarValue) {
        std::vector<AlignmentNode *> newPath_
          = pathToMax(cN, path.back()).second;
        std::vector<AlignmentNode *> newPath;
        newPath.push_back(path.back());
        newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
        q.push(newPath);
      } else {
        std::vector<AlignmentNode *> newPath_
          = pathToMin(cN, path.back()).second;
        std::vector<AlignmentNode *> newPath;
        newPath.push_back(path.back());
        newPath.insert(newPath.end(), newPath_.begin(), newPath_.end());
        q.push(newPath);
      }
    }

    path.back()->branchID = currID;
    currID++;
  }
}

std::pair<float, std::vector<AlignmentNode *>>
  ttk::ContourTreeAlignment::pathToMax(AlignmentNode *root,
                                       AlignmentNode *parent) {

  std::vector<AlignmentNode *> path;
  path.push_back(root);

  if(root->edgeList.size() == 1) {
    return std::make_pair(root->scalarValue, path);
  }

  std::vector<AlignmentNode *> bestPath;
  float bestVal = -FLT_MAX;

  for(AlignmentEdge *cE : root->edgeList) {

    AlignmentNode *nextNode = cE->node1 == root ? cE->node2 : cE->node1;
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

std::pair<float, std::vector<AlignmentNode *>>
  ttk::ContourTreeAlignment::pathToMin(AlignmentNode *root,
                                       AlignmentNode *parent) {

  std::vector<AlignmentNode *> path;
  path.push_back(root);

  if(root->edgeList.size() == 1) {
    return std::make_pair(root->scalarValue, path);
  }

  std::vector<AlignmentNode *> bestPath;
  float bestVal = FLT_MAX;

  for(AlignmentEdge *cE : root->edgeList) {

    AlignmentNode *nextNode = cE->node1 == root ? cE->node2 : cE->node1;
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

///=====================================================================================================================
/// iterated aligning
///=====================================================================================================================

bool ttk::ContourTreeAlignment::initialize(ContourTree *ct) {

  contourtrees.push_back(ct);

  // cancel computation if tree is not binary
  if(!ct->isBinary())
    return false;

  // compute initial Alignment from initial contourtree

  BinaryTree *t = ct->rootAtMax();

  std::queue<std::pair<BinaryTree *, AlignmentNode *>> q;

  AlignmentNode *currNode = new AlignmentNode();
  BinaryTree *currTree;

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

    if(currTree->child1 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();
      childNode->freq = 1;
      childNode->type = currTree->child1->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child1->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child1->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();
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

    if(currTree->child2 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();
      childNode->freq = 1;
      childNode->type = currTree->child2->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child2->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child2->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();
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

bool ttk::ContourTreeAlignment::initialize_consistentRoot(ContourTree *ct,
                                                          int rootIdx) {

  // cancel computation if tree is not binary
  if(!ct->isBinary())
    return false;

  contourtrees.push_back(ct);

  // compute initial Alignment from initial contourtree

  BinaryTree *t = ct->rootAtMax();

  std::queue<std::pair<BinaryTree *, AlignmentNode *>> q;

  AlignmentNode *currNode = new AlignmentNode();
  BinaryTree *currTree;

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

    if(currTree->child1 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();
      childNode->freq = 1;
      childNode->type = currTree->child1->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child1->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child1->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();
      childEdge->area = currTree->child1->area;
      childEdge->scalardistance = currTree->child1->scalardistanceParent;
      childEdge->volume = currTree->child1->volume;
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

    if(currTree->child2 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();
      childNode->freq = 1;
      childNode->type = currTree->child2->type;
      childNode->branchID = -1;
      childNode->scalarValue = currTree->child2->scalarValue;

      childNode->nodeRefs = std::vector<std::pair<int, int>>();
      childNode->nodeRefs.push_back(
        std::make_pair(0, currTree->child2->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();
      childEdge->area = currTree->child2->area;
      childEdge->scalardistance = currTree->child2->scalardistanceParent;
      childEdge->volume = currTree->child2->volume;
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

  ContourTree::deleteBinaryTree(t);

  alignmentRoot = nodes[rootIdx];
  alignmentRootIdx = rootIdx;

  return true;
}

bool ttk::ContourTreeAlignment::alignTree(ContourTree *ct) {

  contourtrees.push_back(ct);

  // cancel computation if tree not binary
  if(!ct->isBinary())
    return false;

  // compute optimal alignment between current alignment and new tree

  BinaryTree *t1, *t2;

  AlignmentTree *res = nullptr;
  float resVal = FLT_MAX;

  std::vector<AlignmentNode *> nodes1 = nodes;
  std::vector<CTNode *> nodes2 = ct->getGraph().first;

  int i = 0;
  int j = 0;

  for(AlignmentNode *node1 : nodes1) {

    t1 = this->rootAtNode(node1);

    j = 0;

    for(CTNode *node2 : nodes2) {

      t2 = ct->rootAtNode(node2);

      // compute matching

      if((node1->type == maxNode && node2->type == maxNode)
         || (node1->type == minNode && node2->type == minNode)) {

        std::pair<float, AlignmentTree *> match = getAlignmentBinary(t1, t2);

        if(match.first < resVal) {
          resVal = match.first;
          res = match.second;
        }
      }

      j++;
    }

    i++;
  }

  computeNewAlignmenttree(res);

  return true;
}

bool ttk::ContourTreeAlignment::alignTree_consistentRoot(ContourTree *ct) {

  // cancel computation if tree not binary
  if(!ct->isBinary())
    return false;

  contourtrees.push_back(ct);

  // compute optimal alignment between current alignment and new tree

  BinaryTree *t1, *t2;

  AlignmentTree *res = NULL;
  float resVal = FLT_MAX;

  std::vector<CTNode *> nodes2 = ct->getGraph().first;

  int i = 0;

  t1 = this->rootAtNode(alignmentRoot);

  for(CTNode *node2 : nodes2) {

    // compute matching

    if((alignmentRoot->type == maxNode && node2->type == maxNode)
       || (alignmentRoot->type == minNode && node2->type == minNode)) {

      t2 = ct->rootAtNode(node2);

      std::pair<float, AlignmentTree *> match = getAlignmentBinary(t1, t2);

      // if(match.second->node1==NULL || match.second->node2==NULL)
      // printWrn("Root not matched...");

      if(match.first < resVal) {
        if(res && res->node2)
          ContourTree::deleteBinaryTree(res->node2);
        if(res)
          deleteAlignmentTree(res);
        resVal = match.first;
        res = match.second;
      } else {
        ContourTree::deleteBinaryTree(t2);
        deleteAlignmentTree(match.second);
      }
    }

    i++;
  }

  computeNewAlignmenttree(res);

  ContourTree::deleteBinaryTree(res->node1);
  ContourTree::deleteBinaryTree(res->node2);
  deleteAlignmentTree(res);

  alignmentVal += resVal;

  alignmentRoot = nodes[0];

  return true;
}

void ttk::ContourTreeAlignment::computeNewAlignmenttree(AlignmentTree *res) {

  for(AlignmentNode *node : nodes) {
    delete node;
  }
  for(AlignmentEdge *edge : arcs) {
    delete edge;
  }

  nodes.clear();
  arcs.clear();

  // std::queue<std::pair<AlignmentTree*,AlignmentNode*>> q;
  std::queue<
    std::tuple<AlignmentTree *, AlignmentNode *, std::vector<AlignmentEdge *>,
               std::vector<AlignmentEdge *>>>
    q;

  AlignmentNode *currNode;
  AlignmentTree *currTree;

  currNode = new AlignmentNode();

  currNode->freq = (res->node1 == NULL ? 0 : res->node1->freq)
                   + (res->node2 == NULL ? 0 : res->node2->freq);
  currNode->type = res->node1 == NULL ? res->node2->type : res->node1->type;
  currNode->branchID = -1;

  if(alignmenttreeType == lastMatchedValue) {

    currNode->scalarValue
      = res->node2 == NULL ? res->node1->scalarValue : res->node2->scalarValue;

    currNode->nodeRefs = std::vector<std::pair<int, int>>();
    if(res->node1 != NULL)
      currNode->nodeRefs.insert(currNode->nodeRefs.end(),
                                res->node1->nodeRefs.begin(),
                                res->node1->nodeRefs.end());
    if(res->node2 != NULL)
      currNode->nodeRefs.push_back(std::make_pair(
        (int)contourtrees.size() - 1, res->node2->nodeRefs[0].second));

  } else if(alignmenttreeType == averageValues) {

    currNode->scalarValue = res->node1 == NULL
                              ? res->node2->scalarValue
                              : res->node2 == NULL
                                  ? res->node1->scalarValue
                                  : (res->node1->scalarValue * res->node1->freq
                                     + res->node2->scalarValue)
                                      / currNode->freq;

  } else {

    if(res->node1 == NULL)
      currNode->scalarValue = res->node2->scalarValue;
    else if(res->node2 == NULL)
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
  if(res->node1 != NULL)
    currNode->nodeRefs.insert(currNode->nodeRefs.end(),
                              res->node1->nodeRefs.begin(),
                              res->node1->nodeRefs.end());
  if(res->node2 != NULL)
    currNode->nodeRefs.push_back(std::make_pair(
      (int)contourtrees.size() - 1, res->node2->nodeRefs[0].second));

  nodes.push_back(currNode);

  // q.push(std::make_pair(res,currNode));
  q.push(std::make_tuple(res, currNode, std::vector<AlignmentEdge *>(),
                         std::vector<AlignmentEdge *>()));

  while(!q.empty()) {

    // currTree = q.front().first;
    // currNode = q.front().second;
    currTree = std::get<0>(q.front());
    currNode = std::get<1>(q.front());
    auto openEdgesOld1 = std::get<2>(q.front());
    auto openEdgesNew1 = std::get<3>(q.front());
    auto openEdgesOld2 = std::get<2>(q.front());
    auto openEdgesNew2 = std::get<3>(q.front());
    q.pop();

    // type 'rootNode' not possible, otherwise has to be catched

    if(currTree->child1 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();

      childNode->freq
        = (currTree->child1->node1 == NULL ? 0 : currTree->child1->node1->freq)
          + (currTree->child1->node2 == NULL ? 0
                                             : currTree->child1->node2->freq);
      childNode->type = currTree->child1->node1 == NULL
                          ? currTree->child1->node2->type
                          : currTree->child1->node1->type;

      childNode->branchID = -1;

      if(alignmenttreeType == lastMatchedValue) {

        childNode->scalarValue = currTree->child1->node2 == NULL
                                   ? currTree->child1->node1->scalarValue
                                   : currTree->child1->node2->scalarValue;

      } else if(alignmenttreeType == averageValues) {

        childNode->scalarValue
          = currTree->child1->node1 == NULL
              ? currTree->child1->node2->scalarValue
              : currTree->child1->node2 == NULL
                  ? currTree->child1->node1->scalarValue
                  : (currTree->child1->node1->scalarValue
                       * currTree->child1->node1->freq
                     + currTree->child1->node2->scalarValue)
                      / childNode->freq;

      } else {

        if(currTree->child1->node1 == NULL)
          childNode->scalarValue = currTree->child1->node2->scalarValue;
        else if(currTree->child1->node2 == NULL)
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
      if(currTree->child1->node1 != NULL)
        childNode->nodeRefs.insert(childNode->nodeRefs.end(),
                                   currTree->child1->node1->nodeRefs.begin(),
                                   currTree->child1->node1->nodeRefs.end());
      if(currTree->child1->node2 != NULL)
        childNode->nodeRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child1->node2->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();

      if(alignmenttreeType == lastMatchedValue) {

        childEdge->area = currTree->child1->node2 == NULL
                            ? currTree->child1->node1->area
                            : currTree->child1->node2->area;

        childEdge->scalardistance
          = currTree->child1->node2 == NULL
              ? currTree->child1->node1->scalardistanceParent
              : currTree->child1->node2->scalardistanceParent;

      } else if(alignmenttreeType == averageValues) {

        childEdge->area = currTree->child1->node1 == NULL
                            ? currTree->child1->node2->area
                            : currTree->child1->node2 == NULL
                                ? currTree->child1->node1->area
                                : (currTree->child1->node1->area
                                     * currTree->child1->node1->freq
                                   + currTree->child1->node2->area)
                                    / childNode->freq;

        childEdge->scalardistance
          = currTree->child1->node1 == NULL
              ? currTree->child1->node2->scalardistanceParent
              : currTree->child1->node2 == NULL
                  ? currTree->child1->node1->scalardistanceParent
                  : (currTree->child1->node1->scalardistanceParent
                       * currTree->child1->node1->freq
                     + currTree->child1->node2->scalardistanceParent)
                      / childNode->freq;

      } else {

        if(currTree->child1->node1 == NULL) {
          childEdge->area = currTree->child1->node2->area;
          childEdge->scalardistance
            = currTree->child1->node2->scalardistanceParent;
        } else if(currTree->child1->node2 == NULL) {
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

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      if(currTree->child1->node1 != NULL) {
        for(AlignmentEdge *e : openEdgesOld1) {
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
      if(currTree->child1->node2 != NULL) {
        for(AlignmentEdge *e : openEdgesNew1) {
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

      // q.push(std::make_pair(currTree->child1,childNode));
      q.push(std::make_tuple(
        currTree->child1, childNode, openEdgesOld1, openEdgesNew1));
    }

    if(currTree->child2 != NULL) {

      AlignmentNode *childNode = new AlignmentNode();

      childNode->freq
        = (currTree->child2->node1 == NULL ? 0 : currTree->child2->node1->freq)
          + (currTree->child2->node2 == NULL ? 0
                                             : currTree->child2->node2->freq);
      childNode->type = currTree->child2->node1 == NULL
                          ? currTree->child2->node2->type
                          : currTree->child2->node1->type;

      childNode->branchID = -1;

      if(alignmenttreeType == lastMatchedValue) {

        childNode->scalarValue = currTree->child2->node2 == NULL
                                   ? currTree->child2->node1->scalarValue
                                   : currTree->child2->node2->scalarValue;

      } else if(alignmenttreeType == averageValues) {

        childNode->scalarValue
          = currTree->child2->node1 == NULL
              ? currTree->child2->node2->scalarValue
              : currTree->child2->node2 == NULL
                  ? currTree->child2->node1->scalarValue
                  : (currTree->child2->node1->scalarValue
                       * currTree->child2->node1->freq
                     + currTree->child2->node2->scalarValue)
                      / childNode->freq;

      } else {

        if(currTree->child2->node1 == NULL)
          childNode->scalarValue = currTree->child2->node2->scalarValue;
        else if(currTree->child2->node2 == NULL)
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
      if(currTree->child2->node1 != NULL)
        childNode->nodeRefs.insert(childNode->nodeRefs.end(),
                                   currTree->child2->node1->nodeRefs.begin(),
                                   currTree->child2->node1->nodeRefs.end());
      if(currTree->child2->node2 != NULL)
        childNode->nodeRefs.push_back(
          std::make_pair((int)contourtrees.size() - 1,
                         currTree->child2->node2->nodeRefs[0].second));

      AlignmentEdge *childEdge = new AlignmentEdge();

      if(alignmenttreeType == lastMatchedValue) {

        childEdge->area = currTree->child2->node2 == NULL
                            ? currTree->child2->node1->area
                            : currTree->child2->node2->area;

        childEdge->scalardistance
          = currTree->child2->node2 == NULL
              ? currTree->child2->node1->scalardistanceParent
              : currTree->child2->node2->scalardistanceParent;

      } else if(alignmenttreeType == averageValues) {

        childEdge->area = currTree->child2->node1 == NULL
                            ? currTree->child2->node2->area
                            : currTree->child2->node2 == NULL
                                ? currTree->child2->node1->area
                                : (currTree->child2->node1->area
                                     * currTree->child2->node1->freq
                                   + currTree->child2->node2->area)
                                    / childNode->freq;

        childEdge->scalardistance
          = currTree->child2->node1 == NULL
              ? currTree->child2->node2->scalardistanceParent
              : currTree->child2->node2 == NULL
                  ? currTree->child2->node1->scalardistanceParent
                  : (currTree->child2->node1->scalardistanceParent
                       * currTree->child2->node1->freq
                     + currTree->child2->node2->scalardistanceParent)
                      / childNode->freq;

      } else {

        if(currTree->child2->node1 == NULL) {
          childEdge->area = currTree->child2->node2->area;
          childEdge->scalardistance
            = currTree->child2->node2->scalardistanceParent;
        } else if(currTree->child2->node2 == NULL) {
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

      childEdge->arcRefs = std::vector<std::pair<int, int>>();
      if(currTree->child2->node1 != NULL) {
        for(AlignmentEdge *e : openEdgesOld2) {
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
      if(currTree->child2->node2 != NULL) {
        for(AlignmentEdge *e : openEdgesNew1) {
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

///=====================================================================================================================
/// getters and setters
///=====================================================================================================================

int ttk::ContourTreeAlignment::getAlignmentRootIdx() {
  int idx = 0;
  for(auto node : nodes) {
    if(node == alignmentRoot)
      break;
    idx++;
  }
  return idx;
}

std::vector<ContourTree *> ttk::ContourTreeAlignment::getContourTrees() {

  /*std::vector<ContourForests*> trees_ttk;

  for(ContourTree* ct : contourtrees){
      trees_ttk.push_back(ct->getContourForest());
  }

  return trees_ttk;*/

  return contourtrees;
}

std::vector<std::pair<std::vector<CTNode *>, std::vector<CTEdge *>>>
  ttk::ContourTreeAlignment::getGraphs() {

  std::vector<std::pair<std::vector<CTNode *>, std::vector<CTEdge *>>>
    trees_simplified;

  for(ContourTree *ct : contourtrees) {
    trees_simplified.push_back(ct->getGraph());
  }

  return trees_simplified;
}

std::pair<std::vector<AlignmentNode *>, std::vector<AlignmentEdge *>>
  ttk::ContourTreeAlignment::getAlignmentGraph() {

  return std::make_pair(nodes, arcs);
}

BinaryTree *ttk::ContourTreeAlignment::getAlignmentGraphRooted() {

  AlignmentNode *root = nodes[0];
  float maxScalar = FLT_MIN;
  for(AlignmentNode *node : nodes) {

    // if(node->type==maxNode && node->edgeList[0]->scalardistance > maxRange){
    if(node->scalarValue > maxScalar) {
      root = node;
      maxScalar = node->scalarValue;
    }
  }

  return rootAtNode(root);
}

///=====================================================================================================================
/// aligning two trees
///=====================================================================================================================

std::pair<float, AlignmentTree *>
  ttk::ContourTreeAlignment::getAlignmentBinary(BinaryTree *t1,
                                                BinaryTree *t2) {

  // initialize memoization tables

  float **memT = new float *[t1->size + 1];
  float **memF = new float *[t1->size + 1];

  for(int i = 0; i <= t1->size; i++) {
    memT[i] = new float[t2->size + 1];
    memF[i] = new float[t2->size + 1];
  }

  for(int i = 0; i <= t1->size; i++) {
    for(int j = 0; j <= t2->size; j++) {
      memT[i][j] = -1;
      memF[i][j] = -1;
    }
  }

  // compute table of distances
  float dist = alignTreeBinary(t1, t2, memT, memF);

  // backtrace through the table to get the alignment
  AlignmentTree *res = traceAlignmentTree(t1, t2, memT, memF);

  for(int i = 0; i <= t1->size; i++) {
    delete[] memT[i];
    delete[] memF[i];
  }

  delete[] memT;
  delete[] memF;

  return std::make_pair(dist, res);
}

float ttk::ContourTreeAlignment::alignTreeBinary(BinaryTree *t1,
                                                 BinaryTree *t2,
                                                 float **memT,
                                                 float **memF) {

  // base cases for matching to empty tree

  if(t1 == NULL && t2 == NULL) {
    if(memT[0][0] < 0) {
      memT[0][0] = 0;
    }
    return memT[0][0];
  }

  else if(t1 == NULL) {
    if(memT[0][t2->id] < 0) {
      memT[0][t2->id]
        = editCost(NULL, t2) + alignForestBinary(NULL, t2, memT, memF);
    }
    return memT[0][t2->id];
  }

  else if(t2 == NULL) {
    if(memT[t1->id][0] < 0) {
      memT[t1->id][0]
        = editCost(t1, NULL) + alignForestBinary(t1, NULL, memT, memF);
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
        memT[t1->id][t2->id] = std::min(
          editCost(t1, NULL) + alignTreeBinary(t1->child2, NULL, memT, memF)
            + alignTreeBinary(t1->child1, t2, memT, memF),
          memT[t1->id][t2->id]);
      if(t1->size > 1)
        memT[t1->id][t2->id] = std::min(
          editCost(t1, NULL) + alignTreeBinary(t1->child1, NULL, memT, memF)
            + alignTreeBinary(t1->child2, t2, memT, memF),
          memT[t1->id][t2->id]);

      // match t2 to blank, one of its children to t1, the other to blank (try
      // both children)
      if(t2->size > 1)
        memT[t1->id][t2->id] = std::min(
          editCost(NULL, t2) + alignTreeBinary(NULL, t2->child2, memT, memF)
            + alignTreeBinary(t1, t2->child1, memT, memF),
          memT[t1->id][t2->id]);
      if(t2->size > 1)
        memT[t1->id][t2->id] = std::min(
          editCost(NULL, t2) + alignTreeBinary(NULL, t2->child1, memT, memF)
            + alignTreeBinary(t1, t2->child2, memT, memF),
          memT[t1->id][t2->id]);
    }
    return memT[t1->id][t2->id];
  }
}

float ttk::ContourTreeAlignment::alignForestBinary(BinaryTree *t1,
                                                   BinaryTree *t2,
                                                   float **memT,
                                                   float **memF) {

  // base cases for matching to empty tree

  if(t1 == NULL && t2 == NULL) {
    if(memF[0][0] < 0) {
      memF[0][0] = 0;
    }
    return memF[0][0];
  }

  else if(t1 == NULL) {
    if(memF[0][t2->id] < 0) {
      memF[0][t2->id] = 0;
      memF[0][t2->id] += alignTreeBinary(NULL, t2->child1, memT, memF);
      memF[0][t2->id] += alignTreeBinary(NULL, t2->child2, memT, memF);
    }
    return memF[0][t2->id];
  }

  else if(t2 == NULL) {
    if(memF[t1->id][0] < 0) {
      memF[t1->id][0] = 0;
      memF[t1->id][0] += alignTreeBinary(t1->child1, NULL, memT, memF);
      memF[t1->id][0] += alignTreeBinary(t1->child2, NULL, memT, memF);
    }
    return memF[t1->id][0];
  }

  // find optimal possible matching in other cases

  else {
    if(memF[t1->id][t2->id] < 0) {

      memF[t1->id][t2->id] = FLT_MAX;

      if(t1->child2 != NULL && t1->child2->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(t1->child2, NULL)
                       + alignForestBinary(t1->child2, t2, memT, memF)
                       + alignTreeBinary(t1->child1, NULL, memT, memF));
      if(t1->child1 != NULL && t1->child1->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(t1->child1, NULL)
                       + alignForestBinary(t1->child1, t2, memT, memF)
                       + alignTreeBinary(t1->child2, NULL, memT, memF));

      if(t2->child2 != NULL && t2->child2->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(NULL, t2->child2)
                       + alignForestBinary(t1, t2->child2, memT, memF)
                       + alignTreeBinary(NULL, t2->child1, memT, memF));
      if(t2->child1 != NULL && t2->child1->size > 1)
        memF[t1->id][t2->id]
          = std::min(memF[t1->id][t2->id],
                     editCost(NULL, t2->child1)
                       + alignForestBinary(t1, t2->child1, memT, memF)
                       + alignTreeBinary(NULL, t2->child2, memT, memF));

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

float ttk::ContourTreeAlignment::editCost(BinaryTree *t1, BinaryTree *t2) {

  float v1 = 0, v2 = 0;
  if(t1 != NULL)
    v1 = arcMatchMode == persistence
           ? t1->scalardistanceParent
           : arcMatchMode == area ? t1->area : t1->volume;
  if(t2 != NULL)
    v2 = arcMatchMode == persistence
           ? t2->scalardistanceParent
           : arcMatchMode == area ? t2->area : t2->volume;

  if(t1 == NULL && t2 == NULL)
    return 0;

  else if(t1 == NULL)
    return weightCombinatorialMatch + weightArcMatch * v2
           + weightScalarValueMatch;

  else if(t2 == NULL)
    return weightCombinatorialMatch + weightArcMatch * v1
           + weightScalarValueMatch;

  else if(t1->type == t2->type)
    return weightArcMatch * std::abs(v1 - v2)
           + weightScalarValueMatch
               * std::abs(t1->scalarValue - t2->scalarValue);

  else
    return FLT_MAX;
}

AlignmentTree *ttk::ContourTreeAlignment::traceAlignmentTree(BinaryTree *t1,
                                                             BinaryTree *t2,
                                                             float **memT,
                                                             float **memF) {

  if(t1 == NULL)
    return traceNullAlignment(t2, false);
  if(t2 == NULL)
    return traceNullAlignment(t1, true);

  auto id = [](BinaryTree *t) { return t == NULL ? 0 : t->id; };

  if(memT[t1->id][t2->id] == editCost(t1, t2) + memF[t1->id][t2->id]) {

    AlignmentTree *resNode = new AlignmentTree;

    resNode->node1 = t1;
    resNode->node2 = t2;
    resNode->child1 = NULL;
    resNode->child2 = NULL;
    resNode->height = 0;
    resNode->size = 1;

    std::vector<AlignmentTree *> resChildren
      = traceAlignmentForest(t1, t2, memT, memF);

    if(resChildren.size() > 0)
      resNode->child1 = resChildren[0];
    if(resChildren.size() > 1)
      resNode->child2 = resChildren[1];

    for(AlignmentTree *child : resChildren) {

      if(child != NULL) {
        resNode->height = std::max(resNode->height, child->height + 1);
        resNode->size += child->size;
      }
    }

    return resNode;
  }

  if(memT[t1->id][t2->id]
     == editCost(t1, NULL) + memT[id(t1->child2)][0]
          + memT[id(t1->child1)]
                [t2->id] /* && t1->type != maxNode && t1->type != minNode */) {

    AlignmentTree *resChild1 = traceAlignmentTree(t1->child1, t2, memT, memF);
    AlignmentTree *resChild2 = traceNullAlignment(t1->child2, true);
    AlignmentTree *res = new AlignmentTree();
    res->node1 = t1;
    res->node2 = NULL;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == NULL ? 0 : resChild1->size;
    res->size += resChild2 == NULL ? 1 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == NULL ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == NULL ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(t1, NULL) + memT[id(t1->child1)][0]
          + memT[id(t1->child2)]
                [t2->id] /* && t1->type != maxNode && t1->type != minNode */) {

    AlignmentTree *resChild1 = traceAlignmentTree(t1->child2, t2, memT, memF);
    AlignmentTree *resChild2 = traceNullAlignment(t1->child1, true);
    AlignmentTree *res = new AlignmentTree();
    res->node1 = t1;
    res->node2 = NULL;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == NULL ? 0 : resChild1->size;
    res->size += resChild2 == NULL ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == NULL ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == NULL ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(NULL, t2) + memT[0][id(t2->child2)]
          + memT[t1->id][id(
            t2->child1)] /* && t2->type != maxNode && t2->type != minNode */) {

    AlignmentTree *resChild1 = traceAlignmentTree(t1, t2->child1, memT, memF);
    AlignmentTree *resChild2 = traceNullAlignment(t2->child2, false);
    AlignmentTree *res = new AlignmentTree();
    res->node1 = NULL;
    res->node2 = t2;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == NULL ? 0 : resChild1->size;
    res->size += resChild2 == NULL ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == NULL ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == NULL ? 0 : resChild2->height + 1);
    return res;
  }

  if(memT[t1->id][t2->id]
     == editCost(NULL, t2) + memT[0][id(t2->child1)]
          + memT[t1->id][id(
            t2->child2)] /* && t2->type != maxNode && t2->type != minNode */) {

    AlignmentTree *resChild1 = traceAlignmentTree(t1, t2->child2, memT, memF);
    AlignmentTree *resChild2 = traceNullAlignment(t2->child1, false);
    AlignmentTree *res = new AlignmentTree();
    res->node1 = NULL;
    res->node2 = t2;
    res->height = 0;
    res->size = 1;
    res->child1 = resChild1;
    res->child2 = resChild2;
    res->size += resChild1 == NULL ? 0 : resChild1->size;
    res->size += resChild2 == NULL ? 0 : resChild2->size;
    res->height
      = std::max(res->height, resChild1 == NULL ? 0 : resChild1->height + 1);
    res->height
      = std::max(res->height, resChild2 == NULL ? 0 : resChild2->height + 1);
    return res;
  }

  return nullptr;
}

std::vector<AlignmentTree *> ttk::ContourTreeAlignment::traceAlignmentForest(
  BinaryTree *t1, BinaryTree *t2, float **memT, float **memF) {

  if(t1 == NULL && t2 == NULL)
    return std::vector<AlignmentTree *>();
  if(t1 == NULL) {
    std::vector<AlignmentTree *> res;
    if(t2->child1 != NULL)
      res.push_back(traceNullAlignment(t2->child1, false));
    if(t2->child2 != NULL)
      res.push_back(traceNullAlignment(t2->child2, false));
  }
  if(t2 == NULL) {
    std::vector<AlignmentTree *> res;
    if(t1->child1 != NULL)
      res.push_back(traceNullAlignment(t1->child1, true));
    if(t1->child2 != NULL)
      res.push_back(traceNullAlignment(t1->child2, true));
  }

  auto id = [](BinaryTree *t) { return t == NULL ? 0 : t->id; };

  if(memF[t1->id][t2->id]
     == memT[id(t1->child1)][id(t2->child1)]
          + memT[id(t1->child2)][id(t2->child2)]) {

    std::vector<AlignmentTree *> res;
    AlignmentTree *res1
      = traceAlignmentTree(t1->child1, t2->child1, memT, memF);
    if(res1 != NULL)
      res.push_back(res1);
    AlignmentTree *res2
      = traceAlignmentTree(t1->child2, t2->child2, memT, memF);
    if(res2 != NULL)
      res.push_back(res2);

    return res;
  }

  if(memF[t1->id][t2->id]
     == memT[id(t1->child1)][id(t2->child2)]
          + memT[id(t1->child2)][id(t2->child1)]) {

    std::vector<AlignmentTree *> res;
    AlignmentTree *res1
      = traceAlignmentTree(t1->child1, t2->child2, memT, memF);
    if(res1 != NULL)
      res.push_back(res1);
    AlignmentTree *res2
      = traceAlignmentTree(t1->child2, t2->child1, memT, memF);
    if(res2 != NULL)
      res.push_back(res2);

    return res;
  }

  if(memF[t1->id][t2->id]
     == editCost(t1->child1, NULL) + memF[id(t1->child1)][t2->id]
          + memT[id(t1->child2)][0]) {

    if(t1->child1 != NULL) {

      std::vector<AlignmentTree *> res;

      AlignmentTree *t = new AlignmentTree;
      t->node1 = t1->child1;
      t->node2 = NULL;

      t->child1 = NULL;
      t->child2 = NULL;

      t->size = 1;
      t->height = 0;

      std::vector<AlignmentTree *> resChildren
        = traceAlignmentForest(t1->child1, t2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(AlignmentTree *c : resChildren) {
        if(c != NULL) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t1->child2 != NULL)
        res.push_back(traceNullAlignment(t1->child2, true));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(t1->child2, NULL) + memF[id(t1->child2)][t2->id]
          + memT[id(t1->child1)][0]) {

    if(t1->child2 != NULL) {

      std::vector<AlignmentTree *> res;

      AlignmentTree *t = new AlignmentTree;
      t->node1 = t1->child2;
      t->node2 = NULL;

      t->child1 = NULL;
      t->child2 = NULL;

      t->size = 1;
      t->height = 0;

      std::vector<AlignmentTree *> resChildren
        = traceAlignmentForest(t1->child2, t2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(AlignmentTree *c : resChildren) {
        if(c != NULL) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t1->child1 != NULL)
        res.push_back(traceNullAlignment(t1->child1, true));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(NULL, t2->child1) + memF[t1->id][id(t2->child1)]
          + memT[0][id(t2->child2)]) {

    if(t2->child1 != NULL) {

      std::vector<AlignmentTree *> res;

      AlignmentTree *t = new AlignmentTree;
      t->node1 = NULL;
      t->node2 = t2->child1;

      t->child1 = NULL;
      t->child2 = NULL;

      t->size = 1;
      t->height = 0;

      std::vector<AlignmentTree *> resChildren
        = traceAlignmentForest(t1, t2->child1, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(AlignmentTree *c : resChildren) {
        if(c != NULL) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t2->child2 != NULL)
        res.push_back(traceNullAlignment(t2->child2, false));

      return res;
    }
  }

  if(memF[t1->id][t2->id]
     == editCost(NULL, t2->child2) + memF[t1->id][id(t2->child2)]
          + memT[0][id(t2->child1)]) {

    if(t2->child2 != NULL) {

      std::vector<AlignmentTree *> res;

      AlignmentTree *t = new AlignmentTree;
      t->node1 = NULL;
      t->node2 = t2->child2;

      t->child1 = NULL;
      t->child2 = NULL;

      t->size = 1;
      t->height = 0;

      std::vector<AlignmentTree *> resChildren
        = traceAlignmentForest(t1, t2->child2, memT, memF);
      if(resChildren.size() > 0)
        t->child1 = resChildren[0];
      if(resChildren.size() > 1)
        t->child2 = resChildren[1];

      for(AlignmentTree *c : resChildren) {
        if(c != NULL) {
          t->height = std::max(t->height, c->height + 1);
          t->size += c->size;
        }
      }

      res.push_back(t);
      if(t2->child1 != NULL)
        res.push_back(traceNullAlignment(t2->child1, false));

      return res;
    }
  }

  return std::vector<AlignmentTree *>();
}

AlignmentTree *ttk::ContourTreeAlignment::traceNullAlignment(BinaryTree *t,
                                                             bool first) {

  if(t == NULL)
    return NULL;
  AlignmentTree *at = new AlignmentTree();
  at->node1 = first ? t : NULL;
  at->node2 = first ? NULL : t;
  at->height = t->height;
  at->size = t->size;
  at->child1 = traceNullAlignment(t->child1, first);
  at->child2 = traceNullAlignment(t->child2, first);
  return at;
}

///=====================================================================================================================
/// helper functions
///=====================================================================================================================

bool ttk::ContourTreeAlignment::isBinary(Tree *t) {

  if(t->children.size() > 2)
    return false;
  else {
    bool res = true;
    for(Tree *c : t->children) {
      res = res & isBinary(c);
    }
    return res;
  }
}

BinaryTree *ttk::ContourTreeAlignment::rootAtNode(AlignmentNode *root) {

  int id = 1;
  BinaryTree *t = computeRootedTree(root, NULL, id);
  return t;
}

BinaryTree *ttk::ContourTreeAlignment::computeRootedTree(AlignmentNode *node,
                                                         AlignmentEdge *parent,
                                                         int &id) {

  BinaryTree *t = new BinaryTree;

  if(parent == NULL) {
    t->scalardistanceParent = 10000;
    t->area = 10000;
    t->volume = 10000;
  } else {
    t->scalardistanceParent = parent->scalardistance;
    t->area = parent->area;
    t->volume = t->area * t->scalardistanceParent;
  }
  t->freq = node->freq;
  t->type = node->type;
  t->child1 = NULL;
  t->child2 = NULL;
  t->id = id;
  t->scalarValue = node->scalarValue;
  t->nodeRefs = node->nodeRefs;
  if(parent != NULL)
    t->arcRefs = parent->arcRefs;
  else
    t->arcRefs = std::vector<std::pair<int, int>>();
  id++;

  t->size = 1;
  t->height = 0;

  std::vector<BinaryTree *> children;

  for(AlignmentEdge *edge : node->edgeList) {

    if(edge != parent) {

      BinaryTree *child = computeRootedTree(
        edge->node1 == node ? edge->node2 : edge->node1, edge, id);
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

BinaryTree *ttk::ContourTreeAlignment::computeRootedDualTree(AlignmentEdge *arc,
                                                             bool parent1,
                                                             int &id) {

  BinaryTree *t = new BinaryTree;

  t->scalardistanceParent = arc->scalardistance;
  t->area = arc->area;
  t->volume = t->area * t->scalardistanceParent;
  t->freq = arc->freq;
  t->type = arc->node1->type == maxNode || arc->node2->type == maxNode
              ? maxNode
              : arc->node1->type == minNode || arc->node2->type == minNode
                  ? minNode
                  : saddleNode;
  t->child1 = NULL;
  t->child2 = NULL;
  t->id = id;
  t->scalarValue = -1;
  t->arcRefs = arc->arcRefs;
  id++;

  t->size = 1;
  t->height = 0;

  std::vector<BinaryTree *> children;

  AlignmentNode *node = parent1 ? arc->node2 : arc->node1;

  for(AlignmentEdge *edge : node->edgeList) {

    if(edge != arc) {

      BinaryTree *child
        = computeRootedDualTree(edge, edge->node1 == node ? true : false, id);
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
