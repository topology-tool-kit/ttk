/*
 * file: ContourForestsTree.cpp
 * description: ContourForestsTree processing package.
 * author: Gueunet Charles
 * date: Juin 2015
 */

#include <iterator>
#include <string>

#include "ContourForestsTree.h"

using namespace std;
using namespace ttk;
using namespace cf;

ContourForestsTree::ContourForestsTree(Params *const params,
                                       Scalars *const scalars,
                                       idPartition part)
  : MergeTree(params, scalars, TreeType::Contour, part),
    jt_(new MergeTree(params, scalars, TreeType::Join, part)),
    st_(new MergeTree(params, scalars, TreeType::Split, part)) {
}

ContourForestsTree::~ContourForestsTree() {
  if(jt_) {
    delete jt_;
    jt_ = nullptr;
  }
  if(st_) {
    delete st_;
    st_ = nullptr;
  }
}

// Process
// {

int ContourForestsTree::combine(const SimplexId &seed0,
                                const SimplexId &seed1) {
  queue<pair<bool, idNode>> growingNodes;
  pair<bool, idNode> head;

  MergeTree *xt = nullptr, *yt = nullptr;
  idNode correspondingNodeId, parentId, node1, node2;
  Node *currentNode, *parentNode;

  const bool DEBUG = false;

  // If a tree add only one leaf, the other tree is the result
  SimplexId nbAddedleavesST = 0, nbAddedleavesJT = 0;

  // Add leves to growing nodes
  // We insert non hidden nodes, only those of the interface or those
  // just beyond, linked to a crossing arc
  for(const idNode &nId : st_->getLeaves()) {
    if(!st_->getNode(nId)->isHidden()
       && st_->getNode(nId)->getNumberOfSuperArcs()) {
      growingNodes.emplace(false, nId);
      ++nbAddedleavesST;
    }
  }

  // filiform tree
  if(nbAddedleavesST == 1) {
    growingNodes.pop();
  }

  // count how many leaves can be added, if more than one : ok!
  for(const idNode &nId : jt_->getLeaves()) {
    if(!jt_->getNode(nId)->isHidden()
       && jt_->getNode(nId)->getNumberOfSuperArcs()) {
      ++nbAddedleavesJT;
    }
    if(nbAddedleavesJT > 1)
      break;
  }

  if(nbAddedleavesJT > 1) {
    for(const idNode &nId : jt_->getLeaves()) {
      if(!jt_->getNode(nId)->isHidden()
         && jt_->getNode(nId)->getNumberOfSuperArcs()) {
        growingNodes.emplace(true, nId);
      }
    }
  }

  if(DEBUG) {
    cout << "growingNodes : " << growingNodes.size() << endl;
  }

  if(nbAddedleavesST == 1 && nbAddedleavesJT == 1) {
    // ultra simplistic case where both tree a filliform
    clone(jt_);
    return 0;
  }

  // Warning, have a reserve here, can't make it at the begnining, need build
  // output
  treeData_.leaves.reserve(jt_->getLeaves().size() + st_->getLeaves().size());

  if(growingNodes.empty()) {
    cout << "[ContourForestsTree::combine ] Nothing to combine" << endl;
  }

  // seed : to keep crossing edges;
  const SimplexId &s0
    = (seed0 == -1) ? nullVertex : scalars_->sortedVertices[seed0];
  const SimplexId &s1
    = (seed1 >= scalars_->size) ? nullVertex : scalars_->sortedVertices[seed1];

  while(!growingNodes.empty()) {
    // i <- Get(Q)
    head = growingNodes.front();

    if(head.first) {
      // node come frome jt
      xt = jt_;
      yt = st_;
    } else {
      // node come from st
      xt = st_;
      yt = jt_;
    }

    currentNode = xt->getNode(head.second);

    if(DEBUG) {
      if(xt == jt_)
        cout << "JT ";
      else
        cout << "ST ";
      cout << "node : " << currentNode->getVertexId() << endl;
    }

    correspondingNodeId
      = yt->getCorrespondingNodeId(currentNode->getVertexId());

    if(isCorrespondingNode(currentNode->getVertexId())) {
      // already a node in the tree
      node1 = getCorrespondingNodeId(currentNode->getVertexId());
    } else {
      // create a new node
      node1 = makeNode(currentNode);

      // check if leaf
      if(!currentNode->getNumberOfDownSuperArcs()
         || !currentNode->getNumberOfUpSuperArcs())
        treeData_.leaves.emplace_back(node1);
    }

    // "choose a non-root leaf that is not a split in ST" so we ignore such
    // nodes
    if(currentNode->getNumberOfUpSuperArcs() == 0
       || yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {

      growingNodes.pop();

      // if(currentNode->getNumberOfUpSuperArcs() == 0){
      // if (DEBUG) {
      // cout << "ignore orphan" << endl;
      //}
      // continue;
      //}

      if(yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
        if(DEBUG) {
          cout << "re-enqueue and ignore " << yt->printNode(correspondingNodeId)
               << endl;
        }

        growingNodes.emplace(head.first, head.second);
        continue;
      }

      if(DEBUG) {
        cout << "  ignore" << endl;
      }

      continue;
    }

    // j <- GetAdj(XT, i)
    idSuperArc curUpArc = currentNode->getUpSuperArcId(0);
    if(xt->getSuperArc(curUpArc)->isMerged())
      curUpArc = xt->getSuperArc(curUpArc)->getReplacantArcId();
    parentId = xt->getSuperArc(curUpArc)->getUpNodeId();
    if(parentId == nullNodes) {
      // security : if not closed arc, close it here
      // can append when 1 vertex is isolated
      parentId = xt->makeNode(
        xt->getSuperArc(currentNode->getUpSuperArcId(0))->getLastVisited());
      // single not si not crossing anything
      xt->closeSuperArc(
        currentNode->getUpSuperArcId(0), parentId, false, false);
    }

    parentNode = xt->getNode(parentId);

    // cout << " parent node :" << parentNode->getVertexId() << endl;

    // HERE parent is null ...
    if(isCorrespondingNode(parentNode->getVertexId())) {
      // already a node in the tree
      node2 = getCorrespondingNodeId(parentNode->getVertexId());
    } else {
      // create a new node
      node2 = makeNode(parentNode);
      if(!parentNode->getNumberOfUpSuperArcs())
        treeData_.leaves.emplace_back(node2);
    }

    // AddArc(CT, ij)
    pair<SimplexId, bool> *arcVertList = nullptr;
    SimplexId arcVertSize = 0;
    {
      arcVertList
        = xt->getSuperArc(currentNode->getUpSuperArcId(0))->getVertList();
      arcVertSize
        = xt->getSuperArc(currentNode->getUpSuperArcId(0))->getVertSize();

      bool overlapB = false, overlapA = false;

      // If the created arc cross tha above or below interface, keep this info
      if(s0 != nullVertex) {
        overlapB = (isHigher(currentNode->getVertexId(), s0)
                    != isHigher(parentNode->getVertexId(), s0));
      }
      if(s1 != nullVertex) {
        overlapA = (isHigher(currentNode->getVertexId(), s1)
                    != isHigher(parentNode->getVertexId(), s1));
      }

      idSuperArc createdArc;
      // create the arc
      if(isLower(currentNode->getVertexId(),
                 parentNode->getVertexId())) { // take care of the order
        createdArc = makeSuperArc(
          node1, node2, overlapB, overlapA, arcVertList, arcVertSize);
      } else {
        createdArc = makeSuperArc(
          node2, node1, overlapB, overlapA, arcVertList, arcVertSize);
      }

      if(overlapB) {
        treeData_.arcsCrossingBelow.emplace_back(createdArc);
      }
      if(overlapA)
        treeData_.arcsCrossingAbove.emplace_back(createdArc);

      // Segmentation, update only if not already set, user vert2tree to check
      SimplexId nbv = 0;
      for(SimplexId vert = 0; vert < arcVertSize; vert++) {
        const SimplexId &v = arcVertList[vert].first;
        if(isCorrespondingNull(v)) {
          updateCorrespondingArc(v, createdArc);
          ++nbv;
        } else {
          getSuperArc(createdArc)->setMasqued(vert);
        }
      }

      if(DEBUG) {
        cout << " arc added : (segm: " << nbv << ") ";
        cout << printArc(createdArc) << endl;
      }
    }

    // DelNode(XT, i)
    {
      if(DEBUG) {
        cout << " delete xt (" << (xt == jt_)
             << ") node :" << xt->getNode(head.second)->getVertexId() << endl;
      }

      xt->delNode(head.second);
    }

    // DelNode(YT, i)
    {
      if(yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() < 2) {
        if(DEBUG) {
          cout << " delete yt (" << head.first << ") node :";
          cout << yt->getNode(correspondingNodeId)->getVertexId();
          cout << " have : ";
          cout << static_cast<unsigned>(
            yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs());
          cout << " down";
          cout << " and : "
               << static_cast<unsigned>(
                    yt->getNode(correspondingNodeId)->getNumberOfUpSuperArcs())
               << " up" << endl;
        }

        yt->delNode(correspondingNodeId, arcVertList, arcVertSize);
      }
    }

    if(parentNode->getNumberOfDownSuperArcs() == 0
       && parentNode->getNumberOfUpSuperArcs()) {
      growingNodes.emplace(head.first, parentId);

      if(DEBUG) {
        cout << "will see : " << parentNode->getVertexId() << endl;
      }
    }

    growingNodes.pop();
  }

  return 0;
}

//}
