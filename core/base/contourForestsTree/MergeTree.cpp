/*
 * file: MergeTree.cpp
 * description: MergeTree processing package.
 * author: Gueunet Charles
 * date: Juin 2015
 */

#include "MergeTree.h"

using namespace std;
using namespace ttk;
using namespace cf;

// Constructors & destructors

MergeTree::MergeTree(Params *const params,
                     Scalars *const scalars,
                     TreeType type,
                     idPartition part)
  : params_(params), scalars_(scalars) {
  if(type == TreeType::Join) {
    this->setDebugMsgPrefix("JoinTree");
  } else if(type == TreeType::Split) {
    this->setDebugMsgPrefix("SplitTree");
  } else if(type == TreeType::Contour) {
    this->setDebugMsgPrefix("ContourTree");
  } else if(type == TreeType::JoinAndSplit) {
    this->setDebugMsgPrefix("SplitJoinTree");
  }
  treeData_.treeType = type;
  treeData_.partition = part;
}

MergeTree::~MergeTree() {
  // all is automatically destroyed in treedata
  // do not touch pointers
}

// }
// Process
// {

// update lately

void MergeTree::updateSegmentation() {
  auto compL
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isLower(a.first, b.first);
      };

  auto compH
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isHigher(a.first, b.first);
      };

  const idSuperArc nbArc = getNumberOfSuperArcs();
  if(treeData_.treeType == TreeType::Split) {
    for(idSuperArc sa = 0; sa < nbArc; sa++) {
      SuperArc *superArc = getSuperArc(sa);
      if(!superArc->isVisible())
        continue;

      auto *segmentation = superArc->getVertList();
      auto &segmSize = superArc->getVertSize();

      // fix
      if(!segmentation)
        continue;

      sort(segmentation, segmentation + segmSize, compH);
      for(SimplexId i = 0; i < segmSize; i++) {
        const SimplexId &vert = segmentation[i].first;
        updateCorrespondingArc(vert, sa);
      }
    }
  } else {
    for(idSuperArc sa = 0; sa < nbArc; sa++) {
      SuperArc *superArc = getSuperArc(sa);
      if(!superArc->isVisible())
        continue;

      auto *segmentation = superArc->getVertList();
      auto &segmSize = superArc->getVertSize();

      // fix
      if(!segmentation)
        continue;

      sort(segmentation, segmentation + segmSize, compL);
      for(SimplexId i = 0; i < segmSize; i++) {
        const SimplexId &vert = segmentation[i].first;
        if(!segmentation[i].second) {
          updateCorrespondingArc(vert, sa);
        }
      }
    }
  }

  const idNode nbNode = getNumberOfNodes();
  for(idNode n = 0; n < nbNode; n++) {
    if(!getNode(n)->isHidden()) {
      updateCorrespondingNode(getNode(n)->getVertexId(), n);
    }
  }
}

void MergeTree::parallelUpdateSegmentation(const bool ttkNotUsed(ct)) {
  // REMOVE THIS FUNCTION AND USE BOOL
  auto compL
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isLower(a.first, b.first);
      };

  auto compH
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isHigher(a.first, b.first);
      };

  const idSuperArc nbArc = getNumberOfSuperArcs();
  if(treeData_.treeType == TreeType::Split) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
    for(idSuperArc sa = 0; sa < nbArc; sa++) {
      SuperArc *superArc = getSuperArc(sa);
      if(!superArc->isVisible())
        continue;

      auto *segmentation = superArc->getVertList();
      auto &segmSize = superArc->getVertSize();

      // fix
      if(!segmentation)
        continue;

      sort(segmentation, segmentation + segmSize, compH);
      for(SimplexId i = 0; i < segmSize; i++) {
        const SimplexId &vert = segmentation[i].first;
        updateCorrespondingArc(vert, sa);
      }
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
    for(idSuperArc sa = 0; sa < nbArc; sa++) {
      SuperArc *superArc = getSuperArc(sa);
      if(!superArc->isVisible())
        continue;

      auto *segmentation = superArc->getVertList();
      auto &segmSize = superArc->getVertSize();

      // fix
      if(!segmentation)
        continue;

      sort(segmentation, segmentation + segmSize, compL);
      for(SimplexId i = 0; i < segmSize; i++) {
        const SimplexId &vert = segmentation[i].first;
        if(!segmentation[i].second) {
          updateCorrespondingArc(vert, sa);
        }
      }
    }
  }

  const idNode nbNode = getNumberOfNodes();
  for(idNode n = 0; n < nbNode; n++) {
    if(!getNode(n)->isHidden()) {
      updateCorrespondingNode(getNode(n)->getVertexId(), n);
    }
  }
}

void MergeTree::parallelInitNodeValence(const int nbThreadValence) {
  // cout << "SENTINEL : Parallel Init Node Valence " << endl;
  const auto &nbNodes = getNumberOfNodes();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nbThreadValence)
#endif
  for(idNode n = 0; n < nbNodes; n++) {
    short downVal = 0, upVal = 0;
    Node *node = getNode(n);

    const auto nbDown = node->getNumberOfDownSuperArcs();
    for(idSuperArc i = 0; i < nbDown; i++) {
      const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
      if(sa->isVisible()) {
        ++downVal;
      }
    }

    const auto nbUp = node->getNumberOfUpSuperArcs();
    for(idSuperArc i = 0; i < nbUp; i++) {
      const SuperArc *sa = getSuperArc(node->getUpSuperArcId(i));
      if(sa->isVisible()) {
        ++upVal;
      }
    }

    node->setDownValence(downVal);
    node->setUpValence(upVal);
  }
  TTK_FORCE_USE(nbThreadValence);
}

// }
// Arcs and node manipulations
// {

// SuperArcs
idSuperArc MergeTree::openSuperArc(const idNode &downNodeId,
                                   const bool overlapB,
                                   const bool overlapA) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(downNodeId >= getNumberOfNodes()) {
    cout << "[Merge Tree] openSuperArc on a inexisting node !" << endl;
    return -2;
  }
#endif

  idSuperArc newSuperArcId = treeData_.superArcs.size();
  treeData_.superArcs.emplace_back(downNodeId, nullNodes, overlapB, overlapA,
                                   treeData_.partition, treeData_.partition);
  treeData_.nodes[downNodeId].addUpSuperArcId(newSuperArcId);
  treeData_.nodes[downNodeId].incUpValence();

  return newSuperArcId;
}

idSuperArc MergeTree::makeSuperArc(const idNode &downNodeId,
                                   const idNode &upNodeId,
                                   const bool overlapB,
                                   const bool overlapA,
                                   pair<SimplexId, bool> *vertexList,
                                   SimplexId vertexSize) {
  idSuperArc newSuperArcId = treeData_.superArcs.size();

  if(downNodeId != upNodeId) {
    treeData_.superArcs.emplace_back(downNodeId, upNodeId, overlapB, overlapA,
                                     treeData_.partition, treeData_.partition);
  } else {
    // arc on a 1 vertex partition ...
    treeData_.superArcs.emplace_back(downNodeId, upNodeId, overlapB, overlapA,
                                     treeData_.partition, treeData_.partition,
                                     ComponentState::Hidden);
  }

  treeData_.superArcs[newSuperArcId].setVertList(vertexList);
  treeData_.superArcs[newSuperArcId].setVertSize(vertexSize);

  treeData_.nodes[downNodeId].addUpSuperArcId(newSuperArcId);
  treeData_.nodes[downNodeId].incUpValence();
  treeData_.nodes[upNodeId].addDownSuperArcId(newSuperArcId);
  treeData_.nodes[upNodeId].incDownValence();

  return newSuperArcId;
}

void MergeTree::closeSuperArc(const idSuperArc &superArcId,
                              const idNode &upNodeId,
                              const bool overlapB,
                              const bool overlapA) {
#ifndef TTK_ENABLE_KAMIKAZE

  if(superArcId >= getNumberOfSuperArcs()) {
    cout << "[Merge Tree] closeSuperArc on a inexisting arc !" << endl;
    return;
  }

  if(upNodeId >= getNumberOfNodes()) {
    cout << "[Merge Tree] closeOpenedArc on a inexisting node !" << endl;
    return;
  }

#endif
  treeData_.superArcs[superArcId].setUpNodeId(upNodeId);
  // TODO why do we need to re-set last-visited ? (maybe Saddle, check it out)
  treeData_.superArcs[superArcId].setLastVisited(
    getNode(upNodeId)->getVertexId());
  treeData_.nodes[upNodeId].addDownSuperArcId(superArcId);
  treeData_.nodes[upNodeId].incDownValence();

  treeData_.superArcs[superArcId].setOverlapBelow(
    treeData_.superArcs[superArcId].getOverlapBelow() != overlapB);

  if(treeData_.superArcs[superArcId].getOverlapBelow()) {
    treeData_.arcsCrossingBelow.emplace_back(superArcId);
  }

  treeData_.superArcs[superArcId].setOverlapAbove(
    treeData_.superArcs[superArcId].getOverlapAbove() != overlapA);

  if(treeData_.superArcs[superArcId].getOverlapAbove()) {
    treeData_.arcsCrossingAbove.emplace_back(superArcId);
  }
}

SimplexId MergeTree::insertNodeAboveSeed(const idSuperArc &arc,
                                         const pair<SimplexId, bool> &seed) {
  auto isLowerComp
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isLower(a.first, b.first);
      };

  SuperArc *crossing = getSuperArc(arc);
  SimplexId stitchVert;
  // get stitching node
  const auto &vertList = crossing->getVertList();
  const auto &vertSize = crossing->getVertSize();

  if(vertSize) {
    auto posVert
      = lower_bound(vertList, vertList + vertSize, seed, isLowerComp);
    while(posVert < vertList + vertSize && posVert->second) {
      ++posVert;
    }

    if(posVert == vertList + vertSize) {
      // up node
      stitchVert = getNode(crossing->getUpNodeId())->getVertexId();
    } else {
      // update segmentation
      crossing->setVertSize(posVert - vertList);
      stitchVert = posVert->first;
      // need to insert node
      const idNode &newNodeId = makeNode(stitchVert);
      Node *newNode = getNode(newNodeId);
      // for the instert node to works
      updateCorrespondingArc(stitchVert, arc);
      insertNode(newNode, false);
      newNode->setUpValence(1);
      newNode->setDownValence(1);
    }
  } else {
    // no segmentation => up node
    stitchVert = getNode(crossing->getUpNodeId())->getVertexId();
  }

  return stitchVert;
}

SimplexId MergeTree::getVertBelowSeed(const idSuperArc &arc,
                                      const pair<SimplexId, bool> &seed,
                                      const vector<idCorresp> &vert2treeOther) {
  auto isLowerComp
    = [&](const pair<SimplexId, bool> &a, const pair<SimplexId, bool> &b) {
        return isLower(a.first, b.first);
      };

  SuperArc *crossing = getSuperArc(arc);
  SimplexId stitchVert;
  // get stitching node
  const auto &vertList = crossing->getVertList();
  const auto &vertSize = crossing->getVertSize();

  if(vertSize) {
    auto posVert
      = lower_bound(vertList, vertList + vertSize, seed, isLowerComp);
    if(posVert == vertList) {
      stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
    } else {
      --posVert; // we want below and not hidden (TODO remove nullCorresp ?)
      while(
        posVert > vertList
        && (posVert->second || vert2treeOther[posVert->first] == nullCorresp)) {
        --posVert;
      }

      if(posVert == vertList) {
        stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
      } else {
        // we only find, do not touch the segmentation
        stitchVert = posVert->first;
      }
    }
  } else {
    stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
  }

  return stitchVert;
}

void MergeTree::removeHiddenDownArcs(const idNode &n) {
  Node *node = getNode(n);

  // need to check number of doan each time as it can change
  for(idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); i++) {
    const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
    if(!sa->isVisible()) {
      node->removeDownSuperArcPos(i);
      // when remove switch position with last
      --i;
    }
  }
}

void MergeTree::removeInternalDownArcs(const idNode &n) {
  Node *node = getNode(n);

  // need to check number of doan each time as it can change
  for(idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); i++) {
    const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
    if(!sa->isExternal()) {
      node->removeDownSuperArcPos(i);
      // when remove switch position with last
      --i;
    }
  }
}

idSuperArc MergeTree::getNumberOfVisibleArcs(const idNode &n) {

  Node *node = getNode(n);

  const auto nbDown = node->getNumberOfDownSuperArcs();
  const auto nbUp = node->getNumberOfUpSuperArcs();
  idSuperArc res = 0;
  for(idSuperArc i = 0; i < nbDown; i++) {
    const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
    if(sa->isVisible()) {
      ++res;
    }
  }
  for(idSuperArc i = 0; i < nbUp; i++) {
    const SuperArc *sa = getSuperArc(node->getUpSuperArcId(i));
    if(sa->isVisible()) {
      ++res;
    }
  }

  return res;
}

idSuperArc MergeTree::getNumberOfUnmergedDownArcs(const idNode &n) {

  Node *node = getNode(n);

  const auto nbDown = node->getNumberOfDownSuperArcs();
  idSuperArc res = 0;
  for(idSuperArc i = 0; i < nbDown; i++) {
    const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
    if(!sa->isMerged()) {
      ++res;
    }
  }

  return res;
}

idSuperArc MergeTree::getNumberOfExternalDownArcs(const idNode &n) {
  Node *node = getNode(n);

  const auto nbDown = node->getNumberOfDownSuperArcs();
  idSuperArc res = 0;
  for(idSuperArc i = 0; i < nbDown; i++) {
    const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
    if(sa->isExternal()) {
      ++res;
    }
  }

  return res;
}

bool MergeTree::alreadyExtLinked(const idNode &node,
                                 const idPartition &tree,
                                 const idNode &treeNode) {
  Node *n = getNode(node);

  const auto nbUp = n->getNumberOfUpSuperArcs();
  for(idSuperArc a = 0; a < nbUp; a++) {
    const idSuperArc sa = n->getUpSuperArcId(a);

    if(getSuperArc(sa)->getUpCT() == tree
       && getSuperArc(sa)->getUpNodeId() == treeNode) {
      return true;
    }
  }

  return false;
}

// state

void MergeTree::hideArc(const idSuperArc &sa) {
  treeData_.superArcs[sa].hide();
  treeData_.nodes[treeData_.superArcs[sa].getUpNodeId()].removeDownSuperArc(sa);
  treeData_.nodes[treeData_.superArcs[sa].getUpNodeId()].decDownValence();
  treeData_.nodes[treeData_.superArcs[sa].getDownNodeId()].removeUpSuperArc(sa);
  treeData_.nodes[treeData_.superArcs[sa].getDownNodeId()].decUpValence();
}

void MergeTree::mergeArc(const idSuperArc &sa,
                         const idSuperArc &recept,
                         const bool changeConnectivity) {
  treeData_.superArcs[sa].merge(recept);

  if(changeConnectivity) {
    treeData_.nodes[treeData_.superArcs[sa].getUpNodeId()].removeDownSuperArc(
      sa);
    treeData_.nodes[treeData_.superArcs[sa].getUpNodeId()].decDownValence();
    treeData_.nodes[treeData_.superArcs[sa].getDownNodeId()].removeUpSuperArc(
      sa);
    treeData_.nodes[treeData_.superArcs[sa].getDownNodeId()].decUpValence();
  }
}

void MergeTree::hideNode(const idNode &node) {
  treeData_.nodes[node].hide();
}

// Nodes

idNode MergeTree::makeNode(const SimplexId &vertexId, const SimplexId &term) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 || vertexId >= scalars_->size) {
    cout << "[Merge Tree] make node, wrong vertex :" << vertexId << " on "
         << scalars_->size << endl;
    return -1;
  }
#endif

  if(isCorrespondingNode(vertexId)) {
    return getCorrespondingNodeId(vertexId);
  }

  SimplexId size_base = (SimplexId)treeData_.nodes.size();
  treeData_.nodes.emplace_back(vertexId, term);
  updateCorrespondingNode(vertexId, size_base);

  return size_base;
}

idNode MergeTree::makeNode(const Node *const n, const SimplexId &term) {
  return makeNode(n->getVertexId(), term);
}

void MergeTree::delNode(const idNode &node,
                        const pair<SimplexId, bool> *markVertices,
                        const SimplexId &nbMark) {
  Node *mainNode = getNode(node);

  if(mainNode->getNumberOfUpSuperArcs() == 0) {

    // ----
    // Root
    // ----

#ifndef TTK_ENABLE_KAMIKAZE
    if(mainNode->getNumberOfDownSuperArcs() != 1) {
      cout << endl << "[MergeTree]:delNode won't delete ";
      cout << mainNode->getVertexId() << " (root) with ";
      cout << static_cast<unsigned>(mainNode->getNumberOfDownSuperArcs())
           << " down ";
      cout << static_cast<unsigned>(mainNode->getNumberOfUpSuperArcs())
           << " up ";
      cout << " partition : " << static_cast<unsigned>(treeData_.partition)
           << endl;
      return;
    }
#endif

    idSuperArc downArc = mainNode->getDownSuperArcId(0);
    Node *downNode = getNode(treeData_.superArcs[downArc].getDownNodeId());
    downNode->removeUpSuperArc(downArc);
    mainNode->clearDownSuperArcs();
    // USELESS : TODO remove
    treeData_.superArcs[downArc].hide();

  } else if(mainNode->getNumberOfDownSuperArcs() < 2) {

    // -------
    // Regular
    // -------

    // We delete the upArc,
    // if there is a down arc, we reattach it to the upNode

    idSuperArc upArc = mainNode->getUpSuperArcId(0);
    idNode upId = treeData_.superArcs[upArc].getUpNodeId();
    Node *upNode = getNode(upId);

    upNode->removeDownSuperArc(upArc);
    treeData_.superArcs[upArc].hide();
    mainNode->clearUpSuperArcs();
    // Have child(s)
    // Should be 0 or 1, verify
    if(mainNode->getNumberOfDownSuperArcs() == 1) {
      const idSuperArc &downArc = mainNode->getDownSuperArcId(0);

      // if have segmentation to process
      if(markVertices != nullptr) {
        // In case the two segmenation are already contiguous,
        // it means we are removing a regular node that was inserted in the
        // tree only for the combinaison.
        if((treeData_.superArcs[downArc].getVertList()
            + treeData_.superArcs[downArc].getVertSize())
           == treeData_.superArcs[upArc].getVertList()) {
          // down arc <- union of the two list
          treeData_.superArcs[downArc].setVertSize(
            treeData_.superArcs[downArc].getVertSize()
            + treeData_.superArcs[upArc].getVertSize());
          // mark removed ones (passed as parameter of this function)
          // the markVertices is sorted in reverse order as it come form the
          // other tree
          SimplexId acc = -1;
          for(SimplexId i = nbMark - 1; i >= 0; --i) {
            if(!treeData_.superArcs[downArc].getVertSize())
              break;
            while(treeData_.superArcs[downArc].getRegularNodeId(++acc)
                  != markVertices[i].first) {
              if(acc == treeData_.superArcs[downArc].getVertSize())
                break;
            }
            if(acc == treeData_.superArcs[downArc].getVertSize())
              break;

            treeData_.superArcs[downArc].setMasqued(acc);
          }

        } else {
          // we are removing a regular node of the tree who is a degenerate node
          // in the final tree. we need to keep : SEGMENTATION DETAILS HERE
          const auto &upSize = treeData_.superArcs[upArc].getVertSize();
          const auto &downSize = treeData_.superArcs[downArc].getVertSize();

          const auto *upSegm = treeData_.superArcs[upArc].getVertList();
          const auto *downSegm = treeData_.superArcs[downArc].getVertList();

          pair<SimplexId, bool> *newSegmentation
            = new pair<SimplexId, bool>[upSize + downSize];

          for(SimplexId i = 0; i < downSize; i++) {
            newSegmentation[i] = downSegm[i];
          }

          for(SimplexId i = 0; i < upSize; i++) {
            newSegmentation[i + downSize] = upSegm[i];
          }

          // avoid some memory leaks
          if(treeData_.superArcs[downArc].getSegmentation().size()) {
            const auto &downVect
              = treeData_.superArcs[downArc].getSegmentation().data();
            if(downSegm < downVect || downSegm >= downVect + downSize) {
              delete[] downSegm;
            }
          }

          if(treeData_.superArcs[upArc].getSegmentation().size()) {
            const auto &upVect
              = treeData_.superArcs[upArc].getSegmentation().data();
            if(upSegm < upVect || upSegm >= upVect + upSize) {
              delete[] upSegm;
            }
          }

          treeData_.superArcs[downArc].setVertList(newSegmentation);
          treeData_.superArcs[downArc].setVertSize(downSize + upSize);
        }
      }

      treeData_.superArcs[downArc].setUpNodeId(upId);
      treeData_.superArcs[upArc].setDownNodeId(
        getSuperArc(downArc)->getDownNodeId());
      upNode->addDownSuperArcId(downArc);
    }

    mainNode->clearDownSuperArcs();
  }
}

// Normal insert : existing arc stay below inserted (JT example)
//  *   - <- upNodeId
//  | \ |   <- newSA
//  |   * <- newNodeId
//  |   |   <- currentSA
//  - - -
idSuperArc MergeTree::insertNode(Node *node, const bool segment) {
  // already present
  if(isCorrespondingNode(node->getVertexId())) {
    Node *myNode = vertex2Node(node->getVertexId());
    // If it has been hidden / replaced we need to re-make it
    if(myNode->isHidden()) {
      SuperArc *sa = getSuperArc(myNode->getUpSuperArcId(0));
      const idSuperArc &correspondingArcId
        = (sa->getReplacantArcId() == nullSuperArc) ? myNode->getUpSuperArcId(0)
                                                    : sa->getReplacantArcId();
      updateCorrespondingArc(myNode->getVertexId(), correspondingArcId);
    } else {
      return nullSuperArc;
    }
  }

  idNode upNodeId, newNodeId;
  idSuperArc currentSA, newSA;
  SimplexId origin;

  // Create new node
  currentSA = getCorrespondingSuperArcId(node->getVertexId());
  upNodeId = treeData_.superArcs[currentSA].getUpNodeId();
  origin = treeData_.nodes[treeData_.superArcs[currentSA].getDownNodeId()]
             .getOrigin();
  newNodeId = makeNode(node, origin);

  // Connectivity
  // TODO use makeSuperArc for the newArc
  // Insert only node inside the partition : created arc don t cross
  newSA = openSuperArc(newNodeId, false, false);

  treeData_.superArcs[newSA].setUpNodeId(upNodeId);
  treeData_.nodes[upNodeId].removeDownSuperArc(currentSA);
  treeData_.nodes[upNodeId].addDownSuperArcId(newSA);

  treeData_.superArcs[currentSA].setUpNodeId(newNodeId);
  treeData_.nodes[newNodeId].addDownSuperArcId(currentSA);

  if(segment) {
    // cut the vertex list at the node position and
    // give each arc its part.
    SuperArc *tmpSA = getSuperArc(currentSA);
    pair<SimplexId, bool> *newNodePosPtr
      = (treeData_.treeType == TreeType::Split)
          ? lower_bound(tmpSA->getVertList(),
                        tmpSA->getVertList() + tmpSA->getVertSize(),
                        make_pair(node->getVertexId(), false),
                        [&](const pair<SimplexId, bool> &a,
                            const pair<SimplexId, bool> &b) {
                          return isHigher(a.first, b.first);
                        })
          : lower_bound(tmpSA->getVertList(),
                        tmpSA->getVertList() + tmpSA->getVertSize(),
                        make_pair(node->getVertexId(), false),
                        [&](const pair<SimplexId, bool> &a,
                            const pair<SimplexId, bool> &b) {
                          return isLower(a.first, b.first);
                        });

    SimplexId newNodePos = newNodePosPtr - tmpSA->getVertList();

    getSuperArc(newSA)->setVertList(newNodePosPtr);
    getSuperArc(newSA)->setVertSize(tmpSA->getVertSize() - newNodePos);

    tmpSA->setVertSize(newNodePos);
  }

  return newSA;
}

// Reverse insert : existing arc stay above inserted (JT example)
//  *   - <- upNodeId
//  |   |   <- currentSa
//  |   * <- newNodeId
//  | / |   <- newSA
//  - - -
idSuperArc MergeTree::reverseInsertNode(Node *node, const bool segment) {
  // already present
  if(isCorrespondingNode(node->getVertexId())) {
    Node *myNode = vertex2Node(node->getVertexId());
    // If it has been hidden / replaced we need to re-make it
    if(myNode->isHidden()) {
      cout << "reverse insert don t  deal with hidden" << endl;
    } else
      return nullSuperArc;
  }

  idNode downNodeId, newNodeId;
  idSuperArc currentSA, newSA;
  SimplexId origin;

  // Create new node
  currentSA = getCorrespondingSuperArcId(node->getVertexId());
  downNodeId = treeData_.superArcs[currentSA].getDownNodeId();
  origin = treeData_.nodes[treeData_.superArcs[currentSA].getDownNodeId()]
             .getOrigin();
  newNodeId = makeNode(node, origin);

  // Connectivity
  // TODO use makeSuperArc for the newArc
  // Insert only node inside the partition : created arc don t cross
  newSA = openSuperArc(downNodeId, false, false);

  treeData_.superArcs[newSA].setUpNodeId(newNodeId);
  treeData_.nodes[downNodeId].removeUpSuperArc(currentSA);
  treeData_.nodes[downNodeId].addUpSuperArcId(newSA);

  treeData_.superArcs[currentSA].setDownNodeId(newNodeId);
  treeData_.nodes[newNodeId].addUpSuperArcId(currentSA);

  treeData_.nodes[newNodeId].addDownSuperArcId(newSA);

  if(segment) {
    // cut the vertex list at the node position and
    // give each arc its part.
    SuperArc *tmpSA = getSuperArc(currentSA);
    pair<SimplexId, bool> *newNodePosPtr
      = (treeData_.treeType == TreeType::Split)
          ? lower_bound(tmpSA->getVertList(),
                        tmpSA->getVertList() + tmpSA->getVertSize(),
                        make_pair(node->getVertexId(), false),
                        [&](const pair<SimplexId, bool> &a,
                            const pair<SimplexId, bool> &b) {
                          return isHigher(a.first, b.first);
                        })
          : lower_bound(tmpSA->getVertList(),
                        tmpSA->getVertList() + tmpSA->getVertSize(),
                        make_pair(node->getVertexId(), false),
                        [&](const pair<SimplexId, bool> &a,
                            const pair<SimplexId, bool> &b) {
                          return isLower(a.first, b.first);
                        });

    SimplexId newNodePos = newNodePosPtr - tmpSA->getVertList();

    getSuperArc(newSA)->setVertList(tmpSA->getVertList());
    getSuperArc(newSA)->setVertSize(newNodePos);

    tmpSA->setVertList(newNodePosPtr);
    tmpSA->setVertSize(tmpSA->getVertSize() - newNodePos);
  }

  return newSA;
}

// traverse

inline Node *MergeTree::getDownNode(const SuperArc *a) {
  return &(treeData_.nodes[a->getDownNodeId()]);
}

inline Node *MergeTree::getUpNode(const SuperArc *a) {
  return &(treeData_.nodes[a->getUpNodeId()]);
}

idNode MergeTree::getParent(const idNode &n) {
  return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
}

// Here the return of the vector use the move constructor
vector<idNode> MergeTree::getNodeNeighbors(const idNode &n) {
  Node *node = getNode(n);
  auto nbUp = node->getNumberOfUpSuperArcs();
  auto nbDown = node->getNumberOfDownSuperArcs();

  vector<idNode> res(nbDown + nbUp);

  SimplexId currentPos = 0;

  // Nodes below
  for(idSuperArc i = 0; i < nbDown; i++) {
    const idSuperArc &corArc = node->getDownSuperArcId(i);
    const idNode &corNode = getSuperArc(corArc)->getDownNodeId();
    res[currentPos++] = corNode;
  }

  // Nodes above
  for(idSuperArc i = 0; i < nbUp; i++) {
    const idSuperArc &corArc = node->getUpSuperArcId(i);
    const idNode &corNode = getSuperArc(corArc)->getUpNodeId();
    res[currentPos++] = corNode;
  }

  return res;
}

vector<idNode> MergeTree::getNodeUpNeighbors(const idNode &n) {
  Node *node = getNode(n);
  auto nbUp = node->getNumberOfUpSuperArcs();

  vector<idNode> res(nbUp);

  // Nodes above
  for(idSuperArc i = 0; i < nbUp; i++) {
    const idSuperArc &corArc = node->getUpSuperArcId(i);
    const idNode &corNode = getSuperArc(corArc)->getUpNodeId();
    res[i] = corNode;
  }

  return res;
}

vector<idNode> MergeTree::getNodeDownNeighbors(const idNode &n) {
  Node *node = getNode(n);
  auto nbDown = node->getNumberOfDownSuperArcs();

  vector<idNode> res(nbDown);

  // Nodes above
  for(idNode i = 0; i < nbDown; i++) {
    const idSuperArc &corArc = node->getDownSuperArcId(i);
    const idNode &corNode = getSuperArc(corArc)->getDownNodeId();
    res[i] = corNode;
  }

  return res;
}

// hide / clear

void MergeTree::hideAndClearArcsAbove(const idNode &baseNode) {
  const idSuperArc nbArc = getNode(baseNode)->getNumberOfUpSuperArcs();
  for(idSuperArc i = 0; i < nbArc; ++i) {
    const idSuperArc &upsaid = getNode(baseNode)->getUpSuperArcId(i);
    SuperArc *curArc = getSuperArc(upsaid);
    if(curArc->getUpCT() != treeData_.partition) {
      continue;
    }

    hideArc(upsaid);
    hideNode(curArc->getUpNodeId());
  }

  getNode(baseNode)->clearUpSuperArcs();
}

void MergeTree::hideAndClearArcsBelow(const idNode &baseNode,
                                      const SimplexId &seed) {
  const idSuperArc nbArc = getNode(baseNode)->getNumberOfDownSuperArcs();
  for(idSuperArc i = 0; i < nbArc; ++i) {
    const idSuperArc &downsaid = getNode(baseNode)->getDownSuperArcId(i);
    SuperArc *curArc = getSuperArc(downsaid);
    if(curArc->getDownCT() != treeData_.partition) {
      continue;
    }

    Node *downNode = getNode(curArc->getDownNodeId());

    if(isLower(downNode->getVertexId(), seed)) {
      hideArc(downsaid);
      hideNode(curArc->getDownNodeId());
    }
  }

  removeHiddenDownArcs(baseNode);
}

// add a from parameter
idSuperArc MergeTree::hideAndClearLeadingTo(const idNode &baseNode,
                                            const SimplexId &v) {
  if(isCorrespondingNode(v)) {
    const auto nbDown = getNode(baseNode)->getNumberOfDownSuperArcs();
    for(idSuperArc aid = 0; aid < nbDown; aid++) {
      const idSuperArc a = getNode(baseNode)->getDownSuperArcId(aid);
      if(getSuperArc(a)->getDownCT() != treeData_.partition
         || !getSuperArc(a)->isVisible()) {
        continue;
      }

      const idNode downNode = getSuperArc(a)->getDownNodeId();

      if(getNode(downNode)->getVertexId() == v) {
        hideArc(a);
        return a;
      }
    }
  } else {
    if(isCorrespondingArc(v)) {
      idSuperArc a = getCorrespondingSuperArcId(v);
      idNode p = getSuperArc(a)->getUpNodeId();
      while(p != baseNode && getNode(p)->getUpValence()) {
        // SHOULD HAVE ONLY ONE ARC
        if(getNode(p)->getUpValence() != 1)
          cout << "Noise with up valence ! (hide&clear Leading to)" << endl;

        a = getNode(p)->getUpSuperArcId(0);
        p = getSuperArc(a)->getUpNodeId();
      }
      hideArc(a);
      return a;
    }
  }

  return nullSuperArc;
}

// }

// Operators : find, print & clone
// {

// Print
void MergeTree::printTree2() {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
  {
    cout << "Partition : " << static_cast<unsigned>(treeData_.partition)
         << endl;

    cout << "Nodes----------" << endl;
    for(idNode nid = 0; nid < getNumberOfNodes(); nid++) {
      const Node &n = treeData_.nodes[nid];
      if(n.isVisible()) {
        cout << printNode(nid) << endl;
      }
    }

    cout << "Arcs-----------" << endl;
    for(idSuperArc said = 0; said < getNumberOfSuperArcs(); ++said) {
      const SuperArc &sa = treeData_.superArcs[said];
      if(sa.isVisible()) {
        cout << printArc(said) << endl;
      }
    }

    cout << "Leaves" << endl;
    for(const auto &l : treeData_.leaves)
      cout << " " << treeData_.nodes[l].getVertexId();
    cout << endl;

    cout << "Roots" << endl;
    for(const auto &r : treeData_.roots)
      cout << " " << treeData_.nodes[r].getVertexId();
    cout << endl;
  }
}

// Clone
MergeTree *MergeTree::clone() const {
  MergeTree *newMT
    = new MergeTree(params_, scalars_, treeData_.treeType, treeData_.partition);

  newMT->treeData_.superArcs = treeData_.superArcs;
  newMT->treeData_.nodes = treeData_.nodes;
  newMT->treeData_.leaves = treeData_.leaves;
  newMT->treeData_.roots = treeData_.roots;
  newMT->treeData_.arcsCrossingBelow = treeData_.arcsCrossingBelow;
  newMT->treeData_.arcsCrossingAbove = treeData_.arcsCrossingAbove;
  newMT->treeData_.vert2tree = treeData_.vert2tree;

  return newMT;
}

void MergeTree::clone(const MergeTree *mt) {
  // we already have common data
  treeData_.superArcs = mt->treeData_.superArcs;
  treeData_.nodes = mt->treeData_.nodes;
  treeData_.leaves = mt->treeData_.leaves;
  treeData_.roots = mt->treeData_.roots;
  treeData_.arcsCrossingBelow = mt->treeData_.arcsCrossingBelow;
  treeData_.arcsCrossingAbove = mt->treeData_.arcsCrossingAbove;
  treeData_.vert2tree = mt->treeData_.vert2tree;
}

void MergeTree::doSwap(MergeTree *mt) {
  // we already have common data
  treeData_.superArcs.swap(mt->treeData_.superArcs);
  treeData_.nodes.swap(mt->treeData_.nodes);
  treeData_.leaves.swap(mt->treeData_.leaves);
  treeData_.roots.swap(mt->treeData_.roots);
  treeData_.arcsCrossingBelow.swap(mt->treeData_.arcsCrossingBelow);
  treeData_.arcsCrossingAbove.swap(mt->treeData_.arcsCrossingAbove);
  treeData_.vert2tree.swap(mt->treeData_.vert2tree);
}

// }
// Simplification
// {

// Hide the basenode and its upArc, merging the segmentation with masterArc
void MergeTree::hideAndMerge(const idSuperArc &mergingArcId,
                             const idSuperArc &receptacleArcId,
                             const bool preserveNode) {
  SuperArc *mergingArc = getSuperArc(mergingArcId);
  SuperArc *receptacleArc = getSuperArc(receptacleArcId);

  // merge     HERE WE LOST THE SORTED SEGM
  if(mergingArc->getVertSize() != -1) {
    receptacleArc->appendSegmentation(mergingArc->getSegmentation());
  }

  if(!preserveNode) {
    hideNode(mergingArc->getDownNodeId());
  }

  mergeArc(mergingArcId, receptacleArcId);
}

void MergeTree::markThisArc(vector<ExtendedUnionFind *> &ufArray,
                            const idNode &curNodeId,
                            const idSuperArc &mergingArcId,
                            const idNode &parentNodeId) {
  // size of this subtree segmentation + segmentation of this arc
  const auto &curSegmenSize = getSuperArc(mergingArcId)->getVertSize()
                              + ufArray[curNodeId]->find()->getOrigin() + 2;
  // +2 for the merging nodes

  // UF propagation
  if(ufArray[parentNodeId] == nullptr) {
    // Parent have never been seen : recopy UF
    ufArray[parentNodeId] = ufArray[curNodeId]->find();
    ufArray[parentNodeId]->find()->setOrigin(curSegmenSize);
    // cout << "will merge " << getNode(curNodeId)->getVertexId() << endl;
  } else {
    // The parent have already been visited : merge UF and segmentation
    const auto &oldSegmentationSize
      = ufArray[parentNodeId]->find()->getOrigin();
    ExtendedUnionFind::makeUnion(
      ufArray[curNodeId]->find(), ufArray[parentNodeId]->find())
      ->setOrigin(oldSegmentationSize + curSegmenSize);
    // cout << "Union on " << getNode(parentNodeId)->getVertexId();
    // cout << " from " << getNode(curNodeId)->getVertexId() << endl;
  }

  // The last parentNode is the root of the subtree
  // cout << "for " << getNode(curNodeId)->getVertexId() << " set root " <<
  // getNode(parentNodeId)->getVertexId() << endl;
  ufArray[parentNodeId]->find()->setData(-((ufDataType)parentNodeId) - 1);
}

idSuperArc MergeTree::newUpArc(const idNode &curNodeId,
                               vector<ExtendedUnionFind *> &ufArray) {

  idSuperArc keepArc = nullSuperArc;
  const auto nbUp = getNode(curNodeId)->getNumberOfUpSuperArcs();
  for(idSuperArc d = 0; d < nbUp; d++) {
    const idSuperArc &curArc = getNode(curNodeId)->getUpSuperArcId(d);
    if(!getSuperArc(curArc)->isVisible())
      continue;

    keepArc = curArc;

    const idNode &newUp = getSuperArc(curArc)->getUpNodeId();
    if(!ufArray[newUp]
       || ufArray[curNodeId]->find() != ufArray[newUp]->find()) {
      return curArc;
    }
  }

  // cout << "node " << getNode(curNodeId)->getVertexId() << " have no up arc to
  // take" << endl;
  return keepArc;
}

idSuperArc MergeTree::newDownArc(const idNode &curNodeId,
                                 vector<ExtendedUnionFind *> &ufArray) {

  idSuperArc keepArc = nullSuperArc;
  const auto nbDown = getNode(curNodeId)->getNumberOfDownSuperArcs();
  for(idSuperArc d = 0; d < nbDown; d++) {
    const idSuperArc &curArc = getNode(curNodeId)->getDownSuperArcId(d);
    if(!getSuperArc(curArc)->isVisible())
      continue;

    keepArc = curArc;

    const idNode &newDown = getSuperArc(curArc)->getDownNodeId();
    if(!ufArray[newDown]
       || ufArray[curNodeId]->find() != ufArray[newDown]->find()) {
      return curArc;
    }
  }

  // cout << "node " << printNode(curNodeId) << " have no down arc to take" <<
  // endl;
  return keepArc;
}

tuple<idNode, idNode, SimplexId> MergeTree::createReceptArc(
  const idNode &root,
  const idSuperArc &receptacleArcId,
  vector<ExtendedUnionFind *> &ufArray,
  const vector<pair<idSuperArc, idSuperArc>> &valenceOffsets) {

  const bool DEBUG = false;

  ExtendedUnionFind *ufRoot = ufArray[root]->find();
  idNode downNode = root;
  idNode upNode = root;

  if(DEBUG) {
    cout << " create receptarc for root : " << printNode(root) << endl;
    cout << " custom valence : " << valenceOffsets[root].first;
    cout << " + " << valenceOffsets[root].second << endl;
  }

  // descend in the tree until valence is not 2
  SimplexId segmentationSize = ufRoot->find()->getOrigin();
  // cout << "init size " << segmentationSize << endl;

  // We need a valence of 2 (we don't want to cross a futur saddle
  // But we want to avoid up && down = root

  while(getNode(downNode)->getUpValence() - valenceOffsets[downNode].second == 1
        && getNode(downNode)->getDownValence() - valenceOffsets[downNode].first
             == 1) {
    // take the down node not leading to the current subtree
    // if have an UF, merge with current subtree
    // (else init it?)
    const idSuperArc &downArc = newDownArc(downNode, ufArray);

    // deal with arc segmentation
    segmentationSize += getSuperArc(downArc)->getVertSize() + 2;
    // + 2 for merging nodes

    downNode = getSuperArc(downArc)->getDownNodeId();
    const idNode tmpUp = getSuperArc(downArc)->getUpNodeId();

    if(DEBUG) {
      cout << "change down to " << getNode(downNode)->getVertexId() << endl;
      cout << " new segmentation : " << segmentationSize << endl;
    }

    // UF
    if(ufArray[downNode]) {
      segmentationSize += ufArray[downNode]->find()->getOrigin();
      // ExtendedUnionFind::makeUnion(ufArray[downNode], ufRoot);
      if(ufArray[downNode]->find()->getData() < 0) {
        ufArray[downNode]->find()->setData(receptacleArcId);
      }
    } else {
      // ufArray[downNode] = ufRoot->find();
    }
    mergeArc(downArc, receptacleArcId, false);
    hideNode(tmpUp);
  }

  if(DEBUG) {
    cout << " continue receptarc for root : " << printNode(root) << endl;
    cout << " custom valence : " << valenceOffsets[root].first;
    cout << " + " << valenceOffsets[root].second << endl;
  }

  // for a node to be regular, it must have a down valence = 1 but
  // but we process the down SO ADAPT HERE
  while(getNode(upNode)->getUpValence() - valenceOffsets[upNode].second == 1
        && getNode(upNode)->getDownValence() - valenceOffsets[upNode].first
             == 1) {

    const idSuperArc &upArc = newUpArc(upNode, ufArray);

    segmentationSize += getSuperArc(upArc)->getVertSize() + 2;

    upNode = getSuperArc(upArc)->getUpNodeId();
    const idNode tmpDown = getSuperArc(upArc)->getDownNodeId();

    if(DEBUG) {
      cout << "change up to " << getNode(upNode)->getVertexId() << endl;
      cout << " new segmentation : " << segmentationSize << endl;
    }

    if(ufArray[upNode]) {
      segmentationSize += ufArray[upNode]->find()->getOrigin();
      // ExtendedUnionFind::makeUnion(ufArray[upNode], ufRoot);
      if(ufArray[upNode]->find()->getData() < 0) {
        ufArray[upNode]->find()->setData(receptacleArcId);
      }
    } else {
      // ufArray[upNode] = ufRoot->find();
    }
    mergeArc(upArc, receptacleArcId, false);
    hideNode(tmpDown);
  }

  // if upNode == downNode, take one none merging arc randomly
  // (this case is possbile if several degenerate node are following)
  if(upNode == downNode) {
    // several degen. nodes adjacent
    // Prefer down for JT / ST
    idSuperArc tmpDown = newDownArc(downNode, ufArray);
    idSuperArc tmpUp = newUpArc(upNode, ufArray);

    if(tmpDown == nullSuperArc) {
      upNode = getSuperArc(tmpUp)->getUpNodeId();
      if(ufArray[upNode]) {
        segmentationSize += ufArray[upNode]->find()->getOrigin();
        // ExtendedUnionFind::makeUnion(ufArray[downNode], ufRoot);
      } else {
        ufArray[upNode] = ufRoot->find();
      }
      getSuperArc(tmpUp)->merge(receptacleArcId);
    } else {
      downNode = getSuperArc(tmpDown)->getDownNodeId();
      if(ufArray[downNode]) {
        segmentationSize += ufArray[downNode]->find()->getOrigin();
        // ExtendedUnionFind::makeUnion(ufArray[upNode], ufRoot);
      } else {
        ufArray[downNode] = ufRoot->find();
      }
      getSuperArc(tmpDown)->merge(receptacleArcId);
    }

    if(tmpDown != nullSuperArc)
      segmentationSize += getSuperArc(tmpDown)->getVertSize() + 2;
    if(tmpUp != nullSuperArc)
      segmentationSize += getSuperArc(tmpUp)->getVertSize() + 2;

    // cout << " special : new segmentation : " << segmentationSize << endl;
  }

  return make_tuple(downNode, upNode, segmentationSize);
}

// }
// ---------------------- Debug
// {

// }

// Operators

std::ostream &ttk::cf::operator<<(std::ostream &o, ttk::cf::SuperArc const &a) {
  o << a.getDownNodeId() << " <>> " << a.getUpNodeId();
  return o;
}

std::ostream &ttk::cf::operator<<(std::ostream &o, ttk::cf::Node const &n) {
  o << n.getNumberOfDownSuperArcs() << " .-. " << n.getNumberOfUpSuperArcs();
  return o;
}
