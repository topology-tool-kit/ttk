/*
 * file:                  ContourForestsTree.cpp
 * description:           ContourForestsTree processing package.
 * author:                Gueunet Charles
 * date:                  Aout 2015
 */

#include <iterator>
#include <list>

#include "ContourForests.h"

using namespace std;
using namespace ttk;
using namespace cf;

Interface::Interface(const SimplexId &seed) : seed_(seed) {
}

// ------------------------- ContourForests

ContourForests::ContourForests()
  : ContourForestsTree(new Params(), nullptr, new Scalars()), parallelParams_(),
    parallelData_() {
  params_->treeType = TreeType::Contour;
  stringstream msg;
  msg << "[ContourForests]: DEPRECATED This module will be removed in a future"
      << "release, please use FTM instead for contour trees"
      << " and FTR for Reeb graphs." << endl;
  dMsg(cerr, msg.str(), timeMsg);
}

ContourForests::~ContourForests() {
  delete params_;
  delete scalars_;
}

// Get
// {

idPartition ContourForests::vertex2partition(const SimplexId &v) {
  const SimplexId &position = scalars_->mirrorVertices[v];
  idPartition partition = 0;
  while(
    partition < parallelParams_.nbInterfaces
    && scalars_->mirrorVertices[parallelData_.interfaces[partition].getSeed()]
         <= position) {
    ++partition;
  }

  return partition;
}

// }
// Init
// {

void ContourForests::initInterfaces() {
  // We have nbThread_ partition of the same size through all vertices
  size_t partitionSize = scalars_->size / parallelParams_.nbPartitions;

  // ------------------
  // Seeds
  // ------------------
  // {

  // We initiate interface with their seed (isovalue) and their adjacent
  // partition
  //  and each partition with it size and bounds.
  for(idInterface i = 0; i < parallelParams_.nbInterfaces; ++i) {
    // interfaces have their first vertex of the sorted array as seed
    parallelData_.interfaces.emplace_back(
      scalars_->sortedVertices[partitionSize * (i + 1)]);
  }

  // }
  // ------------------
  // Print Debug
  // ------------------
  // {

  if(params_->debugLevel >= 4) {
    stringstream partition;
    partition << "seeds :";
    for(const auto &i : parallelData_.interfaces) {
      partition << i.getSeed() << " ";
    }
    partition << endl;
    dMsg(cout, partition.str(), 3);
  }

  // }
}

void ContourForests::initOverlap() {
  const SimplexId nbEdges = mesh_->getNumberOfEdges();

  // if we choose to have less partition, we still want to use all thread for
  // overlap init.

  // ------------------
  // Parallel find border vertices
  // ------------------
  // {

  vector<vector<vector<SimplexId>>> lowers(parallelParams_.nbThreads);
  vector<vector<vector<SimplexId>>> uppers(parallelParams_.nbThreads);

  for(numThread p = 0; p < parallelParams_.nbThreads; p++) {
    lowers[p].resize(parallelParams_.nbInterfaces);
    uppers[p].resize(parallelParams_.nbInterfaces);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(parallelParams_.nbThreads) schedule(static)
#endif
  for(SimplexId e = 0; e < nbEdges; e++) {

#ifdef TTK_ENABLE_OPENMP
    idPartition part = omp_get_thread_num();
#else
    idPartition part = 0;
#endif

    vector<vector<SimplexId>> &localUppers = uppers[part];
    vector<vector<SimplexId>> &localLowers = lowers[part];

    SimplexId v0, v1;
    mesh_->getEdgeVertex(e, 0, v0);
    mesh_->getEdgeVertex(e, 1, v1);

    for(idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
      const bool side0 = isEqHigher(v0, parallelData_.interfaces[i].getSeed());
      const bool side1 = isEqHigher(v1, parallelData_.interfaces[i].getSeed());

      if(side0 != side1) {
        // edge cross this interface, add both extrema in it
        if(side0) {
          // The seed is already in the partition, we do not want to have it
          // twice
          // if (v0 != vect_interfaces_[i].getSeed()) {
          localUppers[i].emplace_back(v0);
          //}
          localLowers[i].emplace_back(v1);
        } else {
          // if (v1 != vect_interfaces_[i].getSeed()) {
          localUppers[i].emplace_back(v1);
          //}
          localLowers[i].emplace_back(v0);
        }
      }
    }
  }

  // }
  // --------------------------
  // Insert in interfaces
  // --------------------------
  // {

  // reserve
  vector<SimplexId> sizeReserveUp(parallelParams_.nbInterfaces, 0);
  vector<SimplexId> sizeReserveLo(parallelParams_.nbInterfaces, 0);
  for(numThread p = 0; p < parallelParams_.nbThreads; p++) {
    for(idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
      sizeReserveUp[i] += uppers[p][i].size();
      sizeReserveLo[i] += lowers[p][i].size();
    }
  }

  for(idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
    parallelData_.interfaces[i].upReserve(sizeReserveUp[i]);
    parallelData_.interfaces[i].loReserve(sizeReserveLo[i]);
  }

  // append
  for(numThread p = 0; p < parallelParams_.nbThreads; p++) {
    for(idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
      parallelData_.interfaces[i].appendUpper(uppers[p][i]);
      parallelData_.interfaces[i].appendLower(lowers[p][i]);
    }
  }

  // }
  // -----------------
  // Sort the overlap
  // ----------------
  // {

  auto vertComp
    = [&](const SimplexId &a, const SimplexId &b) { return isLower(a, b); };

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(parallelParams_.nbThreads) schedule(static)
#endif
  for(idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
    vector<SimplexId> &upOverlap = parallelData_.interfaces[i].getUpper();
    vector<SimplexId> &loOverlap = parallelData_.interfaces[i].getLower();

    // sort & unique via set

    // LESS EFFICIENT IN PARALLEL
    // {
    // set<SimplexId, decltype(vertComp)> setUpOverlap(upOverlap.begin(),
    // upOverlap.end(), vertComp); vector<SimplexId>
    // vectUpOverlap(setUpOverlap.begin(), setUpOverlap.end());
    // parallelData_.interfaces[i].swapUpper(vectUpOverlap);

    // set<SimplexId, decltype(vertComp)> setLoOverlap(loOverlap.begin(),
    // loOverlap.end(), vertComp); vector<SimplexId>
    // vectLoOverlap(setLoOverlap.begin(), setLoOverlap.end());
    // parallelData_.interfaces[i].swapLower(vectLoOverlap);
    // }

    // sort & unique via functions

    sort(upOverlap.begin(), upOverlap.end(), vertComp);
    auto upLast = unique(upOverlap.begin(), upOverlap.end());
    upOverlap.erase(upLast, upOverlap.end());

    sort(loOverlap.begin(), loOverlap.end(), vertComp);
    auto loLast = unique(loOverlap.begin(), loOverlap.end());
    loOverlap.erase(loLast, loOverlap.end());
  }

  // }
  // -----------
  // Debug print
  // -----------
  // {

  // for (idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
  // cout << "interface : " << i << endl;

  // cout << "upper" << endl;
  // for (const SimplexId &v : parallelData_.interfaces[i].getUpper()) {
  // cout << v << ", ";
  //}

  // cout << endl << "lower" << endl;
  // for (const SimplexId &v : parallelData_.interfaces[i].getLower()) {
  // cout << v << ", ";
  //}

  // cout << endl;
  //}

  // }
}

void ContourForests::initNbPartitions() {
  if(parallelParams_.lessPartition && parallelParams_.nbThreads >= 2) {
    parallelParams_.nbPartitions = parallelParams_.nbThreads / 2;
  } else {
    parallelParams_.nbPartitions = parallelParams_.nbThreads;
  }

  parallelParams_.nbInterfaces = parallelParams_.nbPartitions - 1;
}

// }

// Process
// {

void ContourForests::stitch() {
  // We need goods informations here befor starting
  // Get the arc/NODE correspoding to the seed + crossingEdges

  if(params_->treeType == TreeType::Contour) {
    stitchTree(2);
  } else {
    stitchTree(0);
    stitchTree(1);
  }
}

void ContourForests::stitchTree(const char treetype) {

  const bool DEBUG = false;

  // For the instance assume ct :
  //    insert node with true (jt / ct ok)
  auto getTreePart = [&, treetype](const idPartition &i) -> MergeTree * {
    if(treetype == 0) {
      return parallelData_.trees[i].getJoinTree();
    }
    if(treetype == 1) {
      return parallelData_.trees[i].getSplitTree();
    }
    return &parallelData_.trees[i];
  };

  vector<bool> seenSeed(parallelParams_.nbInterfaces, false);

  // For each partition, we stich with above
  for(idPartition i = 0; i < parallelParams_.nbPartitions - 1; i++) {
    MergeTree *curTree = getTreePart(i);
    const auto &seedPair
      = make_pair(parallelData_.interfaces[i].getSeed(), false);

    if(DEBUG) {
      cout << "partition : " << static_cast<unsigned>(i);
      cout << " seed is : " << seedPair.first << endl << endl;
    }

    // For each superarc crossing the upper boundary :
    for(const auto &arc : curTree->treeData_.arcsCrossingAbove) {
      SuperArc *crossing = curTree->getSuperArc(arc);

      // Hidden arc, no need to stich (may have already been stitched)
      if(!crossing->isVisible())
        continue;

      const idNode &downCrossingId = crossing->getDownNodeId();

      // the down node of the arc is not on this partition, It should have
      // already been processed
      if(vertex2partition(curTree->getNode(downCrossingId)->getVertexId()) < i)
        continue;

      // Stitch vertex and insertion in current tree
      const SimplexId &stitchVertex
        = curTree->insertNodeAboveSeed(arc, seedPair);

      // Opposite partition
      const idPartition &otherPartition = vertex2partition(stitchVertex);
      MergeTree *otherTree = getTreePart(otherPartition);

      if(DEBUG) {
        cout << "crossing arc is " << curTree->printArc(arc) << endl;
      }

      if(DEBUG) {
        cout << "stitch vertex : " << stitchVertex << endl;
        cout << "on partition " << static_cast<unsigned>(otherPartition)
             << endl;
        cout << "crossing arc is now " << curTree->printArc(arc) << endl;
      }

      const idNode &curTreeStitchNodeId
        = curTree->getCorrespondingNodeId(stitchVertex);
      Node *curTreeStitchNode = curTree->getNode(curTreeStitchNodeId);

      bool otherTreeAlreadyHide = false;
      if(otherTree->isCorrespondingArc(stitchVertex)) {
        if(DEBUG) {
          const idSuperArc &sa
            = otherTree->getCorrespondingSuperArcId(stitchVertex);
          cout << "other tree arc is : " << otherTree->printArc(sa) << endl;
        }
        const auto &arcToHide
          = otherTree->reverseInsertNode(curTreeStitchNode, true);
        curTreeStitchNode->setDownValence(1);
        curTreeStitchNode->setUpValence(1);
        otherTree->hideArc(arcToHide);
        otherTreeAlreadyHide = true;

        if(DEBUG) {
          cout << "hide arc in other : " << otherTree->printArc(arcToHide)
               << endl;
        }
      }

      const idNode &otherTreeStitchNodeId
        = otherTree->getCorrespondingNodeId(stitchVertex);
      Node *otherTreeStitchNode = otherTree->getNode(otherTreeStitchNodeId);

      if(DEBUG) {
        cout << "Stitch nodes : " << endl;
        cout << "current : " << curTree->printNode(curTreeStitchNodeId) << endl;
        cout << "other   : " << otherTree->printNode(otherTreeStitchNodeId)
             << endl;
      }

      // Now we can remove all arc above the stitch vertex in the current tree
      // (noise) When using debug, all arc are not hidden, they are just
      // disconnected Look at the node information, not all the arcs

      curTreeStitchNode->clearUpSuperArcs();

      // for the other tree we need to replace arc (hide the one that is
      // replaced) Exeption : the seed may contain noise directly above: we
      // clear its down arcs if it is the stitch node.
      if(stitchVertex == seedPair.first) {
        if(DEBUG) {
          cout << "stitch vertex is seed" << endl;
        }
        // we have'nt clear it yet
        if(!seenSeed[otherPartition]) {
          seenSeed[otherPartition] = true;
          otherTree->removeInternalDownArcs(otherTreeStitchNodeId);
          if(DEBUG) {
            cout << "clear below stitch vert " << endl;
          }
        }
      } else if(!otherTreeAlreadyHide && crossing->getDownCT() == i) {
        const SimplexId &currentBelowSeed = curTree->getVertBelowSeed(
          arc, seedPair, otherTree->treeData_.vert2tree);
        const idSuperArc &arcToHide = otherTree->hideAndClearLeadingTo(
          otherTreeStitchNodeId, currentBelowSeed);

        if(DEBUG) {
          cout << "hide arc leading to: ";
          if(arcToHide != nullSuperArc) {
            cout << otherTree->printArc(arcToHide) << endl;
          } else {
            cout << "not found" << endl;
          }
        }
      }

      // Create an external arc to link both nodes on each tree
      // The created arc may be a new crossing arc

      curTree->treeData_.superArcs.emplace_back(curTreeStitchNodeId,
                                                otherTreeStitchNodeId, false,
                                                true, i, otherPartition);
      curTreeStitchNode->addUpSuperArcId(curTree->getNumberOfSuperArcs() - 1);

      // chk if cross the next interface to add it in the vector if crossing
      // above
      bool crossNextInterface = false;
      if(otherPartition < parallelParams_.nbInterfaces) {
        const SimplexId &nextSeed
          = parallelData_.interfaces[otherPartition].getSeed();
        crossNextInterface = isEqLower(nextSeed, stitchVertex);
      }
      otherTree->treeData_.superArcs.emplace_back(
        curTreeStitchNodeId, otherTreeStitchNodeId, true, crossNextInterface, i,
        otherPartition);
      otherTreeStitchNode->addDownSuperArcId(otherTree->getNumberOfSuperArcs()
                                             - 1);
      if(crossNextInterface) {
        otherTree->addCrossingAbove(otherTree->getNumberOfSuperArcs() - 1);
        if(DEBUG) {
          cout << "new crossing above" << endl;
        }
      }

      if(DEBUG) {
        cout << "arc added :" << endl;
        cout << "current : "
             << curTree->printArc(curTree->getNumberOfSuperArcs() - 1) << endl;
        cout << "other   : "
             << otherTree->printArc(otherTree->getNumberOfSuperArcs() - 1)
             << endl;
        cout << endl << endl;
      }

    } // for each arc of this curTree

  } // for each partition

  if(DEBUG) {
    printVectCT();
  }
}

void ContourForests::unify() {

  if(params_->treeType == TreeType::Contour) {
    unifyTree(2);
  } else {
    unifyTree(0);
    unifyTree(1);
  }
}

void ContourForests::unifyTree(const char treetype) {

  const bool DEBUG = false;

  // Get the good tree
  auto getTreePart = [&, treetype](const idPartition &i) -> MergeTree * {
    if(treetype == 0) {
      return parallelData_.trees[i].getJoinTree();
    }
    if(treetype == 1) {
      return parallelData_.trees[i].getSplitTree();
    }
    return &parallelData_.trees[i];
  };

  // this tree will receive the final tree
  // all variables linked to tmpree have a "_tt" suffix
  MergeTree tmpTree(params_, mesh_, scalars_, params_->treeType);
  // for vert2tree
  tmpTree.flush();
  // statistical reserves
  tmpTree.treeData_.nodes.reserve(scalars_->size / 50);
  tmpTree.treeData_.superArcs.reserve(scalars_->size / 50);

  // partition, node in partion, is a leaf
  queue<tuple<idInterface, idNode>> leavesNodes;
  vector<unsigned> nbVisit(scalars_->size, 0);

  // Unify by traversing leaves for each partition starting by the lowest
  for(idPartition partition = 0; partition < parallelParams_.nbPartitions;
      ++partition) {
    MergeTree *currentTree = getTreePart(partition);

    for(auto &l : currentTree->treeData_.leaves) {
      Node *curNode = currentTree->getNode(l);
      SimplexId leafVert = curNode->getVertexId();

      // Condition to keep

      // if not in partition
      if(partition != 0
         && isLower(
              leafVert, parallelData_.interfaces[partition - 1].getSeed()))
        continue;
      if(partition != parallelParams_.nbInterfaces
         && isHigher(leafVert, parallelData_.interfaces[partition].getSeed()))
        continue;

      // if hidden
      if(!curNode->isVisible())
        continue;

      // if max
      if(!curNode->getNumberOfUpSuperArcs())
        continue;

      // this leaf have received an external arc during stitching, not a real
      // leaf anymore
      if(curNode->getNumberOfUpSuperArcs()
         && curNode->getNumberOfDownSuperArcs())
        continue;

      // isolated node
      if(!currentTree->getNumberOfVisibleArcs(l))
        continue;

      // Add the leave

      if(!nbVisit[leafVert]) {
        // will be processed
        leavesNodes.emplace(partition, l);

        if(DEBUG) {
          cout << "will see : partition : " << static_cast<unsigned>(partition);
          cout << " leaf node " << currentTree->printNode(l) << endl;
        }

        // seen once
        ++nbVisit[leafVert];

        // Create node in tmpTree and mark as leaf
        const idNode &newNodeId_tt = tmpTree.makeNode(curNode);
        tmpTree.treeData_.leaves.emplace_back(newNodeId_tt);
      }
    }
  }

  // cross node from min to max to construct the tree

  while(!leavesNodes.empty()) {
    // get the next node
    idPartition currentPartition;
    idNode currentNodeId;
    tie(currentPartition, currentNodeId) = leavesNodes.front();
    leavesNodes.pop();

    // get tree and node
    MergeTree *currentTree = getTreePart(currentPartition);
    Node *currentNode = currentTree->getNode(currentNodeId);

    if(DEBUG) {
      cout << endl;
      cout << "process : partition : "
           << static_cast<unsigned>(currentPartition) << endl;
      cout << " node " << currentTree->printNode(currentNodeId) << endl;
    }

    // create or recover in tmpTree
    const idNode &baseNode_tt = tmpTree.makeNode(currentNode);

    // Cross the ups arc of the node
    const idPartition refPartition = currentPartition;
    const idNode refNodeId = currentNodeId;
    const idSuperArc &nbUpArc = currentNode->getNumberOfUpSuperArcs();

    for(idSuperArc upArcPos = 0; upArcPos < nbUpArc; upArcPos++) {
      // reset partition / node and tree
      currentTree = getTreePart(refPartition);
      currentNodeId = refNodeId;
      currentNode = currentTree->getNode(currentNodeId);

      const idSuperArc &upArcId = currentNode->getUpSuperArcId(upArcPos);
      SuperArc *upArc = currentTree->getSuperArc(upArcId);

      if(DEBUG) {
        cout << " process arc " << currentTree->printArc(upArcId) << endl;
      }

      if(!upArc->isVisible()) {
        if(DEBUG) {
          cout << " - ignore not visible" << endl;
        }
        continue;
      }

      const idSuperArc &newArcId_tt
        = tmpTree.openSuperArc(baseNode_tt, false, false);
      SuperArc *newArc_tt = tmpTree.getSuperArc(newArcId_tt);

      // segmentation related
      list<pair<SimplexId, bool> *> listVertList;
      list<SimplexId> listVertSize;
      int totalSize = 0;

      // add upArc until the upnode is not a regular one
      do {
        listVertList.push_back(upArc->getVertList());
        listVertSize.push_back(upArc->getVertSize());
        totalSize += upArc->getVertSize();

        upArc->hide();

        currentPartition = upArc->getUpCT();
        currentTree = getTreePart(currentPartition);
        currentNodeId = upArc->getUpNodeId();
        currentNode = currentTree->getNode(currentNodeId);

        if(currentNode->getNumberOfUpSuperArcs() == 1
           && currentNode->getNumberOfDownSuperArcs() == 1) {
          upArc = currentTree->getSuperArc(currentNode->getUpSuperArcId(0));
          if(DEBUG) {
            cout << " cross : "
                 << currentTree->printArc(currentNode->getUpSuperArcId(0));
            cout << endl;
          }
        } else {
          // no longer regular : stop here
          if(DEBUG) {
            cout << "stop at " << currentTree->printNode(currentNodeId) << endl;
          }
          break;
        }

      } while(true);

      // Finish the current Arc (segmentation + close)
      if(totalSize) {
        newArc_tt->appendVertLists(listVertList, listVertSize, totalSize);
      }
      const idNode &closingNode_tt = tmpTree.makeNode(currentNode);
      tmpTree.closeSuperArc(newArcId_tt, closingNode_tt, false, false);

      if(DEBUG) {
        cout << " Create arc : " << tmpTree.printArc(newArcId_tt) << endl;
      }

      // push current vertex TODO
      const SimplexId &closingVertex = currentNode->getVertexId();
      ++nbVisit[closingVertex];

      // get the down valence of this node in the partition where it have no
      // noise
      const idPartition &closingPartition = vertex2partition(closingVertex);
      MergeTree *closingTree = getTreePart(closingPartition);
      const idNode &closingNodeId
        = closingTree->getCorrespondingNodeId(closingVertex);
      const unsigned &downVal
        = closingTree->getNode(closingNodeId)->getNumberOfDownSuperArcs();

      if(nbVisit[closingVertex] == downVal) {
        leavesNodes.emplace(closingPartition, closingNodeId);
        if(DEBUG) {
          cout << " push : partition : "
               << static_cast<unsigned>(closingPartition) << endl;
          cout << " push : node : " << closingTree->printNode(closingNodeId)
               << endl;
        }
      } else if(DEBUG) {
        cout << " visit : " << nbVisit[closingVertex] << endl;
        cout << " downVal : " << downVal << endl;
      }
    } // end for each up arc
  } // end while leavesNodes

  tmpTree.treeData_.superArcs.shrink_to_fit();
  tmpTree.treeData_.nodes.shrink_to_fit();

  // could use swap, more efficient
  if(treetype == 0) {
    jt_->clone(&tmpTree);
  } else if(treetype == 1) {
    st_->clone(&tmpTree);
  } else if(treetype == 2) {
    clone(&tmpTree);
  }
}
// }

// Print
// {
void ContourForests::printDebug(DebugTimer &timer, const string &str) {
  stringstream msg;
  msg << "[ContourForests] " << str << " : " << timer.getElapsedTime() << endl;
  dMsg(cout, msg.str(), timeMsg);
}

void ContourForests::printVectCT() {
  int arcCTUp, arcCTDown;

  for(idPartition nb = 0; nb < parallelParams_.nbPartitions; ++nb) {
    cout << "CT " << nb << endl;
    cout << "Nodes" << endl;

    for(const auto &n : parallelData_.trees[nb].getNodes()) {
      if(!n.isHidden()) {
        cout << "Node  " << n.getVertexId();

        if(n.isHidden())
          cout << " X ";

        cout << endl;
        cout << "  arc up : ";

        for(idSuperArc i = 0; i < n.getNumberOfUpSuperArcs(); ++i) {
          cout << n.getUpSuperArcId(i) << " ";
        }

        cout << endl << " arc down : ";

        for(idSuperArc i = 0; i < n.getNumberOfDownSuperArcs(); ++i) {
          cout << n.getDownSuperArcId(i) << " ";
        }

        cout << endl;
      }
    }

    cout << "Arcs" << endl;

    for(const auto &sa : parallelData_.trees[nb].getSuperArc()) {
      if(!sa.isHidden()) {
        arcCTDown = sa.getDownCT();
        arcCTUp = sa.getUpCT();

        if(sa.getDownNodeId() == nullNodes) {
          cout << "||";
        } else {
          cout << static_cast<unsigned>(arcCTDown) << ":";
          cout << parallelData_.trees[arcCTDown]
                    .getNode(sa.getDownNodeId())
                    ->getVertexId();
        }

        if(sa.isHidden())
          cout << " <X> ";
        else if(!sa.isVisible())
          cout << " <-> ";
        else
          cout << " <> ";

        if(sa.getUpNodeId() == nullNodes) {
          cout << "||";
        } else {
          cout << static_cast<unsigned>(arcCTUp) << ":";
          cout << parallelData_.trees[arcCTUp]
                    .getNode(sa.getUpNodeId())
                    ->getVertexId();
        }

        cout << endl;
      }
    }

    if(1) {
      cout << "Leaves" << endl;

      for(const auto &l : parallelData_.trees[nb].getLeaves())
        cout << " " << l;

      cout << endl;

      cout << "Roots" << endl;

      for(const auto &r : parallelData_.trees[nb].getRoots())
        cout << " " << r;

      cout << endl;
    }
  }
}
// }
