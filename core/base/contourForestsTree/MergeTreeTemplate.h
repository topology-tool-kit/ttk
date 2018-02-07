/// \ingroup base
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the
/// sublevel set tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef MERGETREETEMPLATE_H
#define MERGETREETEMPLATE_H

#include "MergeTree.h"

// Init
// {

template <typename scalarType>
void MergeTree::sortInput(void)
{
   const auto &nbVertices = scalars_->size;
   auto &      sortedVect = scalars_->sortedVertices;

   if (!sortedVect.size()) {
      auto indirect_sort = [&](const size_t &a, const size_t &b) {
         return isLower<scalarType>(a, b);
      };

      sortedVect.resize(nbVertices, 0);
      iota(sortedVect.begin(), sortedVect.end(), 0);

#ifdef TTK_ENABLE_OPENMP
      __gnu_parallel::sort(sortedVect.begin(), sortedVect.end(), indirect_sort);
#else
      sort(sortedVect.begin(), sortedVect.end(), indirect_sort);
#endif
   }

   if (!scalars_->mirrorVertices.size()) {
      scalars_->mirrorVertices.resize(nbVertices);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for (idVertex i = 0; i < nbVertices; i++) {
         scalars_->mirrorVertices[sortedVect[i]] = i;
      }
   }
}
// }

// Process
// {

// Simplify

template <typename scalarType>
idEdge MergeTree::localSimplify(const idVertex &posSeed0, const idVertex &posSeed1)
{

    // if null threshold, leave
    if (!params_->simplifyThreshold) {
        return 0;
    }

    const bool DEBUG = false;

    // -----------------
    // Persistance pairs
    // -----------------
    // {

    vector<tuple<idVertex, idVertex, scalarType, bool>> pairs;
    computePersistencePairs<scalarType>(pairs);

    if (DEBUG) {
       cout << "pairs : ( threshold : " << params_->simplifyThreshold << "  )" << endl;
       for (const auto &p : pairs) {
          const idVertex &  thisOriginVert = get<0>(p);
          const idVertex &  thisEndVert    = get<1>(p);
          const scalarType &thisPersist    = get<2>(p);
          const bool        thisNeedUp     = get<3>(p);
          cout << thisOriginVert << " - " << thisEndVert << " : " << thisPersist;
          cout << " ( " << thisNeedUp << " )" << endl;
       }
       cout << endl;
    }

    // }
    // --------------
    // Simplify
    // --------------
    // {

    return simplifyTree<scalarType>(posSeed0, posSeed1, pairs);

    // }
}

template <typename scalarType>
idEdge MergeTree::globalSimplify(const idVertex posSeed0, const idVertex posSeed1)
{

    // if null threshold, leave
    if (!params_->simplifyThreshold) {
        return 0;
    }

    //---------------------
    // Sort Nodes
    //---------------------
    //{

    auto isLowerComp = [&](const idNode &n1, const idNode &n2) {
        return isLower(getNode(n1)->getVertexId(), getNode(n2)->getVertexId());
    };

    const auto nbNode = getNumberOfNodes();

    vector<idNode> sortedNodes(nbNode);
    iota(sortedNodes.begin(), sortedNodes.end(), 0);
// Sort nodes by vertex scalar
//{
#ifdef TTK_ENABLE_OPENMP
    __gnu_parallel::sort(sortedNodes.begin(), sortedNodes.end(), isLowerComp);
#else
    sort(sortedNodes.begin(), sortedNodes.end(), isLowerComp);
#endif
//}

    //}
    //---------------------
    // Make pairs
    //---------------------
    //{

    // origin, end, persistance, needToGoUp
    vector<tuple<idVertex, idVertex, scalarType,bool>> pairsJT;
    vector<tuple<idVertex, idVertex, scalarType,bool>> pairsST;

    recoverMTPairs<scalarType>(sortedNodes, pairsJT, pairsST);

    //}
    //---------------------
    // Fusionne & Sort pairs
    //---------------------
    //{

    auto pairComp = [](const tuple<idVertex, idVertex, scalarType, bool> &a,
                       const tuple<idVertex, idVertex, scalarType, bool> &b) {
       // sort by persistence
       return get<2>(a) < get<2>(b);
    };

    vector<tuple<idVertex, idVertex, scalarType, bool>> sortedPairs;
    size_t sizePJT = pairsJT.size(), sizePST = pairsST.size();

    sortedPairs.reserve(sizePJT + sizePST);

    sortedPairs.insert(sortedPairs.end(), pairsJT.begin(), pairsJT.end());
    sortedPairs.insert(sortedPairs.end(), pairsST.begin(), pairsST.end());

    // Sort pairs by persistence
    //{
    // IS SET STILL BETTER ? (parallel sort) TODO

#ifdef TTK_ENABLE_OPENMP
        __gnu_parallel::sort(sortedPairs.begin(), sortedPairs.end(), pairComp);
#else
        sort(sortedPairs.begin(), sortedPairs.end(), pairComp);
#endif

    auto last = unique(sortedPairs.begin(), sortedPairs.end());
    sortedPairs.erase(last, sortedPairs.end());

    //}
    //---------------------
    // Traverse pairs and merge on the tree
    //---------------------
    //{

        // identify subtrees and merge them in recept'arcs
        return simplifyTree<scalarType>(posSeed0, posSeed1, sortedPairs);
    //}
}

template <typename scalarType>
idEdge MergeTree::simplifyTree(
    const idVertex &posSeed0, const idVertex &posSeed1,
    const vector<tuple<idVertex, idVertex, scalarType, bool>> &sortedPairs)
{
   const auto nbNode = getNumberOfNodes();
   const auto nbArcs = getNumberOfSuperArcs();
   // Retain the relation between merge coming from st, jt
   // also retain info about what we keep
   vector<ExtendedUnionFind *> subtreeUF(nbNode, nullptr);

   // nb arc seen below / above this node
   vector<pair<idSuperArc, idSuperArc>> valenceOffset(nbNode, make_pair(0,0));
   idEdge nbPairMerged = 0;

   const bool DEBUG = false;

   if (DEBUG) {
      cout << "Imapct simplify on tree btwn : " << posSeed0 << " and " << posSeed1 << endl;
   }

   //----------
   // Make subtrees
   //-----------
   //{

   queue<tuple<idNode,bool>> node2see;

   // Add the origin of all pairs that need to merge
   for (const auto & pp : sortedPairs) {
      if (get<2>(pp) < params_->simplifyThreshold ) {
         const idVertex &thisOriginVert = get<0>(pp);
         const idVertex &thisEndVert    = get<1>(pp);
         const idNode &  thisOriginId   = getCorrespondingNodeId(thisOriginVert);

         if (scalars_->mirrorVertices[thisOriginVert] <= posSeed0 ||
             scalars_->mirrorVertices[thisOriginVert] >= posSeed1 ||
             scalars_->mirrorVertices[thisEndVert] <= posSeed0 ||
             scalars_->mirrorVertices[thisEndVert] >= posSeed1) {
            continue;
         }

         node2see.emplace(thisOriginId, get<3>(pp));
         subtreeUF[thisOriginId] = new ExtendedUnionFind(0);
         ++nbPairMerged;
         if(DEBUG){
            cout << "willSee " << printNode(thisOriginId) << endl;
         }
      } else break;
   }

   //---
   // In UF :
   //  Origin is the size of the segmentation
   //  Data is negative : -idNodeRoot-1
   //  Data is positive : Receptacle Arc id
   //--
   // Use the queue to mark Arc that will be merged and UF to identify
   // subtree. When a node have only one way out : enqueue it to continue travresall
   while (!node2see.empty()) {
       idNode curNodeId;
       bool needToGoUp; // identify up/down traversall

       tie(curNodeId, needToGoUp) = node2see.front(); // should have only one arc valid to take
       node2see.pop();

       if (DEBUG) {
           cout << "process : " << printNode(curNodeId) << endl;
       }

       // Here we take the only available arc :
       idSuperArc mergingArcId;
       idNode parentNodeId;
       // continue traversall
       if(needToGoUp){
          mergingArcId = newUpArc(curNodeId, subtreeUF);
          parentNodeId = getSuperArc(mergingArcId)->getUpNodeId();
          ++valenceOffset[curNodeId].second;
          ++valenceOffset[parentNodeId].first;
       } else {
          mergingArcId = newDownArc(curNodeId, subtreeUF);
          parentNodeId = getSuperArc(mergingArcId)->getDownNodeId();
          ++valenceOffset[curNodeId].first;
          ++valenceOffset[parentNodeId].second;
       }

       markThisArc(subtreeUF, curNodeId, mergingArcId, parentNodeId);

       // if we have processed all but one arc of this node, we nee to continue traversall
       // throug it
       if (valenceOffset[parentNodeId].first + valenceOffset[parentNodeId].second + 1 ==
           getNode(parentNodeId)->getValence()) {
           // only one way out, is it up ?
           node2see.emplace(parentNodeId, valenceOffset[parentNodeId].second + 1 ==
                                              getNode(parentNodeId)->getUpValence());
           if(DEBUG){
               cout << " add to see " << printNode(parentNodeId) << endl;
           }
       }
   } // end while node2see

   // for each node valenceOffset is the number of arc attached to this node that will merge

   // Debug print
   if (DEBUG) {
      cout << "node subtrees before creating receptarc " << endl;
      for (idNode nid = 0; nid < nbNode; nid++) {
         if (subtreeUF[nid]) {
            cout << "node " << getNode(nid)->getVertexId() << " is in subtree rooted :";
            const idNode &root = -subtreeUF[nid]->find()->getData() - 1;
            cout << getNode(root)->getVertexId();
            const idVertex &segmSize = subtreeUF[nid]->find()->getOrigin();
            cout << " with segmentation of " << segmSize << endl;
         }
      }
   }

   //}
   //----------
   // Create the recept'arcs
   //-----------
   //{

   // Add the origin of all pairs that need to merge
   for (const auto & pp : sortedPairs) {
      if (get<2>(pp) < params_->simplifyThreshold ) {
         const idVertex &thisOriginVert = get<0>(pp);
         const idVertex &thisEndVert    = get<1>(pp);
         const idNode &  thisOriginId   = getCorrespondingNodeId(thisOriginVert);
         //const idNode &  thisEndId      = getCorrespondingNode(thisEndVert);

         if (scalars_->mirrorVertices[thisOriginVert] <= posSeed0 ||
             scalars_->mirrorVertices[thisOriginVert] >= posSeed1 ||
             scalars_->mirrorVertices[thisEndVert] <= posSeed0 ||
             scalars_->mirrorVertices[thisEndVert] >= posSeed1) {
            continue;
         }

         if (subtreeUF[thisOriginId]->find()->getData() < 0) {
            // create receptarc
            const idNode &subtreeRoot = -subtreeUF[thisOriginId]->find()->getData() - 1;
            // The id of the next arc to be created : NOT PARALLEL
            const idSuperArc receptArcId = treeData_.superArcs.size();
            // down , up, segmentation size
            // create the receptacle arc and merge arc not in sub-tree in it
            const tuple<idNode, idNode, idVertex> &receptArc =
                createReceptArc(subtreeRoot, receptArcId, subtreeUF, valenceOffset);

            // make superArc and do the makeAlloc on it
            const bool overlapB =
                scalars_->mirrorVertices[getNode(get<0>(receptArc))->getVertexId()] < posSeed0;
            const bool overlapA =
                scalars_->mirrorVertices[getNode(get<1>(receptArc))->getVertexId()] >= posSeed1;
            const idSuperArc na = makeSuperArc(get<0>(receptArc), get<1>(receptArc), overlapB, overlapA, nullptr, -1);

            if(overlapB){
                treeData_.arcsCrossingBelow.emplace_back(na);
            }

            if(overlapA){
                treeData_.arcsCrossingAbove.emplace_back(na);
            }

            subtreeUF[thisOriginId]->find()->setData(receptArcId);
            getSuperArc(receptArcId)->makeAllocGlobal(get<2>(receptArc));

            if (DEBUG) {
               cout << "create arc : " << printArc(receptArcId)
                    << " with segm : " << get<2>(receptArc) << endl;
            }
         }
      } else break;
   }

   //}
   //----------
   // Merge these subtrees
   //-----------
   //{

   // The merge is done after because we don't want to forget arcs parallel to the recept'arc
   // but we cannot make the difference before the former are created

   // nbArcs is before the insertion of receptarcs so they will no be crossed here
   for (idSuperArc arc = 0; arc < nbArcs; arc++) {

       if(getSuperArc(arc)->isMerged()){
          const idSuperArc &receptacleArcId = getSuperArc(arc)->getReplacantArcId();
          // take care of connectivity
          mergeArc(arc, receptacleArcId);
          if(DEBUG){
            cout << " parallel merge in " << printArc(receptacleArcId) << endl;
            cout << " arc " << printArc(arc) << " size : " << getSuperArc(arc)->getVertSize()
                 << endl;
          }
          getSuperArc(receptacleArcId)
              ->addSegmentationGlobal(getSuperArc(arc)->getVertList(),
                                      getSuperArc(arc)->getVertSize());
       } else {
          const idNode &downNode = getSuperArc(arc)->getDownNodeId();
          const idNode &upNode   = getSuperArc(arc)->getUpNodeId();

          if (!(subtreeUF[downNode] && subtreeUF[upNode]))
             continue;

          if (subtreeUF[downNode] && subtreeUF[upNode] &&
              subtreeUF[downNode]->find() != subtreeUF[upNode]->find()) {
             if (DEBUG) {
                cout << "Arc between 2 degenerate with mergin " << printArc(arc) << endl;
                cout << "below recept : " << printArc(subtreeUF[downNode]->find()->getData());
                cout << endl;
                cout << "Above recept : " << printArc(subtreeUF[upNode]->find()->getData()) << endl;
                cout << endl;
             }

             continue;
          }

          ExtendedUnionFind *curUF =
              (subtreeUF[upNode]) ? subtreeUF[upNode]->find() : subtreeUF[downNode]->find();

          const idSuperArc &receptacleArcId = curUF->getData();

          if(DEBUG){
             cout << "merge in " << printArc(receptacleArcId) << endl;
             cout << "  arc " << printArc(arc) << " size : " << getSuperArc(arc)->getVertSize();
             cout << endl;
          }

          getSuperArc(arc)->merge(receptacleArcId);
          if (getSuperArc(arc)->getVertSize()) {
             getSuperArc(receptacleArcId)
                 ->addSegmentationGlobal(getSuperArc(arc)->getVertList(),
                                         getSuperArc(arc)->getVertSize());

             getSuperArc(receptacleArcId)->addSegmentationGlobal(getNode(downNode)->getVertexId());
             getSuperArc(receptacleArcId)->addSegmentationGlobal(getNode(upNode)->getVertexId());
          }

          // Tree topology
          getNode(downNode)->removeUpSuperArc(arc);
          getNode(downNode)->decUpValence();
          getNode(upNode)->removeDownSuperArc(arc);
          getNode(upNode)->decDownValence();

          if(!getNode(downNode)->getNumberOfUpSuperArcs())
              hideNode(downNode);
          if(!getNode(upNode)->getNumberOfDownSuperArcs())
              hideNode(upNode);

       }
   } // end for arcs

   //}

        return nbPairMerged;
}

// Persistence

template <typename scalarType>
int MergeTree::computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>> &pairs)
{
#ifndef TTK_ENABLE_KAMIKAZE
   if (!getNumberOfSuperArcs()) {
      return -1;
   }
#endif

   pairs.reserve(treeData_.leaves.size());

   for (const idNode &leave : treeData_.leaves) {
      Node *   curNode  = getNode(leave);
      idVertex curVert  = curNode->getVertexId();
      idVertex termVert = getNode(curNode->getTerminaison())->getVertexId();

      addPair<scalarType>(pairs, curVert, termVert);
   }

   auto pair_sort = [](const tuple<idVertex, idVertex, scalarType> &a,
                       const tuple<idVertex, idVertex, scalarType> &b) {
      return get<2>(a) < get<2>(b);
   };

   sort(pairs.begin(), pairs.end(), pair_sort);

   return 0;
}

template <typename scalarType>
int MergeTree::computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType,bool>> &pairs)
{
    // Need to be called on MergeTree, not ContourTree

#ifndef TTK_ENABLE_KAMIKAZE
   if (!getNumberOfSuperArcs()) {
      return -1;
   }

   if (treeData_.treeType == TreeType::Contour) {
      cout << "WARNING, computePersistencePairs is made to be called on Join or Split Tree" << endl;
   }
#endif

   pairs.reserve(treeData_.leaves.size());

   for (const idNode &leave : treeData_.leaves) {
      Node *   curNode  = getNode(leave);
      idVertex curVert  = curNode->getVertexId();
      idVertex termVert = getNode(curNode->getTerminaison())->getVertexId();

      addPair<scalarType>(pairs, curVert, termVert, treeData_.treeType == TreeType::Join);
   }

   auto pair_sort = [](const tuple<idVertex, idVertex, scalarType, bool> &a,
                       const tuple<idVertex, idVertex, scalarType, bool> &b) {
      return get<2>(a) < get<2>(b);
   };

   sort(pairs.begin(), pairs.end(), pair_sort);

   return 0;
}

template <typename scalarType>
void MergeTree::recoverMTPairs(
    const vector<idNode> &sortedNodes, vector<tuple<idVertex, idVertex, scalarType, bool>> &pairsJT,
    vector<tuple<idVertex, idVertex, scalarType, bool>> &pairsST)
{
    const auto nbNode = getNumberOfNodes();

    vector<ExtendedUnionFind *> vect_JoinUF(nbNode, nullptr);
    vector<ExtendedUnionFind *> vect_SplitUF(nbNode, nullptr);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections num_threads(2)
#endif
    {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
       {
          pairsJT.reserve(treeData_.leaves.size());
          // For the biggest pair of the component
          map<idVertex, idVertex> pendingMinMax;

          for (auto it = sortedNodes.cbegin(); it != sortedNodes.cend(); ++it) {
             const idNode &  n = *it;
             const idVertex &v = getNode(n)->getVertexId();

             const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
             const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

             if (nbDown == 0) {
                // leaf
                vect_JoinUF[n] = new ExtendedUnionFind(getNode(n)->getVertexId());
                vect_JoinUF[n]->setOrigin(v);
                //cout << " jt origin : " << v << endl;
             } else {
                // first descendant
                const idSuperArc   firstSaId = getNode(n)->getDownSuperArcId(0);
                const SuperArc *   firstSA   = getSuperArc(firstSaId);
                const idNode &firstChildNodeId = firstSA->getDownNodeId();

                ExtendedUnionFind *merge     = vect_JoinUF[firstChildNodeId]->find();
                idVertex           further   = merge->getOrigin();
                idSuperArc         furtherI  = 0;

                // Find the most persistant way
                for (idSuperArc ni = 1; ni < nbDown; ++ni) {
                   const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                   const SuperArc * curSA   = getSuperArc(curSaId);
                   const idNode &   neigh   = curSA->getDownNodeId();

                   // fix
                   if (neigh == n)
                      continue;

                   ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                   if (isLower(neighUF->getOrigin(), further)) {
                      further = neighUF->getOrigin();
                      furtherI = ni;
                   }
                }

                if (nbDown > 1) {
                   // close finish pair and make union
                   for (idSuperArc ni = 0; ni < nbDown; ++ni) {
                      const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                      const SuperArc * curSA   = getSuperArc(curSaId);
                      const idNode &   neigh   = curSA->getDownNodeId();

                      // fix
                      if (neigh == n)
                         continue;

                      ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                      if (ni != furtherI) {  // keep the more persitent pair
                         addPair<scalarType>(pairsJT, neighUF->getOrigin(), v, true);
                         pendingMinMax.erase(neighUF->getOrigin());

                         // cout << " jt make pair : " << neighUF->getOrigin() << " - " << v <<
                         // endl;
                      }

                      ExtendedUnionFind::makeUnion(merge, neighUF)->setOrigin(further);
                   }
                }

                merge->find()->setOrigin(further);
                vect_JoinUF[n] = merge->find();

                if(!nbUp){
                   // potential close of the component
                   //cout << "pending for " << further << " is " << v << endl;
                   pendingMinMax[further] = v;
                }
             }
          } // end for each node

          // Add the pending biggest pair of each component
          for (const auto & pair_vert : pendingMinMax) {
              if (isCorrespondingNode(pair_vert.second)) {
                 //cout << " add : " << pair_vert.first << " - " << pair_vert.second << endl;
                 addPair<scalarType>(pairsJT, pair_vert.first, pair_vert.second, true);
              }
          }
       } // end para section

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
       {
          pairsST.reserve(treeData_.leaves.size());
          // For the biggest pair of the component
          map<idVertex, idVertex> pendingMinMax;

          for (auto it = sortedNodes.crbegin(); it != sortedNodes.crend(); ++it) {
             const idNode &  n = *it;
             const idVertex &v = getNode(n)->getVertexId();

             const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
             const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

             if (nbUp  == 0) {
                // leaf
                vect_SplitUF[n] = new ExtendedUnionFind(getNode(n)->getVertexId());
                vect_SplitUF[n]->setOrigin(v);
                //cout << " st origin : " << v << endl;
             } else {
                // first descendant
                const idSuperArc   firstSaId = getNode(n)->getUpSuperArcId(0);
                const SuperArc *   firstSA   = getSuperArc(firstSaId);
                const idNode &firstChildNodeId = firstSA->getUpNodeId();

                ExtendedUnionFind *merge     = vect_SplitUF[firstChildNodeId]->find();
                idVertex           further   = merge->getOrigin();
                idSuperArc         furtherI  = 0;

                for (idSuperArc ni = 1; ni < nbUp; ++ni) {
                   // find the more persistant way
                   const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                   const SuperArc * curSA   = getSuperArc(curSaId);
                   // Ignore hidden / simplified arc
                   if(!curSA->isVisible()) continue;
                   const idNode &   neigh   = curSA->getUpNodeId();

                   // fix
                   if (neigh == n)
                      continue;

                   //cout << "visit neighbor : " << ni << " which is " << getNode(neigh)->getVertexId() << endl;

                   ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                   if (isHigher(neighUF->getOrigin(), further)) {
                      further = neighUF->getOrigin();
                      furtherI = ni;
                   }
                }

                if (nbUp > 1) {
                   // close finsh pair and make union
                   for (idSuperArc ni = 0; ni < nbUp; ++ni) {
                      const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                      const SuperArc * curSA   = getSuperArc(curSaId);
                      const idNode &   neigh   = curSA->getUpNodeId();

                      // fix
                      if (neigh == n)
                         continue;

                      ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                      if(ni != furtherI){
                         addPair<scalarType>(pairsST, neighUF->getOrigin(), v, false);

                         pendingMinMax.erase(neighUF->getOrigin());

                         //cout << " st make pair : " << neighUF->getOrigin() << " - " << v
                              //<< " for neighbor " << getNode(neigh)->getVertexId() << endl;
                      }

                      ExtendedUnionFind::makeUnion(merge, neighUF)->setOrigin(further);
                      // Re-visit after merge lead to add the most persistant pair....
                   }
                }
                merge->find()->setOrigin(further);
                vect_SplitUF[n] = merge->find();

                if (!nbDown) {
                   pendingMinMax[further] = v;
                }

             } // end nbUp == 0 else
          } // end for each node

          // Add the pending biggest pair of each component
          for (const auto & pair_vert : pendingMinMax) {
             addPair<scalarType>(pairsST, pair_vert.first, pair_vert.second, false);
          }
       } // end para section
    } // end para
}

// }
// Tools
// {

template <typename scalarType>
void MergeTree::addPair(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs,
                        const idVertex &orig, const idVertex &term, const bool goUp)
{
   if (params_->simplifyMethod == SimplifMethod::Persist) {
      pairs.emplace_back(orig, term, fabs(getValue<scalarType>(orig) - getValue<scalarType>(term)),
                         goUp);
   } else if (params_->simplifyMethod == SimplifMethod::Span) {
      float coordOrig[3], coordTerm[3], span;
      mesh_->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
      mesh_->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
      span = Geometry::distance(coordOrig, coordTerm);
      pairs.emplace_back(orig, term, span, goUp);
   }
}

template <typename scalarType>
void MergeTree::addPair(vector<tuple<idVertex, idVertex, scalarType>> &pairs,
                        const idVertex &orig, const idVertex &term)
{
   if (params_->simplifyMethod == SimplifMethod::Persist) {
      pairs.emplace_back(orig, term, fabs(getValue<scalarType>(orig) - getValue<scalarType>(term)));
   } else if (params_->simplifyMethod == SimplifMethod::Span) {
      float coordOrig[3], coordTerm[3], span;
      mesh_->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
      mesh_->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
      span = Geometry::distance(coordOrig, coordTerm);
      pairs.emplace_back(orig, term, span);
   }
}

// }

#endif /* end of include guard: MERGETREETEMPLATE_H */

