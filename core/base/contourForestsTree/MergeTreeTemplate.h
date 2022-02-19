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
/// \sa ttkContourForests.cpp %for a usage example.

#pragma once

#include "MergeTree.h"

namespace ttk {
  namespace cf {
    // Init
    // {

    template <typename scalarType>
    void MergeTree::sortInput(void) {
      const auto &nbVertices = scalars_->size;

      if(scalars_->sortedVertices.size() != static_cast<size_t>(nbVertices)) {
        scalars_->sortedVertices.resize(nbVertices);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for(SimplexId i = 0; i < nbVertices; i++) {
        scalars_->sortedVertices[scalars_->sosOffsets[i]] = i;
      }
    }
    // }

    // Process
    // {

    // Simplify

    template <typename scalarType>
    SimplexId MergeTree::localSimplify(const SimplexId &posSeed0,
                                       const SimplexId &posSeed1) {

      // if null threshold, leave
      if(!params_->simplifyThreshold) {
        return 0;
      }

      const bool DEBUG = false;

      // -----------------
      // Persistance pairs
      // -----------------
      // {

      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> pairs;
      computePersistencePairs<scalarType>(pairs);

      if(DEBUG) {
        std::cout << "pairs : ( threshold : " << params_->simplifyThreshold
                  << " )" << std::endl;
        for(const auto &p : pairs) {
          const SimplexId &thisOriginVert = std::get<0>(p);
          const SimplexId &thisEndVert = std::get<1>(p);
          const scalarType &thisPersist = std::get<2>(p);
          const bool thisNeedUp = std::get<3>(p);
          std::cout << thisOriginVert << " - " << thisEndVert << " : "
                    << thisPersist;
          std::cout << " ( " << thisNeedUp << " )" << std::endl;
        }
        std::cout << std::endl;
      }

      // }
      // --------------
      // Simplify
      // --------------
      // {

      return simplifyTree<scalarType>(posSeed0, posSeed1, pairs);

      // }
    }

    template <typename scalarType, typename triangulationType>
    SimplexId MergeTree::globalSimplify(const SimplexId posSeed0,
                                        const SimplexId posSeed1,
                                        const triangulationType &mesh) {

      // if null threshold, leave
      if(!params_->simplifyThreshold) {
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

      std::vector<idNode> sortedNodes(nbNode);
      iota(sortedNodes.begin(), sortedNodes.end(), 0);
      // Sort nodes by vertex scalar
      //{
      TTK_PSORT(this->threadNumber_, sortedNodes.begin(), sortedNodes.end(),
                isLowerComp);
      //}

      //}
      //---------------------
      // Make pairs
      //---------------------
      //{

      // origin, end, persistance, needToGoUp
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> pairsJT;
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> pairsST;

      recoverMTPairs<scalarType>(sortedNodes, pairsJT, pairsST, mesh);

      //}
      //---------------------
      // Fusionne & Sort pairs
      //---------------------
      //{

      auto pairComp
        = [](const std::tuple<SimplexId, SimplexId, scalarType, bool> &a,
             const std::tuple<SimplexId, SimplexId, scalarType, bool> &b) {
            // sort by persistence
            return std::get<2>(a) < std::get<2>(b);
          };

      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>>
        sortedPairs;
      size_t sizePJT = pairsJT.size(), sizePST = pairsST.size();

      sortedPairs.reserve(sizePJT + sizePST);

      sortedPairs.insert(sortedPairs.end(), pairsJT.begin(), pairsJT.end());
      sortedPairs.insert(sortedPairs.end(), pairsST.begin(), pairsST.end());

      // Sort pairs by persistence
      //{
      // IS SET STILL BETTER ? (parallel sort) TODO
      TTK_PSORT(
        this->threadNumber_, sortedPairs.begin(), sortedPairs.end(), pairComp);

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
    SimplexId MergeTree::simplifyTree(
      const SimplexId &posSeed0,
      const SimplexId &posSeed1,
      const std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>>
        &sortedPairs) {
      const auto nbNode = getNumberOfNodes();
      const auto nbArcs = getNumberOfSuperArcs();
      // Retain the relation between merge coming from st, jt
      // also retain info about what we keep
      std::vector<ExtendedUnionFind *> subtreeUF(nbNode, nullptr);

      // nb arc seen below / above this node
      std::vector<std::pair<idSuperArc, idSuperArc>> valenceOffset(
        nbNode, std::make_pair(0, 0));
      SimplexId nbPairMerged = 0;

      const bool DEBUG = false;

      if(DEBUG) {
        std::cout << "Imapct simplify on tree btwn : " << posSeed0 << " and "
                  << posSeed1 << std::endl;
      }

      //----------
      // Make subtrees
      //-----------
      //{

      std::queue<std::tuple<idNode, bool>> node2see;

      // Add the origin of all pairs that need to merge
      for(const auto &pp : sortedPairs) {
        if(std::get<2>(pp) < params_->simplifyThreshold) {
          const SimplexId &thisOriginVert = std::get<0>(pp);
          const SimplexId &thisEndVert = std::get<1>(pp);
          const idNode &thisOriginId = getCorrespondingNodeId(thisOriginVert);

          if(scalars_->sosOffsets[thisOriginVert] <= posSeed0
             || scalars_->sosOffsets[thisOriginVert] >= posSeed1
             || scalars_->sosOffsets[thisEndVert] <= posSeed0
             || scalars_->sosOffsets[thisEndVert] >= posSeed1) {
            continue;
          }

          node2see.emplace(thisOriginId, std::get<3>(pp));
          subtreeUF[thisOriginId] = new ExtendedUnionFind(0);
          ++nbPairMerged;
          if(DEBUG) {
            std::cout << "willSee " << printNode(thisOriginId) << std::endl;
          }
        } else
          break;
      }

      //---
      // In UF :
      //  Origin is the size of the segmentation
      //  Data is negative : -idNodeRoot-1
      //  Data is positive : Receptacle Arc id
      //--
      // Use the queue to mark Arc that will be merged and UF to identify
      // subtree. When a node have only one way out : enqueue it to continue
      // travresall
      while(!node2see.empty()) {
        idNode curNodeId;
        bool needToGoUp; // identify up/down traversall

        std::tie(curNodeId, needToGoUp) = node2see.front();
        // should have only one arc valid to take
        node2see.pop();

        if(DEBUG) {
          std::cout << "process : " << printNode(curNodeId) << std::endl;
        }

        // Here we take the only available arc :
        idSuperArc mergingArcId;
        idNode parentNodeId;
        // continue traversall
        if(needToGoUp) {
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

        // if we have processed all but one arc of this node, we nee to continue
        // traversall
        // throug it
        if(valenceOffset[parentNodeId].first
             + valenceOffset[parentNodeId].second + 1
           == getNode(parentNodeId)->getValence()) {
          // only one way out, is it up ?
          node2see.emplace(
            parentNodeId, valenceOffset[parentNodeId].second + 1 ==

                            getNode(parentNodeId)->getUpValence());
          if(DEBUG) {
            std::cout << " add to see " << printNode(parentNodeId) << std::endl;
          }
        }
      } // end while node2see

      // for each node valenceOffset is the number of arc attached to this node
      // that will merge

      // Debug print
      if(DEBUG) {
        std::cout << "node subtrees before creating receptarc " << std::endl;
        for(idNode nid = 0; nid < nbNode; nid++) {
          if(subtreeUF[nid]) {
            std::cout << "node " << getNode(nid)->getVertexId()
                      << " is in subtree rooted :";
            const idNode &root = -subtreeUF[nid]->find()->getData() - 1;
            std::cout << getNode(root)->getVertexId();
            const SimplexId &segmSize = subtreeUF[nid]->find()->getOrigin();
            std::cout << " with segmentation of " << segmSize << std::endl;
          }
        }
      }

      //}
      //----------
      // Create the recept'arcs
      //-----------
      //{

      // Add the origin of all pairs that need to merge
      for(const auto &pp : sortedPairs) {
        if(std::get<2>(pp) < params_->simplifyThreshold) {
          const SimplexId &thisOriginVert = std::get<0>(pp);
          const SimplexId &thisEndVert = std::get<1>(pp);
          const idNode &thisOriginId = getCorrespondingNodeId(thisOriginVert);
          // const idNode &  thisEndId      = getCorrespondingNode(thisEndVert);

          if(scalars_->sosOffsets[thisOriginVert] <= posSeed0
             || scalars_->sosOffsets[thisOriginVert] >= posSeed1
             || scalars_->sosOffsets[thisEndVert] <= posSeed0
             || scalars_->sosOffsets[thisEndVert] >= posSeed1) {
            continue;
          }

          if(subtreeUF[thisOriginId]->find()->getData() < 0) {
            // create receptarc
            const idNode &subtreeRoot
              = -subtreeUF[thisOriginId]->find()->getData() - 1;
            // The id of the next arc to be created : NOT PARALLEL
            const idSuperArc receptArcId = treeData_.superArcs.size();
            // down , up, segmentation size
            // create the receptacle arc and merge arc not in sub-tree in it
            const std::tuple<idNode, idNode, SimplexId> &receptArc
              = createReceptArc(
                subtreeRoot, receptArcId, subtreeUF, valenceOffset);

            // make superArc and do the makeAlloc on it
            const bool overlapB =

              scalars_
                ->sosOffsets[getNode(std::get<0>(receptArc))->getVertexId()]
              < posSeed0;
            const bool overlapA =

              scalars_
                ->sosOffsets[getNode(std::get<1>(receptArc))->getVertexId()]
              >= posSeed1;
            const idSuperArc na
              = makeSuperArc(std::get<0>(receptArc), std::get<1>(receptArc),
                             overlapB, overlapA, nullptr, -1);

            if(overlapB) {
              treeData_.arcsCrossingBelow.emplace_back(na);
            }

            if(overlapA) {
              treeData_.arcsCrossingAbove.emplace_back(na);
            }

            subtreeUF[thisOriginId]->find()->setData(receptArcId);
            getSuperArc(receptArcId)->makeAllocGlobal(std::get<2>(receptArc));

            if(DEBUG) {
              std::cout << "create arc : " << printArc(receptArcId)
                        << " with segm : " << std::get<2>(receptArc)
                        << std::endl;
            }
          }
        } else
          break;
      }

      //}
      //----------
      // Merge these subtrees
      //-----------
      //{

      // The merge is done after because we don't want to forget arcs parallel
      // to the recept'arc but we cannot make the difference before the former
      // are created

      // nbArcs is before the insertion of receptarcs so they will no be crossed
      // here
      for(idSuperArc arc = 0; arc < nbArcs; arc++) {

        if(getSuperArc(arc)->isMerged()) {
          const idSuperArc &receptacleArcId
            = getSuperArc(arc)->getReplacantArcId();
          // take care of connectivity
          mergeArc(arc, receptacleArcId);
          if(DEBUG) {
            std::cout << " parallel merge in " << printArc(receptacleArcId)
                      << std::endl;
            std::cout << " arc " << printArc(arc)
                      << " size : " << getSuperArc(arc)->getVertSize()
                      << std::endl;
          }
          getSuperArc(receptacleArcId)
            ->addSegmentationGlobal(
              getSuperArc(arc)->getVertList(), getSuperArc(arc)->getVertSize());
        } else {
          const idNode &downNode = getSuperArc(arc)->getDownNodeId();
          const idNode &upNode = getSuperArc(arc)->getUpNodeId();

          if(!(subtreeUF[downNode] && subtreeUF[upNode]))
            continue;

          if(subtreeUF[downNode] && subtreeUF[upNode]
             && subtreeUF[downNode]->find() != subtreeUF[upNode]->find()) {
            if(DEBUG) {
              std::cout << "Arc between 2 degenerate with mergin "
                        << printArc(arc) << std::endl;
              std::cout << "below recept : "
                        << printArc(subtreeUF[downNode]->find()->getData());
              std::cout << std::endl;
              std::cout << "Above recept : "
                        << printArc(subtreeUF[upNode]->find()->getData())
                        << std::endl;
              std::cout << std::endl;
            }

            continue;
          }

          ExtendedUnionFind *curUF = (subtreeUF[upNode])
                                       ? subtreeUF[upNode]->find()
                                       : subtreeUF[downNode]->find();

          const idSuperArc &receptacleArcId = curUF->getData();

          if(DEBUG) {
            std::cout << "merge in " << printArc(receptacleArcId) << std::endl;
            std::cout << "  arc " << printArc(arc)
                      << " size : " << getSuperArc(arc)->getVertSize();
            std::cout << std::endl;
          }

          getSuperArc(arc)->merge(receptacleArcId);
          if(getSuperArc(arc)->getVertSize()) {
            getSuperArc(receptacleArcId)
              ->addSegmentationGlobal(getSuperArc(arc)->getVertList(),
                                      getSuperArc(arc)->getVertSize());

            getSuperArc(receptacleArcId)
              ->addSegmentationGlobal(getNode(downNode)->getVertexId());

            getSuperArc(receptacleArcId)
              ->addSegmentationGlobal(getNode(upNode)->getVertexId());
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

    template <typename scalarType, typename triangulationType>
    int MergeTree::computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
      const triangulationType &mesh) {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!getNumberOfSuperArcs()) {
        return -1;
      }
#endif

      pairs.reserve(treeData_.leaves.size());

      for(const idNode &leave : treeData_.leaves) {
        Node *curNode = getNode(leave);
        SimplexId curVert = curNode->getVertexId();
        SimplexId termVert = getNode(curNode->getTerminaison())->getVertexId();

        addPair<scalarType>(pairs, curVert, termVert, mesh);
      }

      auto pair_sort
        = [](const std::tuple<SimplexId, SimplexId, scalarType> &a,
             const std::tuple<SimplexId, SimplexId, scalarType> &b) {
            return std::get<2>(a) < std::get<2>(b);
          };

      sort(pairs.begin(), pairs.end(), pair_sort);

      return 0;
    }

    template <typename scalarType, typename triangulationType>
    int MergeTree::computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairs,
      const triangulationType &mesh) {
      // Need to be called on MergeTree, not ContourTree

#ifndef TTK_ENABLE_KAMIKAZE
      if(!getNumberOfSuperArcs()) {
        return -1;
      }

      if(treeData_.treeType == TreeType::Contour) {
        std::cout << "WARNING, computePersistencePairs is made to be called on "
                  << "Join or Split Tree" << std::endl;
      }
#endif

      pairs.reserve(treeData_.leaves.size());

      for(const idNode &leave : treeData_.leaves) {
        Node *curNode = getNode(leave);
        SimplexId curVert = curNode->getVertexId();
        SimplexId termVert = getNode(curNode->getTerminaison())->getVertexId();

        addPair<scalarType>(
          pairs, curVert, termVert, mesh, treeData_.treeType == TreeType::Join);
      }

      auto pair_sort
        = [](const std::tuple<SimplexId, SimplexId, scalarType, bool> &a,
             const std::tuple<SimplexId, SimplexId, scalarType, bool> &b) {
            return std::get<2>(a) < std::get<2>(b);
          };

      sort(pairs.begin(), pairs.end(), pair_sort);

      return 0;
    }

    template <typename scalarType, typename triangulationType>
    void MergeTree::recoverMTPairs(
      const std::vector<idNode> &sortedNodes,
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairsJT,
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairsST,
      const triangulationType &mesh) {
      const auto nbNode = getNumberOfNodes();

      std::vector<ExtendedUnionFind *> vect_JoinUF(nbNode, nullptr);
      std::vector<ExtendedUnionFind *> vect_SplitUF(nbNode, nullptr);

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
          std::map<SimplexId, SimplexId> pendingMinMax;

          for(auto it = sortedNodes.cbegin(); it != sortedNodes.cend(); ++it) {
            const idNode &n = *it;
            const SimplexId &v = getNode(n)->getVertexId();

            const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
            const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

            if(nbDown == 0) {
              // leaf
              vect_JoinUF[n] = new ExtendedUnionFind(getNode(n)->getVertexId());
              vect_JoinUF[n]->setOrigin(v);
              // std::cout << " jt origin : " << v << std::endl;
            } else {
              // first descendant
              const idSuperArc firstSaId = getNode(n)->getDownSuperArcId(0);
              const SuperArc *firstSA = getSuperArc(firstSaId);
              const idNode &firstChildNodeId = firstSA->getDownNodeId();

              ExtendedUnionFind *merge = vect_JoinUF[firstChildNodeId]->find();
              SimplexId further = merge->getOrigin();
              idSuperArc furtherI = 0;

              // Find the most persistant way
              for(idSuperArc ni = 1; ni < nbDown; ++ni) {
                const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                const SuperArc *curSA = getSuperArc(curSaId);
                const idNode &neigh = curSA->getDownNodeId();

                // fix
                if(neigh == n)
                  continue;

                ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                if(isLower(neighUF->getOrigin(), further)) {
                  further = neighUF->getOrigin();
                  furtherI = ni;
                }
              }

              if(nbDown > 1) {
                // close finish pair and make union
                for(idSuperArc ni = 0; ni < nbDown; ++ni) {
                  const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                  const SuperArc *curSA = getSuperArc(curSaId);
                  const idNode &neigh = curSA->getDownNodeId();

                  // fix
                  if(neigh == n)
                    continue;

                  ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                  if(ni != furtherI) { // keep the more persitent pair
                    addPair<scalarType>(
                      pairsJT, neighUF->getOrigin(), v, mesh, true);
                    pendingMinMax.erase(neighUF->getOrigin());

                    // std::cout << " jt make pair : " <<
                    // neighUF->getOrigin() << " - " << v <<
                    // std::endl;
                  }

                  ExtendedUnionFind::makeUnion(merge, neighUF)
                    ->setOrigin(further);
                }
              }

              merge->find()->setOrigin(further);
              vect_JoinUF[n] = merge->find();

              if(!nbUp) {
                // potential close of the component
                // std::cout << "pending for " << further << " is " << v <<
                // std::endl;
                pendingMinMax[further] = v;
              }
            }
          } // end for each node

          // Add the pending biggest pair of each component
          for(const auto &pair_vert : pendingMinMax) {
            if(isCorrespondingNode(pair_vert.second)) {
              // std::cout << " add : " << pair_vert.first << " - " <<
              // pair_vert.second << std::endl;
              addPair<scalarType>(
                pairsJT, pair_vert.first, pair_vert.second, mesh, true);
            }
          }
        } // end para section

#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
        {
          pairsST.reserve(treeData_.leaves.size());
          // For the biggest pair of the component
          std::map<SimplexId, SimplexId> pendingMinMax;

          for(auto it = sortedNodes.crbegin(); it != sortedNodes.crend();
              ++it) {
            const idNode &n = *it;
            const SimplexId &v = getNode(n)->getVertexId();

            const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
            const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

            if(nbUp == 0) {
              // leaf
              vect_SplitUF[n]
                = new ExtendedUnionFind(getNode(n)->getVertexId());
              vect_SplitUF[n]->setOrigin(v);
              // std::cout << " st origin : " << v << std::endl;
            } else {
              // first descendant
              const idSuperArc firstSaId = getNode(n)->getUpSuperArcId(0);
              const SuperArc *firstSA = getSuperArc(firstSaId);
              const idNode &firstChildNodeId = firstSA->getUpNodeId();

              ExtendedUnionFind *merge = vect_SplitUF[firstChildNodeId]->find();
              SimplexId further = merge->getOrigin();
              idSuperArc furtherI = 0;

              for(idSuperArc ni = 1; ni < nbUp; ++ni) {
                // find the more persistant way
                const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                const SuperArc *curSA = getSuperArc(curSaId);
                // Ignore hidden / simplified arc
                if(!curSA->isVisible())
                  continue;
                const idNode &neigh = curSA->getUpNodeId();

                // fix
                if(neigh == n)
                  continue;

                // std::cout << "visit neighbor : " << ni << " which is " <<
                // getNode(neigh)->getVertexId() << std::endl;

                ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                if(isHigher(neighUF->getOrigin(), further)) {
                  further = neighUF->getOrigin();
                  furtherI = ni;
                }
              }

              if(nbUp > 1) {
                // close finsh pair and make union
                for(idSuperArc ni = 0; ni < nbUp; ++ni) {
                  const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                  const SuperArc *curSA = getSuperArc(curSaId);
                  const idNode &neigh = curSA->getUpNodeId();

                  // fix
                  if(neigh == n)
                    continue;

                  ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                  if(ni != furtherI) {
                    addPair<scalarType>(
                      pairsST, neighUF->getOrigin(), v, mesh, false);

                    pendingMinMax.erase(neighUF->getOrigin());

                    // std::cout << " st make pair : " <<
                    // neighUF->getOrigin() << " - " << v
                    //<< " for neighbor " <<
                    // getNode(neigh)->getVertexId() << std::endl;
                  }

                  ExtendedUnionFind::makeUnion(merge, neighUF)
                    ->setOrigin(further);
                  // Re-visit after merge lead to add the most persistant
                  // pair....
                }
              }
              merge->find()->setOrigin(further);
              vect_SplitUF[n] = merge->find();

              if(!nbDown) {
                pendingMinMax[further] = v;
              }

            } // end nbUp == 0 else
          } // end for each node

          // Add the pending biggest pair of each component
          for(const auto &pair_vert : pendingMinMax) {
            addPair<scalarType>(
              pairsST, pair_vert.first, pair_vert.second, mesh, false);
          }
        } // end para section
      } // end para
    }

    // }
    // Tools
    // {

    template <typename scalarType, typename triangulationType>
    void MergeTree::addPair(
      std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairs,
      const SimplexId &orig,
      const SimplexId &term,
      const triangulationType &mesh,
      const bool goUp) {
      if(params_->simplifyMethod == SimplifMethod::Persist) {
        pairs.emplace_back(
          orig, term,
          std::abs(static_cast<double>(getValue<scalarType>(orig)
                                       - getValue<scalarType>(term))),
          goUp);
      } else if(params_->simplifyMethod == SimplifMethod::Span) {
        float coordOrig[3], coordTerm[3], span;
        mesh->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
        mesh->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
        span = Geometry::distance(coordOrig, coordTerm);
        pairs.emplace_back(orig, term, span, goUp);
      }
    }

    template <typename scalarType, typename triangulationType>
    void MergeTree::addPair(
      std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
      const SimplexId &orig,
      const SimplexId &term,
      const triangulationType &mesh) {
      if(params_->simplifyMethod == SimplifMethod::Persist) {
        pairs.emplace_back(
          orig, term,
          std::abs(static_cast<double>(getValue<scalarType>(orig)
                                       - getValue<scalarType>(term))));
      } else if(params_->simplifyMethod == SimplifMethod::Span) {
        float coordOrig[3], coordTerm[3], span;
        mesh->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
        mesh->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
        span = Geometry::distance(coordOrig, coordTerm);
        pairs.emplace_back(orig, term, span);
      }
    }

    template <typename triangulationType>
    int MergeTree::build(std::vector<ExtendedUnionFind *> &vect_baseUF,
                         const std::vector<SimplexId> &overlapBefore,
                         const std::vector<SimplexId> &overlapAfter,
                         SimplexId start,
                         SimplexId end,
                         const SimplexId &posSeed0,
                         const SimplexId &posSeed1,
                         const triangulationType &mesh) {
      // idea, work on the neighbohood instead of working on the node itsef.
      // Need lower / higher star construction.
      //
      // at this time, ST have no root except in Adjacency list. These root are
      // isolated vertices.
      //
      // clear and reset tree data (this step should take almost no time)
      flush();

      DebugTimer timerBegin;

      // -----------------
      // Find boundaries
      // -----------------
      // {

      SimplexId sortedNode;
      const SimplexId step = (treeData_.treeType == TreeType::Join) ? 1 : -1;

      // main
      const SimplexId mainStart = start;
      const SimplexId mainEnd = end;

      const bool isJT = treeData_.treeType == TreeType::Join;
      // else Split Tree, can't be called on ContourTree

      // overlap before
      const SimplexId beforeStart = (isJT) ? 0 : overlapBefore.size() - 1;
      const SimplexId beforeEnd = (isJT) ? overlapBefore.size() : -1;

      // overlap after
      const SimplexId afterStart = (isJT) ? 0 : overlapAfter.size() - 1;
      const SimplexId afterEnd = (isJT) ? overlapAfter.size() : -1;

      // print debug
      if(params_->debugLevel >= 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          std::stringstream msg;
          msg << "partition : " << static_cast<unsigned>(treeData_.partition);
          msg << ", isJT : " << isJT;
          msg << ",  size : ";
          msg << "before  : " << std::abs(beforeEnd - beforeStart);
          msg << ", main : " << std::abs(mainEnd - mainStart);
          msg << ", after : " << std::abs(afterEnd - afterStart);
          msg << ", Total : ";
          msg << std::abs(beforeEnd - beforeStart)
                   + std::abs(mainEnd - mainStart)
                   + std::abs(afterEnd - afterStart);
          this->printMsg(msg.str());
        }
      }

      // }
      // --------------
      // Overlap Before
      // --------------
      // {

      // for each vertex of our triangulation
      for(sortedNode = beforeStart; sortedNode != beforeEnd;
          sortedNode += step) {
        const SimplexId currentVertex = overlapBefore[sortedNode];
        const bool overlapB = isJT;
        const bool overlapA = !isJT;
        processVertex(
          currentVertex, vect_baseUF, overlapB, overlapA, mesh, timerBegin);
      } // foreach node

      // }
      // ---------------
      // Partition
      // ---------------
      // {

      // for each vertex of our triangulation
      for(sortedNode = mainStart; sortedNode != mainEnd; sortedNode += step) {
        const SimplexId currentVertex = scalars_->sortedVertices[sortedNode];
        processVertex(
          currentVertex, vect_baseUF, false, false, mesh, timerBegin);
      } // foreach node

      // }
      // ---------------
      // Overlap After
      // ---------------
      // {

      // for each vertex of our triangulation
      for(sortedNode = afterStart; sortedNode != afterEnd; sortedNode += step) {
        const SimplexId currentVertex = overlapAfter[sortedNode];
        const bool overlapB = !isJT;
        const bool overlapA = isJT;
        processVertex(
          currentVertex, vect_baseUF, overlapB, overlapA, mesh, timerBegin);
      } // foreach node

      // }
      // ---------------
      // Close root arcs
      // ---------------
      // {

      // Closing step for openedSuperArc.
      // Not aesthetic but efficient
      // More efficient that using nullity of the neighborhood
      // to detect extrema of the opponent tree
      idNode rootNode;
      SimplexId corrVertex, origin;
      idSuperArc tmp_sa;

      // It can't be more connected component that leaves so test for each
      // leaves (even virtual extrema)
      for(const auto &l : treeData_.leaves) {
        corrVertex = getNode(l)->getVertexId();

        // case of an isolated point
        if(!mesh->getVertexNeighborNumber(corrVertex)) {
          tmp_sa = getNode(l)->getUpSuperArcId(0);
        } else {
          tmp_sa = (idSuperArc)((vect_baseUF[corrVertex])->find()->getData());
          origin = (idSuperArc)((vect_baseUF[corrVertex])->find()->getOrigin());
        }

        if(treeData_.superArcs[tmp_sa].getUpNodeId() == nullNodes) {
          rootNode
            = makeNode(treeData_.superArcs[tmp_sa].getLastVisited(), origin);
          Node *originNode = getNode(
            getNode(getSuperArc(tmp_sa)->getDownNodeId())->getOrigin());
          originNode->setTerminaison(rootNode);

          const bool overlapB
            = scalars_->sosOffsets[getNode(rootNode)->getVertexId()]
              <= posSeed0;
          const bool overlapA
            = scalars_->sosOffsets[getNode(rootNode)->getVertexId()]
              >= posSeed1;

          closeSuperArc(tmp_sa, rootNode, overlapB, overlapA);

          // in the case we have 1 vertex domain,
          // hide the close SuperArc wich is point1 <>> point1
          if(getSuperArc(tmp_sa)->getDownNodeId()
             == getSuperArc(tmp_sa)->getUpNodeId()) {
            hideArc(tmp_sa);
          }

          treeData_.roots.emplace_back(rootNode);
        }
      }

      // }
      // -----------
      // Timer print
      // ------------
      // {

      this->printMsg("Tree " + std::to_string(treeData_.partition)
                       + " computed (" + std::to_string(getNumberOfSuperArcs())
                       + " arcs)",
                     1.0, timerBegin.getElapsedTime(), this->threadNumber_);

      // }

      return 0;
    }

    template <typename triangulationType>
    void MergeTree::processVertex(const SimplexId &currentVertex,
                                  std::vector<ExtendedUnionFind *> &vect_baseUF,
                                  const bool overlapB,
                                  const bool overlapA,
                                  const triangulationType &mesh,
                                  DebugTimer &begin) {
      std::vector<ExtendedUnionFind *> vect_neighUF;
      ExtendedUnionFind *seed = nullptr, *tmpseed;

      SimplexId neighSize;
      const SimplexId neighborNumber
        = mesh->getVertexNeighborNumber(currentVertex);
      const bool isJT = treeData_.treeType == TreeType::Join;

      idSuperArc currentArc;
      idNode closingNode, currentNode;
      SimplexId neighbor;

      // Check UF in neighborhood
      for(SimplexId n = 0; n < neighborNumber; ++n) {
        mesh->getVertexNeighbor(currentVertex, n, neighbor);
        // if the vertex is out: consider it null
        tmpseed = vect_baseUF[neighbor];
        // unvisited vertex, we continue.
        if(tmpseed == nullptr) {
          continue;
        }

        tmpseed = tmpseed->find();

        // get all different UF in neighborhood
        if(find(vect_neighUF.cbegin(), vect_neighUF.cend(), tmpseed)
           == vect_neighUF.end()) {
          vect_neighUF.emplace_back(tmpseed);
          seed = tmpseed;
        }
      }

      (neighSize = vect_neighUF.size());

      // SimplexId test = 1;
      // if (currentVertex == test)
      // cout << test << " : " << vect_neighUF.size() << " " <<
      // vect_interfaceUF.size() << endl; Make output
      if(!neighSize) {
        // we are on a real extrema we have to create a new UNION FIND and a
        // branch a real extrema can't be a virtual extrema

        seed = new ExtendedUnionFind(currentVertex);
        // When creating an extrema we create a pair ending on this node.
        currentNode = makeNode(currentVertex);
        getNode(currentNode)->setOrigin(currentNode);
        currentArc = openSuperArc(currentNode, overlapB, overlapA);
        // if(overlap && partition_ == 1) cout << currentVertex << endl;
        treeData_.leaves.emplace_back(currentNode);

        if(params_->debugLevel >= static_cast<int>(debug::Priority::DETAIL)) {
          if(isJT) {
            this->printMsg("Min node id: " + std::to_string(currentVertex), 1.0,
                           begin.getElapsedTime(), this->threadNumber_);
          } else {
            this->printMsg("Max node id: " + std::to_string(currentVertex), 1.0,
                           begin.getElapsedTime(), this->threadNumber_);
          }
        }
      } else if(neighSize > 1) {
        // Is a saddle if have more than one UF in neighborhood
        // Or is linked to an interface (virtual extrema -> leaves or
        // interface cc representant : regular node )

        // Merge operation => close all arriving SuperArc (ignoring UF from
        // interface) And open a new one
        closingNode = makeNode(currentVertex);
        currentArc = openSuperArc(closingNode, overlapB, overlapA);

        SimplexId farOrigin = vect_neighUF[0]->find()->getOrigin();

        // close each SuperArc finishing here
        for(auto *neigh : vect_neighUF) {
          closeSuperArc((idSuperArc)neigh->find()->getData(), closingNode,
                        overlapB, overlapA);
          // persistance pair closing here.
          // For the one who will continue, it will be overide later
          vertex2Node(neigh->find()->getOrigin())->setTerminaison(closingNode);

          // cout <<
          // getNode(getCorrespondingNode(neigh->find()->getOrigin()))->getVertexId()
          //<< " terminate on " << getNode(closingNode)->getVertexId() << endl;

          if((isJT && isLower(neigh->find()->getOrigin(), farOrigin))
             || (!isJT && isHigher(neigh->find()->getOrigin(), farOrigin))) {
            // here we keep the continuing the most persitant pair.
            // It means a pair end when a parent have another origin thant the
            // current leaf (or is the root) It might be not intuitive but it is
            // more convenient for degenerate cases
            farOrigin = neigh->find()->getOrigin();
            // cout << "find origin  " << farOrigin << " for " << currentVertex
            // << " " << isJT
            //<< endl;
          }
        }

        // Union correspond to the merge
        seed = ExtendedUnionFind::makeUnion(vect_neighUF);
        if(seed == nullptr) {
          return;
        }
        seed->setOrigin(farOrigin);
        getNode(closingNode)->setOrigin(getCorrespondingNodeId(farOrigin));

        // cout << "  " << getNode(closingNode)->getVertexId() << " have origin
        // at "
        //<< getNode(getCorrespondingNode(farOrigin))->getVertexId() << endl;

        this->printMsg("Saddle node id: " + std::to_string(currentVertex),
                       debug::Priority::DETAIL);

      } else {
#ifndef TTK_ENABLE_KAMIKAZE
        if(seed == nullptr) {
          return;
        }
#endif // TTK_ENABLE_KAMIKAZE
       // regular node
        currentArc = (idSuperArc)seed->find()->getData();
        updateCorrespondingArc(currentVertex, currentArc);
      }
      // common
      seed->setData((ufDataType)currentArc);
      getSuperArc(currentArc)->setLastVisited(currentVertex);
      vect_baseUF[currentVertex] = seed;
    }

    template <typename triangulationType>
    bool MergeTree::verifyTree(const triangulationType &mesh) {
      bool res = true;

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
      {
        const idSuperArc &nbArcs = getNumberOfSuperArcs();
        const idSuperArc &nbNodes = getNumberOfNodes();

        std::cout << "Verify Tree : " << std::endl;
        std::cout << "nbNode initial : " << nbNodes << std::endl;
        std::cout << "nbArcs initial : " << nbArcs << std::endl;

        idSuperArc nbArcsVisibles = 0;
        idSuperArc nbNodesVisibles = 0;

        // for each visible arc, verify he is in the node

        for(idSuperArc aid = 0; aid < nbArcs; aid++) {
          const SuperArc &arc = treeData_.superArcs[aid];
          if(arc.isVisible()) {
            ++nbArcsVisibles;
            const idNode &up = arc.getUpNodeId();
            const idNode &down = arc.getDownNodeId();
            if(up == nullNodes || down == nullNodes) {
              res = false;
              std::cout << "[Verif]: arc id : " << aid
                        << "have a null boundary :";
              std::cout << " down :" << down << " up:" << up << std::endl;
            } else {
              bool isIn = false;

              // Arc is present in its upNode
              const Node &nup = treeData_.nodes[up];
              const auto &upNbDown = nup.getNumberOfDownSuperArcs();
              for(idSuperArc d = 0; d < upNbDown; d++) {
                if(nup.getDownSuperArcId(d) == aid) {
                  isIn = true;
                  break;
                }
              }
              if(!isIn) {
                res = false;
                std::cout << "[Verif]: arc " << printArc(aid)
                          << " is not known by its up node :";
                std::cout << treeData_.nodes[up].getVertexId() << std::endl;
              }

              isIn = false;

              // Arc is present in its upNode
              const Node &ndown = treeData_.nodes[arc.getDownNodeId()];
              const auto &upNbUp = ndown.getNumberOfUpSuperArcs();
              for(idSuperArc u = 0; u < upNbUp; u++) {
                if(ndown.getUpSuperArcId(u) == aid) {
                  isIn = true;
                  break;
                }
              }
              if(!isIn) {
                res = false;
                std::cout << "[Verif]: arc " << printArc(aid)
                          << " is not known by its down node :";
                std::cout << treeData_.nodes[down].getVertexId() << std::endl;
              }
            }
          }
        }

        // for each node, verify she is in the arc

        for(idNode nid = 0; nid < nbNodes; nid++) {
          const Node &node = treeData_.nodes[nid];
          if(!node.isHidden()) {
            ++nbNodesVisibles;

            // Verify up arcs
            const auto &nbup = node.getNumberOfUpSuperArcs();
            for(idSuperArc ua = 0; ua < nbup; ua++) {
              const SuperArc &arc
                = treeData_.superArcs[node.getUpSuperArcId(ua)];
              const idNode arcDownNode = arc.getDownNodeId();
              if(arcDownNode != nid || !arc.isVisible()) {
                res = false;
                const idNode upnode = arc.getUpNodeId();
                const idNode downnode = arc.getDownNodeId();
                if(upnode == nullNodes || downnode == nullNodes) {
                  std::cout << "[Verif]: arc id : " << node.getUpSuperArcId(ua);
                  std::cout << "have a null boundary :";
                  std::cout << " down :" << downnode << " up:" << upnode
                            << std::endl;
                } else {
                  std::cout << "[Verif] Node " << node.getVertexId()
                            << " id : " << nid;
                  std::cout << " Problem with up arc : "
                            << printArc(node.getUpSuperArcId(ua)) << std::endl;
                }
              }
            }

            // Verify down arcs
            const auto &nbdown = node.getNumberOfDownSuperArcs();
            for(idSuperArc da = 0; da < nbdown; da++) {
              const SuperArc &arc
                = treeData_.superArcs[node.getDownSuperArcId(da)];
              const idNode arcUpNode = arc.getUpNodeId();
              if(arcUpNode != nid || !arc.isVisible()) {
                res = false;
                const idNode upnode = arc.getUpNodeId();
                const idNode downnode = arc.getDownNodeId();
                if(upnode == nullNodes || downnode == nullNodes) {
                  std::cout << "[Verif]: arc id : "
                            << node.getDownSuperArcId(da);
                  std::cout << "have a null boundary :";
                  std::cout << " down :" << downnode << " up:" << upnode
                            << std::endl;
                } else {
                  std::cout << "[Verif] Node " << node.getVertexId()
                            << " id : " << nid;
                  std::cout << " Problem with down arc : "
                            << printArc(node.getUpSuperArcId(da)) << std::endl;
                }
              }
            }
          }
        }

        // verify segmentation information
        const auto &nbVert = mesh->getNumberOfVertices();
        std::vector<bool> segmSeen(nbVert, false);

        for(idSuperArc aid = 0; aid < nbArcs; aid++) {
          SuperArc &arc = treeData_.superArcs[aid];

          if(!arc.isVisible())
            continue;

          const SimplexId segmSize = arc.getVertSize();
          const std::pair<SimplexId, bool> *segmVect = arc.getVertList();

          if(segmSize && !segmVect) {
            res = false;
            std::cout << "[Verif] Inconsistant segmentation for arc : ";
            std::cout << printArc(aid);
            std::cout << " have size of " << segmSize;
            std::cout << " and a null list" << std::endl;
          }

          if(segmVect != nullptr) {
            for(SimplexId v = 0; v < segmSize; v++) {
              if(!segmVect[v].second) {
                segmSeen.at(segmVect[v].first) = true;
              }
            }
          }
        }

        for(const Node &node : treeData_.nodes) {
          if(node.isHidden())
            continue;

          segmSeen.at(node.getVertexId()) = true;
        }

        std::cout << "Segm missing : ";
        for(SimplexId v = 0; v < nbVert; v++) {
          if(!segmSeen[v]) {
            res = false;
            std::cout << v << ", ";
          }
        }
        std::cout << std::endl;

        std::cout << "Nb visible Node : " << nbNodesVisibles << std::endl;
        std::cout << "Nb visible Arcs : " << nbArcsVisibles << std::endl;
      }
      return res;
    }

    // }
  } // namespace cf
} // namespace ttk
