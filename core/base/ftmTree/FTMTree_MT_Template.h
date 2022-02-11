/// \ingroup base
/// \class ttk::FTMTree_MT
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
/// \sa ttkContourForests.cpp %for a usage example.

#ifndef FTMTREE_MT_TPL_H
#define FTMTREE_MT_TPL_H

#include <functional>

#include "FTMTree_MT.h"

#ifdef TTK_ENABLE_OMP_PRIORITY
#define OPTIONAL_PRIORITY(value) priority(value)
#else
#define OPTIONAL_PRIORITY(value)
#endif

// ----
// Init
// ----

namespace ttk {
  namespace ftm {

    template <class triangulationType>
    void FTMTree_MT::build(const triangulationType *mesh, const bool ct) {
      std::string treeString;
      // Comparator init (template)
      initComp();
      switch(mt_data_.treeType) {
        case TreeType::Join:
          treeString = "JT";
          break;
        case TreeType::Split:
          treeString = "ST";
          break;
        default:
          treeString = "CT";
          break;
      }

      // Build Merge treeString using tasks
      Timer precomputeTime;
      int alreadyDone = leafSearch(mesh);
      printTime(precomputeTime, "leafSearch " + treeString, 3 + alreadyDone);

      Timer buildTime;
      leafGrowth(mesh);
#ifdef TTK_ENABLE_FTM_TREE_PROCESS_SPEED
      // count process
      for(SimplexId i = 0; i < scalars_->size; i++) {
        if((*mt_data_.vert2tree)[i] != nullCorresp)
          ++nbProcessed;
      }
#endif
      printTime(buildTime, "leafGrowth " + treeString, 3);

      Timer bbTime;
      trunk(mesh, ct);
      printTime(bbTime, "trunk " + treeString, 3);

      // Segmentation
      if(ct && params_->segm) {
        Timer segmTime;
        buildSegmentation();
        printTime(segmTime, "segment " + treeString, 3);
      }
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    int FTMTree_MT::leafSearch(const triangulationType *mesh) {
      int ret = 0;
      // if not already computed by CT
      if(getNumberOfNodes() == 0) {
        const auto nbScalars = scalars_->size;
        const auto chunkSize = getChunkSize();
        const auto chunkNb = getChunkCount();

        // Extrema extract and launch tasks
        for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) OPTIONAL_PRIORITY(isPrior())
#endif
          {
            const SimplexId lowerBound = chunkId * chunkSize;
            const SimplexId upperBound
              = std::min(nbScalars, (chunkId + 1) * chunkSize);
            for(SimplexId v = lowerBound; v < upperBound; ++v) {
              const auto &neighNumb = mesh->getVertexNeighborNumber(v);
              valence val = 0;

              for(valence n = 0; n < neighNumb; ++n) {
                SimplexId neigh;
                mesh->getVertexNeighbor(v, n, neigh);
                comp_.vertLower(neigh, v) && ++val;
              }

              (*mt_data_.valences)[v] = val;

              if(!val) {
                makeNode(v);
              }
            }
          }
        }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
      } else {
        ret = 1;
      }

      // fill leaves
      const auto &nbLeaves = mt_data_.nodes->size();
      mt_data_.leaves->resize(nbLeaves);
      std::iota(mt_data_.leaves->begin(), mt_data_.leaves->end(), 0);

      if(debugLevel_ >= 4) {
        this->printMsg("found " + std::to_string(nbLeaves) + " leaves");
      }

      // Reserve Arcs
      mt_data_.superArcs->reserve(nbLeaves * 2 + 1);
#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      createVector<ActiveTask>(mt_data_.activeTasksStats);
      mt_data_.activeTasksStats->resize(nbLeaves * 2 + 1);
#endif

      return ret;
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    void FTMTree_MT::leafGrowth(const triangulationType *mesh) {
      _launchGlobalTime.reStart();

      const auto &nbLeaves = mt_data_.leaves->size();

      // memory allocation here
      initVectStates(nbLeaves + 2);

      // elevation: backbone only
      if(nbLeaves == 1) {
        const SimplexId v = (*mt_data_.nodes)[0].getVertexId();
        (*mt_data_.openedNodes)[v] = 1;
        (*mt_data_.ufs)[v] = new AtomicUF(v);
        return;
      }

      mt_data_.activeTasks = nbLeaves;

      auto comp = [this](const idNode a, const idNode b) {
#ifdef HIGHER
        return this->comp_.vertHigher(
          this->getNode(a)->getVertexId(), this->getNode(b)->getVertexId());
#else
        return this->comp_.vertLower(
          this->getNode(a)->getVertexId(), this->getNode(b)->getVertexId());
#endif
      };
      sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), comp);

      for(idNode n = 0; n < nbLeaves; ++n) {
        const idNode l = (*mt_data_.leaves)[n];
        SimplexId v = getNode(l)->getVertexId();
        // for each node: get vert, create uf and lauch
        (*mt_data_.ufs)[v] = new AtomicUF(v);

#ifdef TTK_ENABLE_OPENMP
#pragma omp task UNTIED() OPTIONAL_PRIORITY(isPrior())
#endif
        arcGrowth(mesh, v, n);
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    void FTMTree_MT::arcGrowth(const triangulationType *mesh,
                               const SimplexId startVert,
                               const SimplexId orig) {
      // current task id / propag

      // local order (ignore non regular verts)
      SimplexId localOrder = -1;
      UF startUF = (*mt_data_.ufs)[startVert]->find();
      // get or recover states
      CurrentState *currentState;
      if(startUF->getNbStates()) {
        currentState = startUF->getFirstState();
      } else {
        const std::size_t currentStateId = mt_data_.states->getNext();
        currentState = &(*mt_data_.states)[currentStateId];
        currentState->setStartVert(startVert);
        startUF->addState(currentState);
      }

      currentState->addNewVertex(startVert);

      // avoid duplicate processing of startVert
      bool seenFirst = false;

      // ARC OPENING
      idNode startNode = getCorrespondingNodeId(startVert);
      idSuperArc currentArc = openSuperArc(startNode);
      startUF->addArcToClose(currentArc);
#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      (*mt_data_.activeTasksStats)[currentArc].begin
        = _launchGlobalTime.getElapsedTime();
      (*mt_data_.activeTasksStats)[currentArc].origin = orig;
#endif

      // TASK PROPAGATION
      while(!currentState->empty()) {
        // Next vertex

        SimplexId currentVert = currentState->getNextMinVertex();

        // ignore duplicate
        if(!isCorrespondingNull(currentVert)
           && !isCorrespondingNode(currentVert)) {
          continue;
        } else {
          // first node can be duplicate, avoid duplicate process
          if(currentVert == startVert) {
            if(!seenFirst) {
              seenFirst = true;
            } else {
              continue;
            }
          }
        }

        // local order to avoid sort
        (*mt_data_.visitOrder)[currentVert] = localOrder++;

        // Saddle & Last detection + propagation
        bool isSaddle, isLast;
        std::tie(isSaddle, isLast) = propage(mesh, *currentState, startUF);

        // regular propagation
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
        (*mt_data_.ufs)[currentVert] = startUF;

        // Saddle case
        if(isSaddle) {

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
          (*mt_data_.activeTasksStats)[currentArc].end
            = _launchGlobalTime.getElapsedTime();
#endif
          // need a node on this vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
          (*mt_data_.openedNodes)[currentVert] = 1;

          // If last close all and merge
          if(isLast) {
            // finish works here
            closeAndMergeOnSaddle(mesh, currentVert);

            // last task detection
            idNode remainingTasks;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
            remainingTasks = mt_data_.activeTasks;
            if(remainingTasks == 1) {
              // only backbone remaining
              return;
            }

            // made a node on this vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
            (*mt_data_.openedNodes)[currentVert] = 0;

            // recursively continue
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
            arcGrowth(mesh, currentVert, orig);
          } else {
            // Active tasks / threads
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
            mt_data_.activeTasks--;
          }

          // stop at saddle
          return;
        }

        if(currentVert != startVert) {
          updateCorrespondingArc(currentVert, currentArc);
        }
        getSuperArc(currentArc)->setLastVisited(currentVert);

      } // end wile propagation

      // close root
      const SimplexId closeVert = getSuperArc(currentArc)->getLastVisited();
      bool existCloseNode = isCorrespondingNode(closeVert);
      idNode closeNode = (existCloseNode) ? getCorrespondingNodeId(closeVert)
                                          : makeNode(closeVert);
      closeSuperArc(currentArc, closeNode);
      getSuperArc(currentArc)->decrNbSeen();
      idNode rootPos = mt_data_.roots->getNext();
      (*mt_data_.roots)[rootPos] = closeNode;

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      (*mt_data_.activeTasksStats)[currentArc].end
        = _launchGlobalTime.getElapsedTime();
#endif
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    std::tuple<bool, bool> FTMTree_MT::propage(const triangulationType *mesh,
                                               CurrentState &currentState,
                                               UF curUF) {
      bool becameSaddle = false, isLast = false;
      const auto nbNeigh = mesh->getVertexNeighborNumber(currentState.vertex);
      valence decr = 0;

      // once for all
      auto *curUFF = curUF->find();

      // propagation / is saddle
      for(valence n = 0; n < nbNeigh; ++n) {
        SimplexId neigh;
        mesh->getVertexNeighbor(currentState.vertex, n, neigh);

        if(comp_.vertLower(neigh, currentState.vertex)) {
          UF neighUF = (*mt_data_.ufs)[neigh];

          // is saddle
          if(!neighUF || neighUF->find() != curUFF) {
            becameSaddle = true;
          } else if(neighUF) {
            ++decr;
          }

        } else {
          if(!(*mt_data_.propagation)[neigh]
             || (*mt_data_.propagation)[neigh]->find() != curUFF) {
            currentState.addNewVertex(neigh);
            (*mt_data_.propagation)[neigh] = curUFF;
          }
        }
      }

      // is last
      valence oldVal;
      valence &tmp = (*mt_data_.valences)[currentState.vertex];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      {
        oldVal = tmp;
        tmp -= decr;
      }
      if(oldVal == decr) {
        isLast = true;
      }

      return std::make_tuple(becameSaddle, isLast);
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    SimplexId FTMTree_MT::trunk(const triangulationType *mesh, const bool ct) {
      Timer bbTimer;

      std::vector<SimplexId> trunkVerts;
      const auto &nbScalars = scalars_->size;

      // trunkVerts
      trunkVerts.reserve(std::max(SimplexId{10}, nbScalars / 500));
      for(SimplexId v = 0; v < nbScalars; ++v) {
        if((*mt_data_.openedNodes)[v]) {
          trunkVerts.emplace_back(v);
        }
      }
      sort(trunkVerts.begin(), trunkVerts.end(), comp_.vertLower);
      for(const SimplexId v : trunkVerts) {
        closeOnBackBone(mesh, v);
      }

      // Arcs
      const auto &nbNodes = trunkVerts.size();
      for(idNode n = 1; n < nbNodes; ++n) {
        idSuperArc na = makeSuperArc(getCorrespondingNodeId(trunkVerts[n - 1]),
                                     getCorrespondingNodeId(trunkVerts[n]));
        getSuperArc(na)->setLastVisited(trunkVerts[n]);
      }

      if(!nbNodes) {
        return 0;
      }
      const idSuperArc lastArc
        = openSuperArc(getCorrespondingNodeId(trunkVerts[nbNodes - 1]));

      // Root (close last arc)
      // if several CC still the backbone is only in one.
      // But the root may not be the max node of the whole dataset: TODO
      const idNode rootNode
        = makeNode(scalars_->sortedVertices[(isJT()) ? scalars_->size - 1 : 0]);
      closeSuperArc(lastArc, rootNode);
      getSuperArc(lastArc)->setLastVisited(getNode(rootNode)->getVertexId());

      printTime(bbTimer, "trunk seq.", 4);
      bbTimer.reStart();

      // Segmentation
      SimplexId begin, stop, processed;
      std::tie(begin, stop) = getBoundsFromVerts(trunkVerts);
      if(ct) {
        processed = trunkCTSegmentation(trunkVerts, begin, stop);
      } else {
        processed = trunkSegmentation(trunkVerts, begin, stop);
      }
      printTime(bbTimer, "trunk para.", 4);

      return processed;
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    void FTMTree_MT::closeAndMergeOnSaddle(const triangulationType *mesh,
                                           SimplexId saddleVert) {
      idNode closeNode = makeNode(saddleVert);

      // Union of the UF coming here (merge propagation and closing arcs)
      const auto &nbNeigh = mesh->getVertexNeighborNumber(saddleVert);
      for(valence n = 0; n < nbNeigh; ++n) {
        SimplexId neigh;
        mesh->getVertexNeighbor(saddleVert, n, neigh);

        if(comp_.vertLower(neigh, saddleVert)) {
          if((*mt_data_.ufs)[neigh]->find()
             != (*mt_data_.ufs)[saddleVert]->find()) {
            (*mt_data_.ufs)[saddleVert] = AtomicUF::makeUnion(
              (*mt_data_.ufs)[saddleVert], (*mt_data_.ufs)[neigh]);
          }
        }
      }

      // close arcs on this node
      closeArcsUF(closeNode, (*mt_data_.ufs)[saddleVert]);

      (*mt_data_.ufs)[saddleVert]->find()->mergeStates();
      (*mt_data_.ufs)[saddleVert]->find()->setExtrema(saddleVert);
    }

    // ------------------------------------------------------------------------

    template <class triangulationType>
    void FTMTree_MT::closeOnBackBone(const triangulationType *mesh,
                                     SimplexId saddleVert) {
      idNode closeNode = makeNode(saddleVert);

      // Union of the UF coming here (merge propagation and closing arcs)
      const auto &nbNeigh = mesh->getVertexNeighborNumber(saddleVert);
      for(valence n = 0; n < nbNeigh; ++n) {
        SimplexId neigh;
        mesh->getVertexNeighbor(saddleVert, n, neigh);

        if(comp_.vertLower(neigh, saddleVert)) {
          if((*mt_data_.ufs)[neigh]
             && (*mt_data_.ufs)[neigh]->find()
                  != (*mt_data_.ufs)[saddleVert]->find()) {
            (*mt_data_.ufs)[saddleVert] = AtomicUF::makeUnion(
              (*mt_data_.ufs)[saddleVert], (*mt_data_.ufs)[neigh]);
          }
        }
      }

      // close arcs on this node
      closeArcsUF(closeNode, (*mt_data_.ufs)[saddleVert]);
    }

    // ------------------------------------------------------------------------

    template <typename scalarType>
    void ftm::FTMTree_MT::sortInput(void) {

      const auto nbVertices = scalars_->size;
      scalars_->sortedVertices.resize(nbVertices);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for(SimplexId i = 0; i < nbVertices; i++) {
        scalars_->sortedVertices[scalars_->offsets[i]] = i;
      }
    }

  } // namespace ftm
} // namespace ttk
// Process

#endif /* end of include guard: FTMTREE_MT_TPL_H */
