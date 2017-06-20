/// \ingroup baseCode
/// \class ttk:FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#include "FTMTree_MT.h"

// -----------
// CONSTRUCT
// -----------
// {

FTMTree_MT::FTMTree_MT(Params *const params, Triangulation *mesh, Scalars *const scalars,
                     TreeType type)
    : params_(params), mesh_(mesh), scalars_(scalars)
{
   treeData_.treeType = type;

   treeData_.superArcs     = nullptr;
   treeData_.nodes         = nullptr;
   treeData_.roots         = nullptr;
   treeData_.leaves        = nullptr;
   treeData_.vert2tree     = nullptr;
   treeData_.trunkSegments = nullptr;
   treeData_.visitOrder    = nullptr;
   treeData_.ufs           = nullptr;
   treeData_.propagation   = nullptr;
   treeData_.valences      = nullptr;
   treeData_.openedNodes   = nullptr;

#ifdef withStatsTime
   treeData_.arcStart = nullptr;
   treeData_.arcEnd   = nullptr;
   treeData_.arcOrig  = nullptr;
   treeData_.arcTasks = nullptr;
#endif
}

FTMTree_MT::~FTMTree_MT()
{
   delete treeData_.superArcs;
   delete treeData_.nodes;
   delete treeData_.roots;
   delete treeData_.leaves;
   delete treeData_.vert2tree;
   delete treeData_.trunkSegments;
   delete treeData_.visitOrder;
   delete treeData_.ufs;
   delete treeData_.propagation;
   delete treeData_.valences;
   delete treeData_.openedNodes;

#ifdef withStatsTime
   delete treeData_.arcStart;
   delete treeData_.arcEnd;
   delete treeData_.arcOrig;
   delete treeData_.arcTasks;
#endif

   treeData_.superArcs     = nullptr;
   treeData_.nodes         = nullptr;
   treeData_.roots         = nullptr;
   treeData_.leaves        = nullptr;
   treeData_.vert2tree     = nullptr;
   treeData_.trunkSegments = nullptr;
   treeData_.visitOrder    = nullptr;
   treeData_.ufs           = nullptr;
   treeData_.propagation   = nullptr;
   treeData_.valences      = nullptr;
   treeData_.openedNodes   = nullptr;

#ifdef withStatsTime
   treeData_.arcStart = nullptr;
   treeData_.arcEnd   = nullptr;
   treeData_.arcOrig  = nullptr;
   treeData_.arcTasks = nullptr;
#endif
}

// }
// -------
// Process
// -------
// {

void FTMTree_MT::build(const bool ct)
{
    string treeString;
   // --------------------------
   // Comparator init (template)
   // --------------------------
   if (treeData_.treeType == TreeType::Join) {
      treeString = "JT";
      comp_.vertLower = [this](idVertex a, idVertex b) -> bool {
         return this->scalars_->isLower(a, b);
      };
      comp_.vertHigher = [this](idVertex a, idVertex b) -> bool {
         return this->scalars_->isHigher(a, b);
      };
   } else {
      treeString = "ST";
      comp_.vertLower = [this](idVertex a, idVertex b) -> bool {
         return this->scalars_->isHigher(a, b);
      };
      comp_.vertHigher = [this](idVertex a, idVertex b) -> bool {
         return this->scalars_->isLower(a, b);
      };
   }

   // ----------------------------
   // Build Merge treeString using tasks
   // ----------------------------
   DebugTimer precomputeTime;
   int alreadyDone = precompute();
   printTime(precomputeTime, "3 precompute " + treeString, scalars_->size, 2 + alreadyDone);

   DebugTimer buildTime;
   leaves();
   int nbProcessed = 0;
#ifdef withProcessSpeed
   // count process
   for (int i = 0; i < scalars_->size; i++) {
       if((*treeData_.vert2tree)[i] != nullCorresp)
           ++nbProcessed;
   }
#endif
   // TODO continue + patch add. mat
   printTime(buildTime, "4 leaves "+treeString, nbProcessed);

   DebugTimer bbTime;
   idVertex bbSize = trunk(ct);
   printTime(bbTime, "5 trunk "+treeString, bbSize);

   // ------------
   // Segmentation
   // ------------
   if (ct) {
      DebugTimer segmTime;
      buildSegmentation();
      printTime(segmTime, "6 segmentation " + treeString, scalars_->size);
   }

   printTree2();
}

// extrema

int FTMTree_MT::precompute()
{
   int ret = 0;
   // if not already computed by CT
   if(getNumberOfNodes() == 0){
      const auto nbScalars = scalars_->size;
      const auto chunkSize = getChunkSize();
      const auto chunkNb   = getChunkCount();

      // --------------------------------
      // Extrema extract and launch tasks
      // --------------------------------
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#pragma omp task firstprivate(chunkId)
         {
            const idVertex lowerBound = chunkId * chunkSize;
            const idVertex upperBound = min(nbScalars, (chunkId + 1) * chunkSize);
            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const auto &neighNumb = mesh_->getVertexNeighborNumber(v);
               valence     val       = 0;

               for (valence n = 0; n < neighNumb; ++n) {
                  idVertex neigh;
                  mesh_->getVertexNeighbor(v, n, neigh);
                  comp_.vertLower(neigh, v) && ++val;
               }

               (*treeData_.valences)[v] = val;

               if (!val) {
                  makeNode(v);
               }
            }
         }
      }

#pragma omp taskwait
   } else {
       ret = 1;
   }

   // fill leaves
   const auto& nbLeaves = treeData_.nodes->size();
   treeData_.leaves->resize(nbLeaves);
   std::iota(treeData_.leaves->begin(), treeData_.leaves->end(), 0);

   if (debugLevel_ >= 3) {
      cout << "nb leaves " << nbLeaves << endl;
   }

   // -------
   // Reserve Arcs
   // -------
   treeData_.superArcs->reserve(nbLeaves * 2 + 1);
#ifdef withStatsTime
   createVector<float>(treeData_.arcStart);
   createVector<float>(treeData_.arcEnd);
   createVector<idVertex>(treeData_.arcOrig);
   createVector<idNode>(treeData_.arcTasks);
   treeData_.arcStart->resize(nbLeaves*2 +1,0);
   treeData_.arcEnd->resize(nbLeaves*2 +1,0);
   treeData_.arcOrig->resize(nbLeaves*2 +1,0);
   treeData_.arcTasks->resize(nbLeaves*2 +1,0);
#endif

   return ret;
}

// skeleton

DebugTimer _launchGlobalTime;

void FTMTree_MT::leaves()
{
   _launchGlobalTime.reStart();

   const auto &nbLeaves = treeData_.leaves->size();

   // elevation: backbone only
   if (nbLeaves == 1) {
      const idVertex v            = (*treeData_.nodes)[0].getVertexId();
      (*treeData_.openedNodes)[v] = 1;
      (*treeData_.ufs)[v]         = new AtomicUF(v);
      return;
   }

   treeData_.activeTasks = nbLeaves;

   // Need testing, simulate priority
   // best with gcc
   auto comp = [this](const idNode a, const idNode b) {
      return this->comp_.vertLower(this->getNode(a)->getVertexId(),
                                   this->getNode(b)->getVertexId());
   };
   sort(treeData_.leaves->begin(), treeData_.leaves->end(), comp);

   for (idNode n = 0; n < nbLeaves; ++n)
   {
      const idNode l = (*treeData_.leaves)[n];
      int          v = getNode(l)->getVertexId();
      // for each node: get vert, create uf and lauch
      (*treeData_.ufs)[v] = new AtomicUF(v);

#pragma omp task untied
      processTask(v, n);
   }

#pragma omp taskwait
}

void FTMTree_MT::processTask(const idVertex startVert, const idVertex orig)
{

   // ------------------------
   // current task id / propag
   // ------------------------
   // local order (ignore non regular verts)
   idVertex localOrder = -1;
   UF startUF = (*treeData_.ufs)[startVert]->find();
   // get or recover states
   CurrentState *currentState;
   if (startUF->getNbStates()) {
      currentState = startUF->getFirstState();
   } else {
      currentState = new CurrentState(startVert, comp_.vertHigher);
      startUF->addState(currentState);
   }

   currentState->addNewVertex(startVert);

   // avoid duplicate processing of startVert
   bool seenFirst = false;

   // -----------
   // ARC OPENING
   // -----------
   idNode     startNode  = getCorrespondingNodeId(startVert);
   idSuperArc currentArc = openSuperArc(startNode);
   startUF->addArcToClose(currentArc);
#ifdef withStatsTime
   (*treeData_.arcStart)[currentArc] = _launchGlobalTime.getElapsedTime();
   (*treeData_.arcOrig)[currentArc] = orig;
#endif

   // ----------------
   // TASK PROPAGATION
   // ----------------
   while (!currentState->empty()) {
      // -----------
      // Next vertex
      // -----------
      idVertex currentVert = currentState->getNextMinVertex();

      // ignore duplicate
      if (!isCorrespondingNull(currentVert) && !isCorrespondingNode(currentVert)) {
         continue;
      } else {
         // first node can be duplicate, avoid duplicate process
         if (currentVert == startVert) {
            if (!seenFirst) {
               seenFirst = true;
            } else {
               continue;
            }
         }
      }

      // local order to avoid sort
      (*treeData_.visitOrder)[currentVert] = localOrder++;

      // -------------------------------------
      // Saddle & Last detection + propagation
      // -------------------------------------
      bool isSaddle, isLast;
      tie(isSaddle, isLast) = propage(*currentState, startUF);

      // regular propagation
#pragma omp atomic write seq_cst
      (*treeData_.ufs)[currentVert] = startUF;

      // -----------
      // Saddle case
      // -----------
      if (isSaddle) {

# ifdef withStatsTime
         (*treeData_.arcEnd)[currentArc] = _launchGlobalTime.getElapsedTime();
         (*treeData_.arcTasks)[currentArc] = treeData_.activeTasks;
# endif
         // need a node on this vertex
         (*treeData_.openedNodes)[currentVert] = 1;

         // ---------------------------
         // If last close all and merge
         // ---------------------------
         if (isLast) {
            // last task detection
            idNode remainingTasks;
#pragma omp atomic read seq_cst
            remainingTasks = treeData_.activeTasks;
            if (remainingTasks == 1) {
                // only backbone remaining
                return;
            }

             // finish works here
            closeAndMergeOnSaddle(currentVert);

            // made a node on this vertex
#pragma omp atomic write seq_cst
            (*treeData_.openedNodes)[currentVert] = 0;

            // recursively continue
#pragma omp taskyield
            processTask(currentVert, orig);
         } else {
            // Active tasks / threads
#pragma omp atomic update seq_cst
            treeData_.activeTasks--;
         }

         // stop at saddle
         return;
      }

      if (currentVert != startVert) {
         updateCorrespondingArc(currentVert, currentArc);
      }
      getSuperArc(currentArc)->setLastVisited(currentVert);

   }  // end wile propagation

// ----------
// close root
// ----------
   const idVertex closeVert      = getSuperArc(currentArc)->getLastVisited();
   bool           existCloseNode = isCorrespondingNode(closeVert);
   idNode closeNode = (existCloseNode) ? getCorrespondingNodeId(closeVert) : makeNode(closeVert);
   closeSuperArc(currentArc, closeNode);
   getSuperArc(currentArc)->decrNbSeen();
   idNode rootPos              = treeData_.roots->getNext();
   (*treeData_.roots)[rootPos] = closeNode;

#ifdef withStatsTime
   (*treeData_.arcEnd)[currentArc] = _launchGlobalTime.getElapsedTime();
#endif
}

tuple<bool, bool> FTMTree_MT::propage(CurrentState &currentState, UF curUF)
{
   bool        becameSaddle = false, isLast = false;
   const auto &nbNeigh = mesh_->getVertexNeighborNumber(currentState.vertex);
   valence decr = 0;

   // once for all
   auto* curUFF = curUF->find();

   // propagation / is saddle
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(currentState.vertex, n, neigh);

      if (comp_.vertLower(neigh, currentState.vertex)) {
         UF neighUF = (*treeData_.ufs)[neigh];

         // is saddle
         if (!neighUF || neighUF->find() != curUFF) {
            becameSaddle = true;
         } else if (neighUF) {
             ++decr;
         }

      } else {
         if (!(*treeData_.propagation)[neigh] ||
             (*treeData_.propagation)[neigh]->find() != curUFF) {
            currentState.addNewVertex(neigh);
            (*treeData_.propagation)[neigh] = curUFF;
         }
      }
   }

   // is last
   valence  oldVal;
#pragma omp atomic capture
   {
      oldVal = (*treeData_.valences)[currentState.vertex];
      (*treeData_.valences)[currentState.vertex] -= decr;
   }
   if (oldVal == decr) {
      isLast = true;
   }

   return make_tuple(becameSaddle, isLast);
}

void FTMTree_MT::closeAndMergeOnSaddle(idVertex saddleVert)
{
   idNode closeNode = makeNode(saddleVert);

   // Union of the UF coming here (merge propagation and closing arcs)
   const auto &nbNeigh = mesh_->getVertexNeighborNumber(saddleVert);
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(saddleVert, n, neigh);

      if (comp_.vertLower(neigh, saddleVert)) {
         if ((*treeData_.ufs)[neigh]->find() != (*treeData_.ufs)[saddleVert]->find()) {

            (*treeData_.ufs)[saddleVert] =
                AtomicUF::makeUnion((*treeData_.ufs)[saddleVert], (*treeData_.ufs)[neigh]);
         }
      }
   }

   // close arcs on this node
   closeArcsUF(closeNode, (*treeData_.ufs)[saddleVert]);

   (*treeData_.ufs)[saddleVert]->find()->mergeStates();
   (*treeData_.ufs)[saddleVert]->find()->setExtrema(saddleVert);
}

void FTMTree_MT::closeOnBackBone(idVertex saddleVert)
{
   idNode closeNode = makeNode(saddleVert);

   // Union of the UF coming here (merge propagation and closing arcs)
   const auto &nbNeigh = mesh_->getVertexNeighborNumber(saddleVert);
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(saddleVert, n, neigh);

      if (comp_.vertLower(neigh, saddleVert)) {
         if ((*treeData_.ufs)[neigh] &&
             (*treeData_.ufs)[neigh]->find() != (*treeData_.ufs)[saddleVert]->find()) {

            (*treeData_.ufs)[saddleVert] =
                AtomicUF::makeUnion((*treeData_.ufs)[saddleVert], (*treeData_.ufs)[neigh]);
         }
      }
   }

   // close arcs on this node
   closeArcsUF(closeNode, (*treeData_.ufs)[saddleVert]);
}

void FTMTree_MT::closeArcsUF(idNode closeNode, UF uf)
{
   for (const auto &sa : uf->find()->getOpenedArcs()) {
      closeSuperArc(sa, closeNode);
   }
   uf->find()->clearOpenedArcs();
}

idVertex FTMTree_MT::trunk(const bool ct)
{
   DebugTimer bbTimer;

   vector<idVertex> pendingNodesVerts;
   const auto &     nbScalars = scalars_->size;

   // -----
   //pendingNodesVerts
   // -----
  pendingNodesVerts.reserve(max(10, nbScalars / 500));
   for (idVertex v = 0; v < nbScalars; ++v) {
       if((*treeData_.openedNodes)[v]){
          pendingNodesVerts.emplace_back(v);
       }
   }
   sort(pendingNodesVerts.begin(), pendingNodesVerts.end(), comp_.vertLower);
   for (const idVertex v : pendingNodesVerts) {
      closeOnBackBone(v);
   }

   // ----
   // Arcs
   // ----

   const auto &nbNodes =pendingNodesVerts.size();
   for (idNode n = 1; n < nbNodes; ++n) {
      idSuperArc na =
          makeSuperArc(getCorrespondingNodeId(pendingNodesVerts[n - 1]), getCorrespondingNodeId(pendingNodesVerts[n]));
      getSuperArc(na)->setLastVisited(pendingNodesVerts[n]);
   }

   if (!nbNodes) {
      return 0;
   }
   const idSuperArc lastArc = openSuperArc(getCorrespondingNodeId(pendingNodesVerts[nbNodes - 1]));

   // debug close Root
   const idNode rootNode = makeNode((*scalars_->sortedVertices)[(isJT())?scalars_->size -1:0]);
   closeSuperArc(lastArc, rootNode);
   getSuperArc(lastArc)->setLastVisited(getNode(rootNode)->getVertexId());

   printTime(bbTimer, "Backbone seq.", -1, 3);
   bbTimer.reStart();

// ------------
// Segmentation
// ------------
   // bounds
   idVertex begin, stop, processed;
   tie(begin, stop) = getBoundsFromVerts(pendingNodesVerts);
   cout << "trunk range " << abs(stop-begin) << endl;
   if(ct){
       processed = trunkCTSegmentation(pendingNodesVerts, begin, stop);
   } else {
       processed = trunkSegmentation(pendingNodesVerts, begin, stop);
   }
   printTime(bbTimer, "Backbone para.", -1, 3);

   // ---------------------
   // Root (close last arc)
   // ---------------------
   // if several CC still the backbone is only in one.
   // But the root may not be the max node of the whole dataset
   return processed;
}

idVertex FTMTree_MT::trunkSegmentation(const vector<idVertex> &pendingNodesVerts,
                                      const idVertex begin, const idVertex stop)
{
   // Assign missing vert to the good arc
   // and also add the corresponding number for
   // futur arc reserve
   const int nbTasksThreads = 40;
   const auto sizeBackBone  = abs(stop - begin);
   const auto chunkSize     = getChunkSize(sizeBackBone, nbTasksThreads);
   const auto chunkNb       = getChunkCount(sizeBackBone, nbTasksThreads);
   // si pas efficace vecteur de la taille de node ici a la place de acc
   idNode   lastVertInRange = 0;
   idVertex acc             = 0;
   idVertex tot             = 0;
   if (isJT()) {
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#pragma omp task firstprivate(chunkId, acc, lastVertInRange) shared(pendingNodesVerts, tot)
         {
            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));
            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = (*scalars_->sortedVertices)[v];
               if (isCorrespondingNull(s)) {
                  const idNode oldVertInRange = lastVertInRange;
                  lastVertInRange          = getVertInRange(pendingNodesVerts, s, lastVertInRange);
                  const idSuperArc thisArc = upArcFromVert(pendingNodesVerts[lastVertInRange]);
                  updateCorrespondingArc(s, thisArc);
                  if (oldVertInRange == lastVertInRange) {
                     ++acc;
                  } else {
                     // accumulated to have only one atomic update when needed
                     const idSuperArc oldArc = upArcFromVert(pendingNodesVerts[oldVertInRange]);
                     getSuperArc(oldArc)->atomicIncVisited(acc);
#ifdef withProcessSpeed
#pragma omp atomic update
                     tot += acc;
#endif
                     acc = 1;
                  }
               }
            }
            // force increment last arc
            const idNode     baseNode = getCorrespondingNodeId(pendingNodesVerts[lastVertInRange]);
            const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
            getSuperArc(upArc)->atomicIncVisited(acc);
#ifdef withProcessSpeed
#pragma omp atomic update
            tot += acc;
#endif
         } // end task
      }
   } else {
      for (idVertex chunkId = chunkNb - 1; chunkId >= 0; --chunkId) {
#pragma omp task firstprivate(chunkId, acc, lastVertInRange) shared(pendingNodesVerts, tot)
         {
            const idVertex upperBound = begin - chunkId * chunkSize;
            const idVertex lowerBound = max(stop, begin - (chunkId + 1) * chunkSize);
            for (idVertex v = upperBound; v > lowerBound; --v) {
               const idVertex s = (*scalars_->sortedVertices)[v];
               if (isCorrespondingNull(s)) {
                  const idNode oldVertInRange = lastVertInRange;
                  lastVertInRange          = getVertInRange(pendingNodesVerts, s, lastVertInRange);
                  const idSuperArc thisArc = upArcFromVert(pendingNodesVerts[lastVertInRange]);
                  updateCorrespondingArc(s, thisArc);
                  if (oldVertInRange == lastVertInRange) {
                     ++acc;
                  } else {
                     // accumulated to have only one atomic update when needed
                     const idSuperArc oldArc = upArcFromVert(pendingNodesVerts[oldVertInRange]);
                     getSuperArc(oldArc)->atomicIncVisited(acc);
#ifdef withProcessSpeed
#pragma omp atomic update
                     tot += acc;
#endif
                     acc = 1;
                  }
               }
            }
            // force increment last arc
            const idNode     baseNode = getCorrespondingNodeId(pendingNodesVerts[lastVertInRange]);
            const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
            getSuperArc(upArc)->atomicIncVisited(acc);
#ifdef withProcessSpeed
#pragma omp atomic update
            tot += acc;
#endif
         }
      }
   }
#pragma omp taskwait
   return tot;
}

idVertex FTMTree_MT::trunkCTSegmentation(const vector<idVertex> &pendingNodesVerts,
                                        const idVertex begin, const idVertex stop)
{
   const int nbTasksThreads = 40;
   const auto sizeBackBone  = abs(stop - begin);
   const auto chunkSize     = getChunkSize(sizeBackBone, nbTasksThreads);
   const auto chunkNb       = getChunkCount(sizeBackBone, nbTasksThreads);
   // si pas efficace vecteur de la taille de node ici a la place de acc
   idNode   lastVertInRange = 0;
   treeData_.trunkSegments->resize(getNumberOfSuperArcs());
   if (isJT()) {
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#pragma omp task firstprivate(chunkId, lastVertInRange) shared(pendingNodesVerts)
         {
            vector<idVertex> regularList;
            regularList.reserve(25);
            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));
            lastVertInRange =
                getVertInRange(pendingNodesVerts, (*scalars_->sortedVertices)[lowerBound], 0);
            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = (*scalars_->sortedVertices)[v];
               if (isCorrespondingNull(s)) {
                  const idNode oldVertInRange = lastVertInRange;
                  lastVertInRange          = getVertInRange(pendingNodesVerts, s, lastVertInRange);
                  const idSuperArc thisArc = upArcFromVert(pendingNodesVerts[lastVertInRange]);
                  updateCorrespondingArc(s, thisArc);
                  if (oldVertInRange == lastVertInRange) {
                     regularList.emplace_back(s);
                  } else {
                     // accumulated to have only one atomic update when needed
                     const idSuperArc oldArc = upArcFromVert(pendingNodesVerts[oldVertInRange]);
                     if (regularList.size()) {
#pragma omp critical
                        {
                           (*treeData_.trunkSegments)[oldArc].emplace_back(regularList);
                           regularList.clear();
                        }
                        regularList.emplace_back(s);
                     }
                  }
               }
            }
            // force increment last arc
            const idNode     baseNode = getCorrespondingNodeId(pendingNodesVerts[lastVertInRange]);
            const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
            if (regularList.size()) {
#pragma omp critical
               {
                  (*treeData_.trunkSegments)[upArc].emplace_back(regularList);
                  regularList.clear();
               }
            }
         }
      }
   } else {
      for (idVertex chunkId = chunkNb - 1; chunkId >= 0; --chunkId) {
#pragma omp task firstprivate(chunkId, lastVertInRange) shared(pendingNodesVerts)
         {
            vector<idVertex> regularList;
            regularList.reserve(25);
            const idVertex upperBound = begin - chunkId * chunkSize;
            const idVertex lowerBound = max(stop, begin - (chunkId + 1) * chunkSize);
            lastVertInRange =
                getVertInRange(pendingNodesVerts, (*scalars_->sortedVertices)[upperBound], 0);
            for (idVertex v = upperBound; v > lowerBound; --v) {
               const idVertex s = (*scalars_->sortedVertices)[v];
               if (isCorrespondingNull(s)) {
                  const idNode oldVertInRange = lastVertInRange;
                  lastVertInRange          = getVertInRange(pendingNodesVerts, s, lastVertInRange);
                  const idSuperArc thisArc = upArcFromVert(pendingNodesVerts[lastVertInRange]);
                  updateCorrespondingArc(s, thisArc);
                  if (oldVertInRange == lastVertInRange) {
                     regularList.emplace_back(s);
                  } else {
                     // accumulated to have only one atomic update when needed
                     const idSuperArc oldArc = upArcFromVert(pendingNodesVerts[oldVertInRange]);
                     if (regularList.size()) {
#pragma omp critical
                        {
                           (*treeData_.trunkSegments)[oldArc].emplace_back(regularList);
                           regularList.clear();
                        }
                     }
                     regularList.emplace_back(s);
                  }
               }
            }
            // force increment last arc
            const idNode baseNode  = getCorrespondingNodeId(pendingNodesVerts[lastVertInRange]);
            const idSuperArc upArc = getNode(baseNode)->getUpSuperArcId(0);
            if (regularList.size()) {
#pragma omp critical
               {
                  (*treeData_.trunkSegments)[upArc].emplace_back(regularList);
                  regularList.clear();
               }
            }
         }
      }
   }
#pragma omp taskwait
   // count added
   idVertex tot = 0;
#ifdef withProcessSpeed
   for (const auto& l : *treeData_.trunkSegments) {
       idVertex arcSize = 0;
       for (const auto& v: l){
          arcSize += v.size();
       }
       tot += arcSize;
   }
#endif
   return tot;
}

// segmentation

void FTMTree_MT::buildSegmentation()
{

   const idSuperArc nbArcs = treeData_.superArcs->size();

   // ------------
   // Make reserve
   // ------------
   // SuperArc i correspond to segment i,
   // one arc correspond to one segment
   vector<idVertex> sizes(nbArcs);

   // get the size of each segment
   const idSuperArc arcChunkSize = getChunkSize(nbArcs);
   const idSuperArc arcChunkNb   = getChunkCount(nbArcs);
   for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId){
       // WHY shared(sizes) is needed ??
#pragma omp task firstprivate(arcChunkId) shared(sizes)
       {
           const idSuperArc lowerBound = arcChunkId*arcChunkSize;
           const idSuperArc upperBound = min(nbArcs, (arcChunkId+1)*arcChunkSize );
           for (idSuperArc a = lowerBound; a < upperBound; ++a) {
              sizes[a] = max(0, (*treeData_.superArcs)[a].getNbVertSeen() - 1);
           }
       }
   }
#pragma omp taskwait

   // change segments size using the created vector
   treeData_.segments_.resize(sizes);

   DebugTimer segmentsSet;
   // -----------------------------
   // Fill segments using vert2tree
   // -----------------------------
   // current status of the segmentation of this arc
   vector<idVertex> posSegm(nbArcs, 0);

   // Segments are connex region of geometrie forming
   // the segmentation (sorted in ascending order)
   const idVertex nbVert = scalars_->size;
   const idVertex chunkSize = getChunkSize();
   const idVertex chunkNb   = getChunkCount();
   for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId)
   {
#pragma omp task firstprivate(chunkId) shared(posSegm)
       {
          const idVertex lowerBound = chunkId * chunkSize;
          const idVertex upperBound = min(nbVert, (chunkId+1)*chunkSize);
          for (idVertex i = lowerBound; i < upperBound; ++i) {
             const auto vert = (*scalars_->sortedVertices)[i];
             if (isCorrespondingArc(vert)) {
                idSuperArc sa = getCorrespondingSuperArcId(vert);
                idVertex   vertToAdd;
                if((*treeData_.visitOrder)[vert] != nullVertex){
                   // Opposite order for Split Tree
                   vertToAdd = (*treeData_.visitOrder)[vert];
                   if(isST()) vertToAdd = getSuperArc(sa)->getNbVertSeen() - vertToAdd -2;
                   treeData_.segments_[sa][vertToAdd] = vert;
                } else if (treeData_.trunkSegments->size() == 0){
                    // MT computation
#pragma omp atomic capture
                   vertToAdd = posSegm[sa]++;
                   treeData_.segments_[sa][vertToAdd] = vert;
                }

             }  // end is arc
          } // end for
       } // end task
   }
#pragma omp taskwait

   printTime(segmentsSet, "segm. set verts", -1, 3);

   if (treeData_.trunkSegments->size() == 0) {
      // sort arc that have been filled by the trunk
      // only for MT
      DebugTimer segmentsSortTime;
      for (idSuperArc a = 0; a < nbArcs; ++a) {
         if (posSegm[a]) {
#pragma omp task firstprivate(a)
            treeData_.segments_[a].sort(scalars_);
         }
      }
#pragma omp taskwait
      printTime(segmentsSortTime, "segm. sort verts", -1, 3);
   } else {
       // Contour tree: we create the arc segmentation for arcs in the trunk
       DebugTimer segmentsArcTime;
       for (idSuperArc a = 0; a < nbArcs; ++a) {
          // CT computation, we have already the vert list
          if ((*treeData_.trunkSegments)[a].size()) {
#pragma omp task firstprivate(a)
             treeData_.segments_[a].createFromList(scalars_, (*treeData_.trunkSegments)[a],
                                                   treeData_.treeType == TreeType::Split);
          }
       }
#pragma omp taskwait
      printTime(segmentsArcTime, "segm. creat trunk verts", -1, 3);
   }

   // ----------------------
   // Update SuperArc region
   // ----------------------
   // ST have a segmentation wich is in the reverse-order of its build
   // ST have a segmentation sorted in ascending order as JT
   for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId){
#pragma omp task firstprivate(arcChunkId)
      {
         const idSuperArc lowerBound = arcChunkId * arcChunkSize;
         const idSuperArc upperBound = min(nbArcs, (arcChunkId + 1) * arcChunkSize);
         for (idSuperArc a = lowerBound; a < upperBound; ++a) {
            // avoid empty region
            if (treeData_.segments_[a].size()) {
               (*treeData_.superArcs)[a].concat(treeData_.segments_[a].begin(),
                                                treeData_.segments_[a].end());
            }
         }

      }
   }
#pragma omp taskwait
}

// }
// ---------------------------
// Arcs and node manipulations
// ---------------------------
// {
// SuperArcs
// .......................{
idSuperArc FTMTree_MT::openSuperArc(idNode downNodeId)
{
#ifndef withKamikaze
   if (downNodeId < 0 || (size_t)downNodeId >= getNumberOfNodes()) {
      cout << "[Merge Tree] openSuperArc on a inexisting node !" << endl;
      return -2;
   }
#endif

   idSuperArc newSuperArcId = treeData_.superArcs->getNext();
   (*treeData_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
   (*treeData_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);

   return newSuperArcId;
}

idSuperArc FTMTree_MT::makeSuperArc(idNode downNodeId, idNode upNodeId)

{
   idSuperArc newSuperArcId = treeData_.superArcs->getNext();
   (*treeData_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
   (*treeData_.superArcs)[newSuperArcId].setUpNodeId(upNodeId);

   (*treeData_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);
   (*treeData_.nodes)[upNodeId].addDownSuperArcId(newSuperArcId);

   return newSuperArcId;
}

void FTMTree_MT::closeSuperArc(idSuperArc superArcId, idNode upNodeId)
{
#ifndef withKamikaze

   if (superArcId < 0 || (size_t)superArcId >= getNumberOfSuperArcs()) {
      cout << "[Merge Tree] closeSuperArc on a inexisting arc !" << endl;
      return;
   }

   if (upNodeId < 0 || (size_t)upNodeId >= getNumberOfNodes()) {
      cout << "[Merge Tree] closeOpenedArc on a inexisting node !" << endl;
      return;
   }

#endif
   (*treeData_.superArcs)[superArcId].setUpNodeId(upNodeId);
   (*treeData_.nodes)[upNodeId].addDownSuperArcId(superArcId);
}


//   }
// nodes
// .....................{

vector<idNode> FTMTree_MT::sortedNodes(const bool para)
{
   vector<idNode> sortedNodes(treeData_.nodes->size());
   std::iota(sortedNodes.begin(), sortedNodes.end(), 0);

   auto indirect_sort = [&](const idNode a, const idNode b) {
      return comp_.vertLower(getNode(a)->getVertexId(), getNode(b)->getVertexId());
   };

   if (para) {
#ifdef __clang__
      std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#else
      __gnu_parallel::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#endif
   } else {
#pragma omp single
      {
         std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
      }
   }

   return sortedNodes;
}

// add

idNode FTMTree_MT::makeNode(idVertex vertexId, idVertex term)
{
#ifndef withKamikaze
   if (vertexId < 0 || vertexId >= scalars_->size) {
      cout << "[Merge Tree] make node, wrong vertex :" << vertexId << " on " << scalars_->size
           << endl;
      return -1;
   }
#endif

   if (isCorrespondingNode(vertexId)) {
      return getCorrespondingNodeId(vertexId);
   }

   idNode newNodeId = treeData_.nodes->getNext();
   (*treeData_.nodes)[newNodeId].setVertexId(vertexId);
   (*treeData_.nodes)[newNodeId].setTerminaison(term);
   updateCorrespondingNode(vertexId, newNodeId);

   return newNodeId;
}

idNode FTMTree_MT::makeNode(const Node *const n, idVertex term)
{
   return makeNode(n->getVertexId());
}

// Normal insert : existing arc stay below inserted (JT example)
//  *   - <- upNodeId
//  | \ |   <- newSA
//  |   * <- newNodeId
//  |   |   <- currentSA
//  - - -
idSuperArc FTMTree_MT::insertNode(Node *node, const bool segm)
{
   // already present
   if (isCorrespondingNode(node->getVertexId())) {
      Node *myNode = vertex2Node(node->getVertexId());
      // If it has been hidden / replaced we need to re-make it
      idSuperArc correspondingArcId = myNode->getUpSuperArcId(0);
      updateCorrespondingArc(myNode->getVertexId(), correspondingArcId);
   }

   idNode     upNodeId, newNodeId;
   idSuperArc currentSA, newSA;
   idVertex   origin;

   // Create new node
   currentSA = getCorrespondingSuperArcId(node->getVertexId());
   upNodeId  = (*treeData_.superArcs)[currentSA].getUpNodeId();
   origin    = (*treeData_.nodes)[(*treeData_.superArcs)[currentSA].getDownNodeId()].getOrigin();
   newNodeId = makeNode(node, origin);

   // Connectivity
   // Insert only node inside the partition : created arc don t cross
   newSA = makeSuperArc(newNodeId, upNodeId);

   (*treeData_.superArcs)[currentSA].setUpNodeId(newNodeId);
   (*treeData_.nodes)[upNodeId].removeDownSuperArc(currentSA);
   (*treeData_.nodes)[newNodeId].addDownSuperArcId(currentSA);

   // cut the vertex list at the node position and
   // give each arc its part.
   if (segm) {
      if (treeData_.treeType == TreeType::Split) {
         (*treeData_.superArcs)[newSA].concat(
             get<1>((*treeData_.superArcs)[currentSA].splitBack(node->getVertexId(), scalars_)));
      } else {
         (*treeData_.superArcs)[newSA].concat(
             get<1>((*treeData_.superArcs)[currentSA].splitFront(node->getVertexId(), scalars_)));
      }
   }

   return newSA;
}

// traverse

Node *FTMTree_MT::getDownNode(const SuperArc *a)
{
   return &((*treeData_.nodes)[a->getDownNodeId()]);
}

Node *FTMTree_MT::getUpNode(const SuperArc *a)
{
   return &((*treeData_.nodes)[a->getUpNodeId()]);
}

idNode FTMTree_MT::getDownNodeId(const SuperArc *a)
{
   return a->getDownNodeId();
}

idNode FTMTree_MT::getUpNodeId(const SuperArc *a)
{
   return a->getUpNodeId();
}

Node *FTMTree_MT::getLowerNode(const SuperArc *a)
{
   if (isST())
      return getUpNode(a);

   return getDownNode(a);
}

Node *FTMTree_MT::getUpperNode(const SuperArc *a)
{
   if (isST())
      return getDownNode(a);

   return getUpNode(a);
}

idNode FTMTree_MT::getLowerNodeId(const SuperArc *a)
{
   if (isST())
      return getUpNodeId(a);

   return getDownNodeId(a);
}

idNode FTMTree_MT::getUpperNodeId(const SuperArc *a)
{
   if (isST())
      return getDownNodeId(a);

   return getUpNodeId(a);
}

// remove

void FTMTree_MT::delNode(idNode node)
{
   Node *mainNode = getNode(node);

   if (mainNode->getNumberOfUpSuperArcs() == 0) {
// -----------------
// Root: No Superarc
// -----------------

#ifndef withKamikaze
      if (mainNode->getNumberOfDownSuperArcs() != 1) {
         // Root with several childs: impossible /\ .
         cout << endl << "[FTMTree_MT]:delNode won't delete ";
         cout << mainNode->getVertexId() << " (root) with ";
         cout << static_cast<unsigned>(mainNode->getNumberOfDownSuperArcs()) << " down ";
         cout << static_cast<unsigned>(mainNode->getNumberOfUpSuperArcs()) << " up ";
         return;
      }
#endif

      idSuperArc downArc  = mainNode->getDownSuperArcId(0);
      Node *     downNode = getNode((*treeData_.superArcs)[downArc].getDownNodeId());

      downNode->removeUpSuperArc(downArc);
      mainNode->clearDownSuperArcs();

   } else if (mainNode->getNumberOfDownSuperArcs() < 2) {
      // ---------------
      // Have one up arc
      // ---------------

      // We delete the upArc of this node,
      // if there is a down arc, we reattach it to the upNode

      idSuperArc upArc  = mainNode->getUpSuperArcId(0);
      idNode     upId   = (*treeData_.superArcs)[upArc].getUpNodeId();
      Node *     upNode = getNode(upId);

      upNode->removeDownSuperArc(upArc);
      mainNode->clearUpSuperArcs();

      if (mainNode->getNumberOfDownSuperArcs()) {
         // -----------------
         // Have one down arc
         // -----------------

         // Reconnect
         idSuperArc downArc = mainNode->getDownSuperArcId(0);
         (*treeData_.superArcs)[downArc].setUpNodeId(upId);
         upNode->addDownSuperArcId(downArc);
         mainNode->clearDownSuperArcs();

         // Segmentation
         (*treeData_.superArcs)[downArc].concat((*treeData_.superArcs)[upArc]);
      }
   }
#ifndef withKamikaze
   else
      cerr << "delete node with multiple childrens " << endl;
#endif
}

// }
// Segmentation
// ...........................{

void FTMTree_MT::finalizeSegmentation(void)
{
   for (auto &arc : *treeData_.superArcs) {
      arc.createSegmentation(scalars_);
   }
}

//    }
// }
// -------------------------------
// Operators : find, print & clone
// -------------------------------
// {

// Clone
FTMTree_MT *FTMTree_MT::clone() const
{
   FTMTree_MT *newMT = new FTMTree_MT(params_, mesh_, scalars_, treeData_.treeType);

   newMT->treeData_.superArcs = treeData_.superArcs;
   newMT->treeData_.nodes     = treeData_.nodes;
   newMT->treeData_.leaves    = treeData_.leaves;
   newMT->treeData_.roots     = treeData_.roots;
   newMT->treeData_.vert2tree = treeData_.vert2tree;

   return newMT;
}

void FTMTree_MT::clone(const FTMTree_MT *mt)
{
   // we already have common data
   treeData_.superArcs = mt->treeData_.superArcs;
   treeData_.nodes     = mt->treeData_.nodes;
   treeData_.leaves    = mt->treeData_.leaves;
   treeData_.roots     = mt->treeData_.roots;
   treeData_.vert2tree = mt->treeData_.vert2tree;
}

// Print
string FTMTree_MT::printArc(idSuperArc a)
{
   const SuperArc *sa = getSuperArc(a);
   stringstream    res;
   res << a;
   res << " : ";
   res << getNode(sa->getDownNodeId())->getVertexId() << " -- ";
   res << getNode(sa->getUpNodeId())->getVertexId();

   res.seekg(0, ios::end);
   while (res.tellg() < 25) {
      res << " ";
      res.seekg(0, ios::end);
   }
   res.seekg(0, ios::beg);

   res << "segm #" << sa->regionSize() << " / " << scalars_->size;  // << " -> ";

   res.seekg(0, ios::end);

   while (res.tellg() < 45) {
      res << " ";
      res.seekg(0, ios::end);
   }
   res.seekg(0, ios::beg);

   res << sa->printReg();
   return res.str();
}

string FTMTree_MT::printNode(idNode n)
{
   const Node * node = getNode(n);
   stringstream res;
   res << n;
   res << " : (";
   res << node->getVertexId() << ") \\ ";

   for (idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); ++i) {
      res << "+";
      res << node->getDownSuperArcId(i) << " ";
   }

   res << " / ";

   for (idSuperArc i = 0; i < node->getNumberOfUpSuperArcs(); ++i) {
      res << "+";
      res << node->getUpSuperArcId(i) << " ";
   }

   return res.str();
}

void FTMTree_MT::printTree2()
{
#ifdef withOpenMP
#pragma omp critical
#endif
   {
      cout << "Nodes----------" << endl;
      for (idNode nid = 0; nid < getNumberOfNodes(); nid++) {
         cout << printNode(nid) << endl;
      }

      cout << "Arcs-----------" << endl;
      for (idSuperArc said = 0; said < getNumberOfSuperArcs(); ++said) {
         cout << printArc(said) << endl;
      }

      cout << "Leaves" << endl;
      for (const auto &l : *treeData_.leaves)
         cout << " " << (*treeData_.nodes)[l].getVertexId();
      cout << endl;

      cout << "Roots" << endl;
      for (const auto &r : *treeData_.roots)
         cout << " " << (*treeData_.nodes)[r].getVertexId();
      cout << endl;
   }
}

void FTMTree_MT::printParams(void) const
{
   if (debugLevel_ > 1) {
      cout << "------------" << endl;
      cout << "nb threads : " << threadNumber_ << endl;
      cout << "debug lvl  : " << debugLevel_ << endl;
      cout << "tree type  : ";
      if (params_->treeType == TreeType::Contour) {
         cout << "Contour";
      } else if (params_->treeType == TreeType::Join) {
         cout << "Join";
      } else if (params_->treeType == TreeType::Split) {
         cout << "Split";
      }
      cout << endl;
      cout << "------------" << endl;
   }
}

int FTMTree_MT::printTime(DebugTimer &t, const string &s, idVertex nbScalars, const int debugLevel) const
{
   if (nbScalars == -1) {
      nbScalars = scalars_->size;
   }

   if (debugLevel_ >= debugLevel) {
      stringstream st;
#ifdef withProcessSpeed
      int          speed = nbScalars / t.getElapsedTime();
#endif
      for (int i = 2; i < debugLevel; i++)
         st << "-";
      st << s << " in ";
      st.seekg(0, ios::end);
      while (st.tellg() < 25) {
         st << " ";
         st.seekg(0, ios::end);
      }
      st.seekg(0, ios::beg);
      st << t.getElapsedTime();

#ifdef withProcessSpeed
      st.seekg(0, ios::end);
      while (st.tellg() < 35) {
         st << " ";
         st.seekg(0, ios::end);
      }
      st.seekg(0, ios::beg);
      st << " at " << speed << " vert/s";
#endif
      cout << st.str() << endl;
   }
   return 1;
}

// }

// ##########
// Protected
// ##########

// -----
// Tools
// -----

idNode FTMTree_MT::getVertInRange(const vector<idVertex> &range, const idVertex v,
                                 const idNode last) const
{
    idNode idRes = last;
    const idNode rangeSize = range.size();
    while (idRes+1 < rangeSize && comp_.vertLower(range[idRes + 1], v)) {
       ++idRes;
    }
    return idRes;
}

tuple<idVertex, idVertex> FTMTree_MT::getBoundsFromVerts(const vector<idVertex> &nodes) const
{
    idVertex begin, stop;

    if(isJT()){
       begin = (*scalars_->mirrorVertices)[nodes[0]];
       stop  = scalars_->size;
    } else {
       begin = (*scalars_->mirrorVertices)[nodes[0]];
       stop  = -1;
    }

    return make_tuple(begin, stop);
}

// ---------
// Operators
// ---------

ostream &ttk::operator<<(ostream &o, SuperArc const &a)
{
   o << a.getDownNodeId() << " <>> " << a.getUpNodeId();
   return o;
}

ostream &ttk::operator<<(ostream &o, Node const &n)
{
   o << n.getNumberOfDownSuperArcs() << " .-. " << n.getNumberOfUpSuperArcs();
   return o;
}

