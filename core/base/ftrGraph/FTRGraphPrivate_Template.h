#ifndef FTRGRAPHPRIVATE_TEMPLATE_H
#define FTRGRAPHPRIVATE_TEMPLATE_H

#include "FTRGraph.h"

// Skeleton + propagation
#ifndef NDEBUG
// #define DEBUG_1(msg) std::cout msg
#define DEBUG_1(msg)
#else
#define DEBUG_1(msg)
#endif

// Dynamic graph structure
#ifndef NDEBUG
// #define DEBUG_2(msg) std::cout msg
#define DEBUG_2(msg)
#else
#define DEBUG_2(msg)
#endif

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localPropagation,
                                                const idSuperArc arcId)
      {
         DEBUG_1(<< "Start " << seed << " go up " << localPropagation->goUp() << std::endl);
#ifndef NDEBUG
         DEBUG_1(<< localPropagation->getRpz() << " " << localPropagation->print() << std::endl);
#endif

         // skeleton
         const idNode     downNode   = graph_.makeNode(seed);
         const idSuperArc currentArc = (arcId != nullSuperArc)? arcId : graph_.openArc(downNode, localPropagation);

#ifndef NDEBUG
         graph_.getArc(currentArc).setFromUp(localPropagation->goUp());
#endif

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle     = false, isSplitSaddle = false;

         // containers
         std::vector<idEdge>                  lowerStarEdges, upperStarEdges;
         std::vector<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localPropagation->empty()) {
            localPropagation->nextVertex();
            const idVertex curVert = localPropagation->getCurVertex();

            // Avoid revisiting things processed by this CC
            if (!graph_.isNode(curVert) && graph_.hasVisited(curVert, localPropagation->getRpz())) {
               DEBUG_1(<< "already seen " << curVert << std::endl);
               continue;
            }

            // Mark this vertex with the current growth
            if (!graph_.isArc(curVert) ) {
               graph_.visit(curVert, currentArc);
               DEBUG_1(<< "visit r: " << curVert << std::endl);
            } else {
               if (graph_.isNode(curVert)) {
                  graph_.visit(curVert, currentArc);
                  DEBUG_1(<< "visit n: " << curVert << std::endl);
               }
               // Caution: crossing tasks can leads to revisit legally
               if (checkSegmentationForArc(seed, curVert, localPropagation)){
                  DEBUG_1(<< "Dismiss current : " << graph_.printArc(currentArc) << std::endl);
                  graph_.getArc(currentArc).hide();
                  return;
               }
            }

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localPropagation);

            lowerComp = lowerComps(lowerStarEdges, localPropagation);
            if(lowerComp.size() > 1){
               isJoinSaddle = true;
            }

            if (isJoinSaddle) {
               isJoinSadlleLast = checkLast(currentArc, localPropagation, lowerStarEdges);
               DEBUG_1(<< ": is join " << isJoinSadlleLast << std::endl);
               // If the current growth reaches a saddle and is not the last
               // reaching this saddle, it just stops here.
               if (!isJoinSadlleLast)
                  break;
            }

            updatePreimage(localPropagation, currentArc);

            upperComp = upperComps(upperStarEdges, localPropagation);
            if (upperComp.size() > 1) {
               DEBUG_1(<< ": is split" << std::endl);
               isSplitSaddle = true;
            }

            if (!isJoinSaddle || isSplitSaddle) {
               // add upper star for futur visit
               localGrowth(localPropagation);
            }
         }

         // if we stop, create/recover the critical point
         const idVertex upVert = localPropagation->getCurVertex();
         const idNode   upNode = updateReebGraph(currentArc, localPropagation);

         // Saddle case

         if (isJoinSadlleLast) {
            localGrowth(localPropagation);
            mergeAtSaddle(upNode, localPropagation);
         }

         if (isSplitSaddle) {
            splitAtSaddle(localPropagation);
         } else if (isJoinSadlleLast) {
            // recursive call
            const idNode     downNode = graph_.getNodeId(upVert);
            const idSuperArc newArc   = graph_.openArc(downNode, localPropagation);
            updateDynGraphCurArc(upVert, newArc, localPropagation);
            graph_.visit(upVert, newArc);
            DEBUG_1(<< "visit m: " << upVert << " with " << newArc << std::endl);
#ifdef TTK_ENABLE_OPENMP
#pragma omp task priority(5)
#endif
            growthFromSeed(upVert, localPropagation, newArc);
         }
      }

      template <typename ScalarType>
      std::pair<std::vector<idEdge>, std::vector<idEdge>> FTRGraph<ScalarType>::visitStar(
          const Propagation* const localPropagation) const
      {
         // TODO re-use the same vectors per thread
         std::vector<idEdge> lowerStar, upperStar;

         const idEdge nbAdjEdges = mesh_.getVertexEdgeNumber(localPropagation->getCurVertex());
         lowerStar.reserve(nbAdjEdges);
         upperStar.reserve(nbAdjEdges);

         for (idEdge e = 0; e < nbAdjEdges; ++e) {
            idEdge edgeId;
            mesh_.getVertexEdge(localPropagation->getCurVertex(), e, edgeId);
            idVertex edgeLowerVert, edgeUpperVert;
            std::tie(edgeLowerVert, edgeUpperVert) = mesh_.getOrderedEdge(edgeId, localPropagation->goUp());
            if (edgeLowerVert == localPropagation->getCurVertex()) {
               upperStar.emplace_back(edgeId);
            } else {
               lowerStar.emplace_back(edgeId);
            }
         }

         return {lowerStar, upperStar};
      }

      template <typename ScalarType>
      std::vector<DynGraphNode<idVertex>*> FTRGraph<ScalarType>::lowerComps(
          const std::vector<idEdge>& finishingEdges, const Propagation* const localProp)
      {
         return dynGraph(localProp).findRoot(finishingEdges);
      }

      template <typename ScalarType>
      std::vector<DynGraphNode<idVertex>*> FTRGraph<ScalarType>::upperComps(
          const std::vector<idEdge>& startingEdges, const Propagation* const localProp)
      {
         return dynGraph(localProp).findRoot(startingEdges);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimage(const Propagation* const localPropagation,
                                                const idSuperArc         curArc)
      {
         const idCell nbAdjTriangles =
             mesh_.getVertexTriangleNumber(localPropagation->getCurVertex());

         DEBUG_1(<< "update preimage " << localPropagation->getCurVertex());

         for (idCell t = 0; t < nbAdjTriangles; ++t) {
            // Classify current cell
            idCell curTriangleid;
            mesh_.getVertexTriangle(localPropagation->getCurVertex(), t, curTriangleid);

            orderedTriangle   oTriangle  = mesh_.getOrderedTriangle(curTriangleid, localPropagation->goUp());
            vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localPropagation);

            // Update DynGraph
            // We can have an end pos on an unvisited triangle
            // in case of saddle points
            switch (curVertPos) {
               case vertPosInTriangle::Start:
                  updatePreimageStartCell(oTriangle, localPropagation, curArc);
                  break;
               case vertPosInTriangle::Middle:
                  updatePreimageMiddleCell(oTriangle, localPropagation, curArc);
                  break;
               case vertPosInTriangle::End:
                  updatePreimageEndCell(oTriangle, localPropagation, curArc);
                  break;
               default:
                  std::cout << "[FTR]: update preimage error, unknown vertPos type" << std::endl;
                  break;
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageStartCell(const orderedTriangle&   oTriangle,
                                                         const Propagation* const localPropagation,
                                                         const idSuperArc         curArc)
      {
         const orderedEdge e0 = mesh_.getOrderedEdge(std::get<0>(oTriangle), localPropagation->goUp());
         const orderedEdge e1 = mesh_.getOrderedEdge(std::get<1>(oTriangle), localPropagation->goUp());
         const idVertex    w  = getWeight(e0, e1, localPropagation);
         bool t = dynGraph(localPropagation).insertEdge(std::get<0>(oTriangle), std::get<1>(oTriangle), w, curArc);

         if (t) {
            DEBUG_2(<< "start add edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< "start no need to create edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageMiddleCell(const orderedTriangle&   oTriangle,
                                                          const Propagation* const localPropagation,
                                                          const idSuperArc         curArc)
      {
         // Check if exist ?
         // If not, the triangle will be visited again once a merge have occured.
         // So we do not add the edge now
         const int t = dynGraph(localPropagation).removeEdge(std::get<0>(oTriangle), std::get<1>(oTriangle));

         // keep history inside the dyngraph structure
         dynGraph(localPropagation).setSubtreeArc(std::get<0>(oTriangle), curArc);
         dynGraph(localPropagation).setSubtreeArc(std::get<1>(oTriangle), curArc);

         if (t) {
            DEBUG_2(<< "mid replace edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
         }
         else {
            DEBUG_2(<< "mid no found edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }

         const orderedEdge e1 = mesh_.getOrderedEdge(std::get<1>(oTriangle), localPropagation->goUp());
         const orderedEdge e2 = mesh_.getOrderedEdge(std::get<2>(oTriangle), localPropagation->goUp());
         const idVertex    w  = getWeight(e1, e2, localPropagation);
         const int u = dynGraph(localPropagation).insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), w, curArc);

         if (u) {
            DEBUG_2(<< " new edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< " mid no need to create edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageEndCell(const orderedTriangle&   oTriangle,
                                                       const Propagation* const localPropagation,
                                                       const idSuperArc         curArc)
      {
         const int t = dynGraph(localPropagation).removeEdge(std::get<1>(oTriangle), std::get<2>(oTriangle));

         dynGraph(localPropagation).setSubtreeArc(std::get<1>(oTriangle), curArc);
         dynGraph(localPropagation).setSubtreeArc(std::get<2>(oTriangle), curArc);

         if (t) {
            DEBUG_2(<< "end remove edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< "end not found edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateDynGraphCurArc(const idVertex seed, const idEdge neigEdge,
                                                      const idSuperArc         curArc,
                                                      const Propagation* const localProp)
      {
         idVertex v0;
         idVertex v1;
         mesh_.getEdgeVertex(neigEdge, 0, v0);
         mesh_.getEdgeVertex(neigEdge, 1, v1);

         const idVertex other = (v0 == seed) ? v1 : v0;

         if (localProp->compare(seed, other)) {
            dynGraph(localProp).setSubtreeArc(neigEdge, curArc);
         } else {
            std::cout << "not updated DG " << seed << " for arc " << graph_.printArc(curArc) << std::endl;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateDynGraphCurArc(const idVertex seed, const idSuperArc curArc,
                                                      const Propagation* const localProp)
      {

         const idVertex nbEdgesNeigh = mesh_.getVertexEdgeNumber(seed);
         for(idVertex nid = 0; nid < nbEdgesNeigh; ++nid) {
            idEdge edgeId;
            mesh_.getVertexEdge(seed, nid, edgeId);

            idVertex v0;
            idVertex v1;
            mesh_.getEdgeVertex(edgeId, 0, v0);
            mesh_.getEdgeVertex(edgeId, 1, v1);

            const idVertex other = (v0 == seed) ? v1 : v0;

            if (localProp->compare(seed, other)) {
               dynGraph(localProp).setSubtreeArc(edgeId, curArc);
            }
         }
      }

      template <typename ScalarType>
      idNode FTRGraph<ScalarType>::updateReebGraph(const idSuperArc         currentArc,
                                                   const Propagation* const localPropagation)
      {
         const idVertex upVert = localPropagation->getCurVertex(); // keep before merge
         const idNode upNode = graph_.makeNode(upVert);
         graph_.closeArc(currentArc, upNode);

         DEBUG_1(<< "close arc " << graph_.printArc(currentArc) << std::endl);

         // To investigate
         if (graph_.getArc(currentArc).getDownNodeId() == graph_.getArc(currentArc).getUpNodeId()) {
            graph_.getArc(currentArc).hide();
         }

         return upNode;
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::localGrowth(Propagation* const localPropagation)
      {
         const idVertex nbNeigh = mesh_.getVertexNeighborNumber(localPropagation->getCurVertex());
         for (idVertex n = 0; n < nbNeigh; ++n) {
            idVertex neighId;
            mesh_.getVertexNeighbor(localPropagation->getCurVertex(), n, neighId);
            if (localPropagation->compare(localPropagation->getCurVertex(), neighId)) {
               if (!toVisit_[neighId] || toVisit_[neighId] != localPropagation->getRpz()) {
                  localPropagation->addNewVertex(neighId);
                  toVisit_[neighId] = localPropagation->getRpz();
               }
            }
         }
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkLast(const idSuperArc           currentArc,
                                           const Propagation* const   localPropagation,
                                           const std::vector<idEdge>& lowerStarEdges)
      {
         const idVertex   curSaddle = localPropagation->getCurVertex();
         valence decr = 0;

         // TODO use LowerStar edge instead of crossing all you dumb foolish stupid moron
         const idVertex nbEdgesNeigh = mesh_.getVertexEdgeNumber(curSaddle);
         for(idVertex nid = 0; nid < nbEdgesNeigh; ++nid) {
            idEdge edgeId;
            mesh_.getVertexEdge(curSaddle, nid, edgeId);

           idVertex v0;
           idVertex v1;
           mesh_.getEdgeVertex(edgeId, 0, v0);
           mesh_.getEdgeVertex(edgeId, 1, v1);

           const idVertex other = (v0 == curSaddle) ? v1 : v0;

           if (localPropagation->compare(other, curSaddle)) {
              const idSuperArc edgeArc = dynGraph(localPropagation).getSubtreeArc(edgeId);
              if (edgeArc == currentArc) {
                 ++decr;
                 DEBUG_1(<< printEdge(edgeId, localPropagation) << " decrement " << static_cast<unsigned>(decr) << " " << curSaddle << std::endl);
              } else {
                 DEBUG_1(<< printEdge(edgeId, localPropagation) << " no decrement " << edgeArc << std::endl);
              }
           }
         }


         valence oldVal = 0;
         if (localPropagation->goUp()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            {
               oldVal = graph_.valDown_[curSaddle];
               graph_.valDown_[curSaddle] -= decr;
            }

         } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            {
               oldVal = graph_.valUp_[curSaddle];
               graph_.valUp_[curSaddle] -= decr;
            }
         }

         if (oldVal == -1) {
            // First task to touch this saddle, compute the valence
            idVertex totalVal = lowerStarEdges.size();
            valence  newVal   = 0;
            if (localPropagation->goUp()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
               {
                  newVal = graph_.valDown_[curSaddle];
                  graph_.valDown_[curSaddle] += (totalVal + 1);
               }
            } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
               {
                  newVal = graph_.valUp_[curSaddle];
                  graph_.valUp_[curSaddle] += (totalVal + 1);
               }
            }
            oldVal = decr + newVal + (totalVal + 1);
         }

         return oldVal == decr;
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkOppositeDGForArc(const idVertex saddle, const idVertex neigh,
                                                       const Propagation* const localProp)
      {
         const idNode checkNode  = graph_.getNodeId(saddle);
         const idEdge nbAdjEdges = mesh_.getVertexEdgeNumber(saddle);

         for (idEdge e = 0; e < nbAdjEdges; ++e) {
            idEdge edgeId;
            mesh_.getVertexEdge(saddle, e, edgeId);
            idVertex edgeLowerVert, edgeUpperVert;
            std::tie(edgeLowerVert, edgeUpperVert) = mesh_.getOrderedEdge(edgeId, localProp->goUp());
            if (edgeUpperVert == neigh || edgeLowerVert == neigh) {
               // curedge is the one between the two vertices
               // WARNING: Check the opposite prop here, In parallel ??
               const idSuperArc rpzArc = dynGraph(!localProp->goUp()).getSubtreeArc(edgeId);
               if (rpzArc != nullSuperArc && graph_.getArc(rpzArc).getUpNodeId() == checkNode) {
                  const idNode downNode = graph_.getArc(rpzArc).getDownNodeId();
                  if (localProp->compare(saddle, graph_.getNode(downNode).getVertexIdentifier())) {
                     return true;
                  }
               }
               // break;
            }
         }

         return false;
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkSegmentationForArc(const idVertex           saddle,
                                                         const idVertex           regular,
                                                         const Propagation* const localProp)
      {
         auto comp = [localProp](const idVertex a, const idVertex b) {
            return localProp->compare(a, b);
         };
         return graph_.hasArcEndingHere(saddle, regular, comp);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::mergeAtSaddle(const idNode saddleId, Propagation* localProp)
      {
         auto comp = [localProp](const idVertex a, const idVertex b) {
            return localProp->compare(a, b);
         };

         graph_.mergeAtSaddle(saddleId, localProp);
         const idVertex saddleVert = graph_.getNode(saddleId).getVertexIdentifier();
         localProp->removeBelow(saddleVert, comp);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::splitAtSaddle(const Propagation* const localProp)
      {
         const idVertex curVert = localProp->getCurVertex();
         const idNode   curNode = graph_.getNodeId(curVert);

         std::vector<std::tuple<idSuperArc, Propagation*>> bfsResults;
         bfsResults.reserve(4);

         const idCell nbTriNeigh = mesh_.getVertexTriangleNumber(curVert);
         for(idCell t = 0; t < nbTriNeigh; ++t) {
            idCell neighTriangle;
            mesh_.getVertexTriangle(curVert, t, neighTriangle);
            const orderedTriangle oNeighTriangle = mesh_.getOrderedTriangle(neighTriangle, localProp->goUp());
            // only if curVert is not the highest point
            if (bfsCells_[neighTriangle] != curVert &&
                getVertPosInTriangle(oNeighTriangle, localProp) != vertPosInTriangle::End) {

               const idVertex endTri          = getEndVertexInTriangle(oNeighTriangle, localProp);
               const bool     alreadyAttached = checkOppositeDGForArc(curVert, endTri, localProp);
               if (alreadyAttached)
                  continue;

               DEBUG_1(<< "split " << curVert << " at " << printTriangle(neighTriangle, localProp) << std::endl);

               // BFS to add vertices in the current propagation for each seed
               // and its corresponfing arc
               Propagation* curProp = newPropagation(curVert, localProp->goUp()/*, localProp->getRpz()*/);
               const idSuperArc newArc = graph_.openArc(curNode, curProp);

               // give this arc to the DG component
               const idEdge crossedEdge = getEdgeFromOTri(oNeighTriangle, curVert, endTri);
               updateDynGraphCurArc(curVert, crossedEdge, newArc, curProp);

               // fill curProp using a BFS on the current seed
               bfsPropagation(curVert, neighTriangle, curProp, newArc);
               // remove the already processed first vertex
               curProp->nextVertex();
               curProp->removeDuplicates(curVert);
               if(!curProp->empty()) {
                  bfsResults.emplace_back(std::make_tuple(newArc, curProp));
               } else {
                  graph_.getArc(newArc).hide();
               }
            }
         }

         for (auto& bfsRes : bfsResults) {
            const auto arc  = std::get<0>(bfsRes);
            const auto prop = std::get<1>(bfsRes);
            graph_.visit(curVert, arc);
            DEBUG_1(<< "visit s: " << curVert << " with " << arc << std::endl);
            // why is the firstprivate required here ?
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(curVert, prop, arc) priority(2)
#endif
            growthFromSeed(curVert, prop, arc);
         }
      }

      /// Tools

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::newPropagation(const idVertex leaf, const bool fromMin)
      {
         Propagation* localPropagation;
         if (fromMin) {
            auto compare_max_fun = [&](idVertex a, idVertex b) { return scalars_->isHigher(a, b); };
            localPropagation = new Propagation(leaf, compare_max_fun, fromMin);
         } else {
            auto compare_min_fun = [&](idVertex a, idVertex b) { return scalars_->isLower(a, b); };
            localPropagation = new Propagation(leaf, compare_min_fun, fromMin);
         }
         const auto   propId   = propagations_.getNext();
         propagations_[propId] = localPropagation;
         return localPropagation;
      }

      template <typename ScalarType>
      idVertex FTRGraph<ScalarType>::getWeight(const orderedEdge& e0, const orderedEdge& e1,
                                               const Propagation* const localPropagation)
      {
         const idVertex end0 = std::get<1>(e0);
         const idVertex end1 = std::get<1>(e1);

         if (localPropagation->compare(end0, end1)) {
            if (localPropagation->goDown()) {
               return -scalars_->getMirror(end0);
            }
            return scalars_->getMirror(end0);
         }

         if (localPropagation->goDown()) {
            return -scalars_->getMirror(end1);
         }
         return scalars_->getMirror(end1);
      }

      template <typename ScalarType>
      vertPosInTriangle FTRGraph<ScalarType>::getVertPosInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localPropagation) const
      {
         orderedEdge firstEdge = mesh_.getOrderedEdge(std::get<0>(oTriangle), localPropagation->goUp());
         if (std::get<0>(firstEdge) == localPropagation->getCurVertex()) {
            return vertPosInTriangle::Start;
         } else if (std::get<1>(firstEdge) == localPropagation->getCurVertex()) {
            return vertPosInTriangle::Middle;
         } else {
            return vertPosInTriangle::End;
         }
      }

      template <typename ScalarType>
      idVertex FTRGraph<ScalarType>::getEndVertexInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localPropagation) const
      {
         const orderedEdge& higherEdge = mesh_.getOrderedEdge(std::get<1>(oTriangle), localPropagation->goUp());
         return std::get<1>(higherEdge);
      }

      template <typename ScalarType>
      idEdge FTRGraph<ScalarType>::getEdgeFromOTri(const orderedTriangle oTri, const idVertex v0, const idVertex v1)
      {
         idVertex edge0Vert, edge1Vert;

         mesh_.getEdgeVertex(std::get<0>(oTri), 0, edge0Vert);
         mesh_.getEdgeVertex(std::get<0>(oTri), 1, edge1Vert);
         if ((edge0Vert == v0 && edge1Vert == v1) || (edge0Vert == v1 && edge1Vert == v0)) {
            return std::get<0>(oTri);
         }

         mesh_.getEdgeVertex(std::get<1>(oTri), 0, edge0Vert);
         mesh_.getEdgeVertex(std::get<1>(oTri), 1, edge1Vert);
         if ((edge0Vert == v0 && edge1Vert == v1) || (edge0Vert == v1 && edge1Vert == v0)) {
            return std::get<1>(oTri);
         }

#ifndef TTK_ENABLE_KAMIKAZE
         mesh_.getEdgeVertex(std::get<2>(oTri), 0, edge0Vert);
         mesh_.getEdgeVertex(std::get<2>(oTri), 1, edge1Vert);
         if ((edge0Vert == v0 && edge1Vert == v1) || (edge0Vert == v1 && edge1Vert == v0)) {
            return std::get<2>(oTri);
         }

         std::cout << "[FTR]: edge not found in triangle " << v0 << " " << v1 << std::endl;
         return nullEdge;
#else
         return std::get<2>(oTri);
#endif
      }
   }
}

#endif /* end of include guard: FTRGRAPHPRIVATE_TEMPLATE_H */
