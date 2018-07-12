#ifndef FTRGRAPHPRIVATE_TEMPLATE_H
#define FTRGRAPHPRIVATE_TEMPLATE_H

#include "FTRGraph.h"
#include "Tasks.h"

// Skeleton + propagation
#ifndef NDEBUG
#define DEBUG_1(msg) std::cout msg
// #define DEBUG_1(msg)
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

#ifdef TTK_ENABLE_KAMIKAZE
#define OPTIONAL_PRIORITY(value) priority(value)
#else
#define OPTIONAL_PRIORITY(value)
#endif

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localProp,
                                                idSuperArc currentArc)
      {
         DEBUG_1(<< "Start " << seed << " go up " << localProp->goUp() << std::endl);
#ifndef NDEBUG
         DEBUG_1(<< localProp->getId() << " " << localProp->print() << std::endl);
#endif

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle     = false, isSplitSaddle = false;

         // containers
         std::vector<idEdge>                  lowerStarEdges, upperStarEdges;
         std::vector<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localProp->empty()) {
            localProp->nextVertex();
            const idVertex curVert = localProp->getCurVertex();

            // Check history for visit (quick test)
            if (!graph_.isNode(curVert) && propagations_.hasVisited(curVert, localProp)) {
               DEBUG_1(<< "already seen " << curVert << " " << localProp->getId() << std::endl);
               continue;
            }

            // Caution: crossing tasks can leads to revisit legally
            // if (graph_.isArc(curVert) && checkSegmentationForArc(seed, curVert, localProp)) {
            //    DEBUG_1(<< "Dismiss current : " << graph_.printArc(currentArc) << std::endl);
            //    graph_.getArc(currentArc).hide();
            //    return;
            // }

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);

            lowerComp = lowerComps(lowerStarEdges, localProp);
            if(lowerComp.size() > 1){
               // This task will stop here
               isJoinSaddle = true;
            } else if (lowerComp.size() /* == 1 */) {
               // recover the arc of this non join saddle. A same propagation
               // can deal with several arc after a split saddle where no BFS
               // where made
               currentArc = lowerComp[0]->getCorArc();
               DEBUG_1(<< "arc: " << lowerComp[0]->getCorArc() << " v " << curVert << std::endl);

               graph_.visit(curVert, currentArc);
               DEBUG_1(<< "visit n: " << curVert << std::endl);
            }

            if (isJoinSaddle) {
               isJoinSadlleLast = checkLast(currentArc, localProp, lowerStarEdges);
               DEBUG_1(<< ": is join " << isJoinSadlleLast << std::endl);
               // We stop here in the join case as we will update preimage
               // only after have processed arcs on the last join
               break;
            }

            updatePreimage(localProp, currentArc);

            upperComp = upperComps(upperStarEdges, localProp);
            if (upperComp.size() > 1) {
               if (!isJoinSaddle) {
                  graph_.visit(curVert, currentArc);
                  DEBUG_1(<< "visit n: " << curVert << std::endl);
               }
               DEBUG_1(<< ": is split" << std::endl);
               DEBUG_1(<< localProp->print() << std::endl);
               isSplitSaddle = true;
            }

            // add upper star for futur visit
            bool seenAbove = localGrowth(localProp);
            if (!seenAbove) {
               // We have reached an opposite leaf
               const idNode upNode = graph_.makeNode(curVert);
               graph_.closeArc(currentArc, upNode);
               DEBUG_1(<< "close arc max " << graph_.printArc(currentArc) << std::endl);
            }
         }

         // get the corresponging critical point on which
         // the propagation has stopped (join, split, max)
         const idVertex upVert = localProp->getCurVertex();
         const idNode   upNode = graph_.makeNode(upVert);

         // At saddle: join or split or both

         idSuperArc joinNewArc;
         // arriving at a join
         if (isJoinSadlleLast) {
            localGrowth(localProp);
            mergeAtSaddle(upNode, localProp, lowerComp);
            const idNode     downNode = graph_.getNodeId(upVert);
            joinNewArc   = graph_.openArc(downNode, localProp);
            updatePreimage(localProp, joinNewArc);
            graph_.visit(upVert, joinNewArc);
            upperComp = upperComps(upperStarEdges, localProp);
            if (upperComp.size() > 1) {
               DEBUG_1(<< ": is split" << std::endl);
               DEBUG_1(<< localProp->print() << std::endl);
               isSplitSaddle = true;
               // will be replaced be new arcs of the split
               graph_.getArc(joinNewArc).hide();
            }
         }

         // starting from the saddle
         if (isSplitSaddle) {
            if (!isJoinSaddle) {
               // only one arc coming here
               graph_.closeArc(currentArc, upNode);
               DEBUG_1(<< "close arc split " << graph_.printArc(currentArc) << std::endl);
            }
            const bool withBFS = false;
            if (withBFS) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task OPTIONAL_PRIORITY(PriorityLevel::Min)
#endif
               splitAtSaddleBFS(localProp);
            } else {
               splitAtSaddle(localProp, upperComp);
            }
         } else if (isJoinSadlleLast) {
            // recursive call
#ifdef TTK_ENABLE_OPENMP
#pragma omp task OPTIONAL_PRIORITY(PriorityLevel::Average)
#endif
            growthFromSeed(upVert, localProp, joinNewArc);
         }
      }

      template <typename ScalarType>
      std::pair<std::vector<idEdge>, std::vector<idEdge>> FTRGraph<ScalarType>::visitStar(
          const Propagation* const localProp) const
      {
         // TODO re-use the same vectors per thread
         std::vector<idEdge> lowerStar, upperStar;

         const idEdge nbAdjEdges = mesh_.getVertexEdgeNumber(localProp->getCurVertex());
         lowerStar.reserve(nbAdjEdges);
         upperStar.reserve(nbAdjEdges);

         for (idEdge e = 0; e < nbAdjEdges; ++e) {
            idEdge edgeId;
            mesh_.getVertexEdge(localProp->getCurVertex(), e, edgeId);
            idVertex edgeLowerVert, edgeUpperVert;
            std::tie(edgeLowerVert, edgeUpperVert) = mesh_.getOrderedEdge(edgeId, localProp->goUp());
            if (edgeLowerVert == localProp->getCurVertex()) {
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
      void FTRGraph<ScalarType>::updatePreimage(const Propagation* const localProp,
                                                const idSuperArc         curArc)
      {
         const idCell nbAdjTriangles =
             mesh_.getVertexTriangleNumber(localProp->getCurVertex());

         DEBUG_1(<< "update preimage " << localProp->getCurVertex() << std::endl);

         orderedTriangle oTriangle;

         for (idCell t = 0; t < nbAdjTriangles; ++t) {
            // Classify current cell
            idCell curTriangleid;
            mesh_.getVertexTriangle(localProp->getCurVertex(), t, curTriangleid);

            mesh_.getOrderedTriangle(curTriangleid, localProp->goUp(), oTriangle);
            vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localProp);

            // Update DynGraph
            // We can have an end pos on an unvisited triangle
            // in case of saddle points
            switch (curVertPos) {
               case vertPosInTriangle::Start:
                  updatePreimageStartCell(oTriangle, localProp, curArc);
                  break;
               case vertPosInTriangle::Middle:
                  updatePreimageMiddleCell(oTriangle, localProp, curArc);
                  break;
               case vertPosInTriangle::End:
                  updatePreimageEndCell(oTriangle, localProp, curArc);
                  break;
               default:
                  std::cout << "[FTR]: update preimage error, unknown vertPos type" << std::endl;
                  break;
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageStartCell(const orderedTriangle&   oTriangle,
                                                         const Propagation* const localProp,
                                                         const idSuperArc         curArc)
      {
         const orderedEdge e0 = mesh_.getOrderedEdge(std::get<0>(oTriangle), localProp->goUp());
         const orderedEdge e1 = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
         const idVertex    w  = getWeight(e0, e1, localProp);

         bool t;
         // if (!dynGraph(localProp).getNode(std::get<0>(oTriangle))->hasParent()) {
         //    t = dynGraph(localProp).insertEdge(std::get<0>(oTriangle), std::get<1>(oTriangle), w,
         //                                       curArc);
         // } else {
            t = dynGraph(localProp).insertEdge(std::get<1>(oTriangle), std::get<0>(oTriangle), w,
                                               curArc);
         // }

         if (t) {
            DEBUG_2(<< "start add edge: " << printEdge(std::get<0>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl);
         } else {
            DEBUG_2(<< "start no need to create edge: " << printEdge(std::get<0>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl);
            DEBUG_2(<< dynGraph(localProp).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageMiddleCell(const orderedTriangle&   oTriangle,
                                                          const Propagation* const localProp,
                                                          const idSuperArc         curArc)
      {
         // Check if exist ?
         // If not, the triangle will be visited again once a merge have occured.
         // So we do not add the edge now
         const int t = dynGraph(localProp).removeEdge(std::get<0>(oTriangle), std::get<1>(oTriangle));

         // keep history inside the dyngraph structure
         dynGraph(localProp).setSubtreeArc(std::get<0>(oTriangle), curArc);
         dynGraph(localProp).setSubtreeArc(std::get<1>(oTriangle), curArc);

         if (t) {
            DEBUG_2(<< "mid replace edge: " << printEdge(std::get<0>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl);
         }
         else {
            DEBUG_2(<< "mid no found edge: " << printEdge(std::get<0>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl);
            DEBUG_2(<< dynGraph(localProp).print() << std::endl);
         }

         const orderedEdge e1 = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
         const orderedEdge e2 = mesh_.getOrderedEdge(std::get<2>(oTriangle), localProp->goUp());
         const idVertex    w  = getWeight(e1, e2, localProp);

         bool u;
         // limit the number of evert
         // if (!dynGraph(localProp).getNode(std::get<2>(oTriangle))->hasParent()) {
         //    u = dynGraph(localProp).insertEdge(std::get<2>(oTriangle), std::get<1>(oTriangle), w,
         //                                       curArc);
         // } else {
            u = dynGraph(localProp).insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), w,
                                               curArc);
         // }

         if (u) {
            DEBUG_2(<< " new edge: " << printEdge(std::get<1>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl);
         } else {
            DEBUG_2(<< " mid no need to create edge: " << printEdge(std::get<1>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl);
            DEBUG_2(<< dynGraph(localProp).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageEndCell(const orderedTriangle&   oTriangle,
                                                       const Propagation* const localProp,
                                                       const idSuperArc         curArc)
      {
         const int t = dynGraph(localProp).removeEdge(std::get<1>(oTriangle), std::get<2>(oTriangle));

         dynGraph(localProp).setSubtreeArc(std::get<1>(oTriangle), curArc);
         dynGraph(localProp).setSubtreeArc(std::get<2>(oTriangle), curArc);

         if (t) {
            DEBUG_2(<< "end remove edge: " << printEdge(std::get<1>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl);
         } else {
            DEBUG_2(<< "end not found edge: " << printEdge(std::get<1>(oTriangle), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl);
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
      bool FTRGraph<ScalarType>::localGrowth(Propagation* const localProp)
      {
         bool seenAbove = false;
         const idVertex curVert = localProp->getCurVertex();
         const idVertex nbNeigh = mesh_.getVertexNeighborNumber(localProp->getCurVertex());
         for (idVertex n = 0; n < nbNeigh; ++n) {
            idVertex neighId;
            mesh_.getVertexNeighbor(curVert, n, neighId);
            if (localProp->compare(curVert, neighId)) {
               seenAbove = true;
               if (!propagations_.willVisit(neighId, localProp)) {
                  localProp->addNewVertex(neighId);
                  propagations_.toVisit(neighId, localProp);
                  DEBUG_1(<< " + " << neighId << std::endl);
               }
            }
         }
         return seenAbove;
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkLast(const idSuperArc           currentArc,
                                           Propagation* const         localProp,
                                           const std::vector<idEdge>& lowerStarEdges)
      {
         const idVertex curSaddle = localProp->getCurVertex();
         UnionFind*     curId     = localProp->getId();
         valence        decr      = 0;

          DEBUG_1(<< "Check last on " << curSaddle << " id " << curId << std::endl);

         // NOTE:
         // Using propagation id allows to decrement by the number of time this propagation
         // has reached the saddle, even if the propagation take care of several of these arcs
         // (after a Hole-split).
         for(idEdge edgeId : lowerStarEdges) {
           const idSuperArc    edgeArc = dynGraph(localProp).getSubtreeArc(edgeId);
           if (edgeArc == nullSuperArc) // ignore unseen
              continue;
           UnionFind* tmpId  = graph_.getArc(edgeArc).getPropagation()->getId();
           if (tmpId == curId) {
              ++decr;
              DEBUG_1(<< printEdge(edgeId, localProp) << " decrement " << static_cast<unsigned>(decr) << " " << curSaddle << std::endl);
           } else {
              DEBUG_1(<< printEdge(edgeId, localProp) << " no decrement " << edgeArc << std::endl);
           }
         }

         valence oldVal = 0;
         if (localProp->goUp()) {
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
            if (localProp->goUp()) {
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
                                                       Propagation* const localProp)
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
      void FTRGraph<ScalarType>::mergeAtSaddle(
          const idNode saddleId, Propagation* localProp,
          const std::vector<DynGraphNode<idVertex>*>& lowerComp)
      {

#ifndef TTK_ENABLE_KAMIKAZE
         if (lowerComp.size() < 2) {
            std::cerr << "[FTR]: merge at saddle with only one lower CC" << std::endl;
         }
#endif

         for(auto* dgNode : lowerComp) {
            // read in the history
            const idSuperArc endingArc = dgNode->getCorArc();
            graph_.closeArc(endingArc, saddleId);
            localProp->merge(*graph_.getArc(endingArc).getPropagation());
            // graph_.getArc(endingArc).setPropagation(localProp);
            // NOTE: when merging update the ID of the old propagation to the new
            // its a bit like union find but better.
            // for other component after a split, update the ID for chaekLast detection
            graph_.getArc(endingArc).getPropagation()->mergeId(localProp->getId());
            DEBUG_1(<< "Close arc join " << graph_.printArc(endingArc) << std::endl);
         }

         // const idVertex saddleVert = graph_.getNode(saddleId).getVertexIdentifier();
         // localProp->removeBelow(saddleVert, comp);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::splitAtSaddleBFS(Propagation* const localProp)
      {
         const idVertex curVert = localProp->getCurVertex();
         const idNode   curNode = graph_.getNodeId(curVert);

         std::vector<std::tuple<idSuperArc, Propagation*>> bfsResults;
         bfsResults.reserve(4);
         const idCell nbTriNeigh = mesh_.getVertexTriangleNumber(curVert);
         for(idCell t = 0; t < nbTriNeigh; ++t) {
            idCell neighTriangle;
            mesh_.getVertexTriangle(curVert, t, neighTriangle);
            const orderedTriangle oNeighTriangle =
                mesh_.getOrderedTriangle(neighTriangle, localProp->goUp());
            // only if curVert is not the highest point
            if (bfsCells_[neighTriangle] != curVert &&
                getVertPosInTriangle(oNeighTriangle, localProp) != vertPosInTriangle::End) {
               const idVertex endTri          = getEndVertexInTriangle(oNeighTriangle, localProp);
               const bool     alreadyAttached = checkOppositeDGForArc(curVert, endTri, localProp);
               if (alreadyAttached)
                  continue;

               DEBUG_1(<< "split " << curVert << " at " << printTriangle(neighTriangle, localProp)
                       << std::endl);

               // BFS to add vertices in the current propagation for each seed
               // and its corresponfing arc
               Propagation*     newProp = newPropagation(curVert, localProp->goUp());
               const idSuperArc newArc  = graph_.openArc(curNode, newProp);

               // give this arc to the DG component
               const idEdge crossedEdge = getEdgeFromOTri(oNeighTriangle, curVert, endTri);
               updateDynGraphCurArc(curVert, crossedEdge, newArc, newProp);

               // fill newProp using a BFS on the current seed
               bfsPropagation(curVert, neighTriangle, newProp, newArc);
               // remove the already processed first vertex
               newProp->nextVertex();
               newProp->removeDuplicates(curVert);
               if (!newProp->empty()) {
                  bfsResults.emplace_back(std::make_tuple(newArc, newProp));
               } else {
                  graph_.getArc(newArc).hide();
               }
            }
         }

         // one growth per connected components
         for (auto& bfsRes : bfsResults) {
            const auto arc  = std::get<0>(bfsRes);
            const auto prop = std::get<1>(bfsRes);
            graph_.visit(curVert, arc);
            DEBUG_1(<< "visit s: " << curVert << " with " << arc << std::endl);
            // why is the firstprivate required here ?
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(curVert, prop, arc) OPTIONAL_PRIORITY(PriorityLevel::Low)
#endif
            growthFromSeed(curVert, prop, arc);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::splitAtSaddle(
          Propagation* const localProp, const std::vector<DynGraphNode<idVertex>*>& upperComp)
      {
         const idVertex curVert = localProp->getCurVertex();
         const idNode   curNode = graph_.getNodeId(curVert);

         for (auto* dgNode : upperComp) {
            const idSuperArc newArc = graph_.openArc(curNode, localProp);
            dgNode->setRootArc(newArc);

            DEBUG_1(<< "set root arc " << newArc << " at ");
            DEBUG_1(<< printEdge(dynGraphs_.up.getNodeId(dgNode), localProp));
            DEBUG_1(<< " id " << localProp->getId() << std::endl);
         }

#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(curVert, localProp) OPTIONAL_PRIORITY(PriorityLevel::Average)
#endif
            growthFromSeed(curVert, localProp);
      }

      /// Tools

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::newPropagation(const idVertex leaf, const bool fromMin)
      {
         VertCompFN comp;
         if (fromMin)
            comp = [&](idVertex a, idVertex b) { return scalars_->isHigher(a, b); };
         else
            comp = [&](idVertex a, idVertex b) { return scalars_->isLower(a, b); };
         return propagations_.newPropagation(leaf, comp, fromMin);
      }

      template <typename ScalarType>
      idVertex FTRGraph<ScalarType>::getWeight(const orderedEdge& e0, const orderedEdge& e1,
                                               const Propagation* const localProp)
      {
         const idVertex end0 = std::get<1>(e0);
         const idVertex end1 = std::get<1>(e1);

         if (localProp->compare(end0, end1)) {
            if (localProp->goDown()) {
               return -scalars_->getMirror(end0);
            }
            return scalars_->getMirror(end0);
         }

         if (localProp->goDown()) {
            return -scalars_->getMirror(end1);
         }
         return scalars_->getMirror(end1);
      }

      template <typename ScalarType>
      vertPosInTriangle FTRGraph<ScalarType>::getVertPosInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localProp) const
      {
         orderedEdge firstEdge = mesh_.getOrderedEdge(std::get<0>(oTriangle), localProp->goUp());
         if (std::get<0>(firstEdge) == localProp->getCurVertex()) {
            return vertPosInTriangle::Start;
         } else if (std::get<1>(firstEdge) == localProp->getCurVertex()) {
            return vertPosInTriangle::Middle;
         } else {
            return vertPosInTriangle::End;
         }
      }

      template <typename ScalarType>
      idVertex FTRGraph<ScalarType>::getEndVertexInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localProp) const
      {
         const orderedEdge& higherEdge = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
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
