#ifndef FTRGRAPHPRIVATE_TEMPLATE_H
#define FTRGRAPHPRIVATE_TEMPLATE_H

// local includes
#include "FTRGraph.h"
#include "Tasks.h"

// c++ incldues
#include <unordered_map>

// Skeleton + propagation
#ifndef NDEBUG
#define DEBUG_1(msg) std::cout msg
// #define DEBUG_1(msg)
#else
#define DEBUG_1(msg)
// #define DEBUG_1(msg) std::cout msg
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
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localProp,
                                                idSuperArc currentArc)
      {
         DEBUG_1(<< seed << " go up " << localProp->goUp() << " arc " << currentArc << " nb "
                 << localProp->getNbArcs() << std::endl);

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle = false, isSplitSaddle = false;

         // containers
         // vitsit
         std::vector<idEdge>                  lowerStarEdges, upperStarEdges;
         std::vector<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localProp->empty()) {
            localProp->nextVertex();
            const idVertex curVert = localProp->getCurVertex();
            idSuperArc mergeIn = nullSuperArc;

            DEBUG_1(<< curVert << std::endl);

            // Check history for visit (quick test)
            if (propagations_.hasVisited(curVert, localProp)) {
               DEBUG_1(<< curVert << "already seen " << localProp->getId() << std::endl);
               continue;
            }

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);

#ifndef TTK_DISABLE_FTR_LAZY
            if (valences_.lower[curVert] < 2 && valences_.upper[curVert] < 2) {

               // not a local min (for local min, currentArc is already set)
               if(lowerStarEdges.size()) {
                  // not a min nor a saddle: 1 CC below
                  currentArc = dynGraph(localProp).getSubtreeArc(lowerStarEdges[0]);
                  if (currentArc == nullSuperArc) {
                     DEBUG_1(<< "- " << curVert << std::endl);
                     // merging component, propagation has been discarded
                     continue;
                  }
                  if(valences_.upper[curVert] && valences_.lower[curVert]){
                     // not saddle neither extrema
                     graph_.getArc(currentArc).visit(curVert);
                  }
               }
               // ensure we will always recover this arc from the upper neighbors
               for (const idEdge dgNode : upperStarEdges) {
                  dynGraph(localProp).setCorArc(dgNode, currentArc);
               }

               mergeIn = visit(localProp, currentArc);

               lazyUpdatePreimage(localProp, currentArc);

            } else {

               // locally apply the lazy one the current growing arc
               for (const idEdge e : lowerStarEdges) {
                  const idSuperArc a = dynGraph(localProp).getNode(e)->findRootArc();
                  if (a != nullSuperArc && graph_.getArc(a).getPropagation()->getId() == localProp->getId()) {
                     lazyApply(localProp, a);
                  }
               }
# else
            {
#endif
               lowerComp = lowerComps(lowerStarEdges, localProp);
               if (lowerComp.size() > 1) {
                  isJoinSaddle = true;
                  isJoinSadlleLast = checkLast(localProp, lowerStarEdges);
                  break;
               } else {
                  if (lowerComp.size()) {
                     currentArc = lowerComp[0]->getCorArc();
                     if (currentArc == nullSuperArc) {
                        // current arc is merging, the lower star has been discontinued
                        continue;
                     }
                  }
                  mergeIn = visit(localProp, currentArc);
               }
               updatePreimage(localProp, currentArc);
               upperComp = upperComps(upperStarEdges, localProp);
               if (upperComp.size() > 1) {
                  isSplitSaddle = true;
               }

               if (!isJoinSaddle && !isSplitSaddle && upperComp.size()){
                  // this arc is not empty (not saddle not max)
                  graph_.getArc(currentArc).visit(curVert);
               }
            }

            // do not propagate on merging arc.
            if (mergeIn != nullSuperArc) {
               if (graph_.getArc(currentArc).isVisible()) {
                  localProp->lessArc();
               }
               graph_.getArc(currentArc).merge(mergeIn);
               DEBUG_1(<< curVert << " arc merging " << graph_.printArc(currentArc) << " in "
                     << graph_.printArc(mergeIn) << std::endl);
            }

            // stop on leaves
            if (!upperStarEdges.size()) {
               // We have reached a local extrema
               const idNode upNode = graph_.makeNode(curVert);
               graph_.closeArc(currentArc, upNode);
               if (graph_.getArc(currentArc).isVisible()) {
                  // do not decrease on merged arcs
                  localProp->lessArc();
               }
               DEBUG_1(<< curVert << " arc max " << graph_.printArc(currentArc) << std::endl);
            }


            // add upper star for futur visit
            localGrowth(localProp, upperStarEdges);
         } // end propagation while

         // get the corresponging critical point on which
         // the propagation has stopped (join, split, max)
         const idVertex upVert = localProp->getCurVertex();
         // reached node id and wether it has been created by this task or already existed
         idNode upNode;

         if (isJoinSaddle) {
            DEBUG_1(<< ": is join " << isJoinSadlleLast << std::endl);
         }
         if (isSplitSaddle) {
            DEBUG_1(<< ": is split " << std::endl);
         }

         // At saddle: join or split or both
         idSuperArc joinNewArc;
         // arriving at a join
         if (isJoinSadlleLast) {
            // ensure we have the good values here, even if other tasks were doing stuff
            {
               // here to solve a 1 over thousands execution bug in parallel
               std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);
               lowerComp = lowerComps(lowerStarEdges, localProp);
            }
            localGrowth(localProp, upperStarEdges);
            upNode = graph_.makeNode(upVert);
            idSuperArc visibleMerged = mergeAtSaddle(upNode, localProp, lowerComp);
            localProp->lessArc(visibleMerged - 1);

            const idNode downNode = graph_.getNodeId(upVert);
            joinNewArc            = graph_.openArc(downNode, localProp);
            visit(localProp, joinNewArc);
            updatePreimage(localProp, joinNewArc);
            upperComp = upperComps(upperStarEdges, localProp);
            if (upperComp.size() > 1) {
               isSplitSaddle = true;
               DEBUG_1(<< ": is join & split : " << localProp->print() << std::endl);
               // will be replaced be new arcs of the split
               graph_.getArc(joinNewArc).hide();
            }
         }


         // this propagation can stop
         if (localProp->getNbArcs() == 0) {
            DEBUG_1(<< "proapgation stop here, no active arcs" << std::endl);
            return;
         }

         // starting from the saddle
         if (isSplitSaddle && (!isJoinSaddle || isJoinSadlleLast)) {
            if (!isJoinSaddle) {
               // only one arc coming here
               upNode = graph_.makeNode(upVert);
               graph_.closeArc(currentArc, upNode);
               DEBUG_1(<< "close arc split " << graph_.printArc(currentArc) << std::endl);
               // current arc stop here
               localProp->lessArc(graph_.getArc(currentArc).isVisible());
            }

#ifdef TTK_ENABLE_FTR_BFS
            Propagation* remainProp = splitAtSaddleBFS(localProp);
            growthFromSeed(upVert, remainProp);
#else
            splitAtSaddle(localProp, upperComp);
            localProp->moreArc(upperComp.size());
            growthFromSeed(upVert, localProp);
            // Can't stop propagation here cause some other arcs thant the current
            // one may nees to growth
#endif
         } else if (isJoinSadlleLast) {
            growthFromSeed(upVert, localProp, joinNewArc);
         }
#ifdef TTK_ENABLE_FTR_STATS
         else {
            // This propagation is dying here
            idVertex curProp;
            float    curTime;
#pragma omp atomic capture
            {
               curProp = nbProp_;
               --nbProp_;
            }
#pragma omp critical(stats)
            {
               curTime = sweepStart_.getElapsedTime();
            }
            propTimes_[curProp - 1] = curTime;
         }
#endif

      }

      template<typename ScalarType>
      void FTRGraph<ScalarType>::growthSequential(const idVertex begin, const idVertex stop)
      {
         std::vector<idEdge>                  lowerStarEdges, upperStarEdges;
         std::vector<DynGraphNode<idVertex>*> lowerComp, upperComp;

         const bool     fromMin         = begin < stop;
         Propagation*   localProp = newPropagation(scalars_->getSortedVert(begin), fromMin);
         const idVertex incr            = fromMin ? 1 : -1;
         for (idVertex idv = begin; idv < stop; idv = idv + incr) {
            const idVertex curVert = scalars_->getSortedVert(idv);
            localProp->setCurvert(curVert);  // never use Fibo Heap here

            idSuperArc currentArc = nullSuperArc;
            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);

            if (valences_.lower[curVert] == 0) { // min
               const idNode minNode = graph_.makeNode(curVert);
               currentArc           = graph_.openArc(minNode, localProp);
               DEBUG_1(<< curVert << " min arc " << currentArc << std::endl);
            }

#ifndef TTK_DISABLE_FTR_LAZY
            if (valences_.lower[curVert] < 2 && valences_.upper[curVert] < 2) {

               // simple reeb regular, lazyness
               if(lowerStarEdges.size()) {
                  // not a min nor a saddle: 1 CC below (need findSubtree)
                  currentArc = dynGraph(localProp).getSubtreeArc(lowerStarEdges[0]);
                  if(valences_.upper[curVert] && valences_.lower[curVert]){
                     // not saddle neither extrema
                     graph_.getArc(currentArc).visit(curVert);
                  }
               }
               if (upperStarEdges.size()) {
                  for (const idEdge dgNode : upperStarEdges) {
                     dynGraph(localProp).setCorArc(dgNode, currentArc);
                  }
               } else { // max
                  const idNode maxNode = graph_.makeNode(curVert);
                  graph_.closeArc(currentArc, maxNode);
                  DEBUG_1(<< curVert << "max arc close " << graph_.printArc(currentArc) << std::endl);
               }

               visit(localProp, currentArc);

               lazyUpdatePreimage(localProp, currentArc);
            } else {

               // locally aply the lazy one the current growing arc
               for (const idEdge e : lowerStarEdges) {
                  const idSuperArc a = dynGraph(localProp).getNode(e)->findRootArc();
                  if (!lazy_.isEmpty(a)) {
                     // process lazy
                     // sort both list
                     // for each:
                     // if del < add: delete in real RG
                     // if add < del: add in real RG
                     // else: drop both, computation avoided
                     lazyApply(localProp, a);
                  }
               }
# else
            {
#endif
               bool isJoin = false;
               lowerComp = lowerComps(lowerStarEdges, localProp);
               if (lowerComp.size() == 1) {  // regular
                  currentArc = lowerComp[0]->getCorArc();
               } else if (lowerComp.size() > 1) {  // join saddle
                  const idNode sadNode = graph_.makeNode(curVert);
                  currentArc           = graph_.openArc(sadNode, localProp);
                  DEBUG_1(<< curVert << " join arc " << currentArc << std::endl);
                  mergeAtSaddle(sadNode, lowerComp);
                  isJoin = true;
               }

               DEBUG_1(<< curVert << " arc " << currentArc << std::endl);
               graph_.visit(curVert, currentArc);
               propagations_.visit(curVert, localProp);
               updatePreimage(localProp, currentArc);

               upperComp = upperComps(upperStarEdges, localProp);
               if (!upperComp.size()) {
                  const idNode maxNode = graph_.makeNode(curVert);
                  graph_.closeArc(currentArc, maxNode);
                  DEBUG_1(<< curVert << "max arc close " << graph_.printArc(currentArc) << std::endl);
               } else if (upperComp.size() < 2) {
                  if (!isJoin) {
                     // this arc is not empty
                     graph_.getArc(currentArc).visit(curVert);
                  }
               } else {
                  if (isJoin) {
                     graph_.getArc(currentArc).hide();
                  } else {
                     const idNode splitNode = graph_.makeNode(curVert);
                     graph_.closeArc(currentArc, splitNode);
                     DEBUG_1(<< curVert << "split arc close " << graph_.printArc(currentArc) << std::endl);
                  }
                  splitAtSaddle(localProp, upperComp);
               }
            }
         } // end for each vertex
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
      std::pair<valence, valence> FTRGraph<ScalarType>::getLinkNbCC(const idVertex curVert,
                                                                    LocalForests&  localForests,
                                                                    VertCompFN     comp)
      {
         // traduce edge id in a local id for the forests
         std::unordered_map<idEdge, std::size_t> mapNeighDown, mapNeighUp;
         std::size_t nextId = 0;

         localForests.up.reset();
         localForests.down.reset();
         const idVertex oldUpCC   = localForests.up.getNbCC();
         const idVertex oldDownCC = localForests.down.getNbCC();
         const idCell   nbTri     = mesh_.getVertexTriangleNumber(curVert);
         for(idCell t = 0; t < nbTri; ++t) {
            idCell curTri;
            mesh_.getVertexTriangle(curVert, t, curTri);
            // edges found, max 2 per triangles.
            std::size_t downSide[2], upSide[2];
            unsigned char curDownSide = 0, curUpSide = 0;
            for(idVertex v = 0; v < 3; ++v) {
               idEdge triEdge;
               mesh_.getTriangleEdge(curTri, v, triEdge);
               idVertex v0, v1;
               mesh_.getEdgeVertex(triEdge, 0, v0);
               mesh_.getEdgeVertex(triEdge, 1, v1);

               // only consider the edges in the star
               if(v0 == curVert) {
                  if (comp(curVert, v1)) {
                     if (mapNeighUp.count(triEdge)) {
                        upSide[curUpSide++] = mapNeighUp[triEdge];
                     } else {
                        upSide[curUpSide++] = nextId;
                        mapNeighUp[triEdge] = nextId++;
                     }
                  } else {
                     if (mapNeighDown.count(triEdge)) {
                        downSide[curDownSide++] = mapNeighDown[triEdge];
                     } else {
                        downSide[curDownSide++] = nextId;
                        mapNeighDown[triEdge]   = nextId++;
                     }
                  }
               } else if (v1 == curVert) {
                  if (comp(curVert, v0)) {
                     if (mapNeighUp.count(triEdge)) {
                        upSide[curUpSide++] = mapNeighUp[triEdge];
                     } else {
                        upSide[curUpSide++] = nextId;
                        mapNeighUp[triEdge] = nextId++;
                     }
                  } else {
                     if (mapNeighDown.count(triEdge)) {
                        downSide[curDownSide++] = mapNeighDown[triEdge];
                     } else {
                        downSide[curDownSide++] = nextId;
                        mapNeighDown[triEdge]   = nextId++;
                     }
                  }
               }
            }

            // if both edges of the triangle were up or down, we add a link btwn them
            // This is how the number of component is reduced
            if(curDownSide == 2) {
               localForests.down.insertEdge(downSide[0], downSide[1], 0, nullSuperArc);
            } else if (curUpSide == 2) {
               localForests.up.insertEdge(upSide[0], upSide[1], 0, nullSuperArc);
            }
         }

         const valence down = mapNeighDown.size() - (oldDownCC - localForests.down.getNbCC());
         const valence up = mapNeighUp.size() - (oldUpCC - localForests.up.getNbCC());
         return {down,up};
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
         // this order for history
         t = dynGraph(localProp).insertEdge(std::get<1>(oTriangle), std::get<0>(oTriangle), w,
                                            curArc);

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
         dynGraph(localProp).setCorArc(std::get<0>(oTriangle), curArc);
         // dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);

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
         // this order for history
         u = dynGraph(localProp).insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), w,
                                            curArc);

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

         dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);
         // dynGraph(localProp).setCorArc(std::get<2>(oTriangle), curArc);

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

#ifndef TTK_DISABLE_FTR_LAZY
      template <typename ScalarType>
      void FTRGraph<ScalarType>::lazyUpdatePreimage(Propagation* const localProp,
                                                    const idSuperArc   curArc)
      {
         const idCell nbAdjTriangles =
             mesh_.getVertexTriangleNumber(localProp->getCurVertex());

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
                  updateLazyStart(oTriangle, localProp, curArc);
                  break;
               case vertPosInTriangle::Middle:
                  updateLazyMiddle(oTriangle, localProp, curArc);
                  break;
               case vertPosInTriangle::End:
                  updateLazyEnd(oTriangle, localProp, curArc);
                  break;
               default:
                  std::cout << "[FTR]: lazy update preimage error, unknown vertPos type" << std::endl;
                  break;
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyStart(const orderedTriangle& oTriangle,
                                                 Propagation* const     localProp,
                                                 const idSuperArc       curArc)
      {
         lazy_.addEmplace(std::get<0>(oTriangle), std::get<1>(oTriangle), curArc);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyMiddle(const orderedTriangle& oTriangle,
                                                  Propagation* const     localProp,
                                                  const idSuperArc       curArc)
      {
         lazy_.delEmplace(std::get<0>(oTriangle), std::get<1>(oTriangle), curArc);
         lazy_.addEmplace(std::get<1>(oTriangle), std::get<2>(oTriangle), curArc);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyEnd(const orderedTriangle& oTriangle,
                                               Propagation* const     localProp,
                                               const idSuperArc       curArc)
      {
         lazy_.delEmplace(std::get<1>(oTriangle), std::get<2>(oTriangle), curArc);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyAdd(const Propagation* const localProp,
                                               const linkEdge edge, const idSuperArc arc)
      {
         const orderedEdge e0      = mesh_.getOrderedEdge(std::get<0>(edge), localProp->goUp());
         const orderedEdge e1      = mesh_.getOrderedEdge(std::get<1>(edge), localProp->goUp());
         const idVertex    w       = getWeight(e0, e1, localProp);
         bool              t;
         t = dynGraph(localProp).insertEdge(std::get<1>(edge), std::get<0>(edge), w, arc);
         if (t) {
            DEBUG_2(<< "start add edge: " << printEdge(std::get<0>(edge), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(edge), localProp) << " w " << w << " arc " << arc << std::endl);
         } else {
            DEBUG_2(<< "start no need to create edge: "
                    << printEdge(std::get<0>(edge), localProp) << std::endl);
         // DEBUG_2(<< " :: " << printEdge(std::get<1>(edge), localProp) << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyDel(const Propagation* const localProp,
                                               const linkEdge edge, const idSuperArc arc)
      {
         const int t = dynGraph(localProp).removeEdge(std::get<0>(edge), std::get<1>(edge));

         // keep history inside the dyngraph structure
         // dynGraph(localProp).setSubtreeArc(std::get<0>(edge), std::get<1>(del));
         // dynGraph(localProp).setSubtreeArc(std::get<1>(edge), std::get<1>(del));

         if (t) {
            DEBUG_2(<< "mid del edge: " << printEdge(std::get<0>(edge), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(edge), localProp) << " : " << arc << std::endl);
         }
         else {
            DEBUG_2(<< "mid no found edge: " << printEdge(std::get<0>(edge), localProp));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(edge), localProp) << " : " << arc << std::endl);
         }
      }

      template<typename ScalarType>
      void FTRGraph<ScalarType>::lazyApply(Propagation* const localProp, const idSuperArc a)
      {
         auto comp = [localProp](const idVertex a, const idVertex b) {
            return localProp->compare(a, b);
         };

         DEBUG_1(<< "lazy apply " << localProp->getCurVertex() << " arc " << graph_.printArc(a)
                 << std::endl);

         auto add = lazy_.addGetNext(a);
         auto del = lazy_.delGetNext(a);
         while (add != nullLink || del != nullLink) {
            if(del == nullLink) {
               updateLazyAdd(localProp, add, a);
               add = lazy_.addGetNext(a);
            } else if (add == nullLink || mesh_.compareLinks(del, add, comp)) {
               updateLazyDel(localProp, del, a);
               del = lazy_.delGetNext(a);
            } else if (mesh_.compareLinks(add, del, comp)) {
               updateLazyAdd(localProp, add, a);
               add = lazy_.addGetNext(a);
            } else {
               // same arc in both list, should be added and removed so we jus ignore it
               // (add and del cant be null both of them at the same time)
               add = lazy_.addGetNext(a);
               del = lazy_.delGetNext(a);
            }
         }
      }

#endif

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateDynGraphCurArc(const idVertex seed, const idSuperArc curArc,
                                                      const Propagation* const localProp)
      {
         const idVertex nbEdgesNeigh = mesh_.getVertexEdgeNumber(seed);
         for(idVertex nid = 0; nid < nbEdgesNeigh; ++nid) {
            idEdge edgeId;
            mesh_.getVertexEdge(seed, nid, edgeId);

            updateDynGraphCurArc(seed, edgeId, curArc, localProp);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::localGrowth(Propagation* const localProp, const std::vector<idEdge>& upperEdges)
      {
         const idVertex curVert = localProp->getCurVertex();
         for (const idEdge e : upperEdges) {
            idVertex v0, v1;
            mesh_.getEdgeVertex(e, 0, v0);
            mesh_.getEdgeVertex(e, 1, v1);
            const idVertex other = (v0 == curVert) ? v1 : v0;
            if (localProp->compare(curVert, other)) {
               if (!propagations_.willVisit(other, localProp)) {
                  localProp->addNewVertex(other);
                  propagations_.toVisit(other, localProp);
                  // DEBUG_1(<< " + " << other << std::endl);
               }
            }
         }
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkLast(Propagation* const         localProp,
                                           const std::vector<idEdge>& lowerStarEdges)
      {
         const idVertex curSaddle = localProp->getCurVertex();
         AtomicUF*      curId     = localProp->getId();
         valence        decr      = 0;

          DEBUG_1(<< "Check last on " << curSaddle <<  std::endl);

          // NOTE:
          // Using propagation id allows to decrement by the number of time this propagation
          // has reached the saddle, even if the propagation take care of several of these arcs
          // (after a Hole-split).
          for (idEdge edgeId : lowerStarEdges) {
             const idSuperArc edgeArc = dynGraph(localProp).getSubtreeArc(edgeId);
             if (edgeArc == nullSuperArc) {
                DEBUG_1(<< printEdge(edgeId, localProp) << " ignore " << std::endl);
                continue;
             }
             AtomicUF* tmpId = graph_.getArc(edgeArc).getPropagation()->getId();
             if (tmpId == curId) {
                ++decr;
                DEBUG_1(<< printEdge(edgeId, localProp) << " decrement "
                        << static_cast<unsigned>(decr) << " " << curSaddle << std::endl);
             } else {
                DEBUG_1(<< printEdge(edgeId, localProp) << " no decrement " << edgeArc
                        << std::endl);
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
      idSuperArc FTRGraph<ScalarType>::mergeAtSaddle(
          const idNode saddleId, Propagation* localProp,
          const std::vector<DynGraphNode<idVertex>*>& lowerComp)
      {

#ifndef TTK_ENABLE_KAMIKAZE
         if (lowerComp.size() < 2) {
            std::cerr << "[FTR]: merge at saddle with only one lower CC" << std::endl;
         }
#endif

         // DEBUG_1(<< " merge in " << graph_.getNode(saddleId).getVertexIdentifier() << std::endl);

         idSuperArc visibleClosed = 0;
         for(auto* dgNode : lowerComp) {
            // read in the history (lower comp already contains roots)
            const idSuperArc endingArc = dgNode->getCorArc();
            graph_.closeArc(endingArc, saddleId);
            if (graph_.getArc(endingArc).isVisible()) {
               ++visibleClosed;
            }
            Propagation * arcProp =  graph_.getArc(endingArc).getPropagation();
            DEBUG_1(<< "merge " << graph_.printArc(endingArc) << std::endl);
            localProp->merge(*arcProp);
         }
         // DEBUG_1(<< " result " << localProp->print() << std::endl);
         return visibleClosed;
      }

      template <typename ScalarType>
      idSuperArc FTRGraph<ScalarType>::mergeAtSaddle(
          const idNode saddleId, const std::vector<DynGraphNode<idVertex>*>& lowerComp)
      {
         // version for the sequential arc growth, do not merge the propagations

#ifndef TTK_ENABLE_KAMIKAZE
         if (lowerComp.size() < 2) {
            std::cerr << "[FTR]: merge at saddle with only one lower CC" << std::endl;
         }
#endif

         idSuperArc visibleClosed = 0;
         for(auto* dgNode : lowerComp) {
            const idSuperArc endingArc = dgNode->getCorArc();
            graph_.closeArc(endingArc, saddleId);
            if (graph_.getArc(endingArc).isVisible()) {
               ++visibleClosed;
            }
            DEBUG_1(<< "close " << graph_.printArc(endingArc) << std::endl);
         }
         return visibleClosed;
      }

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::splitAtSaddleBFS(Propagation* const localProp)
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

         Propagation* remainProp = newPropagation(curVert, localProp->goUp());
         while(!localProp->empty()) {
            localProp->nextVertex();
            const idVertex curVert = localProp->getCurVertex();
            if (bfsVerts_[curVert] != curVert) {
               remainProp->addNewVertex(curVert);
            }
         }

         // one growth per connected components
         for (auto& bfsRes : bfsResults) {
            const auto arc  = std::get<0>(bfsRes);
            const auto prop = std::get<1>(bfsRes);
            graph_.visit(curVert, arc);
            propagations_.visit(curVert, localProp);
            DEBUG_1(<< "visit s: " << curVert << " with " << arc << std::endl);
            // why is the firstprivate required here ?
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(curVert, prop, arc) OPTIONAL_PRIORITY(PriorityLevel::Low)
#endif
            growthFromSeed(curVert, prop, arc);
         }

         return remainProp;
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
            visit(localProp, newArc);

            DEBUG_1(<< "set root arc " << newArc << " at ");
            DEBUG_1(<< printEdge(dynGraph(localProp).getNodeId(dgNode), localProp) << std::endl);
         }
      }

      template<typename ScalarType>
      idSuperArc FTRGraph<ScalarType>::visit(Propagation* const localProp, const idSuperArc curArc)
      {
         const idVertex curVert = localProp->getCurVertex();
         Visit          opposite;
#pragma omp critical
         {
            propagations_.visit(curVert, localProp);
            opposite = propagations_.visitOpposite(curVert, localProp);
         }
         if (!opposite.done) {
            graph_.visit(curVert, curArc);
            DEBUG_1(<< curVert << " visit arc " << curArc << std::endl);
            return nullSuperArc;

         } else {
            return graph_.getArcId(curVert);
         }
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
