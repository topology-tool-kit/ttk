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
         DEBUG_1(<< "Start " << seed << " go up " << localProp->goUp() << std::endl);
         DEBUG_1(<< localProp->print() << " id " << localProp->getId() << std::endl);

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
            if (propagations_.hasVisited(curVert, localProp)) {
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
               graph_.visit(curVert, currentArc);
               DEBUG_1(<< "visit n: " << curVert << " arc " << currentArc << std::endl);
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
               DEBUG_1(<< ": is split : " << localProp->print() << std::endl);
               isSplitSaddle = true;
            }

            // add upper star for futur visit
            bool seenAbove = localGrowth(localProp);
            if (!seenAbove) {
               // We have reached an opposite leaf
               const idNode upNode = graph_.makeNode(curVert);
               graph_.closeArc(currentArc, upNode);
               DEBUG_1(<< "close arc max " << graph_.printArc(currentArc) << std::endl);
#ifdef TTK_ENABLE_FTR_STATS
               if (localProp->empty()) {
                  idVertex curProp;
#pragma omp atomic capture
                  {
                     curProp = nbProp_;
                     nbProp_--;
                  }
                  propTimes_[curProp - 1] = sweepStart_.getElapsedTime();
               }
#endif
            }
         } // end propagation while

         // get the corresponging critical point on which
         // the propagation has stopped (join, split, max)
         const idVertex upVert = localProp->getCurVertex();
         const idNode   upNode = graph_.makeNode(upVert);

         // At saddle: join or split or both

         idSuperArc joinNewArc;
         // arriving at a join
         if (isJoinSadlleLast) {
            // ensure we have the good values here, even if other tasks were doing stuff
            {
               // here to solve a 1 over thousands execution bug
               std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);
               lowerComp = lowerComps(lowerStarEdges, localProp);
            }
            localGrowth(localProp);
            mergeAtSaddle(upNode, localProp, lowerComp);
            const idNode downNode = graph_.getNodeId(upVert);
            joinNewArc            = graph_.openArc(downNode, localProp);
            updatePreimage(localProp, joinNewArc);
            graph_.visit(upVert, joinNewArc);
            upperComp = upperComps(upperStarEdges, localProp);
            if (upperComp.size() > 1) {
               DEBUG_1(<< ": is split : " << localProp->print() << std::endl);
               isSplitSaddle = true;
               // will be replaced be new arcs of the split
               graph_.getArc(joinNewArc).hide();
            }
         }

#ifdef TTK_ENABLE_FTR_STATS
         if (isJoinSaddle && !isJoinSadlleLast) {
            // This propagation is dying here
            idVertex curProp;
#pragma omp atomic capture
            {
               curProp = nbProp_;
               nbProp_--;
            }
            propTimes_[curProp - 1] = sweepStart_.getElapsedTime();
         }
#endif

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
      void FTRGraph<ScalarType>::growthFromSeedWithLazy(const idVertex seed, Propagation* localProp,
                                                        idSuperArc currentArc)
      {

         DEBUG_1(<< "Start " << seed << " go up " << localProp->goUp() << std::endl);
         DEBUG_1(<< localProp->print() << " id " << localProp->getId() << std::endl);

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle     = false, isSplitSaddle = false;

         // containers
         // vitsit
         std::vector<idEdge>                  lowerStarEdges, upperStarEdges;
         std::vector<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localProp->empty()) {
            localProp->nextVertex();
            const idVertex curVert = localProp->getCurVertex();

            // Check history for visit (quick test)
            if (propagations_.hasVisited(curVert, localProp)) {
               DEBUG_1(<< "already seen " << curVert << " " << localProp->getId() << std::endl);
               continue;
            }

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);

            if (valences_.lower[curVert] < 2 && valences_.upper[curVert] < 2) {

               // simple reeb regular, lazyness
               if(lowerStarEdges.size()) {
                  // not a min nor a saddle: 1 CC below
                  currentArc = dynGraph(localProp).getSubtreeArc(lowerStarEdges[0]);
               }
               for (const idEdge dgNode : upperStarEdges) {
                  dynGraph(localProp).setCorArc(dgNode, currentArc);
               }
               graph_.visit(curVert, currentArc);
               DEBUG_1(<< "visit n: " << curVert << " arc " << currentArc << std::endl);
               lazyUpdatePreimage(localProp, currentArc);
            } else {
               // process lazy
               // sort both list
               // for each:
               // if del < add: delete in real RG
               // if add < del: add in real RG
               // else: drop both, computation avoided
               lazyApply(localProp);
               lowerComp = lowerComps(lowerStarEdges, localProp);
               if (lowerComp.size() > 1) {
                  isJoinSaddle = true;
                  isJoinSadlleLast = checkLast(currentArc, localProp, lowerStarEdges);
                  break;
               } else {
                  currentArc = lowerComp[0]->getCorArc();
                  graph_.visit(curVert, currentArc);
               }
               updatePreimage(localProp, currentArc);
               upperComp = upperComps(upperStarEdges, localProp);
               if (upperComp.size() > 1) {
                  isSplitSaddle = true;
               }
            }


            // add upper star for futur visit
            bool seenAbove = localGrowth(localProp);
            if (!seenAbove) {
               // We have reached an opposite leaf
               const idNode upNode = graph_.makeNode(curVert);
               graph_.closeArc(currentArc, upNode);
               DEBUG_1(<< "close arc max " << graph_.printArc(currentArc) << std::endl);
#ifdef TTK_ENABLE_FTR_STATS
               if (localProp->empty()) {
                  idVertex curProp;
#pragma omp atomic capture
                  {
                     curProp = nbProp_;
                     nbProp_--;
                  }
                  propTimes_[curProp - 1] = sweepStart_.getElapsedTime();
               }
#endif
            }
         } // end propagation while

         if (isJoinSaddle) {
            DEBUG_1(<< ": is join " << isJoinSadlleLast << std::endl);
         }
         if (isSplitSaddle) {
            DEBUG_1(<< ": is split " << std::endl);
         }
         // get the corresponging critical point on which
         // the propagation has stopped (join, split, max)
         const idVertex upVert = localProp->getCurVertex();
         const idNode   upNode = graph_.makeNode(upVert);

         // At saddle: join or split or both

         idSuperArc joinNewArc;
         // arriving at a join
         if (isJoinSadlleLast) {
            // ensure we have the good values here, even if other tasks were doing stuff
            {
               // here to solve a 1 over thousands execution bug in parallel
               // std::tie(lowerStarEdges, upperStarEdges) = visitStar(localProp);
               // lowerComp = lowerComps(lowerStarEdges, localProp);
            }
            localGrowth(localProp);
            mergeAtSaddle(upNode, localProp, lowerComp);
            const idNode downNode = graph_.getNodeId(upVert);
            joinNewArc            = graph_.openArc(downNode, localProp);
            graph_.visit(upVert, joinNewArc);
            updatePreimage(localProp, joinNewArc);
            upperComp = upperComps(upperStarEdges, localProp);
            if (upperComp.size() > 1) {
               isSplitSaddle = true;
               DEBUG_1(<< ": is joina & split : " << localProp->print() << std::endl);
               // will be replaced be new arcs of the split
               graph_.getArc(joinNewArc).hide();
            }
         }

#ifdef TTK_ENABLE_FTR_STATS
         if (isJoinSaddle && !isJoinSadlleLast) {
            // This propagation is dying here
            idVertex curProp;
#pragma omp atomic capture
            {
               curProp = nbProp_;
               nbProp_--;
            }
            propTimes_[curProp - 1] = sweepStart_.getElapsedTime();
         }
#endif

         // starting from the saddle
         if (isSplitSaddle && (!isJoinSaddle || isJoinSadlleLast)) {
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
            growthFromSeedWithLazy(upVert, localProp, joinNewArc);
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
         localProp->lazyAdd(std::get<0>(oTriangle), std::get<1>(oTriangle), curArc);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyMiddle(const orderedTriangle& oTriangle,
                                                  Propagation* const     localProp,
                                                  const idSuperArc       curArc)
      {
         localProp->lazyDel(std::get<0>(oTriangle), std::get<1>(oTriangle));
         localProp->lazyAdd(std::get<1>(oTriangle), std::get<2>(oTriangle), curArc);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyEnd(const orderedTriangle& oTriangle,
                                               Propagation* const     localProp,
                                               const idSuperArc       curArc)
      {
         localProp->lazyDel(std::get<1>(oTriangle), std::get<2>(oTriangle));
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyAdd(const Propagation* const                localProp,
                                               const std::tuple<linkEdge, idSuperArc>& add)
      {
         const linkEdge    curLink = std::get<0>(add);
         const orderedEdge e0      = mesh_.getOrderedEdge(std::get<0>(curLink), localProp->goUp());
         const orderedEdge e1      = mesh_.getOrderedEdge(std::get<1>(curLink), localProp->goUp());
         const idVertex    w       = getWeight(e0, e1, localProp);
         bool              t;
         t = dynGraph(localProp).insertEdge(std::get<1>(curLink), std::get<0>(curLink), w,
                                            std::get<1>(add));
         if (t) {
            DEBUG_1(<< "start add edge: " << printEdge(std::get<0>(curLink), localProp));
            DEBUG_1(<< " :: " << printEdge(std::get<1>(curLink), localProp) << " w " << w << " arc " << std::get<1>(add) << std::endl);
         } else {
            DEBUG_1(<< "start no need to create edge: "
                    << printEdge(std::get<0>(curLink), localProp) << std::endl);
            // DEBUG_1(<< " :: " << printEdge(std::get<1>(curLink), localProp) << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateLazyDel(const Propagation* const                localProp,
                                               const std::tuple<linkEdge, idSuperArc>& del)
      {
         const linkEdge curLink = std::get<0>(del);
         const int t = dynGraph(localProp).removeEdge(std::get<0>(curLink), std::get<1>(curLink));

         // keep history inside the dyngraph structure
         // dynGraph(localProp).setSubtreeArc(std::get<0>(curLink), std::get<1>(del));
         // dynGraph(localProp).setSubtreeArc(std::get<1>(curLink), std::get<1>(del));

         if (t) {
            DEBUG_1(<< "mid del edge: " << printEdge(std::get<0>(curLink), localProp));
            DEBUG_1(<< " :: " << printEdge(std::get<1>(curLink), localProp) << " : " << std::get<1>(del) << std::endl);
         }
         else {
            DEBUG_1(<< "mid no found edge: " << printEdge(std::get<0>(curLink), localProp));
            DEBUG_1(<< " :: " << printEdge(std::get<1>(curLink), localProp) << " : " << std::get<1>(del) << std::endl);
         }
      }

      template<typename ScalarType>
      void FTRGraph<ScalarType>::lazyApply(Propagation* const localProp)
      {
         auto comp = [localProp](const idVertex a, const idVertex b) {
            return localProp->compare(a, b);
         };

         auto compEdges = [&, comp](const linkEdge& a, const linkEdge& b) {
            return mesh_.compareLinks(a, b, comp);
         };

         DEBUG_1(<< "lazy apply " << localProp->getCurVertex() << std::endl);

         // use a real sort to avoid problem linked to the tree structure of the DG
         // + also look at the Weight bug
         localProp->sortLazyLists(compEdges);

         auto add = localProp->lazyAddNext();
         auto del = localProp->lazyDelNext();
         while (std::get<0>(add) != nullLink || std::get<0>(del) != nullLink) {
            if(std::get<0>(del) == nullLink) {
               updateLazyAdd(localProp, add);
               add = localProp->lazyAddNext();
            } else if (std::get<0>(add) == nullLink || mesh_.compareLinks(std::get<0>(del), std::get<0>(add), comp)) {
               updateLazyDel(localProp, del);
               del = localProp->lazyDelNext();
            } else if (mesh_.compareLinks(std::get<0>(add), std::get<0>(del), comp)) {
               updateLazyAdd(localProp, add);
               add = localProp->lazyAddNext();
            } else {
               // same arc in both list, should be added and removed so we jus ignore it
               // (add and del cant be null both of them at the same time)
               add = localProp->lazyAddNext();
               del = localProp->lazyDelNext();
            }
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
                  // DEBUG_1(<< " + " << neighId << std::endl);
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
         AtomicUF*      curId     = localProp->getId();
         valence        decr      = 0;

          DEBUG_1(<< "Check last on " << curSaddle << " id " << curId << std::endl);

         // NOTE:
         // Using propagation id allows to decrement by the number of time this propagation
         // has reached the saddle, even if the propagation take care of several of these arcs
         // (after a Hole-split).
         for(idEdge edgeId : lowerStarEdges) {
            // lowerStarEdge already conatins roots
           const idSuperArc    edgeArc = dynGraph(localProp).getCorArc(edgeId);
           if (edgeArc == nullSuperArc) // ignore unseen
              continue;
           AtomicUF* tmpId = graph_.getArc(edgeArc).getPropagation()->getId();
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
               const idSuperArc rpzArc = dynGraph(!localProp->goUp()).getCorArc(edgeId);
               // const idSuperArc rpzArc = dynGraph(!localProp->goUp()).getSubtreeArc(edgeId);
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

         DEBUG_1(<< " merge in " << graph_.getNode(saddleId).getVertexIdentifier() << std::endl);

         for(auto* dgNode : lowerComp) {
            // read in the history (lower comp already contains roots)
            const idSuperArc endingArc = dgNode->getCorArc();
            graph_.closeArc(endingArc, saddleId);
            Propagation * arcProp =  graph_.getArc(endingArc).getPropagation();
            DEBUG_1(<< "merge " << graph_.printArc(endingArc) << " prop " << arcProp << std::endl);
            localProp->merge(*arcProp);
         }
         DEBUG_1(<< " result " << localProp->print() << std::endl);
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
            growthFromSeedWithLazy(curVert, prop, arc);
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
            growthFromSeedWithLazy(curVert, localProp);
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
