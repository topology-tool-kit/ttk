#ifndef FTRGRAPHPRIVATE_TEMPLATE_H
#define FTRGRAPHPRIVATE_TEMPLATE_H

#include "FTRGraph.h"

// Skeleton + propagation
// #define DEBUG_1(msg) std::cout msg
#define DEBUG_1(msg)

// Dynamic graph structure
// #define DEBUG_2(msg) std::cout msg
#define DEBUG_2(msg)

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localPropagation)
      {
         DEBUG_1(<< "Start " << seed << " go up " << localPropagation->goUp() << std::endl);
         DEBUG_1(<< localPropagation->getRpz() << " " << localPropagation->print() << std::endl);

         // Check if next vertex is already visited by an arc coming here
         const bool alreadyAttached =
             checkAlreayAttached(seed, localPropagation->getNextVertex(), localPropagation);
         if (alreadyAttached)
            return;

         // skeleton
         const idNode     downNode   = graph_.makeNode(seed);
         const idSuperArc currentArc = graph_.openArc(downNode, localPropagation);
         graph_.visit(seed, currentArc);

#ifndef NDEBUG
         graph_.getArc(currentArc).setFromUp(localPropagation->goUp());
#endif

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle = false, isSplitSaddle = false;

         // containers
         std::vector<idEdge>               lowerStarEdges, upperStarEdges;
         std::set<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localPropagation->empty()) {
            localPropagation->nextVertex();
            const idVertex curVert = localPropagation->getCurVertex();

            // Avoid revisiting things processed by this CC
            // CAUTION: Diamond case here !
            if (graph_.hasVisited(curVert, localPropagation->getRpz())) {
               continue;
            }

            // Mark this vertex with the current growth
            graph_.visit(curVert, currentArc);

            DEBUG_1(<< "visit: " << curVert << std::endl);

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

            updatePreimage(localPropagation);

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
            mergeAtSaddle(upNode);
         }

         if (isSplitSaddle) {
            std::vector<Propagation*> newProps = splitAtSaddle(localPropagation);
            for (Propagation* newLocalProp : newProps) {
               if (!newLocalProp->empty()) {
                  growthFromSeed(upVert, newLocalProp);
               }
            }
         } else if (isJoinSadlleLast) {
            // recursive call
           growthFromSeed(upVert, localPropagation);
         }
      }

      template <typename ScalarType>
      std::pair<std::vector<idEdge>, std::vector<idEdge>> FTRGraph<ScalarType>::visitStar(
          const Propagation* const localPropagation) const
      {
         // TODO re-use the same vectors per thread
         std::vector<idEdge> lowerStar, upperStar;

         const idEdge nbAdjEdges = mesh_->getVertexEdgeNumber(localPropagation->getCurVertex());
         lowerStar.reserve(nbAdjEdges);
         upperStar.reserve(nbAdjEdges);

         for (idEdge e = 0; e < nbAdjEdges; ++e) {
            idEdge edgeId;
            mesh_->getVertexEdge(localPropagation->getCurVertex(), e, edgeId);
            idVertex edgeLowerVert, edgeUpperVert;
            std::tie(edgeLowerVert, edgeUpperVert, std::ignore) =
                getOrderedEdge(edgeId, localPropagation);
            if (edgeLowerVert == localPropagation->getCurVertex()) {
               upperStar.emplace_back(edgeId);
            } else {
               lowerStar.emplace_back(edgeId);
            }
         }

         return {lowerStar, upperStar};
      }

      template <typename ScalarType>
      std::set<DynGraphNode<idVertex>*> FTRGraph<ScalarType>::lowerComps(
          const std::vector<idEdge>& finishingEdges, const Propagation* const localProp)
      {
         return dynGraph(localProp).findRoot(finishingEdges);
      }

      template <typename ScalarType>
      std::set<DynGraphNode<idVertex>*> FTRGraph<ScalarType>::upperComps(
          const std::vector<idEdge>& startingEdges, const Propagation* const localProp)
      {
         return dynGraph(localProp).findRoot(startingEdges);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimage(const Propagation* const localPropagation)
      {
         const idCell nbAdjTriangles =
             mesh_->getVertexTriangleNumber(localPropagation->getCurVertex());

         for (idCell t = 0; t < nbAdjTriangles; ++t) {
            // Classify current cell
            idCell curTriangleid;
            mesh_->getVertexTriangle(localPropagation->getCurVertex(), t, curTriangleid);

            orderedTriangle   oTriangle  = getOrderedTriangle(curTriangleid, localPropagation);
            vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localPropagation);

            // std::cout << "update preimage at v" << localPropagation->getCurVertex() << " : "
            //           << printTriangle(oTriangle, localPropagation) << std::endl;

            // Update DynGraph
            // We can have an end pos on an unvisited triangle
            // in case of saddle points
            switch (curVertPos) {
               case vertPosInTriangle::Start:
                  updatePreimageStartCell(oTriangle, localPropagation);
                  break;
               case vertPosInTriangle::Middle:
                  updatePreimageMiddleCell(oTriangle, localPropagation);
                  break;
               case vertPosInTriangle::End:
                  updatePreimageEndCell(oTriangle, localPropagation);
                  break;
               default:
                  std::cout << "[FTR]: update preimage error, unknown vertPos type" << std::endl;
                  break;
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageStartCell(const orderedTriangle&   oTriangle,
                                                         const Propagation* const localPropagation)
      {
         const orderedEdge e0 = getOrderedEdge(std::get<0>(oTriangle), localPropagation);
         const orderedEdge e1 = getOrderedEdge(std::get<1>(oTriangle), localPropagation);
         const idVertex    w  = getWeight(e0, e1, localPropagation);
         bool t = dynGraph(localPropagation).insertEdge(std::get<0>(oTriangle), std::get<1>(oTriangle), w);

         if (t) {
            DEBUG_2(<< "start add edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< "no need to create edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageMiddleCell(const orderedTriangle&   oTriangle,
                                                          const Propagation* const localPropagation)
      {
         // Check if exist ?
         // If not, the triangle will be visited again once a merge have occured.
         // So we do not add the edge now
         const int t = dynGraph(localPropagation).removeEdge(std::get<0>(oTriangle), std::get<1>(oTriangle));

         if (t) {
            DEBUG_2(<< "mid replace edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
         }
         else {
            DEBUG_2(<< "mid no found edge: " << printEdge(std::get<0>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }

         const orderedEdge e1 = getOrderedEdge(std::get<1>(oTriangle), localPropagation);
         const orderedEdge e2 = getOrderedEdge(std::get<2>(oTriangle), localPropagation);
         const idVertex    w  = getWeight(e1, e2, localPropagation);
         const int u = dynGraph(localPropagation).insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), w);

         if (u) {
            DEBUG_2(<< " new edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< "no need to create edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
            DEBUG_2(<< dynGraph(localPropagation).print() << std::endl);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageEndCell(const orderedTriangle&   oTriangle,
                                                       const Propagation* const localPropagation)
      {
         const int t = dynGraph(localPropagation).removeEdge(std::get<1>(oTriangle), std::get<2>(oTriangle));

         if (t) {
            DEBUG_2(<< "end remove edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
         } else {
            DEBUG_2(<< "end not found edge: " << printEdge(std::get<1>(oTriangle), localPropagation));
            DEBUG_2(<< " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl);
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

         return upNode;
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::localGrowth(Propagation* const localPropagation)
      {
         const idVertex nbNeigh = mesh_->getVertexNeighborNumber(localPropagation->getCurVertex());
         for (idVertex n = 0; n < nbNeigh; ++n) {
            idVertex neighId;
            mesh_->getVertexNeighbor(localPropagation->getCurVertex(), n, neighId);
            if (localPropagation->compare(localPropagation->getCurVertex(), neighId)) {
               if (!toVisit_[neighId] || toVisit_[neighId]->find() != localPropagation->getRpz()) {
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
         const idVertex   nbNeigh   = mesh_->getVertexNeighborNumber(curSaddle);
         const UnionFind* rpz       = localPropagation->getRpz();

         valence decr = 0;
         for(idVertex nid = 0; nid < nbNeigh; ++nid) {
            idVertex neighId;
            mesh_->getVertexNeighbor(curSaddle, nid, neighId);

            if (localPropagation->compare(neighId, curSaddle)) {
               if (graph_.hasVisited(neighId, rpz)) {
                  DEBUG_1(<< neighId << " decrement " << curSaddle << std::endl);
                  ++decr;
               }
            }
         }

         valence oldVal = 0;
         if (localPropagation->goUp()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            {
               oldVal = graph_.valDown(curSaddle);
               graph_.valDown(curSaddle) -= decr;
            }

         } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            {
               oldVal = graph_.valUp(curSaddle);
               graph_.valUp(curSaddle) -= decr;
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
                  newVal = graph_.valDown(curSaddle);
                  graph_.valDown(curSaddle) += (totalVal + 1);
               }
            } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
               {
                  newVal = graph_.valUp(curSaddle);
                  graph_.valUp(curSaddle) += (totalVal + 1);
               }
            }
            oldVal = decr + newVal + (totalVal + 1);
         }

         return oldVal == decr;
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkAlreayAttached(const idVertex saddle, const idVertex neigh, const Propagation* const localProp)
      {
         const idNode checkNode = graph_.getNodeId(saddle);
         for(const idSegmentation s : graph_.visit(neigh)) {
            if(s >= 0) {
               if (graph_.getArc(s).getUpNodeId() == checkNode) {
                  const idNode downNode = graph_.getArc(s).getDownNodeId();
                  if (localProp->compare(saddle, graph_.getNode(downNode).getVertexIdentifier())) {
                     return true;
                  }
               }
            }
         }
         return false;
      }

      template<typename ScalarType>
      void FTRGraph<ScalarType>::mergeAtSaddle(const idNode saddleId)
      {
         graph_.mergeAtSaddle(saddleId);
      }

      template<typename ScalarType>
      std::vector<Propagation*> FTRGraph<ScalarType>::splitAtSaddle(const Propagation* const localProp)
      {
         std::vector<Propagation*> newLocalProps;
         const idVertex            curVert = localProp->getCurVertex();

         newLocalProps.reserve(4);
         const idCell nbTriNeigh = mesh_->getVertexTriangleNumber(curVert);
         for(idCell t = 0; t < nbTriNeigh; ++t) {
            idCell neighTriangle;
            mesh_->getVertexTriangle(curVert, t, neighTriangle);
            const orderedTriangle oNeighTriangle = getOrderedTriangle(neighTriangle, localProp);
            // only if curVert is not the highest point
            if (bfsCells_[neighTriangle] != curVert &&
                getVertPosInTriangle(oNeighTriangle, localProp) != vertPosInTriangle::End) {

               const idVertex endTri          = getEndVertexInTriangle(oNeighTriangle, localProp);
               const bool     alreadyAttached = checkAlreayAttached(curVert, endTri, localProp);
               if(alreadyAttached) continue;

               // BFS to add vertices in the current propagation for each seed
               Propagation* curProp = newPropagation(curVert, localProp->goUp()/*, localProp->getRpz()*/);
               newLocalProps.emplace_back(curProp);
               // False arc, only used to maintain the history in terms of UF while
               // doing the BFS
               const idSuperArc arcForUf = graph_.makeHiddenArc(curProp);
               // fill curProp using a BFS on the current seed
               bfsPropagation(curVert, neighTriangle, curProp, arcForUf);
               // remove the already processed first vertex
               curProp->nextVertex();
            }
         }

         return newLocalProps;
      }

      /// Tools

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::newPropagation(const idVertex leaf, const bool fromMin,
                                                        UnionFind* rpz)
      {
         Propagation* localPropagation;
         if (fromMin) {
            auto compare_max_fun = [&](idVertex a, idVertex b) { return scalars_->isHigher(a, b); };
            localPropagation = new Propagation(leaf, compare_max_fun, fromMin, rpz);
         } else {
            auto compare_min_fun = [&](idVertex a, idVertex b) { return scalars_->isLower(a, b); };
            localPropagation = new Propagation(leaf, compare_min_fun, fromMin, rpz);
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
      orderedEdge FTRGraph<ScalarType>::getOrderedEdge(const idEdge             edgeId,
                                                       const Propagation* const localPropagation) const
      {
         idVertex edge0Vert, edge1Vert;
         mesh_->getEdgeVertex(edgeId, 0, edge0Vert);
         mesh_->getEdgeVertex(edgeId, 1, edge1Vert);

         if (localPropagation->compare(edge0Vert, edge1Vert)) {
            return {edge0Vert, edge1Vert, edgeId};
         } else {
            return {edge1Vert, edge0Vert, edgeId};
         }
      }

      template <typename ScalarType>
      orderedTriangle FTRGraph<ScalarType>::getOrderedTriangle(
          const idCell cellId, const Propagation* const localPropagation) const
      {
         idEdge edges[3];
         mesh_->getTriangleEdge(cellId, 0, edges[0]);
         mesh_->getTriangleEdge(cellId, 1, edges[1]);
         mesh_->getTriangleEdge(cellId, 2, edges[2]);

         orderedEdge oEdges[3];
         oEdges[0] = getOrderedEdge(edges[0], localPropagation);
         oEdges[1] = getOrderedEdge(edges[1], localPropagation);
         oEdges[2] = getOrderedEdge(edges[2], localPropagation);

         auto compareOEdges = [localPropagation](const orderedEdge& a, const orderedEdge& b) {
            return localPropagation->compare(std::get<0>(a), std::get<0>(b)) ||
                   (std::get<0>(a) == std::get<0>(b) &&
                    localPropagation->compare(std::get<1>(a), std::get<1>(b)));
         };

         if (compareOEdges(oEdges[0], oEdges[1])) {
            // 1 2 3
            // 1 3 2
            // 2 3 1

            if (compareOEdges(oEdges[1], oEdges[2])) {
               // 1 2 3
               return {std::get<2>(oEdges[0]), std::get<2>(oEdges[1]), std::get<2>(oEdges[2]),
                       cellId};
            }
            // 1 3 2
            // 2 3 1

            if (compareOEdges(oEdges[0], oEdges[2])) {
               // 1 3 2
               return {std::get<2>(oEdges[0]), std::get<2>(oEdges[2]), std::get<2>(oEdges[1]), cellId};
            }

            // 2 3 1
            return {std::get<2>(oEdges[2]), std::get<2>(oEdges[0]), std::get<2>(oEdges[1]), cellId};
         }

         // 2 1 3
         // 3 2 1
         // 3 1 2

         if(compareOEdges(oEdges[0], oEdges[2])) {
            // 2 1 3
            return {std::get<2>(oEdges[1]), std::get<2>(oEdges[0]), std::get<2>(oEdges[2]), cellId};
         }

         // 3 2 1
         // 3 1 2

         if(compareOEdges(oEdges[1], oEdges[2])) {
            // 3 1 2
            return {std::get<2>(oEdges[1]), std::get<2>(oEdges[2]), std::get<2>(oEdges[0]), cellId};
         }

         // 3 2 1
         return {std::get<2>(oEdges[2]), std::get<2>(oEdges[1]), std::get<2>(oEdges[0]), cellId};
      }

      template <typename ScalarType>
      vertPosInTriangle FTRGraph<ScalarType>::getVertPosInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localPropagation) const
      {
         orderedEdge firstEdge = getOrderedEdge(std::get<0>(oTriangle), localPropagation);
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
         const orderedEdge& higherEdge = getOrderedEdge(std::get<1>(oTriangle), localPropagation);
         return std::get<1>(higherEdge);
      }
   }
}

#endif /* end of include guard: FTRGRAPHPRIVATE_TEMPLATE_H */
