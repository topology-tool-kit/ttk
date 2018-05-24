#ifndef FTRGRAPHPRIVATE_TEMPLATE_H
#define FTRGRAPHPRIVATE_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localPropagation)
      {
         // skeleton
         const idNode     downNode   = graph_.makeNode(seed);
         const idSuperArc currentArc = graph_.openArc(downNode, localPropagation);
         graph_.visit(seed, currentArc);
         std::cout << "open arc " << graph_.printArc(currentArc) << std::endl;

         // topology
         bool isJoinSadlleLast = false;
         bool isJoinSaddle = false, isSplitSaddle = false;

         // containers
         std::vector<idEdge>               lowerStarEdges, upperStarEdges;
         std::set<DynGraphNode<idVertex>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localPropagation->empty()) {
            localPropagation->getNextVertex();

            if (graph_.hasVisited(localPropagation->getCurVertex(), currentArc)) {
               continue;
            }

            // Mark this vertex with the current growth
            graph_.visit(localPropagation->getCurVertex(), currentArc);

            // Debug print
            std::cout << "visit: " << localPropagation->getCurVertex() << " : "
                      << graph_.printVisit(localPropagation->getCurVertex()) << std::endl;

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localPropagation);

            lowerComp = lowerComps(lowerStarEdges);
            if(lowerComp.size() > 1){
               isJoinSaddle = true;
            }

            if (isJoinSaddle) {
               isJoinSadlleLast = checkLast(currentArc, localPropagation, lowerStarEdges);
               std::cout << "Join " << lowerComp.size() << " isJoinSadlleLast: " << isJoinSadlleLast << std::endl;
               // If the current growth reaches a saddle and is not the last
               // reaching this saddle, it just stops here.
               if (!isJoinSadlleLast)
                  break;
            }

            updatePreimage(localPropagation);

            upperComp = upperComps(upperStarEdges);
            if (upperComp.size() > 1) {
               std::cout << "Split" << std::endl;
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
               growthFromSeed(upVert, newLocalProp);
            }
            // dynGraph_.print();
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
          const std::vector<idEdge>& finishingEdges)
      {
         return dynGraph_.findRoot(finishingEdges);
      }

      template <typename ScalarType>
      std::set<DynGraphNode<idVertex>*> FTRGraph<ScalarType>::upperComps(
          const std::vector<idEdge>& startingEdges)
      {
         return dynGraph_.findRoot(startingEdges);
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
         bool t = dynGraph_.insertEdge(std::get<0>(oTriangle), std::get<1>(oTriangle), w);

         if (t) {
            // std::cout << "start add edge: " << printEdge(std::get<0>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl;
         } else {
            // std::cout << "no need to create edge: " << printEdge(std::get<0>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl;

            // dynGraph_.print();
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageMiddleCell(const orderedTriangle&   oTriangle,
                                                          const Propagation* const localPropagation)
      {
         // Check if exist ?
         // If not, the triangle will be visited again once a merge have occured.
         // So we do not add the edge now
         const int t = dynGraph_.removeEdge(std::get<0>(oTriangle), std::get<1>(oTriangle));
         if (t) {
            // std::cout << "mid replace edge: " << printEdge(std::get<0>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl;

         }
         else {
            // std::cout << "mid no found edge: " << printEdge(std::get<0>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localPropagation) << std::endl;
         }

         const orderedEdge e0 = getOrderedEdge(std::get<0>(oTriangle), localPropagation);
         const orderedEdge e1 = getOrderedEdge(std::get<1>(oTriangle), localPropagation);
         const idVertex    w  = getWeight(e0, e1, localPropagation);
         const int u = dynGraph_.insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), w);

         if (u) {
            // std::cout << " new edge: " << printEdge(std::get<1>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl;
         } else {
            // std::cout << "no need to create edge: " << printEdge(std::get<1>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl;

            // dynGraph_.print();
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageEndCell(const orderedTriangle&   oTriangle,
                                                       const Propagation* const localPropagation)
      {
         const int t =
         dynGraph_.removeEdge(std::get<1>(oTriangle), std::get<2>(oTriangle));
         if (t) {
            // std::cout << "end remove edge: " << printEdge(std::get<1>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl;
         } else {
            // std::cout << "end not found edge: " << printEdge(std::get<1>(oTriangle), localPropagation);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localPropagation) << std::endl;
         }
      }

      template <typename ScalarType>
      idNode FTRGraph<ScalarType>::updateReebGraph(const idSuperArc         currentArc,
                                                   const Propagation* const localPropagation)
      {
         const idVertex upVert = localPropagation->getCurVertex(); // keep before merge
         const idNode upNode = graph_.makeNode(upVert);
         graph_.closeArc(currentArc, upNode);

         std::cout << "close arc " << graph_.printArc(currentArc) << std::endl;

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

         if(curSaddle == 3) {
            volatile int t = 0;
         }

         valence decr = 0;
         for(idVertex nid = 0; nid < nbNeigh; ++nid) {
            idVertex neighId;
            mesh_->getVertexNeighbor(curSaddle, nid, neighId);

            if (localPropagation->compare(neighId, curSaddle)) {
               if (graph_.hasVisited(neighId, rpz)) {
                  ++decr;
               }
            }
         }

         valence oldVal = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
         {
            oldVal = graph_.val(curSaddle);
            graph_.val(curSaddle) -= decr;
         }

         if (oldVal == -1) {
            // First task to touch this saddle, compute the valence
            idVertex totalVal = lowerStarEdges.size();
            valence newVal = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            {
               newVal = graph_.val(curSaddle);
               graph_.val(curSaddle) += (totalVal + 1);
            }
            oldVal = decr + newVal + (totalVal + 1);
         }

         return oldVal == decr;
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

         // find seeds (triangles)
         std::set<idCell> triangleSeeds = upCCtriangleSeeds(curVert, localProp);

         const idCell nbSeed = triangleSeeds.size();
         newLocalProps.reserve(nbSeed);

         // BFS to add vertices in the current propagation for each seed
         for (const idCell curSeed : triangleSeeds) {
            Propagation* curProp = newPropagation(curVert);
            newLocalProps.emplace_back(curProp);
            // fill curProp using a BFS on the current seed
            std::set<idCell>   visitedCells;
            std::set<idVertex> addedVertices;
            bfsPropagation(curVert, curSeed, curProp, visitedCells, addedVertices);
            // remove the already processed first vertex
            curProp->getNextVertex();
         }

         return newLocalProps;
      }

      template <typename ScalarType>
      std::set<idCell> FTRGraph<ScalarType>::upCCtriangleSeeds(const idVertex           v,
                                                               const Propagation* const localProp)
      {
         std::vector<idCell> triangles;
         std::vector<valence> cc;
         // all triangle in eighborhood
         const idCell nbTriNeigh = mesh_->getVertexTriangleNumber(v);
         triangles.reserve(nbTriNeigh);
         for(idCell t = 0; t < nbTriNeigh; ++t) {
            idCell neighTriangle;
            mesh_->getVertexTriangle(v, t, neighTriangle);
            const orderedTriangle oNeighTriangle = getOrderedTriangle(neighTriangle, localProp);
            // only if v is not the highest point
            if (getVertPosInTriangle(oNeighTriangle, localProp) != vertPosInTriangle::End) {
               triangles.emplace_back(neighTriangle);
            }
         }
         const idCell nbTriStar = triangles.size();
         cc.resize(nbTriStar, -1);

         // mark CC
         for (idCell t = 0; t < nbTriStar; ++t) {
            bfsSeed(t, t, triangles, cc, localProp);
         }

         // keep seeds only
         std::set<idCell> seeds;
         std::set<valence> seenCC;
         for (idCell t = 0; t < nbTriStar; ++t) {
            if (seenCC.find(cc[t]) == end(seenCC)) {
               seenCC.emplace(cc[t]);
               seeds.emplace(triangles[t]);
            }
         }
         return seeds;
      }

      /// Tools

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::newPropagation(const idVertex leaf)
      {
         auto compare_fun = [&](idVertex a, idVertex b) { return scalars_->isHigher(a, b); };
         Propagation* localPropagation(new Propagation(leaf, compare_fun));
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
            return scalars_->getMirror(end0);
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

         std::sort(begin(oEdges), end(oEdges), compareOEdges);

         return {std::get<2>(oEdges[0]), std::get<2>(oEdges[1]), std::get<2>(oEdges[2]), cellId};
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
   }
}

#endif /* end of include guard: FTRGRAPHPRIVATE_TEMPLATE_H */
