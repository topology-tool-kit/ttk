#ifndef FTRGRAPH_TEMPLATE_H
#define FTRGRAPH_TEMPLATE_H

#include "FTRGraph.h"
#include "Propagation.h"
#include "Tasks.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      FTRGraph<ScalarType>::FTRGraph()
          : params_(new Params), scalars_(new Scalars<ScalarType>), needDelete_(true)
      {
      }

      template <typename ScalarType>
      FTRGraph<ScalarType>::FTRGraph(Params* const params, Triangulation* mesh,
                                     Scalars<ScalarType>* const scalars)
          : params_(params), mesh_(mesh), scalars_(scalars), needDelete_(false)
      {
      }

      template <typename ScalarType>
      FTRGraph<ScalarType>::~FTRGraph()
      {
         if (needDelete_) {
            delete params_;
            delete scalars_;
         }

         for (Propagation* p : propagations_) {
            delete p;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::build()
      {
         Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
         if (!mesh_) {
            std::cerr << "[FTR Graph]: no mesh input" << std::endl;
            return;
         }
         if (!scalars_) {
            std::cerr << "[FTR Graph]: no scalars given" << std::endl;
            return;
         }
         if (!params_) {
            std::cerr << "[FTR Graph]: no parameters" << std::endl;
            return;
         }
#endif
         // init some values

#ifdef TTK_ENABLE_OPENMP
         omp_set_num_threads(threadNumber_);
         omp_set_nested(1);
#endif

         int vertexNumber = mesh_->getNumberOfVertices();
         scalars_->setSize(vertexNumber);
         graph_.setNumberOfVertices(vertexNumber);
         dynGraph_.setNumberOfNodes(mesh_->getNumberOfEdges());

         params_->printSelf();

         // Precompute
         DebugTimer timeAlloc;
         alloc();
         printTime(timeAlloc, "[FTR Graph]: alloc time: ", infoMsg);

         DebugTimer timeInit;
         init();
         printTime(timeInit, "[FTR Graph]: init time: ", infoMsg);

         DebugTimer timeSort;
         scalars_->sort();
         printTime(timeSort, "[FTR Graph]: sort time: ", infoMsg);

         // Build the graph

         DebugTimer timeBuild;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(params_->threadNumber)
#endif
         {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
             {
               DebugTimer timeLeafSearch;
               leafSearch();
               printTime(timeLeafSearch, "[FTR Graph]: leaf search time: ", timeMsg);

               DebugTimer timeSwipe;
               sweepFrowSeeds();
               printTime(timeSwipe, "[FTR Graph]: sweepFrowSeeds time: ", timeMsg);
            }
         }
         // Debug
         // dynGraph_.test();

         // Message user
         {
            std::stringstream msg;
            msg << "[FTR Graph] Data-set (" << vertexNumber << " points) processed in "
                << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << std::endl;
            dMsg(std::cout, msg.str(), timeMsg);
         }
      }  // namespace ttk

      // protected

      template <typename ScalarType>
      void FTRGraph<ScalarType>::leafSearch()
      {
         TaskChunk leafChunkParams(scalars_->getSize());
         leafChunkParams.grainSize = 10000;
         auto leafChunk            = Tasks::getChunk(leafChunkParams);

         for (idTask leafChunkId = 0; leafChunkId < std::get<1>(leafChunk); ++leafChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(leafChunkId)
#endif
            {
               const idVertex lowerBound = Tasks::getBegin(leafChunkId, std::get<0>(leafChunk));
               const idVertex upperBound =
                   Tasks::getEnd(leafChunkId, std::get<0>(leafChunk), scalars_->getSize());

               for (idVertex v = lowerBound; v < upperBound; ++v) {
                  const valence vNeighNumber = mesh_->getVertexNeighborNumber(v);
                  bool          isMinV       = true;

                  for (valence n = 0; n < vNeighNumber; ++n) {
                     idVertex neigh;
                     mesh_->getVertexNeighbor(v, n, neigh);

                     if (scalars_->isLower(neigh, v)) {
                        isMinV = false;
                        break;
                     }
                  }

                  // v is a minimum, add it to the leaves
                  if (isMinV) {
                     graph_.addLeaf(v);
                     graph_.makeNode(v);
                  }
               }
            }  // end task
         }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
         std::cout << "find: " << graph_.getNumberOfLeaves() << " leaves" << std::endl;
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::sweepFrowSeeds()
      {
         const idNode nbSeed = graph_.getNumberOfLeaves();
         graph_.sortLeaves<ScalarType>(scalars_);

         for (idNode i = 0; i < nbSeed; i++) {
            // TODO
            // #pragma omp task

            // initialize structure
            const idVertex corLeaf          = graph_.getLeaf(i);
            Propagation*   localPropagation = newPropagation(corLeaf);
            // process
            growthFromSeed(corLeaf, localPropagation);
         }

         printGraph(3);
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printEdge(const orderedEdge&       oEdge,
                                                  const Propagation* const localProp) const
      {
         std::stringstream res;
         res << "e" << std::get<2>(oEdge) << ":";
         res << "(";
         res << std::get<0>(oEdge) << " - " << std::get<1>(oEdge);
         res << ")";
         return res.str();
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printEdge(const idEdge             edgeId,
                                                  const Propagation* const localProp) const
      {
         return printEdge(getOrderedEdge(edgeId, localProp), localProp);
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printTriangle(const orderedTriangle&   oTriangle,
                                                      const Propagation* const localProp) const
      {
         std::stringstream res;
         orderedEdge       e0, e1, e2;

         e0 = getOrderedEdge(std::get<0>(oTriangle), localProp);
         e1 = getOrderedEdge(std::get<1>(oTriangle), localProp);
         e2 = getOrderedEdge(std::get<2>(oTriangle), localProp);

         res << "t" << std::get<3>(oTriangle) << ":";
         res << "{";
         res << " " << printEdge(e0, localProp);
         res << " " << printEdge(e1, localProp);
         res << " " << printEdge(e2, localProp);
         res << " }";
         return res.str();
      }

      template <typename ScalarType>
      std::string FTRGraph<ScalarType>::printTriangle(const idCell             cellId,
                                                      const Propagation* const localProp) const
      {
         return printTriangle(getOrderedTriangle(cellId, localProp), localProp);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::printGraph(const int verbosity) const
      {
         graph_.print(verbosity);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::printTime(DebugTimer& timer, const std::string& msg,
                                           const int lvl) const
      {
         std::ostringstream outString(std::string(lvl, ' '));
         outString << msg << timer.getElapsedTime() << std::endl;
         dMsg(std::cout, outString.str(), lvl);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::alloc()
      {
         scalars_->alloc();
         graph_.alloc();
         dynGraph_.alloc();

         propagations_.reserve(scalars_->getSize());
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::init()
      {
         scalars_->removeNaN();
         scalars_->init();
         graph_.init();
         dynGraph_.init();
      }

      // private

      template <typename ScalarType>
      void FTRGraph<ScalarType>::growthFromSeed(const idVertex seed, Propagation* localPropagation)
      {
         // skeleton
         const idNode     downNode   = graph_.makeNode(seed);
         const idSuperArc currentArc = graph_.openArc(downNode, localPropagation);
         std::cout << "open arc " << graph_.printArc(currentArc) << std::endl;

         // topology
         bool isLast       = false;
         bool isJoinSaddle = false, isSplitSaddle = false;

         // containers
         std::vector<idEdge>                 lowerStarEdges, upperStarEdges;
         std::set<DynGraphNode<ScalarType>*> lowerComp, upperComp;

         while (!isJoinSaddle && !isSplitSaddle && !localPropagation->empty()) {
            localPropagation->getNextVertex();

            std::cout << "vert " << localPropagation->getCurVertex() << std::endl;

            // Debug print
            if (graph_.isVisited(localPropagation->getCurVertex())) {
               const auto vv = localPropagation->getCurVertex();
               std::cout << "Revisit: " << vv << " : " << graph_.getFirstVisit(vv) << std::endl;
            }
            // Mark this vertex with the current growth
            graph_.visit(localPropagation->getCurVertex(), currentArc);

            lowerStarEdges.clear();
            upperStarEdges.clear();
            std::tie(lowerStarEdges, upperStarEdges) = visitStar(localPropagation);

            lowerComp = lowerComps(lowerStarEdges);
            if(lowerComp.size() > 1){
               isJoinSaddle = true;
            }

            if (isJoinSaddle) {
               isLast = checkLast(localPropagation, lowerStarEdges);
               std::cout << "isLast: " << isLast << std::endl;
               // If the current growth reaches a saddle and is not the last
               // reaching this saddle, it just stops here.
               if (!isLast)
                  return;
            }

            updatePreimage(localPropagation);

            upperComp = upperComps(upperStarEdges);
            if (upperComp.size() > 1) {
               isSplitSaddle = true;
            }

            // add upper star for futur visit
            localGrowth(localPropagation);
         }

         // // Debug
         // std::cout << "find saddle " << isJoinSaddle << std::endl;
         // const Propagation* const tmpProp = localPropagation;
         // dynGraph_.print([&, tmpProp](std::size_t t) { return printEdge(t, tmpProp); });

         // Skeleton
         const idVertex upVert = localPropagation->getCurVertex(); // keep before merge
         const idNode upNode = graph_.makeNode(upVert);
         graph_.closeArc(currentArc, upNode);
         std::cout << "close arc " << graph_.printArc(currentArc) << std::endl;

         // Data
         if (isJoinSaddle || isSplitSaddle) {
            // if any saddle, we update the skeleton
            updateReebGraph(lowerComp, upperComp, localPropagation);
         }

         if (isJoinSaddle) {  // && isLast already implied
            graph_.mergeAtSaddle(upNode);
         }

         if (isSplitSaddle) {
            std::cout << "Split things" << std::endl;
         } else {
            // recursive call
           growthFromSeed(upVert, localPropagation);
         }
      }

      template <typename ScalarType>
      std::pair<std::vector<idEdge>, std::vector<idEdge>> FTRGraph<ScalarType>::visitStar(
          const Propagation* const localProp) const
      {
         // TODO re-use the same vectors per thread
         std::vector<idEdge> lowerStar, upperStar;

         const idEdge nbAdjEdges = mesh_->getVertexEdgeNumber(localProp->getCurVertex());
         lowerStar.reserve(nbAdjEdges);
         upperStar.reserve(nbAdjEdges);

         for (idEdge e = 0; e < nbAdjEdges; ++e) {
            idEdge edgeId;
            mesh_->getVertexEdge(localProp->getCurVertex(), e, edgeId);
            idVertex edgeLowerVert, edgeUpperVert;
            std::tie(edgeLowerVert, edgeUpperVert, std::ignore) = getOrderedEdge(edgeId, localProp);
            if (edgeLowerVert == localProp->getCurVertex()) {
               upperStar.emplace_back(edgeId);
            } else {
               lowerStar.emplace_back(edgeId);
            }
         }

         return {lowerStar, upperStar};
      }

      template <typename ScalarType>
      std::set<DynGraphNode<ScalarType>*> FTRGraph<ScalarType>::lowerComps(
          const std::vector<idEdge>& finishingEdges)
      {
         return dynGraph_.findRoot(finishingEdges);
      }

      template <typename ScalarType>
      std::set<DynGraphNode<ScalarType>*> FTRGraph<ScalarType>::upperComps(
          const std::vector<idEdge>& startingEdges)
      {
         return dynGraph_.findRoot(startingEdges);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimage(const Propagation* const localProp)
      {
         const idCell nbAdjTriangles = mesh_->getVertexTriangleNumber(localProp->getCurVertex());

         for (idCell t = 0; t < nbAdjTriangles; ++t) {
            // Classify current cell
            idCell curTriangleid;
            mesh_->getVertexTriangle(localProp->getCurVertex(), t, curTriangleid);

            orderedTriangle   oTriangle  = getOrderedTriangle(curTriangleid, localProp);
            vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localProp);

            std::cout << "v" << localProp->getCurVertex() << " : "
                      << printTriangle(oTriangle, localProp) << std::endl;

            // Update DynGraph
            // We can have an end pos on an unvisited triangle
            // in case of saddle points
            switch (curVertPos) {
               case vertPosInTriangle::Start:
                  updatePreimageStartCell(oTriangle, localProp);
                  break;
               case vertPosInTriangle::Middle:
                  updatePreimageMiddleCell(oTriangle, localProp);
                  break;
               case vertPosInTriangle::End:
                  updatePreimageEndCell(oTriangle, localProp);
                  break;
               default:
                  std::cout << "[FTR]: update preimage error, unknown vertPos type" << std::endl;
                  break;
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageStartCell(const orderedTriangle&   oTriangle,
                                                         const Propagation* const localProp)
      {
         dynGraph_.insertEdge(std::get<0>(oTriangle), std::get<1>(oTriangle), 0);
         // std::cout << "start add edge: " << printEdge(std::get<0>(oTriangle), localProp);
         // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl;
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageMiddleCell(const orderedTriangle&   oTriangle,
                                                          const Propagation* const localProp)
      {
         // Check if exist ?
         // If not, the triangle will be visited again once a merge have occured.
         // So we do not add the edge now
         int t = dynGraph_.removeEdge(std::get<0>(oTriangle), std::get<1>(oTriangle));
         if (t) {
            // std::cout << "mid replace edge: " << printEdge(std::get<0>(oTriangle), localProp);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl;

            dynGraph_.insertEdge(std::get<1>(oTriangle), std::get<2>(oTriangle), 0);
            // std::cout << " with edge: " << printEdge(std::get<1>(oTriangle), localProp);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl;
         } else {
            // std::cout << "mid no found edge: " << printEdge(std::get<0>(oTriangle), localProp);
            // std::cout << " :: " << printEdge(std::get<1>(oTriangle), localProp) << std::endl;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updatePreimageEndCell(const orderedTriangle&   oTriangle,
                                                       const Propagation* const localProp)
      {
         int t = dynGraph_.removeEdge(std::get<1>(oTriangle), std::get<2>(oTriangle));
         if (t) {
            // std::cout << "end remove edge: " << printEdge(std::get<1>(oTriangle), localProp);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl;
         } else {
            // std::cout << "end not found edge: " << printEdge(std::get<1>(oTriangle), localProp);
            // std::cout << " :: " << printEdge(std::get<2>(oTriangle), localProp) << std::endl;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::updateReebGraph(
          const std::set<DynGraphNode<ScalarType>*>& lowerComp,
          const std::set<DynGraphNode<ScalarType>*>& upperComp,
          const Propagation* const                   localPropagation)
      {
         graph_.makeNode(localPropagation->getCurVertex());
         // Use lower comp to recover arcs coming here

         if (lowerComp.size() > 1) {
            std::cout << "join" << std::endl;
         }

         if (upperComp.size() > 1) {
            std::cout << "split" << std::endl;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::localGrowth(Propagation* const localProp)
      {
         const idVertex nbNeigh = mesh_->getVertexNeighborNumber(localProp->getCurVertex());
         for (idVertex n = 0; n < nbNeigh; ++n) {
            idVertex neighId;
            mesh_->getVertexNeighbor(localProp->getCurVertex(), n, neighId);
            if (localProp->compare(localProp->getCurVertex(), neighId)) {
               localProp->addNewVertex(neighId);
            }
         }
      }

      template <typename ScalarType>
      bool FTRGraph<ScalarType>::checkLast(const Propagation* const   localPropagation,
                                           const std::vector<idEdge>& lowerStarEdges)
      {
         // WARNING this vertion only works in sequential
         // TODO: find a way to make this work in parallel
         idVertex isLast = true;
         for (const idEdge e : lowerStarEdges) {
            isLast &= !dynGraph_.isDisconnected(e);
         }
         return isLast;
      }

      /// Tools

      template <typename ScalarType>
      Propagation* FTRGraph<ScalarType>::newPropagation(const idVertex leaf)
      {
         auto compare_fun = [&](idVertex a, idVertex b) { return scalars_->isHigher(a, b); };
         Propagation* localPropagation(new Propagation(leaf, compare_fun));
         const auto propId = propagations_.getNext();
         propagations_[propId] = localPropagation;
         return localPropagation;
      }

      template <typename ScalarType>
      orderedEdge FTRGraph<ScalarType>::getOrderedEdge(const idEdge             edgeId,
                                                       const Propagation* const localProp) const
      {
         idVertex edge0Vert, edge1Vert;
         mesh_->getEdgeVertex(edgeId, 0, edge0Vert);
         mesh_->getEdgeVertex(edgeId, 1, edge1Vert);

         if (localProp->compare(edge0Vert, edge1Vert)) {
            return {edge0Vert, edge1Vert, edgeId};
         } else {
            return {edge1Vert, edge0Vert, edgeId};
         }
      }

      template <typename ScalarType>
      orderedTriangle FTRGraph<ScalarType>::getOrderedTriangle(
          const idCell cellId, const Propagation* const localProp) const
      {
         idEdge edges[3];
         mesh_->getTriangleEdge(cellId, 0, edges[0]);
         mesh_->getTriangleEdge(cellId, 1, edges[1]);
         mesh_->getTriangleEdge(cellId, 2, edges[2]);

         orderedEdge oEdges[3];
         oEdges[0] = getOrderedEdge(edges[0], localProp);
         oEdges[1] = getOrderedEdge(edges[1], localProp);
         oEdges[2] = getOrderedEdge(edges[2], localProp);

         auto compareOEdges = [localProp](const orderedEdge& a, const orderedEdge& b) {
            return localProp->compare(std::get<0>(a), std::get<0>(b)) ||
                   (std::get<0>(a) == std::get<0>(b) &&
                    localProp->compare(std::get<1>(a), std::get<1>(b)));
         };

         std::sort(begin(oEdges), end(oEdges), compareOEdges);

         return {std::get<2>(oEdges[0]), std::get<2>(oEdges[1]), std::get<2>(oEdges[2]), cellId};
      }

      template <typename ScalarType>
      vertPosInTriangle FTRGraph<ScalarType>::getVertPosInTriangle(
          const orderedTriangle& oTriangle, const Propagation* const localProp) const
      {
         orderedEdge firstEdge = getOrderedEdge(std::get<0>(oTriangle), localProp);
         if (std::get<0>(firstEdge) == localProp->getCurVertex()) {
            return vertPosInTriangle::Start;
         } else if (std::get<1>(firstEdge) == localProp->getCurVertex()) {
            return vertPosInTriangle::Middle;
         } else {
            return vertPosInTriangle::End;
         }
      }
   }  // namespace ftr
}  // namespace ttk

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
