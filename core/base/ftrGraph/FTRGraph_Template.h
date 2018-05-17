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
          : params_(new Params),
            scalars_(new Scalars<ScalarType>),
            needDelete_(true),
            mesh_(nullptr)
      {
      }

      template <typename ScalarType>
      FTRGraph<ScalarType>::FTRGraph(Params* const params, Triangulation* mesh,
                                     Scalars<ScalarType>* const scalars)
          : params_(params), scalars_(scalars), needDelete_(false), mesh_(mesh)
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

         // set size for allocation / init
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
         printTime(timeBuild, "[FTR Graph]: build time: ", timeMsg);

         // post-process
         graph_.arcs2nodes();

         // Debug print
         printGraph(3);

         // Message user
         {
            std::stringstream msg;
            msg << "[FTR Graph] Data-set (" << vertexNumber << " points) processed in "
                << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << std::endl;
            dMsg(std::cout, msg.str(), timeMsg);
         }
      }

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
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::alloc()
      {
         scalars_->alloc();
         graph_.alloc();
         dynGraph_.alloc();

         propagations_.reserve(scalars_->getSize());
         toVisit_.resize(scalars_->getSize());
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::init()
      {
         scalars_->removeNaN();
         scalars_->init();
         graph_.init();
         dynGraph_.init();

         fillVector<UnionFind*>(toVisit_, nullptr);
      }

   }  // namespace ftr
}  // namespace ttk

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
