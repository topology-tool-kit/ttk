#ifndef FTRGRAPH_TEMPLATE_H
#define FTRGRAPH_TEMPLATE_H

#include "FTRGraph.h"
#include "Tasks.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      FTRGraph<ScalarType>::FTRGraph()
          : params_(new Params),
            mesh_(nullptr),
            scalars_(new Scalars<ScalarType>),
            needDelete_(true)
      {
      }

      template <typename ScalarType>
      FTRGraph<ScalarType>::FTRGraph(Params* const params, Triangulation* mesh,
                                     Scalars<ScalarType>* const scalars)
          : params_(params), mesh_(mesh), scalars_(scalars), needDelete_(false)
      {
      }

      template<typename ScalarType>
      FTRGraph<ScalarType>::~FTRGraph()
      {
         if(needDelete_) {
            delete params_;
            delete scalars_;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::build()
      {
         Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
         if (!mesh_) {
            cerr << "[FTR Graph]: no mesh input" << endl;
            return;
         }
         if (!scalars_) {
            cerr << "[FTR Graph]: no scalars given" << endl;
            return;
         }
         if (!params_) {
            cerr << "[FTR Graph]: no parameters" << endl;
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
         dynGraph_.setNumberOfVertices(vertexNumber);

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
#pragma omp parallel num_threads(params_->threadNumber_)
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
               swipe();
               printTime(timeSwipe, "[FTR Graph]: swipe time: ", timeMsg);
            }
         }

         // Message user
         {
            stringstream msg;
            msg << "[FTR Graph] Data-set (" << vertexNumber << " points) processed in "
                << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::leafSearch()
      {
         TaskChunk leafChunkParams(scalars_->getSize());
         leafChunkParams.grainSize = 10000;
         auto leafChunk            = Tasks::getChunk(leafChunkParams);

         for (idTask leafChunkId = 0; leafChunkId < get<1>(leafChunk); ++leafChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(leafChunkId)
#endif
            {
               const idVertex lowerBound = Tasks::getBegin(leafChunkId, get<0>(leafChunk));
               const idVertex upperBound = Tasks::getEnd(leafChunkId, get<0>(leafChunk), scalars_->getSize());

               for(idVertex v = lowerBound; v < upperBound ; ++v) {
                  const valence vNeighNumber = mesh_->getVertexNeighborNumber(v);
                  bool isMinV = true;

                  for(valence n = 0; n < vNeighNumber; ++n) {
                     idVertex neigh;
                     mesh_->getVertexNeighbor(v, n, neigh);

                     if (scalars_->isLower(neigh, v)){
                        isMinV = false;
                        break;
                     }
                  }

                  // v is a minimum, add it to the leaves
                  if(isMinV) {
                     graph_.addLeaf(v);
                     graph_.makeNode(v);
                  }
               }
            } // end task
         }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
         cout << "find: " << graph_.getNumberOfLeaves() << " leaves" << endl;
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::swipe()
      {
         // TODO
         // DynGraph
         // Propagation FH
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::printTime(DebugTimer& timer, const string& msg, const int lvl) const
      {
         ostringstream outString;
         outString << msg << timer.getElapsedTime()  << endl;
         dMsg(cout, outString.str(), lvl);
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::alloc()
      {
         scalars_->alloc();
         graph_.alloc();
         dynGraph_.alloc();
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::init()
      {
         scalars_->removeNaN();
         scalars_->init();
         graph_.init();
         dynGraph_.init();
      }
   }
}

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
