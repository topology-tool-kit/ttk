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

         params_->printSelf();

         // Precompute

         DebugTimer timeAlloc;
         alloc();
         printTime(timeAlloc, "[FTR Graph]: alloc time: ", infoMsg);

         DebugTimer timeNaN;
         scalars_->removeNaN();
         printTime(timeNaN, "[FTR Graph]: remove NaN time: ", infoMsg);

         DebugTimer timeInit;
         init();
         printTime(timeInit, "[FTR Graph]: init time: ", infoMsg);

         DebugTimer timeSort;
         scalars_->sort();
         printTime(timeSort, "[FTR Graph]: sort time: ", infoMsg);

         // Build the graph

         DebugTimer timeLeafSearch;
         leafSearch();
         printTime(timeLeafSearch, "[FTR Graph]: leaf search time: ", timeMsg);

         DebugTimer timeSwipe;
         swipe();
         printTime(timeSwipe, "[FTR Graph]: swipe time: ", timeMsg);

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

         cout << " grain size " << get<0>(leafChunk) << " nb Tasks " << get<1>(leafChunk) << endl;

         for (idTask leafChunkId = 0; leafChunkId < get<1>(leafChunk); ++leafChunkId) {
            const idVertex lowerBound = Tasks::getBegin(leafChunkId, get<0>(leafChunk));
            const idVertex upperBound = Tasks::getEnd(leafChunkId, get<0>(leafChunk), scalars_->getSize());
            cout << "lower " << lowerBound << endl;
            cout << "upper " << upperBound << endl;
            cout << " " << endl;
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::swipe()
      {
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
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::init()
      {
         scalars_->init();
         graph_.init();
      }
   }
}

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
