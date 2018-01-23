#ifndef FTRGRAPH_TEMPLATE_H
#define FTRGRAPH_TEMPLATE_H

#include "FTRGraph.h"

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
         // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
         if (!mesh_) {
            cerr << "[FTR Graph]: no mesh input" << endl;
            return;
         }
#endif
#ifdef TTK_ENABLE_OPENMP
         omp_set_num_threads(threadNumber_);
         omp_set_nested(1);
#endif

         params_->printSelf();

         int vertexNumber = mesh_->getNumberOfVertices();
         scalars_->setSize(vertexNumber);

         {
            stringstream msg;
            msg << "[FTR Graph] Data-set (" << vertexNumber << " points) processed in "
                << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
         }
      }
   }
}

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
