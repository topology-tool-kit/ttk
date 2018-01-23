#ifndef FTRGRAPH_TEMPLATE_H
#define FTRGRAPH_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk
{
   namespace ftr
   {
      // template functions
      template <class dataType>
      void FTRGraph::build()
      {
         Timer t;

         // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
         if (!mesh_) {
            cerr << "[FTRGraph]: no mesh input" << endl;
            return ;
         }
         if (!inputData_) {
            cerr << "[FTRGraph]: no input data" << endl;
            return ;
         }
         if (!outputData_) {
            cerr << "[FTRGraph]: no ouput data" << endl;
            return ;
         }
#endif

         dataType *outputData = (dataType *)outputData_;
         dataType *inputData  = (dataType *)inputData_;

         int vertexNumber = mesh_->getNumberOfVertices();

         // init the output -- to adapt
         for (int i = 0; i < vertexNumber; i++) {
            outputData[i] = inputData[i];
         }

            // the following open-mp processing is only relevant for embarrassingly
            // parallel algorithms (such as smoothing) -- to adapt
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
         for (int i = 0; i < (int)vertexNumber; i++) {
            // TODO-2
            // processing here!
            // end of TODO-2
         }

         {
            stringstream msg;
            msg << "[FTRGraph] Data-set (" << vertexNumber << " points) processed in "
                << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
            dMsg(cout, msg.str(), timeMsg);
         }
      }
   }
}

#endif /* end of include guard: FTRGRAPH_TEMPLATE_H */
