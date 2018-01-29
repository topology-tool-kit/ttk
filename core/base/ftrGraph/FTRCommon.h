/// \ingroup base
//
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK structures for the reeb graph

#pragma once

#include "DataTypesFTR.h"

#include <iostream>
#include <vector>

namespace ttk
{
   namespace ftr
   {
      // Compute parameters (global)
      struct Params {
         bool segm        = true;
         bool normalize   = true;
         bool advStats    = true;
         int  samplingLvl = 0;

         int threadNumber_ = 1;
         int debugLevel_   = 1;

         void printSelf()
         {
            if (debugLevel_) {
               std::cout << "[FTR Graph]: thread number: " << threadNumber_ << std::endl;
               std::cout << "[FTR Graph]: debug lvl: " << debugLevel_ << std::endl;
               if (debugLevel_ > 2) {
                  std::cout << "[FTR Graph]: segmentation: " << std::boolalpha << segm << std::endl;
                  std::cout << "[FTR Graph]: sampling level: " << samplingLvl << std::endl;
               }
            }
         }
      };

      /// Force the class to use functions alloc and init
      /// alloc is used for all system allocation
      /// init is used to fill arrays that needs to be
      class Allocable
      {
        protected:
         /// Allocation may depends on the number of vertices
         idVertex nbVerts_ = nullVertex;

        public:
         void setNumberOfVertices(const idVertex nbVerts)
         {
            nbVerts_ = nbVerts;
         }

         template<typename type>
         void fillVector(std::vector<type>& vect, const type& elmt)
         {
            const std::size_t nbIt = vect.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t i = 0; i < nbIt; i++) {
               vect[i] = elmt;
            }
         }

         virtual void alloc() = 0;

         virtual void init() = 0;
      };
   }
}

