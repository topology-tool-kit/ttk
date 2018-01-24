/// \ingroup base
//
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-22
///
///\brief TTK structures for the reeb graph

#pragma once

#include<iostream>

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

      void printSelf() {
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
}
}

