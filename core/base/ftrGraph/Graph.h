/// \ingroup base
/// \class ttk::ftr::Graph
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph skeleton
///
/// This class manage nodes and arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#pragma once

#include "DataTypes.h"
#include "FTRCommon.h"

namespace ttk
{
   namespace ftr
   {
      class Graph : public Allocable
      {
        public:
         Graph();
         virtual ~Graph();

         // Accessor on structure
         // ---------------------

         // Initialize functions
         // --------------------

         void alloc() override;

         void init() override;

      };
   }
}

