/// \ingroup base
/// \class ttk::ftr::SuperArc
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph arc
///
/// This class manage  arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#pragma once

#include "DataTypesFTR.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include<iostream>
#endif

namespace ttk
{
   namespace ftr
   {
      class Node;

      class SuperArc
      {
        private:
         idNode upNodeId_;
         idNode downNodeId_;
         // TODO Charles: segmentation

        public:
         SuperArc() : upNodeId_(nullNode), downNodeId_(nullNode)
         {
         }

         SuperArc(const idNode down, const idNode up) : upNodeId_(up), downNodeId_(down)
         {
         }

         idNode getUpNodeId(void) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if(upNodeId_ == nullNode)
            {
               std::cerr << "[FTR Graph]: Arc have null up node" << std::endl;
            }
#endif
            return upNodeId_;
         }

         idNode getDownNodeId(void) const
         {
#ifndef TTK_ENABLE_KAMIKAZE
            if(downNodeId_ == nullNode)
            {
               std::cerr << "[FTR Graph]: Arc have null down node" << std::endl;
            }
#endif
            return downNodeId_;
         }

         void setUpNodeId(const idNode id)
         {
            upNodeId_ = id;
         }

         void seDownNodeId(const idNode id)
         {
            downNodeId_ = id;
         }
      };
   }
}

