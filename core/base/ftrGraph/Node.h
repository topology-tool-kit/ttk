/// \ingroup base
/// \class ttk::ftr::Node
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph node of the Graph
///
/// This class manage nodes of the graph structure
///
/// \sa ttk::FTRGraph

#pragma once

#include "DataTypesFTR.h"

#include <vector>

namespace ttk
{
   namespace ftr
   {
      class Node
      {
         private:
          idVertex vertexIdentifier_;
          std::vector<idSuperArc> upArcsIds_;
          std::vector<idSuperArc> downArcsIds_;

         public:
          Node() : vertexIdentifier_(nullVertex)
          {
          }

          explicit Node(const idVertex vertIdentifier) : vertexIdentifier_(vertIdentifier)
          {
          }

          idVertex getVertexIdentifier(void) const
          {
             return vertexIdentifier_;
          }

          void setVerterIdentifier(const idVertex v)
          {
             vertexIdentifier_ = v;
          }

          idSuperArc getNbUpArcs(void) const
          {
             return upArcsIds_.size();
          }

          idSuperArc getUpArc(const idSuperArc index) const
          {
             return upArcsIds_[index];
          }

          idSuperArc getNbDownArcs(void) const
          {
             return downArcsIds_.size();
          }

          idSuperArc getDownArc(const idSuperArc index) const
          {
             return downArcsIds_[index];
          }
      };

   }
}

