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

#include "FTRDataTypes.h"

#include <vector>

namespace ttk {
  namespace ftr {
    class Node {
    private:
      idVertex vertexIdentifier_;
      std::vector<idSuperArc> upArcsIds_;
      std::vector<idSuperArc> downArcsIds_;

    public:
      Node() : vertexIdentifier_(nullVertex) {
      }

      explicit Node(const idVertex vertIdentifier)
        : vertexIdentifier_(vertIdentifier) {
      }

      idVertex getVertexIdentifier(void) const {
        return vertexIdentifier_;
      }

      void setVerterIdentifier(const idVertex v) {
        vertexIdentifier_ = v;
      }

      ftr::NodeType getType() const {
        const valence upVal = upArcsIds_.size();
        const valence downVal = downArcsIds_.size();

        switch(upVal) {
          case 0:
            if(downVal == 1) {
              return ftr::NodeType::Local_maximum;
            }
            break;
          case 1:
            if(downVal == 1) {
              return ftr::NodeType::Regular;
            }
            if(downVal == 2) {
              return ftr::NodeType::Saddle1;
            }
            break;
        }

        switch(downVal) {
          case 0:
            if(upVal == 1) {
              return ftr::NodeType::Local_minimum;
            }
            break;
          case 1:
            if(upVal == 1) {
              return ftr::NodeType::Regular;
            }
            if(upVal == 2) {
              return ftr::NodeType::Saddle2;
            }
            break;
        }

        return ftr::NodeType::Degenerate;
      }

      void reserveUpArc(const idSuperArc nbUpArc) {
        upArcsIds_.reserve(nbUpArc);
      }

      idSuperArc getNbUpArcs(void) const {
        return upArcsIds_.size();
      }

      idSuperArc getUpArc(const idSuperArc index) const {
        return upArcsIds_[index];
      }

      void addUpArc(const idSuperArc arcId) {
        upArcsIds_.emplace_back(arcId);
      }

      void reserveDownArc(const idSuperArc nbDownArcs) {
        downArcsIds_.reserve(nbDownArcs);
      }

      idSuperArc getNbDownArcs(void) const {
        return downArcsIds_.size();
      }

      idSuperArc getDownArc(const idSuperArc index) const {
        return downArcsIds_[index];
      }

      void addDownArc(const idSuperArc arcId) {
        downArcsIds_.emplace_back(arcId);
      }
    };

  } // namespace ftr
} // namespace ttk
