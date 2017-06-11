/// \ingroup baseCode
//
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK classe representing a SuperArc of a tree,
/// containing regular vertices.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef SUPERARC_H
#define SUPERARC_H

#include <list>
#include <vector>

#include <Debug.h>

#include "DataTypes.h"
#include "Segmentation.h"
#include "Structures.h"

namespace ttk
{
   class SuperArc
   {
     private:
      // Extrema
      idNode downNodeId_, upNodeId_;
      // Stat of this arc (visible, hidden, merged) if merged...
      ComponentState state_;
      // ... use this field to know the replacant arc.
      // Caution, we do not want chained replacant !
      idSuperArc substituteArcId_;

      // Keep th last vertex seen by this arc
      // After the build a a merge tree, a close step is
      // done, using this field to close each root arc
      idVertex lastVisited_;

      // Segmentation related
      ArcRegion region_;
      idVertex  verticesSeen_;

     public:
      // CONSTRUCT
      // -----------------
      // {

      // This arc will needs to receive both ends before being printed
      SuperArc()
          : downNodeId_(nullNodes),
            upNodeId_(nullNodes),
            state_(ComponentState::Visible),
            substituteArcId_(nullSuperArc),
            lastVisited_(nullVertex),
            region_(),
            verticesSeen_(0)
      {
      }

      SuperArc(idNode d, idNode u, const ComponentState &state = ComponentState::Visible)
          : downNodeId_(d),
            upNodeId_(u),
            state_(state),
            substituteArcId_(nullSuperArc),
            lastVisited_(nullVertex),
            region_(),
            verticesSeen_(0)
      {
      }

      // }
      // ------------------
      // ACCESSOR
      // --------------------
      // {
      //
      // node
      // .................................{

      inline idNode getUpNodeId(void) const
      {
         return upNodeId_;
      }

      inline idNode getDownNodeId(void) const
      {
         return downNodeId_;
      }

      inline void setUpNodeId(idNode upId)
      {
         upNodeId_ = upId;
      }

      inline void setDownNodeId(idNode downId)
      {
         downNodeId_ = downId;
      }

      // }
      // last vertex seen & nb vertex seen
      // .................................{

      inline idVertex getLastVisited(void) const
      {
         return lastVisited_;
      }

      inline void setLastVisited(idVertex vertId)
      {
         lastVisited_ = vertId;
         ++verticesSeen_;
      }

      inline void atomicIncVisited(const idVertex nb)
      {
#pragma omp atomic update
         verticesSeen_ += nb;
      }

      inline void decrNbSeen(void)
      {
         --verticesSeen_;
      }

      inline idVertex getNbVertSeen(void) const
      {
         return verticesSeen_;
      }

      // }
      // state
      // .................................{

      inline bool isHidden(void) const
      {
         return state_ == ComponentState::Hidden;
      }

      inline bool isMerged(void) const
      {
         return state_ == ComponentState::Merged;
      }

      inline bool isVisible(void) const
      {
         return state_ == ComponentState::Visible;
      }

      inline void merge(idSuperArc arc)
      {
         substituteArcId_ = arc;
         state_           = ComponentState::Merged;
      }

      // }
      // replacant arc/tree (merge)
      // .................................{

      inline idSuperArc getReplacantArcId(void) const
      {
         return substituteArcId_;
      }

      // }
      // -----------
      // Segmentation
      // ----------
      // {
      // Fonction using ArcRegion

      inline void concat(const segm_it &begin, const segm_it &end)
      {
         region_.concat(begin, end);
      }

      inline void concat(const ArcRegion &r)
      {
         region_.concat(r);
      }

      inline void concat(const SuperArc &s)
      {
         region_.concat(s.region_);
      }

      inline void concat(tuple<segm_it, segm_it> its)
      {
         region_.concat(get<0>(its), get<1>(its));
      }

      // prerequisite for the following segmentation functions
      inline void createSegmentation(const Scalars *s)
      {
         region_.createSegmentation(s);
      }

      // Direct read access to the list of region
      const list<Region> &getRegions(void) const
      {
         return region_.getRegions();
      }

      list<Region> &getRegions(void)
      {
         return region_.getRegions();
      }

      const ArcRegion &getRegion(void) const
      {
         return region_;
      }

      size_t regionSize(void) const
      {
         return region_.count();
      }

      void clearSegmentation(void)
      {
         region_.clear();
      }

      // access segmentation (after createSegmentation)
      // vector-like

      inline size_t size(void) const
      {
         return region_.size();
      }

      vector<idVertex>::iterator begin(void)
      {
         return region_.begin();
      }

      vector<idVertex>::iterator end(void)
      {
         return region_.end();
      }

      idVertex operator[](idVertex v) const
      {
         return region_[v];
      }

      idVertex &operator[](idVertex v)
      {
         return region_[v];
      }

      // Access Segmentation legacy

      idVertex getNumberOfRegularNodes() const
      {
         return region_.size();
      }

      idVertex getRegularNodeId(idVertex id) const
      {
         return region_[id];
      }

      // process segmentation

      // keep the front segmentation, return the back
      tuple<idVertex, ArcRegion> splitFront(idVertex v, const Scalars *s)
      {
         return region_.splitFront(v, s);
      }

      // Keep the back, return the front
      tuple<idVertex, ArcRegion> splitBack(idVertex v, const Scalars *s)
      {
         return region_.splitBack(v, s);
      }

      idVertex findBelow(idVertex v, const Scalars *s,
                         const vector<idCorresp> &vert2treeOther = vector<idCorresp>()) const
      {
         return region_.findBelow(v, s, vert2treeOther);
      }

      bool merge(const SuperArc &s)
      {
         return region_.merge(s.region_);
      }

      string printReg(void) const
      {
         return region_.print();
      }
      // }
   };
}

#endif /* end of include guard: SUPERARC_H */
