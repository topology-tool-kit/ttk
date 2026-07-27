/// \ingroup base
//
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK class representing a SuperArc of a tree,
/// containing regular vertices.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#pragma once

#include <list>
#include <vector>

#include <Debug.h>

#include "FTMDataTypes.h"
#include "FTMSegmentation.h"
#include "FTMStructures.h"

namespace ttk {
  namespace ftm {
    class SuperArc {
    private:
      // Extrema
      idNode downNodeId_, upNodeId_;
      // Stat of this arc (visible, hidden, merged) if merged...
      ComponentState state_;

      // Keep th last vertex seen by this arc
      // After the build a a merge tree, a close step is
      // done, using this field to close each root arc
      SimplexId lastVisited_;

      // Segmentation related
      ArcRegion region_;
      SimplexId verticesSeen_;
      idSuperArc normalizedId_;

    public:
      // -----------------
      // CONSTRUCT
      // -----------------

      // This arc will needs to receive both ends before being printed
      SuperArc()
        : downNodeId_(nullNodes), upNodeId_(nullNodes),
          state_(ComponentState::Visible), lastVisited_(nullVertex),
          verticesSeen_(0), normalizedId_(nullSuperArc) {
      }

      SuperArc(idNode d,
               idNode u,
               const ComponentState &state = ComponentState::Visible)
        : downNodeId_(d), upNodeId_(u), state_(state), lastVisited_(nullVertex),
          verticesSeen_(0), normalizedId_(nullSuperArc) {
      }

      // ------------------
      // ACCESSOR
      // --------------------

      // node

      inline idNode getUpNodeId() const {
        return upNodeId_;
      }

      inline idNode getDownNodeId() const {
        return downNodeId_;
      }

      inline void setUpNodeId(idNode upId) {
        upNodeId_ = upId;
      }

      inline void setDownNodeId(idNode downId) {
        downNodeId_ = downId;
      }

      // last vertex seen, nb vertex seen & ids

      inline SimplexId getLastVisited() const {
        return lastVisited_;
      }

      inline void setLastVisited(SimplexId vertId) {
        lastVisited_ = vertId;
        ++verticesSeen_;
      }

      inline void atomicIncVisited(const SimplexId nb = 1) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
        verticesSeen_ += nb;
      }

      inline void decrNbSeen() {
        --verticesSeen_;
      }

      inline SimplexId getNbVertSeen() const {
        return verticesSeen_;
      }

      inline idSuperArc getNormalizedId() const {
        return normalizedId_;
      }

      inline void setNormalizeIds(const idSuperArc id) {
        normalizedId_ = id;
      }

      // state

      inline bool isHidden() const {
        return state_ == ComponentState::Hidden;
      }

      inline bool isMerged() const {
        return state_ == ComponentState::Merged;
      }

      inline bool isVisible() const {
        return state_ == ComponentState::Visible;
      }

      // ------------
      // Segmentation
      // ------------

      // Function using ArcRegion

      inline void concat(const segm_it &begin, const segm_it &end) {
        region_.concat(begin, end);
      }

      inline void concat(const ArcRegion &r) {
        region_.concat(r);
      }

      inline void concat(const SuperArc &s) {
        region_.concat(s.region_);
      }

      inline void concat(std::tuple<segm_it, segm_it> its) {
        region_.concat(std::get<0>(its), std::get<1>(its));
      }

      // prerequisite for the following segmentation functions
      inline void createSegmentation(const Scalars *s) {
        region_.createSegmentation(s);
      }

      // Direct read access to the list of region
      const std::list<Region> &getRegions() const {
        return region_.getRegions();
      }

      std::list<Region> &getRegions() {
        return region_.getRegions();
      }

      const ArcRegion &getRegion() const {
        return region_;
      }

      size_t regionSize() const {
        return region_.count();
      }

      void clearSegmentation() {
        region_.clear();
      }

      // access segmentation (after createSegmentation)
      // vector-like

      inline size_t size() const {
        return region_.size();
      }

      std::vector<SimplexId>::iterator begin() {
        return region_.begin();
      }

      std::vector<SimplexId>::iterator end() {
        return region_.end();
      }

      SimplexId operator[](SimplexId v) const {
        return region_[v];
      }

      SimplexId &operator[](SimplexId v) {
        return region_[v];
      }

      // Access Segmentation legacy

      SimplexId getNumberOfRegularNodes() const {
        return region_.size();
      }

      SimplexId getRegularNodeId(SimplexId id) const {
        return region_[id];
      }

      // process segmentation

      // keep the front segmentation, return the back
      std::tuple<SimplexId, ArcRegion> splitFront(SimplexId v,
                                                  const Scalars *s) {
        return region_.splitFront(v, s);
      }

      // Keep the back, return the front
      std::tuple<SimplexId, ArcRegion> splitBack(SimplexId v,
                                                 const Scalars *s) {
        return region_.splitBack(v, s);
      }

      SimplexId findBelow(SimplexId v,
                          const Scalars *s,
                          const std::vector<idCorresp> &vert2treeOther
                          = std::vector<idCorresp>()) const {
        return region_.findBelow(v, s, vert2treeOther);
      }

      bool merge(const SuperArc &s) {
        return region_.merge(s.region_);
      }

      std::string printReg() const {
        return region_.print();
      }
    };

  } // namespace ftm
} // namespace ttk
