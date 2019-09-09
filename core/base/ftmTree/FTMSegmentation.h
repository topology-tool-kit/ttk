/// \ingroup base
//
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK classes for the segmentation of the contour tree
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef SEGMENTATION_H_
#define SEGMENTATION_H_

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif
#include <sstream>

#include <list>
#include <tuple>
#include <vector>

#include "FTMDataTypes.h"
#include "FTMStructures.h"

namespace ttk {
  namespace ftm {

    // one segment: like a vector<SimplexId>
    // have a fixed size
    class Segment {
    private:
      std::vector<SimplexId> vertices_;

    public:
      Segment(SimplexId size);

      void sort(const Scalars *s);
      void createFromList(const Scalars *s,
                          std::list<std::vector<SimplexId>> &regularsList,
                          const bool reverse);

      segm_const_it begin(void) const;
      segm_const_it end(void) const;
      segm_it begin(void);
      segm_it end(void);
      SimplexId size(void) const;

      SimplexId operator[](const size_t &idx) const;
      SimplexId &operator[](const size_t &idx);
    };

    // All the segments of the mesh, like a vector<Segment>
    class Segments {
    private:
      std::vector<Segment> segments_;

    public:
      Segments();

      Segments(const Segment &) = delete;

      // add one vertex
      std::tuple<segm_it, segm_it> addLateSimpleSegment(SimplexId v);

      // sort all segment (in parallel?)
      void sortAll(const Scalars *s);

      // callable once
      void resize(const std::vector<SimplexId> &sizes);

      // vector like
      idSegment size(void) const;
      void clear(void);
      Segment &operator[](const size_t &idx);
      const Segment &operator[](const size_t &idx) const;

      // print
      inline std::string print(void) const {
        std::stringstream res;
        res << "{" << std::endl;
        for(const auto &s : segments_) {
          if(s.size()) {
            res << s[0] << " " << s[s.size() - 1] << " : " << s.size()
                << std::endl;
          }
        }
        res << "}" << std::endl;
        return res.str();
      }
    };

    // The segmentation of one arc is a list of segment
    class ArcRegion {
    private:
      // list of segment composing this segmentation and for each segment
      // the begin and the end inside it (as a segment may be subdivided)
      std::list<Region> segmentsIn_;
      // when and how to compact ?
      std::vector<SimplexId> segmentation_;

#ifndef TTK_ENABLE_KAMIKAZE
      // true if the segmentation have been sent to the segmentation_ vector
      bool segmented_;
#endif

    public:
      // create a new arc region with the associated Segment in Segments
      ArcRegion();

      ArcRegion(const segm_it &begin, const segm_it &end);

      // vertex used for the split (not in segmentation anymore), remaining
      // region
      std::tuple<SimplexId, ArcRegion> splitFront(SimplexId v,
                                                  const Scalars *s);

      // vertex used for the split (not in segmentation anymore), remaining
      // region
      std::tuple<SimplexId, ArcRegion> splitBack(SimplexId v, const Scalars *s);

      SimplexId findBelow(SimplexId v,
                          const Scalars *s,
                          const std::vector<idCorresp> &vert2treeOther
                          = std::vector<idCorresp>()) const;

      void concat(const segm_it &begin, const segm_it &end);

      void concat(const ArcRegion &r);

      // if contiguous with current region, merge them
      // else return false
      bool merge(const ArcRegion &r);

      void clear(void) {
        segmentsIn_.clear();
        segmentation_.clear();
      }

      // Put all segments in one vector in the arc
      // Suppose that all segment are already sorted
      // see Segments::sortAll
      // For Contour Tree segmentaion, you have to precise the current arc since
      // a segment can contain vertices for several arcs
      void createSegmentation(const Scalars *s);

      inline SimplexId count(void) const {
        SimplexId res = 0;
        for(const auto &reg : segmentsIn_) {
          res += labs(distance(reg.segmentBegin, reg.segmentEnd));
        }
        return res;
      }

      inline std::string print(void) const {
        std::stringstream res;
        res << "{";
        for(const auto &reg : segmentsIn_) {
          res << " " << *reg.segmentBegin;
          res << "-" << *(reg.segmentEnd - 1);

          // auto it = reg.segmentBegin;
          // for(; it != reg.segmentEnd; ++it){
          //    res << *it << ", ";
          // }
        }
        res << " }";
        return res.str();
      }

      // Direct access to the list of region
      const decltype(segmentsIn_) &getRegions(void) const {
        return segmentsIn_;
      }

      decltype(segmentsIn_) &getRegions(void) {
        return segmentsIn_;
      }

      // vector like manip

      inline SimplexId size(void) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!segmented_)
          std::cerr << "Needs to create segmentation before size" << std::endl;
#endif
        return segmentation_.size();
      }

      SimplexId operator[](SimplexId v) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!segmented_)
          std::cerr
            << "Needs to create segmentation before getting segmentation"
            << std::endl;
#endif
        return segmentation_[v];
      }

      SimplexId &operator[](SimplexId v) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!segmented_)
          std::cerr
            << "Needs to create segmentation before getting segmentation"
            << std::endl;
#endif
        return segmentation_[v];
      }

      decltype(segmentation_)::iterator begin(void) {
        return segmentation_.begin();
      }

      decltype(segmentation_)::iterator end(void) {
        return segmentation_.end();
      }
    };

  } // namespace ftm
} // namespace ttk

#endif /* end of include guard: SEGMENTATION_H_ */
