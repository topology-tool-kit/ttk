/// \ingroup baseCode
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

#ifndef withKamikaze
#include <iostream>
#endif
#include <sstream>

#include <list>
#include <tuple>
#include <vector>

#include "DataTypes.h"
#include "Structures.h"

using namespace std;

namespace ttk
{
namespace ftm
{

   // one segment: like a vector<idVertex>
   // have a fixed size
   class Segment
   {
     private:
      vector<idVertex> vertices_;

     public:
      Segment(idVertex size);

      void sort(const Scalars* s);
      void createFromList(const Scalars* s, list<vector<idVertex>>& regularsList, const bool reverse);

      segm_const_it begin(void) const;
      segm_const_it end(void) const;
      segm_it       begin(void);
      segm_it       end(void);
      idVertex      size(void) const;

      idVertex operator[](const size_t& idx) const;
      idVertex& operator[](const size_t& idx);
   };

   // All the segments of the mesh, like a vector<Segment>
   class Segments
   {
     private:
      vector<Segment> segments_;

     public:
      Segments();

      Segments(const Segment&) = delete;

      // add one vertex
      tuple<segm_it, segm_it> addLateSimpleSegment(idVertex v);

      // sort all segment (in parallel?)
      void sortAll(const Scalars* s);

      // callable once
      void resize(const vector<idVertex>& sizes);

      // vector like
      idSegment size(void) const;
      void clear(void);
      Segment& operator[](const size_t& idx);
      const Segment& operator[](const size_t& idx) const;

      // print
      inline string print(void) const
      {
         stringstream res;
         res << "{" << endl;
         for (const auto& s : segments_) {
            if (s.size()) {
               res << s[0] << " " << s[s.size() - 1] << " : " << s.size() << endl;
            }
         }
         res << "}" << endl;
         return res.str();
      }

   };

   // The segmentation of one arc is a list of segment
   class ArcRegion
   {
     private:
      // list of segment composing this segmentation and for each segment
      // the begin and the end inside it (as a segment may be subdivided)
      list<Region> segmentsIn_;
      // when and how to compact ?
      vector<idVertex> segmentation_;

#ifndef withKamikaze
      // true if the segmentation have been sent to the segmentation_ vector
      bool segmented_;
#endif

     public:
      // create a new arc region with the associated Segment in Segments
      ArcRegion();

      ArcRegion(const segm_it& begin, const segm_it& end);

      // vertex used for the split (not in segmentation anymore), remaining region
      tuple<idVertex, ArcRegion> splitFront(idVertex v, const Scalars* s);

      // vertex used for the split (not in segmentation anymore), remaining region
      tuple<idVertex, ArcRegion> splitBack(idVertex v, const Scalars* s);

      idVertex findBelow(idVertex v, const Scalars* s,
                         const vector<idCorresp>& vert2treeOther = vector<idCorresp>()) const;

      void concat(const segm_it& begin, const segm_it& end);

      void concat(const ArcRegion& r);

      // if contiguous with current region, merge them
      // else return false
      bool merge(const ArcRegion& r);

      void clear(void)
      {
         segmentsIn_.clear();
         segmentation_.clear();
      }

      // Put all segments in one vector in the arc
      // Suppose that all segment are already sorted
      // see Segments::sortAll
      // For Contour Tree segmentaion, you have to precise the current arc since
      // a segment can contain vertices for several arcs
      void createSegmentation(const Scalars* s);

      inline idVertex count(void) const
      {
         size_t res = 0;
         for (const auto& reg : segmentsIn_) {
            res += abs(distance(reg.segmentBegin, reg.segmentEnd));
         }
         return res;
      }

      inline string print(void) const
      {
         stringstream res;
         res << "{";
         for (const auto& reg : segmentsIn_) {
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
      const decltype(segmentsIn_) & getRegions(void) const
      {
         return segmentsIn_;
      }

      decltype(segmentsIn_) & getRegions(void)
      {
         return segmentsIn_;
      }

      // vector like manip

      inline idVertex size(void) const
      {
#ifndef withKamikaze
         if (!segmented_)
            std::cerr << "Needs to create segmentation before size" << std::endl;
#endif
         return segmentation_.size();
      }

      idVertex operator[](idVertex v) const
      {
#ifndef withKamikaze
         if (!segmented_)
            cerr << "Needs to create segmentation before getting segmentation" << endl;
#endif
         return segmentation_[v];
      }

      idVertex& operator[](idVertex v)
      {
#ifndef withKamikaze
         if (!segmented_)
            cerr << "Needs to create segmentation before getting segmentation" << endl;
#endif
         return segmentation_[v];
      }

      decltype(segmentation_)::iterator begin(void)
      {
         return segmentation_.begin();
      }

      decltype(segmentation_)::iterator end(void)
      {
         return segmentation_.end();
      }
   };

}
}

#endif /* end of include guard: SEGMENTATION_H_ */
