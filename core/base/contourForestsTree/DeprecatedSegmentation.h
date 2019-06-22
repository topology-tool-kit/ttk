/// \ingroup base
//
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK classes for the segmentation of the contour tree
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef SEGMENTATION_H_
#define SEGMENTATION_H_

#include <forward_list>
#include <vector>

#include "DeprecatedDataTypes.h"
#include "DeprecatedStructures.h"

namespace ttk {
  namespace cf {

    // one segment: like a std::vector<vertex>
    class Segment {
    public:
      Segment(const bool order = true);

      bool isAscending(void) const;

      void sort(const Scalars *s);

      // std::vector like
      void emplace_back(const SimplexId &v);
      void clear(void);
      void norm_next(segmentIterator &it);

      SimplexId &operator[](size_t idx);
      const SimplexId &operator[](size_t idx) const;

      // custom iterator to cross the segment in sorted order
      sorted_iterator sbegin(void);
      sorted_iterator send(void);

    private:
      std::vector<vertex> vertices_;
      bool ascendingOrder_;
    };

    // All the segments of the mesh, like a std::vector<Segment>
    class Segments {
    public:
      idSegment size(void) const;
      void clear(void);
      Segment &operator[](size_t idx);
      const Segment &operator[](size_t idx) const;

    private:
      std::vector<Segment> segments_;
    };

    // The segmentation of one arc is a list of segment
    class ArcRegion {
    public:
      // create a new arc region with the associated Segment in Segments
      ArcRegion(const segmentIterator &s);

      // During combinaison: concat some segmentations
      void addSegment(const segmentIterator &begin, const segmentIterator &end);

      // Put all segments in one std::vector in the arc
      void createSegmentation(const idSuperArc &thisArc);

    private:
      // list of segment composing this segmentation and for each segment
      // the begin and the end inside it (as a segment may be subdivided)
      std::forward_list<Region> segmentsIn_;
      // when and how to compact ?
      std::vector<vertex> segmentation_; // SimplexId only ?
    };
  } // namespace cf
} // namespace ttk

#endif /* end of include guard: SEGMENTATION_H_ */
