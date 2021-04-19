/// \ingroup base
//
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK structures for the contour tree
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iterator>

#include "DeprecatedDataTypes.h"
#include "DeprecatedNode.h"
#include "DeprecatedSuperArc.h"

namespace ttk {
  namespace cf {
    // Compute parameters (global)
    struct Params {
      int debugLevel;
      TreeType treeType;
      SimplifMethod simplifyMethod;
      double simplifyThreshold;
    };

    // Scalar related containers (global)
    struct Scalars {
      SimplexId size{};
      void *values{};
      const SimplexId *sosOffsets{};
      std::vector<SimplexId> sortedVertices{};

      bool isLower(const SimplexId &a, const SimplexId &b) const {
        return sosOffsets[a] < sosOffsets[b];
      }
      bool isHigher(const SimplexId &a, const SimplexId &b) const {
        return sosOffsets[a] > sosOffsets[b];
      }
    };

    // Tree datas ( 1 per tree )
    struct TreeData {
      TreeType treeType;
      idPartition partition;

      // components : tree / nodes / extrema
      std::vector<SuperArc> superArcs;
      std::vector<Node> nodes;
      std::vector<idNode> leaves, roots;

      // arc crossing an interface (one can be in both)
      std::vector<idSuperArc> arcsCrossingBelow, arcsCrossingAbove;

      // vertex 2 node / superarc
      std::vector<idCorresp> vert2tree;
    };

    // info on one vertex and CT arc in wich it is
    struct vertex {
      SimplexId id;
      idSuperArc ctArc;
    };

    using segmentIterator = std::vector<vertex>::iterator;
    using segmentRevIterator = std::vector<vertex>::reverse_iterator;

    // If we want to cross a Segment in the sorted order,
    // wich is form the end to the beginning in the case of Split Tree,
    // we can do so by using sbegin and send which use this class
    class sorted_iterator : public segmentIterator {
    public:
      sorted_iterator(const segmentIterator &base)
        : segmentIterator(base), forward_(true) {
      }
      sorted_iterator(const segmentRevIterator &base)
        : segmentIterator(base.base()), forward_(false) {
      }

      const sorted_iterator &operator++() {
        if(forward_) {
          std::vector<vertex>::iterator::operator++();
        } else {
          std::vector<vertex>::iterator::operator--();
        }

        return *this;
      }

    private:
      bool forward_;
    };

    // Segmentation data
    struct Region {
      // inverted in case of split tree
      segmentIterator segmentBegin;
      segmentIterator segmentEnd;
      bool forward;

      sorted_iterator sbegin() const {
        if(forward)
          return sorted_iterator(segmentBegin);
        return sorted_iterator(segmentRevIterator(segmentEnd));
      }

      sorted_iterator send() const {
        if(forward)
          return sorted_iterator(segmentEnd);
        return sorted_iterator(segmentRevIterator(segmentBegin));
      }
    };
  } // namespace cf
} // namespace ttk

#endif /* end of include guard: STRUCTURES_H */
