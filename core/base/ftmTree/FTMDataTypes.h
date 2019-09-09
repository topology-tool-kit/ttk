/// \ingroup base
/// \class ttk::ftm::DataTypes
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkContourForests.cpp %for a usage example.

#ifndef DATATYPES_FTM_H
#define DATATYPES_FTM_H

#include <DataTypes.h>
#include <functional>
#include <limits>
#include <set>
#include <tuple>

namespace ttk {
  namespace ftm {
    // Types
    // --------

    /// \brief SuperArc index in vect_superArcs_
    using idSuperArc = long unsigned int;
    /// \brief Node index in vect_nodes_
    using idNode = unsigned int;

    /// \brief type used to recover Node/Arc in vert2tree SIGNED ONLY
    // Warning, in long long int the max super arc is -1, might not be able to
    // deal with too large data
    using idCorresp = long long int;

    /// \brief for the segmentation, we have an array of segment containing area
    /// of the mesh
    using idSegment = idSuperArc;

    /// \brief type use to store threads related numbers
    using numThread = ThreadId;

    /// \brief type stored by UnionFind
    using ufDataType = long int;

    /// \brief manage number of threads
    using idThread = ThreadId;

    /// \brief for task identifiers
    using idTask = TaskId;

    /// \brief for vertex up/down valence
    using valence = SimplexId;

    // For tasks:
    // Set using scalar value comparison
    using VertCompFN = std::function<bool(SimplexId, SimplexId)>;
    using SetPropagation = std::set<SimplexId, VertCompFN>;
    using SetCompFN
      = std::function<bool(const SetPropagation &, const SetPropagation &)>;

    // Special values for types
    // --------------------------

    // QUESTION impact on performance using max (0 would be faster alloacted)
    static const idSuperArc nullSuperArc
      = std::numeric_limits<idSuperArc>::max();
    static const idNode nullNodes = std::numeric_limits<idNode>::max();
    static const SimplexId nullVertex = std::numeric_limits<SimplexId>::max();
    static const SimplexId nullEdge = std::numeric_limits<SimplexId>::max();
    static const SimplexId nullCell = std::numeric_limits<SimplexId>::max();
    static const idCorresp nullCorresp = std::numeric_limits<idCorresp>::max();
    static const idSegment nullSegment = std::numeric_limits<idSegment>::max();
    static const ufDataType nullUfData = std::numeric_limits<ufDataType>::max();
    static constexpr ufDataType specialUfData
      = std::numeric_limits<ufDataType>::max() - 1;
    static const idThread nullThread = std::numeric_limits<idThread>::max();

    // Enum data
    // ----------

    enum TreeType : char { Join = 0, Split = 1, Contour = 2, Join_Split = 3 };

    enum SimplifMethod : char { Persist = 0, Span = 1, NbVert = 2, NbArc = 3 };

    enum ComponentState : char { Visible, Hidden, Pruned, Merged };

    enum class TreeComponent {
      Arc = -1,
      Local_minimum,
      Saddle1,
      Saddle2,
      Local_maximum
    };

    enum class ArcType : char {
      Min_arc = 0,
      Max_arc,
      Saddle1_arc,
      Saddle2_arc,
      Saddle1_saddle2_arc
    };
  } // namespace ftm
} // namespace ttk

// Tests
// Stats time are not compatible with CT, only MT
// #define withStatsTime 1

// #define withProcessSpeed 1

#endif /* end of include guard: DATATYPES_H */
