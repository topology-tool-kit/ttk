/// \ingroup base
/// \class ttk::ftr::DataTypes
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-01-19
///
///\brief TTK DataType management for the reeb graph
///
/// \sa ttkFTRGraph.h %for a usage example.

#pragma once

// core includes
#include <DataTypes.h>

// c++ incldues
#include <functional>
#include <limits>

namespace ttk {
  namespace ftr {
    // Types
    // --------

    /// \brief SuperArc index in vect_superArcs_
    using idSuperArc = long unsigned int;
    /// \brief Node index in vect_nodes_
    using idNode = unsigned int;
    /// \brief Vertex index in scalars_
    using idVertex = SimplexId;
    /// \brief Edge index in vect_edgeList_
    using idEdge = SimplexId;
    /// \brief Cell index in vect_cellList_
    using idCell = SimplexId;
    /// \brief for task identifiers
    using idPropagation = idNode;
    /// \brief for vertex up/down valence
    using valence = SimplexId;
    /// \brief retains history
    using idSegmentation = long int;

    /// \brief Edge represented by its 2 vertices, lower then upper
    using orderedEdge = std::tuple<idVertex, idVertex>;

    /// \brief Triangle represented by its 3 edges
    /// Edges are sorted by their starting vertex (see orderedEdge)
    using orderedTriangle = std::tuple<idEdge, idEdge, idEdge>;

    /// \brief position of a vertex in a trianlge
    enum class vertPosInTriangle : char { Start = 0, Middle, End };

    // For tasks:
    // Set using scalar value comparison
    using VertCompFN = std::function<bool(const idVertex, const idVertex)>;
    using EdgeCompFN = std::function<bool(const idEdge, const idEdge)>;
    using linkEdge = std::pair<idEdge, idEdge>;
    using LinkCompFN = std::function<bool(const linkEdge &, const linkEdge &)>;

    // all null values are dealt with using the maximum of the type
    // WHEN C++14 will be used
    // TODO
    // template <typename T>
    // static const T null = std::numeric_limits<T>::max();
    // nullVertex become: null<idVertex>

    // Special values for types
    // --------------------------

    // QUESTION impact on performance using max (0 would be faster alloacted)
    static const idSuperArc nullSuperArc
      = std::numeric_limits<idSuperArc>::max();
    static const idNode nullNode = std::numeric_limits<idNode>::max();
    static const idVertex nullVertex = std::numeric_limits<idVertex>::max();
    static const idEdge nullEdge = std::numeric_limits<idEdge>::max();
    static const idCell nullCell = std::numeric_limits<idCell>::max();
    static const idPropagation nullProp
      = std::numeric_limits<idPropagation>::max();
    static const idSegmentation nullSegment
      = std::numeric_limits<idSegmentation>::max();

    static const linkEdge nullLink = std::make_pair(nullEdge, nullEdge);

    // Enum data
    // ----------

    enum class GraphComponent {
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

    enum class NodeType {
      Local_minimum = 0,
      Saddle1,
      Saddle2,
      Degenerate,
      Local_maximum,
      Regular
    };
  } // namespace ftr
} // namespace ttk
