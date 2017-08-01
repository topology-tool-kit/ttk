/// \ingroup baseCode
/// \class ttk::DataTypes
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief Define usefull type for TTK (vertices, edges, cells)
//and their null values.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef DATATYPES_H
#define DATATYPES_H

// TODO:---------------------------------------------------------------------------
// This file is the same and have the same sentinel that the ContourForests one
// This is an undesirable situation that can leads to bug if one of this file
// is modified and not the other.
// FIX: We could try to instore a global Type Management file.
// --------------------------------------------------------------------------------

#include <tuple>
#include <limits>

namespace ttk
{

   // Types
   // --------

   /// \brief SuperArc index in vect_superArcs_
   using idSuperArc = long unsigned int;
   /// \brief Node index in vect_nodes_
   using idNode = unsigned int;
   /// \brief Vertex index in scalars_
   using idVertex = int;
   /// \brief Edge index in vect_edgeList_
   using idEdge = int;
   /// \brief Cell index in vect_cellList_
   using idCell = int;

   /// \brief type used to recover Node/Arc in vert2tree SIGNED ONLY
   // Warning, in long long int the max super arc is -1, might not be able to deal with
   // too large data
   using idCorresp = long long int;

   /// \brief for the segmentation, we have an array of segment containing area of the mesh
   using idSegment = idSuperArc;

   /// \brief type use to store threads related numbers
   using numThread = unsigned char;
   /// \brief index of the interface/partition in vect_interfaces_
   using idInterface = numThread;
   using idPartition = numThread;

   /// \brief type stored by UnionFind
   using ufDataType = long int;

   // Special values for types
   // --------------------------

   // QUESTION impact on performance using max (0 would be faster alloacted)
   static const idSuperArc     nullSuperArc       = std::numeric_limits<idSuperArc>::max();
   static const idNode         nullNodes          = std::numeric_limits<idNode>::max();
   static const idVertex       nullVertex         = std::numeric_limits<idVertex>::max();
   static const idCorresp      nullCorresp        = std::numeric_limits<idCorresp>::max();
   static const idSegment      nullSegment        = std::numeric_limits<idSegment>::max();
   static const idInterface    nullInterface      = std::numeric_limits<idInterface>::max();
   static const idPartition    nullPartition      = std::numeric_limits<idPartition>::max();
   static const ufDataType     nullUfData         = std::numeric_limits<ufDataType>::max();
   static constexpr ufDataType specialUfData      = std::numeric_limits<ufDataType>::max() - 1;

   // Enum data
   // ----------

   enum TreeType : char { Join = 0, Split = 1, Contour = 2 , JoinAndSplit = 3};

   enum SimplifMethod : char { Persist = 0, Span = 1, NbVert = 2, NbArc = 3 };

   enum ComponentState : char { Visible, Hidden, Pruned, Merged };

   enum class TreeComponent { Arc = -1, Local_minimum, Saddle1, Saddle2, Local_maximum };

   enum class ArcType { Min_arc = 0, Max_arc, Saddle1_arc, Saddle2_arc, Saddle1_saddle2_arc };

   enum class NodeType { Local_minimum = 0, Saddle1, Saddle2, Local_maximum, Regular, Degenerate };

}

#endif /* end of include guard: DATATYPES_H */
