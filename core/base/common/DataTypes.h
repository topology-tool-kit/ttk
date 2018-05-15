/// \ingroup base
/// \class ttk::DataTypes
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date May 2018.
///
///\brief TTK base package defining the standard types.

#ifndef DATATYPES_H
#define DATATYPES_H

namespace ttk{
  /// \brief Identifier type for vertices.
  using VertexId = long long int;

  /// \brief Identifier type for edges.
  using EdgeId = long long int;

  /// \brief Identifier type for triangles (3D only).
  using TriangleId = long long int;

  /// \brief Identifier type for cells.
  /// Here the notion of cell refers to the simplices of maximal
  /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
  using CellId = long long int;

  /// \brief Identifier type for threads (i.e. with OpenMP).
  using ThreadId = int;

  /// \brief Identifier type for tasks (i.e. with OpenMP).
  using TaskId = int;
}

#endif // DATATYPES_H
