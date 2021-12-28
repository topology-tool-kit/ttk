/// \ingroup base
/// \class ttk::ZeroSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief %ZeroSkeleton processing package.
///
/// %ZeroSkeleton is a processing package that handles the 0-skeleton (vertices)
/// of a triangulation.
/// \sa Triangulation
/// \sa ttkTriangulation

#pragma once

#include <array>
#include <map>

// base code includes
#include <CellArray.h>
#include <FlatJaggedArray.h>
#include <Wrapper.h>

namespace ttk {

  class ZeroSkeleton : public Debug {

  public:
    ZeroSkeleton();

    /// Compute the list of edges connected to each vertex of a triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param edgeList List of edges. Each entry is represented by the
    /// ordered std::pair of identifiers of the entry's edge's vertices.
    /// \param vertexEdges Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be a std::vector listing the identifiers of the edges connected to the
    /// entry's vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexEdges(const SimplexId &vertexNumber,
                         const std::vector<std::array<SimplexId, 2>> &edgeList,
                         FlatJaggedArray &vertexEdges) const;

    /// Compute the link of each vertex of a 2D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param vertexStars List of vertex stars. The size of this std::vector
    /// should be equal to the number of vertices in the triangulation. Each
    /// entry is a std::vector listing the identifiers of triangles.
    /// \param cellEdges List of cell edges. The size of this std::vector
    /// should be equal to the number of triangles. Each entry is a std::array
    /// of identifiers of edges.
    /// \param edgeList List of edges. The size of this std::vector
    /// should be equal to the number of edges. Each entry is a std::array
    /// of vertex identifiers per edge.
    /// \param vertexLinks Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the triangulation. Each
    /// entry will be a std::vector listing the edges in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLinks(const FlatJaggedArray &vertexStars,
                         const std::vector<std::array<SimplexId, 3>> &cellEdges,
                         const std::vector<std::array<SimplexId, 2>> &edgeList,
                         FlatJaggedArray &vertexLinks) const;

    /// Compute the link of each vertex of a 3D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param vertexStars List of vertex stars. The size of this std::vector
    /// should be equal to the number of vertices in the triangulation. Each
    /// entry is a std::vector listing the identifiers of tetrahedra.
    /// \param cellTriangles List of cell triangles. The size of this
    /// std::vector
    /// should be equal to the number of tetrahedra. Each entry is a
    /// std::vector
    /// of identifiers of triangles.
    /// \param triangleList List of triangles. For each triangle, a
    /// std::array identifying its three vertices.
    /// \param vertexLinks Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the triangulation. Each
    /// entry will be a std::vector listing the triangles in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLinks(
      const FlatJaggedArray &vertexStars,
      const std::vector<std::array<SimplexId, 4>> &cellTriangles,
      const std::vector<std::array<SimplexId, 3>> &triangleList,
      FlatJaggedArray &vertexLinks) const;

    /// Compute the list of neighbors of each vertex of a triangulation.
    /// Unspecified behavior if the input mesh is not a valid triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param vertexNeighbors Output neighbor list. The size of this
    /// std::vector will be equal to the number of vertices in the mesh. Each
    /// entry will be std::vector listing the vertex identifiers of the entry's
    /// vertex' neighbors.
    /// \param edgeList List of edges. If this std::vector is not
    /// empty but incorrect, the behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexNeighbors(
      const SimplexId &vertexNumber,
      FlatJaggedArray &vertexNeighbors,
      const std::vector<std::array<SimplexId, 2>> &edgeList) const;

    /// Compute the star of each vertex of a triangulation. Unspecified
    /// behavior if the input mesh is not a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices
    /// of each cell.
    /// \param vertexStars Output vertex stars. The size of this std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be a std::vector listing the identifiers of the maximum-dimensional
    /// cells (3D: tetrahedra, 2D: triangles, etc.) connected to the entry's
    /// vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexStars(const SimplexId &vertexNumber,
                         const CellArray &cellArray,
                         FlatJaggedArray &vertexStars) const;
  };
} // namespace ttk
