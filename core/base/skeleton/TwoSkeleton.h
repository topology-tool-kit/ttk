/// \ingroup base
/// \class ttk::TwoSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief %TwoSkeleton processing package.
///
/// %TwoSkeleton is a processing package that handles the 2-skeleton (triangles)
/// of a triangulation.
/// \sa Triangulation
/// \sa ttkTriangulation

#pragma once

// base code includes
#include <Debug.h>
#include <ZeroSkeleton.h>

#include <algorithm>
#include <array>

namespace ttk {

  class TwoSkeleton : public Debug {

  public:
    TwoSkeleton();

    /// Compute the list of cell-neighbors of each cell of a 2D triangulation
    /// (unspecified behavior if the input mesh is not a triangulation).
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param cellNeighbors Output neighbor list. The size of this std::vector
    /// will be equal to the number of cells in the mesh. Each entry will be a
    /// std::vector listing the cell identifiers of the entry's cell's
    /// neighbors.
    /// \param edgeStars Array of edge stars (list of 2-dimensional
    /// cells connected to each edge).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellNeighborsFromEdges(const CellArray &cellArray,
                                    FlatJaggedArray &cellNeighbors,
                                    const FlatJaggedArray &edgeStars) const;

    /// Compute the list of cell-neighbors of each cell of a 2D triangulation
    /// (unspecified behavior if the input mesh is not a triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param cellNeighbors Output neighbor list. The size of this std::vector
    /// will be equal to the number of cells in the mesh. Each entry will be a
    /// std::vector listing the cell identifiers of the entry's cell's
    /// neighbors.
    /// \param vertexStars Optional list of vertex stars (list of
    /// 2-dimensional cells connected to each vertex). If nullptr, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not nullptr but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not nullptr but pointing
    /// to a non-empty std::vector, this function will use this std::vector as
    /// internal vertex star list. If this std::vector is not empty but
    /// incorrect, the behavior is unspecified.
    ///  \return Returns 0 upon success, negative values otherwise.
    int buildCellNeighborsFromVertices(const SimplexId &vertexNumber,
                                       const CellArray &cellArray,
                                       FlatJaggedArray &cellNeighbors,
                                       FlatJaggedArray *vertexStars
                                       = nullptr) const;

    /// Compute the list of triangles connected to each edge for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param edgeTriangleList Output edge triangle list. The size of this
    /// std::vector will be equal to the number of edges in the triangulation.
    /// Each entry will be a std::vector listing the triangle identifiers for
    /// each triangle connected to the entry's edge.
    /// \param edgeList Edge list (list of std::pairs of vertex
    /// identifiers). If this std::vector is not empty but incorrect,
    /// the behavior is unspecified.
    /// \param triangleEdgeList Optional output triangle edge list (list of
    /// std::vectors of edges identifiers per triangle). If nullptr, the
    /// function will compute this list anyway and free the related memory upon
    /// return. If not nullptr but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs to
    /// be used later on by the calling program). If not nullptr but pointing to
    /// a non-empty std::vector, this function will use this std::vector as
    /// internal edge list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    int
      buildEdgeTriangles(const SimplexId &vertexNumber,
                         const CellArray &cellArray,
                         FlatJaggedArray &edgeTriangleList,
                         const std::vector<std::array<SimplexId, 2>> &edgeList,
                         std::vector<std::array<SimplexId, 3>> *triangleEdgeList
                         = nullptr) const;

    /// Compute the list of triangles of a triangulation represented by a
    /// vtkUnstructuredGrid object. Unspecified behavior if the input mesh is
    /// not a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param triangleList Optional output triangle list (each entry is the
    /// ordered std::vector of the vertex identifiers of the entry's triangle).
    /// \param triangleStars Optional output for triangle tet-adjacency (for
    /// each triangle, list of its adjacent tetrahedra).
    /// \param cellTriangleList Optional list of triangles per
    /// tetrahedron cell.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleList(
      const SimplexId &vertexNumber,
      const CellArray &cellArray,
      std::vector<std::array<SimplexId, 3>> *triangleList = nullptr,
      FlatJaggedArray *triangleStars = nullptr,
      std::vector<std::array<SimplexId, 4>> *cellTriangleList = nullptr) const;

    /// Compute the list of edges connected to each triangle for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// 3D triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param triangleEdgeList Output triangle edge list. The size of this
    /// std::vector will be equal to the number of triangles in the
    /// triangulation. Each entry will be a std::vector listing the edge
    /// identifiers for each edge connected to the entry's triangle.
    /// \param edgeList Edge list (list of std::pairs of vertex
    /// identifiers). If this std::vector is not empty but incorrect,
    /// the behavior is unspecified.
    /// \param vertexEdgeList Optional output vertex edge list (list of edge
    /// identifiers for each vertex). If nullptr, the function will compute this
    /// list anyway and free the related memory upon return. If not nullptr but
    /// pointing to an empty std::vector, the function will fill this empty
    /// std::vector (useful if this list needs to be used later on by the
    /// calling program). If not nullptr but pointing to a non-empty
    /// std::vector, this function will use this std::vector as internal vertex
    /// edge list. If this std::vector is not empty but incorrect, the behavior
    /// is unspecified.
    /// \param triangleList Optional output triangle list (list of std::vectors
    /// of vertex identifiers). If nullptr, the function will compute this list
    /// anyway and free the related memory upon return. If not nullptr but
    /// pointing to an empty std::vector, the function will fill this empty
    /// std::vector (useful if this list needs to be used later on by the
    /// calling program). If not nullptr but pointing to a non-empty
    /// std::vector, this function will use this std::vector as internal
    /// triangle list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \param triangleStarList Optional output triangle star list (list of
    /// tetrahedron identifiers for each triangle). If nullptr, the function
    /// will compute this list anyway and free the related memory upon return.
    /// If not nullptr but pointing to an empty std::vector, the function will
    /// fill this empty std::vector (useful if this list needs to be used later
    /// on by the calling program). If not nullptr but pointing to a non-empty
    /// std::vector, this function will use this std::vector as internal
    /// triangle star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \param cellTriangleList Optional output cell triangle list (list of
    /// triangle identifiers for each tetrahedron). If nullptr, the function
    /// will compute this list anyway and free the related memory upon return.
    /// If not nullptr but pointing to an empty std::vector, the function will
    /// fill this empty std::vector (useful if this list needs to be used later
    /// on by the calling program). If not nullptr but pointing to a non-empty
    /// std::vector, this function will use this std::vector as internal cell
    /// triangle list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleEdgeList(
      const SimplexId &vertexNumber,
      const CellArray &cellArray,
      std::vector<std::array<SimplexId, 3>> &triangleEdgeList,
      const std::vector<std::array<SimplexId, 2>> &edgeList,
      FlatJaggedArray *vertexEdgeList = nullptr,
      std::vector<std::array<SimplexId, 3>> *triangleList = nullptr,
      FlatJaggedArray *triangleStarList = nullptr,
      std::vector<std::array<SimplexId, 4>> *cellTriangleList = nullptr) const;

    /// Compute the links of triangles in a 3D triangulation.
    /// \param triangleList Input triangle list. The number of entries of this
    /// list is equal to the number of triangles in the triangulation. Each
    /// entry lists the vertex identifiers of the corresponding triangle.
    /// \param triangleStars Input triangle star list. The number of entries of
    /// this list is equal to the number of triangles in the triangulation. Each
    /// entry lists the identifiers of the tetrahedra which are the co-faces of
    /// the corresponding triangle.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param triangleLinks Output triangle link list. The number of entries of
    /// this list is equal to the number of triangles in the triangulation. Each
    /// entry lists the identifiers of the vertices in the link of the
    /// corresponding triangle.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleLinks(
      const std::vector<std::array<SimplexId, 3>> &triangleList,
      const FlatJaggedArray &triangleStars,
      const CellArray &cellArray,
      FlatJaggedArray &triangleLinks) const;

    /// Compute the list of triangles connected to each vertex for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param triangleList Input triangle list (list of std::vectors of
    /// vertex identifiers).
    /// \param vertexTriangles Output vertex triangle list (list of
    /// std::vectors of triangle identifiers).
    int buildVertexTriangles(
      const SimplexId &vertexNumber,
      const std::vector<std::array<SimplexId, 3>> &triangleList,
      FlatJaggedArray &vertexTriangles) const;
  };
} // namespace ttk
