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

#ifndef _TWOSKELETON_H
#define _TWOSKELETON_H

// base code includes
#include <OneSkeleton.h>
#include <Wrapper.h>
#include <ZeroSkeleton.h>

#include <algorithm>

namespace ttk {

  class TwoSkeleton : public Debug {

  public:
    TwoSkeleton();

    ~TwoSkeleton();

    /// Compute the list of cell-neighbors of each cell of a 2D triangulation
    /// (unspecified behavior if the input mesh is not a triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param cellNeighbors Output neighbor list. The size of this
    /// std::vector
    /// will be equal to the number of cells in the mesh. Each entry will be
    /// a std::vector listing the cell identifiers of the entry's cell's
    /// neighbors.
    /// \param vertexStars Optional list of vertex stars (list of
    /// 2-dimensional cells connected to each vertex). If NULL, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not NULL but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not NULL but pointing
    /// to a non-empty std::vector, this function will use this std::vector as
    /// internal
    /// vertex star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellNeighborsFromVertices(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &cellNeighbors,
      std::vector<std::vector<SimplexId>> *vertexStars = NULL) const;

    /// Compute the list of triangles connected to each edge for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param edgeTriangleList Output edge triangle list. The size of this
    /// std::vector will be equal to the number of edges in the triangulation.
    /// Each
    /// entry will be a std::vector listing the triangle identifiers for each
    /// triangle connected to the entry's edge.
    /// \param vertexStarList Optional output vertex star list (list of
    /// tetrahedron identifiers for each vertex). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// vertex star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \param edgeList Optional output edge list (list of std::pairs of
    /// vertex
    /// identifiers). If NULL, the function will compute this list anyway and
    /// free the related memory upon return. If not NULL but pointing to an
    /// empty std::vector, the function will fill this empty std::vector
    /// (useful if this
    /// list needs to be used later on by the calling program). If not NULL
    /// but pointing to a non-empty std::vector, this function will use this
    /// std::vector
    /// as internal edge list. If this std::vector is not empty but
    /// incorrect, the
    /// behavior is unspecified.
    /// \param edgeStarList Optional output edge star list (list of
    /// tetrahedron identifiers for each edge). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// edge star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \param triangleList Optional output triangle list (list of
    /// std::vectors of
    /// vertex identifiers). If NULL, the function will compute this list
    /// anyway and free the related memory upon return. If not NULL but
    /// pointing to an empty std::vector, the function will fill this empty
    /// std::vector
    /// (useful if this list needs to be used later on by the calling
    /// program). If not NULL but pointing to a non-empty std::vector, this
    /// function will use this std::vector as internal triangle list. If this
    /// std::vector is not empty but incorrect, the behavior is unspecified.
    /// \param triangleStarList Optional output triangle star list (list of
    /// tetrahedron identifiers for each triangle). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// triangle star list. If this std::vector is not empty but incorrect,
    /// the
    /// behavior is unspecified.
    /// \param cellTriangleList Optional output cell triangle list (list of
    /// triangle identifiers for each tetrahedron). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// cell triangle list. If this std::vector is not empty but incorrect,
    /// the
    /// behavior is unspecified.
    int buildEdgeTriangles(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &edgeTriangleList,
      std::vector<std::vector<SimplexId>> *vertexStarList = NULL,
      std::vector<std::pair<SimplexId, SimplexId>> *edgeList = NULL,
      std::vector<std::vector<SimplexId>> *edgeStarList = NULL,
      std::vector<std::vector<SimplexId>> *triangleList = NULL,
      std::vector<std::vector<SimplexId>> *triangleStarList = NULL,
      std::vector<std::vector<SimplexId>> *cellTriangleList = NULL) const;

    /// Compute the list of triangles of a triangulation represented by a
    /// vtkUnstructuredGrid object. Unspecified behavior if the input mesh is
    /// not a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, only)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param triangleList Optional output triangle list (each entry is the
    /// ordered std::vector of the vertex identifiers of the entry's
    /// triangle).
    /// \param triangleStars Optional output for triangle tet-adjacency (for
    /// each triangle, list of its adjacent tetrahedra).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleList(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> *triangleList = NULL,
      std::vector<std::vector<SimplexId>> *triangleStars = NULL,
      std::vector<std::vector<SimplexId>> *cellTriangleList = NULL) const;

    /// Compute the list of edges connected to each triangle for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// 3D triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param triangleEdgeList Output triangle edge list. The size of this
    /// std::vector will be equal to the number of triangles in the
    /// triangulation.
    /// Each entry will be a std::vector listing the edge identifiers for
    /// each
    /// edge connected to the entry's triangle.
    /// \param vertexEdgeList Optional output vertex edge list (list of
    /// edge identifiers for each vertex). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// vertex edge list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \param edgeList Optional output edge list (list of std::pairs of
    /// vertex
    /// identifiers). If NULL, the function will compute this list anyway and
    /// free the related memory upon return. If not NULL but pointing to an
    /// empty std::vector, the function will fill this empty std::vector
    /// (useful if this
    /// list needs to be used later on by the calling program). If not NULL
    /// but pointing to a non-empty std::vector, this function will use this
    /// std::vector
    /// as internal edge list. If this std::vector is not empty but
    /// incorrect, the
    /// behavior is unspecified.
    /// \param triangleList Optional output triangle list (list of
    /// std::vectors of
    /// vertex identifiers). If NULL, the function will compute this list
    /// anyway and free the related memory upon return. If not NULL but
    /// pointing to an empty std::vector, the function will fill this empty
    /// std::vector
    /// (useful if this list needs to be used later on by the calling
    /// program). If not NULL but pointing to a non-empty std::vector, this
    /// function will use this std::vector as internal triangle list. If this
    /// std::vector is not empty but incorrect, the behavior is unspecified.
    /// \param triangleStarList Optional output triangle star list (list of
    /// tetrahedron identifiers for each triangle). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// triangle star list. If this std::vector is not empty but incorrect,
    /// the
    /// behavior is unspecified.
    /// \param cellTriangleList Optional output cell triangle list (list of
    /// triangle identifiers for each tetrahedron). If NULL, the function
    /// will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty std::vector, the
    /// function
    /// will fill this empty std::vector (useful if this list needs to be
    /// used
    /// later on by the calling program). If not NULL but pointing to a
    /// non-empty std::vector, this function will use this std::vector as
    /// internal
    /// cell triangle list. If this std::vector is not empty but incorrect,
    /// the
    /// behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleEdgeList(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &triangleEdgeList,
      std::vector<std::vector<SimplexId>> *vertexEdgeList = NULL,
      std::vector<std::pair<SimplexId, SimplexId>> *edgeList = NULL,
      std::vector<std::vector<SimplexId>> *triangleList = NULL,
      std::vector<std::vector<SimplexId>> *triangleStarList = NULL,
      std::vector<std::vector<SimplexId>> *cellTriangleList = NULL) const;

    /// Compute the links of triangles in a 3D triangulation.
    /// \param triangleList Input triangle list. The number of entries of this
    /// list is equal to the number of triangles in the triangulation. Each
    /// entry lists the vertex identifiers of the corresponding triangle.
    /// \param triangleStar Input triangle star list. The number of entries of
    /// this list is equal to the number of triangles in the triangulation.
    /// Each entry lists the identifiers of the tetrahedra which are the
    /// co-faces of the corresponding triangle.
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param triangleLinks Output triangle link list. The number of entries
    /// of this list is equal to the number of triangles in the triangulation.
    /// Each entry lists the identifiers of the vertices in the link of the
    /// corresponding triangle.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangleLinks(
      const std::vector<std::vector<SimplexId>> &triangeList,
      const std::vector<std::vector<SimplexId>> &triangleStars,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &triangleLinks) const;

    /// Compute the list of triangles connected to each vertex for 3D
    /// triangulations (unspecified behavior if the input mesh is not a
    /// triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param triangleList Input triangle list (list of std::vectors of
    /// vertex identifiers).
    /// \param vertexTriangleList Output vertex triangle list (list of
    /// std::vectors of triangle identifiers).
    int buildVertexTriangles(
      const SimplexId &vertexNumber,
      const std::vector<std::vector<SimplexId>> &triangleList,
      std::vector<std::vector<SimplexId>> &vertexTriangleList) const;

  protected:
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <TwoSkeleton.cpp>

#endif // TWOSKELETON_H
