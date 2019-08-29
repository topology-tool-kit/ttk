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

#ifndef _ZEROSKELETON_H
#define _ZEROSKELETON_H

#include <map>

// base code includes
#include <Wrapper.h>

namespace ttk {

  class ZeroSkeleton : public Debug {

  public:
    ZeroSkeleton();

    ~ZeroSkeleton();

    /// Compute the list of edges connected to each vertex of a triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param edgeList List of edges. Each entry is represented by the
    /// ordered std::pair of identifiers of the entry's edge's vertices.
    /// \param vertexEdges Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be a std::vector listing the identifiers of the edges connected to the
    /// entry's vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexEdges(
      const SimplexId &vertexNumber,
      const std::vector<std::pair<SimplexId, SimplexId>> &edgeList,
      std::vector<std::vector<SimplexId>> &vertexEdges) const;

    /// Compute the link of a single vertex of a triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param vertexId Input vertex.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param vertexLink Output vertex link. This std::vector contains, for
    /// each
    /// simplex of the link, the number of vertices in the simplex (triangles:
    /// 3, edges: 2) followed by the corresponding vertex identifiers.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLink(const SimplexId &vertexId,
                        const SimplexId &cellNumber,
                        const LongSimplexId *cellArray,
                        std::vector<LongSimplexId> &vertexLink) const;

    /// Compute the link of each vertex of a triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param vertexLinks Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be a std::vector listing the simplices of the link of the entry's
    /// vertex.
    /// In particular, this std::vector contains, for each simplex, the
    /// number of
    /// vertices in the simplex (triangles: 3, edges: 2) followed by the
    /// corresponding vertex identifiers.
    /// \param vertexStars Optional list of vertex stars (list of
    /// 3-dimensional cells connected to each vertex). If NULL, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not NULL but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not NULL but pointing
    /// to a non-empty std::vector, this function will use this std::vector
    /// as internal
    /// vertex star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLinks(const SimplexId &vertexNumber,
                         const SimplexId &cellNumber,
                         const LongSimplexId *cellArray,
                         std::vector<std::vector<LongSimplexId>> &vertexLinks,
                         std::vector<std::vector<SimplexId>> *vertexStars
                         = NULL) const;

    /// Compute the link of each vertex of a 2D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param vertexStars List of vertex stars. The size of this std::vector
    /// should be equal to the number of vertices in the triangulation. Each
    /// entry is a std::vector listing the identifiers of triangles.
    /// \param cellEdges List of cell edges. The size of this std::vector
    /// should be
    /// equal to the number of triangles. Each entry is a std::vector of
    /// identifiers of edges.
    /// \param vertexLinks Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the triangulation. Each
    /// entry will be a std::vector listing the edges in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLinks(
      const std::vector<std::vector<SimplexId>> &vertexStars,
      const std::vector<std::vector<SimplexId>> &cellEdges,
      const std::vector<std::pair<SimplexId, SimplexId>> &edgeList,
      std::vector<std::vector<SimplexId>> &vertexLinks) const;

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
    /// \param vertexLinks Output vertex links. The size of this std::vector
    /// will be equal to the number of vertices in the triangulation. Each
    /// entry will be a std::vector listing the triangles in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int
      buildVertexLinks(const std::vector<std::vector<SimplexId>> &vertexStars,
                       const std::vector<std::vector<SimplexId>> &cellTriangles,
                       const std::vector<std::vector<SimplexId>> &triangleList,
                       std::vector<std::vector<SimplexId>> &vertexLinks) const;

    /// Compute the list of neighbors of each vertex of a triangulation.
    /// Unspecified behavior if the input mesh is not a valid triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param vertexNeighbors Output neighbor list. The size of this
    /// std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be std::vector listing the vertex identifiers of the entry's vertex'
    /// neighbors.
    /// \param edgeList Optional list of edges. If NULL, the function will
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty std::vector, the function will
    /// fill
    /// this empty std::vector (useful if this list needs to be used later on
    /// by
    /// the calling program). If not NULL but pointing to a non-empty
    /// std::vector,
    /// this function will use this std::vector as internal edge list. If
    /// this
    /// std::vector is not empty but incorrect, the behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexNeighbors(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &vertexNeighbors,
      std::vector<std::pair<SimplexId, SimplexId>> *edgeList = NULL) const;

    /// Compute the star of each vertex of a triangulation. Unspecified
    /// behavior if the input mesh is not a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param vertexStars Output vertex stars. The size of this std::vector
    /// will be equal to the number of vertices in the mesh. Each entry will
    /// be a std::vector listing the identifiers of the maximum-dimensional
    /// cells
    /// (3D: tetrahedra, 2D: triangles, etc.) connected to the entry's vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int
      buildVertexStars(const SimplexId &vertexNumber,
                       const SimplexId &cellNumber,
                       const LongSimplexId *cellArray,
                       std::vector<std::vector<SimplexId>> &vertexStars) const;

  protected:
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <ZeroSkeleton.cpp>

#endif // ZEROSKELETON_H
