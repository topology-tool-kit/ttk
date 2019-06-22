/// \ingroup base
/// \class ttk::OneSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief %OneSkeleton processing package.
///
/// %OneSkeleton is a processing package that handles the 1-skeleton (edges)
/// of a triangulation.
/// \sa Triangulation
/// \sa ttkTriangulation

#ifndef _ONESKELETON_H
#define _ONESKELETON_H

#include <map>

// base code includes
#include <Wrapper.h>
#include <ZeroSkeleton.h>

namespace ttk {

  class OneSkeleton : public Debug {

  public:
    OneSkeleton();

    ~OneSkeleton();

    /// Compute the link of each edge of a 2D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param edgeList List of edges. The size of this std::vector
    /// should be equal to the number of edges in the triangulation. Each
    /// entry is a std::pair of vertex identifiers.
    /// \param edgeStars List of edge stars. The size of this std::vector
    /// should be
    /// equal to the number of edges. Each entry is a std::vector of triangle
    /// identifiers.
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param edgeLinks Output edge links. The size of this std::vector
    /// will be equal to the number of edges in the triangulation. Each
    /// entry will be a std::vector listing the vertices in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeLinks(
      const std::vector<std::pair<SimplexId, SimplexId>> &edgeList,
      const std::vector<std::vector<SimplexId>> &edgeStars,
      const LongSimplexId *cellArray,
      std::vector<std::vector<SimplexId>> &edgeLinks) const;

    /// Compute the link of each edge of a 3D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param edgeList List of edges. The size of this std::vector
    /// should be equal to the number of edges in the triangulation. Each
    /// entry is a std::pair of vertex identifiers.
    /// \param edgeStars List of edge stars. The size of this std::vector
    /// should be
    /// equal to the number of edges. Each entry is a std::vector of
    /// tetrahedron
    /// identifiers.
    /// \param cellEdges List of celle edges. The size of this std::vector
    /// should be equal to the number of tetrahedra in the triangulation. Each
    /// entry is a std::vector of edge identifiers.
    /// \param edgeLinks Output edge links. The size of this std::vector
    /// will be equal to the number of edges in the triangulation. Each
    /// entry will be a std::vector listing the vertices in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeLinks(
      const std::vector<std::pair<SimplexId, SimplexId>> &edgeList,
      const std::vector<std::vector<SimplexId>> &edgeStars,
      const std::vector<std::vector<SimplexId>> &cellEdges,
      std::vector<std::vector<SimplexId>> &edgeLinks) const;

    /// Compute the list of edges of a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param edgeList Output edge list (each entry is an ordered std::pair
    /// of
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeList(
      const SimplexId &vertexNumber,
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::pair<SimplexId, SimplexId>> &edgeList) const;

    /// Compute the list of edges of multiple triangulations.
    /// \param cellArrays Vector of cells. For each triangulation, each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param edgeList Output edge list (each entry is an ordered std::pair
    //  of
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int
      buildEdgeLists(const std::vector<std::vector<LongSimplexId>> &cellArrays,
                     std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
                       &edgeLists) const;

    /// Compute the 3-star of all the edges of a triangulation (for each
    /// edge, list of the 3-dimensional cells connected to it).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param starList Output list of 3-stars. The size of this std::vector
    /// will
    /// be equal to the number of edges in the mesh. Each entry stores a
    /// std::vector that lists the identifiers of all 3-dimensional cells
    /// connected to the entry's edge.
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
    /// \param vertexStars Optional list of vertex stars (list of
    /// 3-dimensional cells connected to each vertex). If NULL, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not NULL but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not NULL but pointing
    /// to a non-empty std::vector, this function will use this std::vector as
    /// internal
    /// vertex star list. If this std::vector is not empty but incorrect, the
    /// behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeStars(const SimplexId &vertexNumber,
                       const SimplexId &cellNumber,
                       const LongSimplexId *cellArray,
                       std::vector<std::vector<SimplexId>> &starList,
                       std::vector<std::pair<SimplexId, SimplexId>> *edgeList
                       = NULL,
                       std::vector<std::vector<SimplexId>> *vertexStars
                       = NULL) const;

    /// Compute the list of edges of a sub-portion of a valid triangulation.
    /// \param cellNumber Number of maximum-dimensional cells in the
    /// considered subset of the triangulation (number of tetrahedra in 3D,
    /// triangles in 2D, etc.)
    /// \param cellArray Pointer to a contiguous array of cells. Each entry
    /// starts by the number of vertices in the cell, followed by the vertex
    /// identifiers of the cell.
    /// \param edgeList Output edge list (each entry is an ordered std::pair
    /// of
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeSubList(
      const SimplexId &cellNumber,
      const LongSimplexId *cellArray,
      std::vector<std::pair<SimplexId, SimplexId>> &edgeList) const;

  protected:
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <OneSkeleton.cpp>

#endif // ONESKELETON_H
