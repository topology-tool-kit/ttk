/// \ingroup base
/// \class ttk::ThreeSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief %ThreeSkeleton processing package.
///
/// %ThreeSkeleton is a processing package that handles the 3-skeleton
/// (tetrahedra) of a triangulation.
/// \sa Triangulation
/// \sa ttkTriangulation

#ifndef _THREESKELETON_H
#define _THREESKELETON_H

// base code includes
#include <OneSkeleton.h>
#include <TwoSkeleton.h>
#include <Wrapper.h>
#include <ZeroSkeleton.h>

#include <algorithm>

namespace ttk {

  class ThreeSkeleton : public Debug {

  public:
    ThreeSkeleton();

    ~ThreeSkeleton();

    /// Compute the list of edges of each cell of a triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param cellEdges Output edge lists. The size of this std::vector
    /// will be equal to the number of cells in the mesh. Each entry will be
    /// a std::vector listing the edge identifiers of the entry's cell's
    /// edges.
    /// \param edgeList Optional list of edges. If nullptr, the function will
    /// compute this list anyway and free the related memory upon return. If not
    /// nullptr but pointing to an empty std::vector, the function will fill
    /// this empty std::vector (useful if this list needs to be used later on by
    /// the calling program). If not nullptr but pointing to a non-empty
    /// std::vector, this function will use this std::vector as internal edge
    /// list. If this std::vector is not empty but incorrect, the behavior is
    /// unspecified.
    /// \param vertexEdges Optional list of edges for each vertex.
    /// If nullptr, the function will compute this list anyway and free the
    /// related memory upon return. If not nullptr but pointing to an empty
    /// std::vector, the function will fill this empty std::vector (useful if
    /// this list needs to be used later on by the calling program). If not
    /// nullptr but pointing to a non-empty std::vector, this function will use
    /// this std::vector as internal vertex edge list. If this std::vector is
    /// not empty but incorrect, the behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellEdges(const SimplexId &vertexNumber,
                       const CellArray &cellArray,
                       std::vector<std::vector<SimplexId>> &cellEdges,
                       std::vector<std::pair<SimplexId, SimplexId>> *edgeList
                       = nullptr,
                       std::vector<std::vector<SimplexId>> *vertexEdges
                       = nullptr) const;

    /// Compute the list of cell-neighbors of each cell of a triangulation
    /// (unspecified behavior if the input mesh is not a triangulation).
    /// This implementation is fast only if you already have the triangle
    /// stars computed. Otherwise, please use
    /// ThreeSkeleton::buildCellNeighborsFromVertices instead.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param cellNeighbors Output neighbor list. The size of this std::vector
    /// will be equal to the number of cells in the mesh. Each entry will be a
    /// std::vector listing the cell identifiers of the entry's cell's
    /// neighbors.
    /// \param triangleStars Optional list of triangle stars (list of
    /// 3-dimensional cells connected to each triangle). If nullptr, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not nullptr but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not nullptr but pointing
    /// to a non-empty std::vector, this function will use this std::vector as
    /// internal triangle star list. If this std::vector is not empty but
    /// incorrect, the behavior is unspecified. \return Returns 0 upon success,
    /// negative values otherwise.
    int buildCellNeighborsFromTriangles(
      const SimplexId &vertexNumber,
      const CellArray &cellArray,
      std::vector<std::vector<SimplexId>> &cellNeighbors,
      std::vector<std::vector<SimplexId>> *triangleStars = nullptr) const;

    /// Compute the list of cell-neighbors of each cell of a triangulation
    /// (unspecified behavior if the input mesh is not a triangulation).
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param cellNeighbors Output neighbor list. The size of this
    /// std::vector will be equal to the number of cells in the mesh. Each entry
    /// will be a std::vector listing the cell identifiers of the entry's cell's
    /// neighbors.
    /// \param vertexStars Optional list of vertex stars (list of
    /// 3-dimensional cells connected to each vertex). If nullptr, the
    /// function will compute this list anyway and free the related memory
    /// upon return. If not nullptr but pointing to an empty std::vector, the
    /// function will fill this empty std::vector (useful if this list needs
    /// to be used later on by the calling program). If not nullptr but pointing
    /// to a non-empty std::vector, this function will use this std::vector as
    /// internal vertex star list. If this std::vector is not empty but
    /// incorrect, the behavior is unspecified.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellNeighborsFromVertices(
      const SimplexId &vertexNumber,
      const CellArray &cellArray,
      std::vector<std::vector<SimplexId>> &cellNeighbors,
      std::vector<std::vector<SimplexId>> *vertexStars = nullptr) const;
  };
} // namespace ttk

#endif // THREESKELETON_H
