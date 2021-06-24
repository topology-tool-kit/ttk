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

#pragma once

#include <array>
#include <map>

// base code includes
#include <CellArray.h>
#include <Wrapper.h>
#include <ZeroSkeleton.h>

namespace ttk {

  class OneSkeleton : public Debug {

  public:
    OneSkeleton();

    /// Compute the link of each edge of a 2D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param edgeList List of edges. The size of this std::vector
    /// should be equal to the number of edges in the triangulation. Each
    /// entry is a std::pair of vertex identifiers.
    /// \param edgeStars List of edge stars. The size of this std::vector
    /// should be equal to the number of edges. Each entry is a std::vector of
    /// triangle identifiers.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param edgeLinks Output edge links. The size of this std::vector will be
    /// equal to the number of edges in the triangulation. Each entry will be a
    /// std::vector listing the vertices in the link of the corresponding
    /// vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeLinks(const std::vector<std::array<SimplexId, 2>> &edgeList,
                       const FlatJaggedArray &edgeStars,
                       const CellArray &cellArray,
                       FlatJaggedArray &edgeLinks) const;

    /// Compute the link of each edge of a 3D triangulation (unspecified
    /// behavior if the input mesh is not a valid triangulation).
    /// \param edgeList List of edges. The size of this std::vector
    /// should be equal to the number of edges in the triangulation. Each
    /// entry is a std::pair of vertex identifiers.
    /// \param edgeStars List of edge stars. The size of this std::vector
    /// should be equal to the number of edges. Each entry is a std::vector of
    /// tetrahedron identifiers.
    /// \param cellEdges List of celle edges. The size of this std::vector
    /// should be equal to the number of tetrahedra in the triangulation. Each
    /// entry is a std::vector of edge identifiers.
    /// \param edgeLinks Output edge links. The size of this std::vector
    /// will be equal to the number of edges in the triangulation. Each
    /// entry will be a std::vector listing the vertices in the link of the
    /// corresponding vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeLinks(const std::vector<std::array<SimplexId, 2>> &edgeList,
                       const FlatJaggedArray &edgeStars,
                       const std::vector<std::array<SimplexId, 6>> &cellEdges,
                       FlatJaggedArray &edgeLinks) const;

    /// Compute the list of edges of a valid triangulation.
    /// \param vertexNumber Number of vertices in the triangulation.
    /// \param cellArray Cell container allowing to retrieve the vertices ids
    /// of each cell.
    /// \param edgeList Output edge list (each entry is an ordered
    /// std::array of vertex identifiers).
    /// \param edgeStars Output for edge cell adjacency (for each
    /// edge, a list of adjacent cells)
    /// \param cellEdgeList Output for cell edges: per cell, the list
    /// of its edges identifiers
    /// \return Returns 0 upon success, negative values otherwise.
    template <std::size_t n>
    int
      buildEdgeList(const SimplexId &vertexNumber,
                    const CellArray &cellArray,
                    std::vector<std::array<SimplexId, 2>> &edgeList,
                    FlatJaggedArray &edgeStars,
                    std::vector<std::array<SimplexId, n>> &cellEdgeList) const;
  };
} // namespace ttk
