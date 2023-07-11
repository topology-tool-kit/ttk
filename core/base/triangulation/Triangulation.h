/// \ingroup base
/// \class ttk::Triangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief Triangulation is a class that provides time and memory efficient
/// traversal methods on triangulations of piecewise linear manifolds. It
/// provides the following features:
///   -# Given a vertex, it provides: the list of edges that are connected to
/// it, the list of its neighbors, its link, its star, etc.
///   -# Given an edge, it provides: its vertices, its star, etc.
///   -# Given a triangle, its provides: its vertices, its edges, etc.
///   -# Given a tetrahedron, its provides: its vertices, its edges, its
/// neighbor tetrahedra, etc.
///   -# Given a triangulation, it provides: its list of vertices, edges,
/// triangles and tetrahedra.
///
/// Triangulation supports both explicit and implicit triangulations:
///   -# Explicit triangulations: Given a list of points and a list of cells,
/// Triangulation provides time efficient accesses (requiring adequate
/// preconditioning, see the documentation furtherdown).
///   -# Implicit triangulations: Given a regular grid (origin, spacings and
/// dimensions), Triangulation will perform an implicit triangulation of the
/// grid, enabling both time and memory efficient traversals of triangulations
/// of regular grids.
///
/// Apart from preconditions, Triangulation requires no memory overhead in
/// addition to the input data.
///
/// \note
/// Only precondition the information you need! See the documentation further
/// down.
/// \sa ttkTriangulation

#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CompactTriangulation.h>
#include <ExplicitTriangulation.h>
#include <ImplicitTriangulation.h>
#include <PeriodicImplicitTriangulation.h>

#include <array>

namespace ttk {

  class Triangulation final : public AbstractTriangulation {

  public:
    Triangulation();
    Triangulation(const Triangulation &);
    Triangulation(Triangulation &&) noexcept;
    Triangulation &operator=(const Triangulation &);
    Triangulation &operator=(Triangulation &&) noexcept;
    ~Triangulation() override;

    enum class Type {
      EXPLICIT,
      IMPLICIT,
      HYBRID_IMPLICIT,
      PERIODIC,
      HYBRID_PERIODIC,
      COMPACT
    };

    /**
     * Strategies for implicit & periodic triangulations preconditioning
     */
    enum class STRATEGY {
      /** No preconditions above a number of vertices threshold (given
          by TTK_IMPLICIT_PRECONDITIONS_THRESHOLD, default to
          256^3) */
      DEFAULT = 0,
      /** Always precondition implicit & periodic triangulations */
      WITH_PRECONDITIONS = 1,
      /** Never precondition implicit & periodic triangulations */
      NO_PRECONDITIONS = 2,
    };

    /// Reset the triangulation data-structures.
    /// \return Returns 0 upon success, negative values otherwise.
    inline int clear() {

      if(abstractTriangulation_) {
        return abstractTriangulation_->clear();
      }

      return 0;
    }

    /// Computes and displays the memory footprint of the data-structure.
    /// \return Returns 0 upon success, negative values otherwise.
    inline size_t footprint() const {

      if(abstractTriangulation_) {
        return abstractTriangulation_->footprint();
      }

      return 0;
    }

    /// Get the \p localEdgeId-th edge of the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// In 1D, this function is equivalent to getCellNeighbor().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellEdges() needs to be called on this object prior to any
    /// traversal, in a clearly distinct pre-processing step that involves no
    /// traversal at all. An error will be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getCellEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellEdgeNumber()
    /// \sa getCellNeighbor()
    inline int getCellEdge(const SimplexId &cellId,
                           const int &localEdgeId,
                           SimplexId &edgeId) const override {

      // initialize output variable before early return
      edgeId = -1;

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->getCellEdge(cellId, localEdgeId, edgeId);
    }

    /// Get the number of edges for the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// In 1D, this function is equivalent to getCellNeighborNumber().
    ///
    /// \pre For this function to behave correctly, preconditionCellEdges()
    /// needs to be called on this object prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell edges.
    /// \sa getCellNeighborNumber()
    inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellEdgeNumber(cellId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of edges for all cells.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// The number of entries in this list is equal to the number of cells.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of edges for the corresponding cell.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 1D, this function is equivalent to getCellNeighbors().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell edge list.
    /// \sa getCellNeighbors()
    inline const std::vector<std::vector<SimplexId>> *getCellEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getCellEdges();
    }

    inline int
      getCellIncenter(SimplexId cellId, int dim, float incenter[3]) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellIncenter(cellId, dim, incenter);
    }

    /// Get the \p localNeighborId-th cell neighbor of the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localNeighborId Input local neighbor identifier,
    /// in [0, getCellNeighborNumber()].
    /// \param neighborId Output global neighbor cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellNeighborNumber()
    inline int getCellNeighbor(const SimplexId &cellId,
                               const int &localNeighborId,
                               SimplexId &neighborId) const override {

      // initialize output variable before early return
      neighborId = -1;
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellNeighbor(
        cellId, localNeighborId, neighborId);
    }

    /// Get the number of cell neighbors for the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell neighbors.
    inline SimplexId
      getCellNeighborNumber(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellNeighborNumber(cellId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of cell neighbors for all cells.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// The number of entries in this list is equal to the number of cells.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of neighbor cells for the corresponding cell.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell neighbor list.
    inline const std::vector<std::vector<SimplexId>> *
      getCellNeighbors() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getCellNeighbors();
    }

    /// Get the \p localTriangleId-th triangle id of the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// In 2D, this function is equivalent to getCellNeighbor().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getCellTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellTriangleNumber()
    /// \sa getCellNeighbor()
    inline int getCellTriangle(const SimplexId &cellId,
                               const int &localTriangleId,
                               SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellTriangle(
        cellId, localTriangleId, triangleId);
    }

    /// Get the number of triangles for the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// In 2D, this function is equivalent to getCellNeighborNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell triangles.
    /// \sa getCellNeighborNumber()
    inline SimplexId
      getCellTriangleNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellTriangleNumber(cellId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of triangles for all cells.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// Also, the notion of triangle only makes sense if the triangulation
    /// has a dimension greater than 2 (otherwise, use the cell information).
    ///
    /// The number of entries in this list is equal to the number of cells.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of triangles for the corresponding cell.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 2D, this function is equivalent to getCellNeighbors().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionCellTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell triangle list.
    /// \sa getCellNeighbors()
    inline const std::vector<std::vector<SimplexId>> *
      getCellTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getCellTriangles();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p cellId-th
    /// cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \param cellId Input global cell identifier.
    /// \param localVertexId Input local vertex identifier,
    /// in [0, getCellVertexNumber()].
    /// \param vertexId Output global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertexNumber()
    inline int getCellVertex(const SimplexId &cellId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->getCellVertex(
        cellId, localVertexId, vertexId);
    }

    /// Get the number of vertices in a cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \param cellId Input global cell identifier.
    /// \returns Number of vertices in the cell.
    inline SimplexId
      getCellVertexNumber(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getCellVertexNumber(cellId);
    }

    /// Get the dimensionality of the triangulation (this value is equal to
    /// the dimension of the simplex with largest dimensionality).
    /// \return Returns the dimensionality of the triangulation.
    inline int getDimensionality() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getDimensionality();
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of edges of the triangulation.
    ///
    /// Here the notion of edge only makes sense if the triangulation has a
    /// dimension greater than 1 (otherwise, use the cell information).
    ///
    /// The number of entries in this list is equal to the number of edges.
    /// Each entry is a std::pair of vertex identifiers.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge list.
    inline const std::vector<std::array<SimplexId, 2>> *getEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getEdges();
    }

    /// Compute the barycenter of the points of the given edge identifier.
    /// \pre For this function to behave correctly,
    /// preconditionEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \param edgeId Input global edge identifier.
    /// \param incenter Output barycenter.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleIncenter()
    /// \sa getCellIncenter()
    inline int getEdgeIncenter(SimplexId edgeId, float incenter[3]) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeIncenter(edgeId, incenter);
    }

    /// Compute the barycenter of the points of the given triangle identifier.
    /// \pre For this function to behave correctly,
    /// preconditionTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \param triangleId Input global triangle identifier.
    /// \param incenter Output barycenter.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeIncenter()
    /// \sa getCellIncenter()
    inline int getTriangleIncenter(SimplexId triangleId,
                                   float incenter[3]) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleIncenter(triangleId, incenter);
    }

    /// Get the \p localLinkId-th simplex of the link of the \p edgeId-th
    /// edge.
    ///
    /// The output \p linkId refers in 2D to a vertex identifier and in 3D
    /// to an edge identifier. It returns a negative value in 1D.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getEdgeLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeLinkNumber()
    inline int getEdgeLink(const SimplexId &edgeId,
                           const int &localLinkId,
                           SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeLink(edgeId, localLinkId, linkId);
    }

    /// Get the number of simplicies in the link of the \p edgeId-th edge.
    ///
    /// In 2D, this will return the number of vertices in the link, in 3D the
    /// number of edges. It returns a negative value in 1D.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of cells in the link of the edge.
    inline SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeLinkNumber(edgeId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of link simplices for all edges.
    ///
    /// The number of entries in this list is equal to the number of edges.
    /// Each entry is a std::vector of identifiers representing vertices in
    /// 2D and
    /// edges in 3D. It returns NULL in 1D.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge link list.
    inline const std::vector<std::vector<SimplexId>> *getEdgeLinks() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getEdgeLinks();
    }

    /// Get the \p localStarId-th cell of the star of the \p edgeId-th edge.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// Also, the notion of edge only makes sense if the triangulation has a
    /// dimension greater than 1 (otherwise, use the cell information).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param edgeId Input global edge identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getEdgeStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeStarNumber()
    inline int getEdgeStar(const SimplexId &edgeId,
                           const int &localStarId,
                           SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeStar(edgeId, localStarId, starId);
    }

    /// Get the number of star cells for the \p edgeId-th edge.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// Also, the notion of edge only makes sense if the triangulation has a
    /// dimension greater than 1 (otherwise, use the cell information).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier
    /// \return Returns the number of star cells.
    inline SimplexId getEdgeStarNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeStarNumber(edgeId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of star cell identifiers for all edges.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// Also, the notion of edge only makes sense if the triangulation has a
    /// dimension greater than 1 (otherwise, use the cell information).
    ///
    /// The number of entries in this list is equal to the number of edges.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of star cells for the corresponding edge.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge star list.
    inline const std::vector<std::vector<SimplexId>> *getEdgeStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getEdgeStars();
    }

    /// Get the \p localTriangleId-th triangle id of the \p edgeId-th edge.
    ///
    /// In 2D, this function is equivalent to getEdgeStar().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param edgeId Input global edge identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getEdgeTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeTriangleNumber()
    /// \sa getEdgeStar()
    inline int getEdgeTriangle(const SimplexId &edgeId,
                               const int &localTriangleId,
                               SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeTriangle(
        edgeId, localTriangleId, triangleId);
    }

    /// Get the number of triangles for the \p edgeId-th edge.
    ///
    /// In 2D, this function is equivalent to getEdgeStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of edge triangles.
    /// \sa getEdgeStarNumber
    inline SimplexId
      getEdgeTriangleNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeTriangleNumber(edgeId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of triangles for all edges.
    ///
    /// The number of entries in this list is equal to the number of edges.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of triangles for the corresponding edge.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 2D, this function is equivalent to getEdgeStars().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdgeTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge triangle list.
    /// \sa getEdgeStars
    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getEdgeTriangles();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p edgeId-th
    /// edge.
    ///
    /// In 1D, this function is equivalent to getCellVertex().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \param localVertexId Input local vertex identifier (0 or 1).
    /// \param vertexId Output global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertex()
    inline int getEdgeVertex(const SimplexId &edgeId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getEdgeVertex(
        edgeId, localVertexId, vertexId);
    }

    /// Get the internal abstract triangulation object.
    /// \return Returns a pointer to the internal abstract triangulation object.
    inline AbstractTriangulation *getData() {
      return abstractTriangulation_;
    }

    /// Get the number of cells in the triangulation.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \return Returns the number of cells.
    inline SimplexId getNumberOfCells() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->getNumberOfCells();
    }

    /// Get the number of edges in the triangulation.
    ///
    /// In 1D, this function is equivalent to getNumberOfCells().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns the number of edges.
    /// \sa getNumberOfCells()
    inline SimplexId getNumberOfEdges() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getNumberOfEdges();
    }

    /// Get the number of triangles in the triangulation.
    ///
    /// In 2D, this function is equivalent to getNumberOfCells().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns the number of triangles.
    /// \sa getNumberOfCells()
    inline SimplexId getNumberOfTriangles() const override {
#ifndef TTK_ENABLE_KAMIKAZE

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getNumberOfTriangles();
    }

    /// Get the number of vertices in the triangulation.
    /// \return Returns the number of vertices.
    inline SimplexId getNumberOfVertices() const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getNumberOfVertices();
    }

    /// Compute the barycenter of the points of the given tet identifier.
    /// \param tetraId Input global tet identifier.
    /// \param incenter Output barycenter.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleIncenter()
    /// \sa getEdgeIncenter()
    int getTetraIncenter(SimplexId tetraId, float incenter[3]) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTetraIncenter(tetraId, incenter);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of triangles of the triangulation.
    ///
    /// Here the notion of triangle only makes sense if the triangulation has
    /// a dimension greater than 2 (otherwise, use the cell information).
    ///
    /// The number of entries in this list is equal to the number of
    /// triangles.
    /// Each entry is a std::vector of vertex identifiers.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    /// \pre For this function to behave correctly,
    /// preconditionTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle list.
    inline const std::vector<std::array<SimplexId, 3>> *
      getTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getTriangles();
    }

    /// Get the \p localEdgeId-th edge of the \p triangleId-th triangle.
    ///
    /// In 2D, this function is equivalent to getCellEdge().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getTriangleEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleEdgeNumber()
    /// \sa getCellEdge()
    inline int getTriangleEdge(const SimplexId &triangleId,
                               const int &localEdgeId,
                               SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleEdge(
        triangleId, localEdgeId, edgeId);
    }

    /// Get the number of edges of the \p triangleId-th triangle.
    ///
    /// In 2D, this function is equivalent to getCellEdgeNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of cells in the link of the triangle.
    /// \sa getCellEdgeNumber()
    inline SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleEdgeNumber(triangleId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of edges for all triangles.
    ///
    /// The number of entries in this list is equal to the number of
    /// triangles. Each entry is a std::vector of identifiers representing the
    /// edges connected to the triangle (3).
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 2D, this function is equivalent to getCellEdges().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle edge list.
    /// \sa getCellEdges()
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getTriangleEdges();
    }

    /// Get the \p localLinkId-th simplex of the link of the \p triangleId-th
    /// triangle.
    ///
    /// The notion of triangle link only makes sense in 3D, where the output
    /// \p linkId refers to a vertex identifier.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getTriangleLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleLinkNumber()
    inline int getTriangleLink(const SimplexId &triangleId,
                               const int &localLinkId,
                               SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleLink(
        triangleId, localLinkId, linkId);
    }

    /// Get the number of simplices in the link of the \p triangleId-th
    /// triangle.
    ///
    /// The notion of triangle link only makes sense in 3D, where the number
    /// of vertices in the link will be returned.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of simplices in the link of the triangle.
    inline SimplexId
      getTriangleLinkNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleLinkNumber(triangleId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of link simplices for all triangles.
    ///
    /// The number of entries in this list is equal to the number of
    /// triangles.
    /// Each entry is a std::vector of identifiers representing a vertex.
    ///
    /// The notion of triangle link only makes sense in 3D.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle link list.
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinks() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getTriangleLinks();
    }

    /// Get the \p localStarId-th cell of the star of the \p triangleId-th
    /// triangle.
    ///
    /// The notion of triangle star only makes sense in 3D, where the output
    /// \p starId refers to a tetrahedron identifier.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param triangleId Input global triangle identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getTriangleStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleStarNumber()
    inline int getTriangleStar(const SimplexId &triangleId,
                               const int &localStarId,
                               SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleStar(
        triangleId, localStarId, starId);
    }

    /// Get the number of star cells for the \p triangleId-th triangle.
    ///
    /// The notion of triangle star only makes sense in 3D, where the number
    /// of tetrahedra in the star will be returned.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of star cells.
    inline SimplexId
      getTriangleStarNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleStarNumber(triangleId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of star cell identifiers for all triangles.
    ///
    /// The number of entries in this list is equal to the number of
    /// triangles.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of star cells for the corresponding triangle.
    ///
    /// The notion of triangle star only makes sense in 3D.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangleStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle star list.
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getTriangleStars();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p triangleId-th
    /// triangle.
    ///
    /// In 2D, this function is equivalent to getCellVertex().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global edge identifier.
    /// \param localVertexId Input local vertex identifier (in [0, 2]).
    /// \param vertexId Output global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertex()
    inline int getTriangleVertex(const SimplexId &triangleId,
                                 const int &localVertexId,
                                 SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleVertex(
        triangleId, localVertexId, vertexId);
    }

    /// Get the type of internal representation for the triangulation
    /// (explicit, implicit, periodic).
    ///
    /// \return Returns the current type of the triangulation.
    /// \sa setPeriodicBoundaryConditions()
    inline Triangulation::Type getType() const {
      if(abstractTriangulation_ == &explicitTriangulation_)
        return Triangulation::Type::EXPLICIT;
      else if(abstractTriangulation_ == &implicitTriangulation_)
        return Triangulation::Type::IMPLICIT;
      else if(abstractTriangulation_ == &implicitPreconditionsTriangulation_)
        return Triangulation::Type::HYBRID_IMPLICIT;
      else if(abstractTriangulation_ == &compactTriangulation_)
        return Triangulation::Type::COMPACT;
      else if(abstractTriangulation_ == &periodicImplicitTriangulation_)
        return Triangulation::Type::PERIODIC;
      else
        return Triangulation::Type::HYBRID_PERIODIC;
    }

    /// Get the \p localEdgeId-th edge identifier connected to the
    /// \p vertexId-th
    /// vertex.
    ///
    /// In 1D, this function is equivalent to getVertexStar().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getVertexEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexEdgeNumber()
    /// \sa getVertexStar()
    inline int getVertexEdge(const SimplexId &vertexId,
                             const int &localEdgeId,
                             SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexEdge(
        vertexId, localEdgeId, edgeId);
    }

    /// Get the number of edges connected to the \p vertexId-th vertex.
    ///
    /// In 1D, this function is equivalent to getVertexStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of edges connected to the vertex.
    /// \sa getVertexStarNumber()
    inline SimplexId
      getVertexEdgeNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexEdgeNumber(vertexId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of edge identifiers for all vertices.
    ///
    /// The number of entries in this list is equal to the number of
    /// vertices.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of edges connected to the corresponding vertex.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 1D, this function is equivalent to getVertexStars()
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex edge list.
    /// \sa getVertexStars()
    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getVertexStars();

      return abstractTriangulation_->getVertexEdges();
    }

#ifdef TTK_ENABLE_MPI

    /// Get the corresponding global id for a given local id of a vertex.
    ///    ///
    /// \pre For this function to behave correctly,
    /// preconditionDistributedVertices() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param leid Input local vertex identifier.
    /// \return vertexId Input global vertex identifier.
    inline SimplexId getVertexGlobalId(const SimplexId lvid) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexGlobalId(lvid);
    }

    /// Get the global id to local id map for the triangulation.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionDistributedVertices() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param map the std::unordered_map<SimplexId, SimplexId> in which we want
    /// our GidToLidMap. \return 0 if successful, -1 else.
    inline std::unordered_map<SimplexId, SimplexId> &getVertexGlobalIdMap() {
      return this->explicitTriangulation_.getVertexGlobalIdMap();
    }

    /// Set the flag for precondtioning of distributed vertices of the
    /// triangulation.
    inline void setHasPreconditionedDistributedVertices(bool flag) override {
      abstractTriangulation_->setHasPreconditionedDistributedVertices(flag);
    }

    inline bool hasPreconditionedDistributedCells() const override {
      return abstractTriangulation_->hasPreconditionedDistributedCells();
    }
    inline bool hasPreconditionedDistributedVertices() const override {
      return abstractTriangulation_->hasPreconditionedDistributedVertices();
    }

    inline const std::vector<int> &getNeighborRanks() const override {
      return abstractTriangulation_->getNeighborRanks();
    }
    inline std::vector<int> &getNeighborRanks() override {
      return abstractTriangulation_->getNeighborRanks();
    }

    inline const std::vector<std::array<ttk::SimplexId, 6>> &
      getNeighborVertexBBoxes() const override {
      return abstractTriangulation_->getNeighborVertexBBoxes();
    }

    inline const std::vector<std::vector<SimplexId>> &
      getGhostCellsPerOwner() const override {
      return abstractTriangulation_->getGhostCellsPerOwner();
    }

    inline const std::vector<std::vector<SimplexId>> &
      getRemoteGhostCells() const override {
      return abstractTriangulation_->getRemoteGhostCells();
    }

    inline int getVertexRank(const SimplexId lvid) const override {
      return this->abstractTriangulation_->getVertexRank(lvid);
    }

    /// Get the corresponding local id for a given global id of a vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionDistributedVertices() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param geid Input global vertex identifier.
    /// \return vertexId Input local vertex identifier.
    inline SimplexId getVertexLocalId(const SimplexId geid) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexLocalId(geid);
    }

    /**
     * @brief Create a meta grid for implicit triangulations
     *
     * In an MPI context, input domains are split into separate,
     * overlapping local grids. This methods makes each of the local
     * implicit triangulations aware of the dimensions of the
     * original, global grid.
     *
     * @param[in] dimensions Global grid dimensions
     */
    inline void createMetaGrid(const double *const bounds) {
      this->implicitPreconditionsTriangulation_.createMetaGrid(bounds);
      this->periodicImplicitTriangulation_.createMetaGrid(bounds);
      this->implicitTriangulation_.createMetaGrid(bounds);
      this->periodicPreconditionsTriangulation_.createMetaGrid(bounds);
      // also pass bounding box to ExplicitTriangulation...
      this->explicitTriangulation_.setBoundingBox(bounds);
    }

    inline void setIsBoundaryPeriodic(std::array<unsigned char, 6> boundary) {
      this->periodicImplicitTriangulation_.setIsBoundaryPeriodic(boundary);
      this->periodicPreconditionsTriangulation_.setIsBoundaryPeriodic(boundary);
    }

    /**
     * @brief  Get the Global Id of the simplex by calling the appropriate
     * global id retrieval function based on the simplex dimension cellDim
     *
     * @param localCellId: local id of the simplex
     * @param cellDim: dimension of the simplex
     * @param globalCellId global id of the simplex
     */
    inline int getDistributedGlobalCellId(const ttk::SimplexId &localCellId,
                                          const int &cellDim,
                                          ttk::SimplexId &globalCellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return this->abstractTriangulation_->getDistributedGlobalCellId(
        localCellId, cellDim, globalCellId);
    }

#endif // TTK_ENABLE_MPI

    /// Get the \p localLinkId-th simplex of the link of the \p vertexId-th
    /// vertex.
    ///
    /// The output \p linkId refers in 2D to an edge identifier and in 3D to
    /// a triangle identifier.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getVertexLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexLinkNumber()
    inline int getVertexLink(const SimplexId &vertexId,
                             const int &localLinkId,
                             SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexLink(
        vertexId, localLinkId, linkId);
    }

    /// Get the number of simplices in the link of the \p vertexId-th vertex.
    ///
    /// In 2D, this will return the number of edges in the link, in 3D the
    /// number of triangles.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of cells in the link of the vertex.
    inline SimplexId
      getVertexLinkNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexLinkNumber(vertexId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of link simplices for all vertices.
    ///
    /// The number of entries in this list is equal to the number of
    /// vertices.
    /// Each entry is a std::vector of identifiers representing edges in 2D
    /// and
    /// triangles in 3D.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex link list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getVertexLinks();
    }

    /// Get the \p localNeighborId-th vertex neighbor of the \p vertexId-th
    /// vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localNeighborId Input local neighbor identifier,
    /// in [0, getVertexNeighborNumber()].
    /// \param neighborId Output global neighbor vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexNeighborNumber()
    inline int getVertexNeighbor(const SimplexId &vertexId,
                                 const int &localNeighborId,
                                 SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      neighborId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexNeighbor(
        vertexId, localNeighborId, neighborId);
    }

    /// Get the number of vertex neighbors for the \p vertexId-th vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number vertex neighbors.
    inline SimplexId
      getVertexNeighborNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexNeighborNumber(vertexId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of vertex neighbor identifiers for all vertices.
    ///
    /// The number of entries in this list is equal to the number of
    /// vertices.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of vertex neighbors for the corresponding vertex.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex neighbor list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighbors() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getVertexNeighbors();
    }

    /// Get the point (3D coordinates) for the \p vertexId-th vertex.
    /// \param vertexId Input global vertex identifier.
    /// \param x Output x coordinate.
    /// \param y Output y coordinate.
    /// \param z Output z coordinate.
    /// \return Returns 0 upon success, negative values otherwise.
    inline int getVertexPoint(const SimplexId &vertexId,
                              float &x,
                              float &y,
                              float &z) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variables before early return
      x = NAN;
      y = NAN;
      z = NAN;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexPoint(vertexId, x, y, z);
    }

    /// Get the \p localStarId-th cell of the star of the \p vertexId-th
    /// vertex.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param vertexId Input global vertex identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getVertexStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexStarNumber()
    inline int getVertexStar(const SimplexId &vertexId,
                             const int &localStarId,
                             SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexStar(
        vertexId, localStarId, starId);
    }

    /// Get the number of star cells for the \p vertexId-th vertex.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier
    /// \return Returns the number of star cells.
    inline SimplexId
      getVertexStarNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexStarNumber(vertexId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of star cell identifiers for all vertices.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// The number of entries in this list is equal to the number of vertices.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of star cells for the corresponding vertex.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex star list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getVertexStars();
    }

    /// Get the \p localTriangleId-th triangle id of the
    /// \p vertexId-th vertex.
    ///
    /// In 2D, this function is equivalent to getVertexStar().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    ///
    /// \param vertexId Input global vertex identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getVertexTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexTriangleNumber()
    /// \sa getVertexStar()
    inline int getVertexTriangle(const SimplexId &vertexId,
                                 const int &localTriangleId,
                                 SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexTriangle(
        vertexId, localTriangleId, triangleId);
    }

    /// Get the number of triangles for the \p vertexId-th vertex.
    ///
    /// In 2D, this function is equivalent to getVertexStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of vertex triangles.
    /// \sa getVertexStarNumber()
    inline SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getVertexTriangleNumber(vertexId);
    }

    /// \warning
    /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
    /// DOING.
    ///
    /// Get the list of triangles for all vertices.
    ///
    /// The number of entries in this list is equal to the number of vertices.
    /// Each entry is a std::vector of identifiers whose size is equal to the
    /// number of triangles for the corresponding vertex.
    ///
    /// In implicit mode, this function will force the creation of such a
    /// list (which will be time and memory consuming).
    /// THIS IS USUALLY A BAD IDEA.
    ///
    /// In 2D, this function is equivalent to getVertexStars().
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex triangle list.
    /// \sa getVertexStars()
    inline const std::vector<std::vector<SimplexId>> *
      getVertexTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return nullptr;
#endif
      return abstractTriangulation_->getVertexTriangles();
    }

    /// Check if the edge with global identifier \p edgeId is on the boundary
    /// of the domain.
    ///
    /// For 2D triangulations, this function will return true if the edge is
    /// a boundary edge. For 3D triangulations, this function will return
    /// true if the edge belongs to a boundary triangle.
    ///
    /// Here the notion of edge only makes sense if the triangulation
    /// has a dimension greater than 1 (otherwise, use the cell information).
    ///
    /// \pre For this function to behave correctly,
    /// preconditionBoundaryEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns true if the edge is on the boundary, false otherwise.
    inline bool isEdgeOnBoundary(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;
#endif
      return abstractTriangulation_->isEdgeOnBoundary(edgeId);
    }

    /// Check if the data structure is empty or not.
    /// \return Returns true if empty, false otherwise.
    inline bool isEmpty() const override {
      return !abstractTriangulation_;
    }

    /// Check if the triangulation is manifold or not (Rips Complexes
    /// are not manifold)
    /// \return True if the triangulation is manifold
    inline bool isManifold() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(this->isEmptyCheck()) {
        return true;
      }
#endif // TTK_ENABLE_KAMIKAZE
      return this->abstractTriangulation_->isManifold();
    }

    /// Check if the triangulation is manifold or not.
    ///
    /// \ref ttk::ExplicitTriangulation (and maybe \ref
    /// ttk::CompactTriangulation too) can be generated from
    /// non-manifold datasets (such as a Rips Complex). Some TTK
    /// modules may be valid only for manifold triangulations, other
    /// may have alternatives for non-manifold data-sets (\see
    /// ttk::PersistenceDiagram::checkManifold).
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following function(s):
    ///   - isManifold()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isManifold()
    inline int preconditionManifold() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(this->isEmptyCheck())
        return false;
#endif // TTK_ENABLE_KAMIKAZE
      return this->abstractTriangulation_->preconditionManifold();
    }

    /// Check if the triangle with global identifier \p triangleId is on the
    /// boundary of the domain.
    ///
    /// Here the notion of triangle only makes sense if the triangulation
    /// has a dimension greater than 2 (otherwise, use the cell information).
    ///
    /// For 2D triangulations, this function will return false all the time.
    /// For 3D triangulations, this function will return true if the triangle
    /// is a boundary triangle.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionBoundaryTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns true if the triangle is on the boundary, false
    /// otherwise.
    inline bool
      isTriangleOnBoundary(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;
#endif
      return abstractTriangulation_->isTriangleOnBoundary(triangleId);
    }

    /// Check if the vertex with global identifier \p vertexId is on the
    /// boundary of the domain.
    ///
    /// For 2D triangulations, this function will return true if the vertex
    /// belongs to a boundary edge. For 3D triangulations, this function will
    /// return true if the vertex belongs to a boundary triangle.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionBoundaryVertices() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a preconditioning step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns true if the vertex is on the boundary, false
    /// otherwise.
    inline bool isVertexOnBoundary(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;
#endif
      return abstractTriangulation_->isVertexOnBoundary(vertexId);
    }

    /// Pre-process the boundary edges.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following function(s):
    ///   - isEdgeOnBoundary()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isEdgeOnBoundary()
    inline int preconditionBoundaryEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionBoundaryEdges();
    }

    /// Pre-process the boundary triangles.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following function(s):
    ///   - isTriangleOnBoundary()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isTriangleOnBoundary()
    inline int preconditionBoundaryTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionBoundaryTriangles();
    }

    /// Pre-process the boundary vertices.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following function(s):
    ///   - isVertexOnBoundary()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isVertexOnBoundary()
    inline int preconditionBoundaryVertices() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionBoundaryVertices();
    }

    /// Pre-process the cell edges.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getCellEdge()
    ///   - getCellEdgeNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellEdge()
    /// \sa getCellEdgeNumber()
    inline int preconditionCellEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionCellEdges();
    }

    /// Pre-process the cell neighbors.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getCellNeighbor()
    ///   - getCellNeighbors()
    ///   - getCellNeighborNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellNeighbor()
    /// \sa getCellNeighbors()
    /// \sa getCellNeighborNumber()
    inline int preconditionCellNeighbors() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionCellNeighbors();
    }

    /// Pre-process the cell triangles.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getCellTriangle()
    ///   - getCellTriangles()
    ///   - getCellTriangleNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellTriangle()
    /// \sa getCellTriangles()
    /// \sa getCellTriangleNumber()
    inline int preconditionCellTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionCellTriangles();
    }

    /// Pre-process the edges.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getEdges()
    ///   - getEdgeVertex()
    ///   - getNumberOfEdges()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdges()
    /// \sa getEdgeVertex()
    /// \sa getNumberOfEdges()
    inline int preconditionEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionEdges();
    }

    /// Pre-process the edge links.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getEdgeLink()
    ///   - getEdgeLinks()
    ///   - getEdgeLinkNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeLink()
    /// \sa getEdgeLinks()
    /// \sa getEdgeLinkNumber()
    inline int preconditionEdgeLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionEdgeLinks();
    }

    /// Pre-process the edge stars.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getEdgeStar()
    ///   - getEdgeStars()
    ///   - getEdgeStarNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeStar()
    /// \sa getEdgeStars()
    /// \sa getEdgeStarNumber()
    inline int preconditionEdgeStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionEdgeStars();
    }

    /// Pre-process the edge triangles.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getEdgeTriangle()
    ///   - getEdgeTriangles()
    ///   - getEdgeTriangleNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeTriangle()
    /// \sa getEdgeTriangles()
    /// \sa getEdgeTriangleNumber()
    inline int preconditionEdgeTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionEdgeTriangles();
    }

    /// Pre-process the triangles.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getNumberOfTriangles()
    ///   - getTriangles()
    ///   - getTriangleVertex()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getNumberOfTriangles()
    /// \sa getTriangles()
    /// \sa getTriangleVertex()
    inline int preconditionTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionTriangles();
    }

    /// Pre-process the triangle edges.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getTriangleEdge()
    ///   - getTriangleEdges()
    ///   - getTriangleEdgeNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleEdge()
    /// \sa getTriangleEdges()
    /// \sa getTriangleEdgeNumber()
    inline int preconditionTriangleEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionTriangleEdges();
    }

    /// Pre-process the triangle links.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getTriangleLink()
    ///   - getTriangleLinks()
    ///   - getTriangleLinkNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleLink()
    /// \sa getTriangleLinks()
    /// \sa getTriangleLinkNumber()
    inline int preconditionTriangleLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionTriangleLinks();
    }

    /// Pre-process the triangle stars.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getTriangleStar()
    ///   - getTriangleStars()
    ///   - getTriangleStarNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleStar()
    /// \sa getTriangleStars()
    /// \sa getTriangleStarNumber()
    inline int preconditionTriangleStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionTriangleStars();
    }

    /// Pre-process the vertex edges.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexEdge()
    ///   - getVertexEdges()
    ///   - getVertexEdgeNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexEdge()
    /// \sa getVertexEdges()
    /// \sa getVertexEdgeNumber()
    inline int preconditionVertexEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionVertexEdges();
    }

    ///
    /// \pre Get the VTK id of a cell using the TTK id of a cell.
    /// For implicit triangulations, VTK cells are squares or cubes,
    /// whereas TTK cells are triangles and tetrahedron.
    /// For other triangulations, the TTK and VTK cells coincide.
    ///
    /// \param ttkId TTK cell id
    /// \param vtkId Output VTK cell id
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    inline int getCellVTKID(const int &ttkId, int &vtkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->getCellVTKID(ttkId, vtkId);
    }

#ifdef TTK_ENABLE_MPI
    /// Pre-process the distributed vertex ids.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexGlobalId()
    ///   - getVertexGlobalIdMap()
    ///   - getVertexLocalId()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexGlobalId()
    /// \sa getVertexGlobalIdMap()
    /// \sa getVertexLocalId()
    inline int preconditionDistributedVertices() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionDistributedVertices();
    }

    inline int preconditionEdgeRankArray() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionEdgeRankArray();
    }

    inline int preconditionTriangleRankArray() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionTriangleRankArray();
    }

    inline int setVertexRankArray(const int *rankArray) override {
      return abstractTriangulation_->setVertexRankArray(rankArray);
    }

    inline int setCellRankArray(const int *rankArray) override {
      return abstractTriangulation_->setCellRankArray(rankArray);
    }

    /// Pre-process the global boundaries when using MPI. Local bounds should
    /// be set prior to using this function.
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa globalBounds_

    inline int preconditionGlobalBoundary() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionGlobalBoundary();
    }
    /// Pre-process the distributed ghost cells .
    ///
    /// This function should ONLY be called as a pre-condition for
    /// handling cells when executing with MPI.
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    inline int preconditionDistributedCells() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionDistributedCells();
    }

    /// Pre-process the distributed ghost cells .
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getGhostCellsPerOwner()
    ///   - getRemoteGhostCells()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getGhostCellsPerOwner()
    /// \sa getRemoteGhostCells()
    inline int preconditionExchangeGhostCells() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionExchangeGhostCells();
    }

    /// Pre-process the distributed ghost vertices .
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getGhostVerticesPerOwner()
    ///   - getRemoteGhostVertices()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getGhostVerticesPerOwner()
    /// \sa getRemoteGhostVertices()
    inline int preconditionExchangeGhostVertices() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionExchangeGhostVertices();
    }
#endif // TTK_ENABLE_MPI

    /// Pre-process the vertex links.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexLink()
    ///   - getVertexLinks()
    ///   - getVertexLinkNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexLink()
    /// \sa getVertexLinks()
    /// \sa getVertexLinkNumber()
    inline int preconditionVertexLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionVertexLinks();
    }

    /// Pre-process the vertex neighbors.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexNeighbor()
    ///   - getVertexNeighbors()
    ///   - getVertexNeighborNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexNeighbor()
    /// \sa getVertexNeighbors()
    /// \sa getVertexNeighborNumber()
    inline int preconditionVertexNeighbors() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionVertexNeighbors();
    }

    /// Pre-process the vertex stars.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexStar()
    ///   - getVertexStars()
    ///   - getVertexStarNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexStar()
    /// \sa getVertexStars()
    /// \sa getVertexStarNumber()
    inline int preconditionVertexStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->preconditionVertexStars();
    }

    /// Pre-process the vertex triangles.
    ///
    /// This function should ONLY be called as a pre-condition to the
    /// following functions:
    ///   - getVertexTriangle()
    ///   - getVertexTriangles()
    ///   - getVertexTriangleNumber()
    ///
    /// \pre This function should be called prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude this preconditioning function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexTriangle()
    /// \sa getVertexTriangles()
    /// \sa getVertexTriangleNumber()
    inline int preconditionVertexTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->preconditionVertexTriangles();
    }

    /// Tune the debug level (default: 0)
    inline int setDebugLevel(const int &debugLevel) override {
      explicitTriangulation_.setDebugLevel(debugLevel);
      compactTriangulation_.setDebugLevel(debugLevel);
      implicitTriangulation_.setDebugLevel(debugLevel);
      implicitPreconditionsTriangulation_.setDebugLevel(debugLevel);
      periodicImplicitTriangulation_.setDebugLevel(debugLevel);
      periodicPreconditionsTriangulation_.setDebugLevel(debugLevel);
      debugLevel_ = debugLevel;
      return 0;
    }

    // Set the cache size
    inline int setCacheSize(const float &ratio) {
      if(abstractTriangulation_ == &compactTriangulation_) {
        compactTriangulation_.initCache(ratio);
      }
      return 0;
    }

#ifdef TTK_CELL_ARRAY_NEW
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \param cellNumber Number of input cells.
    /// \param connectivity Pointer to an array of long long int. It contains
    /// the list of point ids of each cell.
    /// \param offsets Pointer to an array of long long int. It has a size of
    /// cellNumber+1 and each cell contains the position of the first vertex id
    /// of this cell in the connectivity array.
    /// This corresponds to the default cell array representation in VTK 9.
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that preconditioning functions should be called again).
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *connectivity,
                             const LongSimplexId *offset) {
      abstractTriangulation_ = &explicitTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return explicitTriangulation_.setInputCells(
        cellNumber, connectivity, offset);
    }

    inline int setStellarInputCells(const SimplexId &cellNumber,
                                    const LongSimplexId *connectivity,
                                    const LongSimplexId *offset) {
      abstractTriangulation_ = &compactTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return compactTriangulation_.setInputCells(
        cellNumber, connectivity, offset);
    }
#else
    /// Set the input cells for the triangulation.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \param cellNumber Number of input cells.
    /// \param cellArray Pointer to the input cells. This pointer should point
    /// to an array of long long int where cells are stored one after the
    /// other. In particular, each cell starts by the number of vertices in
    /// it, followed by the identifiers of its vertices. This corresponds to
    /// the default cell array representation in VTK < 9.
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that preconditioning functions should be called again).
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {
      abstractTriangulation_ = &explicitTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return explicitTriangulation_.setInputCells(cellNumber, cellArray);
    }

    inline int setStellarInputCells(const SimplexId &cellNumber,
                                    const LongSimplexId *cellArray) {

      abstractTriangulation_ = &compactTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;

      return compactTriangulation_.setInputCells(cellNumber, cellArray);
    }
#endif
    /// Set the specifications of the input grid to implicitly represent as a
    /// triangulation.
    /// \param xOrigin Input x coordinate of the grid origin.
    /// \param yOrigin Input y coordinate of the grid origin.
    /// \param zOrigin Input z coordinate of the grid origin.
    /// \param xSpacing Input spacing along the x dimension.
    /// \param ySpacing Input spacing along the y dimension.
    /// \param zSpacing Input spacing along the z dimension.
    /// \param xDim Input number of vertices along the x dimension.
    /// \param yDim Input number of vertices along the y dimension.
    /// \param zDim Input number of vertices along the z dimension.
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that preconditioning functions should be called again).
    inline int setInputGrid(const float &xOrigin,
                            const float &yOrigin,
                            const float &zOrigin,
                            const float &xSpacing,
                            const float &ySpacing,
                            const float &zSpacing,
                            const SimplexId &xDim,
                            const SimplexId &yDim,
                            const SimplexId &zDim) {
      gridDimensions_[0] = xDim;
      gridDimensions_[1] = yDim;
      gridDimensions_[2] = zDim;

      int ret{};

      ret |= periodicImplicitTriangulation_.setInputGrid(
        xOrigin, yOrigin, zOrigin, xSpacing, ySpacing, zSpacing, xDim, yDim,
        zDim);
      ret |= periodicPreconditionsTriangulation_.setInputGrid(
        xOrigin, yOrigin, zOrigin, xSpacing, ySpacing, zSpacing, xDim, yDim,
        zDim);
      ret |= implicitTriangulation_.setInputGrid(xOrigin, yOrigin, zOrigin,
                                                 xSpacing, ySpacing, zSpacing,
                                                 xDim, yDim, zDim);
      ret |= implicitPreconditionsTriangulation_.setInputGrid(
        xOrigin, yOrigin, zOrigin, xSpacing, ySpacing, zSpacing, xDim, yDim,
        zDim);
      const auto useImplicitPreconditions = this->processImplicitStrategy();

      this->switchGrid(this->hasPeriodicBoundaries_, useImplicitPreconditions);
      return ret;
    }

    /// Set the input grid to use period boundary conditions.
    ///
    /// \param usePeriodicBoundaries If this set to true then a triangulation
    /// with periodic boundaries will be used.
    inline void
      setPeriodicBoundaryConditions(const bool &usePeriodicBoundaries) {

      if(abstractTriangulation_ == &implicitTriangulation_
         || abstractTriangulation_ == &periodicImplicitTriangulation_
         || abstractTriangulation_ == &implicitPreconditionsTriangulation_
         || abstractTriangulation_ == &periodicPreconditionsTriangulation_) {
        if(usePeriodicBoundaries == hasPeriodicBoundaries_) {
          return;
        }
        const auto hasPreconditions{this->hasImplicitPreconditions()};

        this->switchGrid(usePeriodicBoundaries, hasPreconditions);

        // reset hasPreconditioned boolean
        AbstractTriangulation::clear();
        // but don't forget to set hasPeriodicBoundaries_
        hasPeriodicBoundaries_ = usePeriodicBoundaries;
      }
    }

    /**
     * @brief Set the input grid preconditioning strategy.
     * @param[in] strategy Strategy to implement.
     */
    inline void setImplicitPreconditions(const STRATEGY strategy) {
      if(abstractTriangulation_ == &implicitTriangulation_
         || abstractTriangulation_ == &periodicImplicitTriangulation_
         || abstractTriangulation_ == &implicitPreconditionsTriangulation_
         || abstractTriangulation_ == &periodicPreconditionsTriangulation_) {

        const auto useImplicitPreconditions
          = this->processImplicitStrategy(strategy);

        if(useImplicitPreconditions == this->hasImplicitPreconditions()) {
          return;
        }
        const auto isPeriodic{this->hasPeriodicBoundaries_};

        this->switchGrid(isPeriodic, useImplicitPreconditions);

        // reset hasPreconditioned boolean
        AbstractTriangulation::clear();
        // but don't forget to set hasImplicitPreconditions_
        hasPeriodicBoundaries_ = isPeriodic;
      }
    }

    /// Set the input 3D points of the triangulation.
    /// \param pointNumber Number of input vertices.
    /// \param pointSet Pointer to the 3D points. This pointer should point to
    /// an array of float where points are stored one after the other.
    /// In particular, each point is represented by X-Y-Z coordinates (one
    /// after the other). This corresponds to the default point set
    /// representation in VTK.
    /// \param doublePrecision Should we use double precision or stay
    /// with simple?
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that preconditioning functions should be called again).
    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const bool &doublePrecision = false) {

      abstractTriangulation_ = &explicitTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return explicitTriangulation_.setInputPoints(
        pointNumber, pointSet, doublePrecision);
    }

    inline int setStellarInputPoints(const SimplexId &pointNumber,
                                     const void *pointSet,
                                     const int *indexArray,
                                     const bool &doublePrecision = false) {

      abstractTriangulation_ = &compactTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return compactTriangulation_.setInputPoints(
        pointNumber, pointSet, indexArray, doublePrecision);
    }

    /// Tune the number of active threads (default: number of logical cores)
    inline int setThreadNumber(const ThreadId threadNumber) override {
      explicitTriangulation_.setThreadNumber(threadNumber);
      implicitTriangulation_.setThreadNumber(threadNumber);
      implicitPreconditionsTriangulation_.setThreadNumber(threadNumber);
      periodicImplicitTriangulation_.setThreadNumber(threadNumber);
      periodicPreconditionsTriangulation_.setThreadNumber(threadNumber);
      compactTriangulation_.setThreadNumber(threadNumber);
      threadNumber_ = threadNumber;
      return 0;
    }

    /// Internal usage. Pass the execution context (debug level, number of
    /// threads, etc.) to the implementing classes.
    inline int setWrapper(const Wrapper *wrapper) override {
      explicitTriangulation_.setWrapper(wrapper);
      implicitTriangulation_.setWrapper(wrapper);
      implicitPreconditionsTriangulation_.setWrapper(wrapper);
      periodicImplicitTriangulation_.setWrapper(wrapper);
      periodicPreconditionsTriangulation_.setWrapper(wrapper);
      compactTriangulation_.setWrapper(wrapper);
      return 0;
    }

    /// Returns true if the grid uses preconditions.
    inline bool hasImplicitPreconditions() const {
      return abstractTriangulation_ == &implicitPreconditionsTriangulation_
             || abstractTriangulation_ == &periodicPreconditionsTriangulation_;
    }

#ifdef TTK_ENABLE_MPI

    // GlobalPointIds, GlobalCellIds (only for ExplicitTriangulation)

    inline void setVertsGlobalIds(const LongSimplexId *data) {
      this->explicitTriangulation_.setVertsGlobalIds(data);
    }
    inline void setCellsGlobalIds(const LongSimplexId *const data) {
      this->explicitTriangulation_.setCellsGlobalIds(data);
    }

    // "vtkGhostType" on vertices & cells

    inline void setVertexGhostArray(const unsigned char *data) {
      this->abstractTriangulation_->setVertexGhostArray(data);
    }
    inline void setCellGhostArray(const unsigned char *data) {
      this->abstractTriangulation_->setCellGhostArray(data);
    }

    inline bool getIsMPIValid() const {
      return isMPIValid_;
    }

    inline void setIsMPIValid(bool flag) {
      isMPIValid_ = flag;
    }
#endif // TTK_ENABLE_MPI

  protected:
#ifdef TTK_ENABLE_MPI
    bool isMPIValid_{true};
#endif
    inline bool isEmptyCheck() const {
      if(!abstractTriangulation_) {
        printErr("Trying to access an empty data-structure!");
        return true;
      }
      return false;
    }

    /**
     * @brief Should we precondition the implicit/periodic triangulations?
     *
     * @param[in] strategy The strategy to follow
     * @return True to perform preconditioning.
     */
    bool processImplicitStrategy(const STRATEGY strategy
                                 = STRATEGY::DEFAULT) const;

    /**
     * @brief Switch regular grid triangulation type
     *
     * @param[in] usePeriodic Use a periodic triangulation
     * @param[in] usePreconditions Perform triangulation preconditions
     */
    void switchGrid(const bool usePeriodic, const bool usePreconditions);

    AbstractTriangulation *abstractTriangulation_;
    ExplicitTriangulation explicitTriangulation_;
    ImplicitNoPreconditions implicitTriangulation_;
    ImplicitWithPreconditions implicitPreconditionsTriangulation_;
    PeriodicNoPreconditions periodicImplicitTriangulation_;
    PeriodicWithPreconditions periodicPreconditionsTriangulation_;
    CompactTriangulation compactTriangulation_;
  };
} // namespace ttk
