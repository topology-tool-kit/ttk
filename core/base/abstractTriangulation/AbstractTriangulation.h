/// \ingroup base
/// \class ttk::AbstractTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief AbstractTriangulation is an interface class that defines an interface
/// for efficient traversal methods on triangulations of piecewise linear
/// manifolds.
///
/// \sa Triangulation

#ifndef _ABSTRACTTRIANGULATION_H
#define _ABSTRACTTRIANGULATION_H

// base code includes
#include <Geometry.h>
#include <Wrapper.h>

#include <array>
#include <ostream>

#ifdef TTK_ENABLE_KAMIKAZE
#define TTK_TRIANGULATION_INTERNAL(NAME) NAME
#else
#define TTK_TRIANGULATION_INTERNAL(NAME) NAME##Internal
#endif // TTK_ENABLE_KAMIKAZE

#define ttkTemplateMacroCase(triangulationType, triangulationClass, call) \
  case triangulationType: {                                               \
    using TTK_TT = triangulationClass;                                    \
    call;                                                                 \
  }; break

#define ttkTemplateMacro(triangulationType, call)                            \
  switch(triangulationType) {                                                \
    ttkTemplateMacroCase(                                                    \
      ttk::Triangulation::Type::EXPLICIT, ttk::ExplicitTriangulation, call); \
    ttkTemplateMacroCase(                                                    \
      ttk::Triangulation::Type::IMPLICIT, ttk::ImplicitTriangulation, call); \
    ttkTemplateMacroCase(ttk::Triangulation::Type::PERIODIC,                 \
                         ttk::PeriodicImplicitTriangulation, call);          \
    ttkTemplateMacroCase(                                                    \
      ttk::Triangulation::Type::COMPACT, ttk::CompactTriangulation, call);   \
  }

namespace ttk {

  class AbstractTriangulation : public Wrapper {

  public:
    AbstractTriangulation();

    virtual ~AbstractTriangulation();

    AbstractTriangulation(const AbstractTriangulation &) = default;
    AbstractTriangulation(AbstractTriangulation &&) = default;
    AbstractTriangulation &operator=(const AbstractTriangulation &) = default;
    AbstractTriangulation &operator=(AbstractTriangulation &&) = default;

    /// Reset the triangulation data-structures.
    /// \return Returns 0 upon success, negative values otherwise.
    int clear();

    /// Computes and displays the memory footprint of the data-structure.
    /// \return Returns 0 upon success, negative values otherwise.
    size_t footprint(size_t size = 0) const;

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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getCellEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellEdgeNumber()
    /// \sa getCellNeighbor()
    virtual inline int getCellEdge(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(!hasPreconditionedCellEdges())
        return -1;
#endif

      if(getDimensionality() == 1)
        return getCellNeighbor(cellId, localEdgeId, edgeId);

      else if(getDimensionality() == 2)
        return getTriangleEdgeInternal(cellId, localEdgeId, edgeId);

      return getCellEdgeInternal(cellId, localEdgeId, edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell edges.
    /// \sa getCellNeighborNumber()
    virtual inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellEdges())
        return -1;
#endif
      if(getDimensionality() == 1)
        return getCellNeighborNumber(cellId);

      else if(getDimensionality() == 2)
        return getTriangleEdgeNumber(cellId);

      return getCellEdgeNumberInternal(cellId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell edge list.
    /// \sa getCellNeighbors()
    virtual inline const std::vector<std::vector<SimplexId>> *getCellEdges() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellEdges())
        return NULL;
#endif
      if(getDimensionality() == 1)
        return getCellNeighbors();

      else if(getDimensionality() == 2)
        return getTriangleEdgesInternal();

      return getCellEdgesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localNeighborId Input local neighbor identifier,
    /// in [0, getCellNeighborNumber()].
    /// \param neighborId Output global neighbor cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellNeighborNumber()
    virtual inline int getCellNeighbor(const SimplexId &cellId,
                                       const int &localNeighborId,
                                       SimplexId &neighborId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      neighborId = -1;

      if(!hasPreconditionedCellNeighbors())
        return -1;
#endif
      return getCellNeighborInternal(cellId, localNeighborId, neighborId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell neighbors.
    virtual inline SimplexId
      getCellNeighborNumber(const SimplexId &cellId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellNeighbors())
        return -1;
#endif
      return getCellNeighborNumberInternal(cellId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell neighbor list.
    virtual inline const std::vector<std::vector<SimplexId>> *
      getCellNeighbors() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellNeighbors())
        return NULL;
#endif
      return getCellNeighborsInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param cellId Input global cell identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getCellTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellTriangleNumber()
    /// \sa getCellNeighbor()
    virtual inline int getCellTriangle(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedCellTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getCellNeighbor(cellId, localTriangleId, triangleId);

      return getCellTriangleInternal(cellId, localTriangleId, triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell triangles.
    /// \sa getCellNeighborNumber()
    virtual inline SimplexId
      getCellTriangleNumber(const SimplexId &cellId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedCellTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getCellNeighborNumber(cellId);

      return getCellTriangleNumberInternal(cellId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell triangle list.
    /// \sa getCellNeighbors()
    virtual inline const std::vector<std::vector<SimplexId>> *
      getCellTriangles() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedCellTriangles())
        return NULL;
#endif
      if(getDimensionality() == 2)
        return getCellNeighbors();

      return getCellTrianglesInternal();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p cellId-th
    /// cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \param cellId Input global cell identifier.
    /// \param localVertexId Input local vertex identifier,
    /// in [0, getCellVertexNumber()].
    /// \param vertexId Ouput global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertexNumber()
    virtual inline int getCellVertex(const SimplexId &cellId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const {

      return getCellVertexInternal(cellId, localVertexId, vertexId);
    }

    /// Get the number of vertices in a cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \param cellId Input global cell identifier.
    /// \returns Number of vertices in the cell.
    virtual inline SimplexId
      getCellVertexNumber(const SimplexId &cellId) const {
      return getCellVertexNumberInternal(cellId);
    }

    /// Get the dimensionality of the triangulation (this value is equal to
    /// the dimension of the simplex with largest dimensionality).
    /// \return Returns the dimensionality of the triangulation.
    virtual inline int getDimensionality() const {
      return getDimensionalityInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge list.
    virtual inline const std::vector<std::array<SimplexId, 2>> *getEdges() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedEdges())
        return NULL;
#endif
      return getEdgesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getEdgeLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeLinkNumber()
    virtual inline int getEdgeLink(const SimplexId &edgeId,
                                   const int &localLinkId,
                                   SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeLinks())
        return -2;
#endif
      return getEdgeLinkInternal(edgeId, localLinkId, linkId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of cells in the link of the edge.
    virtual inline SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeLinks())
        return -2;
#endif
      return getEdgeLinkNumberInternal(edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge link list.
    virtual inline const std::vector<std::vector<SimplexId>> *getEdgeLinks() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedEdgeLinks())
        return NULL;
#endif
      return getEdgeLinksInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param edgeId Input global edge identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getEdgeStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeStarNumber()
    virtual inline int getEdgeStar(const SimplexId &edgeId,
                                   const int &localStarId,
                                   SimplexId &starId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeStars())
        return -2;
#endif
      return getEdgeStarInternal(edgeId, localStarId, starId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier
    /// \return Returns the number of star cells.
    virtual inline SimplexId getEdgeStarNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeStars())
        return -2;
#endif
      return getEdgeStarNumberInternal(edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge star list.
    virtual inline const std::vector<std::vector<SimplexId>> *getEdgeStars() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedEdgeStars())
        return NULL;
#endif
      return getEdgeStarsInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param edgeId Input global edge identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getEdgeTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeTriangleNumber()
    /// \sa getEdgeStar()
    virtual inline int getEdgeTriangle(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getEdgeStar(edgeId, localTriangleId, triangleId);

      return getEdgeTriangleInternal(edgeId, localTriangleId, triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of edge triangles.
    /// \sa getEdgeStarNumber
    virtual inline SimplexId
      getEdgeTriangleNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedEdgeTriangles())
        return -2;
#endif

      if(getDimensionality() == 2)
        return getEdgeStarNumber(edgeId);

      return getEdgeTriangleNumberInternal(edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge triangle list.
    /// \sa getEdgeStars
    virtual inline const std::vector<std::vector<SimplexId>> *
      getEdgeTriangles() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedEdgeTriangles())
        return NULL;
#endif

      if(getDimensionality() == 2)
        return getEdgeStars();

      return getEdgeTrianglesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \param localVertexId Input local vertex identifier (0 or 1).
    /// \param vertexId Output global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertex()
    virtual inline int getEdgeVertex(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(!hasPreconditionedEdges())
        return -2;
#endif
      if(getDimensionality() == 1)
        return getCellVertex(edgeId, localVertexId, vertexId);

      return getEdgeVertexInternal(edgeId, localVertexId, vertexId);
    }

    /// Get the number of vertices of a particular edge.
    /// Always returns 2.
    inline int getEdgeVertexNumber(const SimplexId ttkNotUsed(edgeId)) const {
      return 2;
    }

    /// Get the dimensions of the grid if the current object is the implicit
    /// triangulation of a regular grid.
    /// \param dimensions Vector that will be filled with the dimensions of
    /// the grid. This std::vector has 3 entries (first: x, second: y,
    /// third: z).
    /// \return Returns 0 upon success, negative values otherwise (for
    /// instance, if the object is not representing a regular grid).
    virtual inline int getGridDimensions(std::vector<int> &dimensions) {

      if((gridDimensions_[0] == -1) && (gridDimensions_[1] == -1)
         && (gridDimensions_[2] == -1)) {
        return -1;
      }

      dimensions.resize(3);
      dimensions[0] = gridDimensions_[0];
      dimensions[1] = gridDimensions_[1];
      dimensions[2] = gridDimensions_[2];

      return 0;
    }

    /// Get the number of cells in the triangulation.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    /// \return Returns the number of cells.
    virtual inline SimplexId getNumberOfCells() const {
      return getNumberOfCellsInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns the number of edges.
    /// \sa getNumberOfCells()
    virtual inline SimplexId getNumberOfEdges() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdges())
        return -2;
#endif
      if(getDimensionality() == 1)
        return getNumberOfCells();

      return getNumberOfEdgesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns the number of triangles.
    /// \sa getNumberOfCells()
    virtual inline SimplexId getNumberOfTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getNumberOfCells();

      return getNumberOfTrianglesInternal();
    }

    /// Get the number of vertices in the triangulation.
    /// \return Returns the number of vertices.
    virtual inline SimplexId getNumberOfVertices() const {
      return getNumberOfVerticesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle list.
    virtual inline const std::vector<std::array<SimplexId, 3>> *getTriangles() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedTriangles())
        return NULL;
#endif
      return getTrianglesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getTriangleEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleEdgeNumber()
    /// \sa getCellEdge()
    virtual inline int getTriangleEdge(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedTriangleEdges())
        return -2;
#endif

      return getTriangleEdgeInternal(triangleId, localEdgeId, edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of cells in the link of the triangle.
    /// \sa getCellEdgeNumber()
    virtual inline SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedTriangleEdges())
        return -2;
#endif

      return getTriangleEdgeNumberInternal(triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle edge list.
    /// \sa getCellEdges()
    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdges() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedTriangleEdges())
        return NULL;
#endif

      return getTriangleEdgesInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getTriangleLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleLinkNumber()
    virtual inline int getTriangleLink(const SimplexId &triangleId,
                                       const int &localLinkId,
                                       SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(getDimensionality() != 3)
        return -1;

      if(!hasPreconditionedTriangleLinks())
        return -2;
#endif
      return getTriangleLinkInternal(triangleId, localLinkId, linkId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of simplices in the link of the triangle.
    virtual inline SimplexId
      getTriangleLinkNumber(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() != 3)
        return -1;

      if(!hasPreconditionedTriangleLinks())
        return -2;
#endif
      return getTriangleLinkNumberInternal(triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle link list.
    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinks() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() != 3)
        return NULL;

      if(!hasPreconditionedTriangleLinks())
        return NULL;
#endif
      return getTriangleLinksInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param triangleId Input global triangle identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getTriangleStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleStarNumber()
    virtual inline int getTriangleStar(const SimplexId &triangleId,
                                       const int &localStarId,
                                       SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(getDimensionality() != 3)
        return -1;

      if(!hasPreconditionedTriangleStars())
        return -2;
#endif
      return getTriangleStarInternal(triangleId, localStarId, starId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of star cells.
    virtual inline SimplexId
      getTriangleStarNumber(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() != 3)
        return -1;

      if(!hasPreconditionedTriangleStars())
        return -2;
#endif
      return getTriangleStarNumberInternal(triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle star list.
    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleStars() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() != 3)
        return NULL;

      if(!hasPreconditionedTriangleStars())
        return NULL;
#endif
      return getTriangleStarsInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global edge identifier.
    /// \param localVertexId Input local vertex identifier (in [0, 2]).
    /// \param vertexId Output global vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellVertex()
    virtual inline int getTriangleVertex(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getCellVertex(triangleId, localVertexId, vertexId);

      return getTriangleVertexInternal(triangleId, localVertexId, vertexId);
    }

    /// Get the number of vertices of a particular triangle.
    /// Always returns 3.
    inline int
      getTriangleVertexNumber(const SimplexId ttkNotUsed(triangleId)) const {
      return 3;
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localEdgeId Input local edge identifier,
    /// in [0, getVertexEdgeNumber()].
    /// \param edgeId Output global edge identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexEdgeNumber()
    /// \sa getVertexStar()
    virtual inline int getVertexEdge(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(!hasPreconditionedVertexEdges())
        return -1;
#endif
      if(getDimensionality() == 1)
        return getVertexStar(vertexId, localEdgeId, edgeId);

      return getVertexEdgeInternal(vertexId, localEdgeId, edgeId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of edges connected to the vertex.
    /// \sa getVertexStarNumber()
    virtual inline SimplexId
      getVertexEdgeNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexEdges())
        return -1;
#endif
      if(getDimensionality() == 1)
        return getVertexStarNumber(vertexId);

      return getVertexEdgeNumberInternal(vertexId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex edge list.
    /// \sa getVertexStars()
    virtual inline const std::vector<std::vector<SimplexId>> *getVertexEdges() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexEdges())
        return NULL;
#endif
      if(getDimensionality() == 1)
        return getVertexStars();

      return getVertexEdgesInternal();
    }

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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localLinkId Input local link simplex identifier,
    /// in [0, getVertexLinkNumber()].
    /// \param linkId Output link simplex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexLinkNumber()
    virtual inline int getVertexLink(const SimplexId &vertexId,
                                     const int &localLinkId,
                                     SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(!hasPreconditionedVertexLinks())
        return -1;
#endif
      return getVertexLinkInternal(vertexId, localLinkId, linkId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of cells in the link of the vertex.
    virtual inline SimplexId
      getVertexLinkNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexLinks())
        return -1;
#endif
      return getVertexLinkNumberInternal(vertexId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex link list.
    virtual inline const std::vector<std::vector<SimplexId>> *getVertexLinks() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexLinks())
        return NULL;
#endif
      return getVertexLinksInternal();
    }

    /// Get the \p localNeighborId-th vertex neighbor of the \p vertexId-th
    /// vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \param localNeighborId Input local neighbor identifier,
    /// in [0, getVertexNeighborNumber()].
    /// \param neighborId Output global neighbor vertex identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexNeighborNumber()
    virtual inline int getVertexNeighbor(const SimplexId &vertexId,
                                         const int &localNeighborId,
                                         SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      neighborId = -1;

      if(!hasPreconditionedVertexNeighbors())
        return -1;
#endif
      return getVertexNeighborInternal(vertexId, localNeighborId, neighborId);
    }

    /// Get the number of vertex neighbors for the \p vertexId-th vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preconditionVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number vertex neighbors.
    virtual inline SimplexId
      getVertexNeighborNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexNeighbors())
        return -1;
#endif
      return getVertexNeighborNumberInternal(vertexId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex neighbor list.
    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighbors() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexNeighbors())
        return NULL;
#endif
      return getVertexNeighborsInternal();
    }

    /// Get the point (3D coordinates) for the \p vertexId-th vertex.
    /// \param vertexId Input global vertex identifier.
    /// \param x Output x coordinate.
    /// \param y Output y coordinate.
    /// \param z Output z coordinate.
    /// \return Returns 0 upon success, negative values otherwise.
    virtual inline int getVertexPoint(const SimplexId &vertexId,
                                      float &x,
                                      float &y,
                                      float &z) const {
      return getVertexPointInternal(vertexId, x, y, z);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param vertexId Input global vertex identifier.
    /// \param localStarId Input local star cell identifier,
    /// in [0, getVertexStarNumber()].
    /// \param starId Output global star cell identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexStarNumber()
    virtual inline int getVertexStar(const SimplexId &vertexId,
                                     const int &localStarId,
                                     SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(!hasPreconditionedVertexStars())
        return -1;
#endif
      return getVertexStarInternal(vertexId, localStarId, starId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier
    /// \return Returns the number of star cells.
    virtual inline SimplexId
      getVertexStarNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexStars())
        return -1;
#endif
      return getVertexStarNumberInternal(vertexId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex star list.
    virtual inline const std::vector<std::vector<SimplexId>> *getVertexStars() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexStars())
        return NULL;
#endif
      return getVertexStarsInternal();
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    ///
    /// \param vertexId Input global vertex identifier.
    /// \param localTriangleId Input local triangle identifier,
    /// in [0, getVertexTriangleNumber()].
    /// \param triangleId Output global triangle identifier.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexTriangleNumber()
    /// \sa getVertexStar()
    virtual inline int getVertexTriangle(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedVertexTriangles())
        return -2;
#endif
      if(getDimensionality() == 2)
        return getVertexStar(vertexId, localTriangleId, triangleId);

      return getVertexTriangleInternal(vertexId, localTriangleId, triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of vertex triangles.
    /// \sa getVertexStarNumber()
    virtual inline SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return -1;

      if(!hasPreconditionedVertexTriangles())
        return -2;
#endif

      if(getDimensionality() == 2)
        return getVertexStarNumber(vertexId);

      return getVertexTriangleNumberInternal(vertexId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex triangle list.
    /// \sa getVertexStars()
    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexTriangles() {
#ifndef TTK_ENABLE_KAMIKAZE
      if(getDimensionality() == 1)
        return NULL;

      if(!hasPreconditionedVertexTriangles())
        return NULL;
#endif
      if(getDimensionality() == 2)
        return getVertexStars();

      return getVertexTrianglesInternal();
    }

    /// Returns true if the grid uses period boundary conditions.
    inline bool hasPeriodicBoundaries() const {
      return hasPeriodicBoundaries_;
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns true if the edge is on the boundary, false otherwise.
    virtual inline bool isEdgeOnBoundary(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryEdges())
        return false;
#endif
      return isEdgeOnBoundaryInternal(edgeId);
    }

    /// Check if the data structure is empty or not.
    /// \return Returns true if empty, false otherwise.
    virtual inline bool isEmpty() const {
      return true;
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns true if the triangle is on the boundary, false
    /// otherwise.
    virtual inline bool
      isTriangleOnBoundary(const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryTriangles())
        return false;
#endif
      return isTriangleOnBoundaryInternal(triangleId);
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
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns true if the vertex is on the boundary, false
    /// otherwise.
    virtual inline bool isVertexOnBoundary(const SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryVertices())
        return false;
#endif
      return isVertexOnBoundaryInternal(vertexId);
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isEdgeOnBoundary()
    virtual inline int preconditionBoundaryEdges() {

      if(!hasPreconditionedBoundaryEdges_) {
        preconditionEdges();
        preconditionBoundaryEdgesInternal();
        hasPreconditionedBoundaryEdges_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isTriangleOnBoundary()
    virtual inline int preconditionBoundaryTriangles() {

      if(!hasPreconditionedBoundaryTriangles_) {
        preconditionTriangles();
        preconditionBoundaryTrianglesInternal();
        hasPreconditionedBoundaryTriangles_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isVertexOnBoundary()
    virtual inline int preconditionBoundaryVertices() {

      if(!hasPreconditionedBoundaryVertices_) {
        preconditionBoundaryVerticesInternal();
        hasPreconditionedBoundaryVertices_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellEdge()
    /// \sa getCellEdgeNumber()
    virtual inline int preconditionCellEdges() {

      if(!hasPreconditionedCellEdges_) {
        hasPreconditionedCellEdges_ = true;
        if(getDimensionality() == 1)
          return preconditionCellNeighbors();

        preconditionEdges();
        preconditionCellEdgesInternal();
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellNeighbor()
    /// \sa getCellNeighbors()
    /// \sa getCellNeighborNumber()
    virtual inline int preconditionCellNeighbors() {

      if(!hasPreconditionedCellNeighbors_) {
        preconditionCellNeighborsInternal();
        hasPreconditionedCellNeighbors_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getCellTriangle()
    /// \sa getCellTriangles()
    /// \sa getCellTriangleNumber()
    virtual inline int preconditionCellTriangles() {

      if(!hasPreconditionedCellTriangles_) {
        hasPreconditionedCellTriangles_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif
        if(getDimensionality() == 2)
          return preconditionCellNeighbors();

        preconditionTriangles();
        preconditionCellTrianglesInternal();
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdges()
    /// \sa getEdgeVertex()
    /// \sa getNumberOfEdges()
    virtual inline int preconditionEdges() {

      if(!hasPreconditionedEdges_) {
        preconditionEdgesInternal();
        hasPreconditionedEdges_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeLink()
    /// \sa getEdgeLinks()
    /// \sa getEdgeLinkNumber()
    virtual inline int preconditionEdgeLinks() {

      if(!hasPreconditionedEdgeLinks_) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif
        preconditionEdges();
        preconditionEdgeLinksInternal();
        hasPreconditionedEdgeLinks_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeStar()
    /// \sa getEdgeStars()
    /// \sa getEdgeStarNumber()
    virtual inline int preconditionEdgeStars() {

      if(!hasPreconditionedEdgeStars_) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif
        preconditionEdges();
        preconditionEdgeStarsInternal();
        hasPreconditionedEdgeStars_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getEdgeTriangle()
    /// \sa getEdgeTriangles()
    /// \sa getEdgeTriangleNumber()
    virtual inline int preconditionEdgeTriangles() {

      if(!hasPreconditionedEdgeTriangles_) {
        hasPreconditionedEdgeTriangles_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif

        if(getDimensionality() == 2)
          return preconditionEdgeStars();

        preconditionEdges();
        preconditionTriangles();
        preconditionEdgeTrianglesInternal();
      }

      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getNumberOfTriangles()
    /// \sa getTriangles()
    /// \sa getTriangleVertex()
    virtual inline int preconditionTriangles() {

      if(!hasPreconditionedTriangles_) {
        hasPreconditionedTriangles_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif
        if(getDimensionality() == 2)
          return 0;

        preconditionTrianglesInternal();
      }

      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleEdge()
    /// \sa getTriangleEdges()
    /// \sa getTriangleEdgeNumber()
    virtual inline int preconditionTriangleEdges() {

      if(!hasPreconditionedTriangleEdges_) {
        hasPreconditionedTriangleEdges_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif
        if(getDimensionality() == 2)
          return preconditionCellEdges();

        preconditionEdges();
        preconditionTriangles();
        preconditionTriangleEdgesInternal();
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleLink()
    /// \sa getTriangleLinks()
    /// \sa getTriangleLinkNumber()
    virtual inline int preconditionTriangleLinks() {

      if(!hasPreconditionedTriangleLinks_) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() != 3)
          return -2;
#endif
        preconditionTriangles();
        preconditionTriangleLinksInternal();
        hasPreconditionedTriangleLinks_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getTriangleStar()
    /// \sa getTriangleStars()
    /// \sa getTriangleStarNumber()
    virtual inline int preconditionTriangleStars() {

      if(!hasPreconditionedTriangleStars_) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() != 3)
          return -1;
#endif
        preconditionTriangles();
        preconditionTriangleStarsInternal();
        hasPreconditionedTriangleStars_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexEdge()
    /// \sa getVertexEdges()
    /// \sa getVertexEdgeNumber()
    virtual inline int preconditionVertexEdges() {

      if(!hasPreconditionedVertexEdges_) {
        hasPreconditionedVertexEdges_ = true;

        if(getDimensionality() == 1)
          return preconditionVertexStars();

        preconditionEdges();
        preconditionVertexEdgesInternal();
      }
      return 0;
    }

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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexLink()
    /// \sa getVertexLinks()
    /// \sa getVertexLinkNumber()
    virtual inline int preconditionVertexLinks() {

      if(!hasPreconditionedVertexLinks_) {
        preconditionVertexLinksInternal();
        hasPreconditionedVertexLinks_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexNeighbor()
    /// \sa getVertexNeighbors()
    /// \sa getVertexNeighborNumber()
    virtual inline int preconditionVertexNeighbors() {

      if(!hasPreconditionedVertexNeighbors_) {
        preconditionVertexNeighborsInternal();
        hasPreconditionedVertexNeighbors_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexStar()
    /// \sa getVertexStars()
    /// \sa getVertexStarNumber()
    virtual inline int preconditionVertexStars() {

      if(!hasPreconditionedVertexStars_) {
        preconditionVertexStarsInternal();
        hasPreconditionedVertexStars_ = true;
      }
      return 0;
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa getVertexTriangle()
    /// \sa getVertexTriangles()
    /// \sa getVertexTriangleNumber()
    virtual inline int preconditionVertexTriangles() {

      if(!hasPreconditionedVertexTriangles_) {
        hasPreconditionedVertexTriangles_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
        if(getDimensionality() == 1)
          return -1;
#endif

        if(getDimensionality() == 2)
          return preconditionVertexStars();

        preconditionTriangles();
        preconditionVertexTrianglesInternal();
      }
      return 0;
    }

    /**
     * Compute the barycenter of the points of the given edge identifier.
     */
    inline int getEdgeIncenter(const SimplexId edgeId,
                               float incenter[3]) const {
      std::array<SimplexId, 2> vertexId{};
      getEdgeVertex(edgeId, 0, vertexId[0]);
      getEdgeVertex(edgeId, 1, vertexId[1]);

      std::array<float, 6> p{};
      getVertexPoint(vertexId[0], p[0], p[1], p[2]);
      getVertexPoint(vertexId[1], p[3], p[4], p[5]);

      incenter[0] = 0.5 * p[0] + 0.5 * p[3];
      incenter[1] = 0.5 * p[1] + 0.5 * p[4];
      incenter[2] = 0.5 * p[2] + 0.5 * p[5];

      return 0;
    }

    /**
     * Compute the incenter of the points of the given triangle identifier.
     */
    inline int getTriangleIncenter(const SimplexId triangleId,
                                   float incenter[3]) const {

      std::array<SimplexId, 3> vertexId{};
      if(getDimensionality() == 2) {
        getCellVertex(triangleId, 0, vertexId[0]);
        getCellVertex(triangleId, 1, vertexId[1]);
        getCellVertex(triangleId, 2, vertexId[2]);
      } else if(getDimensionality() == 3) {
        getTriangleVertex(triangleId, 0, vertexId[0]);
        getTriangleVertex(triangleId, 1, vertexId[1]);
        getTriangleVertex(triangleId, 2, vertexId[2]);
      }

      std::array<float, 9> p{};
      getVertexPoint(vertexId[0], p[0], p[1], p[2]);
      getVertexPoint(vertexId[1], p[3], p[4], p[5]);
      getVertexPoint(vertexId[2], p[6], p[7], p[8]);

      std::array<float, 3> d{};
      d[0] = Geometry::distance(&p[3], &p[6]);
      d[1] = Geometry::distance(&p[0], &p[6]);
      d[2] = Geometry::distance(&p[0], &p[3]);
      const float sum = d[0] + d[1] + d[2];

      d[0] = d[0] / sum;
      d[1] = d[1] / sum;
      d[2] = d[2] / sum;

      incenter[0] = d[0] * p[0] + d[1] * p[3] + d[2] * p[6];
      incenter[1] = d[0] * p[1] + d[1] * p[4] + d[2] * p[7];
      incenter[2] = d[0] * p[2] + d[1] * p[5] + d[2] * p[8];

      return 0;
    }

    /**
     * Compute the barycenter of the incenters of the triangles of the given
       tetra identifier.
     */
    inline int getTetraIncenter(const SimplexId tetraId,
                                float incenter[3]) const {
      incenter[0] = 0.0f;
      incenter[1] = 0.0f;
      incenter[2] = 0.0f;

      std::array<float, 3> p{};
      for(int i = 0; i < 4; ++i) {
        SimplexId triangleId;
        getCellTriangle(tetraId, i, triangleId);
        getTriangleIncenter(triangleId, p.data());
        incenter[0] += p[0];
        incenter[1] += p[1];
        incenter[2] += p[2];
      }

      incenter[0] /= 4.0f;
      incenter[1] /= 4.0f;
      incenter[2] /= 4.0f;

      return 0;
    }

    /**
     * Compute the geometric barycenter of a given cell.
     */
    inline int getCellIncenter(const SimplexId cellid,
                               const int dim,
                               float incenter[3]) const {
      switch(dim) {
        case 0:
          getVertexPoint(cellid, incenter[0], incenter[1], incenter[2]);
          break;
        case 1:
          getEdgeIncenter(cellid, incenter);
          break;
        case 2:
          getTriangleIncenter(cellid, incenter);
          break;
        case 3:
          getTetraIncenter(cellid, incenter);
          break;
      }
      return 0;
    }

  protected:
    virtual int getCellEdgeInternal(const SimplexId &ttkNotUsed(cellId),
                                    const int &ttkNotUsed(localEdgeId),
                                    SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &ttkNotUsed(cellId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() {
      return NULL;
    }

    virtual inline int
      getCellNeighborInternal(const SimplexId &ttkNotUsed(cellId),
                              const int &ttkNotUsed(localNeighborId),
                              SimplexId &ttkNotUsed(neighborId)) const {
      return 0;
    }

    virtual inline SimplexId
      getCellNeighborNumberInternal(const SimplexId &ttkNotUsed(cellId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getCellNeighborsInternal() {
      return NULL;
    }

    virtual inline int
      getCellTriangleInternal(const SimplexId &ttkNotUsed(cellId),
                              const int &ttkNotUsed(localTriangleId),
                              SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &ttkNotUsed(cellId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() {
      return NULL;
    }

    virtual inline int
      getCellVertexInternal(const SimplexId &ttkNotUsed(cellId),
                            const int &ttkNotUsed(localVertexId),
                            SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline SimplexId
      getCellVertexNumberInternal(const SimplexId &ttkNotUsed(cellId)) const {
      return 0;
    }

    virtual inline int getDimensionalityInternal() const {
      return 0;
    }

    virtual inline const std::vector<std::array<SimplexId, 2>> *
      getEdgesInternal() {
      return NULL;
    }

    virtual inline int
      getEdgeLinkInternal(const SimplexId &ttkNotUsed(edgeId),
                          const int &ttkNotUsed(localLinkId),
                          SimplexId &ttkNotUsed(linkId)) const {
      return 0;
    }

    virtual inline SimplexId
      getEdgeLinkNumberInternal(const SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getEdgeLinksInternal() {
      return NULL;
    }

    virtual inline int
      getEdgeStarInternal(const SimplexId &ttkNotUsed(edgeId),
                          const int &ttkNotUsed(localStarId),
                          SimplexId &ttkNotUsed(starId)) const {
      return 0;
    }

    virtual inline SimplexId
      getEdgeStarNumberInternal(const SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getEdgeStarsInternal() {
      return NULL;
    }

    virtual inline int
      getEdgeTriangleInternal(const SimplexId &ttkNotUsed(edgeId),
                              const int &ttkNotUsed(localTriangleId),
                              SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() {
      return NULL;
    }

    virtual inline int
      getEdgeVertexInternal(const SimplexId &ttkNotUsed(edgeId),
                            const int &ttkNotUsed(localVertexId),
                            SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline SimplexId getNumberOfCellsInternal() const {
      return 0;
    }

    virtual inline SimplexId getNumberOfEdgesInternal() const {
      return 0;
    }

    virtual inline SimplexId getNumberOfTrianglesInternal() const {
      return 0;
    }

    virtual inline SimplexId getNumberOfVerticesInternal() const {
      return 0;
    }

    virtual inline const std::vector<std::array<SimplexId, 3>> *
      getTrianglesInternal() {
      return NULL;
    }

    virtual inline int
      getTriangleEdgeInternal(const SimplexId &ttkNotUsed(triangleId),
                              const int &ttkNotUsed(localEdgeId),
                              SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() {
      return NULL;
    }

    virtual inline int
      getTriangleLinkInternal(const SimplexId &ttkNotUsed(triangleId),
                              const int &ttkNotUsed(localLinkId),
                              SimplexId &ttkNotUsed(linkId)) const {
      return 0;
    }

    virtual inline SimplexId getTriangleLinkNumberInternal(
      const SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinksInternal() {
      return NULL;
    }

    virtual inline int
      getTriangleStarInternal(const SimplexId &ttkNotUsed(triangleId),
                              const int &ttkNotUsed(localStarId),
                              SimplexId &ttkNotUsed(starId)) const {
      return 0;
    }

    virtual inline SimplexId getTriangleStarNumberInternal(
      const SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getTriangleStarsInternal() {
      return NULL;
    }

    virtual inline int
      getTriangleVertexInternal(const SimplexId &ttkNotUsed(triangleId),
                                const int &ttkNotUsed(localVertexId),
                                SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline int
      getVertexEdgeInternal(const SimplexId &ttkNotUsed(vertexId),
                            const int &ttkNotUsed(localEdgeId),
                            SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() {
      return NULL;
    }

    virtual inline int
      getVertexLinkInternal(const SimplexId &ttkNotUsed(vertexId),
                            const int &ttkNotUsed(localLinkId),
                            SimplexId &ttkNotUsed(linkId)) const {
      return 0;
    }

    virtual inline SimplexId
      getVertexLinkNumberInternal(const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexLinksInternal() {
      return NULL;
    }

    virtual inline int
      getVertexNeighborInternal(const SimplexId &ttkNotUsed(vertexId),
                                const int &ttkNotUsed(localNeighborId),
                                SimplexId &ttkNotUsed(neighborId)) const {
      return 0;
    }

    virtual inline SimplexId getVertexNeighborNumberInternal(
      const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighborsInternal() {
      return NULL;
    }

    virtual inline int
      getVertexPointInternal(const SimplexId &ttkNotUsed(vertexId),
                             float &ttkNotUsed(x),
                             float &ttkNotUsed(y),
                             float &ttkNotUsed(z)) const {
      return 0;
    }

    virtual inline int
      getVertexStarInternal(const SimplexId &ttkNotUsed(vertexId),
                            const int &ttkNotUsed(localStarId),
                            SimplexId &ttkNotUsed(starId)) const {
      return 0;
    }

    virtual inline SimplexId
      getVertexStarNumberInternal(const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexStarsInternal() {
      return NULL;
    }

    virtual inline int
      getVertexTriangleInternal(const SimplexId &ttkNotUsed(vertexId),
                                const int &ttkNotUsed(localTriangleId),
                                SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    virtual inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() {
      return NULL;
    }

    virtual inline bool
      isEdgeOnBoundaryInternal(const SimplexId &ttkNotUsed(edgeId)) const {
      return 0;
    }

    virtual inline bool isTriangleOnBoundaryInternal(
      const SimplexId &ttkNotUsed(triangleId)) const {
      return 0;
    }

    virtual inline bool
      isVertexOnBoundaryInternal(const SimplexId &ttkNotUsed(vertexId)) const {
      return 0;
    }

    inline bool hasPreconditionedBoundaryEdges() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryEdges_) {

        printMsg("BoundaryEdge query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionBoundaryEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif

      return hasPreconditionedBoundaryEdges_;
    }

    inline bool hasPreconditionedBoundaryTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryTriangles_) {

        printMsg("BoundaryTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg(
          "Please call preconditionBoundaryTriangles() in a pre-process.",
          debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedBoundaryTriangles_;
    }

    inline bool hasPreconditionedBoundaryVertices() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedBoundaryVertices_) {

        printMsg("BoundaryVertex query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionBoundaryVertices() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedBoundaryVertices_;
    }

    inline bool hasPreconditionedCellEdges() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 1) && (!hasPreconditionedCellNeighbors_))
         || ((getDimensionality() > 1) && (!hasPreconditionedCellEdges_))) {

        printMsg("CellEdge query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedCellNeighbors() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellNeighbors_) {

        printMsg("CellNeighbor query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellNeighbors() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedCellNeighbors_;
    }

    inline bool hasPreconditionedCellTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedCellNeighbors_))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedCellTriangles_))) {

        printMsg("CellTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedEdgeLinks() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdgeLinks_) {

        printMsg("EdgeLink query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeLinks() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedEdgeLinks_;
    }

    inline bool hasPreconditionedEdgeStars() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdgeStars_) {

        printMsg("EdgeStar query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeStars() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedEdgeStars_;
    }

    inline bool hasPreconditionedEdgeTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedEdgeStars_))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedEdgeTriangles_))) {

        printMsg("EdgeTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedEdges() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((getDimensionality() != 1) && (!hasPreconditionedEdges_)) {

        printMsg("Edge query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if((getDimensionality() == 3) && (!hasPreconditionedTriangles_)) {

        printMsg("Triangle query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedTriangleEdges() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedCellEdges_))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedTriangleEdges_))) {

        printMsg("TriangleEdge query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedTriangleLinks() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedTriangleLinks_) {

        printMsg("TriangleLink query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleLinks() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedTriangleLinks_;
    }

    inline bool hasPreconditionedTriangleStars() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedTriangleStars_) {

        printMsg("TriangleStar query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleStars() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedTriangleStars_;
    }

    inline bool hasPreconditionedVertexEdges() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 1) && (!hasPreconditionedVertexStars_))
         || ((getDimensionality() > 1) && (!hasPreconditionedVertexEdges_))) {

        printMsg("VertexEdge query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionVertexEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool hasPreconditionedVertexLinks() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexLinks_) {

        printMsg("VertexLink query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionVertexLinks() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedVertexLinks_;
    }

    inline bool hasPreconditionedVertexNeighbors() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexNeighbors_) {

        printMsg("VertexNeighbor query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionVertexNeighbors() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedVertexNeighbors_;
    }

    inline bool hasPreconditionedVertexStars() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexStars_) {

        printMsg("VertexStar query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionVertexStars() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
      }
#endif
      return hasPreconditionedVertexStars_;
    }

    inline bool hasPreconditionedVertexTriangles() const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedVertexStars_))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedVertexTriangles_))) {

        printMsg("VertexTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionVertexTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    // empty wrapping to VTK for now
    bool needsToAbort() {
      return false;
    }

    virtual inline int preconditionBoundaryEdgesInternal() {
      return 0;
    }

    virtual inline int preconditionBoundaryTrianglesInternal() {
      return 0;
    }

    virtual inline int preconditionBoundaryVerticesInternal() {
      return 0;
    }

    virtual inline int preconditionCellEdgesInternal() {
      return 0;
    }

    virtual inline int preconditionCellNeighborsInternal() {
      return 0;
    }

    virtual inline int preconditionCellTrianglesInternal() {
      return 0;
    }

    virtual inline int preconditionEdgesInternal() {
      return 0;
    }

    virtual inline int preconditionEdgeLinksInternal() {
      return 0;
    }

    virtual inline int preconditionEdgeStarsInternal() {
      return 0;
    }

    virtual inline int preconditionEdgeTrianglesInternal() {
      return 0;
    }

    virtual inline int preconditionTrianglesInternal() {
      return 0;
    }

    virtual inline int preconditionTriangleEdgesInternal() {
      return 0;
    }

    virtual inline int preconditionTriangleLinksInternal() {
      return 0;
    }

    virtual inline int preconditionTriangleStarsInternal() {
      return 0;
    }

    virtual inline int preconditionVertexEdgesInternal() {
      return 0;
    }

    virtual inline int preconditionVertexLinksInternal() {
      return 0;
    }

    virtual inline int preconditionVertexNeighborsInternal() {
      return 0;
    }

    virtual inline int preconditionVertexStarsInternal() {
      return 0;
    }

    virtual inline int preconditionVertexTrianglesInternal() {
      return 0;
    }

    template <class itemType>
    size_t tableFootprint(const std::vector<itemType> &table,
                          const std::string &tableName = "",
                          std::ostream &stream = std::cout) const {

      std::stringstream msg;
      if((table.size()) && (tableName.length())) {
        msg << tableName << ": " << table.size() * sizeof(itemType) << " bytes";
        printMsg(
          msg.str(), debug::Priority::INFO, debug::LineMode::NEW, stream);
      }

      return table.size() * sizeof(itemType);
    }

    template <class itemType>
    size_t tableTableFootprint(const std::vector<std::vector<itemType>> &table,
                               const std::string &tableName = "",
                               std::ostream &stream = std::cout) const;

    int updateProgress(const float &ttkNotUsed(progress)) {
      return 0;
    }

    bool hasPeriodicBoundaries_, hasPreconditionedBoundaryEdges_,
      hasPreconditionedBoundaryTriangles_, hasPreconditionedBoundaryVertices_,
      hasPreconditionedCellEdges_, hasPreconditionedCellNeighbors_,
      hasPreconditionedCellTriangles_, hasPreconditionedEdges_,
      hasPreconditionedEdgeLinks_, hasPreconditionedEdgeStars_,
      hasPreconditionedEdgeTriangles_, hasPreconditionedTriangles_,
      hasPreconditionedTriangleEdges_, hasPreconditionedTriangleLinks_,
      hasPreconditionedTriangleStars_, hasPreconditionedVertexEdges_,
      hasPreconditionedVertexLinks_, hasPreconditionedVertexNeighbors_,
      hasPreconditionedVertexStars_, hasPreconditionedVertexTriangles_;

    std::array<int, 3> gridDimensions_;

    std::vector<bool> boundaryEdges_, boundaryTriangles_, boundaryVertices_;
    std::vector<std::array<SimplexId, 6>> tetraEdgeList_;
    std::vector<std::vector<SimplexId>> cellNeighborList_;
    std::vector<std::array<SimplexId, 4>> tetraTriangleList_;
    std::vector<std::vector<SimplexId>> edgeLinkList_;
    std::vector<std::array<SimplexId, 2>> edgeList_;
    std::vector<std::vector<SimplexId>> edgeStarList_;
    std::vector<std::vector<SimplexId>> edgeTriangleList_;
    std::vector<std::array<SimplexId, 3>> triangleList_;
    std::vector<std::array<SimplexId, 3>> triangleEdgeList_;
    std::vector<std::vector<SimplexId>> triangleLinkList_;
    std::vector<std::vector<SimplexId>> triangleStarList_;
    std::vector<std::vector<SimplexId>> vertexEdgeList_;
    std::vector<std::vector<SimplexId>> vertexLinkList_;
    std::vector<std::vector<SimplexId>> vertexNeighborList_;
    std::vector<std::vector<SimplexId>> vertexStarList_;
    std::vector<std::vector<SimplexId>> vertexTriangleList_;

    // keep compatibility between getCellEdges(), getCellTriangles(),
    // getCellNeighbors() and getTriangleEdges()
    std::vector<std::vector<SimplexId>> cellEdgeVector_{};
    std::vector<std::vector<SimplexId>> cellTriangleVector_{};
    std::vector<std::vector<SimplexId>> triangleEdgeVector_{};
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <AbstractTriangulation.cpp>

#endif // _ABSTRACTTRIANGULATION_H
