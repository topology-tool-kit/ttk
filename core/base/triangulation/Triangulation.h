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
/// pre-processing, see the documentation furtherdown).
///   -# Implicit triangulations: Given a regular grid (origin, spacings and
/// dimensions), Triangulation will perform an implicit triangulation of the
/// grid, enabling both time and memory efficient traversals of triangulations
/// of regular grids.
///
/// Apart from pre-processes, Triangulation requires no memory overhead in
/// addition to the input data.
///
/// \note
/// Only pre-process the information you need! See the documentation further
/// down.
/// \sa ttkTriangulation

#ifndef _TRIANGULATION_H
#define _TRIANGULATION_H

// base code includes
#include <AbstractTriangulation.h>
#include <ExplicitTriangulation.h>
#include <ImplicitTriangulation.h>
#include <PeriodicImplicitTriangulation.h>

#include <array>

namespace ttk {

  class Triangulation final : public AbstractTriangulation {

  public:
    Triangulation();
    Triangulation(const Triangulation &);
    Triangulation(Triangulation &&);
    Triangulation &operator=(const Triangulation &);
    Triangulation &operator=(Triangulation &&);
    ~Triangulation();

    /// Reset the triangulation data-structures.
    /// \return Returns 0 upon success, negative values otherwise.
    inline int clear() override {

      if(abstractTriangulation_) {
        return abstractTriangulation_->clear();
      }

      return 0;
    }

    /// Computes and displays the memory footprint of the data-structure.
    /// \return Returns 0 upon success, negative values otherwise.
    inline size_t footprint() const override {

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
    /// preprocessCellEdges() needs to be called on this object prior to any
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
    inline int getCellEdge(const SimplexId &cellId,
                           const int &localEdgeId,
                           SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(isEmptyCheck())
        return -1;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedCellEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellEdge query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getCellNeighbor(
          cellId, localEdgeId, edgeId);

      return abstractTriangulation_->getCellEdge(cellId, localEdgeId, edgeId);
    }

    /// Get the number of edges for the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// In 1D, this function is equivalent to getCellNeighborNumber().
    ///
    /// \pre For this function to behave correctly, preprocessCellEdges()
    /// needs to be called on this object prior to any traversal, in a
    /// clearly distinct pre-processing step that involves no traversal at
    /// all. An error will be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell edges.
    /// \sa getCellNeighborNumber()
    inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedCellEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellEdgeNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getCellNeighborNumber(cellId);

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
    /// preprocessCellEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell edge list.
    /// \sa getCellNeighbors()
    inline const std::vector<std::vector<SimplexId>> *getCellEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedCellEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellEdges query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getCellNeighbors();

      return abstractTriangulation_->getCellEdges();
    }

    /// Get the \p localNeighborId-th cell neighbor of the \p cellId-th cell.
    ///
    /// Here the notion of cell refers to the simplicices of maximal
    /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
    ///
    /// \pre For this function to behave correctly,
    /// preprocessCellNeighbors() needs to be called
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
    inline int getCellNeighbor(const SimplexId &cellId,
                               const int &localNeighborId,
                               SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      neighborId = -1;

      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedCellNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellNeighbor query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessCellNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell neighbors.
    inline SimplexId
      getCellNeighborNumber(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedCellNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellNeighborNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessCellNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell neighbor list.
    inline const std::vector<std::vector<SimplexId>> *
      getCellNeighbors() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedCellNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellNeighbors query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
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
    /// preprocessCellTriangles() needs to be called
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
    inline int getCellTriangle(const SimplexId &cellId,
                               const int &localTriangleId,
                               SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedCellTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellTriangle query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellNeighbor(
          cellId, localTriangleId, triangleId);

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
    /// preprocessCellTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param cellId Input global cell identifier.
    /// \return Returns the number of cell triangles.
    /// \sa getCellNeighborNumber()
    inline SimplexId
      getCellTriangleNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedCellTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellTriangleNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellNeighborNumber(cellId);

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
    /// preprocessCellTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the cell triangle list.
    /// \sa getCellNeighbors()
    inline const std::vector<std::vector<SimplexId>> *
      getCellTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellNeighbors()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedCellTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "CellTriangles query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessCellTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellNeighbors();

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
    /// \param vertexId Ouput global vertex identifier.
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
    /// preprocessEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge list.
    inline const std::vector<std::pair<SimplexId, SimplexId>> *
      getEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedEdges()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "Edges query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getEdges();
    }

    /// Get the \p localLinkId-th simplex of the link of the \p edgeId-th
    /// edge.
    ///
    /// The output \p linkId refers in 2D to a vertex identifier and in 3D
    /// to an edge identifier. It returns a negative value in 1D.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessEdgeLinks() needs to be called
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
    inline int getEdgeLink(const SimplexId &edgeId,
                           const int &localLinkId,
                           SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedEdgeLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeLink query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      return abstractTriangulation_->getEdgeLink(edgeId, localLinkId, linkId);
    }

    /// Get the number of simplicies in the link of the \p edgeId-th edge.
    ///
    /// In 2D, this will return the number of vertices in the link, in 3D the
    /// number of edges. It returns a negative value in 1D.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessEdgeLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of cells in the link of the edge.
    inline SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedEdgeLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeLinkNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessEdgeLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge link list.
    inline const std::vector<std::vector<SimplexId>> *getEdgeLinks() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedEdgeLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeLinks query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
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
    /// preprocessEdgeStars() needs to be called
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
    inline int getEdgeStar(const SimplexId &edgeId,
                           const int &localStarId,
                           SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedEdgeStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeStar query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessEdgeStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier
    /// \return Returns the number of star cells.
    inline SimplexId getEdgeStarNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedEdgeStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeStarNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessEdgeStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge star list.
    inline const std::vector<std::vector<SimplexId>> *getEdgeStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedEdgeStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeStars query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getEdgeStars();
    }

    /// Get the \p localTriangleId-th triangle id of the \p edgeId-th edge.
    ///
    /// In 2D, this function is equivalent to getEdgeStar().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessEdgeTriangles() needs to be called
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
    inline int getEdgeTriangle(const SimplexId &edgeId,
                               const int &localTriangleId,
                               SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedEdgeStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedEdgeTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeTriangle query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getEdgeStar(
          edgeId, localTriangleId, triangleId);

      return abstractTriangulation_->getEdgeTriangle(
        edgeId, localTriangleId, triangleId);
    }

    /// Get the number of triangles for the \p edgeId-th edge.
    ///
    /// In 2D, this function is equivalent to getEdgeStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessEdgeTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns the number of edge triangles.
    /// \sa getEdgeStarNumber
    inline SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedEdgeStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedEdgeTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeTriangleNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif

      if(getDimensionality() == 2)
        return abstractTriangulation_->getEdgeStarNumber(edgeId);

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
    /// preprocessEdgeTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the edge triangle list.
    /// \sa getEdgeStars
    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedEdgeStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedEdgeTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeTriangles query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdgeTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif

      if(getDimensionality() == 2)
        return abstractTriangulation_->getEdgeStars();

      return abstractTriangulation_->getEdgeTriangles();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p edgeId-th
    /// edge.
    ///
    /// In 1D, this function is equivalent to getCellVertex().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessEdges() needs to be called
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
    inline int getEdgeVertex(const SimplexId &edgeId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(isEmptyCheck())
        return -1;

      if((getDimensionality() != 1)
         && (!abstractTriangulation_->hasPreprocessedEdges())) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "EdgeVertex query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getCellVertex(
          edgeId, localVertexId, vertexId);

      return abstractTriangulation_->getEdgeVertex(
        edgeId, localVertexId, vertexId);
    }

    /// Get the dimensions of the grid if the current object is the implicit
    /// triangulation of a regular grid.
    /// \param dimensions Vector that will be filled with the dimensions of
    /// the grid. This std::vector has 3 entries (first: x, second: y,
    /// third: z).
    /// \return Returns 0 upon success, negative values otherwise (for
    /// instance, if the object is not representing a regular grid).
    inline int getGridDimensions(std::vector<int> &dimensions) {

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
    /// preprocessEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns the number of edges.
    /// \sa getNumberOfCells()
    inline SimplexId getNumberOfEdges() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if((getDimensionality() != 1)
         && (!abstractTriangulation_->hasPreprocessedEdges())) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "NumberOfEdges query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getNumberOfCells();

      return abstractTriangulation_->getNumberOfEdges();
    }

    /// Get the number of triangles in the triangulation.
    ///
    /// In 2D, this function is equivalent to getNumberOfCells().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns the number of triangles.
    /// \sa getNumberOfCells()
    inline SimplexId getNumberOfTriangles() const override {
#ifndef TTK_ENABLE_KAMIKAZE

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if((getDimensionality() == 3)
         && (!abstractTriangulation_->hasPreprocessedTriangles())) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "NumberOfTriangles query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getNumberOfCells();

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
    /// preprocessTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle list.
    inline const std::vector<std::vector<SimplexId>> *getTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedTriangles()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "Triangles query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getTriangles();
    }

    /// Get the \p localEdgeId-th edge of the \p triangleId-th triangle.
    ///
    /// In 2D, this function is equivalent to getCellEdge().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessTriangleEdges() needs to be called
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
    inline int getTriangleEdge(const SimplexId &triangleId,
                               const int &localEdgeId,
                               SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellEdges()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedTriangleEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleEdge query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellEdge(
          triangleId, localEdgeId, edgeId);

      return abstractTriangulation_->getTriangleEdge(
        triangleId, localEdgeId, edgeId);
    }

    /// Get the number of edges of the \p triangleId-th triangle.
    ///
    /// In 2D, this function is equivalent to getCellEdgeNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessTriangleEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of cells in the link of the triangle.
    /// \sa getCellEdgeNumber()
    inline SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellEdges()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedTriangleEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleEdgeNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellEdgeNumber(triangleId);

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
    /// preprocessTriangleEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle edge list.
    /// \sa getCellEdges()
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedCellEdges()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedTriangleEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleEdges query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellEdges();

      return abstractTriangulation_->getTriangleEdges();
    }

    /// Get the \p localLinkId-th simplex of the link of the \p triangleId-th
    /// triangle.
    ///
    /// The notion of triangle link only makes sense in 3D, where the output
    /// \p linkId refers to a vertex identifier.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessTriangleLinks() needs to be called
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
    inline int getTriangleLink(const SimplexId &triangleId,
                               const int &localLinkId,
                               SimplexId &linkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedTriangleLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleLink query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -1;
      }
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
    /// preprocessTriangleLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of simplices in the link of the triangle.
    inline SimplexId
      getTriangleLinkNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedTriangleLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleLinkNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessTriangleLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle link list.
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleLinks() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() != 3)
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedTriangleLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleLinks query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
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
    /// preprocessTriangleStars() needs to be called
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
    inline int getTriangleStar(const SimplexId &triangleId,
                               const int &localStarId,
                               SimplexId &starId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedTriangleStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleStar query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessTriangleStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns the number of star cells.
    inline SimplexId
      getTriangleStarNumber(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;

      if(!abstractTriangulation_->hasPreprocessedTriangleStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleStarNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessTriangleStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the triangle star list.
    inline const std::vector<std::vector<SimplexId>> *
      getTriangleStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() != 3)
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedTriangleStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleStars query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangleStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getTriangleStars();
    }

    /// Get the \p localVertexId-th vertex identifier of the \p triangleId-th
    /// triangle.
    ///
    /// In 2D, this function is equivalent to getCellVertex().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessTriangles() needs to be called
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
    inline int getTriangleVertex(const SimplexId &triangleId,
                                 const int &localVertexId,
                                 SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if((getDimensionality() == 3)
         && (!abstractTriangulation_->hasPreprocessedTriangles())) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "TriangleVertex query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getCellVertex(
          triangleId, localVertexId, vertexId);

      return abstractTriangulation_->getTriangleVertex(
        triangleId, localVertexId, vertexId);
    }

    /// Get the \p localEdgeId-th edge identifier connected to the
    /// \p vertexId-th
    /// vertex.
    ///
    /// In 1D, this function is equivalent to getVertexStar().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexEdges() needs to be called
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
    inline int getVertexEdge(const SimplexId &vertexId,
                             const int &localEdgeId,
                             SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(isEmptyCheck())
        return -1;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedVertexEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexEdge query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getVertexStar(
          vertexId, localEdgeId, edgeId);

      return abstractTriangulation_->getVertexEdge(
        vertexId, localEdgeId, edgeId);
    }

    /// Get the number of edges connected to the \p vertexId-th vertex.
    ///
    /// In 1D, this function is equivalent to getVertexStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of edges connected to the vertex.
    /// \sa getVertexStarNumber()
    inline SimplexId
      getVertexEdgeNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedVertexEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexEdgeNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getVertexStarNumber(vertexId);

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
    /// preprocessVertexEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex edge list.
    /// \sa getVertexStars()
    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(((getDimensionality() == 1)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() > 1)
             && (!abstractTriangulation_->hasPreprocessedVertexEdges()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexEdges query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      if(getDimensionality() == 1)
        return abstractTriangulation_->getVertexStars();

      return abstractTriangulation_->getVertexEdges();
    }

    /// Get the \p localLinkId-th simplex of the link of the \p vertexId-th
    /// vertex.
    ///
    /// The output \p linkId refers in 2D to an edge identifier and in 3D to
    /// a triangle identifier.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexLinks() needs to be called
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
    inline int getVertexLink(const SimplexId &vertexId,
                             const int &localLinkId,
                             SimplexId &linkId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      linkId = -1;

      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexLink query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessVertexLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of cells in the link of the vertex.
    inline SimplexId
      getVertexLinkNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexLinkNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessVertexLinks() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex link list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedVertexLinks()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexLinks query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexLinks() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getVertexLinks();
    }

    /// Get the \p localNeighborId-th vertex neighbor of the \p vertexId-th
    /// vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexNeighbors() needs to be called
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
    inline int getVertexNeighbor(const SimplexId &vertexId,
                                 const int &localNeighborId,
                                 SimplexId &neighborId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      neighborId = -1;

      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexNeighbor query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      return abstractTriangulation_->getVertexNeighbor(
        vertexId, localNeighborId, neighborId);
    }

    /// Get the number of vertex neighbors for the \p vertexId-th vertex.
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number vertex neighbors.
    inline SimplexId
      getVertexNeighborNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexNeighborNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessVertexNeighbors() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex neighbor list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexNeighbors() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedVertexNeighbors()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexNeighbors query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexNeighbors() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
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
    /// preprocessVertexStars() needs to be called
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
    inline int getVertexStar(const SimplexId &vertexId,
                             const int &localStarId,
                             SimplexId &starId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      starId = -1;

      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexStar query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexStar() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessVertexStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier
    /// \return Returns the number of star cells.
    inline SimplexId
      getVertexStarNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(!abstractTriangulation_->hasPreprocessedVertexStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexStarNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
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
    /// preprocessVertexStars() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex star list.
    inline const std::vector<std::vector<SimplexId>> *
      getVertexStars() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(!abstractTriangulation_->hasPreprocessedVertexStars()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexStars query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexStars() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      return abstractTriangulation_->getVertexStars();
    }

    /// Get the \p localTriangleId-th triangle id of the
    /// \p vertexId-th vertex.
    ///
    /// In 2D, this function is equivalent to getVertexStar().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexTriangles() needs to be called
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
    inline int getVertexTriangle(const SimplexId &vertexId,
                                 const int &localTriangleId,
                                 SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedVertexTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexTriangle query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getVertexStar(
          vertexId, localTriangleId, triangleId);

      return abstractTriangulation_->getVertexTriangle(
        vertexId, localTriangleId, triangleId);
    }

    /// Get the number of triangles for the \p vertexId-th vertex.
    ///
    /// In 2D, this function is equivalent to getVertexStarNumber().
    ///
    /// \pre For this function to behave correctly,
    /// preprocessVertexTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns the number of vertex triangles.
    /// \sa getVertexStarNumber()
    inline SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedVertexTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexTriangleNumber query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return -2;
      }
#endif

      if(getDimensionality() == 2)
        return abstractTriangulation_->getVertexStarNumber(vertexId);

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
    /// preprocessVertexTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \return Returns a pointer to the vertex triangle list.
    /// \sa getVertexStars()
    inline const std::vector<std::vector<SimplexId>> *
      getVertexTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return NULL;

      if(getDimensionality() == 1)
        return NULL;

      if(((getDimensionality() == 2)
          && (!abstractTriangulation_->hasPreprocessedVertexStars()))
         || ((getDimensionality() == 3)
             && (!abstractTriangulation_->hasPreprocessedVertexTriangles()))) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "VertexTriangles query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessVertexTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return NULL;
      }
#endif
      if(getDimensionality() == 2)
        return abstractTriangulation_->getVertexStars();

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
    /// preprocessBoundaryEdges() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param edgeId Input global edge identifier.
    /// \return Returns true if the edge is on the boundary, false otherwise.
    inline bool isEdgeOnBoundary(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;

      if(!abstractTriangulation_->hasPreprocessedBoundaryEdges()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "BoundaryEdge query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessBoundaryEdges() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return false;
      }
#endif
      return abstractTriangulation_->isEdgeOnBoundary(edgeId);
    }

    /// Check if the data structure is empty or not.
    /// \return Returns true if empty, false otherwise.
    inline bool isEmpty() const override {
      return !abstractTriangulation_;
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
    /// preprocessBoundaryTriangles() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param triangleId Input global triangle identifier.
    /// \return Returns true if the triangle is on the boundary, false
    /// otherwise.
    inline bool
      isTriangleOnBoundary(const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;

      if(!abstractTriangulation_->hasPreprocessedBoundaryTriangles()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "BoundaryTriangle query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessBoundaryTriangles() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return false;
      }
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
    /// preprocessBoundaryVertices() needs to be called
    /// on this object prior to any traversal, in a clearly distinct
    /// pre-processing step that involves no traversal at all. An error will
    /// be returned otherwise.
    /// \note It is recommended to exclude such a pre-processing step
    /// from any time performance measurement.
    /// \param vertexId Input global vertex identifier.
    /// \return Returns true if the vertex is on the boundary, false
    /// otherwise.
    inline bool isVertexOnBoundary(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return false;

      if(!abstractTriangulation_->hasPreprocessedBoundaryVertices()) {
        std::stringstream msg;
        msg << "[Triangulation] "
            << "BoundaryVertex query without pre-process!" << std::endl;
        msg << "[Triangulation] "
            << "Please call preprocessBoundaryVertices() in a"
            << " pre-process." << std::endl;
        dMsg(std::cerr, msg.str(), Debug::fatalMsg);
        return false;
      }
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
    /// \note It is recommended to exclude this pre-processing function from
    /// any time performance measurement.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa isEdgeOnBoundary()
    inline int preprocessBoundaryEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessBoundaryEdges())
               && (hasPreprocessedBoundaryEdges_ = true));
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
    inline int preprocessBoundaryTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessBoundaryTriangles())
               && (hasPreprocessedBoundaryTriangles_ = true));
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
    inline int preprocessBoundaryVertices() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessBoundaryVertices())
               && (hasPreprocessedBoundaryVertices_ = true));
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
    inline int preprocessCellEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      if(getDimensionality() == 1)
        return !((!abstractTriangulation_->preprocessCellNeighbors())
                 && (hasPreprocessedCellEdges_ = true)
                 && (hasPreprocessedCellNeighbors_ = true));

      return !((!abstractTriangulation_->preprocessCellEdges())
               && (hasPreprocessedCellEdges_ = true));
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
    inline int preprocessCellNeighbors() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessCellNeighbors())
               && (hasPreprocessedCellNeighbors_ = true));
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
    inline int preprocessCellTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;
#endif
      if(getDimensionality() == 2)
        return !((!abstractTriangulation_->preprocessCellNeighbors())
                 && (hasPreprocessedCellTriangles_ = true)
                 && (hasPreprocessedCellNeighbors_ = true));

      return !((!abstractTriangulation_->preprocessCellTriangles())
               && (hasPreprocessedCellTriangles_ = true));
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
    inline int preprocessEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessEdges())
               && (hasPreprocessedBoundaryEdges_ = true));
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
    inline int preprocessEdgeLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;
#endif

      return !((!abstractTriangulation_->preprocessEdgeLinks())
               && (hasPreprocessedEdgeLinks_ = true));
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
    inline int preprocessEdgeStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;
#endif

      return !((!abstractTriangulation_->preprocessEdgeStars())
               && (hasPreprocessedEdgeStars_ = true));
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
    inline int preprocessEdgeTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;
#endif

      if(getDimensionality() == 2) {
        return !((!abstractTriangulation_->preprocessEdgeStars())
                 && (hasPreprocessedEdgeTriangles_ = true)
                 && (hasPreprocessedEdgeStars_ = true));
      }

      return !((!abstractTriangulation_->preprocessEdgeTriangles())
               && (hasPreprocessedEdgeTriangles_ = true));
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
    inline int preprocessTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -1;
#endif

      if(getDimensionality() == 2)
        return 0;

      return !((!abstractTriangulation_->preprocessTriangles())
               && (hasPreprocessedTriangles_ = true));
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
    inline int preprocessTriangleEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() == 1)
        return -2;
#endif

      if(getDimensionality() == 2)
        return !((!abstractTriangulation_->preprocessCellEdges())
                 && (hasPreprocessedTriangleEdges_ = true)
                 && (hasPreprocessedCellEdges_ = true));

      return !((!abstractTriangulation_->preprocessTriangleEdges())
               && (hasPreprocessedTriangleEdges_ = true));
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
    inline int preprocessTriangleLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;
#endif

      return !((!abstractTriangulation_->preprocessTriangleLinks())
               && (hasPreprocessedTriangleLinks_ = true));
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
    inline int preprocessTriangleStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;

      if(getDimensionality() != 3)
        return -2;
#endif

      return !((!abstractTriangulation_->preprocessTriangleStars())
               && (hasPreprocessedTriangleStars_ = true));
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
    inline int preprocessVertexEdges() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      if(getDimensionality() == 1)
        return !((!abstractTriangulation_->preprocessVertexStars())
                 && (hasPreprocessedVertexEdges_ = true)
                 && (hasPreprocessedVertexStars_ = true));

      return !((!abstractTriangulation_->preprocessVertexEdges())
               && (hasPreprocessedVertexEdges_ = true));
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
    inline int preprocessVertexLinks() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessVertexLinks())
               && (hasPreprocessedVertexLinks_ = true));
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
    inline int preprocessVertexNeighbors() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessVertexNeighbors())
               && (hasPreprocessedVertexNeighbors_ = true));
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
    inline int preprocessVertexStars() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return !((!abstractTriangulation_->preprocessVertexStars())
               && (hasPreprocessedVertexStars_ = true));
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
    inline int preprocessVertexTriangles() override {

#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
      if(getDimensionality() == 1)
        return -2;
#endif
      if(getDimensionality() == 2) {
        return !((!abstractTriangulation_->preprocessVertexStars())
                 && (hasPreprocessedVertexTriangles_ = true)
                 && (hasPreprocessedVertexStars_ = true));
      }

      return !((!abstractTriangulation_->preprocessVertexTriangles())
               && (hasPreprocessedVertexTriangles_ = true));
    }

    int getEdgeIncenter(SimplexId edgeId, float incenter[3]) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif

      return abstractTriangulation_->getEdgeIncenter(edgeId, incenter);
    };

    int getTriangleIncenter(SimplexId triangleId,
                            float incenter[3]) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTriangleIncenter(triangleId, incenter);
    };

    int getTetraIncenter(SimplexId tetraId, float incenter[3]) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(isEmptyCheck())
        return -1;
#endif
      return abstractTriangulation_->getTetraIncenter(tetraId, incenter);
    };

    /// Tune the debug level (default: 0)
    inline int setDebugLevel(const int &debugLevel) override {
      explicitTriangulation_.setDebugLevel(debugLevel);
      implicitTriangulation_.setDebugLevel(debugLevel);
      periodicImplicitTriangulation_.setDebugLevel(debugLevel);
      debugLevel_ = debugLevel;
      return 0;
    }

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
    /// the default cell array representation in VTK.
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that pre-processing functions should be called again).
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {

      abstractTriangulation_ = &explicitTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;

      return explicitTriangulation_.setInputCells(cellNumber, cellArray);
    }

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
    /// means that pre-processing functions should be called again).
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

      int retPeriodic = periodicImplicitTriangulation_.setInputGrid(
        xOrigin, yOrigin, zOrigin, xSpacing, ySpacing, zSpacing, xDim, yDim,
        zDim);
      int ret = implicitTriangulation_.setInputGrid(xOrigin, yOrigin, zOrigin,
                                                    xSpacing, ySpacing,
                                                    zSpacing, xDim, yDim, zDim);

      if(usePeriodicBoundaries_) {
        abstractTriangulation_ = &periodicImplicitTriangulation_;
        return retPeriodic;
      } else {
        abstractTriangulation_ = &implicitTriangulation_;
        return ret;
      }
      return 0;
    }

    /// Set the input grid to use period boundary conditions.
    ///
    /// \param usePeriodicBoundaries If this set to true then a triangulation
    /// with periodic boundaries will be used.
    inline void
      setPeriodicBoundaryConditions(const bool &usePeriodicBoundaries) {

      if((abstractTriangulation_ == &implicitTriangulation_)
         || (abstractTriangulation_ == &periodicImplicitTriangulation_)) {
        if(usePeriodicBoundaries == usePeriodicBoundaries_) {
          return;
        }
        usePeriodicBoundaries_ = usePeriodicBoundaries;
        if(usePeriodicBoundaries_) {
          abstractTriangulation_ = &periodicImplicitTriangulation_;
        } else {
          abstractTriangulation_ = &implicitTriangulation_;
        }
      }
    }

    /// Returns true if the grid uses period boundary conditions.
    /// \see setPeriodicBoundaryConditions
    bool usesPeriodicBoundaryConditions() const {
      return usePeriodicBoundaries_;
    }

    /// Set the input 3D points of the triangulation.
    /// \param pointNumber Number of input vertices.
    /// \param pointSet Pointer to the 3D points. This pointer should point to
    /// an array of float where points are stored one after the other.
    /// In particular, each point is represented by X-Y-Z coordinates (one
    /// after the other). This corresponds to the default point set
    /// representation in VTK.
    /// \return Returns 0 upon success, negative values otherwise.
    ///
    /// \note This function does not need to be called if the current object
    /// is a vtkTriangulation (this function is automatically called
    /// if needed through vtkTriangulation::setInputData()).
    ///
    /// \warning If this ttk::Triangulation object is already representing a
    /// valid triangulation, this information will be over-written (which
    /// means that pre-processing functions should be called again).
    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const bool &doublePrecision = false) {

      abstractTriangulation_ = &explicitTriangulation_;
      gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
      return explicitTriangulation_.setInputPoints(
        pointNumber, pointSet, doublePrecision);
    }

    /// Tune the number of active threads (default: number of logical cores)
    inline int setThreadNumber(const ThreadId &threadNumber) {
      explicitTriangulation_.setThreadNumber(threadNumber);
      implicitTriangulation_.setThreadNumber(threadNumber);
      periodicImplicitTriangulation_.setThreadNumber(threadNumber);
      threadNumber_ = threadNumber;
      return 0;
    }

    /// Internal usage. Pass the execution context (debug level, number of
    /// threads, etc.) to the implementing classes.
    inline int setWrapper(const Wrapper *wrapper) override {
      explicitTriangulation_.setWrapper(wrapper);
      implicitTriangulation_.setWrapper(wrapper);
      periodicImplicitTriangulation_.setWrapper(wrapper);
      return 0;
    }

  protected:
    inline bool isEmptyCheck() const {
      if(!abstractTriangulation_) {
        std::stringstream msg;
        msg << "[Triangulation] Trying to access an empty data-structure!"
            << std::endl;
        dMsg(std::cerr, msg.str(), fatalMsg);
        return true;
      }
      return false;
    }

    std::array<int, 3> gridDimensions_;

    AbstractTriangulation *abstractTriangulation_;
    ExplicitTriangulation explicitTriangulation_;
    ImplicitTriangulation implicitTriangulation_;
    PeriodicImplicitTriangulation periodicImplicitTriangulation_;
    bool usePeriodicBoundaries_;
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <Triangulation.cpp>

#endif // _TRIANGULATION_H
