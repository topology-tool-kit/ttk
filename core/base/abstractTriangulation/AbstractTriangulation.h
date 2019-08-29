/// \ingroup base
/// \class ttk::AbstractTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief AbstractTriangulation is a virtual class that defines an interface
/// for efficient traversal methods on triangulations of piecewise linear
/// manifolds.
///
/// \sa Triangulation

#ifndef _ABSTRACTTRIANGULATION_H
#define _ABSTRACTTRIANGULATION_H

// base code includes
#include <Wrapper.h>

namespace ttk {

  class AbstractTriangulation : public Wrapper {

  public:
    AbstractTriangulation();

    ~AbstractTriangulation();

    virtual int clear();

    virtual size_t footprint() const;

    virtual int getCellEdge(const SimplexId &cellId,
                            const int &localEdgeId,
                            SimplexId &edgeId) const = 0;

    virtual SimplexId getCellEdgeNumber(const SimplexId &cellId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getCellEdges() = 0;

    virtual int getCellNeighbor(const SimplexId &cellId,
                                const int &localNeighborId,
                                SimplexId &neighborId) const = 0;

    virtual SimplexId getCellNeighborNumber(const SimplexId &cellId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getCellNeighbors() = 0;

    virtual int getCellTriangle(const SimplexId &cellId,
                                const int &localTriangleId,
                                SimplexId &triangleId) const = 0;

    virtual SimplexId getCellTriangleNumber(const SimplexId &cellId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getCellTriangles() = 0;

    virtual int getCellVertex(const SimplexId &cellId,
                              const int &localVertexId,
                              SimplexId &vertexId) const = 0;

    virtual SimplexId getCellVertexNumber(const SimplexId &cellId) const = 0;

    virtual int getDimensionality() const = 0;

    virtual const std::vector<std::pair<SimplexId, SimplexId>> *getEdges() = 0;

    virtual int getEdgeLink(const SimplexId &edgeId,
                            const int &localLinkId,
                            SimplexId &linkId) const = 0;

    virtual SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getEdgeLinks() = 0;

    virtual int getEdgeStar(const SimplexId &edgeId,
                            const int &localStarId,
                            SimplexId &starId) const = 0;

    virtual SimplexId getEdgeStarNumber(const SimplexId &edgeId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getEdgeStars() = 0;

    virtual int getEdgeTriangle(const SimplexId &edgeId,
                                const int &localTriangleId,
                                SimplexId &triangleId) const = 0;

    virtual SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getEdgeTriangles() = 0;

    virtual int getEdgeVertex(const SimplexId &edgeId,
                              const int &localVertexId,
                              SimplexId &vertexId) const = 0;

    virtual SimplexId getNumberOfCells() const = 0;

    virtual SimplexId getNumberOfEdges() const = 0;

    virtual SimplexId getNumberOfTriangles() const = 0;

    virtual SimplexId getNumberOfVertices() const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getTriangles() = 0;

    virtual int getTriangleEdge(const SimplexId &triangleId,
                                const int &localEdgeId,
                                SimplexId &edgeId) const = 0;

    virtual SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getTriangleEdges() = 0;

    virtual int getTriangleLink(const SimplexId &triangleId,
                                const int &localLinkId,
                                SimplexId &linkId) const = 0;

    virtual SimplexId
      getTriangleLinkNumber(const SimplexId &triangleId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getTriangleLinks() = 0;

    virtual int getTriangleStar(const SimplexId &triangleId,
                                const int &localStarId,
                                SimplexId &starId) const = 0;

    virtual SimplexId
      getTriangleStarNumber(const SimplexId &triangleId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getTriangleStars() = 0;

    virtual int getTriangleVertex(const SimplexId &triangleId,
                                  const int &localVertexId,
                                  SimplexId &vertexId) const = 0;

    virtual int getVertexEdge(const SimplexId &vertexId,
                              const int &localEdgeId,
                              SimplexId &edgeId) const = 0;

    virtual SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getVertexEdges() = 0;

    virtual int getVertexLink(const SimplexId &vertexId,
                              const int &localLinkId,
                              SimplexId &linkId) const = 0;

    virtual SimplexId getVertexLinkNumber(const SimplexId &vertexId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getVertexLinks() = 0;

    virtual int getVertexNeighbor(const SimplexId &vertexId,
                                  const int &localNeighborId,
                                  SimplexId &neighborId) const = 0;

    virtual SimplexId
      getVertexNeighborNumber(const SimplexId &vertexId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getVertexNeighbors() = 0;

    virtual int getVertexPoint(const SimplexId &vertexId,
                               float &x,
                               float &y,
                               float &z) const = 0;

    virtual int getVertexStar(const SimplexId &vertexId,
                              const int &localStarId,
                              SimplexId &starId) const = 0;

    virtual SimplexId getVertexStarNumber(const SimplexId &vertexId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getVertexStars() = 0;

    virtual int getVertexTriangle(const SimplexId &vertexId,
                                  const int &localTriangleId,
                                  SimplexId &triangleId) const = 0;

    virtual SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const = 0;

    virtual const std::vector<std::vector<SimplexId>> *getVertexTriangles() = 0;

    virtual inline bool hasPreprocessedBoundaryEdges() const {
      return hasPreprocessedBoundaryEdges_;
    }

    virtual inline bool hasPreprocessedBoundaryTriangles() const {
      return hasPreprocessedBoundaryTriangles_;
    }

    virtual inline bool hasPreprocessedBoundaryVertices() const {
      return hasPreprocessedBoundaryVertices_;
    }

    virtual inline bool hasPreprocessedCellEdges() const {
      return hasPreprocessedCellEdges_;
    }

    virtual inline bool hasPreprocessedCellNeighbors() const {
      return hasPreprocessedCellNeighbors_;
    }

    virtual inline bool hasPreprocessedCellTriangles() const {
      return hasPreprocessedCellTriangles_;
    }

    virtual inline bool hasPreprocessedEdgeLinks() const {
      return hasPreprocessedEdgeLinks_;
    }

    virtual inline bool hasPreprocessedEdgeStars() const {
      return hasPreprocessedEdgeStars_;
    }

    virtual inline bool hasPreprocessedEdgeTriangles() const {
      return hasPreprocessedEdgeTriangles_;
    }

    virtual inline bool hasPreprocessedEdges() const {
      return hasPreprocessedEdges_;
    }

    virtual inline bool hasPreprocessedTriangles() const {
      return hasPreprocessedTriangles_;
    }

    virtual inline bool hasPreprocessedTriangleEdges() const {
      return hasPreprocessedTriangleEdges_;
    }

    virtual inline bool hasPreprocessedTriangleLinks() const {
      return hasPreprocessedTriangleLinks_;
    }

    virtual inline bool hasPreprocessedTriangleStars() const {
      return hasPreprocessedTriangleStars_;
    }

    virtual inline bool hasPreprocessedVertexEdges() const {
      return hasPreprocessedVertexEdges_;
    }

    virtual inline bool hasPreprocessedVertexLinks() const {
      return hasPreprocessedVertexLinks_;
    }

    virtual inline bool hasPreprocessedVertexNeighbors() const {
      return hasPreprocessedVertexNeighbors_;
    }

    virtual inline bool hasPreprocessedVertexStars() const {
      return hasPreprocessedVertexStars_;
    }

    virtual inline bool hasPreprocessedVertexTriangles() const {
      return hasPreprocessedVertexTriangles_;
    }

    virtual bool isEdgeOnBoundary(const SimplexId &edgeId) const = 0;

    virtual bool isEmpty() const = 0;

    virtual bool isTriangleOnBoundary(const SimplexId &triangleId) const = 0;

    virtual bool isVertexOnBoundary(const SimplexId &vertexId) const = 0;

    virtual int preprocessBoundaryEdges() {
      preprocessEdges();
      hasPreprocessedBoundaryEdges_ = true;
      return 0;
    }

    virtual int preprocessBoundaryTriangles() {
      preprocessTriangles();
      hasPreprocessedBoundaryTriangles_ = true;
      return 0;
    }

    virtual int preprocessBoundaryVertices() {
      hasPreprocessedBoundaryVertices_ = true;
      return 0;
    }

    virtual int preprocessCellEdges() {
      preprocessEdges();
      hasPreprocessedCellEdges_ = true;
      return 0;
    }

    virtual int preprocessCellNeighbors() {
      hasPreprocessedCellNeighbors_ = true;
      return 0;
    }

    virtual int preprocessCellTriangles() {
      preprocessTriangles();
      hasPreprocessedCellTriangles_ = true;
      return 0;
    }

    virtual int preprocessEdges() {
      hasPreprocessedEdges_ = true;
      return 0;
    }

    virtual int preprocessEdgeLinks() {
      preprocessEdges();
      hasPreprocessedEdgeLinks_ = true;
      return 0;
    }

    virtual int preprocessEdgeStars() {
      preprocessEdges();
      hasPreprocessedEdgeStars_ = true;
      return 0;
    }

    virtual int preprocessEdgeTriangles() {
      preprocessEdges();
      preprocessTriangles();
      hasPreprocessedEdgeTriangles_ = true;
      return 0;
    }

    virtual int preprocessTriangles() {
      hasPreprocessedTriangles_ = true;
      return 0;
    }

    virtual int preprocessTriangleEdges() {
      preprocessEdges();
      preprocessTriangles();
      hasPreprocessedTriangleEdges_ = true;
      return 0;
    }

    virtual int preprocessTriangleLinks() {
      preprocessTriangles();
      hasPreprocessedTriangleLinks_ = true;
      return 0;
    }

    virtual int preprocessTriangleStars() {
      preprocessTriangles();
      hasPreprocessedTriangleStars_ = true;
      return 0;
    }

    virtual int preprocessVertexEdges() {
      preprocessEdges();
      hasPreprocessedVertexEdges_ = true;
      return 0;
    }

    virtual int preprocessVertexLinks() {
      hasPreprocessedVertexLinks_ = true;
      return 0;
    }

    virtual int preprocessVertexNeighbors() {
      hasPreprocessedVertexNeighbors_ = true;
      return 0;
    }

    virtual int preprocessVertexStars() {
      hasPreprocessedVertexStars_ = true;
      return 0;
    }

    virtual int preprocessVertexTriangles() {
      preprocessTriangles();
      hasPreprocessedVertexTriangles_ = true;
      return 0;
    }

  protected:
    // empty wrapping to VTK for now
    bool needsToAbort() {
      return false;
    };

    template <class itemType>
    size_t tableFootprint(const std::vector<itemType> &table,
                          const std::string tableName = "",
                          std::stringstream *msg = NULL) const {

      if((table.size()) && (tableName.length()) && (msg)) {
        (*msg) << "[AbstractTriangulation] " << tableName << ": "
               << table.size() * sizeof(itemType) << " bytes" << std::endl;
      }

      return table.size() * sizeof(itemType);
    }

    template <class itemType>
    size_t tableTableFootprint(const std::vector<std::vector<itemType>> &table,
                               const std::string tableName = "",
                               std::stringstream *msg = NULL) const;

    int updateProgress(const float &progress) {
      return 0;
    };

    bool hasPreprocessedBoundaryEdges_, hasPreprocessedBoundaryTriangles_,
      hasPreprocessedBoundaryVertices_, hasPreprocessedCellEdges_,
      hasPreprocessedCellNeighbors_, hasPreprocessedCellTriangles_,
      hasPreprocessedEdges_, hasPreprocessedEdgeLinks_,
      hasPreprocessedEdgeStars_, hasPreprocessedEdgeTriangles_,
      hasPreprocessedTriangles_, hasPreprocessedTriangleEdges_,
      hasPreprocessedTriangleLinks_, hasPreprocessedTriangleStars_,
      hasPreprocessedVertexEdges_, hasPreprocessedVertexLinks_,
      hasPreprocessedVertexNeighbors_, hasPreprocessedVertexStars_,
      hasPreprocessedVertexTriangles_;

    std::vector<bool> boundaryEdges_, boundaryTriangles_, boundaryVertices_;
    std::vector<std::vector<SimplexId>> cellEdgeList_;
    std::vector<std::vector<SimplexId>> cellNeighborList_;
    std::vector<std::vector<SimplexId>> cellTriangleList_;
    std::vector<std::vector<SimplexId>> edgeLinkList_;
    std::vector<std::pair<SimplexId, SimplexId>> edgeList_;
    std::vector<std::vector<SimplexId>> edgeStarList_;
    std::vector<std::vector<SimplexId>> edgeTriangleList_;
    std::vector<std::vector<SimplexId>> triangleList_;
    std::vector<std::vector<SimplexId>> triangleEdgeList_;
    std::vector<std::vector<SimplexId>> triangleLinkList_;
    std::vector<std::vector<SimplexId>> triangleStarList_;
    std::vector<std::vector<SimplexId>> vertexEdgeList_;
    std::vector<std::vector<SimplexId>> vertexLinkList_;
    std::vector<std::vector<SimplexId>> vertexNeighborList_;
    std::vector<std::vector<SimplexId>> vertexStarList_;
    std::vector<std::vector<SimplexId>> vertexTriangleList_;
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include                  <AbstractTriangulation.cpp>

#endif // _ABSTRACTTRIANGULATION_H
