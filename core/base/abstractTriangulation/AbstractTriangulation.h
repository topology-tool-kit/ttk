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
#include <Geometry.h>
#include <Wrapper.h>

namespace ttk {

  class AbstractTriangulation : public Wrapper {

  public:
    AbstractTriangulation();

    virtual ~AbstractTriangulation();

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

    virtual inline bool hasPreconditionedBoundaryEdges() const {
      return hasPreconditionedBoundaryEdges_;
    }

    virtual inline bool hasPreconditionedBoundaryTriangles() const {
      return hasPreconditionedBoundaryTriangles_;
    }

    virtual inline bool hasPreconditionedBoundaryVertices() const {
      return hasPreconditionedBoundaryVertices_;
    }

    virtual inline bool hasPreconditionedCellEdges() const {
      return hasPreconditionedCellEdges_;
    }

    virtual inline bool hasPreconditionedCellNeighbors() const {
      return hasPreconditionedCellNeighbors_;
    }

    virtual inline bool hasPreconditionedCellTriangles() const {
      return hasPreconditionedCellTriangles_;
    }

    virtual inline bool hasPreconditionedEdgeLinks() const {
      return hasPreconditionedEdgeLinks_;
    }

    virtual inline bool hasPreconditionedEdgeStars() const {
      return hasPreconditionedEdgeStars_;
    }

    virtual inline bool hasPreconditionedEdgeTriangles() const {
      return hasPreconditionedEdgeTriangles_;
    }

    virtual inline bool hasPreconditionedEdges() const {
      return hasPreconditionedEdges_;
    }

    virtual inline bool hasPreconditionedTriangles() const {
      return hasPreconditionedTriangles_;
    }

    virtual inline bool hasPreconditionedTriangleEdges() const {
      return hasPreconditionedTriangleEdges_;
    }

    virtual inline bool hasPreconditionedTriangleLinks() const {
      return hasPreconditionedTriangleLinks_;
    }

    virtual inline bool hasPreconditionedTriangleStars() const {
      return hasPreconditionedTriangleStars_;
    }

    virtual inline bool hasPreconditionedVertexEdges() const {
      return hasPreconditionedVertexEdges_;
    }

    virtual inline bool hasPreconditionedVertexLinks() const {
      return hasPreconditionedVertexLinks_;
    }

    virtual inline bool hasPreconditionedVertexNeighbors() const {
      return hasPreconditionedVertexNeighbors_;
    }

    virtual inline bool hasPreconditionedVertexStars() const {
      return hasPreconditionedVertexStars_;
    }

    virtual inline bool hasPreconditionedVertexTriangles() const {
      return hasPreconditionedVertexTriangles_;
    }

    virtual bool isEdgeOnBoundary(const SimplexId &edgeId) const = 0;

    virtual bool isEmpty() const = 0;

    virtual bool isTriangleOnBoundary(const SimplexId &triangleId) const = 0;

    virtual bool isVertexOnBoundary(const SimplexId &vertexId) const = 0;

    virtual int preconditionBoundaryEdges() {
      preconditionEdges();
      hasPreconditionedBoundaryEdges_ = true;
      return 0;
    }

    virtual int preconditionBoundaryTriangles() {
      preconditionTriangles();
      hasPreconditionedBoundaryTriangles_ = true;
      return 0;
    }

    virtual int preconditionBoundaryVertices() {
      hasPreconditionedBoundaryVertices_ = true;
      return 0;
    }

    virtual int preconditionCellEdges() {
      preconditionEdges();
      hasPreconditionedCellEdges_ = true;
      return 0;
    }

    virtual int preconditionCellNeighbors() {
      hasPreconditionedCellNeighbors_ = true;
      return 0;
    }

    virtual int preconditionCellTriangles() {
      preconditionTriangles();
      hasPreconditionedCellTriangles_ = true;
      return 0;
    }

    virtual int preconditionEdges() {
      hasPreconditionedEdges_ = true;
      return 0;
    }

    virtual int preconditionEdgeLinks() {
      preconditionEdges();
      hasPreconditionedEdgeLinks_ = true;
      return 0;
    }

    virtual int preconditionEdgeStars() {
      preconditionEdges();
      hasPreconditionedEdgeStars_ = true;
      return 0;
    }

    virtual int preconditionEdgeTriangles() {
      preconditionEdges();
      preconditionTriangles();
      hasPreconditionedEdgeTriangles_ = true;
      return 0;
    }

    virtual int preconditionTriangles() {
      hasPreconditionedTriangles_ = true;
      return 0;
    }

    virtual int preconditionTriangleEdges() {
      preconditionEdges();
      preconditionTriangles();
      hasPreconditionedTriangleEdges_ = true;
      return 0;
    }

    virtual int preconditionTriangleLinks() {
      preconditionTriangles();
      hasPreconditionedTriangleLinks_ = true;
      return 0;
    }

    virtual int preconditionTriangleStars() {
      preconditionTriangles();
      hasPreconditionedTriangleStars_ = true;
      return 0;
    }

    virtual int preconditionVertexEdges() {
      preconditionEdges();
      hasPreconditionedVertexEdges_ = true;
      return 0;
    }

    virtual int preconditionVertexLinks() {
      hasPreconditionedVertexLinks_ = true;
      return 0;
    }

    virtual int preconditionVertexNeighbors() {
      hasPreconditionedVertexNeighbors_ = true;
      return 0;
    }

    virtual int preconditionVertexStars() {
      hasPreconditionedVertexStars_ = true;
      return 0;
    }

    virtual int preconditionVertexTriangles() {
      preconditionTriangles();
      hasPreconditionedVertexTriangles_ = true;
      return 0;
    }

    /**
     * Compute the barycenter of the points of the given edge identifier.
     */
    virtual int getEdgeIncenter(SimplexId edgeId, float incenter[3]) const {
      SimplexId vertexId[2];
      getEdgeVertex(edgeId, 0, vertexId[0]);
      getEdgeVertex(edgeId, 1, vertexId[1]);

      float p[6];
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
    virtual int getTriangleIncenter(SimplexId triangleId,
                                    float incenter[3]) const {

      SimplexId vertexId[3];
      if(getDimensionality() == 2) {
        getCellVertex(triangleId, 0, vertexId[0]);
        getCellVertex(triangleId, 1, vertexId[1]);
        getCellVertex(triangleId, 2, vertexId[2]);
      } else if(getDimensionality() == 3) {
        getTriangleVertex(triangleId, 0, vertexId[0]);
        getTriangleVertex(triangleId, 1, vertexId[1]);
        getTriangleVertex(triangleId, 2, vertexId[2]);
      }

      float p[9];
      getVertexPoint(vertexId[0], p[0], p[1], p[2]);
      getVertexPoint(vertexId[1], p[3], p[4], p[5]);
      getVertexPoint(vertexId[2], p[6], p[7], p[8]);

      float d[3];
      d[0] = Geometry::distance(p + 3, p + 6);
      d[1] = Geometry::distance(p, p + 6);
      d[2] = Geometry::distance(p, p + 3);
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
    virtual int getTetraIncenter(SimplexId tetraId, float incenter[3]) const {
      incenter[0] = 0.0f;
      incenter[1] = 0.0f;
      incenter[2] = 0.0f;

      float p[3];
      for(int i = 0; i < 4; ++i) {
        SimplexId triangleId;
        getCellTriangle(tetraId, i, triangleId);
        getTriangleIncenter(triangleId, p);
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
    int getCellIncenter(SimplexId cellid, int dim, float incenter[3]) const {
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
    inline bool isCellEdgePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 1) && (!hasPreconditionedCellNeighbors()))
         || ((getDimensionality() > 1) && (!hasPreconditionedCellEdges()))) {

        printMsg("CellEdge query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isCellNeighborPreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellNeighbors()) {

        printMsg("CellNeighbor query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellNeighbors() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isCellTrianglePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedCellNeighbors()))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedCellTriangles()))) {

        printMsg("CellTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionCellTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isEdgePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((getDimensionality() != 1) && (!hasPreconditionedEdges())) {

        printMsg("Edge query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isEdgeLinkPreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdgeLinks()) {

        printMsg("EdgeLink query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeLinks() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isEdgeStarPreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdgeStars()) {

        printMsg("EdgeStar query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeStars() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isEdgeTrianglePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedEdgeStars()))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedEdgeTriangles()))) {

        printMsg("EdgeTriangle query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionEdgeTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isTrianglePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if((getDimensionality() == 3) && (!hasPreconditionedTriangles())) {

        printMsg("Triangle query without pre-process!", debug::Priority::ERROR,
                 debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangles() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isTriangleEdgePreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(((getDimensionality() == 2) && (!hasPreconditionedCellEdges()))
         || ((getDimensionality() == 3)
             && (!hasPreconditionedTriangleEdges()))) {

        printMsg("TriangleEdge query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleEdges() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isTriangleLinkPreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedTriangleLinks()) {

        printMsg("TriangleLink query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleLinks() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

    inline bool isTriangleStarPreconditioned() const {

#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedTriangleStars()) {

        printMsg("TriangleStar query without pre-process!",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);
        printMsg("Please call preconditionTriangleStars() in a pre-process.",
                 debug::Priority::ERROR, debug::LineMode::NEW, std::cerr);

        return false;
      }
#endif
      return true;
    }

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

    bool hasPreconditionedBoundaryEdges_, hasPreconditionedBoundaryTriangles_,
      hasPreconditionedBoundaryVertices_, hasPreconditionedCellEdges_,
      hasPreconditionedCellNeighbors_, hasPreconditionedCellTriangles_,
      hasPreconditionedEdges_, hasPreconditionedEdgeLinks_,
      hasPreconditionedEdgeStars_, hasPreconditionedEdgeTriangles_,
      hasPreconditionedTriangles_, hasPreconditionedTriangleEdges_,
      hasPreconditionedTriangleLinks_, hasPreconditionedTriangleStars_,
      hasPreconditionedVertexEdges_, hasPreconditionedVertexLinks_,
      hasPreconditionedVertexNeighbors_, hasPreconditionedVertexStars_,
      hasPreconditionedVertexTriangles_;

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
