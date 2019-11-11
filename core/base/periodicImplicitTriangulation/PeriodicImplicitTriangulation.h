/// \ingroup base
/// \class ttk::PeriodicImplicitTriangulation
/// \author Talha Bin Masood <talha.bin.masood@liu.se>
/// \date September 2019
///
/// \brief TTK triangulation class for grids with periodic boundary conditions
/// implemented in all directions.
///
/// \sa ttk::Triangulation
/// \sa ttk::Triangulation::setPeriodicBoundaryConditions
///

#ifndef _PERIODICIMPLICITTRIANGULATION_H
#define _PERIODICIMPLICITTRIANGULATION_H

// base code includes
#include <AbstractTriangulation.h>

#ifdef _WIN32
#include <ciso646>
#endif

#include <array>

namespace ttk {

  class PeriodicImplicitTriangulation final : public AbstractTriangulation {

  public:
    PeriodicImplicitTriangulation();
    ~PeriodicImplicitTriangulation();

    int getCellEdge(const SimplexId &cellId,
                    const int &id,
                    SimplexId &edgeId) const override;

    SimplexId getCellEdgeNumber(const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *getCellEdges() override;

    int getCellNeighbor(const SimplexId &cellId,
                        const int &localNeighborId,
                        SimplexId &neighborId) const override;

    SimplexId getCellNeighborNumber(const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *getCellNeighbors() override;

    int getCellTriangle(const SimplexId &cellId,
                        const int &id,
                        SimplexId &triangleId) const override;

    SimplexId getCellTriangleNumber(const SimplexId &cellId) const override {
      // NOTE: the output is always 4 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 4;
    };

    const std::vector<std::vector<SimplexId>> *getCellTriangles() override;

    int getCellVertex(const SimplexId &cellId,
                      const int &localVertexId,
                      SimplexId &vertexId) const override;

    SimplexId getCellVertexNumber(const SimplexId &cellId) const override;

    int getDimensionality() const override {
      return dimensionality_;
    };

    int getEdgeLink(const SimplexId &edgeId,
                    const int &localLinkId,
                    SimplexId &linkId) const override;

    SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const override;

    const std::vector<std::vector<SimplexId>> *getEdgeLinks() override;

    int getEdgeStar(const SimplexId &edgeId,
                    const int &localStarId,
                    SimplexId &starId) const override;

    SimplexId getEdgeStarNumber(const SimplexId &edgeId) const override;

    const std::vector<std::vector<SimplexId>> *getEdgeStars() override;

    int getEdgeTriangle(const SimplexId &edgeId,
                        const int &id,
                        SimplexId &triangleId) const override;

    SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const override;

    const std::vector<std::vector<SimplexId>> *getEdgeTriangles() override;

    int getEdgeVertex(const SimplexId &edgeId,
                      const int &localVertexId,
                      SimplexId &vertexId) const override;

    const std::vector<std::pair<SimplexId, SimplexId>> *getEdges() override;

    SimplexId getNumberOfCells() const override {
      return cellNumber_;
    };

    SimplexId getNumberOfEdges() const override {
      return edgeNumber_;
    };

    SimplexId getNumberOfTriangles() const override {
      return triangleNumber_;
    };

    SimplexId getNumberOfVertices() const override {
      return vertexNumber_;
    };

    int getTetrahedronEdge(const SimplexId &tetId,
                           const int &id,
                           SimplexId &edgeId) const;

    int getTetrahedronEdges(std::vector<std::vector<SimplexId>> &edges) const;

    int getTetrahedronTriangle(const SimplexId &tetId,
                               const int &id,
                               SimplexId &triangleId) const;

    int getTetrahedronTriangles(
      std::vector<std::vector<SimplexId>> &triangles) const;

    int getTetrahedronNeighbor(const SimplexId &tetId,
                               const int &localNeighborId,
                               SimplexId &neighborId) const;

    SimplexId getTetrahedronNeighborNumber(const SimplexId &tetId) const;

    int getTetrahedronNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    int getTetrahedronVertex(const SimplexId &tetId,
                             const int &localVertexId,
                             SimplexId &vertexId) const;

    int getTriangleEdge(const SimplexId &triangleId,
                        const int &id,
                        SimplexId &edgeId) const override;

    SimplexId
      getTriangleEdgeNumber(const SimplexId &triangleId) const override {
      // NOTE: the output is always 3 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 3;
    }

    const std::vector<std::vector<SimplexId>> *getTriangleEdges() override;

    int getTriangleEdges(std::vector<std::vector<SimplexId>> &edges) const;

    int getTriangleLink(const SimplexId &triangleId,
                        const int &localLinkId,
                        SimplexId &linkId) const override;

    SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const override;

    const std::vector<std::vector<SimplexId>> *getTriangleLinks() override;

    int getTriangleNeighbor(const SimplexId &triangleId,
                            const int &localNeighborId,
                            SimplexId &neighborId) const;

    SimplexId getTriangleNeighborNumber(const SimplexId &triangleId) const;

    int getTriangleNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    int getTriangleStar(const SimplexId &triangleId,
                        const int &localStarId,
                        SimplexId &starId) const override;

    SimplexId getTriangleStarNumber(const SimplexId &triangleId) const override;

    const std::vector<std::vector<SimplexId>> *getTriangleStars() override;

    int getTriangleVertex(const SimplexId &triangleId,
                          const int &localVertexId,
                          SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getTriangles() override;

    int getVertexEdge(const SimplexId &vertexId,
                      const int &id,
                      SimplexId &edgeId) const override;

    SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getVertexEdges() override;

    int getVertexLink(const SimplexId &vertexId,
                      const int &localLinkId,
                      SimplexId &linkId) const override;

    SimplexId getVertexLinkNumber(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getVertexLinks() override;

    int getVertexNeighbor(const SimplexId &vertexId,
                          const int &localNeighborId,
                          SimplexId &neighborId) const override;

    SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getVertexNeighbors() override;

    int getVertexPoint(const SimplexId &vertexId,
                       float &x,
                       float &y,
                       float &z) const override;

    int getVertexStar(const SimplexId &vertexId,
                      const int &localStarId,
                      SimplexId &starId) const override;

    SimplexId getVertexStarNumber(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getVertexStars() override;

    int getVertexTriangle(const SimplexId &vertexId,
                          const int &id,
                          SimplexId &triangleId) const override;

    SimplexId getVertexTriangleNumber(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *getVertexTriangles() override;

    bool isEdgeOnBoundary(const SimplexId &edgeId) const override;

    bool isEmpty() const override {
      return !vertexNumber_;
    };

    bool isTriangleOnBoundary(const SimplexId &triangleId) const override;

    bool isVertexOnBoundary(const SimplexId &vertexId) const override;

    int setInputGrid(const float &xOrigin,
                     const float &yOrigin,
                     const float &zOrigin,
                     const float &xSpacing,
                     const float &ySpacing,
                     const float &zSpacing,
                     const int &xDim,
                     const int &yDim,
                     const int &zDim);

    std::array<SimplexId, 3>
      vertexToPositionNd(const SimplexId vertexId) const {
      std::array<SimplexId, 3> p{};
      if(dimensionality_ == 1) {
        p[0] = vertexId;
      } else if(dimensionality_ == 2) {
        vertexToPosition2d(vertexId, p.data());
      } else if(dimensionality_ == 3) {
        vertexToPosition(vertexId, p.data());
      }
      return p;
    }

    /**
     * Compute the barycenter of the points of the given edge identifier.
     */
    virtual int getEdgeIncenter(SimplexId edgeId,
                                float incenter[3]) const override {
      std::array<SimplexId, 2> vertexId;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getEdgeVertex(edgeId, i, vertexId[i]);
      }

      std::array<std::array<float, 3>, vertexId.size()> p;
      std::array<std::array<SimplexId, 3>, vertexId.size()> ind;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getVertexPoint(vertexId[i], p[i][0], p[i][1], p[i][2]);
        ind[i] = vertexToPositionNd(vertexId[i]);
      }

      for(int i = 0; i < dimensionality_; ++i) {
        if(ind[1][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[0][i] == nbvoxels_[i]) {
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
        }
      }

      for(int i = 0; i < 3; ++i) {
        incenter[i] = 0.5f * (p[0][i] + p[1][i]);
      }

      return 0;
    }

    /**
     * Compute the incenter of the points of the given triangle
     * identifier.
     */
    virtual int getTriangleIncenter(SimplexId triangleId,
                                    float incenter[3]) const override {

      std::array<SimplexId, 3> vertexId;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getTriangleVertex(triangleId, i, vertexId[i]);
      }

      std::array<std::array<float, 3>, vertexId.size()> p;
      std::array<std::array<SimplexId, 3>, vertexId.size()> ind;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getVertexPoint(vertexId[i], p[i][0], p[i][1], p[i][2]);
        ind[i] = vertexToPositionNd(vertexId[i]);
      }

      for(int i = 0; i < dimensionality_; ++i) {
        if(ind[0][i] == nbvoxels_[i]) {
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
          p[2][i] += (ind[2][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[1][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
          p[2][i] += (ind[2][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[2][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
        }
      }

      std::array<float, p.size()> d;
      for(int i = 0; i < (int)d.size(); ++i) {
        d[i] = Geometry::distance(p[(i + 1) % 3].data(), p[(i + 2) % 3].data());
      }
      const float sum = d[0] + d[1] + d[2];
      for(int i = 0; i < (int)d.size(); ++i) {
        d[i] = d[i] / sum;
      }

      for(int i = 0; i < 3; ++i) {
        incenter[i] = d[0] * p[0][i] + d[1] * p[1][i] + d[2] * p[2][i];
      }

      return 0;
    }

    /**
     * Compute the barycenter of the incenters of the triangles of the
     * given tetra identifier.
     */
    virtual int getTetraIncenter(SimplexId tetraId,
                                 float incenter[3]) const override {

      std::array<SimplexId, 4> vertexId;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getCellVertex(tetraId, i, vertexId[i]);
      }

      std::array<std::array<float, 3>, vertexId.size()> p;
      std::array<std::array<SimplexId, 3>, vertexId.size()> ind;
      for(int i = 0; i < (int)vertexId.size(); ++i) {
        getVertexPoint(vertexId[i], p[i][0], p[i][1], p[i][2]);
        ind[i] = vertexToPositionNd(vertexId[i]);
      }

      for(int i = 0; i < dimensionality_; ++i) {
        if(ind[0][i] == nbvoxels_[i]) {
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
          p[2][i] += (ind[2][i] == 0) * dimensions_[i] * spacing_[i];
          p[3][i] += (ind[3][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[1][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
          p[2][i] += (ind[2][i] == 0) * dimensions_[i] * spacing_[i];
          p[3][i] += (ind[3][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[2][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
          p[3][i] += (ind[3][i] == 0) * dimensions_[i] * spacing_[i];
        } else if(ind[3][i] == nbvoxels_[i]) {
          p[0][i] += (ind[0][i] == 0) * dimensions_[i] * spacing_[i];
          p[1][i] += (ind[1][i] == 0) * dimensions_[i] * spacing_[i];
          p[2][i] += (ind[2][i] == 0) * dimensions_[i] * spacing_[i];
        }
      }

      for(int i = 0; i < 3; ++i) {
        incenter[i] = 0.25f * (p[0][i] + p[1][i] + p[2][i] + p[3][i]);
      }
      return 0;
    }

  protected:
    int dimensionality_; //
    float origin_[3]; //
    float spacing_[3]; //
    SimplexId dimensions_[3]; // dimensions
    SimplexId nbvoxels_[3]; // nombre de voxels par axe
    SimplexId wrap_[3];

    // Vertex helper //
    SimplexId vshift_[2]; // VertexShift

    // Edge helper //
    SimplexId esetdims_[7]; // EdgeSetDimensions
    SimplexId esetshift_[7]; // EdgeSetShift
    SimplexId eshift_[14]; // EdgeShift

    // Triangle helper //
    SimplexId tsetdims_[6]; // TriangleSetDimensions
    SimplexId tsetshift_[6]; // TriangleSetShift
    SimplexId tshift_[12]; // TriangleShift

    // Tetrahedron helper //
    SimplexId tetshift_[2]; // TetrahedronShift

    SimplexId cellNumber_; // number of cells
    SimplexId vertexNumber_; // number of vertices
    SimplexId edgeNumber_; // number of edges
    SimplexId triangleNumber_; // number of triangles
    SimplexId tetrahedronNumber_; // number of tetrahedra

    // 2d helpers
    SimplexId Di_;
    SimplexId Dj_;

    // acceleration variables
    bool isAccelerated_;
    SimplexId mod_[2];
    SimplexId div_[2];

    // acceleration functions
    int checkAcceleration();
    bool isPowerOfTwo(unsigned long long int v, unsigned long long int &r);

    //\cond
    // 2D //
    void vertexToPosition2d(const SimplexId vertex, SimplexId p[2]) const;
    void
      edgeToPosition2d(const SimplexId edge, const int k, SimplexId p[2]) const;
    void triangleToPosition2d(const SimplexId triangle, SimplexId p[2]) const;

    SimplexId getVertexNeighbor2d(const SimplexId p[2],
                                  const SimplexId v,
                                  const int id) const;
    SimplexId getVertexEdge2d(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2d(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2d(const SimplexId p[2], const int id) const;

    SimplexId getEdgeTriangle2dL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle2dH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle2dD1(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLink2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dH(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dD1(const SimplexId p[2], const int id) const;

    SimplexId getEdgeStar2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeStar2dH(const SimplexId p[2], const int id) const;

    // 3D //
    void vertexToPosition(const SimplexId vertex, SimplexId p[3]) const;
    void
      edgeToPosition(const SimplexId edge, const int k, SimplexId p[3]) const;
    void triangleToPosition(const SimplexId triangle,
                            const int k,
                            SimplexId p[3]) const;
    void tetrahedronToPosition(const SimplexId tetrahedron,
                               SimplexId p[3]) const;

    SimplexId getVertexNeighbor3d(const SimplexId p[3],
                                  const SimplexId v,
                                  const int id) const;
    SimplexId getVertexEdge3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangle3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexLink3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexStar3d(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangle3dL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD3(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD4(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLinkL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD3(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD4(const SimplexId p[3], const int id) const;

    SimplexId getEdgeStarL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleVertexF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleEdgeF_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeF_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_1(const SimplexId p[3], const int id) const;

    SimplexId getTriangleLinkF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleStarF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronVertexABCG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBCDG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexABEG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBEFG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBFGH(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBDGH(const SimplexId p[3],
                                       const int id) const;

    SimplexId getTetrahedronEdgeABCG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBCDG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeABEG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBEFG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBFGH(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBDGH(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronTriangleABCG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBCDG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleABEG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBEFG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBFGH(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBDGH(const SimplexId p[3],
                                         const int id) const;

    SimplexId getTetrahedronNeighborABCG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBCDG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborABEG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBEFG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBFGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBDGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    //\endcond
  };
} // namespace ttk

inline void
  ttk::PeriodicImplicitTriangulation::vertexToPosition2d(const SimplexId vertex,
                                                         SimplexId p[2]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = vertex >> div_[0];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = vertex / vshift_[0];
  }
}

inline void ttk::PeriodicImplicitTriangulation::edgeToPosition2d(
  const SimplexId edge, const int k, SimplexId p[2]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = e / eshift_[2 * k];
}

inline void ttk::PeriodicImplicitTriangulation::triangleToPosition2d(
  const SimplexId triangle, SimplexId p[2]) const {
  p[0] = triangle % tshift_[0];
  p[1] = triangle / tshift_[0];
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2d(
  const SimplexId p[2], const SimplexId v, const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[Di_])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[Dj_])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return v - 1 + wrapXLeft;
    case 1:
      return v - vshift_[0] + wrapYTop;
    case 2:
      return v - vshift_[0] + 1 + wrapXRight + wrapYTop;
    case 3:
      return v + 1 + wrapXRight;
    case 4:
      return v + vshift_[0] + wrapYBottom;
    case 5:
      return v + vshift_[0] - 1 + wrapXLeft + wrapYBottom;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexEdge2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + wrapYTop;
    case 1:
      return p[0] + p[1] * eshift_[0] - 1 + wrapXLeft;
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] + wrapYTop;
    case 3:
      return p[0] + p[1] * eshift_[0];
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2];
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1 + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexStar2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + wrapXLeft;
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1 + wrapXLeft;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + wrapYTop;
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1 + wrapYTop;
    case 5:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1 + wrapXLeft
             + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexLink2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[Di_])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[Dj_])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1 + wrapXLeft;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1 + wrapXLeft + wrapYBottom;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4];
    case 3:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1 + wrapXRight
             + wrapYTop;
    case 4:
      return p[0] + (p[1] - 1) * eshift_[0] + wrapYTop;
    case 5:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1 + wrapXLeft
             + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle2dL(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1 + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle2dH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1 + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle2dD1(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLink2dL(const SimplexId p[2],
                                                     const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]));
    }
    return -1;
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]))
               + wrap_[1];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0] - wrap_[1];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]));
    }
    return -1;
  }
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLink2dH(const SimplexId p[2],
                                                     const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]));
    }
    return -1;
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1 + wrap_[0];
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]));
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]))
               - wrap_[0];
    }
    return -1;
  }
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLink2dD1(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapX = (p[0] < nbvoxels_[Di_]) ? 0 : wrap_[0];
  SimplexId wrapY = (p[1] < nbvoxels_[Dj_]) ? 0 : wrap_[1];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0];
    case 1:
      return p[0] + (p[1] + 1) * vshift_[0] + 1 - wrapX - wrapY;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStar2dL(const SimplexId p[2],
                                                     const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    }
    return -1;
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1 + 2 * wrap_[1];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    }
    return -1;
  }
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStar2dH(const SimplexId p[2],
                                                     const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    }
    return -1;
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1 + 2 * wrap_[0];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    }
    return -1;
  }
}

inline void
  ttk::PeriodicImplicitTriangulation::vertexToPosition(const SimplexId vertex,
                                                       SimplexId p[3]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = (vertex & mod_[1]) >> div_[0];
    p[2] = vertex >> div_[1];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = (vertex % vshift_[1]) / vshift_[0];
    p[2] = vertex / vshift_[1];
  }
}

inline void ttk::PeriodicImplicitTriangulation::edgeToPosition(
  const SimplexId edge, const int k, SimplexId p[3]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = (e % eshift_[2 * k + 1]) / eshift_[2 * k];
  p[2] = e / eshift_[2 * k + 1];
}

inline void ttk::PeriodicImplicitTriangulation::triangleToPosition(
  const SimplexId triangle, const int k, SimplexId p[3]) const {
  const SimplexId t = (k) ? triangle - tsetshift_[k - 1] : triangle;
  p[0] = t % tshift_[2 * k];
  p[1] = (t % tshift_[2 * k + 1]) / tshift_[2 * k];
  p[2] = t / tshift_[2 * k + 1];
}

inline void ttk::PeriodicImplicitTriangulation::tetrahedronToPosition(
  const SimplexId tetrahedron, SimplexId p[3]) const {
  p[0] = (tetrahedron % tetshift_[0]) / 6;
  p[1] = (tetrahedron % tetshift_[1]) / tetshift_[0];
  p[2] = tetrahedron / tetshift_[1];
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor3d(
  const SimplexId p[3], const SimplexId v, const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return v - vshift_[0] - vshift_[1] + wrapYTop + wrapZBack;
    case 1:
      return v + 1 - vshift_[0] - vshift_[1] + wrapXRight + wrapYTop
             + wrapZBack;
    case 2:
      return v - vshift_[1] + wrapZBack;
    case 3:
      return v + 1 - vshift_[1] + wrapXRight + wrapZBack;
    case 4:
      return v - vshift_[0] + wrapYTop;
    case 5:
      return v + 1 - vshift_[0] + wrapXRight + wrapYTop;
    case 6:
      return v + 1 + wrapXRight;
    case 7:
      return v - 1 + vshift_[1] + wrapXLeft + wrapZFront;
    case 8:
      return v + vshift_[1] + wrapZFront;
    case 9:
      return v - 1 + wrapXLeft;
    case 10:
      return v - 1 + vshift_[0] + wrapXLeft + wrapYBottom;
    case 11:
      return v + vshift_[0] + wrapYBottom;
    case 12:
      return v - 1 + vshift_[0] + vshift_[1] + wrapXLeft + wrapYBottom
             + wrapZFront;
    case 13:
      return v + vshift_[0] + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexEdge3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9] + wrapYTop + wrapZBack;
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapYTop + wrapZBack;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5]
             + wrapZBack;
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11] + wrapZBack;
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + p[2] * eshift_[3]
             + wrapYTop;
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6] + p[2] * eshift_[7]
             + wrapYTop;
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 7:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11] - 1
             + wrapXLeft;
    case 8:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1 + wrapXLeft;
    case 10:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7] - 1
             + wrapXLeft;
    case 11:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 12:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13] - 1
             + wrapXLeft;
    case 13:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexLink3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -2 * wrap_[2];
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + wrapXLeft;
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1 + wrapXLeft;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + wrapXLeft + wrapYBottom;
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft + wrapYBottom;
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + wrapXLeft + wrapZFront;
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
             + wrapXLeft + wrapZFront;
    case 8:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + wrapYTop;
    case 9:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapYTop;
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapXLeft + wrapYTop;
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft + wrapYTop;
    case 12:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + wrapZBack;
    case 13:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapZBack;
    case 14:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapXLeft + wrapZBack;
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapXLeft + wrapZBack;
    case 16:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + wrapXRight + wrapYTop + wrapZBack;
    case 17:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapXRight + wrapYTop + wrapZBack;
    case 18:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + wrapYTop + wrapZBack;
    case 19:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapYTop + wrapZBack;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + wrapYTop + wrapZBack;
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1] + 1
             + wrapYTop + wrapZBack;
    case 22:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapXLeft + wrapYTop + wrapZBack;
    case 23:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + wrapXLeft + wrapYTop + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexTriangle3d(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + wrapXLeft;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1 + wrapXLeft;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1 + wrapXLeft;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + wrapXLeft;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + wrapXLeft;
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + wrapXLeft;
    case 12:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + wrapYTop + wrapZBack;
    case 13:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapYTop + wrapZBack;
    case 14:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + wrapYTop + wrapZBack;
    case 15:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapYTop + wrapZBack;
    case 16:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + wrapYTop + wrapZBack;
    case 17:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapYTop + wrapZBack;
    case 18:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + wrapYTop + wrapZBack;
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapYTop + wrapZBack;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + wrapYTop;
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapYTop;
    case 22:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + wrapZBack;
    case 23:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapZBack;
    case 24:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft + wrapYTop;
    case 25:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft + wrapYTop;
    case 26:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5] + wrapYTop;
    case 27:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapZBack;
    case 28:
      return p[0] * 2 + tsetshift_[1] + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapZBack;
    case 29:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 30:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapXLeft + wrapZBack;
    case 31:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapXLeft + wrapYTop + wrapZBack;
    case 32:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 33:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 34:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapXLeft + wrapZBack;
    case 35:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getVertexStar3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + wrapXLeft;
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1
             + wrapXLeft;
    case 2:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2
             + wrapXLeft;
    case 3:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3
             + wrapXLeft;
    case 4:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4
             + wrapXLeft;
    case 5:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5
             + wrapXLeft;
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + wrapYTop + wrapZBack;
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 1 + wrapYTop + wrapZBack;
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2 + wrapYTop + wrapZBack;
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3 + wrapYTop + wrapZBack;
    case 10:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4 + wrapYTop + wrapZBack;
    case 11:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 5 + wrapYTop + wrapZBack;
    case 12:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
    case 13:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
    case 14:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + wrapYTop;
    case 15:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1] + 1
             + wrapYTop;
    case 16:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1 + wrapXLeft + wrapYTop;
    case 17:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5 + wrapXLeft + wrapYTop;
    case 18:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 2
             + wrapZBack;
    case 19:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 3
             + wrapZBack;
    case 20:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3 + wrapXLeft + wrapZBack;
    case 21:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4 + wrapXLeft + wrapZBack;
    case 22:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4 + wrapXLeft + wrapYTop + wrapZBack;
    case 23:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5 + wrapXLeft + wrapYTop + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dL(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapYTop;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 2:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapZBack;
    case 3:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapYTop + wrapZBack;
    case 4:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 5:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapZBack;
    case 4:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapZBack;
    case 5:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dP(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1 + wrapXLeft;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 4:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapYTop;
    case 5:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5] + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dD1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dD2(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dD3(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeTriangle3dD4(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 4:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 5:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkL(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + p[2] * eshift_[13] + wrapYTop;
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] - 1) * eshift_[4]
             + (p[2] - 1) * eshift_[5] + wrapXRight + wrapYTop + wrapZBack;
    case 4:
      return esetshift_[0] + p[0] + 1 + (p[1] - 1) * eshift_[2]
             + (p[2] - 1) * eshift_[3] + wrapXRight + wrapYTop + wrapZBack;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkH(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1] + wrapZBack;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXLeft + wrapYBottom;
    case 3:
      return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5] + wrapXRight + wrapZBack;
    case 4:
      return esetshift_[5] + (p[0] - 1) + p[1] * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapXLeft + wrapZBack;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkP(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1] + wrapYTop;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    case 3:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + p[2] * eshift_[13] - 1 + wrapXLeft + wrapYTop;
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             - 1 + wrapXLeft + wrapZFront;
    case 5:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkD1(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11] + wrapZBack;
    case 1:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 2:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + (p[2] - 1) * eshift_[9]
             + 1 + wrapXRight + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkD2(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 1:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] - 1 + wrapXLeft + wrapYBottom;
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             - 1 + wrapXLeft + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkD3(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6] + p[2] * eshift_[7]
             + wrapYTop;
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 2:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 3:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8] + p[2] * eshift_[9]
             + 1 + wrapXRight + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeLinkD4(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1] + wrapYBottom;
    case 1:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1] + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXRight + wrapYBottom;
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 5:
      return esetshift_[0] + p[0] + 1 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapXRight + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarL(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 4:
      return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + 4 + wrapYTop + wrapZBack;
    case 5:
      return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + 5 + wrapYTop + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarH(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2
             + wrapZBack;
    case 2:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
             + 4 + wrapXLeft + wrapZBack;
    case 4:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1
             + wrapXLeft;
    case 5:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5
             + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarP(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
             + 5 + wrapXLeft + wrapYTop;
    case 1:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + wrapYTop;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 4:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 3
             + wrapXLeft;
    case 5:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 4
             + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarD1(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 1;
    case 2:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 4
             + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarD2(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 2:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5
             + wrapXLeft;
    case 3:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 4
             + wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getEdgeStarD3(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 5
             + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexF(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
               + wrapXRight + wrapYBottom;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
               + wrapXRight + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexC(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
      case 2:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexD1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
               + wrapXRight + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexD2(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleVertexD3(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
               + wrapXRight + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeF_0(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeF_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1]
             + wrapYBottom;
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeH_0(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeH_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapZFront;
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeC_0(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapYBottom;
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeC_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapZFront;
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD1_0(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight;
    case 1:
      return esetshift_[4] + p[0] / 2 + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD1_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapZFront;
    case 1:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD2_0(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD2_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD3_0(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapYBottom;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleEdgeD3_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1 + wrapXRight;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7] + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkF(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1
             + wrapXRight + wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkH(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1
             + wrapXRight + wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkC(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == 0)
    wrapXLeft = wrap_[0];
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkD1(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1]
               + wrapYBottom;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkD2(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight + wrapYBottom;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1]
               + wrapYBottom;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1]
               + wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleLinkD3(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1]
               + wrapZFront;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight + wrapYBottom;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarF(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
               + 4 + wrapZBack;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 3
               + wrapZBack;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarH(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
      case 1:
        return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
               + 5 + wrapYTop;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1] + 1
               + wrapYTop;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarC(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] < 2)
    wrapXLeft = 6 * wrap_[0];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return (p[0] / 2 - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4 + wrapXLeft;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return (p[0] / 2 - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 5 + wrapXLeft;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarD1(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarD2(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTriangleStarD3(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexABCG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + wrapYBottom;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBCDG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + wrapYBottom;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
             + wrapXRight + wrapYBottom;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexABEG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBEFG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + wrapZFront;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
             + wrapXRight + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBFGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
             + wrapXRight + wrapZFront;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBDGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
             + wrapXRight + wrapYBottom;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeABCG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 4:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBCDG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1] + wrapYBottom;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2] + p[2] * eshift_[3]
             + wrapXRight;
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeABEG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBEFG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1] + wrapZFront;
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4] + p[2] * eshift_[5]
             + wrapXRight;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBFGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapXRight + wrapZFront;
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4] + p[2] * eshift_[5]
             + wrapXRight;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 4:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8] + p[2] * eshift_[9]
             + wrapXRight;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBDGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2] + p[2] * eshift_[3]
             + wrapXRight;
    case 2:
      return esetshift_[1] + (p[0] + 1) + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXRight + wrapYBottom;
    case 3:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8] + p[2] * eshift_[9]
             + wrapXRight;
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleABCG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBCDG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 3:
      return (p[1] < nbvoxels_[1])
               ? tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3]
               : tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] - wrap_[1] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleABEG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBEFG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return (p[2] < nbvoxels_[2])
               ? p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
               : p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
                   - wrap_[2] * 2;
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBFGH(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
    case 2:
      return (p[0] < nbvoxels_[0])
               ? tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 3
               : tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 3 - wrap_[0] * 2;
    case 3:
      return (p[2] < nbvoxels_[2])
               ? p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
               : p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
                   - wrap_[2] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBDGH(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return (p[0] < nbvoxels_[0])
               ? tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 2
               : tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 2 - wrap_[0] * 2;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 3:
      return (p[1] < nbvoxels_[1])
               ? tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] + 1
               : tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] + 1 - wrap_[1] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborABCG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t + 1;
    case 1:
      return t + 2;
    case 2:
      return p[0] > 0 ? t - 1 : t - 1 + wrap_[0] * 6;
    case 3:
      return p[2] > 0 ? t - tetshift_[1] + 3
                      : t - tetshift_[1] + 3 + wrap_[2] * 6;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBCDG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 4;
    case 2:
      return p[2] > 0 ? t - tetshift_[1] + 3
                      : t - tetshift_[1] + 3 + wrap_[2] * 6;
    case 3:
      return p[1] < nbvoxels_[1] ? t + tetshift_[0] + 1
                                 : t + tetshift_[0] + 1 - wrap_[1] * 6;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborABEG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 2;
    case 1:
      return t + 1;
    case 2:
      return p[0] > 0 ? t - 4 : t - 4 + wrap_[0] * 6;
    case 3:
      return p[1] > 0 ? t - tetshift_[0] - 1
                      : t - tetshift_[0] - 1 + wrap_[1] * 6;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBEFG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      return p[1] > 0 ? t - tetshift_[0] + 2
                      : t - tetshift_[0] + 2 + wrap_[1] * 6;
    case 3:
      return p[2] < nbvoxels_[2] ? t + tetshift_[1] - 3
                                 : t + tetshift_[1] - 3 - wrap_[2] * 6;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBFGH(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      return p[0] < nbvoxels_[0] ? t + 4 : t + 4 - wrap_[0] * 6;
    case 3:
      return p[2] < nbvoxels_[2] ? t + tetshift_[1] - 3
                                 : t + tetshift_[1] - 3 - wrap_[2] * 6;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBDGH(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t - 4;
    case 2:
      return p[0] < nbvoxels_[0] ? t + 1 : t + 1 - wrap_[0] * 6;
    case 3:
      return p[1] < nbvoxels_[1] ? t + tetshift_[0] - 2
                                 : t + tetshift_[0] - 2 - wrap_[1] * 6;
  }
  return -1;
}

#endif // _PERIODICIMPLICITTRIANGULATION_H
