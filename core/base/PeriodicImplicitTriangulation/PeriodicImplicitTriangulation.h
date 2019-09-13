#ifndef _PERIODICIMPLICITTRIANGULATION_H
#define _PERIODICIMPLICITTRIANGULATION_H

// base code includes
#include                  <AbstractTriangulation.h>

#ifdef _WIN32
#include <ciso646>
#endif

namespace ttk{

  class PeriodicImplicitTriangulation : public AbstractTriangulation{

    public:
      PeriodicImplicitTriangulation();
      ~PeriodicImplicitTriangulation();

      int getCellEdge(const SimplexId &cellId, const int &id, SimplexId &edgeId) const;

      SimplexId getCellEdgeNumber(const SimplexId &cellId) const;

      const std::vector<std::vector<SimplexId>>* getCellEdges();

      int getCellNeighbor(const SimplexId &cellId,
        const int &localNeighborId, SimplexId &neighborId) const;

      SimplexId getCellNeighborNumber(const SimplexId &cellId) const;

      const std::vector<std::vector<SimplexId>>* getCellNeighbors();

      int getCellTriangle(const SimplexId &cellId, const int &id,
        SimplexId &triangleId) const;

       SimplexId getCellTriangleNumber(const SimplexId &cellId) const{
        // NOTE: the output is always 4 here. let's keep the function in there
        // in case of further generalization to CW-complexes
        return 4;
      };

      const std::vector<std::vector<SimplexId>>* getCellTriangles();

      int getCellVertex(const SimplexId &cellId, const int &localVertexId,
        SimplexId &vertexId) const;

      SimplexId getCellVertexNumber(const SimplexId &cellId) const;

      int getDimensionality() const {return dimensionality_;};

      int getEdgeLink(const SimplexId &edgeId, const int &localLinkId,
        SimplexId &linkId) const;

       SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const;

      const std::vector<std::vector<SimplexId>>* getEdgeLinks();

      int getEdgeStar(const SimplexId &edgeId, const int &localStarId,
        SimplexId &starId) const;

      SimplexId getEdgeStarNumber(const SimplexId &edgeId) const;

      const std::vector<std::vector<SimplexId>>* getEdgeStars();

      int getEdgeTriangle(const SimplexId &edgeId, const int &id,
        SimplexId &triangleId) const;

      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const;

      const std::vector<std::vector<SimplexId>>* getEdgeTriangles();

      int getEdgeVertex(const SimplexId &edgeId, const int &localVertexId,
        SimplexId &vertexId) const;

      const std::vector<std::pair<SimplexId,SimplexId>>* getEdges();

      SimplexId getNumberOfCells() const{return cellNumber_;};

      SimplexId getNumberOfEdges() const{return edgeNumber_;};

      SimplexId getNumberOfTriangles() const{return triangleNumber_;};

      SimplexId getNumberOfVertices() const{return vertexNumber_;};

      int getTetrahedronEdge(const SimplexId &tetId, const int &id,
        SimplexId &edgeId) const;

      int getTetrahedronEdges(std::vector<std::vector<SimplexId>>& edges) const;

      int getTetrahedronTriangle(const SimplexId &tetId, const int &id,
SimplexId &triangleId) const;

      int getTetrahedronTriangles(std::vector<std::vector<SimplexId>>& triangles)
const;

      int getTetrahedronNeighbor(const SimplexId &tetId, const int &localNeighborId,
SimplexId &neighborId) const;

      SimplexId getTetrahedronNeighborNumber(const SimplexId &tetId) const;

      int getTetrahedronNeighbors(std::vector<std::vector<SimplexId>>& neighbors);

      int getTetrahedronVertex(const SimplexId &tetId, const int& localVertexId,
SimplexId &vertexId) const;

      int getTriangleEdge(const SimplexId &triangleId, const int &id,
SimplexId &edgeId) const;

      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        // NOTE: the output is always 3 here. let's keep the function in there
        // in case of further generalization to CW-complexes
        return 3;
      }

      const std::vector<std::vector<SimplexId>>* getTriangleEdges();

      int getTriangleEdges(std::vector<std::vector<SimplexId>>& edges) const;

      int getTriangleLink(const SimplexId &triangleId, const int &localLinkId,
SimplexId &linkId) const;

      SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const;

      const std::vector<std::vector<SimplexId> >* getTriangleLinks();

      int getTriangleNeighbor(const SimplexId &triangleId, const int &localNeighborId,
SimplexId &neighborId) const;

      SimplexId getTriangleNeighborNumber(const SimplexId &triangleId) const;

      int getTriangleNeighbors(std::vector<std::vector<SimplexId>>& neighbors);

      int getTriangleStar(const SimplexId &triangleId, const int &localStarId,
SimplexId &starId) const;

      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const;

      const std::vector<std::vector<SimplexId>>* getTriangleStars();

      int getTriangleVertex(const SimplexId &triangleId, const int &localVertexId,
SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId>>* getTriangles();

      int getVertexEdge(const SimplexId &vertexId, const int &id,
SimplexId &edgeId) const;

      SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId>>* getVertexEdges();

      int getVertexLink(const SimplexId &vertexId, const int& localLinkId,
SimplexId &linkId) const;

      SimplexId getVertexLinkNumber(const SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId> >* getVertexLinks();

      int getVertexNeighbor(const SimplexId &vertexId, const int &localNeighborId,
SimplexId &neighborId) const;

      SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId>>* getVertexNeighbors();

      int getVertexPoint(const SimplexId &vertexId, float &x, float &y, float &z)
const;

      int getVertexStar(const SimplexId &vertexId, const int &localStarId,
SimplexId &starId) const;

      SimplexId getVertexStarNumber(const SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId>>* getVertexStars();

      int getVertexTriangle(const SimplexId &vertexId, const int &id,
SimplexId &triangleId) const;

      SimplexId getVertexTriangleNumber(const SimplexId &vertexId) const;

      const std::vector<std::vector<SimplexId>>* getVertexTriangles();

      bool isEdgeOnBoundary(const SimplexId &edgeId) const;

      bool isEmpty() const{return !vertexNumber_;};

      bool isTriangleOnBoundary(const SimplexId &triangleId) const;

      bool isVertexOnBoundary(const SimplexId &vertexId) const;

      int setInputGrid(const float &xOrigin, const float &yOrigin, const float &zOrigin,
          const float &xSpacing, const float &ySpacing, const float &zSpacing,
          const int &xDim, const int &yDim, const int &zDim);

    protected:
      int dimensionality_;//
      float origin_[3];//
      float spacing_[3];//
      SimplexId dimensions_[3];// dimensions
      SimplexId nbvoxels_[3];// nombre de voxels par axe
      SimplexId wrap_[3];

      // Vertex helper //
      SimplexId vshift_[2];// VertexShift

      // Edge helper //
      SimplexId esetdims_[7];// EdgeSetDimensions
      SimplexId esetshift_[7];// EdgeSetShift
      SimplexId eshift_[14];// EdgeShift

      // Triangle helper //
      SimplexId tsetdims_[6];// TriangleSetDimensions
      SimplexId tsetshift_[6];// TriangleSetShift
      SimplexId tshift_[12];// TriangleShift

      // Tetrahedron helper //
      SimplexId tetshift_[2];// TetrahedronShift

      SimplexId cellNumber_;// number of cells
      SimplexId vertexNumber_;// number of vertices
      SimplexId edgeNumber_;// number of edges
      SimplexId triangleNumber_;// number of triangles
      SimplexId tetrahedronNumber_;// number of tetrahedra

      // 2d helpers
      SimplexId Di_;
      SimplexId Dj_;

      // acceleration variables
      bool isAccelerated_;
      SimplexId mod_[2];
      SimplexId div_[2];

      // acceleration functions
      int checkAcceleration();
      bool isPowerOfTwo(unsigned long long int v, unsigned long long int& r);

      //\cond
      // 2D //
      void vertexToPosition2d(const SimplexId vertex, SimplexId p[2]) const;
      void edgeToPosition2d(const SimplexId edge, const int k, SimplexId p[2]) const;
      void triangleToPosition2d(const SimplexId triangle, SimplexId p[2]) const;

      SimplexId getVertexNeighbor2dA(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dB(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dAB(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dCD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dAC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dBD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighbor2dABCD(const SimplexId v,const int id) const;

      SimplexId getVertexEdge2dA(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dB(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dC(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dD(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dAB(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dCD(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dAC(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dBD(const SimplexId p[2],const int id) const;
      SimplexId getVertexEdge2dABCD(const SimplexId p[2],const int id) const;

      SimplexId getVertexStar2dA(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dB(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dC(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dD(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dAB(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dCD(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dAC(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dBD(const SimplexId p[2], const int id) const;
      SimplexId getVertexStar2dABCD(const SimplexId p[2], const int id) const;

      SimplexId getVertexLink2dA(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dB(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dC(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dD(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dAB(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dCD(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dAC(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dBD(const SimplexId p[2],const int id) const;
      SimplexId getVertexLink2dABCD(const SimplexId p[2],const int id) const;

      SimplexId getEdgeTriangleL_x0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xN(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_0y(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_ny(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_Ny(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD1_xy(const SimplexId p[3], const int id) const;

      SimplexId getEdgeLink2dL(const SimplexId p[2], const int id) const;
      SimplexId getEdgeLink2dH(const SimplexId p[2], const int id) const;
      SimplexId getEdgeLink2dD1(const SimplexId p[2], const int id) const;

      SimplexId getEdgeStar2dL(const SimplexId p[2], const int id) const;
      SimplexId getEdgeStar2dH(const SimplexId p[2], const int id) const;

      // 3D //
      void vertexToPosition(const SimplexId vertex, SimplexId p[3]) const;
      void edgeToPosition(const SimplexId edge, const int k, SimplexId p[3]) const;
      void triangleToPosition(const SimplexId triangle, const int k, SimplexId p[3]) const;
      void tetrahedronToPosition(const SimplexId tetrahedron, SimplexId p[3]) const;

      SimplexId getVertexNeighborA(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborB(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborE(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborF(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborG(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborH(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborAB(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborCD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborEF(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborGH(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborAC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborBD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborEG(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborFH(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborAE(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborBF(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborCG(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborDH(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborABDC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborEFHG(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborAEGC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborBFHD(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborAEFB(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborGHDC(const SimplexId v,const int id) const;
      SimplexId getVertexNeighborABCDEFGH(const SimplexId v,const int id) const;

      SimplexId getVertexEdgeA(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeB(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeC(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeD(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeE(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeF(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeG(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeH(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeAB(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeCD(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeEF(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeGH(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeAC(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeBD(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeEG(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeFH(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeAE(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeBF(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeCG(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeDH(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeABDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeEFHG(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeAEGC(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeBFHD(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeAEFB(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeGHDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexEdgeABCDEFGH(const SimplexId p[3],const int id) const;

      SimplexId getVertexTriangleA(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleB(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleC(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleD(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleE(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleF(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleG(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleH(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleAB(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleCD(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleEF(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleGH(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleAC(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleBD(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleEG(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleFH(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleAE(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleBF(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleCG(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleDH(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleABDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleEFHG(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleAEGC(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleBFHD(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleAEFB(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleGHDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexTriangleABCDEFGH(const SimplexId p[3],const int id) const;

      SimplexId getVertexLinkA(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkB(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkC(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkD(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkE(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkF(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkG(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkH(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkAB(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkCD(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkEF(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkGH(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkAC(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkBD(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkEG(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkFH(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkAE(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkBF(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkCG(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkDH(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkABDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkEFHG(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkAEGC(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkBFHD(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkAEFB(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkGHDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexLinkABCDEFGH(const SimplexId p[3],const int id) const;

      SimplexId getVertexStarA(const SimplexId p[3], const int id) const;
      SimplexId getVertexStarB(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarC(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarD(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarE(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarF(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarG(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarH(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarAB(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarCD(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarEF(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarGH(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarAC(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarBD(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarEG(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarFH(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarAE(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarBF(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarCG(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarDH(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarABDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarEFHG(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarAEGC(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarBFHD(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarAEFB(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarGHDC(const SimplexId p[3],const int id) const;
      SimplexId getVertexStarABCDEFGH(const SimplexId p[3],const int id) const;

      SimplexId getEdgeTriangleL_x00(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_x0n(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_x0N(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xn0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xnn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xnN(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xN0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xNn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleL_xNN(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleH_0y0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_0yn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_0yN(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_ny0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_nyn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_nyN(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_Ny0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_Nyn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleH_NyN(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleP_00z(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_0nz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_0Nz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_n0z(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_nnz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_nNz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_N0z(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_Nnz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleP_NNz(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleD1_xy0(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD1_xyn(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD1_xyN(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleD2_0yz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD2_nyz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD2_Nyz(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleD3_x0z(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD3_xnz(const SimplexId p[3], const int id) const;
      SimplexId getEdgeTriangleD3_xNz(const SimplexId p[3], const int id) const;

      SimplexId getEdgeTriangleD4_xyz(const SimplexId p[3], const int id) const;

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

      SimplexId getTetrahedronVertexABCG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronVertexBCDG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronVertexABEG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronVertexBEFG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronVertexBFGH(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronVertexBDGH(const SimplexId p[3], const int id) const;

      SimplexId getTetrahedronEdgeABCG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronEdgeBCDG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronEdgeABEG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronEdgeBEFG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronEdgeBFGH(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronEdgeBDGH(const SimplexId p[3], const int id) const;

      SimplexId getTetrahedronTriangleABCG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronTriangleBCDG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronTriangleABEG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronTriangleBEFG(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronTriangleBFGH(const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronTriangleBDGH(const SimplexId p[3], const int id) const;

      SimplexId getTetrahedronNeighborABCG(const SimplexId t, const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronNeighborBCDG(const SimplexId t, const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronNeighborABEG(const SimplexId t, const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronNeighborBEFG(const SimplexId t, const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronNeighborBFGH(const SimplexId t, const SimplexId p[3], const int id) const;
      SimplexId getTetrahedronNeighborBDGH(const SimplexId t, const SimplexId p[3], const int id) const;
      //\endcond
  };
}

inline void ttk::PeriodicImplicitTriangulation::vertexToPosition2d(const SimplexId vertex,
SimplexId p[2]) const{
  if(isAccelerated_){
    p[0]=vertex&mod_[0];
    p[1]=vertex>>div_[0];
  }
  else{
    p[0]=vertex%vshift_[0];
    p[1]=vertex/vshift_[0];
  }
}

inline void ttk::PeriodicImplicitTriangulation::edgeToPosition2d(const SimplexId edge,
const int k, SimplexId p[2]) const{
  const int e=(k)?edge-esetshift_[k-1]:edge;
  p[0]=e%eshift_[2*k];
  p[1]=e/eshift_[2*k];
}

inline void ttk::PeriodicImplicitTriangulation::triangleToPosition2d(const SimplexId triangle,
SimplexId p[2]) const{
  p[0]=triangle%tshift_[0];
  p[1]=triangle/tshift_[0];
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dA(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1+wrap_[0];
    case 1: return v-vshift_[0]+wrap_[1];
    case 2: return v-vshift_[0]+1+wrap_[1];
    case 3: return v+1;
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dAB(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0]+wrap_[1];
    case 2: return v-vshift_[0]+1+wrap_[1];
    case 3: return v+1;
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dB(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0]+wrap_[1];
    case 2: return v-vshift_[0]+1-wrap_[0]+wrap_[1];
    case 3: return v+1-wrap_[0];
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dAC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1+wrap_[0];
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1;
    case 3: return v+1;
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dABCD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1;
    case 3: return v+1;
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dBD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1-wrap_[0];
    case 3: return v+1-wrap_[0];
    case 4: return v+vshift_[0];
    case 5: return v+vshift_[0]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1+wrap_[0];
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1;
    case 3: return v+1;
    case 4: return v+vshift_[0]-wrap_[1];
    case 5: return v+vshift_[0]-1+wrap_[0]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dCD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1;
    case 3: return v+1;
    case 4: return v+vshift_[0]-wrap_[1];
    case 5: return v+vshift_[0]-1-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighbor2dD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-1;
    case 1: return v-vshift_[0];
    case 2: return v-vshift_[0]+1-wrap_[0];
    case 3: return v+1-wrap_[0];
    case 4: return v+vshift_[0]-wrap_[1];
    case 5: return v+vshift_[0]-1-wrap_[1];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dA(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+wrap_[1];
    case 1: return p[0]+p[1]*eshift_[0]-1+wrap_[0];
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+wrap_[1];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dAB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+wrap_[1];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+wrap_[1];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+wrap_[1];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+wrap_[1];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dAC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1+wrap_[0];
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dABCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dBD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1+wrap_[0];
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdge2dD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];
    case 1: return p[0]+p[1]*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];
    case 3: return p[0]+p[1]*eshift_[0];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dA(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+2*wrap_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1+2*wrap_[0];
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0]+2*wrap_[1];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1+2*wrap_[0]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dAB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0]+2*wrap_[1];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0]+2*wrap_[1];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dAC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+2*wrap_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1+2*wrap_[0];
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dABCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dBD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+2*wrap_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1+2*wrap_[0];
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStar2dD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0];
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 5: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dA(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1+wrap_[0];
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1+wrap_[0];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1+wrap_[1];
    case 4: return p[0]+(p[1]-1)*eshift_[0]+wrap_[1];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1+wrap_[0]+wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dAB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1+wrap_[1];
    case 4: return p[0]+(p[1]-1)*eshift_[0]+wrap_[1];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1+wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dB(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1-wrap_[0]+wrap_[1];
    case 4: return p[0]+(p[1]-1)*eshift_[0]+wrap_[1];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1+wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dAC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1+wrap_[0];
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1+wrap_[0];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dABCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dBD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1-wrap_[0];
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dC(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1+wrap_[0];
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1+wrap_[0]-wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1+wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dCD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1-wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLink2dD(const SimplexId p[2], const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1-wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1-wrap_[0];
    case 4: return p[0]+(p[1]-1)*eshift_[0];
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_x0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+(p[Dj_]-1)*tshift_[0]+1+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+(p[Dj_]-1)*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+(p[Dj_]-1)*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_0y(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return (p[Di_]-1)*2+p[Dj_]*tshift_[0]+1+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_ny(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return (p[Di_]-1)*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_Ny(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return (p[Di_]-1)*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD1_xy(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLink2dL(const SimplexId p[2], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[Dj_]){
    switch(id){
      case 0: return p[0]+(p[1]+1)*vshift_[0];
      case 1: return ((p[0] < nbvoxels_[Di_])? (p[0]+(p[1]-1)*vshift_[0]+1) : (p[0]+(p[1]-1)*vshift_[0]+1-wrap_[0]));
    }
    return -1;
  }
  else if(p[1]==0){
    switch(id){
      case 0: return p[0]+(p[1]+1)*vshift_[0];
      case 1: return ((p[0] < nbvoxels_[Di_])? (p[0]+(p[1]-1)*vshift_[0]+1) : (p[0]+(p[1]-1)*vshift_[0]+1-wrap_[0]))+wrap_[1];
    }
    return -1;
  }
  else{
    switch(id){
      case 0: return p[0]+(p[1]+1)*vshift_[0]-wrap_[1];
      case 1: return ((p[0] < nbvoxels_[Di_])? (p[0]+(p[1]-1)*vshift_[0]+1) : (p[0]+(p[1]-1)*vshift_[0]+1-wrap_[0]));
    }
    return -1;
  }
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLink2dH(const SimplexId p[2], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[Di_]){
    switch(id){
      case 0: return p[0]+p[1]*vshift_[0]+1;
      case 1: return ((p[1] < nbvoxels_[Dj_])? (p[0]+(p[1]+1)*vshift_[0]-1) : (p[0]+(p[1]+1)*vshift_[0]-1-wrap_[1]));
    }
    return -1;
  }
  else if(p[0]==0){
    switch(id){
      case 0: return p[0]+p[1]*vshift_[0]+1+wrap_[0];
      case 1: return ((p[1] < nbvoxels_[Dj_])? (p[0]+(p[1]+1)*vshift_[0]-1) : (p[0]+(p[1]+1)*vshift_[0]-1-wrap_[1]));
    }
    return -1;
  }
  else{
    switch(id){
      case 0: return p[0]+p[1]*vshift_[0]+1;
      case 1: return ((p[1] < nbvoxels_[Dj_])? (p[0]+(p[1]+1)*vshift_[0]-1) : (p[0]+(p[1]+1)*vshift_[0]-1-wrap_[1]))-wrap_[0];
    }
    return -1;
  }
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLink2dD1(const SimplexId p[2], const int id) const{
  SimplexId wrapX = (p[0] < nbvoxels_[Di_]) ? 0 : wrap_[0];
  SimplexId wrapY = (p[1] < nbvoxels_[Dj_]) ? 0 : wrap_[1];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0];
    case 1: return p[0]+(p[1]+1)*vshift_[0]+1-wrapX-wrapY;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStar2dL(const SimplexId p[2], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[Dj_]){
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    }
    return -1;
  }
  else if(p[1]==0){
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1+2*wrap_[1];
    }
    return -1;
  }
  else{
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    }
    return -1;
  }
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStar2dH(const SimplexId p[2], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[Di_]){
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    }
    return -1;
  }
  else if(p[0]==0){
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1+2*wrap_[0];
    }
    return -1;
  }
  else{
    switch(id){
      case 0: return p[0]*2+p[1]*tshift_[0];
      case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    }
    return -1;
  }
}

inline void ttk::PeriodicImplicitTriangulation::vertexToPosition(const SimplexId vertex,
SimplexId p[3]) const{
  if(isAccelerated_){
    p[0]=vertex&mod_[0];
    p[1]=(vertex&mod_[1])>>div_[0];
    p[2]=vertex>>div_[1];
  }
  else{
    p[0]=vertex%vshift_[0];
    p[1]=(vertex%vshift_[1])/vshift_[0];
    p[2]=vertex/vshift_[1];
  }
}

inline void ttk::PeriodicImplicitTriangulation::edgeToPosition(const SimplexId edge,
const int k, SimplexId p[3]) const{
  const int e=(k)?edge-esetshift_[k-1]:edge;
  p[0]=e%eshift_[2*k];
  p[1]=(e%eshift_[2*k+1])/eshift_[2*k];
  p[2]=e/eshift_[2*k+1];
}

inline void ttk::PeriodicImplicitTriangulation::triangleToPosition(const SimplexId triangle,
const int k, SimplexId p[3]) const{
  const SimplexId t=(k)?triangle-tsetshift_[k-1]:triangle;
  p[0]=t%tshift_[2*k];
  p[1]=(t%tshift_[2*k+1])/tshift_[2*k];
  p[2]=t/tshift_[2*k+1];
}

inline void ttk::PeriodicImplicitTriangulation::tetrahedronToPosition(
const SimplexId tetrahedron, SimplexId p[3]) const{
  p[0]=(tetrahedron%tetshift_[0])/6;
  p[1]=(tetrahedron%tetshift_[1])/tetshift_[0];
  p[2]=tetrahedron/tetshift_[1];
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborA(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborAB(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborB(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0]+wrap_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]-wrap_[0]+wrap_[2];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]-wrap_[0]+wrap_[1];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborAC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborABDC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborBD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]-wrap_[0]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborCD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[2];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0]+wrap_[2];
    case 2: return v-vshift_[1]+wrap_[2];
    case 3: return v+1-vshift_[1]-wrap_[0]+wrap_[2];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborAE(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborAEFB(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborBF(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]-wrap_[0]+wrap_[1];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborAEGC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborABCDEFGH(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborBFHD(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1];
    case 13: return v+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborCG(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0];
    case 8: return v+vshift_[1];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborGHDC(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborDH(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1];
    case 8: return v+vshift_[1];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborE(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborEF(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]+wrap_[1];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborF(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1]+wrap_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0]+wrap_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0]+wrap_[1];
    case 5: return v+1-vshift_[0]-wrap_[0]+wrap_[1];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborEG(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborEFHG(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborFH(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0];
    case 11: return v+vshift_[0];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborG(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]+wrap_[0]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1+wrap_[0];
    case 10: return v-1+vshift_[0]+wrap_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]+wrap_[0]-wrap_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborGH(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0];
    case 6: return v+1;
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1]-wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexNeighborH(const SimplexId v, const int id) const{
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];
    case 1: return v+1-vshift_[0]-vshift_[1]-wrap_[0];
    case 2: return v-vshift_[1];
    case 3: return v+1-vshift_[1]-wrap_[0];
    case 4: return v-vshift_[0];
    case 5: return v+1-vshift_[0]-wrap_[0];
    case 6: return v+1-wrap_[0];
    case 7: return v-1+vshift_[1]-wrap_[2];
    case 8: return v+vshift_[1]-wrap_[2];
    case 9: return v-1;
    case 10: return v-1+vshift_[0]-wrap_[1];
    case 11: return v+vshift_[0]-wrap_[1];
    case 12: return v-1+vshift_[0]+vshift_[1]-wrap_[1]-wrap_[2];
    case 13: return v+vshift_[0]+vshift_[1]-wrap_[1]-wrap_[2];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeA(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeAB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeAC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeABDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeBD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeCD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[2];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[2];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrap_[2];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrap_[2];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeAE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeAEFB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeBF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeAEGC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeABCDEFGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeBFHD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeCG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeGHDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeDH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeEF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9]+wrap_[1];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13]+wrap_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+wrap_[1];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrap_[1];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeEG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeEFHG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeFH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1+wrap_[0];
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1+wrap_[0];
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1+wrap_[0];
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1+wrap_[0];
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexEdgeH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkA(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0]+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0]+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1]+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1]+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[1]+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0]+2*wrap_[1]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkAB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1]+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1]+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1]+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0]+2*wrap_[1]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0]+2*wrap_[1]+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1]+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1]+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkAC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0]+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkABDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkBD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0]+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0]-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0]+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkCD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0]+2*wrap_[2];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0]+2*wrap_[2];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[2];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[2];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkAE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0]+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkAEFB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkBF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0]+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkAEGC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkABCDEFGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkBFHD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkCG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0]-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkGHDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkDH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0]-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0]+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkEF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+2*wrap_[1];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[1];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0]+2*wrap_[1];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0]+2*wrap_[1];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[1];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1+2*wrap_[1];
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkEG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0]-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkEFHG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkFH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2*wrap_[0];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1+2*wrap_[0];
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0]-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+2*wrap_[0]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1+2*wrap_[0]-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[0];
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexLinkH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-2*wrap_[1];
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-2*wrap_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-2*wrap_[2];
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-2*wrap_[2];
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]-2*wrap_[0];
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1-2*wrap_[0];
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleA(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1]+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1]+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1]+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0]+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0]+2*wrap_[1]+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleAB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1]+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1]+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1]+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1]+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1]+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1]+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1]+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1]+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1]+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleAC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0]+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleABDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleBD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0]+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0]+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0]+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleCD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[2];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[2];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[2];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[2];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+2*wrap_[2];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[2];
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[2];
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[2];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[2];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[2];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleAE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0]+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0]+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleAEFB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleBF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleAEGC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleABCDEFGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleBFHD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleCG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleGHDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleDH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0]+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0]+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleEF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+2*wrap_[1];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1+2*wrap_[1];
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+2*wrap_[1];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1+2*wrap_[1];
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+2*wrap_[1];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[1];
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+2*wrap_[1];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[1];
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+2*wrap_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[1];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[1];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+2*wrap_[1];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[1];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+2*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleEG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleEFHG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleFH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+2*wrap_[0];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+2*wrap_[0];
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+2*wrap_[0];
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+2*wrap_[0];
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1+2*wrap_[0];
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+2*wrap_[0];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+2*wrap_[0];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+2*wrap_[0];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+2*wrap_[0];
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+2*wrap_[0];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+2*wrap_[0];
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+2*wrap_[0];
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+2*wrap_[0];
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexTriangleH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 13: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 14: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 15: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 16: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 17: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 18: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 22: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 23: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 24: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 25: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 26: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 27: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 28: return p[0]*2+tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 29: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 30: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 31: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 32: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 33: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 34: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 35: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarA(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1]+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1]+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1]+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1]+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1]+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0]+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0]+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[1]+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1]+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarAB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1]+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1]+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1]+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1]+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1]+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1]+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1]+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1]+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1]+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1]+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1]+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1]+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1]+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1]+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarAC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0]+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0]+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarABDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarBD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0]+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0]+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarCD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[2];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[2];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[2];
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[2];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[2];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[2];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarAE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0]+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarAEFB(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarBF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarAEGC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarABCDEFGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarBFHD(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarCG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarGHDC(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarDH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarE(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0]+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0]+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0]+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarEF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarF(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+6*wrap_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1+6*wrap_[1];
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2+6*wrap_[1];
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[1];
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[1];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[1];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[1];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarEG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarEFHG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarFH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+6*wrap_[0];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2+6*wrap_[0];
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3+6*wrap_[0];
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+6*wrap_[0];
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+6*wrap_[0];
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+6*wrap_[0];
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+6*wrap_[0];
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+6*wrap_[0];
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5+6*wrap_[0];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getVertexStarH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    case 4: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    case 5: return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1];
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+1;
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 10: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 11: return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_x00(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+wrap_[1]*2;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+wrap_[2]*2;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+wrap_[1]*2+wrap_[2]*2;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xn0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+wrap_[2]*2;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+wrap_[2]*2;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xN0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1+wrap_[2]*2;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+wrap_[2]*2;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_x0n(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+wrap_[1]*2;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+wrap_[1]*2;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xnn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xNn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_x0N(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1+wrap_[1]*2;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1+wrap_[1]*2;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xnN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleL_xNN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_0y0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+wrap_[0]*2;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+wrap_[2]*2;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+wrap_[2]*2;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_ny0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+wrap_[2]*2;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+wrap_[2]*2;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_Ny0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1+wrap_[2]*2;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1+wrap_[2]*2;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_0yn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+wrap_[0]*2;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_nyn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_Nyn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_0yN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1+wrap_[0]*2;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7]+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_nyN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleH_NyN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_00z(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+wrap_[0]*2;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+wrap_[1]*2;
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+wrap_[1]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_n0z(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+wrap_[1]*2;
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+wrap_[1]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_N0z(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11]+wrap_[1]*2;
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5]+wrap_[1]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_0nz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+wrap_[0]*2;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_nnz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_Nnz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_0Nz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1+wrap_[0]*2;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1+wrap_[0]*2;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_nNz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleP_NNz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD1_xy0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1+wrap_[2]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD1_xyn(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD1_xyN(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD2_0yz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1+wrap_[0]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD2_nyz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD2_Nyz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD3_x0z(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 3: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7]+wrap_[1]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD3_xnz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 3: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD3_xNz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 3: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeTriangleD4_xyz(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkL(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==0) wrapYTop = wrap_[1];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==0) wrapZBack = wrap_[2];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapYBottom;
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapZFront;
    case 2: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]+wrapYTop;
    case 3: return esetshift_[1]+p[0]+1+(p[1]-1)*eshift_[4]+(p[2]-1)*eshift_[5]+wrapXRight+wrapYTop+wrapZBack;
    case 4: return esetshift_[0]+p[0]+1+(p[1]-1)*eshift_[2]+(p[2]-1)*eshift_[3]+wrapXRight+wrapYTop+wrapZBack;
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+(p[2]-1)*eshift_[13]+wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkH(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==0) wrapXLeft = wrap_[0];
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==0) wrapZBack = wrap_[2];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+(p[2]-1)*eshift_[1]+wrapZBack;
    case 1: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1+wrapXLeft+wrapYBottom+wrapZFront;
    case 2: return esetshift_[1]+p[0]-1+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapXLeft+wrapYBottom;
    case 3: return esetshift_[1]+p[0]+1+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+wrapXRight+wrapZBack;
    case 4: return esetshift_[5]+(p[0]-1)+p[1]*eshift_[12]+(p[2]-1)*eshift_[13]+wrapXLeft+wrapZBack;
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkP(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==0) wrapXLeft = wrap_[0];
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==0) wrapYTop = wrap_[1];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+(p[1]-1)*eshift_[0]+p[2]*eshift_[1]+wrapYTop;
    case 1: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1+wrapXLeft+wrapYBottom+wrapZFront;
    case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
    case 3: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]-1+wrapXLeft+wrapYTop;
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]-1+wrapXLeft+wrapZFront;
    case 5: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+1+wrapXRight+wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkD1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==0) wrapZBack = wrap_[2];
  switch(id){
    case 0: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11]+wrapZBack;
    case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]+wrapYBottom;
    case 2: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+(p[2]-1)*eshift_[9]+1+wrapXRight+wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkD2(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==0) wrapXLeft = wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];
    case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]-1+wrapXLeft+wrapYBottom;
    case 2: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]-1+wrapXLeft+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkD3(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==0) wrapYTop = wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7]+wrapYTop;
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]+wrapZFront;
    case 2: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 3: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+p[2]*eshift_[9]+1+wrapXRight+wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeLinkD4(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+p[2]*eshift_[1]+wrapYBottom;
    case 1: return p[0]+p[1]*eshift_[0]+(p[2]+1)*eshift_[1]+wrapZFront;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 3: return esetshift_[1]+p[0]+1+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapXRight+wrapYBottom;
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[0]+p[0]+1+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapXRight+wrapZFront;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarL(const SimplexId p[3], const int id) const{
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[1]==0) wrapYTop = 6*wrap_[1];
  if(p[2]==0) wrapZBack = 6*wrap_[2];
  switch(id){
    case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6;
    case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2;
    case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1+wrapYTop;
    case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3+wrapZBack;
    case 4: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+4+wrapYTop+wrapZBack;
    case 5: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5+wrapYTop+wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarH(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapZBack = 0;
  if(p[0]==0) wrapXLeft = 6*wrap_[0];
  if(p[2]==0) wrapZBack = 6*wrap_[2];
  switch(id){
    case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6;
    case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2+wrapZBack;
    case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3+wrapZBack;
    case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4+wrapXLeft+wrapZBack;
    case 4: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+1+wrapXLeft;
    case 5: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5+wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarP(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0]==0) wrapXLeft = 6*wrap_[0];
  if(p[1]==0) wrapYTop = 6*wrap_[1];
  switch(id){
    case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+(p[0]-1)*6+5+wrapXLeft+wrapYTop;
    case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+wrapYTop;
    case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1+wrapYTop;
    case 3: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2;
    case 4: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+3+wrapXLeft;
    case 5: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4+wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarD1(const SimplexId p[3], const int id) const{
  SimplexId wrapZBack = 0;
  if(p[2]==0) wrapZBack = 6*wrap_[2];
  switch(id){
    case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6;
    case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+1;
    case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3+wrapZBack;
    case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+4+wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarD2(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  if(p[0]==0) wrapXLeft = 6*wrap_[0];
  switch(id){
    case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6;
    case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2;
    case 2: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5+wrapXLeft;
    case 3: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4+wrapXLeft;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getEdgeStarD3(const SimplexId p[3], const int id) const{
  SimplexId wrapYTop = 0;
  if(p[1]==0) wrapYTop = 6*wrap_[1];
  switch(id){
    case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2;
    case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3;
    case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1+wrapYTop;
    case 3: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5+wrapYTop;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexF(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1+wrapXRight+wrapYBottom;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1+wrapXRight+wrapZFront;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexC(const SimplexId p[3], const int id) const{
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
      case 2: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0]+wrapYBottom+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
      case 2: return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0]+wrapYBottom+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexD1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0]+wrapYBottom+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1+wrapXRight+wrapYBottom;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0]+wrapYBottom+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexD2(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1+wrapXRight+wrapYBottom+wrapZFront;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleVertexD3(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1+wrapXRight+wrapZFront;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
      case 1: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
      case 2: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    }
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeF_0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 2: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeF_1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return p[0]/2+(p[1]+1)*eshift_[0]+p[2]*eshift_[1]+wrapYBottom;
    case 1: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3]+1+wrapXRight;
    case 2: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeH_0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 2: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeH_1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  SimplexId wrapZFront = 0;
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]/2+p[1]*eshift_[0]+(p[2]+1)*eshift_[1]+wrapZFront;
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5]+1+wrapXRight;
    case 2: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeC_0(const SimplexId p[3], const int id) const{
  SimplexId wrapYBottom = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 1: return esetshift_[1]+p[0]/2+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapYBottom;
    case 2: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeC_1(const SimplexId p[3], const int id) const{
  SimplexId wrapZFront = 0;
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapZFront;
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 2: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD1_0(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3]+1+wrapXRight;
    case 1: return esetshift_[4]+p[0]/2+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]+wrapYBottom;
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD1_1(const SimplexId p[3], const int id) const{
  SimplexId wrapZFront = 0;
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapZFront;
    case 1: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD2_0(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD2_1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]/2+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]+wrapYBottom+wrapZFront;
    case 1: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9]+1+wrapXRight;
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD3_0(const SimplexId p[3], const int id) const{
  SimplexId wrapYBottom = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return esetshift_[1]+p[0]/2+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapYBottom;
    case 1: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleEdgeD3_1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5]+1+wrapXRight;
    case 1: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]+wrapZFront;
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkF(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==0) wrapZBack = wrap_[2];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+wrapYBottom+wrapZFront;
    case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]-1)*vshift_[1]+1+wrapXRight+wrapZBack;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==0) wrapYTop = wrap_[1];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+wrapYBottom+wrapZFront;
    case 1: return p[0]/2+(p[1]-1)*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight+wrapYTop;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkC(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==0) wrapXLeft = wrap_[0];
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]-1+wrapXLeft+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkD1(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+1+wrapXRight+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+wrapYBottom;
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+1+wrapXRight+wrapYBottom+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkD2(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight+wrapYBottom;
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+1+wrapXRight+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+wrapYBottom;
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+wrapZFront;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleLinkD3(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]/2==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+wrapZFront;
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+1+wrapXRight+wrapYBottom+wrapZFront;
    }
  }
  else {
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight+wrapYBottom;
    }
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarF(const SimplexId p[3], const int id) const{
  SimplexId wrapZBack = 0;
  if(p[2]==0) wrapZBack = 6*wrap_[2];
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4+wrapZBack;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];
      case 1: return p[0]*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3+wrapZBack;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarH(const SimplexId p[3], const int id) const{
  SimplexId wrapYTop = 0;
  if(p[1]==0) wrapYTop = 6*wrap_[1];
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
      case 1: return (p[0]-1)*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5+wrapYTop;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
      case 1: return p[0]*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1+wrapYTop;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarC(const SimplexId p[3], const int id) const{
  SimplexId wrapXLeft = 0;
  if(p[0]<2) wrapXLeft = 6*wrap_[0];
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
      case 1: return (p[0]/2-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4+wrapXLeft;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];
      case 1: return (p[0]/2-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5+wrapXLeft;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarD1(const SimplexId p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarD2(const SimplexId p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;
    }
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTriangleStarD3(const SimplexId p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;
    }
  }
  else {
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;
    }
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexABCG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1];
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBCDG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+wrapYBottom;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1+wrapXRight+wrapYBottom;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexABEG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1];
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBEFG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+wrapZFront;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1+wrapXRight+wrapZFront;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBFGH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1+wrapXRight+wrapZFront;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1+wrapXRight+wrapYBottom+wrapZFront;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronVertexBDGH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1+wrapXRight;
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1+wrapXRight+wrapYBottom;
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+wrapYBottom+wrapZFront;
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1+wrapXRight+wrapYBottom+wrapZFront;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeABCG(const SimplexId p[3], const int id) const{
  SimplexId wrapYBottom = 0;
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 2: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapYBottom;
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];
    case 4: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBCDG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+p[2]*eshift_[1]+wrapYBottom;
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+p[2]*eshift_[3]+wrapXRight;
    case 2: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapYBottom;
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];
    case 4: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]+wrapYBottom;
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeABEG(const SimplexId p[3], const int id) const{
  SimplexId wrapZFront = 0;
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapZFront;
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 4: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBEFG(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+(p[2]+1)*eshift_[1]+wrapZFront;
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapZFront;
    case 2: return esetshift_[1]+(p[0]+1)+p[1]*eshift_[4]+p[2]*eshift_[5]+wrapXRight;
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]+wrapZFront;
    case 4: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBFGH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]+wrapYBottom+wrapZFront;
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]+wrapXRight+wrapZFront;
    case 2: return esetshift_[1]+(p[0]+1)+p[1]*eshift_[4]+p[2]*eshift_[5]+wrapXRight;
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]+wrapZFront;
    case 4: return esetshift_[3]+(p[0]+1)+p[1]*eshift_[8]+p[2]*eshift_[9]+wrapXRight;
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronEdgeBDGH(const SimplexId p[3], const int id) const{
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0]==nbvoxels_[0]) wrapXRight = -wrap_[0];
  if(p[1]==nbvoxels_[1]) wrapYBottom = -wrap_[1];
  if(p[2]==nbvoxels_[2]) wrapZFront = -wrap_[2];
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]+wrapYBottom+wrapZFront;
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+p[2]*eshift_[3]+wrapXRight;
    case 2: return esetshift_[1]+(p[0]+1)+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]+wrapXRight+wrapYBottom;
    case 3: return esetshift_[3]+(p[0]+1)+p[1]*eshift_[8]+p[2]*eshift_[9]+wrapXRight;
    case 4: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]+wrapYBottom;
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleABCG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBCDG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return (p[1]<nbvoxels_[1]) ? tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3] : tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]-wrap_[1]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleABEG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBEFG(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return (p[2]<nbvoxels_[2]) ? p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1] : p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]-wrap_[2]*2;
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBFGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return (p[0]<nbvoxels_[0]) ? tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+3 : tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+3-wrap_[0]*2;
    case 3: return (p[2]<nbvoxels_[2]) ? p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1 : p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1-wrap_[2]*2;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronTriangleBDGH(const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return (p[0]<nbvoxels_[0]) ? tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2 : tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2-wrap_[0]*2;
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 3: return (p[1]<nbvoxels_[1]) ? tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1 : tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1-wrap_[1]*2;
  }
  return -1;
}


inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborABCG(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t+1;
    case 1: return t+2;
    case 2: return p[0]>0 ? t-1 : t-1+wrap_[0]*6;
    case 3: return p[2]>0 ? t-tetshift_[1]+3 : t-tetshift_[1]+3+wrap_[2]*6;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBCDG(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+4;
    case 2: return p[2]>0 ? t-tetshift_[1]+3 : t-tetshift_[1]+3+wrap_[2]*6;
    case 3: return p[1]<nbvoxels_[1] ? t+tetshift_[0]+1 : t+tetshift_[0]+1-wrap_[1]*6;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborABEG(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t-2;
    case 1: return t+1;
    case 2: return p[0]>0 ? t-4 : t-4+wrap_[0]*6;
    case 3: return p[1]>0 ? t-tetshift_[0]-1 : t-tetshift_[0]-1+wrap_[1]*6;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBEFG(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+1;
    case 2: return p[1]>0 ? t-tetshift_[0]+2 : t-tetshift_[0]+2+wrap_[1]*6;
    case 3: return p[2]<nbvoxels_[2] ? t+tetshift_[1]-3 : t+tetshift_[1]-3-wrap_[2]*6;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBFGH(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+1;
    case 2: return p[0]<nbvoxels_[0] ? t+4 : t+4-wrap_[0]*6;
    case 3: return p[2]<nbvoxels_[2] ? t+tetshift_[1]-3 : t+tetshift_[1]-3-wrap_[2]*6;
  }
  return -1;
}

inline ttk::SimplexId ttk::PeriodicImplicitTriangulation::getTetrahedronNeighborBDGH(const SimplexId t, const SimplexId p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t-4;
    case 2: return p[0]<nbvoxels_[0] ? t+1 : t+1-wrap_[0]*6;
    case 3: return p[1]<nbvoxels_[1] ? t+tetshift_[0]-2 : t+tetshift_[0]-2-wrap_[1]*6;
  }
  return -1;
}

#endif // _PERIODICIMPLICITTRIANGULATION_H
