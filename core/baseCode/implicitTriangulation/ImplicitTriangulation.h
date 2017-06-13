/// \ingroup baseCode
/// \class ttk::ImplicitTriangulation
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date January 2016.
///
/// \brief ImplicitTriangulation is a class that provides time and memory
/// efficient traversal methods on triangulations of piecewise linear
/// manifolds represented by regular grids.
/// \sa Triangulation

#ifndef _IMPLICITTRIANGULATION_H
#define _IMPLICITTRIANGULATION_H

// base code includes
#include                  <AbstractTriangulation.h>

namespace ttk{

  class ImplicitTriangulation : public AbstractTriangulation{

    public:
      ImplicitTriangulation();
      ~ImplicitTriangulation();

      int getCellEdge(const int &cellId, const int &localEdgeId, int &edgeId) const;

      int getCellEdgeNumber(const int &cellId) const;

      const vector<vector<int>>* getCellEdges();

      int getCellNeighbor(const int &cellId, const int &localNeighborId, int &neighborId) const;

      int getCellNeighborNumber(const int &cellId) const;

      const vector<vector<int>>* getCellNeighbors();

      int getCellTriangle(const int &cellId, const int &localTriangleId, int &triangleId) const;

      int getCellTriangleNumber(const int &cellId) const{
        // NOTE: the output is always 4 here. let's keep the function in there
        // in case of further generalization to CW-complexes
        return 4;
      };

      const vector<vector<int>>* getCellTriangles();

      int getCellVertex(const int &cellId, const int &localVertexId, int &vertexId) const;

      int getCellVertexNumber(const int &cellId) const;

      int getDimensionality() const {return dimensionality_;};

      int getEdgeLink(const int &edgeId, const int &localLinkId, int &linkId) const;

      int getEdgeLinkNumber(const int &edgeId) const;

      const vector<vector<int> >* getEdgeLinks();

      int getEdgeStar(const int &edgeId, const int &localStarId, int &starId) const;

      int getEdgeStarNumber(const int &edgeId) const;

      const vector<vector<int>>* getEdgeStars();

      int getEdgeTriangle(const int &edgeId, const int &localTriangleId, int &triangleId) const;

      int getEdgeTriangleNumber(const int &edgeId) const;

      const vector<vector<int>>* getEdgeTriangles();

      int getEdgeVertex(const int &edgeId, const int &localVertexId, int &vertexId) const;

      const vector<pair<int,int>>* getEdges();

      int getNumberOfCells() const{return cellNumber_;};

      int getNumberOfEdges() const{return edgeNumber_;};

      int getNumberOfTriangles() const{return triangleNumber_;};

      int getNumberOfVertices() const{return vertexNumber_;};

      int getTetrahedronEdge(const int &tetId, const int &localEdgeId, int &edgeId) const;

      int getTetrahedronEdges(vector<vector<int>>& edges) const;

      int getTetrahedronTriangle(const int &tetId, const int &localTriangleId, int &triangleId) const;

      int getTetrahedronTriangles(vector<vector<int>>& triangles) const;

      int getTetrahedronNeighbor(const int &tetId, const int &localNeighborId, int &neighborId) const;

      int getTetrahedronNeighborNumber(const int &tetId) const;

      int getTetrahedronNeighbors(vector<vector<int>>& neighbors);

      int getTetrahedronVertex(const int& tetId, const int& localVertexId, int& vertexId) const;

      int getTriangleEdge(const int &triangleId, const int &localEdgeId, int &edgeId) const;

      int getTriangleEdgeNumber(const int &triangleId) const{
        // NOTE: the output is always 3 here. let's keep the function in there
        // in case of further generalization to CW-complexes
        return 3;
      }

      const vector<vector<int>>* getTriangleEdges();

      int getTriangleEdges(vector<vector<int>>& edges) const;

      int getTriangleLink(const int &triangleId, const int &localLinkId, int &linkId) const;

      int getTriangleLinkNumber(const int &triangleId) const;

      const vector<vector<int> >* getTriangleLinks();

      int getTriangleNeighbor(const int &triangleId, const int &localNeighborId, int &neighborId) const;

      int getTriangleNeighborNumber(const int &triangleId) const;

      int getTriangleNeighbors(vector<vector<int>>& neighbors);

      int getTriangleStar(const int &triangleId, const int &localStarId, int &starId) const;

      int getTriangleStarNumber(const int &triangleId) const;

      const vector<vector<int>>* getTriangleStars();

      int getTriangleVertex(const int &triangleId, const int &localVertexId, int &vertexId) const;

      const vector<vector<int>>* getTriangles();

      int getVertexEdge(const int &vertexId, const int &localEdgeId, int &edgeId) const;

      int getVertexEdgeNumber(const int &vertexId) const;

      const vector<vector<int>>* getVertexEdges();

      int getVertexLink(const int& vertexId, const int& localLinkId, int &linkId) const;

      int getVertexLinkNumber(const int &vertexId) const;

      const vector<vector<int> >* getVertexLinks();

      int getVertexNeighbor(const int &vertexId, const int &localNeighborId, int &neighborId) const;

      int getVertexNeighborNumber(const int &vertexId) const;

      const vector<vector<int>>* getVertexNeighbors();

      int getVertexPoint(const int &vertexId, float &x, float &y, float &z) const;

      int getVertexStar(const int &vertexId, const int &localStarId, int &starId) const;

      int getVertexStarNumber(const int &vertexId) const;

      const vector<vector<int>>* getVertexStars();

      int getVertexTriangle(const int &vertexId, const int &localTriangleId, int &triangleId) const;

      int getVertexTriangleNumber(const int &vertexId) const;

      const vector<vector<int>>* getVertexTriangles();

      bool isEdgeOnBoundary(const int &edgeId) const;

      bool isEmpty() const{return !vertexNumber_;};

      bool isTriangleOnBoundary(const int &triangleId) const;

      bool isVertexOnBoundary(const int &vertexId) const;

      int setInputGrid(const float &xOrigin, const float &yOrigin, const float &zOrigin,
          const float &xSpacing, const float &ySpacing, const float &zSpacing,
          const int &xDim, const int &yDim, const int &zDim);

    protected:
      int dimensionality_;//
      float origin_[3];//
      float spacing_[3];//
      int dimensions_[3];// dimensions
      int nbvoxels_[3];// nombre de voxels par axe

      // Vertex helper //
      int vshift_[2];// VertexShift

      // Edge helper //
      int esetdims_[7];// EdgeSetDimensions
      int esetshift_[7];// EdgeSetShift
      int eshift_[14];// EdgeShift

      // Triangle helper //
      int tsetdims_[6];// TriangleSetDimensions
      int tsetshift_[6];// TriangleSetShift
      int tshift_[12];// TriangleShift

      // Tetrahedron helper //
      int tetshift_[2];// TetrahedronShift

      int cellNumber_;// number of cells
      int vertexNumber_;// number of vertices
      int edgeNumber_;// number of edges
      int triangleNumber_;// number of triangles
      int tetrahedronNumber_;// number of tetrahedra

      // 2d helpers
      int Di_;
      int Dj_;

      // acceleration variables
      bool isAccelerated_;
      int mod_[2];
      int div_[2];

      // acceleration functions
      int checkAcceleration();
      bool isPowerOfTwo(unsigned int v, unsigned int& r);

      //\cond
      // 2D //
      void vertexToPosition2d(const int vertex, int p[2]) const;
      void edgeToPosition2d(const int edge, const int k, int p[2]) const;
      void triangleToPosition2d(const int triangle, int p[2]) const;

      int getVertexNeighbor2dA(const int v,const int id) const;
      int getVertexNeighbor2dB(const int v,const int id) const;
      int getVertexNeighbor2dC(const int v,const int id) const;
      int getVertexNeighbor2dD(const int v,const int id) const;
      int getVertexNeighbor2dAB(const int v,const int id) const;
      int getVertexNeighbor2dCD(const int v,const int id) const;
      int getVertexNeighbor2dAC(const int v,const int id) const;
      int getVertexNeighbor2dBD(const int v,const int id) const;
      int getVertexNeighbor2dABCD(const int v,const int id) const;

      int getVertexEdge2dA(const int p[2],const int id) const;
      int getVertexEdge2dB(const int p[2],const int id) const;
      int getVertexEdge2dC(const int p[2],const int id) const;
      int getVertexEdge2dD(const int p[2],const int id) const;
      int getVertexEdge2dAB(const int p[2],const int id) const;
      int getVertexEdge2dCD(const int p[2],const int id) const;
      int getVertexEdge2dAC(const int p[2],const int id) const;
      int getVertexEdge2dBD(const int p[2],const int id) const;
      int getVertexEdge2dABCD(const int p[2],const int id) const;

      int getVertexStar2dA(const int p[2], const int id) const;
      int getVertexStar2dB(const int p[2], const int id) const;
      int getVertexStar2dC(const int p[2], const int id) const;
      int getVertexStar2dD(const int p[2], const int id) const;
      int getVertexStar2dAB(const int p[2], const int id) const;
      int getVertexStar2dCD(const int p[2], const int id) const;
      int getVertexStar2dAC(const int p[2], const int id) const;
      int getVertexStar2dBD(const int p[2], const int id) const;
      int getVertexStar2dABCD(const int p[2], const int id) const;

      int getVertexLink2dA(const int p[2],const int id) const;
      int getVertexLink2dB(const int p[2],const int id) const;
      int getVertexLink2dC(const int p[2],const int id) const;
      int getVertexLink2dD(const int p[2],const int id) const;
      int getVertexLink2dAB(const int p[2],const int id) const;
      int getVertexLink2dCD(const int p[2],const int id) const;
      int getVertexLink2dAC(const int p[2],const int id) const;
      int getVertexLink2dBD(const int p[2],const int id) const;
      int getVertexLink2dABCD(const int p[2],const int id) const;

      int getEdgeTriangleL_x0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xN(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_0y(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_ny(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_Ny(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD1_xy(const int p[3], const int localTriangleId) const;

      int getEdgeLink2dL(const int p[2], const int id) const;
      int getEdgeLink2dH(const int p[2], const int id) const;
      int getEdgeLink2dD1(const int p[2], const int id) const;

      int getEdgeStar2dL(const int p[2], const int id) const;
      int getEdgeStar2dH(const int p[2], const int id) const;

      // 3D //
      void vertexToPosition(const int vertex, int p[3]) const;
      void edgeToPosition(const int edge, const int k, int p[3]) const;
      void triangleToPosition(const int triangle, const int k, int p[3]) const;
      void tetrahedronToPosition(const int tetrahedron, int p[3]) const;

      int getVertexNeighborA(const int v,const int id) const;
      int getVertexNeighborB(const int v,const int id) const;
      int getVertexNeighborC(const int v,const int id) const;
      int getVertexNeighborD(const int v,const int id) const;
      int getVertexNeighborE(const int v,const int id) const;
      int getVertexNeighborF(const int v,const int id) const;
      int getVertexNeighborG(const int v,const int id) const;
      int getVertexNeighborH(const int v,const int id) const;
      int getVertexNeighborAB(const int v,const int id) const;
      int getVertexNeighborCD(const int v,const int id) const;
      int getVertexNeighborEF(const int v,const int id) const;
      int getVertexNeighborGH(const int v,const int id) const;
      int getVertexNeighborAC(const int v,const int id) const;
      int getVertexNeighborBD(const int v,const int id) const;
      int getVertexNeighborEG(const int v,const int id) const;
      int getVertexNeighborFH(const int v,const int id) const;
      int getVertexNeighborAE(const int v,const int id) const;
      int getVertexNeighborBF(const int v,const int id) const;
      int getVertexNeighborCG(const int v,const int id) const;
      int getVertexNeighborDH(const int v,const int id) const;
      int getVertexNeighborABDC(const int v,const int id) const;
      int getVertexNeighborEFHG(const int v,const int id) const;
      int getVertexNeighborAEGC(const int v,const int id) const;
      int getVertexNeighborBFHD(const int v,const int id) const;
      int getVertexNeighborAEFB(const int v,const int id) const;
      int getVertexNeighborGHDC(const int v,const int id) const;
      int getVertexNeighborABCDEFGH(const int v,const int id) const;

      int getVertexEdgeA(const int p[3],const int id) const;
      int getVertexEdgeB(const int p[3],const int id) const;
      int getVertexEdgeC(const int p[3],const int id) const;
      int getVertexEdgeD(const int p[3],const int id) const;
      int getVertexEdgeE(const int p[3],const int id) const;
      int getVertexEdgeF(const int p[3],const int id) const;
      int getVertexEdgeG(const int p[3],const int id) const;
      int getVertexEdgeH(const int p[3],const int id) const;
      int getVertexEdgeAB(const int p[3],const int id) const;
      int getVertexEdgeCD(const int p[3],const int id) const;
      int getVertexEdgeEF(const int p[3],const int id) const;
      int getVertexEdgeGH(const int p[3],const int id) const;
      int getVertexEdgeAC(const int p[3],const int id) const;
      int getVertexEdgeBD(const int p[3],const int id) const;
      int getVertexEdgeEG(const int p[3],const int id) const;
      int getVertexEdgeFH(const int p[3],const int id) const;
      int getVertexEdgeAE(const int p[3],const int id) const;
      int getVertexEdgeBF(const int p[3],const int id) const;
      int getVertexEdgeCG(const int p[3],const int id) const;
      int getVertexEdgeDH(const int p[3],const int id) const;
      int getVertexEdgeABDC(const int p[3],const int id) const;
      int getVertexEdgeEFHG(const int p[3],const int id) const;
      int getVertexEdgeAEGC(const int p[3],const int id) const;
      int getVertexEdgeBFHD(const int p[3],const int id) const;
      int getVertexEdgeAEFB(const int p[3],const int id) const;
      int getVertexEdgeGHDC(const int p[3],const int id) const;
      int getVertexEdgeABCDEFGH(const int p[3],const int id) const;

      int getVertexTriangleA(const int p[3],const int id) const;
      int getVertexTriangleB(const int p[3],const int id) const;
      int getVertexTriangleC(const int p[3],const int id) const;
      int getVertexTriangleD(const int p[3],const int id) const;
      int getVertexTriangleE(const int p[3],const int id) const;
      int getVertexTriangleF(const int p[3],const int id) const;
      int getVertexTriangleG(const int p[3],const int id) const;
      int getVertexTriangleH(const int p[3],const int id) const;
      int getVertexTriangleAB(const int p[3],const int id) const;
      int getVertexTriangleCD(const int p[3],const int id) const;
      int getVertexTriangleEF(const int p[3],const int id) const;
      int getVertexTriangleGH(const int p[3],const int id) const;
      int getVertexTriangleAC(const int p[3],const int id) const;
      int getVertexTriangleBD(const int p[3],const int id) const;
      int getVertexTriangleEG(const int p[3],const int id) const;
      int getVertexTriangleFH(const int p[3],const int id) const;
      int getVertexTriangleAE(const int p[3],const int id) const;
      int getVertexTriangleBF(const int p[3],const int id) const;
      int getVertexTriangleCG(const int p[3],const int id) const;
      int getVertexTriangleDH(const int p[3],const int id) const;
      int getVertexTriangleABDC(const int p[3],const int id) const;
      int getVertexTriangleEFHG(const int p[3],const int id) const;
      int getVertexTriangleAEGC(const int p[3],const int id) const;
      int getVertexTriangleBFHD(const int p[3],const int id) const;
      int getVertexTriangleAEFB(const int p[3],const int id) const;
      int getVertexTriangleGHDC(const int p[3],const int id) const;
      int getVertexTriangleABCDEFGH(const int p[3],const int id) const;

      int getVertexLinkA(const int p[3],const int id) const;
      int getVertexLinkB(const int p[3],const int id) const;
      int getVertexLinkC(const int p[3],const int id) const;
      int getVertexLinkD(const int p[3],const int id) const;
      int getVertexLinkE(const int p[3],const int id) const;
      int getVertexLinkF(const int p[3],const int id) const;
      int getVertexLinkG(const int p[3],const int id) const;
      int getVertexLinkH(const int p[3],const int id) const;
      int getVertexLinkAB(const int p[3],const int id) const;
      int getVertexLinkCD(const int p[3],const int id) const;
      int getVertexLinkEF(const int p[3],const int id) const;
      int getVertexLinkGH(const int p[3],const int id) const;
      int getVertexLinkAC(const int p[3],const int id) const;
      int getVertexLinkBD(const int p[3],const int id) const;
      int getVertexLinkEG(const int p[3],const int id) const;
      int getVertexLinkFH(const int p[3],const int id) const;
      int getVertexLinkAE(const int p[3],const int id) const;
      int getVertexLinkBF(const int p[3],const int id) const;
      int getVertexLinkCG(const int p[3],const int id) const;
      int getVertexLinkDH(const int p[3],const int id) const;
      int getVertexLinkABDC(const int p[3],const int id) const;
      int getVertexLinkEFHG(const int p[3],const int id) const;
      int getVertexLinkAEGC(const int p[3],const int id) const;
      int getVertexLinkBFHD(const int p[3],const int id) const;
      int getVertexLinkAEFB(const int p[3],const int id) const;
      int getVertexLinkGHDC(const int p[3],const int id) const;
      int getVertexLinkABCDEFGH(const int p[3],const int id) const;

      int getVertexStarA(const int p[3], const int id) const;
      int getVertexStarB(const int p[3],const int id) const;
      int getVertexStarC(const int p[3],const int id) const;
      int getVertexStarD(const int p[3],const int id) const;
      int getVertexStarE(const int p[3],const int id) const;
      int getVertexStarF(const int p[3],const int id) const;
      int getVertexStarG(const int p[3],const int id) const;
      int getVertexStarH(const int p[3],const int id) const;
      int getVertexStarAB(const int p[3],const int id) const;
      int getVertexStarCD(const int p[3],const int id) const;
      int getVertexStarEF(const int p[3],const int id) const;
      int getVertexStarGH(const int p[3],const int id) const;
      int getVertexStarAC(const int p[3],const int id) const;
      int getVertexStarBD(const int p[3],const int id) const;
      int getVertexStarEG(const int p[3],const int id) const;
      int getVertexStarFH(const int p[3],const int id) const;
      int getVertexStarAE(const int p[3],const int id) const;
      int getVertexStarBF(const int p[3],const int id) const;
      int getVertexStarCG(const int p[3],const int id) const;
      int getVertexStarDH(const int p[3],const int id) const;
      int getVertexStarABDC(const int p[3],const int id) const;
      int getVertexStarEFHG(const int p[3],const int id) const;
      int getVertexStarAEGC(const int p[3],const int id) const;
      int getVertexStarBFHD(const int p[3],const int id) const;
      int getVertexStarAEFB(const int p[3],const int id) const;
      int getVertexStarGHDC(const int p[3],const int id) const;
      int getVertexStarABCDEFGH(const int p[3],const int id) const;

      int getEdgeTriangleL_x00(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_x0n(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_x0N(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xn0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xnn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xnN(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xN0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xNn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleL_xNN(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleH_0y0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_0yn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_0yN(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_ny0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_nyn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_nyN(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_Ny0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_Nyn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleH_NyN(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleP_00z(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_0nz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_0Nz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_n0z(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_nnz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_nNz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_N0z(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_Nnz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleP_NNz(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleD1_xy0(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD1_xyn(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD1_xyN(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleD2_0yz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD2_nyz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD2_Nyz(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleD3_x0z(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD3_xnz(const int p[3], const int localTriangleId) const;
      int getEdgeTriangleD3_xNz(const int p[3], const int localTriangleId) const;

      int getEdgeTriangleD4_xyz(const int p[3], const int localTriangleId) const;

      int getEdgeLinkL(const int p[3], const int id) const;
      int getEdgeLinkH(const int p[3], const int id) const;
      int getEdgeLinkP(const int p[3], const int id) const;
      int getEdgeLinkD1(const int p[3], const int id) const;
      int getEdgeLinkD2(const int p[3], const int id) const;
      int getEdgeLinkD3(const int p[3], const int id) const;
      int getEdgeLinkD4(const int p[3], const int id) const;

      int getEdgeStarL(const int p[3], const int id) const;
      int getEdgeStarH(const int p[3], const int id) const;
      int getEdgeStarP(const int p[3], const int id) const;
      int getEdgeStarD1(const int p[3], const int id) const;
      int getEdgeStarD2(const int p[3], const int id) const;
      int getEdgeStarD3(const int p[3], const int id) const;

      int getTriangleVertexF(const int p[3], const int id) const;
      int getTriangleVertexH(const int p[3], const int id) const;
      int getTriangleVertexC(const int p[3], const int id) const;
      int getTriangleVertexD1(const int p[3], const int id) const;
      int getTriangleVertexD2(const int p[3], const int id) const;
      int getTriangleVertexD3(const int p[3], const int id) const;

      int getTriangleEdgeF_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeF_1(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeH_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeH_1(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeC_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeC_1(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD1_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD1_1(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD2_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD2_1(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD3_0(const int p[3], const int localEdgeId) const;
      int getTriangleEdgeD3_1(const int p[3], const int localEdgeId) const;

      int getTriangleLinkF(const int p[3], const int id) const;
      int getTriangleLinkH(const int p[3], const int id) const;
      int getTriangleLinkC(const int p[3], const int id) const;
      int getTriangleLinkD1(const int p[3], const int id) const;
      int getTriangleLinkD2(const int p[3], const int id) const;
      int getTriangleLinkD3(const int p[3], const int id) const;

      int getTriangleStarF(const int p[3], const int id) const;
      int getTriangleStarH(const int p[3], const int id) const;
      int getTriangleStarC(const int p[3], const int id) const;
      int getTriangleStarD1(const int p[3], const int id) const;
      int getTriangleStarD2(const int p[3], const int id) const;
      int getTriangleStarD3(const int p[3], const int id) const;

      int getTetrahedronVertexABCG(const int p[3], const int id) const;
      int getTetrahedronVertexBCDG(const int p[3], const int id) const;
      int getTetrahedronVertexABEG(const int p[3], const int id) const;
      int getTetrahedronVertexBEFG(const int p[3], const int id) const;
      int getTetrahedronVertexBFGH(const int p[3], const int id) const;
      int getTetrahedronVertexBDGH(const int p[3], const int id) const;

      int getTetrahedronEdgeABCG(const int p[3], const int id) const;
      int getTetrahedronEdgeBCDG(const int p[3], const int id) const;
      int getTetrahedronEdgeABEG(const int p[3], const int id) const;
      int getTetrahedronEdgeBEFG(const int p[3], const int id) const;
      int getTetrahedronEdgeBFGH(const int p[3], const int id) const;
      int getTetrahedronEdgeBDGH(const int p[3], const int id) const;

      int getTetrahedronTriangleABCG(const int p[3], const int id) const;
      int getTetrahedronTriangleBCDG(const int p[3], const int id) const;
      int getTetrahedronTriangleABEG(const int p[3], const int id) const;
      int getTetrahedronTriangleBEFG(const int p[3], const int id) const;
      int getTetrahedronTriangleBFGH(const int p[3], const int id) const;
      int getTetrahedronTriangleBDGH(const int p[3], const int id) const;

      int getTetrahedronNeighborABCG(const int t, const int p[3], const int id) const;
      int getTetrahedronNeighborBCDG(const int t, const int p[3], const int id) const;
      int getTetrahedronNeighborABEG(const int t, const int p[3], const int id) const;
      int getTetrahedronNeighborBEFG(const int t, const int p[3], const int id) const;
      int getTetrahedronNeighborBFGH(const int t, const int p[3], const int id) const;
      int getTetrahedronNeighborBDGH(const int t, const int p[3], const int id) const;
      //\endcond
  };
}

inline void ImplicitTriangulation::vertexToPosition2d(const int vertex, int p[2]) const{
  if(isAccelerated_){
    p[0]=vertex&mod_[0];
    p[1]=vertex>>div_[0];
  }
  else{
    p[0]=vertex%vshift_[0];
    p[1]=vertex/vshift_[0];
  }
}

inline void ImplicitTriangulation::edgeToPosition2d(const int edge, const int k, int p[2]) const{
  const int e=(k)?edge-esetshift_[k-1]:edge;
  p[0]=e%eshift_[2*k];
  p[1]=e/eshift_[2*k];
}

inline void ImplicitTriangulation::triangleToPosition2d(const int triangle, int p[2]) const{
  p[0]=triangle%tshift_[0];
  p[1]=triangle/tshift_[0];
}

inline int ImplicitTriangulation::getVertexNeighbor2dA(const int v,const int id) const{
  //V(a)={b,c}
  switch(id){
    case 0: return v+1;//b
    case 1: return v+vshift_[0];//c
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dB(const int v,const int id) const{
  //V(b)={a,c,d}
  switch(id){
    case 0: return v-1;//a
    case 1: return v+vshift_[0];//d
    case 2: return v+vshift_[0]-1;//c
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dC(const int v,const int id) const{
  //V(c)={a,b,d}
  switch(id){
    case 0: return v+1;//d
    case 1: return v-vshift_[0];//a
    case 2: return v-vshift_[0]+1;//b
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dD(const int v,const int id) const{
  //V(d)={c,b}
  switch(id){
    case 0: return v-1;//c
    case 1: return v-vshift_[0];//b
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dAB(const int v,const int id) const{
  //V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id){
    case 0: return v-1;//V(b)::a
    case 1: return v+vshift_[0]-1;//V(b)::c
    case 2: return v+vshift_[0];//V(b)::d
    case 3: return v+1;//V(a)::b
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dCD(const int v,const int id) const{
  //V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id){
    case 0: return v-1;//V(d)::c
    case 1: return v-vshift_[0];//V(c)::a
    case 2: return v-vshift_[0]+1;//V(c)::b
    case 3: return v+1;//V(c)::d
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dAC(const int v,const int id) const{
  //V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id){
    case 0: return v-vshift_[0];//V(c)::{a}
    case 1: return v-vshift_[0]+1;//V(c)::{b}
    case 2: return v+1;//V(c)::{d}
    case 3: return v+vshift_[0];//V(a)::{c}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dBD(const int v,const int id) const{
  //V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id){
    case 0: return v+vshift_[0]-1;//V(b)::{c}
    case 1: return v+vshift_[0];//V(b)::{d}
    case 2: return v-vshift_[0];//V(d)::{b}
    case 3: return v-1;//V(d)::{c}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighbor2dABCD(const int v,const int id) const{
  //V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
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

inline int ImplicitTriangulation::getVertexEdge2dA(const int p[2],const int id) const{
  //V(a)={b,c}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2];//ac-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dB(const int p[2],const int id) const{
  //V(b)={a,c,d}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]-1;//ba-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2];//bd-H
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;//bc-D1
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dC(const int p[2],const int id) const{
  //V(c)={a,b,d}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//ca-H
    case 1: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];//cb-D1
    case 2: return p[0]+p[1]*eshift_[0];//cd-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dD(const int p[2],const int id) const{
  //V(d)={c,b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]-1;//dc-L
    case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//db-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dAB(const int p[2],const int id) const{
  //V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]-1;//ba-L
    case 1: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2];//bd-H
    case 3: return p[0]+p[1]*eshift_[0];//ab-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dCD(const int p[2],const int id) const{
  //V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//ca-H
    case 1: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];//cb-D1
    case 2: return p[0]+p[1]*eshift_[0];//cd-L
    case 3: return p[0]+p[1]*eshift_[0]-1;//dc-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dAC(const int p[2],const int id) const{
  //V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//ca-H
    case 1: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];//cb-D1
    case 2: return p[0]+p[1]*eshift_[0];//cd-L
    case 3: return esetshift_[0]+p[0]+p[1]*eshift_[2];//ac-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dBD(const int p[2],const int id) const{
  //V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id){
    case 0: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;//bc-D1
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2];//bd-H
    case 2: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//db-H
    case 3: return p[0]+p[1]*eshift_[0]-1;//dc-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdge2dABCD(const int p[2],const int id) const{
  //V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{c}+V(b)::{c}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2];//db-H
    case 1: return p[0]+p[1]*eshift_[0]-1;//dc-L
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4];//cb-D1
    case 3: return p[0]+p[1]*eshift_[0];//cd-L
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2];//ac-H
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]-1;//bc-D1
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dA(const int p[2], const int id) const{
  return p[0]*2+p[1]*tshift_[0];
}

inline int ImplicitTriangulation::getVertexStar2dB(const int p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dC(const int p[2], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0];
    case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dD(const int p[2], const int id) const{
  return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
}

inline int ImplicitTriangulation::getVertexStar2dAB(const int p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dCD(const int p[2], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0];
    case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 2: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dAC(const int p[2], const int id) const{
  switch(id){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0];
    case 1: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 2: return p[0]*2+p[1]*tshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dBD(const int p[2], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0];
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 2: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStar2dABCD(const int p[2], const int id) const{
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

inline int ImplicitTriangulation::getVertexLink2dA(const int p[2],const int id) const{
  return esetshift_[1]+p[0]+p[1]*eshift_[4];//D1::bc
}

inline int ImplicitTriangulation::getVertexLink2dB(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;//H::ac
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;//L::ab
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dC(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;//H::ac
    case 1: return p[0]+(p[1]-1)*eshift_[0];//L::ab
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dD(const int p[2],const int id) const{
  return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;//D1::bc
}

inline int ImplicitTriangulation::getVertexLink2dAB(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;//H::ac
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;//L::ab
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];//D1::bc
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dCD(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;//H::ac
    case 1: return p[0]+(p[1]-1)*eshift_[0];//L::ab
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;//D1::bc
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dAC(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;//H::ac
    case 1: return p[0]+(p[1]-1)*eshift_[0];//L::ab
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];//D1::bc
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dBD(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;//H::ac
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;//L::ab
    case 2: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;//D1::bc
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLink2dABCD(const int p[2],const int id) const{
  switch(id){
    case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]-1;//H::ac
    case 1: return p[0]+(p[1]+1)*eshift_[0]-1;//L::ab
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4];//D1::bc
    case 3: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+1;//H::ac
    case 4: return p[0]+(p[1]-1)*eshift_[0];//L::ab
    case 5: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]-1;//D1::bc
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_x0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Di_]*2;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+(p[Dj_]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Di_]*2+(p[Dj_]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_0y(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Dj_]*tshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_ny(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return (p[Di_]-1)*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_Ny(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return (p[Di_]-1)*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD1_xy(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[Di_]*2+p[Dj_]*tshift_[0];
    case 1: return p[Di_]*2+p[Dj_]*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLink2dL(const int p[2], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[Dj_]){
    switch(id){
      case 0: return p[0]+(p[1]+1)*vshift_[0];
      case 1: return p[0]+(p[1]-1)*vshift_[0]+1;
    }
  }
  else if(p[1]==0) return p[0]+vshift_[0];
  else return p[0]+(p[1]-1)*vshift_[0]+1;
  return -1;
}

inline int ImplicitTriangulation::getEdgeLink2dH(const int p[2], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[Di_]){
    switch(id){
      case 0: return p[0]+p[1]*vshift_[0]+1;
      case 1: return p[0]+(p[1]+1)*vshift_[0]-1;
    }
  }
  else if(p[0]==0) return p[1]*vshift_[0]+1;
  else return p[0]+(p[1]+1)*vshift_[0]-1;
  return -1;
}

inline int ImplicitTriangulation::getEdgeLink2dD1(const int p[2], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0];
    case 1: return p[0]+(p[1]+1)*vshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStar2dL(const int p[2], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[Dj_]){
    if(id==0) return p[0]*2+p[1]*tshift_[0];
    else return p[0]*2+(p[1]-1)*tshift_[0]+1;
  }
  else if(p[1]==0) return p[0]*2+p[1]*tshift_[0];
  else return p[0]*2+(p[1]-1)*tshift_[0]+1;

  return -1;
}

inline int ImplicitTriangulation::getEdgeStar2dH(const int p[2], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[Di_]){
    if(id==0) return p[0]*2+p[1]*tshift_[0];
    else return (p[0]-1)*2+p[1]*tshift_[0]+1;
  }
  else if(p[0]==0) return p[0]*2+p[1]*tshift_[0];
  else return (p[0]-1)*2+p[1]*tshift_[0]+1;
  return -1;
}

inline void ImplicitTriangulation::vertexToPosition(const int vertex, int p[3]) const{
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

inline void ImplicitTriangulation::edgeToPosition(const int edge, const int k, int p[3]) const{
  const int e=(k)?edge-esetshift_[k-1]:edge;
  p[0]=e%eshift_[2*k];
  p[1]=(e%eshift_[2*k+1])/eshift_[2*k];
  p[2]=e/eshift_[2*k+1];
}

inline void ImplicitTriangulation::triangleToPosition(const int triangle, const int k, int p[3]) const{
  const int t=(k)?triangle-tsetshift_[k-1]:triangle;
  p[0]=t%tshift_[2*k];
  p[1]=(t%tshift_[2*k+1])/tshift_[2*k];
  p[2]=t/tshift_[2*k+1];
}

inline void ImplicitTriangulation::tetrahedronToPosition(const int tetrahedron, int p[3]) const{
  p[0]=(tetrahedron%tetshift_[0])/6;
  p[1]=(tetrahedron%tetshift_[1])/tetshift_[0];
  p[2]=tetrahedron/tetshift_[1];
}

inline int ImplicitTriangulation::getVertexNeighborA(const int v,const int id) const{
  // V(a)={b,c,e,g}
  switch(id){
    case 0: return v+1;//b
    case 1: return v+vshift_[0];//c
    case 2: return v+vshift_[1];//e
    case 3: return v+vshift_[0]+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborB(const int v,const int id) const{
  // V(b)={a,c,d,e,f,g,h}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborC(const int v,const int id) const{
  // V(c)={a,b,d,g}
  switch(id){
    case 0: return v-vshift_[0];//a
    case 1: return v+1-vshift_[0];//b
    case 2: return v+1;//d
    case 3: return v+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborD(const int v,const int id) const{
  // V(d)={b,c,g,h}
  switch(id){
    case 0: return v-vshift_[0];//b
    case 1: return v-1;//c
    case 2: return v-1+vshift_[1];//g
    case 3: return v+vshift_[1];//h
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborE(const int v,const int id) const{
  // V(e)={a,b,f,g}
  switch(id){
    case 0: return v-vshift_[1];//a
    case 1: return v+1-vshift_[1];//b
    case 2: return v+1;//f
    case 3: return v+vshift_[0];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborF(const int v,const int id) const{
  // V(f)={b,e,g,h}
  switch(id){
    case 0: return v-vshift_[1];//b
    case 1: return v-1;//e
    case 2: return v-1+vshift_[0];//g
    case 3: return v+vshift_[0];//h
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborG(const int v,const int id) const{
  // V(g)={a,b,c,d,e,f,h}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborH(const int v,const int id) const{
  // V(h)={b,d,f,g}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//b
    case 1: return v-vshift_[1];//d
    case 2: return v-vshift_[0];//f
    case 3: return v-1;//g
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborAB(const int v,const int id) const{
  // V(ab)=V(b)+V(a)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v+1;//V(a)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborCD(const int v,const int id) const{
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id){
    case 0: return v-vshift_[0];//b
    case 1: return v-1;//c
    case 2: return v-1+vshift_[1];//g
    case 3: return v+vshift_[1];//h
    case 4: return v+1-vshift_[0];//V(c)::{b}
    case 5: return v+1;//V(c)::{d}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborEF(const int v,const int id) const{
  // V(ef)=V(f)+V(e)::{b,f}
  switch(id){
    case 0: return v-vshift_[1];//b
    case 1: return v-1;//e
    case 2: return v-1+vshift_[0];//g
    case 3: return v+vshift_[0];//h
    case 4: return v+1-vshift_[1];//V(e)::{b}
    case 5: return v+1;//V(e)::{f}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborGH(const int v,const int id) const{
  // V(gh)=V(g)+V(h)::{g}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v-1;//V(h)::{g}
  }

  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborAC(const int v,const int id) const{
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id){
    case 0: return v-vshift_[0];//a
    case 1: return v+1-vshift_[0];//b
    case 2: return v+1;//d
    case 3: return v+vshift_[1];//g
    case 4: return v+vshift_[0];//V(a)::{c}
    case 5: return v+vshift_[0]+vshift_[1];//V(a)::{c}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborBD(const int v,const int id) const{
  // V(bd)=V(b)+V(d)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v-vshift_[0];//V(d)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborEG(const int v,const int id) const{
  //V(eg)=V(g)+V(e)::{g}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v+vshift_[0];//V(e)::{g}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborFH(const int v,const int id) const{
  //V(fh)=V(f)+V(h)::{b,f}
  switch(id){
    case 0: return v-vshift_[1];//b
    case 1: return v-1;//e
    case 2: return v-1+vshift_[0];//g
    case 3: return v+vshift_[0];//h
    case 4: return v-vshift_[0]-vshift_[1];//V(h)::{b}
    case 5: return v-vshift_[0];//V(h)::{f}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborAE(const int v,const int id) const{
  //V(ae)=V(a)+V(e)::{a,b}
  switch(id){
    case 0: return v+1;//b
    case 1: return v+vshift_[0];//c
    case 2: return v+vshift_[1];//e
    case 3: return v+vshift_[0]+vshift_[1];//g
    case 4: return v-vshift_[1];//V(e)::{a}
    case 5: return v+1-vshift_[1];//V(e)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborBF(const int v,const int id) const{
  //V(bf)=V(b)+V(f)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v-vshift_[1];//V(f)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborCG(const int v,const int id) const{
  //V(cg)=V(g)+V(c)::{g}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v+vshift_[1];//V(c)::{g}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborDH(const int v,const int id) const{
  //V(dh)=V(d)+V(h)::{b,d}
  switch(id){
    case 0: return v-vshift_[0];//b
    case 1: return v-1;//c
    case 2: return v-1+vshift_[1];//g
    case 3: return v+vshift_[1];//h
    case 4: return v-vshift_[0]-vshift_[1];//V(h)::{b}
    case 5: return v-vshift_[1];//V(h)::{d}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborABDC(const int v,const int id) const{
  //V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v-vshift_[0];//V(d)::{b}
    case 8: return v+1-vshift_[0];//V(c)::{b}
    case 9: return v+1;//V(a)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborEFHG(const int v,const int id) const{
  //V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v-1;//V(h)::{g}
    case 8: return v-1+vshift_[0];//V(f)::{g}
    case 9: return v+vshift_[0];//V(f)::{h}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborAEGC(const int v,const int id) const{
  //V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v+vshift_[0];//V(a)::{c}
    case 8: return v+vshift_[0]+vshift_[1];//V(a)::{g}
    case 9: return v+vshift_[1];//V(c)::{g}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborBFHD(const int v,const int id) const{
  //V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v-vshift_[1];//V(f)::{b}
    case 8: return v-vshift_[0]-vshift_[1];//V(h)::{b}
    case 9: return v-vshift_[0];//V(d)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborAEFB(const int v,const int id) const{
  //V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id){
    case 0: return v-1;//a
    case 1: return v-1+vshift_[0];//c
    case 2: return v+vshift_[0];//d
    case 3: return v-1+vshift_[1];//e
    case 4: return v+vshift_[1];  //f
    case 5: return v-1+vshift_[0]+vshift_[1];  //g
    case 6: return v+vshift_[0]+vshift_[1];//h
    case 7: return v+1;//V(a)::{b}
    case 8: return v+1-vshift_[1];//V(e)::{b}
    case 9: return v-vshift_[1];//V(f)::{b}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborGHDC(const int v,const int id) const{
  //V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v-1;//V(h)::{g}
    case 8: return v-1+vshift_[1];//V(d)::{g}
    case 9: return v+vshift_[1];//V(d)::{h}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexNeighborABCDEFGH(const int v,const int id) const{
  //V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id){
    case 0: return v-vshift_[0]-vshift_[1];//a
    case 1: return v+1-vshift_[0]-vshift_[1];//b
    case 2: return v-vshift_[1];//c
    case 3: return v+1-vshift_[1];//d
    case 4: return v-vshift_[0];//e
    case 5: return v+1-vshift_[0];//f
    case 6: return v+1;//h
    case 7: return v-1+vshift_[1];//V(d)::{g}
    case 8: return v+vshift_[1];//V(d)::{h}
    case 9: return v-1;//V(h)::{g}
    case 10: return v-1+vshift_[0];//V(b)::{c}
    case 11: return v+vshift_[0];//V(b)::{d}
    case 12: return v-1+vshift_[0]+vshift_[1];//V(b)::{g}
    case 13: return v+vshift_[0]+vshift_[1];//V(b)::{h}
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeA(const int p[3],const int id) const{
  // V(a)={b,c,e,g}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeB(const int p[3],const int id) const{
  // V(b)={a,c,d,e,f,g,h}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeC(const int p[3],const int id) const{
  // V(c)={a,b,d,g}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ca-H
    case 1: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//cb-D1
    case 2: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//cd-L
    case 3: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//cg-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeD(const int p[3],const int id) const{
  // V(d)={b,c,g,h}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//dc-L
    case 2: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//dg-D3
    case 3: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//dh-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeE(const int p[3],const int id) const{
  // V(e)={a,b,f,g}
  switch(id){
    case 0: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//ea-P
    case 1: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//eb-D3
    case 2: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ef-L
    case 3: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//eg-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeF(const int p[3],const int id) const{
  // V(f)={b,e,g,h}
  switch(id){
    case 0: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//fe-L
    case 2: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//fg-D1
    case 3: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//fh-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeG(const int p[3],const int id) const{
  // V(g)={a,b,c,d,e,f,h}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeH(const int p[3],const int id) const{
  // V(h)={b,d,f,g}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//hb-D2
    case 1: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//hd-P
    case 2: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//hf-H
    case 3: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//hg-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeAB(const int p[3],const int id) const{
  // V(ab)=V(b)+V(a)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[11]+p[2]*eshift_[12]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeCD(const int p[3],const int id) const{
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//dc-L
    case 2: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//dg-D3
    case 3: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//dh-P
    case 4: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//cb-D1
    case 5: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//cd-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeEF(const int p[3],const int id) const{
  // V(fe)=V(f)+V(e)::{b,f}
  switch(id){
    case 0: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//fe-L
    case 2: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//fg-D1
    case 3: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//fh-H
    case 4: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//eb-D3
    case 5: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ef-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeGH(const int p[3],const int id) const{
  // V(gh)=V(g)+V(h)::{g}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//hg-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeAC(const int p[3],const int id) const{
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ca-H
    case 1: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//cb-D1
    case 2: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//cd-L
    case 3: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//cg-P
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 5: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeBD(const int p[3],const int id) const{
  // V(bd)=V(b)+V(d)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeEG(const int p[3],const int id) const{
  //V(eg)=V(g)+V(e)::{g}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//eg-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeFH(const int p[3],const int id) const{
  //V(fh)=V(f)+V(h)::{b,f}
  switch(id){
    case 0: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//fe-L
    case 2: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//fg-D1
    case 3: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//fh-H
    case 4: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//hb-D2
    case 5: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//hf-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeAE(const int p[3],const int id) const{
  //V(ae)=V(a)+V(e)::{a,b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//ea-P
    case 5: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//eb-D3
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeBF(const int p[3],const int id) const{
  //V(bf)=V(b)+V(f)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeCG(const int p[3],const int id) const{
  //V(cg)=V(g)+V(c)::{g}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//cg-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeDH(const int p[3],const int id) const{
  //V(dh)=V(d)+V(h)::{b,d}
  switch(id){
    case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
    case 1: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//dc-L
    case 2: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//dg-D3
    case 3: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//dh-P
    case 4: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//hb-D2
    case 5: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//hd-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeABDC(const int p[3],const int id) const{
  //V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
    case 8: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//cb-D1
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeEFHG(const int p[3],const int id) const{
  //V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//hg-L
    case 8: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//fg-D1
    case 9: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//fh-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeAEGC(const int p[3],const int id) const{
  //V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 8: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 9: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//cg-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeBFHD(const int p[3],const int id) const{
  //V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
    case 8: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//hb-D2
    case 9: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//db-H
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeAEFB(const int p[3],const int id) const{
  //V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//ba-L
    case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 2: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//be-D3
    case 4: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//bf-P
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 6: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
    case 7: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 8: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//eb-D3
    case 9: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//fb-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeGHDC(const int p[3],const int id) const{
  //V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//hg-L
    case 8: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//dg-D3
    case 9: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//dh-P
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexEdgeABCDEFGH(const int p[3],const int id) const{
  //V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id){
    case 0: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+(p[2]-1)*eshift_[9];//ga-D2
    case 1: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+(p[2]-1)*eshift_[13];//gb-D4
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];//gc-P
    case 3: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];//gd-D3
    case 4: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3];//ge-H
    case 5: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];//gf-D1
    case 6: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//gh-L
    case 7: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11]-1;//dg-D3
    case 8: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//dh-P
    case 9: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1]-1;//hg-L
    case 10: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7]-1;//bc-D1
    case 11: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//bd-H
    case 12: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13]-1;//bg-D4
    case 13: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//bh-D2
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleA(const int p[3],const int id) const{
  switch(id){
    case 0: return 0;
    case 1: return tsetshift_[0];
    case 2: return tsetshift_[1];
    case 3: return tsetshift_[3];
    case 4: return tsetshift_[1]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleB(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2;
    case 2: return tsetshift_[2]+(p[0]-1)*2;
    case 3: return tsetshift_[3]+(p[0]-1)*2+1;
    case 4: return tsetshift_[1]+p[0]*2;
    case 5: return tsetshift_[1]+p[0]*2+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2;
    case 10: return tsetshift_[0]+(p[0]-1)*2;
    case 11: return (p[0]-1)*2;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleC(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[1]-1)*tshift_[0];
    case 1: return (p[1]-1)*tshift_[0]+1;
    case 2: return tsetshift_[4]+(p[1]-1)*tshift_[10];
    case 3: return tsetshift_[0]+p[1]*tshift_[2];
    case 4: return tsetshift_[1]+(p[1]-1)*tshift_[4];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleD(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6];
    case 2: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4];
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleE(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[2]-1)*tshift_[3];
    case 1: return tsetshift_[2]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[2]-1)*tshift_[5]+1;
    case 3: return p[2]*tshift_[1];
    case 4: return tsetshift_[0]+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+(p[2]-1)*tshift_[3]+1;
    case 1: return (p[0]-1)*2+p[2]*tshift_[1];
    case 2: return (p[0]-1)*2+p[2]*tshift_[1]+1;
    case 3: return tsetshift_[4]+(p[0]-1)*2+(p[2]-1)*tshift_[11]+1;
    case 4: return tsetshift_[1]+p[0]*2+(p[2]-1)*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 1: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleAB(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2;
    case 2: return tsetshift_[2]+(p[0]-1)*2;
    case 3: return tsetshift_[3]+(p[0]-1)*2+1;
    case 4: return tsetshift_[1]+p[0]*2;
    case 5: return tsetshift_[1]+p[0]*2+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2;
    case 10: return tsetshift_[0]+(p[0]-1)*2;
    case 11: return (p[0]-1)*2;
    case 12: return p[0]*2;
    case 13: return tsetshift_[0]+p[0]*2;
    case 14: return tsetshift_[3]+p[0]*2;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleCD(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6];
    case 2: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4];
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+1;
    case 5: return p[0]*2+(p[1]-1)*tshift_[0];
    case 6: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 7: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10];
    case 8: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleEF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+(p[2]-1)*tshift_[3]+1;
    case 1: return (p[0]-1)*2+p[2]*tshift_[1];
    case 2: return (p[0]-1)*2+p[2]*tshift_[1]+1;
    case 3: return tsetshift_[4]+(p[0]-1)*2+(p[2]-1)*tshift_[11]+1;
    case 4: return tsetshift_[1]+p[0]*2+(p[2]-1)*tshift_[5]+1;
    case 5: return p[0]*2+tsetshift_[0]+(p[2]-1)*tshift_[3];
    case 6: return p[0]*2+tsetshift_[2]+(p[2]-1)*tshift_[7]+1;
    case 7: return p[0]*2+p[2]*tshift_[1];
    case 8: return p[0]*2+tsetshift_[0]+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleGH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 1: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 5: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 6: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 7: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 8: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 9: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 10: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 11: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 12: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 13: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 14: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleAC(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[1]-1)*tshift_[0];
    case 1: return (p[1]-1)*tshift_[0]+1;
    case 2: return tsetshift_[4]+(p[1]-1)*tshift_[10];
    case 3: return tsetshift_[0]+p[1]*tshift_[2];
    case 4: return tsetshift_[1]+(p[1]-1)*tshift_[4];
    case 5: return p[1]*tshift_[0];
    case 6: return tsetshift_[1]+p[1]*tshift_[4];
    case 7: return tsetshift_[3]+p[1]*tshift_[8];
    case 8: return tsetshift_[1]+p[1]*tshift_[4]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleBD(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6];
    case 2: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4];
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+1;
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10];
    case 7: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6];
    case 8: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+1;
    case 9: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4];
    case 10: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+1;
    case 11: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+1;
    case 12: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+1;
    case 13: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8];
    case 14: return (p[0]-1)*2+p[1]*tshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleEG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 12: return tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 13: return tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 14: return p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleFH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 1: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 7: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleAE(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[2]-1)*tshift_[3];
    case 1: return tsetshift_[2]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[2]-1)*tshift_[5]+1;
    case 3: return p[2]*tshift_[1];
    case 4: return tsetshift_[0]+(p[2]-1)*tshift_[3]+1;
    case 5: return tsetshift_[0]+p[2]*tshift_[3];
    case 6: return tsetshift_[1]+p[2]*tshift_[5];
    case 7: return tsetshift_[3]+p[2]*tshift_[9];
    case 8: return tsetshift_[1]+p[2]*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleBF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+(p[0]-1)*2+(p[2]-1)*tshift_[3]+1;
    case 1: return (p[0]-1)*2+p[2]*tshift_[1];
    case 2: return (p[0]-1)*2+p[2]*tshift_[1]+1;
    case 3: return tsetshift_[4]+(p[0]-1)*2+(p[2]-1)*tshift_[11]+1;
    case 4: return tsetshift_[1]+p[0]*2+(p[2]-1)*tshift_[5]+1;
    case 5: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11];
    case 6: return tsetshift_[2]+(p[0]-1)*2+p[2]*tshift_[7];
    case 7: return tsetshift_[3]+(p[0]-1)*2+p[2]*tshift_[9]+1;
    case 8: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5];
    case 9: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5]+1;
    case 10: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11]+1;
    case 11: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3]+1;
    case 12: return tsetshift_[2]+(p[0]-1)*2+p[2]*tshift_[7]+1;
    case 13: return tsetshift_[3]+(p[0]-1)*2+p[2]*tshift_[9];
    case 14: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleCG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 12: return tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 13: return tsetshift_[0]+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 14: return tsetshift_[1]+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleDH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 1: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 5: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 6: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 7: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 8: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleABDC(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4];
    case 5: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2];
    case 11: return (p[0]-1)*2+p[1]*tshift_[0];
    case 12: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+1;
    case 13: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6];
    case 14: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4];
    case 15: return p[0]*2+(p[1]-1)*tshift_[0];
    case 16: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 17: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10];
    case 18: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2];
    case 19: return p[0]*2+p[1]*tshift_[0];
    case 20: return p[0]*2+tsetshift_[3]+p[1]*tshift_[8];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleEFHG(const int p[3],const int id) const{
  switch(id){
    case 0: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 12: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 13: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 14: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 15: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 16: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 17: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 18: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 19: return p[0]*2+tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 20: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleAEGC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return (p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 12: return tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 13: return tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 14: return p[1]*tshift_[0]+p[2]*tshift_[1];
    case 15: return tsetshift_[0]+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 16: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 17: return tsetshift_[3]+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 18: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 19: return tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 20: return tsetshift_[1]+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleBFHD(const int p[3],const int id) const{
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
    case 12: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 13: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 14: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 15: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 16: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 17: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 18: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 19: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 20: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleAEFB(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*2+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[2]*tshift_[7];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5];
    case 5: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5]+1;
    case 6: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11]+1;
    case 7: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3]+1;
    case 8: return tsetshift_[2]+(p[0]-1)*2+p[2]*tshift_[7]+1;
    case 9: return tsetshift_[3]+(p[0]-1)*2+p[2]*tshift_[9];
    case 10: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3];
    case 11: return (p[0]-1)*2+p[2]*tshift_[1];
    case 12: return tsetshift_[0]+(p[0]-1)*2+(p[2]-1)*tshift_[3]+1;
    case 13: return tsetshift_[4]+(p[0]-1)*2+(p[2]-1)*tshift_[11]+1;
    case 14: return tsetshift_[1]+p[0]*2+(p[2]-1)*tshift_[5]+1;
    case 15: return p[0]*2+tsetshift_[0]+(p[2]-1)*tshift_[3];
    case 16: return p[0]*2+tsetshift_[2]+(p[2]-1)*tshift_[7]+1;
    case 17: return p[0]*2+p[2]*tshift_[1];
    case 18: return p[0]*2+tsetshift_[0]+(p[2]-1)*tshift_[3]+1;
    case 19: return p[0]*2+tsetshift_[0]+p[2]*tshift_[3];
    case 20: return p[0]*2+tsetshift_[3]+p[2]*tshift_[9];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleGHDC(const int p[3],const int id) const{
  switch(id){
    case 0: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];
    case 1: return p[0]*2+tsetshift_[2]+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];
    case 3: return p[0]*2+tsetshift_[1]+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9];
    case 5: return p[0]*2+tsetshift_[3]+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 6: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11];
    case 7: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;
    case 8: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1];
    case 9: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 10: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3];
    case 11: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 12: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 13: return (p[0]-1)*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 14: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 15: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
    case 16: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 17: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 18: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 19: return p[0]*2+tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 20: return p[0]*2+tsetshift_[0]+p[1]*tshift_[2]+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexTriangleABCDEFGH(const int p[3],const int id) const{
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

inline int ImplicitTriangulation::getVertexLinkA(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkB(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 4: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 1: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkD(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkE(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 1: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 5: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 1: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkAB(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkCD(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 1: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
    case 2: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 3: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkEF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkGH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 5: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
    case 6: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 7: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkAC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 3: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkBD(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 4: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
    case 6: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 7: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkEG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 2: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 3: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 4: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 5: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 6: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 7: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkFH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 1: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
    case 2: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 3: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkAE(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkBF(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 4: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
    case 6: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 7: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkCG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 1: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
    case 2: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 3: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 4: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 5: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 6: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 7: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkDH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
    case 2: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 3: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkABDC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
    case 8: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 9: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
    case 10: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 11: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkEFHG(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 2: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
    case 4: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 5: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
    case 6: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 7: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 8: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 9: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 10: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 11: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkAEGC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 4: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 5: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 6: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 7: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 8: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 9: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
    case 10: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 11: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkBFHD(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 4: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 5: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
    case 6: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 7: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
    case 8: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 9: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkAEFB(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 4: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 5: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
    case 6: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 7: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 8: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 9: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 10: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 11: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkGHDC(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 1: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 2: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 3: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 4: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 5: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
    case 6: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 7: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
    case 8: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 9: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
    case 10: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 11: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexLinkABCDEFGH(const int p[3],const int id) const{
  switch(id){
    case 0: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;//D1::beg
    case 2: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5];//C::acg
    case 3: return tsetshift_[1]+(p[0]-1)*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;//C::aeg
    case 4: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];//H::abe
    case 5: return tsetshift_[0]+(p[0]-1)*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;//H::bef
    case 6: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];//F::abc
    case 7: return (p[0]-1)*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;//F::bcd
    case 8: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9];//D2::abg
    case 9: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];//D1::bdg
    case 10: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];//D3::bcg
    case 11: return tsetshift_[3]+(p[0]-1)*2+(p[1]-1)*tshift_[8]+p[2]*tshift_[9]+1;//D2::bgh
    case 12: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9];//D2::abg
    case 13: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 14: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;//D1::beg
    case 15: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+(p[2]-1)*tshift_[9]+1;//D2::bgh
    case 16: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5];//C::acg
    case 17: return tsetshift_[1]+(p[0]+1)*2+(p[1]-1)*tshift_[4]+(p[2]-1)*tshift_[5]+1;//C::aeg
    case 18: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3];//H::abe
    case 19: return tsetshift_[0]+p[0]*2+(p[1]-1)*tshift_[2]+(p[2]-1)*tshift_[3]+1;//H::bef
    case 20: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1];//F::abc
    case 21: return p[0]*2+(p[1]-1)*tshift_[0]+(p[2]-1)*tshift_[1]+1;//F::bcd
    case 22: return tsetshift_[4]+(p[0]-1)*2+(p[1]-1)*tshift_[10]+(p[2]-1)*tshift_[11]+1;//D3::bfg
    case 23: return tsetshift_[2]+(p[0]-1)*2+(p[1]-1)*tshift_[6]+(p[2]-1)*tshift_[7];//D1::bdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarA(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//abcg
    case 1: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//abeg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarB(const int p[3], const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;
  return -1;
}

inline int ImplicitTriangulation::getVertexStarC(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//abcg
    case 1: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//bcdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarD(const int p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//bcdg
    case 1: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarE(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//abeg
    case 1: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarF(const int p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//befg
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//bfgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarG(const int p[3], const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;
  return -1;
}

inline int ImplicitTriangulation::getVertexStarH(const int p[3], const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//bfgh
    case 1: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarAB(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 7: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarCD(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 1: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
    case 2: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 3: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarEF(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 1: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
    case 2: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 3: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarGH(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 7: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarAC(const int p[3],const int id) const{
  switch(id){
    case 0: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 1: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
    case 2: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 3: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarBD(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 7: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarEG(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 7: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarFH(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 1: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
    case 2: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 3: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarAE(const int p[3],const int id) const{
  switch(id){
    case 0: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 1: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
    case 2: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 3: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarBF(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 7: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarCG(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarDH(const int p[3],const int id) const{
  switch(id){
    case 0: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 1: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
    case 2: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 3: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarABDC(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 7: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
    case 10: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 11: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarEFHG(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 7: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
    case 8: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 9: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
    case 10: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 11: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarAEGC(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 7: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
    case 8: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 9: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
    case 10: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 11: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarBFHD(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 7: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
    case 8: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 9: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
    case 10: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 11: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarAEFB(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  switch(id){
    case 6: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 7: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
    case 8: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 9: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
    case 10: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 11: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarGHDC(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id;//tet(g)
  switch(id){
    case 6: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 7: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
    case 8: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 9: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
    case 10: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 11: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getVertexStarABCDEFGH(const int p[3],const int id) const{
  if(id>=0 && id<=5)
    return (p[0]-1)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+id;//tet(b)
  if(id>=6 && id<=11)
    return p[0]*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+id-6;//tet(g)
  switch(id){
    case 12: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1];//tet(a)::abcg
    case 13: return p[0]*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//tet(a)::abeg
    case 14: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1];//tet(c)::abcg
    case 15: return p[0]*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(c)::bcdg
    case 16: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//tet(d)::bcdg
    case 17: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//tet(d)::bdgh
    case 18: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+2;//tet(e)::abeg
    case 19: return p[0]*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(e)::befg
    case 20: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//tet(f)::befg
    case 21: return (p[0]-1)*6+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(f)::bfgh
    case 22: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//tet(h)::bfgh
    case 23: return (p[0]-1)*6+(p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+5;//tet(h)::bdgh
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_x00(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2;
    case 1: return tsetshift_[0]+p[0]*2;
    case 2: return tsetshift_[3]+p[0]*2;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_x0n(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[2]*tshift_[1];
    case 1: return tsetshift_[0]+p[0]*2+p[2]*tshift_[3];
    case 2: return tsetshift_[3]+p[0]*2+p[2]*tshift_[9];
    case 3: return tsetshift_[0]+p[0]*2+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_x0N(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[2]*tshift_[1];
    case 1: return tsetshift_[0]+p[0]*2+(p[2]-1)*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xn0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[1]*tshift_[0];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2];
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8];
    case 3: return p[0]*2+(p[1]-1)*tshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xnn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 2: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 3: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 4: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 5: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xnN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 2: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 3: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xN0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xNn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 2: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleL_xNN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+(p[1]-1)*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+(p[2]-1)*tshift_[3]+1;
    case 2: return tsetshift_[3]+p[0]*2+(p[1]-1)*tshift_[8]+(p[2]-1)*tshift_[9]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_0y0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[1]*tshift_[0];
    case 1: return tsetshift_[1]+p[1]*tshift_[4];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_0yn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 1: return tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return p[1]*tshift_[0]+p[2]*tshift_[1];
    case 3: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_0yN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 1: return tsetshift_[2]+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_ny0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4];
    case 3: return p[0]*2+p[1]*tshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_nyn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 4: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 5: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_nyN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+(p[2]-1)*tshift_[7]+1;
    case 2: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 3: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_Ny0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_Nyn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[2]+(p[0]-1)*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 3: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleH_NyN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+(p[2]-1)*tshift_[5]+1;
    case 1: return (p[0]-1)*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_00z(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+p[2]*tshift_[3];
    case 1: return tsetshift_[1]+p[2]*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_0nz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[0]+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 3: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_0Nz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[4]+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[0]+p[1]*tshift_[2]+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_n0z(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_nnz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 4: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_nNz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
    case 2: return tsetshift_[4]+p[0]*2+(p[1]-1)*tshift_[10]+p[2]*tshift_[11];
    case 3: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_N0z(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[2]*tshift_[5]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_Nnz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[4]+(p[0]-1)*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleP_NNz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+(p[0]-1)*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 1: return tsetshift_[1]+p[0]*2+(p[1]-1)*tshift_[4]+p[2]*tshift_[5];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD1_xy0(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[1]*tshift_[0];
    case 1: return p[0]*2+p[1]*tshift_[0]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD1_xyn(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD1_xyN(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 2: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+(p[2]-1)*tshift_[11]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD2_0yz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+p[1]*tshift_[8]+p[2]*tshift_[9];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD2_nyz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD2_Nyz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 2: return tsetshift_[3]+(p[0]-1)*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD3_x0z(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+p[0]*2+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+p[2]*tshift_[7]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD3_xnz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 3: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD3_xNz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return tsetshift_[2]+p[0]*2+(p[1]-1)*tshift_[6]+p[2]*tshift_[7];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeTriangleD4_xyz(const int p[3], const int localTriangleId) const{
  switch(localTriangleId){
    case 0: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 4: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 5: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkL(const int p[3], const int id) const{
  if(p[2]==0 and p[1]==0){
    switch(id){
      case 0: return esetshift_[1]+p[0]+eshift_[4]; //CG
      case 1: return esetshift_[0]+p[0]+eshift_[3]; //EG
    }
  }
  else if(p[2]==0 and p[1]==nbvoxels_[1]) return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]; //BG
  else if(p[2]==nbvoxels_[2] and p[1]==0) return esetshift_[5]+p[0]+(p[2]-1)*eshift_[13]; //BG
  else if(p[2]==nbvoxels_[2] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+(p[2]-1)*eshift_[5]+1; //BF
      case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+(p[2]-1)*eshift_[3]+1; //BD
    }
  }
  else if(p[2]==0 and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]; //CG
      case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+eshift_[3]; //EG
      case 2: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]; //BG
    }
  }
  else if(p[2]==nbvoxels_[2] and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+(p[2]-1)*eshift_[5]+1; //BF
      case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+(p[2]-1)*eshift_[3]+1; //BD
      case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12]+(p[2]-1)*eshift_[13]; //BG
   }
  }
  else if(p[2]>0 and p[2]<nbvoxels_[2] and p[1]==0){
    switch(id){
      case 0: return esetshift_[1]+p[0]+p[2]*eshift_[5]+eshift_[4]; //CG
      case 1: return esetshift_[0]+p[0]+(p[2]+1)*eshift_[3]; //EG
      case 2: return esetshift_[5]+p[0]+(p[2]-1)*eshift_[13]; //BG
    }
  }
  else if(p[2]>0 and p[2]<nbvoxels_[2] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]-1)*eshift_[4]+(p[2]-1)*eshift_[5]+1; //BF
      case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+(p[2]-1)*eshift_[3]+1; //BD
      case 2: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]; //BG
    }
  }
  else{
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5]; //CG
      case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]; //EG
      case 2: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13];
      case 3: return esetshift_[1]+p[0]+1+(p[1]-1)*eshift_[4]+(p[2]-1)*eshift_[5]; //CG
      case 4: return esetshift_[0]+p[0]+1+(p[1]-1)*eshift_[2]+(p[2]-1)*eshift_[3]; //EG
      case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+(p[2]-1)*eshift_[13];
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkH(const int p[3], const int id) const{
  if(p[0]==0 and p[2]==0) return esetshift_[5]+p[1]*eshift_[12];
  else if(p[0]==nbvoxels_[0] and p[2]==0){
    switch(id){
      case 0: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]-1;
      case 1: return p[0]+(p[1]+1)*eshift_[0]+eshift_[1]-1;
    }
  }
  else if(p[0]==0 and p[2]==nbvoxels_[2]){
    switch(id){
      case 0: return p[1]*eshift_[0]+(p[2]-1)*eshift_[1];
      case 1: return esetshift_[1]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+1;
    }
  }
  else if(p[0]==nbvoxels_[0] and p[2]==nbvoxels_[2]) return esetshift_[5]+p[0]-1+p[1]*eshift_[12]+(p[2]-1)*eshift_[13];
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[2]==0){
    switch(id){
      case 0: return p[0]-1+(p[1]+1)*eshift_[0]+eshift_[1];
      case 1: return esetshift_[1]+p[0]-1+(p[1]+1)*eshift_[4];
      case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12];
    }
  }
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[2]==nbvoxels_[2]){
    switch(id){
      case 0: return esetshift_[5]+p[0]-1+p[1]*eshift_[12]+(p[2]-1)*eshift_[13];
      case 1: return p[0]+p[1]*eshift_[0]+(p[2]-1)*eshift_[1];
      case 2: return esetshift_[1]+p[0]+1+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
    }
  }
  else if(p[0]==0 and p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return esetshift_[1]+p[1]*eshift_[4]+(p[2]-1)*eshift_[5]+1;
      case 1: return esetshift_[5]+p[1]*eshift_[12]+p[2]*eshift_[13];
      case 2: return p[1]*eshift_[0]+(p[2]-1)*eshift_[1];
    }
  }
  else if(p[0]==nbvoxels_[0] and p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return esetshift_[5]+p[0]-1+p[1]*eshift_[12]+(p[2]-1)*eshift_[13];
      case 1: return esetshift_[1]+p[0]-1+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];
      case 2: return p[0]-1+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1];
    }
  }
  else{
    switch(id){
      case 0: return p[0]+p[1]*eshift_[0]+(p[2]-1)*eshift_[1];
      case 1: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1;
      case 2: return esetshift_[1]+p[0]-1+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];
      case 3: return esetshift_[1]+p[0]+1+p[1]*eshift_[4]+(p[2]-1)*eshift_[5];
      case 4: return esetshift_[5]+(p[0]-1)+p[1]*eshift_[12]+(p[2]-1)*eshift_[13];
      case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkP(const int p[3], const int id) const{
  if(p[0]==0 and p[1]==0) return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
  else if(p[0]==0 and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+1;
      case 1: return p[0]+(p[1]-1)*eshift_[0]+p[2]*eshift_[1];
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]==0){
    switch(id){
      case 0: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]-1;
      case 1: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1;
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]==nbvoxels_[1]) return esetshift_[5]+p[0]-1+(p[1]-1)*eshift_[12]+p[2]*eshift_[13];
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[1]==0){
    switch(id){
      case 0: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1;
      case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]-1;
      case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
    }
  }
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return p[0]+(p[1]-1)*eshift_[0]+p[2]*eshift_[1];
      case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+1;
      case 2: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]-1;
    }
  }
  else if(p[0]==0 and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[0]+(p[1]-1)*eshift_[0]+p[2]*eshift_[1];
      case 1: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+1;
      case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]-1;
      case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]-1;
      case 2: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1;
    }
  }
  else{
    switch(id){
      case 0: return p[0]+(p[1]-1)*eshift_[0]+p[2]*eshift_[1];
      case 1: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1]-1;
      case 2: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];
      case 3: return esetshift_[5]+p[0]+(p[1]-1)*eshift_[12]+p[2]*eshift_[13]-1;
      case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3]-1;
      case 5: return esetshift_[0]+p[0]+(p[1]-1)*eshift_[2]+p[2]*eshift_[3]+1;
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkD1(const int p[3], const int id) const{
  if(p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
      case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11];
      case 2: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
      case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+(p[2]-1)*eshift_[9]+1;
    }
  }
  else if(p[2]==0){
    switch(id){
      case 0: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
      case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11];
    }
  }
  else{
    switch(id){
      case 0: return esetshift_[3]+p[0]+p[1]*eshift_[8]+(p[2]-1)*eshift_[9]+1;
      case 1: return esetshift_[4]+p[0]+p[1]*eshift_[10]+(p[2]-1)*eshift_[11];
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkD2(const int p[3], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[0]){
    switch(id){
      case 0: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];
      case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]-1;
      case 2: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];
      case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]-1;
    }
  }
  else if(p[0]==0){
    switch(id){
      case 0: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];
      case 1: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];
    }
  }
  else{
    switch(id){
      case 0: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7]-1;
      case 1: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11]-1;
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkD3(const int p[3], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
      case 1: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7];
      case 2: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
      case 3: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+p[2]*eshift_[9]+1;
    }
  }
  else if(p[1]==0){
    switch(id){
      case 0: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7];
      case 1: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];
    }
  }
  else{
    switch(id){
      case 0: return esetshift_[2]+p[0]+(p[1]-1)*eshift_[6]+p[2]*eshift_[7];
      case 1: return esetshift_[3]+p[0]+(p[1]-1)*eshift_[8]+p[2]*eshift_[9]+1;
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeLinkD4(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+p[2]*eshift_[1];
    case 1: return p[0]+p[1]*eshift_[0]+(p[2]+1)*eshift_[1];
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 3: return esetshift_[1]+p[0]+1+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];
    case 4: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 5: return esetshift_[0]+p[0]+1+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStarL(const int p[3], const int id) const{
  if(p[2]==0 and p[1]==0){
    switch(id){
      case 0: return p[0]*6; //ABCG
      case 1: return p[0]*6+2; //ABEG
    }
  }
  else if(p[2]==0 and p[1]==nbvoxels_[1]) return (p[1]-1)*tetshift_[0]+p[0]*6+1;//BCDG
  else if(p[2]==nbvoxels_[2] and p[1]==0) return (p[2]-1)*tetshift_[1]+p[0]*6+3;//BEFG
  else if(p[2]==nbvoxels_[2] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return (p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+p[0]*6+4; //BFGH
      case 1: return (p[1]-1)*tetshift_[0]+(p[2]-1)*tetshift_[1]+p[0]*6+5; //BDGH
    }
  }
  else if(p[2]==0 and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 2: return (p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
    }
  }
  else if(p[2]==nbvoxels_[2] and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 1: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+4; //BFGH
      case 2: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5; //BDGH
    }
  }
  else if(p[2]>0 and p[2]<nbvoxels_[2] and p[1]==0){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[0]*6; //ABCG
      case 1: return p[2]*tetshift_[1]+p[0]*6+2; //ABEG
      case 2: return (p[2]-1)*tetshift_[1]+p[0]*6+3; //BEFG
    }
  }
  else if(p[2]>0 and p[2]<nbvoxels_[2] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 1: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+4; //BFGH
      case 2: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5; //BDGH
    }
  }
  else{
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 4: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+4; //BFGH
      case 5: return (p[2]-1)*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5; //BDGH
    }
  }

  return -1;
}

inline int ImplicitTriangulation::getEdgeStarH(const int p[3], const int id) const{
  if(p[0]==0 and p[2]==0) return p[1]*tetshift_[0]; //ABCG
  else if(p[0]==nbvoxels_[0] and p[2]==0){
    switch(id){
      case 0: return p[1]*tetshift_[0]+(p[0]-1)*6+1; //BCDG
      case 1: return p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
    }
  }
  else if(p[0]==0 and p[2]==nbvoxels_[2]){
    switch(id){
      case 0: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+2; //ABEG
      case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+3; //BEFG
    }
  }
  else if(p[0]==nbvoxels_[0] and p[2]==nbvoxels_[2]) return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4;//BFGH
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[2]==0){
    switch(id){
      case 0: return p[1]*tetshift_[0]+(p[0]-1)*6+1; //BCDG
      case 1: return p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
      case 2: return p[1]*tetshift_[0]+p[0]*6; //ABCG
    }
  }
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[2]==nbvoxels_[2]){
    switch(id){
      case 0: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
    }
  }
  else if(p[0]==0 and p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+2; //ABEG
      case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+3; //BEFG
      case 2: return p[2]*tetshift_[1]+p[1]*tetshift_[0]; //ABCG
    }
  }
  else if(p[0]==nbvoxels_[0] and p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+1; //BCDG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
      case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
    }
  }
  else{
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
      case 4: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+1; //BCDG
      case 5: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStarP(const int p[3], const int id) const{
  if(p[0]==0 and p[1]==0) return p[2]*tetshift_[1]+2;//ABEG
  else if(p[0]==0 and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]; //ABCG
      case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+1; //BCDG
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]==0){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[0]-1)*6+3; //BEFG
      case 1: return p[2]*tetshift_[1]+(p[0]-1)*6+4; //BFGH
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]==nbvoxels_[1]) return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+(p[0]-1)*6+5;//BDGH
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[1]==0){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[0]-1)*6+3; //BEFG
      case 1: return p[2]*tetshift_[1]+(p[0]-1)*6+4; //BFGH
      case 2: return p[2]*tetshift_[1]+p[0]*6+2; //ABEG
    }
  }
  else if(p[0]>0 and p[0]<nbvoxels_[0] and p[1]==nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+(p[0]-1)*6+5; //BDGH
    }
  }
  else if(p[0]==0 and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]; //ABCG
      case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+1; //BCDG
      case 2: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+2; //ABEG
    }
  }
  else if(p[0]==nbvoxels_[0] and p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+3; //BEFG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
      case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+(p[0]-1)*6+5; //BDGH
    }
  }
  else{
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+(p[0]-1)*6+5; //BDGH
      case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6; //ABCG
      case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 3: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 4: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+3; //BEFG
      case 5: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStarD1(const int p[3], const int id) const{
  if(p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+1; //BCDG
      case 2: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 3: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+4; //BFGH
    }
  }
  else if(p[2]==0){
    switch(id){
      case 0: return p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[1]*tetshift_[0]+p[0]*6+1; //BCDG
    }
  }
  else{
    switch(id){
      case 0: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 1: return (p[2]-1)*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+4; //BFGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStarD2(const int p[3], const int id) const{
  if(p[0]>0 and p[0]<nbvoxels_[0]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6; //ABCG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 2: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
      case 3: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
    }
  }
  else if(p[0]==0){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]; //ABCG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+2; //ABEG
    }
  }
  else{
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+5; //BDGH
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+(p[0]-1)*6+4; //BFGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getEdgeStarD3(const int p[3], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
      case 2: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 3: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5; //BDGH
    }
  }
  else if(p[1]==0){
    switch(id){
      case 0: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+2; //ABEG
      case 1: return p[2]*tetshift_[1]+p[1]*tetshift_[0]+p[0]*6+3; //BEFG
    }
  }
  else{
    switch(id){
      case 0: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+1; //BCDG
      case 1: return p[2]*tetshift_[1]+(p[1]-1)*tetshift_[0]+p[0]*6+5; //BDGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexF(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1;
  }
  else{
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexH(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];
  }
  else{
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexC(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1];
    else if(id==1) return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];
    else return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0];
  }
  else{
    if(id==0) return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1];
    else if(id==1) return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];
    else return (p[0]/2)+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexD1(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0];
  }
  else{
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+vshift_[0];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexD2(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];
  }
  else{
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleVertexD3(const int p[3], const int id) const{
  if(p[0]%2){
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1;
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];
  }
  else{
    if(id==0) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
    else if(id==1) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];
    else return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeF_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 2: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeF_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+(p[1]+1)*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3]+1;
    case 2: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeH_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 2: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeH_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+p[1]*eshift_[0]+(p[2]+1)*eshift_[1];
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5]+1;
    case 2: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeC_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3];
    case 1: return esetshift_[1]+p[0]/2+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];
    case 2: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeC_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];
    case 1: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5];
    case 2: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD1_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+p[2]*eshift_[3]+1;
    case 1: return esetshift_[4]+p[0]/2+(p[1]+1)*eshift_[10]+p[2]*eshift_[11];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD1_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[0]+p[0]/2+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];
    case 1: return esetshift_[4]+p[0]/2+p[1]*eshift_[10]+p[2]*eshift_[11];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD2_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+p[1]*eshift_[0]+p[2]*eshift_[1];
    case 1: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD2_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return p[0]/2+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1];
    case 1: return esetshift_[3]+p[0]/2+p[1]*eshift_[8]+p[2]*eshift_[9]+1;
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD3_0(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[1]+p[0]/2+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];
    case 1: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+p[2]*eshift_[7];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleEdgeD3_1(const int p[3], const int localEdgeId) const{
  switch(localEdgeId){
    case 0: return esetshift_[1]+p[0]/2+p[1]*eshift_[4]+p[2]*eshift_[5]+1;
    case 1: return esetshift_[2]+p[0]/2+p[1]*eshift_[6]+(p[2]+1)*eshift_[7];
    case 2: return esetshift_[5]+p[0]/2+p[1]*eshift_[12]+p[2]*eshift_[13];
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkF(const int p[3], const int id) const{
  if(p[2]>0 and p[2]<nbvoxels_[2]){
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]-1)*vshift_[1]+1;
    }
  }
  else if(p[2]==0) return p[0]/2+(p[1]+1)*vshift_[0]+vshift_[1];
  else return p[0]/2+p[1]*vshift_[0]+(p[2]-1)*vshift_[1]+1;

  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkH(const int p[3], const int id) const{
  if(p[1]>0 and p[1]<nbvoxels_[1]){
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1];
      case 1: return p[0]/2+(p[1]-1)*vshift_[0]+p[2]*vshift_[1]+1;
    }
  }
  else if(p[1]==0) return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1];
  else return p[0]/2+(p[1]-1)*vshift_[0]+p[2]*vshift_[1]+1;

  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkC(const int p[3], const int id) const{
  if(p[0]>1 and p[0]<(dimensions_[0]*2-2)){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]-1;
    }
  }
  else if(p[0]<2) return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1]+1;
  else return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]-1;

  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkD1(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+1;
    }
  }
  else{
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+1;
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkD2(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+1;
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1]+1;
    }
  }
  else{
    switch(id){
      case 0: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1];
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleLinkD3(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+(p[2]+1)*vshift_[1];
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+(p[2]+1)*vshift_[1]+1;
    }
  }
  else{
    switch(id){
      case 0: return p[0]/2+p[1]*vshift_[0]+p[2]*vshift_[1];
      case 1: return p[0]/2+(p[1]+1)*vshift_[0]+p[2]*vshift_[1]+1;
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleStarF(const int p[3], const int id) const{
  if(p[0]%2){
    if(p[2]>0 and p[2]<nbvoxels_[2]){
      switch(id){
        case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
        case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//BFGH
      }
    }
    else if(p[2]==0) return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
    else return (p[0]-1)*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+4;//BFGH
  }
  else{
    if(p[2]>0 and p[2]<nbvoxels_[2]){
      switch(id){
        case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
        case 1: return p[0]*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//BEFG
      }
    }
    else if(p[2]==0) return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
    else return p[0]*3+p[1]*tetshift_[0]+(p[2]-1)*tetshift_[1]+3;//BEFG
  }

  return -1;
}

inline int ImplicitTriangulation::getTriangleStarH(const int p[3], const int id) const{
  if(p[0]%2){
    if(p[1]>0 and p[1]<nbvoxels_[1]){
      switch(id){
        case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;//BEFG
        case 1: return (p[0]-1)*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
      }
    }
    else if(p[1]==0) return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;//BEFG
    else return (p[0]-1)*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
  }
  else{
    if(p[1]>0 and p[1]<nbvoxels_[1]){
      switch(id){
        case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
        case 1: return p[0]*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
      }
    }
    else if(p[1]==0) return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
    else return p[0]*3+(p[1]-1)*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleStarC(const int p[3], const int id) const{
  if(p[0]%2){
    if(p[0]>1 and p[0]<(dimensions_[0]*2-2)){
      switch(id){
        case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
        case 1: return ((p[0]-2)/2)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;//BFGH
      }
    }
    else if(p[0]<2) return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
    else return ((p[0]-2)/2)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;//BFGH
  }
  else{
    if(p[0]>1 and p[0]<(dimensions_[0]*2-2)){
      switch(id){
        case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
        case 1: return ((p[0]-1)/2)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
      }
    }
    else if(p[0]<2) return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
    else return ((p[0]-1)/2)*6+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
  }

  return -1;
}

inline int ImplicitTriangulation::getTriangleStarD1(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;//BEFG
    }
  }
  else{
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleStarD2(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+5;//BDGH
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;//BFGH
    }
  }
  else{
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+2;//ABEG
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTriangleStarD3(const int p[3], const int id) const{
  if(p[0]%2){
    switch(id){
      case 0: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+3;//BEFG
      case 1: return (p[0]-1)*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+4;//BFGH
    }
  }
  else{
    switch(id){
      case 0: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1];//ABCG
      case 1: return p[0]*3+p[1]*tetshift_[0]+p[2]*tetshift_[1]+1;//BCDG
    }
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexABCG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1];//a
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];//c
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexBCDG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0];//c
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1;//d
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexABEG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1];//a
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];//e
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexBEFG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1];//e
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1;//f
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexBFGH(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[1]+1;//f
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1;//h
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronVertexBDGH(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+1;//b
    case 1: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+1;//d
    case 2: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1];//g
    case 3: return p[0]+p[1]*vshift_[0]+p[2]*vshift_[1]+vshift_[0]+vshift_[1]+1;//h
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeABCG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 2: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];//bc-D1
    case 4: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeBCDG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 2: return esetshift_[1]+p[0]+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+p[2]*eshift_[7];//bc-D1
    case 4: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11];//be-D3
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeABEG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+p[2]*eshift_[1];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];//ac-H
    case 2: return esetshift_[1]+p[0]+p[1]*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[3]+p[0]+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 4: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];//be-D3
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeBEFG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+p[1]*eshift_[0]+(p[2]+1)*eshift_[1];//ab-L
    case 1: return esetshift_[0]+p[0]+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];//ac-H
    case 2: return esetshift_[1]+(p[0]+1)+p[1]*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7];//bc-D1
    case 4: return esetshift_[4]+p[0]+p[1]*eshift_[10]+p[2]*eshift_[11];//be-D3
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeBFGH(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1];//ab-L
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+(p[2]+1)*eshift_[3];//ac-H
    case 2: return esetshift_[1]+(p[0]+1)+p[1]*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[2]+p[0]+p[1]*eshift_[6]+(p[2]+1)*eshift_[7];//bc-D1
    case 4: return esetshift_[3]+(p[0]+1)+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronEdgeBDGH(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]+(p[1]+1)*eshift_[0]+(p[2]+1)*eshift_[1];//ab-L
    case 1: return esetshift_[0]+(p[0]+1)+p[1]*eshift_[2]+p[2]*eshift_[3];//ac-H
    case 2: return esetshift_[1]+(p[0]+1)+(p[1]+1)*eshift_[4]+p[2]*eshift_[5];//ae-P
    case 3: return esetshift_[3]+(p[0]+1)+p[1]*eshift_[8]+p[2]*eshift_[9];//ag-D2
    case 4: return esetshift_[4]+p[0]+(p[1]+1)*eshift_[10]+p[2]*eshift_[11];//be-D3
    case 5: return esetshift_[5]+p[0]+p[1]*eshift_[12]+p[2]*eshift_[13];//bg-D4
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleABCG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1];
    case 1: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5];
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleBCDG(const int p[3], const int id) const{
  switch(id){
    case 0: return p[0]*2+p[1]*tshift_[0]+p[2]*tshift_[1]+1;
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11];
    case 2: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 3: return tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3];
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleABEG(const int p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3];
    case 1: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9];
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+1;
    case 3: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleBEFG(const int p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7]+1;
    case 1: return tsetshift_[0]+p[0]*2+p[1]*tshift_[2]+p[2]*tshift_[3]+1;
    case 2: return p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1];
    case 3: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleBFGH(const int p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 1: return tsetshift_[4]+p[0]*2+p[1]*tshift_[10]+p[2]*tshift_[11]+1;
    case 2: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+3;
    case 3: return p[0]*2+p[1]*tshift_[0]+(p[2]+1)*tshift_[1]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronTriangleBDGH(const int p[3], const int id) const{
  switch(id){
    case 0: return tsetshift_[1]+p[0]*2+p[1]*tshift_[4]+p[2]*tshift_[5]+2;
    case 1: return tsetshift_[2]+p[0]*2+p[1]*tshift_[6]+p[2]*tshift_[7];
    case 2: return tsetshift_[3]+p[0]*2+p[1]*tshift_[8]+p[2]*tshift_[9]+1;
    case 3: return tsetshift_[0]+p[0]*2+(p[1]+1)*tshift_[2]+p[2]*tshift_[3]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborABCG(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t+1;
    case 1: return t+2;
    case 2:
            if(p[0]>0) return t-1;
            else return t-tetshift_[1]+3;
    case 3: return t-tetshift_[1]+3;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborBCDG(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+4;
    case 2:
            if(p[2]>0) return t-tetshift_[1]+3;
            else return t+tetshift_[0]+1;
    case 3: return t+tetshift_[0]+1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborABEG(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t-2;
    case 1: return t+1;
    case 2:
            if(p[0]>0) return t-4;
            else return t-tetshift_[0]-1;
    case 3: return t-tetshift_[0]-1;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborBEFG(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+1;
    case 2:
            if(p[1]>0) return t-tetshift_[0]+2;
            else return t+tetshift_[1]-3;
    case 3: return t+tetshift_[1]-3;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborBFGH(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t+1;
    case 2:
            if(p[0]<nbvoxels_[0]-1) return t+4;
            else return t+tetshift_[1]-3;
    case 3: return t+tetshift_[1]-3;
  }
  return -1;
}

inline int ImplicitTriangulation::getTetrahedronNeighborBDGH(const int t, const int p[3], const int id) const{
  switch(id){
    case 0: return t-1;
    case 1: return t-4;
    case 2:
            if(p[0]<nbvoxels_[0]-1) return t+1;
            else return t+tetshift_[0]-2;
    case 3: return t+tetshift_[0]-2;
  }
  return -1;
}

#endif // _IMPLICITTRIANGULATION_H
