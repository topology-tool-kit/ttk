/// \ingroup base
/// \class ttk::ImplicitTriangulation
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date January 2016.
///
/// \brief ImplicitTriangulation is a class that provides time and memory
/// efficient traversal methods on triangulations of piecewise linear
/// manifolds represented by regular grids.
/// \sa Triangulation

#pragma once

#include <array>

// base code includes
#include <RegularGridTriangulation.h>

#define CASE_EDGE_POSITION_L_3D \
  case EdgePosition::L_xnn_3D:  \
  case EdgePosition::L_xn0_3D:  \
  case EdgePosition::L_xnN_3D:  \
  case EdgePosition::L_x0n_3D:  \
  case EdgePosition::L_x00_3D:  \
  case EdgePosition::L_x0N_3D:  \
  case EdgePosition::L_xNn_3D:  \
  case EdgePosition::L_xN0_3D:  \
  case EdgePosition::L_xNN_3D
#define CASE_EDGE_POSITION_H_3D \
  case EdgePosition::H_nyn_3D:  \
  case EdgePosition::H_ny0_3D:  \
  case EdgePosition::H_nyN_3D:  \
  case EdgePosition::H_0yn_3D:  \
  case EdgePosition::H_0y0_3D:  \
  case EdgePosition::H_0yN_3D:  \
  case EdgePosition::H_Nyn_3D:  \
  case EdgePosition::H_Ny0_3D:  \
  case EdgePosition::H_NyN_3D
#define CASE_EDGE_POSITION_P_3D \
  case EdgePosition::P_nnz_3D:  \
  case EdgePosition::P_n0z_3D:  \
  case EdgePosition::P_nNz_3D:  \
  case EdgePosition::P_0nz_3D:  \
  case EdgePosition::P_00z_3D:  \
  case EdgePosition::P_0Nz_3D:  \
  case EdgePosition::P_Nnz_3D:  \
  case EdgePosition::P_N0z_3D:  \
  case EdgePosition::P_NNz_3D
#define CASE_EDGE_POSITION_D1_3D \
  case EdgePosition::D1_xyn_3D:  \
  case EdgePosition::D1_xy0_3D:  \
  case EdgePosition::D1_xyN_3D
#define CASE_EDGE_POSITION_D2_3D \
  case EdgePosition::D2_nyz_3D:  \
  case EdgePosition::D2_0yz_3D:  \
  case EdgePosition::D2_Nyz_3D
#define CASE_EDGE_POSITION_D3_3D \
  case EdgePosition::D3_xnz_3D:  \
  case EdgePosition::D3_x0z_3D:  \
  case EdgePosition::D3_xNz_3D
#define CASE_EDGE_POSITION_L_2D \
  case EdgePosition::L_xn_2D:   \
  case EdgePosition::L_x0_2D:   \
  case EdgePosition::L_xN_2D
#define CASE_EDGE_POSITION_H_2D \
  case EdgePosition::H_ny_2D:   \
  case EdgePosition::H_0y_2D:   \
  case EdgePosition::H_Ny_2D

namespace ttk {

  class ImplicitTriangulation : public RegularGridTriangulation {

  public:
    ImplicitTriangulation();
    ~ImplicitTriangulation() override;

    ImplicitTriangulation(const ImplicitTriangulation &) = default;
    ImplicitTriangulation(ImplicitTriangulation &&) = default;
    ImplicitTriangulation &operator=(const ImplicitTriangulation &) = default;
    ImplicitTriangulation &operator=(ImplicitTriangulation &&) = default;

    inline const std::array<SimplexId, 3> &getGridDimensions() const final {
      return this->dimensions_;
    }

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const final {

      if(dimensionality_ == 3)
        getTetrahedronEdge(cellId, localEdgeId, edgeId);
      else if(dimensionality_ == 2)
        getTriangleEdgeInternal(cellId, localEdgeId, edgeId);
      else if(dimensionality_ == 1)
        getCellNeighbor(cellId, localEdgeId, edgeId);

      return 0;
    }

    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId & /*cellId*/) const final {

      if(dimensionality_ == 3)
        return 6;
      else if(dimensionality_ == 2)
        return 3;

      return 0;
    }

    const std::vector<std::vector<SimplexId>> *getCellEdgesInternal() final;

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const final {

      if(dimensionality_ == 3)
        getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
      else if(dimensionality_ == 2)
        getTriangleNeighbor(cellId, localNeighborId, neighborId);
      else if(dimensionality_ == 1) {
        printErr("getCellNeighbor() not implemented in 1D! (TODO)");
        return -1;
      }

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const final {

      if(dimensionality_ == 3)
        return getTetrahedronNeighborNumber(cellId);
      else if(dimensionality_ == 2)
        return getTriangleNeighborNumber(cellId);
      else if(dimensionality_ == 1) {
        printErr("getCellNeighborNumber() not implemented in 1D! (TODO)");
        return -1;
      }

      return 0;
    }

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() final;

    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const final {

      if(dimensionality_ == 3)
        getTetrahedronTriangle(cellId, localTriangleId, triangleId);

      return 0;
    }

    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId & /*cellId*/) const final {
      // NOTE: the output is always 4 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 4;
    }

    const std::vector<std::vector<SimplexId>> *getCellTrianglesInternal() final;

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const final {

      if(dimensionality_ == 3)
        getTetrahedronVertex(cellId, localVertexId, vertexId);
      else if(dimensionality_ == 2)
        getTriangleVertexInternal(cellId, localVertexId, vertexId);
      else if(dimensionality_ == 1)
        getEdgeVertexInternal(cellId, localVertexId, vertexId);

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId & /*cellId*/) const final {

      return dimensionality_ + 1;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const final {
      return dimensionality_;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const final {

      return TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(edgeId);
    }

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() final;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() final;

    const std::vector<std::vector<SimplexId>> *getEdgeTrianglesInternal() final;

    const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() final;

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const final {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdgesInternal() const final {
      return edgeNumber_;
    }

    inline SimplexId getNumberOfTrianglesInternal() const final {
      return triangleNumber_;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const final {
      return vertexNumber_;
    }

    virtual int getTetrahedronEdge(const SimplexId &tetId,
                                   const int &id,
                                   SimplexId &edgeId) const = 0;

    inline int
      getTetrahedronEdges(std::vector<std::vector<SimplexId>> &edges) const {
      edges.resize(tetrahedronNumber_);
      for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
        edges[i].resize(6);
        for(int j = 0; j < 6; ++j)
          getTetrahedronEdge(i, j, edges[i][j]);
      }

      return 0;
    }

    virtual int getTetrahedronTriangle(const SimplexId &tetId,
                                       const int &id,
                                       SimplexId &triangleId) const = 0;

    inline int getTetrahedronTriangles(
      std::vector<std::vector<SimplexId>> &triangles) const {

      triangles.resize(tetrahedronNumber_);
      for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
        triangles[i].resize(4);
        for(int j = 0; j < 4; ++j)
          getTetrahedronTriangle(i, j, triangles[i][j]);
      }

      return 0;
    }

    virtual int getTetrahedronNeighbor(const SimplexId &tetId,
                                       const int &localNeighborId,
                                       SimplexId &neighborId) const = 0;

    virtual SimplexId
      getTetrahedronNeighborNumber(const SimplexId &tetId) const = 0;

    inline int
      getTetrahedronNeighbors(std::vector<std::vector<SimplexId>> &neighbors) {

      neighbors.resize(tetrahedronNumber_);
      for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
        neighbors[i].resize(getTetrahedronNeighborNumber(i));
        for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
          getTetrahedronNeighbor(i, j, neighbors[i][j]);
      }

      return 0;
    }

    virtual int getTetrahedronVertex(const SimplexId &tetId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const = 0;

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId & /*triangleId*/) const final {
      // NOTE: the output is always 3 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 3;
    }

    const std::vector<std::vector<SimplexId>> *getTriangleEdgesInternal() final;

    inline int getTriangleEdgesInternal(
      std::vector<std::vector<SimplexId>> &edges) const {

      edges.resize(triangleNumber_);
      for(SimplexId i = 0; i < triangleNumber_; ++i) {
        edges[i].resize(3);
        for(int j = 0; j < 3; ++j)
          getTriangleEdgeInternal(i, j, edges[i][j]);
      }
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const final {

      return TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(triangleId);
    }

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() final;

    virtual int getTriangleNeighbor(const SimplexId &triangleId,
                                    const int &localNeighborId,
                                    SimplexId &neighborId) const = 0;

    virtual SimplexId
      getTriangleNeighborNumber(const SimplexId &triangleId) const = 0;

    inline int
      getTriangleNeighbors(std::vector<std::vector<SimplexId>> &neighbors) {

      neighbors.resize(triangleNumber_);
      for(SimplexId i = 0; i < triangleNumber_; ++i) {
        neighbors[i].resize(getTriangleNeighborNumber(i));
        for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
          getTriangleNeighbor(i, j, neighbors[i][j]);
      }
      return 0;
    }

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() final;

    const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() final;

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const final {
      return TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(vertexId);
    }

    const std::vector<std::vector<SimplexId>> *getVertexEdgesInternal() final;

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const final {
      return TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(vertexId);
    }

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() final;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() final;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() final;

    const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() final;

    inline bool isEmpty() const final {
      return !vertexNumber_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const final {

#ifdef TTK_ENABLE_MPI
      if(this->metaGrid_ != nullptr) {
        return this->isTriangleOnGlobalBoundaryInternal(triangleId);
      }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 or triangleId >= triangleNumber_)
        return false;
#endif // !TTK_ENABLE_KAMIKAZE

      if(dimensionality_ == 3)
        return (TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(triangleId)
                == 1);

      return false;
    }

    int setInputGrid(const float &xOrigin,
                     const float &yOrigin,
                     const float &zOrigin,
                     const float &xSpacing,
                     const float &ySpacing,
                     const float &zSpacing,
                     const SimplexId &xDim,
                     const SimplexId &yDim,
                     const SimplexId &zDim) final;

    virtual int preconditionVerticesInternal() = 0;
    int preconditionVertexNeighborsInternal() final;
    int preconditionEdgesInternal() override = 0;
    int preconditionTrianglesInternal() override = 0;
    virtual int preconditionTetrahedronsInternal() = 0;

    inline int preconditionCellsInternal() {
      if(dimensionality_ == 3) {
        return this->preconditionTetrahedronsInternal();
      } else if(dimensionality_ == 2 && !hasPreconditionedTriangles_) {
        hasPreconditionedTriangles_ = true;
        return this->preconditionTrianglesInternal();
      }
      return 0;
    }

    inline int preconditionVerticesAndCells() {
      if(!this->hasPreconditionedVerticesAndCells_) {
        this->preconditionVerticesInternal();
        this->preconditionCellsInternal();
        this->hasPreconditionedVerticesAndCells_ = true;
      }
      return 0;
    }

    inline int getCellVTKIDInternal(const int &ttkId, int &vtkId) const final {
#ifndef TTK_ENABLE_KAMIKAZE
      if(ttkId < 0) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      const SimplexId nSimplexPerCell{this->getDimensionality() == 3 ? 6 : 2};
      vtkId = ttkId / nSimplexPerCell;
      return 0;
    }

#ifdef TTK_ENABLE_MPI

  protected:
    int preconditionDistributedCells() final;

  public:
    void createMetaGrid(const double *const bounds) final;
    inline int getCellRankInternal(const SimplexId lcid) const final {
      const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
      const auto locCubeId{lcid / nTetraPerCube};

      if(this->cellGhost_[locCubeId] == 0) {
        return ttk::MPIrank_;
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(this->neighborRanks_.empty() && ttk::MPIsize_ > 1) {
        this->printErr("Empty neighborsRanks_!");
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE

      const auto nVertsCell{this->getCellVertexNumber(lcid)};
      std::vector<bool> inRank(nVertsCell);
      for(const auto neigh : this->neighborRanks_) {
        std::fill(inRank.begin(), inRank.end(), false);
        const auto &bbox{this->neighborCellBBoxes_[neigh]};
        for(SimplexId i = 0; i < nVertsCell; ++i) {
          SimplexId v{};
          this->getCellVertex(lcid, i, v);
          if(this->vertexGhost_[v] == 0) {
            inRank[i] = true;
          } else {
            const auto p{this->getVertGlobalCoords(v)};
            if(p[0] >= bbox[0] && p[0] <= bbox[1] && p[1] >= bbox[2]
               && p[1] <= bbox[3] && p[2] >= bbox[4] && p[2] <= bbox[5]) {
              inRank[i] = true;
            }
          }
        }
        if(std::all_of(
             inRank.begin(), inRank.end(), [](const bool v) { return v; })) {
          return neigh;
        }
      }

      return -1;
    }

  protected:
    inline bool
      isVertexOnGlobalBoundaryInternal(const SimplexId lvid) const final {

      if(!ttk::isRunningWithMPI()) {
        // NOTE: perf killer
        return this->isVertexOnBoundary(lvid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(lvid > this->TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() - 1
         || lvid < 0) {
        return false;
      }
      if(this->metaGrid_ == nullptr) {
        return false;
      }
#endif // TTK_ENABLE_KAMIKAZE

      const auto gvid{this->getVertexGlobalIdInternal(lvid)};
      if(gvid == -1) {
        return false;
      }
      return this->metaGrid_->isVertexOnBoundary(gvid);
    }

    inline bool
      isEdgeOnGlobalBoundaryInternal(const SimplexId leid) const final {

      if(!ttk::isRunningWithMPI()) {
        // NOTE: perf killer
        return this->isEdgeOnBoundary(leid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(leid > this->getNumberOfEdgesInternal() - 1 || leid < 0) {
        return false;
      }
      if(this->metaGrid_ == nullptr) {
        return false;
      }
#endif // TTK_ENABLE_KAMIKAZE

      const auto geid{this->getEdgeGlobalIdInternal(leid)};
      if(geid == -1) {
        return false;
      }
      return this->metaGrid_->isEdgeOnBoundary(geid);
    }

    inline bool
      isTriangleOnGlobalBoundaryInternal(const SimplexId ltid) const final {

      if(!ttk::isRunningWithMPI()) {
        // NOTE: perf killer
        return this->isTriangleOnBoundary(ltid);
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if(ltid > this->getNumberOfTrianglesInternal() - 1 || ltid < 0) {
        return false;
      }
      if(this->metaGrid_ == nullptr) {
        return false;
      }
#endif // TTK_ENABLE_KAMIKAZE

      const auto gtid{this->getTriangleGlobalIdInternal(ltid)};
      if(gtid == -1) {
        return false;
      }
      return this->metaGrid_->isTriangleOnBoundary(gtid);
    }

  private:
    inline std::array<SimplexId, 3>
      getVertGlobalCoords(const SimplexId lvid) const final {

      // local vertex coordinates
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 3) {
        this->vertexToPosition(lvid, p.data());
      } else if(this->dimensionality_ == 2) {
        this->vertexToPosition2d(lvid, p.data());
      }

      // global vertex coordinates
      p[0] += this->localGridOffset_[0];
      p[1] += this->localGridOffset_[1];
      p[2] += this->localGridOffset_[2];

      return p;
    }

    inline std::array<SimplexId, 3>
      getVertLocalCoords(const SimplexId gvid) const final {

      // global vertex coordinates
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 3) {
        this->metaGrid_->vertexToPosition(gvid, p.data());
      } else if(this->dimensionality_ == 2) {
        this->metaGrid_->vertexToPosition2d(gvid, p.data());
      }

      // local vertex coordinates
      p[0] -= this->localGridOffset_[0];
      p[1] -= this->localGridOffset_[1];
      p[2] -= this->localGridOffset_[2];

      return p;
    }

#endif // TTK_ENABLE_MPI

  protected:
    enum class VertexPosition : char {
      // a--------b

      LEFT_CORNER_1D, // a
      RIGHT_CORNER_1D, // b
      CENTER_1D,
      // total: 3 1D cases

      // a--------b
      // |        |
      // |        |
      // |        |
      // c--------d

      // 2D corners
      TOP_LEFT_CORNER_2D, // a
      TOP_RIGHT_CORNER_2D, // b
      BOTTOM_LEFT_CORNER_2D, // c
      BOTTOM_RIGHT_CORNER_2D, // d
      // 2D edges
      TOP_EDGE_2D, // ab
      BOTTOM_EDGE_2D, // cd
      LEFT_EDGE_2D, // ac
      RIGHT_EDGE_2D, // bd
      // 2D central strip
      CENTER_2D,
      // total: 9 2D cases

      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // 3D corners
      TOP_LEFT_FRONT_CORNER_3D, // a
      TOP_RIGHT_FRONT_CORNER_3D, // b
      BOTTOM_LEFT_FRONT_CORNER_3D, // c
      BOTTOM_RIGHT_FRONT_CORNER_3D, // d
      TOP_LEFT_BACK_CORNER_3D, // e
      TOP_RIGHT_BACK_CORNER_3D, // f
      BOTTOM_LEFT_BACK_CORNER_3D, // g
      BOTTOM_RIGHT_BACK_CORNER_3D, // h
      // 3D edges
      TOP_FRONT_EDGE_3D, // ab
      BOTTOM_FRONT_EDGE_3D, // cd
      LEFT_FRONT_EDGE_3D, // ac
      RIGHT_FRONT_EDGE_3D, // bd
      TOP_BACK_EDGE_3D, // ef
      BOTTOM_BACK_EDGE_3D, // gh
      LEFT_BACK_EDGE_3D, // eg
      RIGHT_BACK_EDGE_3D, // fh
      TOP_LEFT_EDGE_3D, // ae
      TOP_RIGHT_EDGE_3D, // bf
      BOTTOM_LEFT_EDGE_3D, // cg
      BOTTOM_RIGHT_EDGE_3D, // dh
      // 3D faces
      FRONT_FACE_3D, // abcd
      BACK_FACE_3D, // efgh
      TOP_FACE_3D, // abef
      BOTTOM_FACE_3D, // cdgh
      LEFT_FACE_3D, // aceg
      RIGHT_FACE_3D, // bdfh
      // 3D central part
      CENTER_3D,
      // total: 27 3D cases
    };

    // vertex neighbor shifts
    std::array<SimplexId, 14> vertexNeighborABCDEFGH_{};

    std::array<SimplexId, 10> vertexNeighborABCD_{};
    std::array<SimplexId, 10> vertexNeighborEFGH_{};
    std::array<SimplexId, 10> vertexNeighborAEFB_{};
    std::array<SimplexId, 10> vertexNeighborGHDC_{};
    std::array<SimplexId, 10> vertexNeighborAEGC_{};
    std::array<SimplexId, 10> vertexNeighborBFHD_{};

    std::array<SimplexId, 8> vertexNeighborAB_{};
    std::array<SimplexId, 8> vertexNeighborBD_{};
    std::array<SimplexId, 8> vertexNeighborGH_{};
    std::array<SimplexId, 8> vertexNeighborEG_{};
    std::array<SimplexId, 8> vertexNeighborCG_{};
    std::array<SimplexId, 8> vertexNeighborBF_{};

    std::array<SimplexId, 7> vertexNeighborB_{};
    std::array<SimplexId, 7> vertexNeighborG_{};

    std::array<SimplexId, 6> vertexNeighborEF_{};
    std::array<SimplexId, 6> vertexNeighborCD_{};
    std::array<SimplexId, 6> vertexNeighborAC_{};
    std::array<SimplexId, 6> vertexNeighborAE_{};
    std::array<SimplexId, 6> vertexNeighborFH_{};
    std::array<SimplexId, 6> vertexNeighborDH_{};

    std::array<SimplexId, 4> vertexNeighborA_{};
    std::array<SimplexId, 4> vertexNeighborC_{};
    std::array<SimplexId, 4> vertexNeighborD_{};
    std::array<SimplexId, 4> vertexNeighborE_{};
    std::array<SimplexId, 4> vertexNeighborF_{};
    std::array<SimplexId, 4> vertexNeighborH_{};

    std::array<SimplexId, 6> vertexNeighbor2dABCD_{};
    std::array<SimplexId, 4> vertexNeighbor2dAB_{};
    std::array<SimplexId, 4> vertexNeighbor2dCD_{};
    std::array<SimplexId, 4> vertexNeighbor2dAC_{};
    std::array<SimplexId, 4> vertexNeighbor2dBD_{};
    std::array<SimplexId, 3> vertexNeighbor2dB_{};
    std::array<SimplexId, 3> vertexNeighbor2dC_{};
    std::array<SimplexId, 2> vertexNeighbor2dA_{};
    std::array<SimplexId, 2> vertexNeighbor2dD_{};

    enum class EdgePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // length (ab)
      L_xnn_3D,
      L_xn0_3D,
      L_xnN_3D,
      L_x0n_3D,
      L_x00_3D,
      L_x0N_3D,
      L_xNn_3D,
      L_xN0_3D,
      L_xNN_3D,
      // height (ac)
      H_nyn_3D,
      H_ny0_3D,
      H_nyN_3D,
      H_0yn_3D,
      H_0y0_3D,
      H_0yN_3D,
      H_Nyn_3D,
      H_Ny0_3D,
      H_NyN_3D,
      // depth (ae)
      P_nnz_3D,
      P_n0z_3D,
      P_nNz_3D,
      P_0nz_3D,
      P_00z_3D,
      P_0Nz_3D,
      P_Nnz_3D,
      P_N0z_3D,
      P_NNz_3D,
      // diagonal1 (bc)
      D1_xyn_3D,
      D1_xy0_3D,
      D1_xyN_3D,
      // diagonal2 (ag)
      D2_nyz_3D,
      D2_0yz_3D,
      D2_Nyz_3D,
      // diagonal3 (be)
      D3_xnz_3D,
      D3_x0z_3D,
      D3_xNz_3D,
      // diagonal4 (bg)
      D4_3D,

      // length (ab)
      L_xn_2D,
      L_x0_2D,
      L_xN_2D,
      // height (ac)
      H_ny_2D,
      H_0y_2D,
      H_Ny_2D,
      // diagonal1 (bc)
      D1_2D,

      FIRST_EDGE_1D,
      LAST_EDGE_1D,
      CENTER_1D,
    };

    enum class TrianglePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      F_3D, // face (abc, bcd)
      C_3D, // side (abe, bef)
      H_3D, // top (acg, aeg)
      D1_3D, // diagonal1 (bdg, beg)
      D2_3D, // diagonal2 (abg, bgh)
      D3_3D, // diagonal3 (bcg, bfg)

      TOP_2D, // abc
      BOTTOM_2D, // bcd
    };

    bool hasPreconditionedVerticesAndCells_{false};

    float origin_[3]; //
    float spacing_[3]; //
    SimplexId nbvoxels_[3]; // nombre de voxels par axe

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
    inline bool isPowerOfTwo(unsigned long long int v,
                             unsigned long long int &r) const {
      if(v && !(v & (v - 1))) {
        r = 0;
        while(v >>= 1)
          r++;
        return true;
      }
      return false;
    }

    //\cond
    // 2D //
    void vertexToPosition2d(const SimplexId vertex, SimplexId p[2]) const final;
    void
      edgeToPosition2d(const SimplexId edge, const int k, SimplexId p[2]) const;
    void triangleToPosition2d(const SimplexId triangle,
                              SimplexId p[2]) const final;

    inline SimplexId getVertexEdge2dA(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexEdge2dB(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexEdge2dC(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexEdge2dD(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexEdge2dAB(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexEdge2dCD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexEdge2dAC(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexEdge2dBD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexEdge2dABCD(const SimplexId p[2],
                                         const int id) const;

    inline SimplexId getVertexStar2dA(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexStar2dB(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexStar2dC(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexStar2dD(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexStar2dAB(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexStar2dCD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexStar2dAC(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexStar2dBD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexStar2dABCD(const SimplexId p[2],
                                         const int id) const;

    inline SimplexId getVertexLink2dA(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexLink2dB(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexLink2dC(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexLink2dD(const SimplexId p[2], const int id) const;
    inline SimplexId getVertexLink2dAB(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexLink2dCD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexLink2dAC(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexLink2dBD(const SimplexId p[2],
                                       const int id) const;
    inline SimplexId getVertexLink2dABCD(const SimplexId p[2],
                                         const int id) const;

    inline SimplexId getEdgeTriangleL_x0(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleL_xn(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleL_xN(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleH_0y(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleH_ny(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleH_Ny(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getEdgeTriangleD1_xy(const SimplexId p[3],
                                          const int id) const;

    inline SimplexId getEdgeLink2dL(const SimplexId p[2], const int id) const;
    inline SimplexId getEdgeLink2dH(const SimplexId p[2], const int id) const;
    inline SimplexId getEdgeLink2dD1(const SimplexId p[2], const int id) const;

    inline SimplexId getEdgeStar2dL(const SimplexId p[2], const int id) const;
    inline SimplexId getEdgeStar2dH(const SimplexId p[2], const int id) const;

    // 3D //
    void vertexToPosition(const SimplexId vertex, SimplexId p[3]) const final;
    void
      edgeToPosition(const SimplexId edge, const int k, SimplexId p[3]) const;
    void triangleToPosition(const SimplexId triangle,
                            const int k,
                            SimplexId p[3]) const final;
    void tetrahedronToPosition(const SimplexId tetrahedron,
                               SimplexId p[3]) const final;

    inline SimplexId getVertexEdgeA(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeAB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeCD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeEF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeGH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeAC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeBD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeEG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeFH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeAE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeBF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeCG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeDH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexEdgeABDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeEFHG(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeAEGC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeBFHD(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeAEFB(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeGHDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexEdgeABCDEFGH(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getVertexTriangleA(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleB(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleC(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleD(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleE(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleF(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleG(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleH(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getVertexTriangleAB(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleCD(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleEF(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleGH(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleAC(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleBD(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleEG(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleFH(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleAE(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleBF(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleCG(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleDH(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getVertexTriangleABDC(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleEFHG(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleAEGC(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleBFHD(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleAEFB(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleGHDC(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getVertexTriangleABCDEFGH(const SimplexId p[3],
                                               const int id) const;

    inline SimplexId getVertexLinkA(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkAB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkCD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkEF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkGH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkAC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkBD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkEG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkFH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkAE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkBF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkCG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkDH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexLinkABDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkEFHG(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkAEGC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkBFHD(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkAEFB(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkGHDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexLinkABCDEFGH(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getVertexStarA(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarAB(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarCD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarEF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarGH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarAC(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarBD(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarEG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarFH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarAE(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarBF(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarCG(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarDH(const SimplexId p[3], const int id) const;
    inline SimplexId getVertexStarABDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarEFHG(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarAEGC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarBFHD(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarAEFB(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarGHDC(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getVertexStarABCDEFGH(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getEdgeTriangleL_x00(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_x0n(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_x0N(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xn0(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xnn(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xnN(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xN0(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xNn(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleL_xNN(const SimplexId p[3],
                                          const int id) const;

    inline SimplexId getEdgeTriangleH_0y0(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_0yn(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_0yN(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_ny0(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_nyn(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_nyN(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_Ny0(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_Nyn(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleH_NyN(const SimplexId p[3],
                                          const int id) const;

    inline SimplexId getEdgeTriangleP_00z(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_0nz(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_0Nz(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_n0z(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_nnz(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_nNz(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_N0z(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_Nnz(const SimplexId p[3],
                                          const int id) const;
    inline SimplexId getEdgeTriangleP_NNz(const SimplexId p[3],
                                          const int id) const;

    inline SimplexId getEdgeTriangleD1_xy0(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD1_xyn(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD1_xyN(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getEdgeTriangleD2_0yz(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD2_nyz(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD2_Nyz(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getEdgeTriangleD3_x0z(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD3_xnz(const SimplexId p[3],
                                           const int id) const;
    inline SimplexId getEdgeTriangleD3_xNz(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getEdgeTriangleD4_xyz(const SimplexId p[3],
                                           const int id) const;

    inline SimplexId getEdgeLinkL(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkH(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkP(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkD1(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkD2(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkD3(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeLinkD4(const SimplexId p[3], const int id) const;

    inline SimplexId getEdgeStarL(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeStarH(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeStarP(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeStarD1(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeStarD2(const SimplexId p[3], const int id) const;
    inline SimplexId getEdgeStarD3(const SimplexId p[3], const int id) const;

    inline SimplexId getTriangleVertexF(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleVertexH(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleVertexC(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleVertexD1(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleVertexD2(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleVertexD3(const SimplexId p[3],
                                         const int id) const;

    inline SimplexId getTriangleEdgeF_0(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeF_1(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeH_0(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeH_1(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeC_0(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeC_1(const SimplexId p[3],
                                        const int id) const;
    inline SimplexId getTriangleEdgeD1_0(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleEdgeD1_1(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleEdgeD2_0(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleEdgeD2_1(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleEdgeD3_0(const SimplexId p[3],
                                         const int id) const;
    inline SimplexId getTriangleEdgeD3_1(const SimplexId p[3],
                                         const int id) const;

    inline SimplexId getTriangleLinkF(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleLinkH(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleLinkC(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleLinkD1(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getTriangleLinkD2(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getTriangleLinkD3(const SimplexId p[3],
                                       const int id) const;

    inline SimplexId getTriangleStarF(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleStarH(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleStarC(const SimplexId p[3], const int id) const;
    inline SimplexId getTriangleStarD1(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getTriangleStarD2(const SimplexId p[3],
                                       const int id) const;
    inline SimplexId getTriangleStarD3(const SimplexId p[3],
                                       const int id) const;

    inline SimplexId getTetrahedronVertexABCG(const SimplexId p[3],
                                              const int id) const;
    inline SimplexId getTetrahedronVertexBCDG(const SimplexId p[3],
                                              const int id) const;
    inline SimplexId getTetrahedronVertexABEG(const SimplexId p[3],
                                              const int id) const;
    inline SimplexId getTetrahedronVertexBEFG(const SimplexId p[3],
                                              const int id) const;
    inline SimplexId getTetrahedronVertexBFGH(const SimplexId p[3],
                                              const int id) const;
    inline SimplexId getTetrahedronVertexBDGH(const SimplexId p[3],
                                              const int id) const;

    inline SimplexId getTetrahedronEdgeABCG(const SimplexId p[3],
                                            const int id) const;
    inline SimplexId getTetrahedronEdgeBCDG(const SimplexId p[3],
                                            const int id) const;
    inline SimplexId getTetrahedronEdgeABEG(const SimplexId p[3],
                                            const int id) const;
    inline SimplexId getTetrahedronEdgeBEFG(const SimplexId p[3],
                                            const int id) const;
    inline SimplexId getTetrahedronEdgeBFGH(const SimplexId p[3],
                                            const int id) const;
    inline SimplexId getTetrahedronEdgeBDGH(const SimplexId p[3],
                                            const int id) const;

    inline SimplexId getTetrahedronTriangleABCG(const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronTriangleBCDG(const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronTriangleABEG(const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronTriangleBEFG(const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronTriangleBFGH(const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronTriangleBDGH(const SimplexId p[3],
                                                const int id) const;

    inline SimplexId getTetrahedronNeighborABCG(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronNeighborBCDG(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronNeighborABEG(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronNeighborBEFG(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronNeighborBFGH(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    inline SimplexId getTetrahedronNeighborBDGH(const SimplexId t,
                                                const SimplexId p[3],
                                                const int id) const;
    //\endcond
  };

  template <typename Derived>
  class ImplicitTriangulationCRTP : public ImplicitTriangulation {
    inline Derived &underlying() {
      return static_cast<Derived &>(*this);
    }
    inline Derived const &underlying() const {
      return static_cast<Derived const &>(*this);
    }

  public:
    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexId < 0 or vertexId >= vertexNumber_)
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          return 14;
        case VertexPosition::FRONT_FACE_3D:
        case VertexPosition::BACK_FACE_3D:
        case VertexPosition::TOP_FACE_3D:
        case VertexPosition::BOTTOM_FACE_3D:
        case VertexPosition::LEFT_FACE_3D:
        case VertexPosition::RIGHT_FACE_3D:
          return 10;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          return 8;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          return 7;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
        case VertexPosition::CENTER_2D:
          return 6;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
        case VertexPosition::TOP_EDGE_2D:
        case VertexPosition::BOTTOM_EDGE_2D:
        case VertexPosition::LEFT_EDGE_2D:
        case VertexPosition::RIGHT_EDGE_2D:
          return 4;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          return 3;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
        case VertexPosition::CENTER_1D:
          return 2;
        case VertexPosition::LEFT_CORNER_1D:
        case VertexPosition::RIGHT_CORNER_1D:
          return 1;
      }

      return -1;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const final {
#ifdef TTK_ENABLE_MPI
      if(this->metaGrid_ != nullptr) {
        return this->isVertexOnGlobalBoundaryInternal(vertexId);
      }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexId < 0 or vertexId >= vertexNumber_)
        return false;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
        case VertexPosition::CENTER_2D:
        case VertexPosition::CENTER_1D:
          return false;
        default:
          return true;
      }
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const final {

#ifdef TTK_ENABLE_MPI
      if(this->metaGrid_ != nullptr) {
        return this->isEdgeOnGlobalBoundaryInternal(edgeId);
      }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
      if(edgeId < 0 or edgeId >= edgeNumber_)
        return false;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getEdgePosition(edgeId)) {
        case EdgePosition::L_xnn_3D:
        case EdgePosition::H_nyn_3D:
        case EdgePosition::P_nnz_3D:
        case EdgePosition::D1_xyn_3D:
        case EdgePosition::D2_nyz_3D:
        case EdgePosition::D3_xnz_3D:
        case EdgePosition::D4_3D:
        case EdgePosition::L_xn_2D:
        case EdgePosition::H_ny_2D:
        case EdgePosition::D1_2D:
          return false;
        default:
          break;
      }
      return true;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId < 0
         or localNeighborId >= getVertexNeighborNumber(vertexId))
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          neighborId
            = vertexId + this->vertexNeighborABCDEFGH_[localNeighborId];
          break;
        case VertexPosition::FRONT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborABCD_[localNeighborId];
          break;
        case VertexPosition::BACK_FACE_3D:
          neighborId = vertexId + this->vertexNeighborEFGH_[localNeighborId];
          break;
        case VertexPosition::TOP_FACE_3D:
          neighborId = vertexId + this->vertexNeighborAEFB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          neighborId = vertexId + this->vertexNeighborGHDC_[localNeighborId];
          break;
        case VertexPosition::LEFT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborAEGC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborBFHD_[localNeighborId];
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          neighborId = vertexId + this->vertexNeighborAB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          neighborId = vertexId + this->vertexNeighborCD_[localNeighborId];
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          neighborId = vertexId + this->vertexNeighborAC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          neighborId = vertexId + this->vertexNeighborBD_[localNeighborId];
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          neighborId = vertexId + this->vertexNeighborEF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          neighborId = vertexId + this->vertexNeighborGH_[localNeighborId];
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          neighborId = vertexId + this->vertexNeighborEG_[localNeighborId];
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          neighborId = vertexId + this->vertexNeighborFH_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          neighborId = vertexId + this->vertexNeighborAE_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          neighborId = vertexId + this->vertexNeighborBF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          neighborId = vertexId + this->vertexNeighborCG_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          neighborId = vertexId + this->vertexNeighborDH_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          neighborId = vertexId + this->vertexNeighborA_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          neighborId = vertexId + this->vertexNeighborB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          neighborId = vertexId + this->vertexNeighborC_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          neighborId = vertexId + this->vertexNeighborD_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          neighborId = vertexId + this->vertexNeighborE_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          neighborId = vertexId + this->vertexNeighborF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          neighborId = vertexId + this->vertexNeighborG_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          neighborId = vertexId + this->vertexNeighborH_[localNeighborId];
          break;
        case VertexPosition::CENTER_2D:
          neighborId = vertexId + this->vertexNeighbor2dABCD_[localNeighborId];
          break;
        case VertexPosition::TOP_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dAB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dCD_[localNeighborId];
          break;
        case VertexPosition::LEFT_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dAC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dBD_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
          neighborId = vertexId + this->vertexNeighbor2dA_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
          neighborId = vertexId + this->vertexNeighbor2dB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          neighborId = vertexId + this->vertexNeighbor2dC_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          neighborId = vertexId + this->vertexNeighbor2dD_[localNeighborId];
          break;
        case VertexPosition::CENTER_1D:
          neighborId = (localNeighborId == 0 ? vertexId + 1 : vertexId - 1);
          break;
        case VertexPosition::LEFT_CORNER_1D:
          neighborId = vertexId + 1;
          break;
        case VertexPosition::RIGHT_CORNER_1D:
          neighborId = vertexId - 1;
          break;
        default:
          neighborId = -1;
          break;
      }

      return 0;
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localEdgeId < 0
         or localEdgeId >= getVertexEdgeNumberInternal(vertexId))
        return -1;
#endif
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--g-----b--h
      // | /      | /
      // |/       |/
      // c--------d
      //
      // Classement des "Edges" et dans cet ordre:
      // L: largeur (type ab)
      // H: hauteur (type ac)
      // P: profondeur (type ae)
      // D1: diagonale1 (type bc)
      // D2: diagonale2 (type ag)
      // D3: diagonale3 (type be)
      // D4: diagonale4 (type bg)

      const auto &p = this->underlying().getVertexCoords(vertexId);

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          edgeId = getVertexEdgeABCDEFGH(p.data(), localEdgeId);
          break;
        case VertexPosition::FRONT_FACE_3D:
          edgeId = getVertexEdgeABDC(p.data(), localEdgeId);
          break;
        case VertexPosition::BACK_FACE_3D:
          edgeId = getVertexEdgeEFHG(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_FACE_3D:
          edgeId = getVertexEdgeAEFB(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          edgeId = getVertexEdgeGHDC(p.data(), localEdgeId);
          break;
        case VertexPosition::LEFT_FACE_3D:
          edgeId = getVertexEdgeAEGC(p.data(), localEdgeId);
          break;
        case VertexPosition::RIGHT_FACE_3D:
          edgeId = getVertexEdgeBFHD(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          edgeId = getVertexEdgeAB(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          edgeId = getVertexEdgeCD(p.data(), localEdgeId);
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          edgeId = getVertexEdgeAC(p.data(), localEdgeId);
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          edgeId = getVertexEdgeBD(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          edgeId = getVertexEdgeEF(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          edgeId = getVertexEdgeGH(p.data(), localEdgeId);
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          edgeId = getVertexEdgeEG(p.data(), localEdgeId);
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          edgeId = getVertexEdgeFH(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          edgeId = getVertexEdgeAE(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          edgeId = getVertexEdgeBF(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          edgeId = getVertexEdgeCG(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          edgeId = getVertexEdgeDH(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          edgeId = getVertexEdgeA(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          edgeId = getVertexEdgeB(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          edgeId = getVertexEdgeC(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          edgeId = getVertexEdgeD(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          edgeId = getVertexEdgeE(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          edgeId = getVertexEdgeF(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          edgeId = getVertexEdgeG(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          edgeId = getVertexEdgeH(p.data(), localEdgeId);
          break;
        case VertexPosition::CENTER_2D:
          edgeId = getVertexEdge2dABCD(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_EDGE_2D:
          edgeId = getVertexEdge2dAB(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_EDGE_2D:
          edgeId = getVertexEdge2dCD(p.data(), localEdgeId);
          break;
        case VertexPosition::LEFT_EDGE_2D:
          edgeId = getVertexEdge2dAC(p.data(), localEdgeId);
          break;
        case VertexPosition::RIGHT_EDGE_2D:
          edgeId = getVertexEdge2dBD(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
          edgeId = getVertexEdge2dA(p.data(), localEdgeId);
          break;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
          edgeId = getVertexEdge2dB(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          edgeId = getVertexEdge2dC(p.data(), localEdgeId);
          break;
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          edgeId = getVertexEdge2dD(p.data(), localEdgeId);
          break;
        case VertexPosition::CENTER_1D:
          edgeId = (localEdgeId == 0 ? vertexId : vertexId - 1);
          break;
        case VertexPosition::LEFT_CORNER_1D:
          edgeId = vertexId;
          break;
        case VertexPosition::RIGHT_CORNER_1D:
          edgeId = vertexId - 1;
          break;
        default:
          edgeId = -1;
      }

      return 0;
    }

    inline SimplexId
      getVertexTriangleNumberInternal(const SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexId < 0 or vertexId >= vertexNumber_)
        return -1;
#endif

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          return 36;
        case VertexPosition::FRONT_FACE_3D:
        case VertexPosition::BACK_FACE_3D:
        case VertexPosition::TOP_FACE_3D:
        case VertexPosition::BOTTOM_FACE_3D:
        case VertexPosition::LEFT_FACE_3D:
        case VertexPosition::RIGHT_FACE_3D:
          return 21;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          return 15;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          return 12;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          return 9;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          return 5;
        default: // 1D + 2D
          break;
      }

      return 0;
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const final {
#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId < 0
         or localTriangleId >= getVertexTriangleNumberInternal(vertexId))
        return -1;
#endif

      const auto &p = this->underlying().getVertexCoords(vertexId);

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          triangleId = getVertexTriangleABCDEFGH(p.data(), localTriangleId);
          break;
        case VertexPosition::FRONT_FACE_3D:
          triangleId = getVertexTriangleABDC(p.data(), localTriangleId);
          break;
        case VertexPosition::BACK_FACE_3D:
          triangleId = getVertexTriangleEFHG(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_FACE_3D:
          triangleId = getVertexTriangleAEFB(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          triangleId = getVertexTriangleGHDC(p.data(), localTriangleId);
          break;
        case VertexPosition::LEFT_FACE_3D:
          triangleId = getVertexTriangleAEGC(p.data(), localTriangleId);
          break;
        case VertexPosition::RIGHT_FACE_3D:
          triangleId = getVertexTriangleBFHD(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          triangleId = getVertexTriangleAB(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          triangleId = getVertexTriangleCD(p.data(), localTriangleId);
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          triangleId = getVertexTriangleAC(p.data(), localTriangleId);
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          triangleId = getVertexTriangleBD(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          triangleId = getVertexTriangleEF(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          triangleId = getVertexTriangleGH(p.data(), localTriangleId);
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          triangleId = getVertexTriangleEG(p.data(), localTriangleId);
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          triangleId = getVertexTriangleFH(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          triangleId = getVertexTriangleAE(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          triangleId = getVertexTriangleBF(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          triangleId = getVertexTriangleCG(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          triangleId = getVertexTriangleDH(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          triangleId = getVertexTriangleA(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          triangleId = getVertexTriangleB(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          triangleId = getVertexTriangleC(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          triangleId = getVertexTriangleD(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          triangleId = getVertexTriangleE(p.data(), localTriangleId);
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          triangleId = getVertexTriangleF(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          triangleId = getVertexTriangleG(p.data(), localTriangleId);
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          triangleId = getVertexTriangleH(p.data(), localTriangleId);
          break;
        default: // 1D + 2D
          triangleId = -1;
          break;
      }

      return 0;
    }

    inline int
      TTK_TRIANGULATION_INTERNAL(getVertexLink)(const SimplexId &vertexId,
                                                const int &localLinkId,
                                                SimplexId &linkId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      const auto &p = this->underlying().getVertexCoords(vertexId);

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          linkId = getVertexLinkABCDEFGH(p.data(), localLinkId);
          break;
        case VertexPosition::FRONT_FACE_3D:
          linkId = getVertexLinkABDC(p.data(), localLinkId);
          break;
        case VertexPosition::BACK_FACE_3D:
          linkId = getVertexLinkEFHG(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_FACE_3D:
          linkId = getVertexLinkAEFB(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          linkId = getVertexLinkGHDC(p.data(), localLinkId);
          break;
        case VertexPosition::LEFT_FACE_3D:
          linkId = getVertexLinkAEGC(p.data(), localLinkId);
          break;
        case VertexPosition::RIGHT_FACE_3D:
          linkId = getVertexLinkBFHD(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          linkId = getVertexLinkAB(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          linkId = getVertexLinkCD(p.data(), localLinkId);
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          linkId = getVertexLinkAC(p.data(), localLinkId);
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          linkId = getVertexLinkBD(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          linkId = getVertexLinkEF(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          linkId = getVertexLinkGH(p.data(), localLinkId);
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          linkId = getVertexLinkEG(p.data(), localLinkId);
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          linkId = getVertexLinkFH(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          linkId = getVertexLinkAE(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          linkId = getVertexLinkBF(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          linkId = getVertexLinkCG(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          linkId = getVertexLinkDH(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          linkId = getVertexLinkA(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          linkId = getVertexLinkB(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          linkId = getVertexLinkC(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          linkId = getVertexLinkD(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          linkId = getVertexLinkE(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          linkId = getVertexLinkF(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          linkId = getVertexLinkG(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          linkId = getVertexLinkH(p.data(), localLinkId);
          break;
        case VertexPosition::CENTER_2D:
          linkId = getVertexLink2dABCD(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_EDGE_2D:
          linkId = getVertexLink2dAB(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_EDGE_2D:
          linkId = getVertexLink2dCD(p.data(), localLinkId);
          break;
        case VertexPosition::LEFT_EDGE_2D:
          linkId = getVertexLink2dAC(p.data(), localLinkId);
          break;
        case VertexPosition::RIGHT_EDGE_2D:
          linkId = getVertexLink2dBD(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
          linkId = getVertexLink2dA(p.data(), localLinkId);
          break;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
          linkId = getVertexLink2dB(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          linkId = getVertexLink2dC(p.data(), localLinkId);
          break;
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          linkId = getVertexLink2dD(p.data(), localLinkId);
          break;
        default: // 1D
          linkId = -1;
          break;
      };

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(vertexId < 0 or vertexId >= vertexNumber_)
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          return 24;
        case VertexPosition::FRONT_FACE_3D:
        case VertexPosition::BACK_FACE_3D:
        case VertexPosition::TOP_FACE_3D:
        case VertexPosition::BOTTOM_FACE_3D:
        case VertexPosition::LEFT_FACE_3D:
        case VertexPosition::RIGHT_FACE_3D:
          return 12;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          return 8;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
        case VertexPosition::CENTER_2D:
          return 6;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          return 4;
        case VertexPosition::TOP_EDGE_2D: // ab
        case VertexPosition::BOTTOM_EDGE_2D: // cd
        case VertexPosition::LEFT_EDGE_2D: // ac
        case VertexPosition::RIGHT_EDGE_2D: // bd
          return 3;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          return 2;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          return 1;
        default: // 1D
          break;
      }

      return 0;
    }

    inline int
      TTK_TRIANGULATION_INTERNAL(getVertexStar)(const SimplexId &vertexId,
                                                const int &localStarId,
                                                SimplexId &starId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      const auto &p = this->underlying().getVertexCoords(vertexId);

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          starId = getVertexStarABCDEFGH(p.data(), localStarId);
          break;
        case VertexPosition::FRONT_FACE_3D:
          starId = getVertexStarABDC(p.data(), localStarId);
          break;
        case VertexPosition::BACK_FACE_3D:
          starId = getVertexStarEFHG(p.data(), localStarId);
          break;
        case VertexPosition::TOP_FACE_3D:
          starId = getVertexStarAEFB(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          starId = getVertexStarGHDC(p.data(), localStarId);
          break;
        case VertexPosition::LEFT_FACE_3D:
          starId = getVertexStarAEGC(p.data(), localStarId);
          break;
        case VertexPosition::RIGHT_FACE_3D:
          starId = getVertexStarBFHD(p.data(), localStarId);
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          starId = getVertexStarAB(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          starId = getVertexStarCD(p.data(), localStarId);
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          starId = getVertexStarAC(p.data(), localStarId);
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          starId = getVertexStarBD(p.data(), localStarId);
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          starId = getVertexStarEF(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          starId = getVertexStarGH(p.data(), localStarId);
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          starId = getVertexStarEG(p.data(), localStarId);
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          starId = getVertexStarFH(p.data(), localStarId);
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          starId = getVertexStarAE(p.data(), localStarId);
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          starId = getVertexStarBF(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          starId = getVertexStarCG(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          starId = getVertexStarDH(p.data(), localStarId);
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          starId = getVertexStarA(p.data(), localStarId);
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          starId = getVertexStarB(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          starId = getVertexStarC(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          starId = getVertexStarD(p.data(), localStarId);
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          starId = getVertexStarE(p.data(), localStarId);
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          starId = getVertexStarF(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          starId = getVertexStarG(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          starId = getVertexStarH(p.data(), localStarId);
          break;
        case VertexPosition::CENTER_2D:
          starId = getVertexStar2dABCD(p.data(), localStarId);
          break;
        case VertexPosition::TOP_EDGE_2D:
          starId = getVertexStar2dAB(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_EDGE_2D:
          starId = getVertexStar2dCD(p.data(), localStarId);
          break;
        case VertexPosition::LEFT_EDGE_2D:
          starId = getVertexStar2dAC(p.data(), localStarId);
          break;
        case VertexPosition::RIGHT_EDGE_2D:
          starId = getVertexStar2dBD(p.data(), localStarId);
          break;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
          starId = getVertexStar2dA(p.data(), localStarId);
          break;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
          starId = getVertexStar2dB(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          starId = getVertexStar2dC(p.data(), localStarId);
          break;
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          starId = getVertexStar2dD(p.data(), localStarId);
          break;
        default: // 1D
          starId = -1;
          break;
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
      const SimplexId &vertexId, float &x, float &y, float &z) const final {
      if(dimensionality_ == 3) {
        const auto &p = this->underlying().getVertexCoords(vertexId);

        x = origin_[0] + spacing_[0] * p[0];
        y = origin_[1] + spacing_[1] * p[1];
        z = origin_[2] + spacing_[2] * p[2];
      } else if(dimensionality_ == 2) {
        const auto &p = this->underlying().getVertexCoords(vertexId);

        if(dimensions_[0] > 1 and dimensions_[1] > 1) {
          x = origin_[0] + spacing_[0] * p[0];
          y = origin_[1] + spacing_[1] * p[1];
          z = origin_[2];
        } else if(dimensions_[1] > 1 and dimensions_[2] > 1) {
          x = origin_[0];
          y = origin_[1] + spacing_[1] * p[0];
          z = origin_[2] + spacing_[2] * p[1];
        } else if(dimensions_[0] > 1 and dimensions_[2] > 1) {
          x = origin_[0] + spacing_[0] * p[0];
          y = origin_[1];
          z = origin_[2] + spacing_[2] * p[1];
        }
      } else if(dimensionality_ == 1) {
        if(dimensions_[0] > 1) {
          x = origin_[0] + spacing_[0] * vertexId;
          y = origin_[1];
          z = origin_[2];
        } else if(dimensions_[1] > 1) {
          x = origin_[0];
          y = origin_[1] + spacing_[1] * vertexId;
          z = origin_[2];
        } else if(dimensions_[2] > 1) {
          x = origin_[0];
          y = origin_[1];
          z = origin_[2] + spacing_[2] * vertexId;
        }
      }

      return 0;
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(edgeId < 0 or edgeId >= edgeNumber_)
        return -1;
      if(localVertexId < 0 or localVertexId >= 2)
        return -2;
#endif

      const auto &p = this->underlying().getEdgeCoords(edgeId);

      const auto helper3d
        = [&](const SimplexId a, const SimplexId b) -> SimplexId {
        if(isAccelerated_) {
          const auto tmp = p[0] + (p[1] << div_[0]) + (p[2] << div_[1]);
          return (localVertexId == 0) ? tmp + a : tmp + b;
        } else {
          const auto tmp = p[0] + (p[1] * vshift_[0]) + (p[2] * vshift_[1]);
          return (localVertexId == 0) ? tmp + a : tmp + b;
        }
      };

      const auto helper2d
        = [&](const SimplexId a, const SimplexId b) -> SimplexId {
        if(isAccelerated_) {
          const auto tmp = p[0] + (p[1] << div_[0]);
          return localVertexId == 0 ? tmp + a : tmp + b;
        } else {
          const auto tmp = p[0] + (p[1] * vshift_[0]);
          return localVertexId == 0 ? tmp + a : tmp + b;
        }
      };

      switch(this->underlying().getEdgePosition(edgeId)) {
      CASE_EDGE_POSITION_L_3D:
        vertexId = helper3d(0, 1);
        break;
      CASE_EDGE_POSITION_H_3D:
        vertexId = helper3d(0, vshift_[0]);
        break;
      CASE_EDGE_POSITION_P_3D:
        vertexId = helper3d(0, vshift_[1]);
        break;
      CASE_EDGE_POSITION_D1_3D:
        vertexId = helper3d(1, vshift_[0]);
        break;
      CASE_EDGE_POSITION_D2_3D:
        vertexId = helper3d(0, vshift_[0] + vshift_[1]);
        break;
      CASE_EDGE_POSITION_D3_3D:
        vertexId = helper3d(1, vshift_[1]);
        break;
        case EdgePosition::D4_3D:
          vertexId = helper3d(1, vshift_[0] + vshift_[1]);
          break;
        CASE_EDGE_POSITION_L_2D:
          vertexId = helper2d(0, 1);
          break;
        CASE_EDGE_POSITION_H_2D:
          vertexId = helper2d(0, vshift_[0]);
          break;
        case EdgePosition::D1_2D:
          vertexId = helper2d(1, vshift_[0]);
          break;

        case EdgePosition::FIRST_EDGE_1D:
          vertexId = localVertexId == 0 ? 0 : 1;
          break;
        case EdgePosition::LAST_EDGE_1D:
          vertexId = localVertexId == 0 ? edgeNumber_ - 1 : edgeNumber_;
          break;
        case EdgePosition::CENTER_1D:
          vertexId = localVertexId == 0 ? edgeId : edgeId + 1;
          break;
      }

      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(edgeId < 0 or edgeId >= edgeNumber_)
        return -1;
#endif

      switch(this->underlying().getEdgePosition(edgeId)) {
        case EdgePosition::L_xnn_3D:
        case EdgePosition::H_nyn_3D:
        case EdgePosition::P_nnz_3D:
        case EdgePosition::D4_3D:
          return 6;
        case EdgePosition::L_x0n_3D:
        case EdgePosition::L_xNn_3D:
        case EdgePosition::L_xn0_3D:
        case EdgePosition::L_xnN_3D:
        case EdgePosition::H_ny0_3D:
        case EdgePosition::H_nyN_3D:
        case EdgePosition::H_0yn_3D:
        case EdgePosition::H_Nyn_3D:
        case EdgePosition::P_n0z_3D:
        case EdgePosition::P_nNz_3D:
        case EdgePosition::P_0nz_3D:
        case EdgePosition::P_Nnz_3D:
        case EdgePosition::D1_xyn_3D:
        case EdgePosition::D2_nyz_3D:
        case EdgePosition::D3_xnz_3D:
          return 4;
        case EdgePosition::L_x00_3D:
        case EdgePosition::L_xNN_3D:
        case EdgePosition::H_0yN_3D:
        case EdgePosition::H_Ny0_3D:
        case EdgePosition::P_0Nz_3D:
        case EdgePosition::P_N0z_3D:
        case EdgePosition::D1_xy0_3D:
        case EdgePosition::D1_xyN_3D:
        case EdgePosition::D2_0yz_3D:
        case EdgePosition::D2_Nyz_3D:
        case EdgePosition::D3_x0z_3D:
        case EdgePosition::D3_xNz_3D:
          return 3;
        case EdgePosition::L_xN0_3D:
        case EdgePosition::L_x0N_3D:
        case EdgePosition::H_0y0_3D:
        case EdgePosition::H_NyN_3D:
        case EdgePosition::P_00z_3D:
        case EdgePosition::P_NNz_3D:
        case EdgePosition::L_xn_2D:
        case EdgePosition::H_ny_2D:
        case EdgePosition::D1_2D:
          return 2;
        case EdgePosition::L_x0_2D:
        case EdgePosition::L_xN_2D:
        case EdgePosition::H_0y_2D:
        case EdgePosition::H_Ny_2D:
          return 1;

        default: // 1D
          break;
      }

      return 0;
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const final {
#ifndef TTK_ENABLE_KAMIKAZE
      if(localTriangleId < 0
         or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
        return -1;
#endif

      const auto &p = this->underlying().getEdgeCoords(edgeId);

      switch(this->underlying().getEdgePosition(edgeId)) {
        case EdgePosition::L_xnn_3D:
          triangleId = getEdgeTriangleL_xnn(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xn0_3D:
          triangleId = getEdgeTriangleL_xn0(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xnN_3D:
          triangleId = getEdgeTriangleL_xnN(p.data(), localTriangleId);
          break;
        case EdgePosition::L_x0n_3D:
          triangleId = getEdgeTriangleL_x0n(p.data(), localTriangleId);
          break;
        case EdgePosition::L_x00_3D:
          triangleId = getEdgeTriangleL_x00(p.data(), localTriangleId);
          break;
        case EdgePosition::L_x0N_3D:
          triangleId = getEdgeTriangleL_x0N(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xNn_3D:
          triangleId = getEdgeTriangleL_xNn(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xN0_3D:
          triangleId = getEdgeTriangleL_xN0(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xNN_3D:
          triangleId = getEdgeTriangleL_xNN(p.data(), localTriangleId);
          break;
        case EdgePosition::H_nyn_3D:
          triangleId = getEdgeTriangleH_nyn(p.data(), localTriangleId);
          break;
        case EdgePosition::H_ny0_3D:
          triangleId = getEdgeTriangleH_ny0(p.data(), localTriangleId);
          break;
        case EdgePosition::H_nyN_3D:
          triangleId = getEdgeTriangleH_nyN(p.data(), localTriangleId);
          break;
        case EdgePosition::H_0yn_3D:
          triangleId = getEdgeTriangleH_0yn(p.data(), localTriangleId);
          break;
        case EdgePosition::H_0y0_3D:
          triangleId = getEdgeTriangleH_0y0(p.data(), localTriangleId);
          break;
        case EdgePosition::H_0yN_3D:
          triangleId = getEdgeTriangleH_0yN(p.data(), localTriangleId);
          break;
        case EdgePosition::H_Nyn_3D:
          triangleId = getEdgeTriangleH_Nyn(p.data(), localTriangleId);
          break;
        case EdgePosition::H_Ny0_3D:
          triangleId = getEdgeTriangleH_Ny0(p.data(), localTriangleId);
          break;
        case EdgePosition::H_NyN_3D:
          triangleId = getEdgeTriangleH_NyN(p.data(), localTriangleId);
          break;
        case EdgePosition::P_nnz_3D:
          triangleId = getEdgeTriangleP_nnz(p.data(), localTriangleId);
          break;
        case EdgePosition::P_n0z_3D:
          triangleId = getEdgeTriangleP_n0z(p.data(), localTriangleId);
          break;
        case EdgePosition::P_nNz_3D:
          triangleId = getEdgeTriangleP_nNz(p.data(), localTriangleId);
          break;
        case EdgePosition::P_0nz_3D:
          triangleId = getEdgeTriangleP_0nz(p.data(), localTriangleId);
          break;
        case EdgePosition::P_00z_3D:
          triangleId = getEdgeTriangleP_00z(p.data(), localTriangleId);
          break;
        case EdgePosition::P_0Nz_3D:
          triangleId = getEdgeTriangleP_0Nz(p.data(), localTriangleId);
          break;
        case EdgePosition::P_Nnz_3D:
          triangleId = getEdgeTriangleP_Nnz(p.data(), localTriangleId);
          break;
        case EdgePosition::P_N0z_3D:
          triangleId = getEdgeTriangleP_N0z(p.data(), localTriangleId);
          break;
        case EdgePosition::P_NNz_3D:
          triangleId = getEdgeTriangleP_NNz(p.data(), localTriangleId);
          break;
        case EdgePosition::D1_xyn_3D:
          triangleId = getEdgeTriangleD1_xyn(p.data(), localTriangleId);
          break;
        case EdgePosition::D1_xy0_3D:
          triangleId = getEdgeTriangleD1_xy0(p.data(), localTriangleId);
          break;
        case EdgePosition::D1_xyN_3D:
          triangleId = getEdgeTriangleD1_xyN(p.data(), localTriangleId);
          break;
        case EdgePosition::D2_nyz_3D:
          triangleId = getEdgeTriangleD2_nyz(p.data(), localTriangleId);
          break;
        case EdgePosition::D2_0yz_3D:
          triangleId = getEdgeTriangleD2_0yz(p.data(), localTriangleId);
          break;
        case EdgePosition::D2_Nyz_3D:
          triangleId = getEdgeTriangleD2_Nyz(p.data(), localTriangleId);
          break;
        case EdgePosition::D3_xnz_3D:
          triangleId = getEdgeTriangleD3_xnz(p.data(), localTriangleId);
          break;
        case EdgePosition::D3_x0z_3D:
          triangleId = getEdgeTriangleD3_x0z(p.data(), localTriangleId);
          break;
        case EdgePosition::D3_xNz_3D:
          triangleId = getEdgeTriangleD3_xNz(p.data(), localTriangleId);
          break;
        case EdgePosition::D4_3D:
          triangleId = getEdgeTriangleD4_xyz(p.data(), localTriangleId);
          break;

        case EdgePosition::L_xn_2D:
          triangleId = getEdgeTriangleL_xn(p.data(), localTriangleId);
          break;
        case EdgePosition::L_x0_2D:
          triangleId = getEdgeTriangleL_x0(p.data(), localTriangleId);
          break;
        case EdgePosition::L_xN_2D:
          triangleId = getEdgeTriangleL_xN(p.data(), localTriangleId);
          break;
        case EdgePosition::H_ny_2D:
          triangleId = getEdgeTriangleH_ny(p.data(), localTriangleId);
          break;
        case EdgePosition::H_0y_2D:
          triangleId = getEdgeTriangleH_0y(p.data(), localTriangleId);
          break;
        case EdgePosition::H_Ny_2D:
          triangleId = getEdgeTriangleH_Ny(p.data(), localTriangleId);
          break;
        case EdgePosition::D1_2D:
          triangleId = getEdgeTriangleD1_xy(p.data(), localTriangleId);
          break;

        default: // 1D
          triangleId = -1;
          break;
      }

      return 0;
    }

    inline int
      TTK_TRIANGULATION_INTERNAL(getEdgeLink)(const SimplexId &edgeId,
                                              const int &localLinkId,
                                              SimplexId &linkId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
        return -1;
#endif

      const auto &p = this->underlying().getEdgeCoords(edgeId);

      switch(this->underlying().getEdgePosition(edgeId)) {
      CASE_EDGE_POSITION_L_3D:
        linkId = getEdgeLinkL(p.data(), localLinkId);
        break;
      CASE_EDGE_POSITION_H_3D:
        linkId = getEdgeLinkH(p.data(), localLinkId);
        break;
      CASE_EDGE_POSITION_P_3D:
        linkId = getEdgeLinkP(p.data(), localLinkId);
        break;
      CASE_EDGE_POSITION_D1_3D:
        linkId = getEdgeLinkD1(p.data(), localLinkId);
        break;
      CASE_EDGE_POSITION_D2_3D:
        linkId = getEdgeLinkD2(p.data(), localLinkId);
        break;
      CASE_EDGE_POSITION_D3_3D:
        linkId = getEdgeLinkD3(p.data(), localLinkId);
        break;
        case EdgePosition::D4_3D:
          linkId = getEdgeLinkD4(p.data(), localLinkId);
          break;

        CASE_EDGE_POSITION_L_2D:
          linkId = getEdgeLink2dL(p.data(), localLinkId);
          break;
        CASE_EDGE_POSITION_H_2D:
          linkId = getEdgeLink2dH(p.data(), localLinkId);
          break;
        case EdgePosition::D1_2D:
          linkId = getEdgeLink2dD1(p.data(), localLinkId);
          break;

        default: // 1D
          linkId = -1;
          break;
      }

      return 0;
    }

    inline int
      TTK_TRIANGULATION_INTERNAL(getEdgeStar)(const SimplexId &edgeId,
                                              const int &localStarId,
                                              SimplexId &starId) const final {
#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
        return -1;
#endif

      const auto &p = this->underlying().getEdgeCoords(edgeId);

      switch(this->underlying().getEdgePosition(edgeId)) {
      CASE_EDGE_POSITION_L_3D:
        starId = getEdgeStarL(p.data(), localStarId);
        break;
      CASE_EDGE_POSITION_H_3D:
        starId = getEdgeStarH(p.data(), localStarId);
        break;
      CASE_EDGE_POSITION_P_3D:
        starId = getEdgeStarP(p.data(), localStarId);
        break;
      CASE_EDGE_POSITION_D1_3D:
        starId = getEdgeStarD1(p.data(), localStarId);
        break;
      CASE_EDGE_POSITION_D2_3D:
        starId = getEdgeStarD2(p.data(), localStarId);
        break;
      CASE_EDGE_POSITION_D3_3D:
        starId = getEdgeStarD3(p.data(), localStarId);
        break;
        case EdgePosition::D4_3D:
          starId = p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
                   + localStarId;
          break;

        CASE_EDGE_POSITION_L_2D:
          starId = getEdgeStar2dL(p.data(), localStarId);
          break;
        CASE_EDGE_POSITION_H_2D:
          starId = getEdgeStar2dH(p.data(), localStarId);
          break;
        case EdgePosition::D1_2D:
          starId = p[0] * 2 + p[1] * tshift_[0] + localStarId;
          break;

        default: // 1D
          starId = -1;
          break;
      }

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(edgeId < 0 or edgeId >= edgeNumber_)
        return -1;
#endif

      switch(this->underlying().getEdgePosition(edgeId)) {
        case EdgePosition::L_xnn_3D:
        case EdgePosition::H_nyn_3D:
        case EdgePosition::P_nnz_3D:
        case EdgePosition::D4_3D:
          return 6;
        case EdgePosition::D1_xyn_3D:
        case EdgePosition::D2_nyz_3D:
        case EdgePosition::D3_xnz_3D:
          return 4;
        case EdgePosition::L_x0n_3D:
        case EdgePosition::L_xNn_3D:
        case EdgePosition::L_xn0_3D:
        case EdgePosition::L_xnN_3D:
        case EdgePosition::H_ny0_3D:
        case EdgePosition::H_nyN_3D:
        case EdgePosition::H_0yn_3D:
        case EdgePosition::H_Nyn_3D:
        case EdgePosition::P_n0z_3D:
        case EdgePosition::P_nNz_3D:
        case EdgePosition::P_0nz_3D:
        case EdgePosition::P_Nnz_3D:
          return 3;
        case EdgePosition::L_x00_3D:
        case EdgePosition::L_xNN_3D:
        case EdgePosition::H_0yN_3D:
        case EdgePosition::H_Ny0_3D:
        case EdgePosition::P_0Nz_3D:
        case EdgePosition::P_N0z_3D:
        case EdgePosition::D1_xy0_3D:
        case EdgePosition::D1_xyN_3D:
        case EdgePosition::D2_0yz_3D:
        case EdgePosition::D2_Nyz_3D:
        case EdgePosition::D3_x0z_3D:
        case EdgePosition::D3_xNz_3D:
        case EdgePosition::L_xn_2D:
        case EdgePosition::H_ny_2D:
        case EdgePosition::D1_2D:
          return 2;
        case EdgePosition::L_xN0_3D:
        case EdgePosition::L_x0N_3D:
        case EdgePosition::H_0y0_3D:
        case EdgePosition::H_NyN_3D:
        case EdgePosition::P_00z_3D:
        case EdgePosition::P_NNz_3D:
        case EdgePosition::L_x0_2D:
        case EdgePosition::L_xN_2D:
        case EdgePosition::H_0y_2D:
        case EdgePosition::H_Ny_2D:
          return 1;

        default: // 1D
          break;
      }

      return 0;
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 or triangleId >= triangleNumber_)
        return -1;
      if(localVertexId < 0 or localVertexId >= 3)
        return -2;
#endif

      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--g-----b--h
      // | /      | /
      // |/       |/
      // c--------d
      //
      // Classement des "Triangles" et dans cet ordre:
      // F: face (type abc/bcd)
      // C: cote (type abe/bef)
      // H: haut (type acg/aeg)
      // D1: diagonale1 (type bdg/beg)
      // D2: diagonale2 (type abg/bgh)
      // D3: diagonale3 (type bcg/bfg)

      const auto &p = this->underlying().getTriangleCoords(triangleId);
      vertexId = -1;

      switch(this->underlying().getTrianglePosition(triangleId)) {
        case TrianglePosition::F_3D:
          vertexId = getTriangleVertexF(p.data(), localVertexId);
          break;
        case TrianglePosition::H_3D:
          vertexId = getTriangleVertexH(p.data(), localVertexId);
          break;
        case TrianglePosition::C_3D:
          vertexId = getTriangleVertexC(p.data(), localVertexId);
          break;
        case TrianglePosition::D1_3D:
          vertexId = getTriangleVertexD1(p.data(), localVertexId);
          break;
        case TrianglePosition::D2_3D:
          vertexId = getTriangleVertexD2(p.data(), localVertexId);
          break;
        case TrianglePosition::D3_3D:
          vertexId = getTriangleVertexD3(p.data(), localVertexId);
          break;
        case TrianglePosition::TOP_2D:
          switch(localVertexId) {
            break;
            case 0:
              vertexId = p[0] / 2 + p[1] * vshift_[0];
              break;
            case 1:
              vertexId = p[0] / 2 + p[1] * vshift_[0] + 1;
              break;
            case 2:
              vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0];
              break;
          }
          break;
        case TrianglePosition::BOTTOM_2D:
          switch(localVertexId) {
            break;
            case 0:
              vertexId = p[0] / 2 + p[1] * vshift_[0] + 1;
              break;
            case 1:
              vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0] + 1;
              break;
            case 2:
              vertexId = p[0] / 2 + p[1] * vshift_[0] + vshift_[0];
              break;
          }
      }

      return 0;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 or triangleId >= triangleNumber_)
        return -1;
      if(localEdgeId < 0 or localEdgeId >= 3)
        return -2;
#endif

      const auto &p = this->underlying().getTriangleCoords(triangleId);
      const auto par = triangleId % 2;
      edgeId = -1;

      switch(this->underlying().getTrianglePosition(triangleId)) {
        case TrianglePosition::F_3D:
          edgeId = (par == 1) ? getTriangleEdgeF_1(p.data(), localEdgeId)
                              : getTriangleEdgeF_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::H_3D:
          edgeId = (par == 1) ? getTriangleEdgeH_1(p.data(), localEdgeId)
                              : getTriangleEdgeH_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::C_3D:
          edgeId = (par == 1) ? getTriangleEdgeC_1(p.data(), localEdgeId)
                              : getTriangleEdgeC_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::D1_3D:
          edgeId = (par == 1) ? getTriangleEdgeD1_1(p.data(), localEdgeId)
                              : getTriangleEdgeD1_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::D2_3D:
          edgeId = (par == 1) ? getTriangleEdgeD2_1(p.data(), localEdgeId)
                              : getTriangleEdgeD2_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::D3_3D:
          edgeId = (par == 1) ? getTriangleEdgeD3_1(p.data(), localEdgeId)
                              : getTriangleEdgeD3_0(p.data(), localEdgeId);
          break;
        case TrianglePosition::TOP_2D:
          switch(localEdgeId) {
            break;
            case 0:
              edgeId = p[0] / 2 + p[1] * eshift_[0];
              break;
            case 1:
              edgeId = esetshift_[0] + p[0] / 2 + p[1] * eshift_[2];
              break;
            case 2:
              edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
              break;
          }
          break;
        case TrianglePosition::BOTTOM_2D:
          switch(localEdgeId) {
            break;
            case 0:
              edgeId = p[0] / 2 + (p[1] + 1) * eshift_[0];
              break;
            case 1:
              edgeId = esetshift_[0] + (p[0] + 1) / 2 + p[1] * eshift_[2];
              break;
            case 2:
              edgeId = esetshift_[1] + p[0] / 2 + p[1] * eshift_[4];
              break;
          }
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
        return -1;
#endif

      const auto &p = this->underlying().getTriangleCoords(triangleId);

      switch(this->underlying().getTrianglePosition(triangleId)) {
        case TrianglePosition::F_3D:
          linkId = getTriangleLinkF(p.data(), localLinkId);
          break;
        case TrianglePosition::H_3D:
          linkId = getTriangleLinkH(p.data(), localLinkId);
          break;
        case TrianglePosition::C_3D:
          linkId = getTriangleLinkC(p.data(), localLinkId);
          break;
        case TrianglePosition::D1_3D:
          linkId = getTriangleLinkD1(p.data(), localLinkId);
          break;
        case TrianglePosition::D2_3D:
          linkId = getTriangleLinkD2(p.data(), localLinkId);
          break;
        case TrianglePosition::D3_3D:
          linkId = getTriangleLinkD3(p.data(), localLinkId);
          break;
        default: // 2D
          linkId = -1;
          break;
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
        return -1;
#endif

      const auto &p = this->underlying().getTriangleCoords(triangleId);

      switch(this->underlying().getTrianglePosition(triangleId)) {
        case TrianglePosition::F_3D:
          starId = getTriangleStarF(p.data(), localStarId);
          break;
        case TrianglePosition::H_3D:
          starId = getTriangleStarH(p.data(), localStarId);
          break;
        case TrianglePosition::C_3D:
          starId = getTriangleStarC(p.data(), localStarId);
          break;
        case TrianglePosition::D1_3D:
          starId = getTriangleStarD1(p.data(), localStarId);
          break;
        case TrianglePosition::D2_3D:
          starId = getTriangleStarD2(p.data(), localStarId);
          break;
        case TrianglePosition::D3_3D:
          starId = getTriangleStarD3(p.data(), localStarId);
          break;
        default: // 2D
          starId = -1;
          break;
      }

      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 or triangleId >= triangleNumber_)
        return -1;
#endif

      const auto &p = this->underlying().getTriangleCoords(triangleId);

      switch(this->underlying().getTrianglePosition(triangleId)) {
        case TrianglePosition::F_3D:
          return (p[2] > 0 and p[2] < nbvoxels_[2]) ? 2 : 1;
        case TrianglePosition::H_3D:
          return (p[1] > 0 and p[1] < nbvoxels_[1]) ? 2 : 1;
        case TrianglePosition::C_3D:
          return (p[0] < 2 or p[0] >= (dimensions_[0] * 2 - 2)) ? 1 : 2;

        case TrianglePosition::D1_3D:
        case TrianglePosition::D2_3D:
        case TrianglePosition::D3_3D:
          return 2;
        default: // 2D
          break;
      }
      return 0;
    }

    inline int getTriangleNeighbor(const SimplexId &triangleId,
                                   const int &localNeighborId,
                                   SimplexId &neighborId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId < 0
         or localNeighborId >= getTriangleNeighborNumber(triangleId))
        return -1;
#endif

      neighborId = -1;

      if(dimensionality_ == 2) {
        const auto &p = this->underlying().getTriangleCoords(triangleId);
        const SimplexId id = triangleId % 2;

        if(id) {
          if(p[0] / 2 == nbvoxels_[Di_] - 1 and p[1] == nbvoxels_[Dj_] - 1)
            neighborId = triangleId - 1;
          else if(p[0] / 2 == nbvoxels_[Di_] - 1) {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId - 1;
                break;
              case 1:
                neighborId = triangleId + tshift_[0] - 1;
                break;
            }
          } else if(p[1] == nbvoxels_[Dj_] - 1) {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId - 1;
                break;
              case 1:
                neighborId = triangleId + 1;
                break;
            }
          } else {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId - 1;
                break;
              case 1:
                neighborId = triangleId + 1;
                break;
              case 2:
                neighborId = triangleId + tshift_[0] - 1;
                break;
            }
          }
        } else {
          if(p[0] == 0 and p[1] == 0)
            neighborId = triangleId + 1;
          else if(p[0] == 0) {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId + 1;
                break;
              case 1:
                neighborId = triangleId - tshift_[0] + 1;
                break;
            }
          } else if(p[1] == 0) {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId + 1;
                break;
              case 1:
                neighborId = triangleId - 1;
                break;
            }
          } else {
            switch(localNeighborId) {
              case 0:
                neighborId = triangleId + 1;
                break;
              case 1:
                neighborId = triangleId - 1;
                break;
              case 2:
                neighborId = triangleId - tshift_[0] + 1;
                break;
            }
          }
        }
      }

      return 0;
    }

    inline SimplexId
      getTriangleNeighborNumber(const SimplexId &triangleId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(triangleId < 0 or triangleId >= triangleNumber_)
        return -1;
#endif

      if(dimensionality_ == 2) {
        const auto &p = this->underlying().getTriangleCoords(triangleId);
        const SimplexId id = triangleId % 2;

        if(id) {
          if(p[0] / 2 == nbvoxels_[Di_] - 1 and p[1] == nbvoxels_[Dj_] - 1)
            return 1;
          else if(p[0] / 2 == nbvoxels_[Di_] - 1 or p[1] == nbvoxels_[Dj_] - 1)
            return 2;
          else
            return 3;
        } else {
          if(p[0] == 0 and p[1] == 0)
            return 1;
          else if(p[0] == 0 or p[1] == 0)
            return 2;
          else
            return 3;
        }
      }

      return 0;
    }

    inline int getTetrahedronVertex(const SimplexId &tetId,
                                    const int &localVertexId,
                                    SimplexId &vertexId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(tetId < 0 or tetId >= tetrahedronNumber_)
        return -1;
      if(localVertexId < 0 or localVertexId >= 4)
        return -2;
#endif

      vertexId = -1;

      if(dimensionality_ == 3) {
        const SimplexId id = tetId % 6;
        const auto &c = this->underlying().getTetrahedronCoords(tetId);
        const auto p{c.data()};

        switch(id) {
          case 0:
            vertexId = getTetrahedronVertexABCG(p, localVertexId);
            break;
          case 1:
            vertexId = getTetrahedronVertexBCDG(p, localVertexId);
            break;
          case 2:
            vertexId = getTetrahedronVertexABEG(p, localVertexId);
            break;
          case 3:
            vertexId = getTetrahedronVertexBEFG(p, localVertexId);
            break;
          case 4:
            vertexId = getTetrahedronVertexBFGH(p, localVertexId);
            break;
          case 5:
            vertexId = getTetrahedronVertexBDGH(p, localVertexId);
            break;
        }
      }
      return 0;
    }

    inline int getTetrahedronEdge(const SimplexId &tetId,
                                  const int &localEdgeId,
                                  SimplexId &edgeId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(tetId < 0 or tetId >= tetrahedronNumber_)
        return -1;
      if(localEdgeId < 0 or localEdgeId >= 6)
        return -2;
#endif

      edgeId = -1;

      if(dimensionality_ == 3) {
        const SimplexId id = tetId % 6;
        const auto &c = this->underlying().getTetrahedronCoords(tetId);
        const auto p{c.data()};

        switch(id) {
          case 0:
            edgeId = getTetrahedronEdgeABCG(p, localEdgeId);
            break;
          case 1:
            edgeId = getTetrahedronEdgeBCDG(p, localEdgeId);
            break;
          case 2:
            edgeId = getTetrahedronEdgeABEG(p, localEdgeId);
            break;
          case 3:
            edgeId = getTetrahedronEdgeBEFG(p, localEdgeId);
            break;
          case 4:
            edgeId = getTetrahedronEdgeBFGH(p, localEdgeId);
            break;
          case 5:
            edgeId = getTetrahedronEdgeBDGH(p, localEdgeId);
            break;
        }
      }

      return 0;
    }

    inline int getTetrahedronTriangle(const SimplexId &tetId,
                                      const int &localTriangleId,
                                      SimplexId &triangleId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(tetId < 0 or tetId >= tetrahedronNumber_)
        return -1;
      if(localTriangleId < 0 or localTriangleId >= 4)
        return -2;
#endif

      triangleId = -1;

      if(dimensionality_ == 3) {
        const SimplexId id = tetId % 6;
        const auto &c = this->underlying().getTetrahedronCoords(tetId);
        const auto p{c.data()};

        switch(id) {
          case 0:
            triangleId = getTetrahedronTriangleABCG(p, localTriangleId);
            break;
          case 1:
            triangleId = getTetrahedronTriangleBCDG(p, localTriangleId);
            break;
          case 2:
            triangleId = getTetrahedronTriangleABEG(p, localTriangleId);
            break;
          case 3:
            triangleId = getTetrahedronTriangleBEFG(p, localTriangleId);
            break;
          case 4:
            triangleId = getTetrahedronTriangleBFGH(p, localTriangleId);
            break;
          case 5:
            triangleId = getTetrahedronTriangleBDGH(p, localTriangleId);
            break;
        }
      }

      return 0;
    }

    inline SimplexId
      getTetrahedronNeighborNumber(const SimplexId &tetId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(tetId < 0 or tetId >= tetrahedronNumber_)
        return -1;
#endif

      if(dimensionality_ == 3) {
        const SimplexId id = tetId % 6;
        const auto &c = this->underlying().getTetrahedronCoords(tetId);
        const auto p{c.data()};

        switch(id) {
          case 0: // ABCG
            if(p[0] == 0 and p[2] == 0)
              return 2;
            else if(p[0] == 0 or p[2] == 0)
              return 3;
            else
              return 4;
            break;
          case 1: // BCDG
            if(p[1] == nbvoxels_[1] - 1 and p[2] == 0)
              return 2;
            else if(p[1] == nbvoxels_[1] - 1 or p[2] == 0)
              return 3;
            else
              return 4;
            break;
          case 2: // ABEG
            if(p[0] == 0 and p[1] == 0)
              return 2;
            else if(p[0] == 0 or p[1] == 0)
              return 3;
            else
              return 4;
            break;
          case 3: // BEFG
            if(p[1] == 0 and p[2] == nbvoxels_[2] - 1)
              return 2;
            else if(p[1] == 0 or p[2] == nbvoxels_[2] - 1)
              return 3;
            else
              return 4;
            break;
          case 4: // BFGH
            if(p[0] == nbvoxels_[0] - 1 and p[2] == nbvoxels_[2] - 1)
              return 2;
            else if(p[0] == nbvoxels_[0] - 1 or p[2] == nbvoxels_[2] - 1)
              return 3;
            else
              return 4;
            break;
          case 5: // BDGH
            if(p[0] == nbvoxels_[0] - 1 and p[1] == nbvoxels_[1] - 1)
              return 2;
            else if(p[0] == nbvoxels_[0] - 1 or p[1] == nbvoxels_[1] - 1)
              return 3;
            else
              return 4;
            break;
        }
      }

      return 0;
    }

    inline int getTetrahedronNeighbor(const SimplexId &tetId,
                                      const int &localNeighborId,
                                      SimplexId &neighborId) const final {

#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId < 0
         or localNeighborId >= getTetrahedronNeighborNumber(tetId))
        return -1;
#endif

      neighborId = -1;

      if(dimensionality_ == 3) {
        const SimplexId id = tetId % 6;
        const auto &c = this->underlying().getTetrahedronCoords(tetId);
        const auto p{c.data()};

        switch(id) {
          case 0:
            neighborId = getTetrahedronNeighborABCG(tetId, p, localNeighborId);
            break;
          case 1:
            neighborId = getTetrahedronNeighborBCDG(tetId, p, localNeighborId);
            break;
          case 2:
            neighborId = getTetrahedronNeighborABEG(tetId, p, localNeighborId);
            break;
          case 3:
            neighborId = getTetrahedronNeighborBEFG(tetId, p, localNeighborId);
            break;
          case 4:
            neighborId = getTetrahedronNeighborBFGH(tetId, p, localNeighborId);
            break;
          case 5:
            neighborId = getTetrahedronNeighborBDGH(tetId, p, localNeighborId);
            break;
        }
      }

      return 0;
    }
  };
} // namespace ttk

/// @cond

inline void
  ttk::ImplicitTriangulation::vertexToPosition2d(const SimplexId vertex,
                                                 SimplexId p[2]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = vertex >> div_[0];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = vertex / vshift_[0];
  }
}

inline void ttk::ImplicitTriangulation::edgeToPosition2d(const SimplexId edge,
                                                         const int k,
                                                         SimplexId p[2]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = e / eshift_[2 * k];
}

inline void
  ttk::ImplicitTriangulation::triangleToPosition2d(const SimplexId triangle,
                                                   SimplexId p[2]) const {
  p[0] = triangle % tshift_[0];
  p[1] = triangle / tshift_[0];
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dA(const SimplexId p[2],
                                               const int id) const {
  // V(a)={b,c}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dB(const SimplexId p[2],
                                               const int id) const {
  // V(b)={a,c,d}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // ba-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dC(const SimplexId p[2],
                                               const int id) const {
  // V(c)={a,b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dD(const SimplexId p[2],
                                               const int id) const {
  // V(d)={c,b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
    case 1:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dAB(const SimplexId p[2],
                                                const int id) const {
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // ba-L
    case 1:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 3:
      return p[0] + p[1] * eshift_[0]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dCD(const SimplexId p[2],
                                                const int id) const {
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 3:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dAC(const SimplexId p[2],
                                                const int id) const {
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dBD(const SimplexId p[2],
                                                const int id) const {
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 2:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
    case 3:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdge2dABCD(const SimplexId p[2],
                                                  const int id) const {
  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{c}+V(b)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 3:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dA(const SimplexId p[2],
                                               const int /*id*/) const {
  return p[0] * 2 + p[1] * tshift_[0];
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dB(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dC(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dD(const SimplexId p[2],
                                               const int /*id*/) const {
  return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dAB(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dCD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dAC(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dBD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStar2dABCD(const SimplexId p[2],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 5:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dA(const SimplexId p[2],
                                               const int /*id*/) const {
  return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dB(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dC(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dD(const SimplexId p[2],
                                               const int /*id*/) const {
  return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dAB(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dCD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dAC(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dBD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLink2dABCD(const SimplexId p[2],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
    case 3:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 4:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 5:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_x0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xn(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xN(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_0y(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Dj_] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_ny(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_Ny(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD1_xy(const SimplexId p[3],
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
  ttk::ImplicitTriangulation::getEdgeLink2dL(const SimplexId p[2],
                                             const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return p[0] + (p[1] - 1) * vshift_[0] + 1;
    }
  } else if(p[1] == 0)
    return p[0] + vshift_[0];
  else
    return p[0] + (p[1] - 1) * vshift_[0] + 1;
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLink2dH(const SimplexId p[2],
                                             const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return p[0] + (p[1] + 1) * vshift_[0] - 1;
    }
  } else if(p[0] == 0)
    return p[1] * vshift_[0] + 1;
  else
    return p[0] + (p[1] + 1) * vshift_[0] - 1;
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLink2dD1(const SimplexId p[2],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0];
    case 1:
      return p[0] + (p[1] + 1) * vshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStar2dL(const SimplexId p[2],
                                             const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    if(id == 0)
      return p[0] * 2 + p[1] * tshift_[0];
    else
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  } else if(p[1] == 0)
    return p[0] * 2 + p[1] * tshift_[0];
  else
    return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStar2dH(const SimplexId p[2],
                                             const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    if(id == 0)
      return p[0] * 2 + p[1] * tshift_[0];
    else
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
  } else if(p[0] == 0)
    return p[0] * 2 + p[1] * tshift_[0];
  else
    return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
}

inline void ttk::ImplicitTriangulation::vertexToPosition(const SimplexId vertex,
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

inline void ttk::ImplicitTriangulation::edgeToPosition(const SimplexId edge,
                                                       const int k,
                                                       SimplexId p[3]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = (e % eshift_[2 * k + 1]) / eshift_[2 * k];
  p[2] = e / eshift_[2 * k + 1];
}

inline void ttk::ImplicitTriangulation::triangleToPosition(
  const SimplexId triangle, const int k, SimplexId p[3]) const {
  const SimplexId t = (k) ? triangle - tsetshift_[k - 1] : triangle;
  p[0] = t % tshift_[2 * k];
  p[1] = (t % tshift_[2 * k + 1]) / tshift_[2 * k];
  p[2] = t / tshift_[2 * k + 1];
}

inline void
  ttk::ImplicitTriangulation::tetrahedronToPosition(const SimplexId tetrahedron,
                                                    SimplexId p[3]) const {
  p[0] = (tetrahedron % tetshift_[0]) / 6;
  p[1] = (tetrahedron % tetshift_[1]) / tetshift_[0];
  p[2] = tetrahedron / tetshift_[1];
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeA(const SimplexId p[3],
                                             const int id) const {
  // V(a)={b,c,e,g}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeB(const SimplexId p[3],
                                             const int id) const {
  // V(b)={a,c,d,e,f,g,h}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeC(const SimplexId p[3],
                                             const int id) const {
  // V(c)={a,b,d,g}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ca-H
    case 1:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeD(const SimplexId p[3],
                                             const int id) const {
  // V(d)={b,c,g,h}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeE(const SimplexId p[3],
                                             const int id) const {
  // V(e)={a,b,f,g}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // ea-P
    case 1:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ef-L
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // eg-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeF(const SimplexId p[3],
                                             const int id) const {
  // V(f)={b,e,g,h}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeG(const SimplexId p[3],
                                             const int id) const {
  // V(g)={a,b,c,d,e,f,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeH(const SimplexId p[3],
                                             const int id) const {
  // V(h)={b,d,f,g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 1:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // hd-P
    case 2:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // hf-H
    case 3:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeAB(const SimplexId p[3],
                                              const int id) const {
  // V(ab)=V(b)+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[11] + p[2] * eshift_[12]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeCD(const SimplexId p[3],
                                              const int id) const {
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 4:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 5:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeEF(const SimplexId p[3],
                                              const int id) const {
  // V(fe)=V(f)+V(e)::{b,f}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 5:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ef-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeGH(const SimplexId p[3],
                                              const int id) const {
  // V(gh)=V(g)+V(h)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeAC(const SimplexId p[3],
                                              const int id) const {
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ca-H
    case 1:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 5:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeBD(const SimplexId p[3],
                                              const int id) const {
  // V(bd)=V(b)+V(d)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeEG(const SimplexId p[3],
                                              const int id) const {
  // V(eg)=V(g)+V(e)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // eg-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeFH(const SimplexId p[3],
                                              const int id) const {
  // V(fh)=V(f)+V(h)::{b,f}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
    case 4:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 5:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // hf-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeAE(const SimplexId p[3],
                                              const int id) const {
  // V(ae)=V(a)+V(e)::{a,b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // ea-P
    case 5:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeBF(const SimplexId p[3],
                                              const int id) const {
  // V(bf)=V(b)+V(f)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeCG(const SimplexId p[3],
                                              const int id) const {
  // V(cg)=V(g)+V(c)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeDH(const SimplexId p[3],
                                              const int id) const {
  // V(dh)=V(d)+V(h)::{b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 4:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // hd-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeABDC(const SimplexId p[3],
                                                const int id) const {
  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 8:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeEFHG(const SimplexId p[3],
                                                const int id) const {
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 8:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 9:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeAEGC(const SimplexId p[3],
                                                const int id) const {
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 8:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeBFHD(const SimplexId p[3],
                                                const int id) const {
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 8:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 9:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeAEFB(const SimplexId p[3],
                                                const int id) const {
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 8:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeGHDC(const SimplexId p[3],
                                                const int id) const {
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 8:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexEdgeABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 8:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 10:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 11:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 12:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 13:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleA(const SimplexId /*p*/[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return 0;
    case 1:
      return tsetshift_[0];
    case 2:
      return tsetshift_[1];
    case 3:
      return tsetshift_[3];
    case 4:
      return tsetshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleB(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2;
    case 5:
      return tsetshift_[1] + p[0] * 2 + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2;
    case 11:
      return (p[0] - 1) * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleC(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[1] - 1) * tshift_[0];
    case 1:
      return (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 3:
      return tsetshift_[0] + p[1] * tshift_[2];
    case 4:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleD(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleE(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 1:
      return tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return tsetshift_[1] + (p[2] - 1) * tshift_[5] + 1;
    case 3:
      return p[2] * tshift_[1];
    case 4:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleF(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleG(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleH(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleAB(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2;
    case 5:
      return tsetshift_[1] + p[0] * 2 + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2;
    case 11:
      return (p[0] - 1) * 2;
    case 12:
      return p[0] * 2;
    case 13:
      return tsetshift_[0] + p[0] * 2;
    case 14:
      return tsetshift_[3] + p[0] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleCD(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 8:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleEF(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 5:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 6:
      return p[0] * 2 + tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 7:
      return p[0] * 2 + p[2] * tshift_[1];
    case 8:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleGH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 6:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 7:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 8:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 9:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 10:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 11:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 13:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 14:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleAC(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[1] - 1) * tshift_[0];
    case 1:
      return (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 3:
      return tsetshift_[0] + p[1] * tshift_[2];
    case 4:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4];
    case 5:
      return p[1] * tshift_[0];
    case 6:
      return tsetshift_[1] + p[1] * tshift_[4];
    case 7:
      return tsetshift_[3] + p[1] * tshift_[8];
    case 8:
      return tsetshift_[1] + p[1] * tshift_[4] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleBD(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10];
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 8:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8] + 1;
    case 9:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 10:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + 1;
    case 11:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10] + 1;
    case 12:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6] + 1;
    case 13:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8];
    case 14:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleEG(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 13:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 14:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleFH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 7:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleAE(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 1:
      return tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return tsetshift_[1] + (p[2] - 1) * tshift_[5] + 1;
    case 3:
      return p[2] * tshift_[1];
    case 4:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return tsetshift_[0] + p[2] * tshift_[3];
    case 6:
      return tsetshift_[1] + p[2] * tshift_[5];
    case 7:
      return tsetshift_[3] + p[2] * tshift_[9];
    case 8:
      return tsetshift_[1] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleBF(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 5:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11];
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7];
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9] + 1;
    case 8:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5];
    case 9:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 11:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7] + 1;
    case 13:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9];
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleCG(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 13:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 14:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleDH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 6:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 7:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 8:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleABDC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 12:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 13:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 15:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 16:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 17:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 18:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2];
    case 19:
      return p[0] * 2 + p[1] * tshift_[0];
    case 20:
      return p[0] * 2 + tsetshift_[3] + p[1] * tshift_[8];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleEFHG(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 1:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 5:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 6:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 11:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 13:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 15:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 16:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 17:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 18:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 20:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleAEGC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 13:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 14:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
    case 15:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 16:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 17:
      return tsetshift_[3] + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 18:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
    case 19:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 20:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleBFHD(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 13:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 16:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 17:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 18:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 19:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 20:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleAEFB(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 12:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 13:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 15:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 16:
      return p[0] * 2 + tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 17:
      return p[0] * 2 + p[2] * tshift_[1];
    case 18:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[0] + p[2] * tshift_[3];
    case 20:
      return p[0] * 2 + tsetshift_[3] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleGHDC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 1:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 5:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 6:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 11:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 13:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 15:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 16:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 17:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 18:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 20:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexTriangleABCDEFGH(const SimplexId p[3],
                                                        const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 13:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 14:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 15:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 16:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 17:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 18:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 22:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 23:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 24:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 25:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 26:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 27:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 28:
      return p[0] * 2 + tsetshift_[1] + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 29:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 30:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 31:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 32:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 33:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 34:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 35:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkA(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkB(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkC(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkD(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkE(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkF(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkG(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkH(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkAB(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkCD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkEF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkGH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkAC(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkBD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkEG(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 7:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkFH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 3:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkAE(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkBF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkCG(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 2:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 7:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkDH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 3:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkABDC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkEFHG(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 4:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 5:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 6:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 7:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 8:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 9:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 10:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 11:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkAEGC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 4:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 5:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 6:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 7:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkBFHD(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 9:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkAEFB(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 4:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 5:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 6:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 7:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 8:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 9:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 10:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkGHDC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexLinkABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 8:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 9:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 12:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 13:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 14:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 16:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 17:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 18:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 19:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 22:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 23:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarA(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // abcg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarB(const SimplexId p[3],
                                             const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + id;
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarC(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]; // abcg
    case 1:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarD(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // bcdg
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarE(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // abeg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarF(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // befg
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarG(const SimplexId p[3],
                                             const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id;
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarH(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarAB(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarCD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 2:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 3:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarEF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarGH(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarAC(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 1:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarBD(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarEG(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarFH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 2:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 3:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarAE(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarBF(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 7:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarCG(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarDH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 2:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 3:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarABDC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarEFHG(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 8:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 9:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarAEGC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 8:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 9:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarBFHD(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 7:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 8:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 9:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 10:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 11:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarAEFB(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 8:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 9:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarGHDC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 10:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 11:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getVertexStarABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  if(id >= 6 && id <= 11)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + id
           - 6; // tet(g)
  switch(id) {
    case 12:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 13:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 14:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 15:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 16:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 17:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 18:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 19:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
    case 20:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 21:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 22:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 23:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_x00(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2;
    case 1:
      return tsetshift_[0] + p[0] * 2;
    case 2:
      return tsetshift_[3] + p[0] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_x0n(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[2] * tshift_[1];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[2] * tshift_[9];
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_x0N(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[2] * tshift_[1];
    case 1:
      return tsetshift_[0] + p[0] * 2 + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xn0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xnn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 2:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 3:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 4:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 5:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xnN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 3:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xN0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xNn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleL_xNN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_0y0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[1] * tshift_[0];
    case 1:
      return tsetshift_[1] + p[1] * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_0yn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
    case 3:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_0yN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_ny0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 3:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_nyn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 5:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_nyN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 3:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_Ny0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_Nyn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleH_NyN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_00z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[1] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_0nz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 3:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_0Nz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_n0z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_nnz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 4:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 5:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_nNz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 2:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_N0z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_Nnz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleP_NNz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD1_xy0(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD1_xyn(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD1_xyN(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD2_0yz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
    case 2:
      return tsetshift_[3] + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD2_nyz(const SimplexId p[3],
                                                    const int id) const {
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
             + p[2] * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD2_Nyz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 2:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD3_x0z(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3] + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[2] * tshift_[7] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD3_xnz(const SimplexId p[3],
                                                    const int id) const {
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
             + p[2] * tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD3_xNz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeTriangleD4_xyz(const SimplexId p[3],
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
  ttk::ImplicitTriangulation::getEdgeLinkL(const SimplexId p[3],
                                           const int id) const {
  if(p[2] == 0 and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + eshift_[3]; // EG
    }
  } else if(p[2] == 0 and p[1] == nbvoxels_[1])
    return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]; // BG
  else if(p[2] == nbvoxels_[2] and p[1] == 0)
    return esetshift_[5] + p[0] + (p[2] - 1) * eshift_[13]; // BG
  else if(p[2] == nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
    }
  } else if(p[2] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2] + eshift_[3]; // EG
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]; // BG
    }
  } else if(p[2] == nbvoxels_[2] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
        // BG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + p[2] * eshift_[5] + eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + (p[2] + 1) * eshift_[3]; // EG
      case 2:
        return esetshift_[5] + p[0] + (p[2] - 1) * eshift_[13]; // BG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13];
        // BG
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
        // CG
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3];
        // EG
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13];
      case 3:
        return esetshift_[1] + p[0] + 1 + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5]; // CG
      case 4:
        return esetshift_[0] + p[0] + 1 + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3]; // EG
      case 5:
        return esetshift_[5] + p[0] + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkH(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[2] == 0)
    return esetshift_[5] + p[1] * eshift_[12];
  else if(p[0] == nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] - 1;
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + eshift_[1] - 1;
    }
  } else if(p[0] == 0 and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 1:
        return esetshift_[1] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5] + 1;
    }
  } else if(p[0] == nbvoxels_[0] and p[2] == nbvoxels_[2])
    return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
           + (p[2] - 1) * eshift_[13];
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[0] - 1 + (p[1] + 1) * eshift_[0] + eshift_[1];
      case 1:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4];
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12];
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 1:
        return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 2:
        return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
               + (p[2] - 1) * eshift_[5];
    }
  } else if(p[0] == 0 and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5] + 1;
      case 1:
        return esetshift_[5] + p[1] * eshift_[12] + p[2] * eshift_[13];
      case 2:
        return p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
    }
  } else if(p[0] == nbvoxels_[0] and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 1:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
      case 2:
        return p[0] - 1 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1];
    }
  } else {
    switch(id) {
      case 0:
        return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 2:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
      case 3:
        return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
               + (p[2] - 1) * eshift_[5];
      case 4:
        return esetshift_[5] + (p[0] - 1) + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 5:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkP(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[1] == 0)
    return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  else if(p[0] == 0 and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 1:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == nbvoxels_[1])
    return esetshift_[5] + p[0] - 1 + (p[1] - 1) * eshift_[12]
           + p[2] * eshift_[13];
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
    }
  } else if(p[0] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  } else if(p[0] == nbvoxels_[0] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 2:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
      case 3:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
      case 4:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 5:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkD1(const SimplexId p[3],
                                            const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[4] + p[0] + p[1] * eshift_[10]
               + (p[2] - 1) * eshift_[11];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11];
      case 2:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 3:
        return esetshift_[3] + p[0] + p[1] * eshift_[8]
               + (p[2] - 1) * eshift_[9] + 1;
    }
  } else if(p[2] == 0) {
    switch(id) {
      case 0:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[3] + p[0] + p[1] * eshift_[8]
               + (p[2] - 1) * eshift_[9] + 1;
      case 1:
        return esetshift_[4] + p[0] + p[1] * eshift_[10]
               + (p[2] - 1) * eshift_[11];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkD2(const SimplexId p[3],
                                            const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[0]) {
    switch(id) {
      case 0:
        return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11] - 1;
      case 2:
        return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
      case 3:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7] - 1;
    }
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
      case 1:
        return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7] - 1;
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11] - 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkD3(const SimplexId p[3],
                                            const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
               + p[2] * eshift_[7];
      case 1:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7];
      case 2:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 3:
        return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
               + p[2] * eshift_[9] + 1;
    }
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7];
      case 1:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
               + p[2] * eshift_[7];
      case 1:
        return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
               + p[2] * eshift_[9] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeLinkD4(const SimplexId p[3],
                                            const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 5:
      return esetshift_[0] + p[0] + 1 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarL(const SimplexId p[3],
                                           const int id) const {
  if(p[2] == 0 and p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] * 6; // ABCG
      case 1:
        return p[0] * 6 + 2; // ABEG
    }
  } else if(p[2] == 0 and p[1] == nbvoxels_[1])
    return (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1; // BCDG
  else if(p[2] == nbvoxels_[2] and p[1] == 0)
    return (p[2] - 1) * tetshift_[1] + p[0] * 6 + 3; // BEFG
  else if(p[2] == nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + p[0] * 6
               + 4;
        // BFGH
      case 1:
        return (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + p[0] * 6
               + 5;
        // BDGH
    }
  } else if(p[2] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1; // BCDG
    }
  } else if(p[2] == nbvoxels_[2] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 1:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[0] * 6 + 2; // ABEG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[0] * 6 + 3; // BEFG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 1:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 4:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 5:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarH(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[2] == 0)
    return p[1] * tetshift_[0]; // ABCG
  else if(p[0] == nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1; // BCDG
      case 1:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5; // BDGH
    }
  } else if(p[0] == 0 and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 3; // BEFG
    }
  } else if(p[0] == nbvoxels_[0] and p[2] == nbvoxels_[2])
    return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
           + 4; // BFGH
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1; // BCDG
      case 1:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5; // BDGH
      case 2:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
    }
  } else if(p[0] == 0 and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 3; // BEFG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0]; // ABCG
    }
  } else if(p[0] == nbvoxels_[0] and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 1; // BCDG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 2; // ABEG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
      case 4:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 1; // BCDG
      case 5:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarP(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[1] == 0)
    return p[2] * tetshift_[1] + 2; // ABEG
  else if(p[0] == 0 and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + 1; // BCDG
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 4; // BFGH
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == nbvoxels_[1])
    return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
           + 5; // BDGH
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 4; // BFGH
      case 2:
        return p[2] * tetshift_[1] + p[0] * 6 + 2; // ABEG
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]
               + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
    }
  } else if(p[0] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + 1; // BCDG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
    }
  } else if(p[0] == nbvoxels_[0] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]
               + p[0] * 6; // ABCG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 4:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 3; // BEFG
      case 5:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarD1(const SimplexId p[3],
                                            const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 1; // BCDG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 4; // BFGH
    }
  } else if(p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[1] * tetshift_[0] + p[0] * 6 + 1; // BCDG
    }
  } else {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarD2(const SimplexId p[3],
                                            const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[0]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 3:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getEdgeStarD3(const SimplexId p[3],
                                            const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3; // BEFG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5; // BDGH
    }
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3; // BEFG
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexF(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1;
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexH(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexC(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
    else
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  } else {
    if(id == 0)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexD1(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexD2(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleVertexD3(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeF_0(const SimplexId p[3],
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
  ttk::ImplicitTriangulation::getTriangleEdgeF_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1;
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeH_0(const SimplexId p[3],
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
  ttk::ImplicitTriangulation::getTriangleEdgeH_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1;
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeC_0(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeC_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeD1_0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1;
    case 1:
      return esetshift_[4] + p[0] / 2 + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeD1_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
    case 1:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeD2_0(const SimplexId p[3],
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
  ttk::ImplicitTriangulation::getTriangleEdgeD2_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9]
             + 1;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeD3_0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleEdgeD3_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkF(const SimplexId p[3],
                                               const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1;
    }
  } else if(p[2] == 0)
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + vshift_[1];
  else
    return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkH(const SimplexId p[3],
                                               const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1;
    }
  } else if(p[1] == 0)
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
  else
    return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkC(const SimplexId p[3],
                                               const int id) const {
  if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1;
    }
  } else if(p[0] < 2)
    return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
  else
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkD1(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkD2(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleLinkD3(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarF(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[2] > 0 and p[2] < nbvoxels_[2]) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 1; // BCDG
        case 1:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0]
                 + (p[2] - 1) * tetshift_[1] + 4; // BFGH
      }
    } else if(p[2] == 0)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // BCDG
    else
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // BFGH
  } else {
    if(p[2] > 0 and p[2] < nbvoxels_[2]) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
        case 1:
          return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
                 + 3; // BEFG
      }
    } else if(p[2] == 0)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
    else
      return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // BEFG
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarH(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[1] > 0 and p[1] < nbvoxels_[1]) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 3; // BEFG
        case 1:
          return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0]
                 + p[2] * tetshift_[1] + 5; // BDGH
      }
    } else if(p[1] == 0)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 3; // BEFG
    else
      return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // BDGH
  } else {
    if(p[1] > 0 and p[1] < nbvoxels_[1]) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 2; // ABEG
        case 1:
          return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
                 + 1; // BCDG
      }
    } else if(p[1] == 0)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // ABEG
    else
      return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // BCDG
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarC(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 2; // ABEG
        case 1:
          return ((p[0] - 2) / 2) * 6 + p[1] * tetshift_[0]
                 + p[2] * tetshift_[1] + 4; // BFGH
      }
    } else if(p[0] < 2)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // ABEG
    else
      return ((p[0] - 2) / 2) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 4; // BFGH
  } else {
    if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
        case 1:
          return ((p[0] - 1) / 2) * 6 + p[1] * tetshift_[0]
                 + p[2] * tetshift_[1] + 5; // BDGH
      }
    } else if(p[0] < 2)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
    else
      return ((p[0] - 1) / 2) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // BDGH
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarD1(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 2; // ABEG
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 3; // BEFG
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1; // BCDG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarD2(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 5; // BDGH
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4; // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // ABEG
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTriangleStarD3(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 3; // BEFG
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4; // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1; // BCDG
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexABCG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1]; // a
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]; // c
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexBCDG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]; // c
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1; // d
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexABEG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1]; // a
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]; // e
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexBEFG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]; // e
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1; // f
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexBFGH(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1; // f
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1; // h
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronVertexBDGH(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1; // d
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1; // h
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeABCG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + p[2] * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeBCDG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + p[2] * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeABEG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeBEFG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeBFGH(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronEdgeBDGH(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronTriangleABCG(const SimplexId p[3],
                                                         const int id) const {
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
  ttk::ImplicitTriangulation::getTetrahedronTriangleBCDG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronTriangleABEG(const SimplexId p[3],
                                                         const int id) const {
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
  ttk::ImplicitTriangulation::getTetrahedronTriangleBEFG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronTriangleBFGH(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 3;
    case 3:
      return p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ImplicitTriangulation::getTetrahedronTriangleBDGH(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 2;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborABCG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t + 1;
    case 1:
      return t + 2;
    case 2:
      if(p[0] > 0)
        return t - 1;
      else
        return t - tetshift_[1] + 3;
    case 3:
      return t - tetshift_[1] + 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborBCDG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 4;
    case 2:
      if(p[2] > 0)
        return t - tetshift_[1] + 3;
      else
        return t + tetshift_[0] + 1;
    case 3:
      return t + tetshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborABEG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 2;
    case 1:
      return t + 1;
    case 2:
      if(p[0] > 0)
        return t - 4;
      else
        return t - tetshift_[0] - 1;
    case 3:
      return t - tetshift_[0] - 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborBEFG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      if(p[1] > 0)
        return t - tetshift_[0] + 2;
      else
        return t + tetshift_[1] - 3;
    case 3:
      return t + tetshift_[1] - 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborBFGH(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      if(p[0] < nbvoxels_[0] - 1)
        return t + 4;
      else
        return t + tetshift_[1] - 3;
    case 3:
      return t + tetshift_[1] - 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ImplicitTriangulation::getTetrahedronNeighborBDGH(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t - 4;
    case 2:
      if(p[0] < nbvoxels_[0] - 1)
        return t + 1;
      else
        return t + tetshift_[0] - 2;
    case 3:
      return t + tetshift_[0] - 2;
  }
  return -1;
}

#include <ImplicitPreconditions.h>

/// @endcond
