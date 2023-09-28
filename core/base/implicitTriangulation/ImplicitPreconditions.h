#pragma once

#include <ImplicitTriangulation.h>

namespace ttk {
  /**
   * @brief Implicit Triangulation class with preconditioning
   */
  class ImplicitWithPreconditions final
    : public ImplicitTriangulationCRTP<ImplicitWithPreconditions> {
  public:
    ImplicitWithPreconditions() {
      this->setDebugMsgPrefix("ImplicitTriangulationWithPreconditions");
    }

    int preconditionVerticesInternal() override;
    int preconditionEdgesInternal() override;
    int preconditionTrianglesInternal() override;
    int preconditionTetrahedronsInternal() override;

    inline VertexPosition getVertexPosition(const SimplexId v) const {
      return this->vertexPositions_[v];
    }
    inline std::array<SimplexId, 3> const &
      getVertexCoords(const SimplexId v) const {
      return this->vertexCoords_[v];
    }
    inline EdgePosition getEdgePosition(const SimplexId e) const {
      return this->edgePositions_[e];
    }
    inline std::array<SimplexId, 3> const &
      getEdgeCoords(const SimplexId e) const {
      return this->edgeCoords_[e];
    }
    inline TrianglePosition getTrianglePosition(const SimplexId t) const {
      return this->trianglePositions_[t];
    }
    inline std::array<SimplexId, 3> const &
      getTriangleCoords(const SimplexId t) const {
      return this->triangleCoords_[t];
    }
    inline const std::array<SimplexId, 3> &
      getTetrahedronCoords(const SimplexId t) const {
      return this->tetrahedronCoords_[t];
    }

    inline int clear() {
      vertexPositions_ = std::vector<VertexPosition>{};
      vertexCoords_ = std::vector<std::array<SimplexId, 3>>{};
      edgePositions_ = std::vector<EdgePosition>{};
      edgeCoords_ = std::vector<std::array<SimplexId, 3>>{};
      trianglePositions_ = std::vector<TrianglePosition>{};
      triangleCoords_ = std::vector<std::array<SimplexId, 3>>{};
      tetrahedronCoords_ = std::vector<std::array<SimplexId, 3>>{};
      hasPreconditionedVerticesAndCells_ = false;
      return AbstractTriangulation::clear();
    }

  private:
    // for every vertex, its position on the grid
    std::vector<VertexPosition> vertexPositions_{};
    // for  every vertex, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> vertexCoords_{};
    // for every edge, its position on the grid
    std::vector<EdgePosition> edgePositions_{};
    // for every edge, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> edgeCoords_{};
    // for every triangle, its position on the grid
    std::vector<TrianglePosition> trianglePositions_{};
    // for every triangle, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> triangleCoords_{};
    // for every tetrahedron, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> tetrahedronCoords_{};
  };

  /**
   * @brief Implicit Triangulation class without preconditioning
   */
  class ImplicitNoPreconditions final
    : public ImplicitTriangulationCRTP<ImplicitNoPreconditions> {
  public:
    ImplicitNoPreconditions() {
      this->setDebugMsgPrefix("ImplicitTriangulationNoPreconditions");
    }

    inline int preconditionVerticesInternal() override {
      return 0;
    }
    inline int preconditionEdgesInternal() override {
      return 0;
    }
    inline int preconditionTrianglesInternal() override {
      return 0;
    }
    inline int preconditionTetrahedronsInternal() override {
      return 0;
    }

    inline VertexPosition getVertexPosition(const SimplexId v) const {
      if(this->dimensionality_ == 1) {
        if(v == 0) {
          return VertexPosition::LEFT_CORNER_1D;
        } else if(v == this->vertexNumber_ - 1) {
          return VertexPosition::RIGHT_CORNER_1D;
        }
        return VertexPosition::CENTER_1D;
      } else if(this->dimensionality_ == 2) {
        const auto p{this->getVertexCoords(v)};
        if(0 < p[0] and p[0] < this->nbvoxels_[this->Di_]) {
          if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
            return VertexPosition::CENTER_2D;
          else if(p[1] == 0)
            return VertexPosition::TOP_EDGE_2D; // ab
          else
            return VertexPosition::BOTTOM_EDGE_2D; // cd
        } else if(p[0] == 0) {
          if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
            return VertexPosition::LEFT_EDGE_2D; // ac
          else if(p[1] == 0)
            return VertexPosition::TOP_LEFT_CORNER_2D; // a
          else
            return VertexPosition::BOTTOM_LEFT_CORNER_2D; // c
        } else {
          if(0 < p[1] and p[1] < this->nbvoxels_[this->Dj_])
            return VertexPosition::RIGHT_EDGE_2D; // bd
          else if(p[1] == 0)
            return VertexPosition::TOP_RIGHT_CORNER_2D; // b
          else
            return VertexPosition::BOTTOM_RIGHT_CORNER_2D; // d
        }

      } else if(this->dimensionality_ == 3) {
        const auto p{this->getVertexCoords(v)};
        if(0 < p[0] and p[0] < this->nbvoxels_[0]) {
          if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::CENTER_3D;
            else if(p[2] == 0)
              return VertexPosition::FRONT_FACE_3D; // abcd
            else
              return VertexPosition::BACK_FACE_3D; // efgh
          } else if(p[1] == 0) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::TOP_FACE_3D; // abef
            else if(p[2] == 0)
              return VertexPosition::TOP_FRONT_EDGE_3D; // ab
            else
              return VertexPosition::TOP_BACK_EDGE_3D; // ef
          } else {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::BOTTOM_FACE_3D; // cdgh
            else if(p[2] == 0)
              return VertexPosition::BOTTOM_FRONT_EDGE_3D; // cd
            else
              return VertexPosition::BOTTOM_BACK_EDGE_3D; // gh
          }
        } else if(p[0] == 0) {
          if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::LEFT_FACE_3D; // aceg
            else if(p[2] == 0)
              return VertexPosition::LEFT_FRONT_EDGE_3D; // ac
            else
              return VertexPosition::LEFT_BACK_EDGE_3D; // eg
          } else if(p[1] == 0) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::TOP_LEFT_EDGE_3D; // ae
            else if(p[2] == 0)
              return VertexPosition::TOP_LEFT_FRONT_CORNER_3D; // a
            else
              return VertexPosition::TOP_LEFT_BACK_CORNER_3D; // e
          } else {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::BOTTOM_LEFT_EDGE_3D; // cg
            else if(p[2] == 0)
              return VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D; // c
            else
              return VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D; // g
          }
        } else {
          if(0 < p[1] and p[1] < this->nbvoxels_[1]) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::RIGHT_FACE_3D; // bdfh
            else if(p[2] == 0)
              return VertexPosition::RIGHT_FRONT_EDGE_3D; // bd
            else
              return VertexPosition::RIGHT_BACK_EDGE_3D; // fh
          } else if(p[1] == 0) {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::TOP_RIGHT_EDGE_3D; // bf
            else if(p[2] == 0)
              return VertexPosition::TOP_RIGHT_FRONT_CORNER_3D; // b
            else
              return VertexPosition::TOP_RIGHT_BACK_CORNER_3D; // f
          } else {
            if(0 < p[2] and p[2] < this->nbvoxels_[2])
              return VertexPosition::BOTTOM_RIGHT_EDGE_3D; // dh
            else if(p[2] == 0)
              return VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D; // d
            else
              return VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D; // h
          }
        }
      }
      return VertexPosition::CENTER_3D;
    }

    inline std::array<SimplexId, 3> getVertexCoords(const SimplexId v) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 2) {
        this->vertexToPosition2d(v, p.data());
      } else if(this->dimensionality_ == 3) {
        this->vertexToPosition(v, p.data());
      }
      return p;
    }

    inline EdgePosition getEdgePosition(const SimplexId e) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 3) {

        if(e < this->esetshift_[0]) {
          this->edgeToPosition(e, 0, p.data());
          if(p[1] > 0 and p[1] < this->nbvoxels_[1]) {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::L_xnn_3D;
            else if(p[2] == 0)
              return EdgePosition::L_xn0_3D;
            else
              return EdgePosition::L_xnN_3D;
          } else if(p[1] == 0) {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::L_x0n_3D;
            else if(p[2] == 0)
              return EdgePosition::L_x00_3D;
            else
              return EdgePosition::L_x0N_3D;
          } else {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::L_xNn_3D;
            else if(p[2] == 0)
              return EdgePosition::L_xN0_3D;
            else
              return EdgePosition::L_xNN_3D;
          }
        } else if(e < this->esetshift_[1]) {
          this->edgeToPosition(e, 1, p.data());
          if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::H_nyn_3D;
            else if(p[2] == 0)
              return EdgePosition::H_ny0_3D;
            else
              return EdgePosition::H_nyN_3D;
          } else if(p[0] == 0) {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::H_0yn_3D;
            else if(p[2] == 0)
              return EdgePosition::H_0y0_3D;
            else
              return EdgePosition::H_0yN_3D;
          } else {
            if(p[2] > 0 and p[2] < this->nbvoxels_[2])
              return EdgePosition::H_Nyn_3D;
            else if(p[2] == 0)
              return EdgePosition::H_Ny0_3D;
            else
              return EdgePosition::H_NyN_3D;
          }
        } else if(e < this->esetshift_[2]) {
          this->edgeToPosition(e, 2, p.data());
          if(p[0] > 0 and p[0] < this->nbvoxels_[0]) {
            if(p[1] > 0 and p[1] < this->nbvoxels_[1])
              return EdgePosition::P_nnz_3D;
            else if(p[1] == 0)
              return EdgePosition::P_n0z_3D;
            else
              return EdgePosition::P_nNz_3D;
          } else if(p[0] == 0) {
            if(p[1] > 0 and p[1] < this->nbvoxels_[1])
              return EdgePosition::P_0nz_3D;
            else if(p[1] == 0)
              return EdgePosition::P_00z_3D;
            else
              return EdgePosition::P_0Nz_3D;
          } else {
            if(p[1] > 0 and p[1] < this->nbvoxels_[1])
              return EdgePosition::P_Nnz_3D;
            else if(p[1] == 0)
              return EdgePosition::P_N0z_3D;
            else
              return EdgePosition::P_NNz_3D;
          }
        } else if(e < this->esetshift_[3]) {
          this->edgeToPosition(e, 3, p.data());
          if(p[2] > 0 and p[2] < this->nbvoxels_[2])
            return EdgePosition::D1_xyn_3D;
          else if(p[2] == 0)
            return EdgePosition::D1_xy0_3D;
          else
            return EdgePosition::D1_xyN_3D;
        } else if(e < this->esetshift_[4]) {
          this->edgeToPosition(e, 4, p.data());
          if(p[0] > 0 and p[0] < this->nbvoxels_[0])
            return EdgePosition::D2_nyz_3D;
          else if(p[0] == 0)
            return EdgePosition::D2_0yz_3D;
          else
            return EdgePosition::D2_Nyz_3D;
        } else if(e < this->esetshift_[5]) {
          this->edgeToPosition(e, 5, p.data());
          if(p[1] > 0 and p[1] < this->nbvoxels_[1])
            return EdgePosition::D3_xnz_3D;
          else if(p[1] == 0)
            return EdgePosition::D3_x0z_3D;
          else
            return EdgePosition::D3_xNz_3D;
        } else if(e < this->esetshift_[6]) {
          return EdgePosition::D4_3D;
        }

      } else if(this->dimensionality_ == 2) {
        if(e < this->esetshift_[0]) {
          this->edgeToPosition2d(e, 0, p.data());
          if(p[1] > 0 and p[1] < this->nbvoxels_[this->Dj_])
            return EdgePosition::L_xn_2D;
          else if(p[1] == 0)
            return EdgePosition::L_x0_2D;
          else
            return EdgePosition::L_xN_2D;
        } else if(e < this->esetshift_[1]) {
          this->edgeToPosition2d(e, 1, p.data());
          if(p[0] > 0 and p[0] < this->nbvoxels_[this->Di_])
            return EdgePosition::H_ny_2D;
          else if(p[0] == 0)
            return EdgePosition::H_0y_2D;
          else
            return EdgePosition::H_Ny_2D;
        } else if(e < this->esetshift_[2]) {
          return EdgePosition::D1_2D;
        }

      } else if(this->dimensionality_ == 1) {
        if(e == 0) {
          return EdgePosition::FIRST_EDGE_1D;
        } else if(e == this->edgeNumber_ - 1)
          return EdgePosition::CENTER_1D;
      } else {
        return EdgePosition::LAST_EDGE_1D;
      }

      return EdgePosition::CENTER_1D;
    }

    inline std::array<SimplexId, 3> getEdgeCoords(const SimplexId e) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 3) {
        if(e < this->esetshift_[0]) {
          this->edgeToPosition(e, 0, p.data());
        } else if(e < this->esetshift_[1]) {
          this->edgeToPosition(e, 1, p.data());
        } else if(e < this->esetshift_[2]) {
          this->edgeToPosition(e, 2, p.data());
        } else if(e < this->esetshift_[3]) {
          this->edgeToPosition(e, 3, p.data());
        } else if(e < this->esetshift_[4]) {
          this->edgeToPosition(e, 4, p.data());
        } else if(e < this->esetshift_[5]) {
          this->edgeToPosition(e, 5, p.data());
        } else if(e < this->esetshift_[6]) {
          this->edgeToPosition(e, 6, p.data());
        }

      } else if(this->dimensionality_ == 2) {
        if(e < this->esetshift_[0]) {
          this->edgeToPosition2d(e, 0, p.data());
        } else if(e < this->esetshift_[1]) {
          this->edgeToPosition2d(e, 1, p.data());
        } else if(e < this->esetshift_[2]) {
          this->edgeToPosition2d(e, 2, p.data());
        }
      }
      return p;
    }

    inline TrianglePosition getTrianglePosition(const SimplexId t) const {
      if(this->dimensionality_ == 2) {
        if(t % 2 == 0) {
          return TrianglePosition::TOP_2D;
        } else {
          return TrianglePosition::BOTTOM_2D;
        }
      } else if(this->dimensionality_ == 3) {
        if(t < this->tsetshift_[0]) {
          return TrianglePosition::F_3D;
        } else if(t < this->tsetshift_[1]) {
          return TrianglePosition::H_3D;
        } else if(t < this->tsetshift_[2]) {
          return TrianglePosition::C_3D;
        } else if(t < this->tsetshift_[3]) {
          return TrianglePosition::D1_3D;
        } else if(t < this->tsetshift_[4]) {
          return TrianglePosition::D2_3D;
        } else if(t < this->tsetshift_[5]) {
          return TrianglePosition::D3_3D;
        }
      }
      return TrianglePosition::C_3D;
    }

    inline std::array<SimplexId, 3> getTriangleCoords(const SimplexId t) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 2) {
        this->triangleToPosition2d(t, p.data());
      } else if(this->dimensionality_ == 3) {
        if(t < this->tsetshift_[0]) {
          this->triangleToPosition(t, 0, p.data());
        } else if(t < this->tsetshift_[1]) {
          this->triangleToPosition(t, 1, p.data());
        } else if(t < this->tsetshift_[2]) {
          this->triangleToPosition(t, 2, p.data());
        } else if(t < this->tsetshift_[3]) {
          this->triangleToPosition(t, 3, p.data());
        } else if(t < this->tsetshift_[4]) {
          this->triangleToPosition(t, 4, p.data());
        } else if(t < this->tsetshift_[5]) {
          this->triangleToPosition(t, 5, p.data());
        }
      }
      return p;
    }

    inline std::array<SimplexId, 3>
      getTetrahedronCoords(const SimplexId t) const {
      std::array<SimplexId, 3> p{};
      this->tetrahedronToPosition(t, p.data());
      return p;
    }
  };

} // namespace ttk
