#pragma once

#include <ImplicitTriangulation.h>

namespace ttk {
  /**
   * @brief Implicit Triangulation class with preconditioning
   */
  class ImplicitWithPreconditions final : public ImplicitTriangulation {
  public:
    int preconditionVerticesInternal() override;
    int preconditionEdgesInternal() override;
    int preconditionTrianglesInternal() override;
    int preconditionTetrahedronsInternal() override;

    inline VertexPosition getVertexPosition(const SimplexId v) const override {
      return this->vertexPositions_[v];
    }
    inline std::array<SimplexId, 3>
      getVertexCoords(const SimplexId v) const override {
      return this->vertexCoords_[v];
    }
    inline EdgePosition getEdgePosition(const SimplexId e) const override {
      return this->edgePositions_[e];
    }
    inline std::array<SimplexId, 3>
      getEdgeCoords(const SimplexId e) const override {
      return this->edgeCoords_[e];
    }
    inline TrianglePosition
      getTrianglePosition(const SimplexId t) const override {
      return this->trianglePositions_[t];
    }
    inline std::array<SimplexId, 3>
      getTriangleCoords(const SimplexId t) const override {
      return this->triangleCoords_[t];
    }
    inline std::array<SimplexId, 3>
      getTetrahedronCoords(const SimplexId t) const override {
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
  class ImplicitNoPreconditions final : public ImplicitTriangulation {
  public:
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

    VertexPosition getVertexPosition(const SimplexId v) const override;

    inline std::array<SimplexId, 3>
      getVertexCoords(const SimplexId v) const override {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 2) {
        this->vertexToPosition2d(v, p.data());
      } else if(this->dimensionality_ == 3) {
        this->vertexToPosition(v, p.data());
      }
      return p;
    }

    EdgePosition getEdgePosition(const SimplexId e) const override;
    std::array<SimplexId, 3> getEdgeCoords(const SimplexId e) const override;
    TrianglePosition getTrianglePosition(const SimplexId t) const override;
    std::array<SimplexId, 3>
      getTriangleCoords(const SimplexId t) const override;

    inline std::array<SimplexId, 3>
      getTetrahedronCoords(const SimplexId t) const override {
      std::array<SimplexId, 3> p{};
      this->tetrahedronToPosition(t, p.data());
      return p;
    }
  };

} // namespace ttk
