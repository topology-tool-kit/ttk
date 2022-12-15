#pragma once

#include <PeriodicImplicitTriangulation.h>

namespace ttk {
  /**
   * @brief Periodic implicit Triangulation class with preconditioning
   */
  class PeriodicWithPreconditions final
    : public PeriodicImplicitTriangulationCRTP<PeriodicWithPreconditions> {
  public:
    PeriodicWithPreconditions() {
      this->setDebugMsgPrefix("PeriodicTriangulationWithPreconditions");
    }

    int preconditionVerticesInternal() override;
    int preconditionEdgesInternal() override;
    int preconditionTrianglesInternal() override;
    int preconditionTetrahedronsInternal() override;

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
    inline std::array<SimplexId, 3> const &
      getTetrahedronCoords(const SimplexId t) const {
      return this->tetrahedronCoords_[t];
    }
    inline SimplexId getEdgeVertexAccelerated(const SimplexId e) const {
      return this->edgeVertexAccelerated_[e];
    }

    inline int clear() {
      vertexCoords_ = std::vector<std::array<SimplexId, 3>>{};
      edgePositions_ = std::vector<EdgePosition>{};
      edgeCoords_ = std::vector<std::array<SimplexId, 3>>{};
      trianglePositions_ = std::vector<TrianglePosition>{};
      triangleCoords_ = std::vector<std::array<SimplexId, 3>>{};
      tetrahedronCoords_ = std::vector<std::array<SimplexId, 3>>{};
      edgeVertexAccelerated_ = std::vector<SimplexId>{};
      hasPreconditionedVerticesAndCells_ = false;
      return AbstractTriangulation::clear();
    }

  private:
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

    // cache some edge vertex computation wrt acceleration
    std::vector<SimplexId> edgeVertexAccelerated_{};
  };

  /**
   * @brief Periodic implicit Triangulation class without preconditioning
   */
  class PeriodicNoPreconditions final
    : public PeriodicImplicitTriangulationCRTP<PeriodicNoPreconditions> {
  public:
    PeriodicNoPreconditions() {
      this->setDebugMsgPrefix("PeriodicTriangulationNoPreconditions");
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

    inline std::array<SimplexId, 3> getVertexCoords(const SimplexId v) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 2) {
        this->vertexToPosition2d(v, p.data());
      } else if(this->dimensionality_ == 3) {
        this->vertexToPosition(v, p.data());
      }
      return p;
    }

    EdgePosition getEdgePosition(const SimplexId e) const;
    std::array<SimplexId, 3> getEdgeCoords(const SimplexId e) const;
    TrianglePosition getTrianglePosition(const SimplexId t) const;
    std::array<SimplexId, 3> getTriangleCoords(const SimplexId t) const;

    inline std::array<SimplexId, 3>
      getTetrahedronCoords(const SimplexId t) const {
      std::array<SimplexId, 3> p{};
      this->tetrahedronToPosition(t, p.data());
      return p;
    }

    inline SimplexId getEdgeVertexAccelerated(const SimplexId e) const {
      const auto p{this->getEdgeCoords(e)};
      if(this->isAccelerated_) {
        return (p[1] << this->div_[0]) + (p[2] << this->div_[1]);
      } else {
        return p[1] * this->vshift_[0] + p[2] * this->vshift_[1];
      }
    }
  };

} // namespace ttk
