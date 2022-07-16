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
#ifdef TTK_ENABLE_MPI
      return this->preconditionDistributedVertices();
#endif // TTK_ENABLE_MPI
      return 0;
    }
    inline int preconditionEdgesInternal() override {
#ifdef TTK_ENABLE_MPI
      return this->preconditionDistributedEdges();
#endif // TTK_ENABLE_MPI
      return 0;
    }
    inline int preconditionTrianglesInternal() override {
#ifdef TTK_ENABLE_MPI
      return this->preconditionDistributedTriangles();
#endif // TTK_ENABLE_MPI
      return 0;
    }
    inline int preconditionTetrahedronsInternal() override {
      return 0;
    }

    VertexPosition getVertexPosition(const SimplexId v) const;

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
  };

} // namespace ttk
