#include <SurfaceGeometrySmoother.h>

#include <array>
#include <iostream>

namespace ttk {

  class Quadrangulation : virtual public Debug {
  public:
    Quadrangulation();
    ~Quadrangulation() override = default;

    int preconditionVertexNeighbors();
    int preconditionVertexStars();
    int preconditionEdges();

    inline void setInputCells(const SimplexId cellNumber,
                              void *const quadCells) {
      this->nCells_ = cellNumber;
      this->cells_ = static_cast<const Quad *>(quadCells);
    }

    inline void setInputPoints(const SimplexId pointNumber,
                               void *const pointCoords) {
      this->nVerts_ = pointNumber;
      this->vertCoords_ = static_cast<const Point *>(pointCoords);
    }

    inline int isVertexExtraordinary(const SimplexId v) {
      const size_t ordinary_valence{4};
      return this->vertexNeighbors_[v].size() != ordinary_valence;
    }

    inline void
      getVertexPoint(const SimplexId v, float &x, float &y, float &z) const {
      const auto &c{this->vertCoords_[v]};
      x = c[0];
      y = c[1];
      z = c[2];
    }

    inline void
      getCellVertex(const SimplexId c, const int l, SimplexId &v) const {
      v = this->cells_[c][l];
    }

    inline SimplexId getVertexNeighborNumber(const SimplexId v) const {
      return this->vertexNeighbors_[v].size();
    }
    inline void
      getVertexNeighbor(const SimplexId v, const int l, SimplexId &n) const {
      n = this->vertexNeighbors_[v][l];
    }

    inline SimplexId getVertexStarNumber(const SimplexId v) const {
      return this->vertexStars_[v].size();
    }
    inline SimplexId getVertexStar(const SimplexId v, const int l) const {
      return this->vertexStars_[v][l];
    }
    inline void
      getVertexStar(const SimplexId v, const int l, SimplexId &vs) const {
      vs = this->vertexStars_[v][l];
    }

    inline SimplexId getEdgeStarNumber(const SimplexId e) const {
      return this->edgeStars_[e].size();
    }
    inline SimplexId getEdgeStar(const SimplexId e, const int l) const {
      return this->edgeStars_[e][l];
    }

    inline SimplexId getCellEdge(const SimplexId c, const int l) {
      return this->quadEdges_[c][l];
    }

    const std::array<SimplexId, 2> &getEdge(const SimplexId e) const {
      return this->edges_[e];
    }

    inline SimplexId getCellVertexNumber(const SimplexId ttkNotUsed(c)) const {
      return 4;
    }
    inline int getDimensionality() const {
      return 2;
    }
    inline const SimplexId &getNumberOfVertices() const {
      return this->nVerts_;
    }
    inline const SimplexId &getNumberOfCells() const {
      return this->nCells_;
    }
    inline SimplexId getNumberOfEdges() const {
      return this->edges_.size();
    }

    inline void computeDensityAndDeformity(const SimplexId v,
                                           float &density,
                                           float &deformity) const {
      auto minDist{std::numeric_limits<float>::max()};
      auto maxDist{std::numeric_limits<float>::lowest()};
      const auto &pi{this->vertCoords_[v]};
      for(const auto j : this->vertexNeighbors_[v]) {
        const auto &pj{this->vertCoords_[j]};
        const auto dist{Geometry::distance(pi.data(), pj.data())};
        minDist = std::min(minDist, dist);
        maxDist = std::max(maxDist, dist);
      }
      density = std::exp(-minDist);
      deformity = std::exp(-minDist / maxDist);
    }

    void computeStatistics(std::vector<SimplexId> &vertsValence,
                           std::vector<float> &vertsDensity,
                           std::vector<float> &vertsDifformity,
                           std::vector<float> &quadArea,
                           std::vector<float> &quadDiagsRatio,
                           std::vector<float> &quadEdgesRatio,
                           std::vector<float> &quadAnglesRatio) const;

    /**
     * @brief Ad-hoc quad data structure (4 vertex ids)
     */
    using Quad = std::array<ttk::LongSimplexId, 4>;
    /**
     * @brief Ad-hoc vertex coordinates data structure (3 floats)
     */
    using Point = ttk::SurfaceGeometrySmoother::Point;

  private:
    CellArray buildQuadOffsets();

    const Point *vertCoords_{};
    const Quad *cells_{};
    SimplexId nVerts_{};
    SimplexId nCells_{};
    FlatJaggedArray vertexNeighbors_{};
    FlatJaggedArray vertexStars_{};
    std::vector<std::array<SimplexId, 2>> edges_{};
    FlatJaggedArray edgeStars_{};
    std::vector<std::array<SimplexId, 4>> quadEdges_{};
    std::vector<LongSimplexId> quadOffsets_{};
  };

} // namespace ttk
