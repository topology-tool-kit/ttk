#include <Quadrangulation.h>

#include <OneSkeleton.h>
#include <ZeroSkeleton.h>

ttk::Quadrangulation::Quadrangulation() {
  this->setDebugMsgPrefix("Quadrangulation");
}

int ttk::Quadrangulation::preconditionVertexNeighbors() {

  this->preconditionEdges();

  ZeroSkeleton zsk{};
  zsk.setDebugLevel(this->debugLevel_);
  zsk.setThreadNumber(this->threadNumber_);

  return zsk.buildVertexNeighbors(
    this->nVerts_, this->vertexNeighbors_, this->edges_);
}

int ttk::Quadrangulation::preconditionVertexStars() {

  ZeroSkeleton zsk{};
  zsk.setDebugLevel(this->debugLevel_);
  zsk.setThreadNumber(this->threadNumber_);

  return zsk.buildVertexStars(
    this->nVerts_, this->buildQuadOffsets(), this->vertexStars_);
}

int ttk::Quadrangulation::preconditionEdges() {

  OneSkeleton osk{};
  osk.setDebugLevel(this->debugLevel_);
  osk.setThreadNumber(this->threadNumber_);

  return osk.buildEdgeList(this->nVerts_, this->buildQuadOffsets(),
                           this->edges_, this->edgeStars_, this->quadEdges_);
}

ttk::CellArray ttk::Quadrangulation::buildQuadOffsets() {

  if(static_cast<SimplexId>(this->quadOffsets_.size()) != this->nCells_ + 1) {
    this->quadOffsets_.resize(this->nCells_ + 1);
    this->quadOffsets_[0] = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->nCells_; ++i) {
      this->quadOffsets_[i + 1] = this->cells_[i].size() * (i + 1);
    }
  }

  return CellArray{this->cells_[0].data(), this->quadOffsets_.data(),
                   static_cast<LongSimplexId>(this->nCells_)};
}

void ttk::Quadrangulation::computeStatistics(
  std::vector<SimplexId> &vertsValence,
  std::vector<float> &vertsDensity,
  std::vector<float> &vertsDeformity,
  std::vector<float> &quadArea,
  std::vector<float> &quadDiagsRatio,
  std::vector<float> &quadEdgesRatio,
  std::vector<float> &quadAnglesRatio) const {

  Timer tm;

  vertsValence.resize(this->nVerts_);
  vertsDensity.resize(this->nVerts_);
  vertsDeformity.resize(this->nVerts_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    vertsValence[i] = this->getVertexNeighborNumber(i);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nVerts_; ++i) {
    this->computeDensityAndDeformity(i, vertsDensity[i], vertsDeformity[i]);
  }

  quadArea.resize(this->nCells_);
  quadDiagsRatio.resize(this->nCells_);
  quadEdgesRatio.resize(this->nCells_);
  quadAnglesRatio.resize(this->nCells_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->nCells_; ++i) {
    const auto &q = this->cells_[i];
    const auto &pi = this->vertCoords_[q[0]];
    const auto &pj = this->vertCoords_[q[1]];
    const auto &pk = this->vertCoords_[q[2]];
    const auto &pl = this->vertCoords_[q[3]];

    // quadrangle area
    float area0{}, area1{};
    Geometry::computeTriangleArea(pi.data(), pj.data(), pk.data(), area0);
    Geometry::computeTriangleArea(pi.data(), pk.data(), pl.data(), area1);
    quadArea[i] = area0 + area1;

    // diagonals ratio
    const auto diag0 = Geometry::distance(pi.data(), pk.data());
    const auto diag1 = Geometry::distance(pj.data(), pl.data());
    quadDiagsRatio[i] = std::min(diag0, diag1) / std::max(diag0, diag1);

    // edges ratio
    const std::array<float, 4> edges{
      Geometry::distance(pi.data(), pj.data()), // ij
      Geometry::distance(pj.data(), pk.data()), // jk
      Geometry::distance(pk.data(), pl.data()), // kl
      Geometry::distance(pl.data(), pi.data()), // li
    };
    quadEdgesRatio[i] = *std::min_element(edges.begin(), edges.end())
                        / *std::max_element(edges.begin(), edges.end());

    // angles ratio
    const std::array<float, 4> angles{
      Geometry::angle(pi.data(), pl.data(), pi.data(), pj.data()), // lij
      Geometry::angle(pj.data(), pi.data(), pj.data(), pk.data()), // ijk
      Geometry::angle(pk.data(), pj.data(), pk.data(), pl.data()), // jkl
      Geometry::angle(pl.data(), pk.data(), pl.data(), pi.data()), // kli
    };

    const auto min_max{std::minmax_element(angles.begin(), angles.end())};
    quadAnglesRatio[i] = *min_max.first / *min_max.second;
  }

  // compute ratio between quad area and mean quad area

  // global surface area
  float sumArea{};
  for(const auto a : quadArea) {
    sumArea += a;
  }
  for(auto &a : quadArea) {
    a *= quadArea.size() / sumArea;
  }

  this->printMsg("Computed statistics", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
}
