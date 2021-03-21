#include <QuadrangulationSubdivision.h>

ttk::SimplexId ttk::QuadrangulationSubdivision::findQuadBary(
  const std::vector<size_t> &quadVertices) const {

  std::vector<float> sum(vertexDistance_[*quadVertices.begin()].size(),
                         std::numeric_limits<float>::infinity());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = false;

    for(const auto vert : quadVertices) {
      if(vertexDistance_[vert][i] == std::numeric_limits<float>::infinity()) {
        skip = true;
        break;
      }
    }

    if(skip) {
      continue;
    }

    float m = vertexDistance_[quadVertices[0]][i];
    float n = vertexDistance_[quadVertices[1]][i];
    float o = vertexDistance_[quadVertices[2]][i];
    float p = vertexDistance_[quadVertices[3]][i];

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

ttk::QuadrangulationSubdivision::Point
  ttk::QuadrangulationSubdivision::getQuadNormal(const size_t a) const {
  Point quadNormal{};

  // current vertex 3d coordinates
  Point pa = outputPoints_[a];

  // find all quads that have a as vertex
  std::vector<Quad> quads{};
  for(auto &q : outputQuads_) {
    auto _a = static_cast<LongSimplexId>(a);
    if(q[0] == _a || q[1] == _a || q[2] == _a || q[3] == _a) {
      quads.emplace_back(q);
    }
  }

  // store for each quad around a its two direct neighbors
  std::set<std::array<size_t, 2>> couples{};

  // find couple of neighbors of a sharing a quad
  for(auto &q : quads) {
    auto _a = static_cast<LongSimplexId>(a);
    std::array<size_t, 2> tmp{};
    if(q[0] == _a) {
      tmp[0] = static_cast<size_t>(q[3]);
      tmp[1] = static_cast<size_t>(q[1]);
    } else if(q[1] == _a) {
      tmp[0] = static_cast<size_t>(q[0]);
      tmp[1] = static_cast<size_t>(q[2]);
    } else if(q[2] == _a) {
      tmp[0] = static_cast<size_t>(q[1]);
      tmp[1] = static_cast<size_t>(q[3]);
    } else if(q[3] == _a) {
      tmp[0] = static_cast<size_t>(q[2]);
      tmp[1] = static_cast<size_t>(q[0]);
    }
    couples.insert(tmp);
  }

  // triangle normals around current quadrangulation vertex
  std::vector<Point> normals{};

  for(auto &t : couples) {
    Point pb = outputPoints_[t[0]];
    Point pc = outputPoints_[t[1]];

    // triangle normal: cross product of two edges
    Point crossP{};
    // ab, ac vectors
    Point ab = pb - pa;
    Point ac = pc - pa;
    // compute ab ^ ac
    Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
    // magnitude
    auto mag = Geometry::magnitude(&crossP.x);
    // ensure normal not null
    if(mag > powf(10, -FLT_DIG)) {
      // unitary normal vector
      normals.emplace_back(crossP / mag);
    }
  }

  // ensure normals have same direction
  for(size_t i = 1; i < normals.size(); ++i) {
    auto dotProd = Geometry::dotProduct(&normals[0].x, &normals[i].x);
    if(dotProd < 0.0F) {
      normals[i] = normals[i] * -1.0F;
    }
  }

  if(!normals.empty()) {
    // compute mean of normals
    quadNormal = std::accumulate(
      normals.begin(), normals.end(), Point{}, std::plus<Point>());
    quadNormal = quadNormal / static_cast<float>(normals.size());
  } else {
    // set error value directly in output variable...
    quadNormal.x = NAN;
  }

  return quadNormal;
}

int ttk::QuadrangulationSubdivision::getQuadNeighbors(
  const std::vector<Quad> &quads,
  std::vector<std::set<size_t>> &neighbors,
  const bool secondNeighbors) const {
  Timer tm;

  for(auto &q : quads) {
    auto i = static_cast<size_t>(q[0]);
    auto j = static_cast<size_t>(q[1]);
    auto k = static_cast<size_t>(q[2]);
    auto l = static_cast<size_t>(q[3]);
    if(secondNeighbors) {
      neighbors[i].insert(j);
      neighbors[i].insert(k);
      neighbors[i].insert(l);
      neighbors[j].insert(i);
      neighbors[j].insert(k);
      neighbors[j].insert(l);
      neighbors[k].insert(i);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[l].insert(i);
      neighbors[l].insert(j);
      neighbors[l].insert(k);
    } else {
      neighbors[i].insert(j);
      neighbors[i].insert(l);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[j].insert(k);
      neighbors[j].insert(i);
      neighbors[l].insert(k);
      neighbors[l].insert(i);
    }
  }

  this->printMsg("Computed neighbors mapping of "
                   + std::to_string(outputPoints_.size()) + " points",
                 1.0, tm.getElapsedTime(), debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

int ttk::QuadrangulationSubdivision::relax(const std::set<size_t> &filtered) {
  Timer tm;

  // temp storage for relaxed points
  std::vector<Point> tmp(outputPoints_.size());

  // loop over output points, do not touch input MSC critical points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      tmp[i] = outputPoints_[i];
      continue;
    }

    // barycenter of curr neighbors
    Point relax{};
    for(auto &neigh : quadNeighbors_[i]) {
      relax = relax + outputPoints_[neigh];
    }
    relax = relax * (1.0F / static_cast<float>(quadNeighbors_[i].size()));

    tmp[i] = relax;
  }

  outputPoints_ = std::move(tmp);

  this->printMsg(
    "Relaxed " + std::to_string(outputPoints_.size() - inputVertexNumber_)
      + " points",
    1.0, tm.getElapsedTime(), debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0;
}

int ttk::QuadrangulationSubdivision::findExtraordinaryVertices(
  std::set<size_t> &output) const {

  // clear output
  output.clear();

  // hold input quads in a vector
  std::vector<Quad> inputQuads(inputQuadNumber_);

  std::vector<std::set<size_t>> neighbors(inputVertexNumber_);

  // use outputQuads_ here because it contains the input quadrangles before
  // subdivision wrt function call position in execute()
  getQuadNeighbors(outputQuads_, neighbors);

  const size_t NORMAL_VALENCE = 4;

  for(size_t i = 0; i < neighbors.size(); ++i) {
    if(neighbors[i].size() != NORMAL_VALENCE) {
      output.insert(i);
    }
  }

  return 0;
}

void ttk::QuadrangulationSubdivision::clearData() {
  outputQuads_.clear();
  outputPoints_.clear();
  outputValences_.clear();
  outputVertType_.clear();
  outputSubdivision_.clear();
  quadNeighbors_.clear();
  vertexDistance_.clear();
  trianglesChecked_.clear();
  projSucceeded_.clear();
}
