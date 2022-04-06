#include <QuadrangulationSubdivision.h>

ttk::SimplexId
  ttk::QuadrangulationSubdivision::findQuadBary(std::vector<float> &sum,
                                                const Quad &quad) const {

  sum.resize(vertexDistance_[quad[0]].size());
  std::fill(sum.begin(), sum.end(), std::numeric_limits<float>::infinity());

  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = false;

    for(const auto vert : quad) {
      if(vertexDistance_[vert][i] == std::numeric_limits<float>::infinity()) {
        skip = true;
        break;
      }
    }

    if(skip) {
      continue;
    }

    const auto &m{vertexDistance_[quad[0]][i]};
    const auto &n{vertexDistance_[quad[1]][i]};
    const auto &o{vertexDistance_[quad[2]][i]};
    const auto &p{vertexDistance_[quad[3]][i]};

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

void ttk::QuadrangulationSubdivision::clearData() {
  outputQuads_.clear();
  outputPoints_.clear();
  outputValences_.clear();
  outputVertType_.clear();
  outputSubdivision_.clear();
  quadNeighbors_.clear();
  vertexDistance_.clear();
}
