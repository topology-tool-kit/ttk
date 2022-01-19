#include <RipsComplex.h>

ttk::RipsComplex::RipsComplex() {
  this->setDebugMsgPrefix("RipsComplex");
}

int ttk::RipsComplex::execute(
  std::vector<LongSimplexId> &connectivity,
  const std::vector<std::vector<double>> &inputMatrix) const {

  Timer tm{};

  for(size_t i = 0; i < inputMatrix.size() - 1; ++i) {
    std::vector<size_t> candidates{};
    for(size_t j = i + 1; j < inputMatrix.size(); ++j) {
      if(inputMatrix[i][j] < this->Epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(static_cast<SimplexId>(candidates.size()) < this->OutputDimension) {
      continue;
    }

    if(this->OutputDimension == 1) {
      std::array<LongSimplexId, 2> verts{static_cast<LongSimplexId>(i), -1};
      for(const auto c : candidates) {
        verts[1] = c;
        connectivity.insert(connectivity.end(), verts.begin(), verts.end());
      }
    } else if(this->OutputDimension == 2) {
      std::array<LongSimplexId, 3> verts{static_cast<LongSimplexId>(i), -1};
      for(size_t j = 0; j < candidates.size(); ++j) {
        verts[1] = candidates[j];
        for(size_t k = j + 1; k < candidates.size(); ++k) {
          verts[2] = candidates[k];
          connectivity.insert(connectivity.end(), verts.begin(), verts.end());
        }
      }
    } else if(this->OutputDimension == 3) {
      std::array<LongSimplexId, 4> verts{static_cast<LongSimplexId>(i), -1};
      for(size_t j = 0; j < candidates.size(); ++j) {
        verts[1] = candidates[j];
        for(size_t k = j + 1; k < candidates.size(); ++k) {
          verts[2] = candidates[k];
          for(size_t l = k + 1; l < candidates.size(); ++l) {
            verts[3] = candidates[l];
            connectivity.insert(connectivity.end(), verts.begin(), verts.end());
          }
        }
      }
    }
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);
  return 0;
}
