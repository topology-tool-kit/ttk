#include <RipsComplex.h>

ttk::RipsComplex::RipsComplex() {
  this->setDebugMsgPrefix("RipsComplex");
}

void computeEdges(std::vector<ttk::LongSimplexId> &connectivity,
                  std::vector<double> &diameters,
                  const double epsilon,
                  const std::vector<std::vector<double>> &inputMatrix) {
  for(size_t i = 0; i < inputMatrix.size() - 1; ++i) {
    for(size_t j = i + 1; j < inputMatrix.size(); ++j) {
      if(inputMatrix[i][j] < epsilon) {
        connectivity.emplace_back(i);
        connectivity.emplace_back(j);
        diameters.emplace_back(inputMatrix[i][j]);
      }
    }
  }
}

void computeTriangles(std::vector<ttk::LongSimplexId> &connectivity,
                      std::vector<double> &diameters,
                      const double epsilon,
                      const std::vector<std::vector<double>> &inputMatrix) {
  using ttk::LongSimplexId;

  for(size_t i = 0; i < inputMatrix.size() - 1; ++i) {
    std::vector<size_t> candidates{};
    for(size_t j = i + 1; j < inputMatrix.size(); ++j) {
      if(inputMatrix[i][j] < epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(candidates.size() < 2) {
      continue;
    }
    std::array<LongSimplexId, 3> verts{static_cast<LongSimplexId>(i), -1, -1};
    for(size_t j = 0; j < candidates.size(); ++j) {
      verts[1] = candidates[j];
      double diam = inputMatrix[verts[0]][verts[1]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        verts[2] = candidates[k];
        diam = std::max(diam, inputMatrix[verts[1]][verts[2]]);
        if(diam < epsilon) {
          connectivity.insert(connectivity.end(), verts.begin(), verts.end());
          diameters.emplace_back(diam);
        }
      }
    }
  }
}

void computeTetras(std::vector<ttk::LongSimplexId> &connectivity,
                   std::vector<double> &diameters,
                   const double epsilon,
                   const std::vector<std::vector<double>> &inputMatrix) {
  using ttk::LongSimplexId;

  for(size_t i = 0; i < inputMatrix.size() - 1; ++i) {
    std::vector<size_t> candidates{};
    for(size_t j = i + 1; j < inputMatrix.size(); ++j) {
      if(inputMatrix[i][j] < epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(candidates.size() < 3) {
      continue;
    }
    std::array<LongSimplexId, 4> verts{
      static_cast<LongSimplexId>(i), -1, -1, -1};
    for(size_t j = 0; j < candidates.size(); ++j) {
      verts[1] = candidates[j];
      double diam = inputMatrix[verts[0]][verts[1]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        verts[2] = candidates[k];
        diam = std::max(diam, inputMatrix[verts[1]][verts[2]]);
        if(diam > epsilon) {
          continue;
        }
        for(size_t l = k + 1; l < candidates.size(); ++l) {
          verts[3] = candidates[l];
          diam = std::max(diam, inputMatrix[verts[0]][verts[3]]);
          diam = std::max(diam, inputMatrix[verts[1]][verts[3]]);
          diam = std::max(diam, inputMatrix[verts[2]][verts[3]]);
          if(diam > epsilon) {
            continue;
          }
          connectivity.insert(connectivity.end(), verts.begin(), verts.end());
          diameters.emplace_back(diam);
        }
      }
    }
  }
}

int ttk::RipsComplex::execute(
  std::vector<LongSimplexId> &connectivity,
  std::vector<double> &diameters,
  const std::vector<std::vector<double>> &inputMatrix) const {

  Timer tm{};

  if(this->OutputDimension == 1) {
    computeEdges(connectivity, diameters, this->Epsilon, inputMatrix);
  } else if(this->OutputDimension == 2) {
    computeTriangles(connectivity, diameters, this->Epsilon, inputMatrix);
  } else if(this->OutputDimension == 3) {
    computeTetras(connectivity, diameters, this->Epsilon, inputMatrix);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);
  return 0;
}
