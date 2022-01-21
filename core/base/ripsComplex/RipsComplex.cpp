#include <LDistanceMatrix.h>
#include <RipsComplex.h>

ttk::RipsComplex::RipsComplex() {
  this->setDebugMsgPrefix("RipsComplex");
}

int ttk::RipsComplex::computeDistanceMatrix(
  std::vector<std::vector<double>> &distanceMatrix,
  const std::vector<std::vector<double>> &inputMatrix) const {

  Timer tm{};

  std::vector<const double *> inputPtrs(inputMatrix.size());
  for(size_t i = 0; i < inputMatrix.size(); ++i) {
    const auto &vec{inputMatrix[i]};
    inputPtrs[i] = vec.data();
  }

  ttk::LDistanceMatrix worker{};
  worker.setThreadNumber(this->threadNumber_);
  worker.setDebugLevel(this->debugLevel_);
  worker.execute(distanceMatrix, inputPtrs, inputMatrix[0].size());

  this->printMsg(
    "Computed distance matrix", 1.0, tm.getElapsedTime(), this->threadNumber_);
  return 0;
}

void computeEdges(std::vector<ttk::LongSimplexId> &connectivity,
                  std::vector<double> &diameters,
                  const double epsilon,
                  const std::vector<std::vector<double>> &distanceMatrix) {
  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
        connectivity.emplace_back(i);
        connectivity.emplace_back(j);
        diameters.emplace_back(distanceMatrix[i][j]);
      }
    }
  }
}

void computeTriangles(std::vector<ttk::LongSimplexId> &connectivity,
                      std::vector<double> &diameters,
                      const double epsilon,
                      const std::vector<std::vector<double>> &distanceMatrix) {
  using ttk::LongSimplexId;

  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    std::vector<size_t> candidates{};
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(candidates.size() < 2) {
      continue;
    }
    std::array<LongSimplexId, 3> verts{static_cast<LongSimplexId>(i), -1, -1};
    for(size_t j = 0; j < candidates.size(); ++j) {
      verts[1] = candidates[j];
      double diam = distanceMatrix[verts[0]][verts[1]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        verts[2] = candidates[k];
        diam = std::max(diam, distanceMatrix[verts[1]][verts[2]]);
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
                   const std::vector<std::vector<double>> &distanceMatrix) {
  using ttk::LongSimplexId;

  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    std::vector<size_t> candidates{};
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
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
      double diam = distanceMatrix[verts[0]][verts[1]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        verts[2] = candidates[k];
        diam = std::max(diam, distanceMatrix[verts[1]][verts[2]]);
        if(diam > epsilon) {
          continue;
        }
        for(size_t l = k + 1; l < candidates.size(); ++l) {
          verts[3] = candidates[l];
          diam = std::max(diam, distanceMatrix[verts[0]][verts[3]]);
          diam = std::max(diam, distanceMatrix[verts[1]][verts[3]]);
          diam = std::max(diam, distanceMatrix[verts[2]][verts[3]]);
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

  std::vector<std::vector<double>> distanceMatrix_{};
  if(!this->InputIsADistanceMatrix) {
    computeDistanceMatrix(distanceMatrix_, inputMatrix);
  }
  const auto &distanceMatrix
    = this->InputIsADistanceMatrix ? inputMatrix : distanceMatrix_;

  if(this->OutputDimension == 1) {
    computeEdges(connectivity, diameters, this->Epsilon, distanceMatrix);
  } else if(this->OutputDimension == 2) {
    computeTriangles(connectivity, diameters, this->Epsilon, distanceMatrix);
  } else if(this->OutputDimension == 3) {
    computeTetras(connectivity, diameters, this->Epsilon, distanceMatrix);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);
  return 0;
}
