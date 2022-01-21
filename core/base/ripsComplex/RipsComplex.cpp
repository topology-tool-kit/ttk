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

void computeEdges(std::vector<ttk::SimplexId> &connectivity,
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

void computeTriangles(std::vector<ttk::SimplexId> &connectivity,
                      std::vector<double> &diameters,
                      const double epsilon,
                      const std::vector<std::vector<double>> &distanceMatrix) {

  std::vector<std::array<ttk::SimplexId, 3>> triangles{};

  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    std::vector<ttk::SimplexId> candidates{};
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(candidates.size() < 2) {
      continue;
    }
    for(size_t j = 0; j < candidates.size(); ++j) {
      double diam = distanceMatrix[i][candidates[j]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        const auto jkd{distanceMatrix[candidates[j]][candidates[k]]};
        if(diam < jkd) {
          diam = jkd;
        }
        if(diam < epsilon) {
          connectivity.emplace_back(static_cast<ttk::SimplexId>(i));
          connectivity.emplace_back(candidates[j]);
          connectivity.emplace_back(candidates[k]);
          diameters.emplace_back(diam);
        }
      }
    }
  }
}

void computeTetras(std::vector<ttk::SimplexId> &connectivity,
                   std::vector<double> &diameters,
                   const double epsilon,
                   const std::vector<std::vector<double>> &distanceMatrix) {

  std::vector<std::array<ttk::SimplexId, 4>> tetras{};

  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    std::vector<ttk::SimplexId> candidates{};
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
        candidates.emplace_back(j);
      }
    }
    if(candidates.size() < 3) {
      continue;
    }
    for(size_t j = 0; j < candidates.size(); ++j) {
      double diam = distanceMatrix[i][candidates[j]];
      for(size_t k = j + 1; k < candidates.size(); ++k) {
        const auto jkd{distanceMatrix[candidates[j]][candidates[k]]};
        if(diam < jkd) {
          diam = jkd;
        }
        if(diam > epsilon) {
          continue;
        }
        for(size_t l = k + 1; l < candidates.size(); ++l) {
          const auto jld{distanceMatrix[candidates[j]][candidates[l]]};
          if(diam < jld) {
            diam = jld;
          }
          const auto kld{distanceMatrix[candidates[k]][candidates[l]]};
          if(diam < kld) {
            diam = kld;
          }
          if(diam > epsilon) {
            continue;
          }
          connectivity.emplace_back(static_cast<ttk::SimplexId>(i));
          connectivity.emplace_back(candidates[j]);
          connectivity.emplace_back(candidates[k]);
          connectivity.emplace_back(candidates[l]);
          diameters.emplace_back(diam);
        }
      }
    }
  }
}

int ttk::RipsComplex::execute(
  std::vector<SimplexId> &connectivity,
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
