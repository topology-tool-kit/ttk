#include <LDistanceMatrix.h>
#include <RipsComplex.h>

#include <limits>

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

int ttk::RipsComplex::computeGaussianDensity(
  double *const density,
  const std::vector<std::vector<double>> &distanceMatrix) const {

  const auto sq = [](const double a) -> double { return a * a; };

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    density[i] = 0.0;
    for(size_t j = 0; j < distanceMatrix.size(); ++j) {
      if(i == j) {
        density[i] += 1.0;
      }
      density[i]
        += std::exp(-sq(distanceMatrix[i][j]) / (2.0 * sq(this->StdDev)));
    }
  }

  return 0;
}

int ttk::RipsComplex::computeDiameterStats(
  const SimplexId nPoints,
  std::array<double *const, 3> diamStats,
  const std::vector<SimplexId> &connectivity,
  const std::vector<double> &cellDiameters) const {

  const auto nCells{cellDiameters.size()};
  if(nCells != connectivity.size() / (this->OutputDimension + 1)) {
    this->printErr("Cell number mismatch");
    return -1;
  }

  std::vector<size_t> nCellsAroundVert(nPoints, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nPoints; ++i) {
    diamStats[0][i] = this->Epsilon; // min
    diamStats[1][i] = 0.0; // mean
    diamStats[2][i] = 0.0; // max
  }

  for(size_t i = 0; i < nCells; ++i) {
    for(int j = 0; j < this->OutputDimension + 1; ++j) {
      const auto p{connectivity[i * (this->OutputDimension + 1) + j]};
      nCellsAroundVert[p]++;
      diamStats[0][p] = std::min(diamStats[0][p], cellDiameters[i]);
      diamStats[1][p] += cellDiameters[i];
      diamStats[2][p] = std::max(diamStats[2][p], cellDiameters[i]);
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nPoints; ++i) {
    if(nCellsAroundVert[i] == 0) {
      diamStats[0][i] = 0.0; // min
      diamStats[1][i] = 0.0; // mean
      diamStats[2][i] = 0.0; // max
    } else {
      diamStats[1][i] /= nCellsAroundVert[i];
    }
  }

  return 0;
}

int ttk::RipsComplex::execute(
  std::vector<SimplexId> &connectivity,
  std::vector<double> &diameters,
  std::array<double *const, 3> diamStats,
  const std::vector<std::vector<double>> &inputMatrix,
  double *const density) const {

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

  this->computeDiameterStats(
    distanceMatrix.size(), diamStats, connectivity, diameters);

  if(this->ComputeGaussianDensity) {
    this->computeGaussianDensity(density, distanceMatrix);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), 1);
  return 0;
}
