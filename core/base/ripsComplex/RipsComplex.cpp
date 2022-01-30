#include <RipsComplex.h>

#include <array>
#include <limits>

ttk::RipsComplex::RipsComplex() {
  this->setDebugMsgPrefix("RipsComplex");
}

template <size_t n>
struct LocCell {
  double diam;
  std::array<ttk::SimplexId, n> verts;
};

void computeEdges(std::vector<ttk::SimplexId> &connectivity,
                  std::vector<double> &diameters,
                  const double epsilon,
                  const std::vector<std::vector<double>> &distanceMatrix,
                  const int nThreads) {

  TTK_FORCE_USE(nThreads);

  std::vector<std::vector<LocCell<1>>> edges(distanceMatrix.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] < epsilon) {
        edges[i].emplace_back(LocCell<1>{
          distanceMatrix[i][j],
          std::array<ttk::SimplexId, 1>{static_cast<ttk::SimplexId>(j)}});
      }
    }
  }

  std::vector<size_t> psum(edges.size() + 1);
  for(size_t i = 0; i < edges.size(); ++i) {
    psum[i + 1] = psum[i] + edges[i].size();
  }

  const auto nCells{psum.back()};
  diameters.resize(nCells);
  connectivity.resize(2 * nCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < edges.size(); ++i) {
    for(size_t j = 0; j < edges[i].size(); ++j) {
      connectivity[2 * (psum[i] + j) + 0] = static_cast<ttk::SimplexId>(i);
      connectivity[2 * (psum[i] + j) + 1] = edges[i][j].verts[0];
      diameters[psum[i] + j] = edges[i][j].diam;
    }
  }
}

static inline void maxAssign(double &a, const double b) {
  if(a < b) {
    a = b;
  }
}

void computeTriangles(std::vector<ttk::SimplexId> &connectivity,
                      std::vector<double> &diameters,
                      const double epsilon,
                      const std::vector<std::vector<double>> &distanceMatrix,
                      const int nThreads) {

  TTK_FORCE_USE(nThreads);

  std::vector<std::vector<LocCell<2>>> triangles(distanceMatrix.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] > epsilon) {
        continue;
      }
      for(size_t k = j + 1; k < distanceMatrix.size(); ++k) {
        if(distanceMatrix[i][k] > epsilon || distanceMatrix[j][k] > epsilon) {
          continue;
        }
        auto diam{distanceMatrix[i][j]};
        maxAssign(diam, distanceMatrix[i][k]);
        maxAssign(diam, distanceMatrix[j][k]);
        triangles[i].emplace_back(
          LocCell<2>{diam, std::array<ttk::SimplexId, 2>{
                             static_cast<ttk::SimplexId>(j),
                             static_cast<ttk::SimplexId>(k),
                           }});
      }
    }
  }

  std::vector<size_t> psum(triangles.size() + 1);
  for(size_t i = 0; i < triangles.size(); ++i) {
    psum[i + 1] = psum[i] + triangles[i].size();
  }

  const auto nCells{psum.back()};
  diameters.resize(nCells);
  connectivity.resize(3 * nCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < triangles.size(); ++i) {
    for(size_t j = 0; j < triangles[i].size(); ++j) {
      connectivity[3 * (psum[i] + j) + 0] = static_cast<ttk::SimplexId>(i);
      connectivity[3 * (psum[i] + j) + 1] = triangles[i][j].verts[0];
      connectivity[3 * (psum[i] + j) + 2] = triangles[i][j].verts[1];
      diameters[psum[i] + j] = triangles[i][j].diam;
    }
  }
}

void computeTetras(std::vector<ttk::SimplexId> &connectivity,
                   std::vector<double> &diameters,
                   const double epsilon,
                   const std::vector<std::vector<double>> &distanceMatrix,
                   const int nThreads) {

  TTK_FORCE_USE(nThreads);

  std::vector<std::vector<LocCell<3>>> tetras(distanceMatrix.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < distanceMatrix.size(); ++i) {
    for(size_t j = i + 1; j < distanceMatrix.size(); ++j) {
      if(distanceMatrix[i][j] > epsilon) {
        continue;
      }
      for(size_t k = j + 1; k < distanceMatrix.size(); ++k) {
        if(distanceMatrix[i][k] > epsilon || distanceMatrix[j][k] > epsilon) {
          continue;
        }
        for(size_t l = k + 1; l < distanceMatrix.size(); ++l) {
          if(distanceMatrix[i][l] > epsilon || distanceMatrix[j][l] > epsilon
             || distanceMatrix[k][l] > epsilon) {
            continue;
          }
          auto diam{distanceMatrix[i][j]};
          maxAssign(diam, distanceMatrix[i][k]);
          maxAssign(diam, distanceMatrix[i][l]);
          maxAssign(diam, distanceMatrix[j][k]);
          maxAssign(diam, distanceMatrix[j][l]);
          maxAssign(diam, distanceMatrix[k][l]);
          tetras[i].emplace_back(
            LocCell<3>{diam, std::array<ttk::SimplexId, 3>{
                               static_cast<ttk::SimplexId>(j),
                               static_cast<ttk::SimplexId>(k),
                               static_cast<ttk::SimplexId>(l),
                             }});
        }
      }
    }
  }

  std::vector<size_t> psum(tetras.size() + 1);
  for(size_t i = 0; i < tetras.size(); ++i) {
    psum[i + 1] = psum[i] + tetras[i].size();
  }

  const auto nCells{psum.back()};
  diameters.resize(nCells);
  connectivity.resize(4 * nCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < tetras.size(); ++i) {
    for(size_t j = 0; j < tetras[i].size(); ++j) {
      connectivity[4 * (psum[i] + j) + 0] = static_cast<ttk::SimplexId>(i);
      connectivity[4 * (psum[i] + j) + 1] = tetras[i][j].verts[0];
      connectivity[4 * (psum[i] + j) + 2] = tetras[i][j].verts[1];
      connectivity[4 * (psum[i] + j) + 3] = tetras[i][j].verts[2];
      diameters[psum[i] + j] = tetras[i][j].diam;
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
  const std::vector<std::vector<double>> &distanceMatrix,
  double *const density) const {

  Timer tm{};

  if(distanceMatrix.empty()
     || distanceMatrix.size() != distanceMatrix[0].size()) {
    this->printErr("Invalid distance matrix");
    return 1;
  }

  Timer tm_rips{};

  if(this->OutputDimension == 1) {
    computeEdges(connectivity, diameters, this->Epsilon, distanceMatrix,
                 this->threadNumber_);
  } else if(this->OutputDimension == 2) {
    computeTriangles(connectivity, diameters, this->Epsilon, distanceMatrix,
                     this->threadNumber_);
  } else if(this->OutputDimension == 3) {
    computeTetras(connectivity, diameters, this->Epsilon, distanceMatrix,
                  this->threadNumber_);
  }

  this->printMsg("Generated Rips complex from distance matrix", 1.0,
                 tm_rips.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  this->computeDiameterStats(
    distanceMatrix.size(), diamStats, connectivity, diameters);

  if(this->ComputeGaussianDensity) {
    this->computeGaussianDensity(density, distanceMatrix);
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);
  return 0;
}
