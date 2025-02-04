/// \ingroup baseCode
/// \class ttk::dcg::DiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date November 2016.
///
/// \brief TTK %discreteGradient processing package.
///
/// %DiscreteGradient is a TTK processing package that handles discrete gradient
/// (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation

#pragma once

#include <DiscreteGradient.h>
#include <random>

using ttk::SimplexId;
using ttk::dcg::Cell;
using ttk::dcg::CellExt;
using ttk::dcg::DiscreteGradient;

template <typename dataType, typename triangulationType>
dataType DiscreteGradient::getPersistence(
  const Cell &up,
  const Cell &down,
  const dataType *const scalars,
  const triangulationType &triangulation) const {

  return scalars[getCellGreaterVertex(up, triangulation)]
         - scalars[getCellGreaterVertex(down, triangulation)];
}

template <typename triangulationType>
int DiscreteGradient::buildGradient(const triangulationType &triangulation,
                                    bool bypassCache,
                                    const std::vector<bool> *updateMask) {



  auto &cacheHandler = *triangulation.getGradientCacheHandler();
  const auto findGradient
    = [this, &cacheHandler]() -> AbstractTriangulation::gradientType * {
    if(this->inputScalarField_.first == nullptr) {
      return {};
    }
    return cacheHandler.get(this->inputScalarField_);
  };

#ifdef TTK_ENABLE_OPENMP
  if(!bypassCache && omp_in_parallel()) {
    this->printWrn(
      "buildGradient() called inside a parallel region, disabling cache...");
    bypassCache = true;
  }
#endif // TTK_ENABLE_OPENMP

  // set member variables at each buildGradient() call
  this->dimensionality_ = triangulation.getCellVertexNumber(0) - 1;
  this->numberOfVertices_ = triangulation.getNumberOfVertices();

  this->gradient_ = bypassCache ? &this->localGradient_ : findGradient();
  if(this->gradient_ == nullptr || bypassCache) {

    if(!bypassCache) {
      // add new cache entry
      cacheHandler.insert(this->inputScalarField_, {});
      this->gradient_ = cacheHandler.get(this->inputScalarField_);
    }

    // allocate gradient memory
    this->initMemory(triangulation);
    Timer tm{};
    if(updateMask) {
     this->processLowerStarsWithMask(
     this->inputOffsets_, triangulation, updateMask);
     this->printMsg("Update cached discrete gradient", 1.0,
                tm.getElapsedTime(), this->threadNumber_);
    }else if(this->BackEnd == BACKEND::STOCHASTIC_BACKEND){
      this->processLowerStarsStochastic (this->inputOffsets_, triangulation);
      this->printMsg("Build stochastic discrete gradient", 1.0, tm.getElapsedTime(), 
                    this->threadNumber_);

    }else if(this->BackEnd == BACKEND::CLASSIC_BACKEND){
      this->processLowerStars(this->inputOffsets_, triangulation);
      this->printMsg("Built discrete gradient", 1.0, tm.getElapsedTime(),
                     this->threadNumber_);
    }
  } else {
    this->printMsg("Fetched cached discrete gradient");
    if(updateMask) {
      Timer tm{};
      this->processLowerStarsWithMask(
        this->inputOffsets_, triangulation, updateMask);
      this->printMsg("Update cached discrete gradient", 1.0,
                     tm.getElapsedTime(), this->threadNumber_);
    }
  }
  return 0;
}

template <typename triangulationType>
int DiscreteGradient::setCriticalPoints(
  const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
  std::vector<std::array<float, 3>> &points,
  std::vector<char> &cellDimensions,
  std::vector<SimplexId> &cellIds,
  std::vector<char> &isOnBoundary,
  std::vector<SimplexId> &PLVertexIdentifiers,
  const triangulationType &triangulation) const {

  std::array<size_t, 5> partSums{};
  for(size_t i = 0; i < criticalCellsByDim.size(); ++i) {
    partSums[i + 1] = partSums[i] + criticalCellsByDim[i].size();
  }

  const auto nCritPoints = partSums.back();

  points.resize(nCritPoints);
  cellDimensions.resize(nCritPoints);
  cellIds.resize(nCritPoints);
  isOnBoundary.resize(nCritPoints);
  PLVertexIdentifiers.resize(nCritPoints);

  for(size_t i = 0; i < criticalCellsByDim.size(); ++i) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < criticalCellsByDim[i].size(); ++j) {
      const SimplexId cellId = criticalCellsByDim[i][j];
      const int cellDim = i;
      const auto o{partSums[i] + j};

      triangulation.getCellIncenter(cellId, i, points[o].data());
      cellDimensions[o] = cellDim;
#ifdef TTK_ENABLE_MPI
      ttk::SimplexId globalId{-1};
      triangulation.getDistributedGlobalCellId(cellId, cellDim, globalId);
      cellIds[o] = globalId;
#else
      cellIds[o] = cellId;
#endif // TTK_ENABLE_MPI
      const Cell cell{static_cast<int>(i), cellId};
      isOnBoundary[o] = this->isBoundary(cell, triangulation);
      PLVertexIdentifiers[o] = this->getCellGreaterVertex(cell, triangulation);
    }
  }

  std::vector<std::vector<std::string>> rows(this->dimensionality_ + 1);
  for(int i = 0; i < this->dimensionality_ + 1; ++i) {
    rows[i]
      = std::vector<std::string>{"#" + std::to_string(i) + "-cell(s)",
                                 std::to_string(criticalCellsByDim[i].size())};
  }
  this->printMsg(rows);

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::setCriticalPoints(
  std::vector<std::array<float, 3>> &points,
  std::vector<char> &cellDimensions,
  std::vector<SimplexId> &cellIds,
  std::vector<char> &isOnBoundary,
  std::vector<SimplexId> &PLVertexIdentifiers,
  const triangulationType &triangulation) const {

  std::array<std::vector<SimplexId>, 4> criticalCellsByDim;
  getCriticalPoints(criticalCellsByDim, triangulation);
  setCriticalPoints(criticalCellsByDim, points, cellDimensions, cellIds,
                    isOnBoundary, PLVertexIdentifiers, triangulation);

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::getCriticalPoints(
  std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
  const triangulationType &triangulation) const {

  const auto dims{this->getNumberOfDimensions()};
  for(int i = 0; i < dims; ++i) {

    // map: store critical cell per dimension per thread
    std::vector<std::vector<SimplexId>> critCellsPerThread(this->threadNumber_);
    const auto numberOfCells{this->getNumberOfCells(i, triangulation)};

    // use static scheduling to ensure that critical cells
    // are sorted by id

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) schedule(static)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId j = 0; j < numberOfCells; ++j) {
#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
      if(this->isCellCritical(i, j)) {
        critCellsPerThread[tid].emplace_back(j);
      }
    }

    // reduce: aggregate critical cells per thread
    criticalCellsByDim[i] = std::move(critCellsPerThread[0]);
    for(size_t j = 1; j < critCellsPerThread.size(); ++j) {
      const auto &vec{critCellsPerThread[j]};
      criticalCellsByDim[i].insert(
        criticalCellsByDim[i].end(), vec.begin(), vec.end());
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::getCriticalPoints(
  std::vector<Cell> &criticalPoints,
  const triangulationType &triangulation) const {

  // foreach dimension
  const int numberOfDimensions = getNumberOfDimensions();
  for(int i = 0; i < numberOfDimensions; ++i) {

    // foreach cell of that dimension
    const SimplexId numberOfCells = getNumberOfCells(i, triangulation);
    for(SimplexId j = 0; j < numberOfCells; ++j) {
      const Cell cell(i, j);

      if(isCellCritical(cell)) {
        criticalPoints.push_back(cell);
      }
    }
  }

  return 0;
}

template <typename triangulationType>
SimplexId DiscreteGradient::getNumberOfCells(
  const int dimension, const triangulationType &triangulation) const {

  if(dimension > this->dimensionality_ || dimension < 0) {
    return -1;
  }

  switch(dimension) {
    case 0:
      return triangulation.getNumberOfVertices();
      break;

    case 1:
      return triangulation.getNumberOfEdges();
      break;

    case 2:
      return triangulation.getNumberOfTriangles();
      break;

    case 3:
      return triangulation.getNumberOfCells();
      break;
  }

  return -1;
}

template <typename triangulationType>
inline void DiscreteGradient::lowerStarWithMask(
  lowerStarType &ls,
  const SimplexId a,
  const SimplexId *const offsets,
  const triangulationType &triangulation,
  const std::vector<bool> *updateMask) const {

  // make sure that ls is cleared
  for(auto &vec : ls) {
    vec.clear();
  }

  // a belongs to its lower star
  ls[0].emplace_back(CellExt{0, a});

  if(updateMask != nullptr) {
    for(auto &vec : ls[0]) {
      int vertexId = vec.id_;
      int edgeId = (*gradient_)[0][vertexId];

      (*gradient_)[0][vertexId] = -1;
      if(edgeId != -1) {
        (*gradient_)[1][edgeId] = -1;
      }
    }
  }

  // store lower edges
  const auto nedges = triangulation.getVertexEdgeNumber(a);
  ls[1].reserve(nedges);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    triangulation.getVertexEdge(a, i, edgeId);
    SimplexId vertexId;
    triangulation.getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == a) {
      triangulation.getEdgeVertex(edgeId, 1, vertexId);
    }
    if(offsets[vertexId] < offsets[a]) {
      ls[1].emplace_back(CellExt{1, edgeId, {offsets[vertexId], -1, -1}, {}});
    }
  }

  if(updateMask != nullptr) {
    for(auto &vec : ls[1]) {
      int edgeId = vec.id_;
      int vertexId = (*gradient_)[1][edgeId];
      int triangleId = (*gradient_)[2][edgeId];

      (*gradient_)[1][edgeId] = -1;
      if(vertexId != -1) {
        (*gradient_)[0][vertexId] = -1;
      }

      (*gradient_)[2][edgeId] = -1;
      if(triangleId != -1) {
        (*gradient_)[3][triangleId] = -1;
      }
    }
  }

  if(ls[1].size() < 2) {
    // at least two edges in the lower star for one triangle
    return;
  }

  const auto processTriangle
    = [&](const SimplexId triangleId, const SimplexId v0, const SimplexId v1,
          const SimplexId v2) {
        std::array<SimplexId, 3> lowVerts{-1, -1, -1};
        if(v0 == a) {
          lowVerts[0] = offsets[v1];
          lowVerts[1] = offsets[v2];
        } else if(v1 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v2];
        } else if(v2 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
        }
        // higher order vertex first
        if(lowVerts[0] < lowVerts[1]) {
          std::swap(lowVerts[0], lowVerts[1]);
        }
        if(offsets[a] > lowVerts[0]) { // triangle in lowerStar
          uint8_t j{}, k{};
          // store edges indices of current triangle
          std::array<uint8_t, 3> faces{};
          for(const auto &e : ls[1]) {
            if(e.lowVerts_[0] == lowVerts[0] || e.lowVerts_[0] == lowVerts[1]) {
              faces[k++] = j;
            }
            j++;
          }
          ls[2].emplace_back(CellExt{2, triangleId, lowVerts, faces});
        }
      };

  if(dimensionality_ == 2) {
    // store lower triangles

    // use optimised triangulation methods:
    // getVertexStar instead of getVertexTriangle
    // getCellVertex instead of getTriangleVertex
    const auto ncells = triangulation.getVertexStarNumber(a);
    ls[2].reserve(ncells);
    for(SimplexId i = 0; i < ncells; ++i) {
      SimplexId cellId;
      triangulation.getVertexStar(a, i, cellId);
      SimplexId v0{}, v1{}, v2{};
      triangulation.getCellVertex(cellId, 0, v0);
      triangulation.getCellVertex(cellId, 1, v1);
      triangulation.getCellVertex(cellId, 2, v2);
      processTriangle(cellId, v0, v1, v2);
    }

    if(updateMask != nullptr) {
      for(auto &vec : ls[2]) {
        int triangleId = vec.id_;
        int edgeId = (*gradient_)[3][triangleId];

        (*gradient_)[3][triangleId] = -1;
        if(edgeId != -1) {
          (*gradient_)[2][edgeId] = -1;
        }
      }
    }

  } else if(dimensionality_ == 3) {
    // store lower triangles
    const auto ntri = triangulation.getVertexTriangleNumber(a);
    ls[2].reserve(ntri);
    for(SimplexId i = 0; i < ntri; i++) {
      SimplexId triangleId;
      triangulation.getVertexTriangle(a, i, triangleId);
      SimplexId v0{}, v1{}, v2{};
      triangulation.getTriangleVertex(triangleId, 0, v0);
      triangulation.getTriangleVertex(triangleId, 1, v1);
      triangulation.getTriangleVertex(triangleId, 2, v2);
      processTriangle(triangleId, v0, v1, v2);
    }

    if(updateMask != nullptr) {
      for(auto &vec : ls[2]) {
        int triangleId = vec.id_;
        int edgeId = (*gradient_)[3][triangleId];
        int tetraId = (*gradient_)[4][triangleId];

        (*gradient_)[3][triangleId] = -1;
        if(edgeId != -1) {
          (*gradient_)[2][edgeId] = -1;
        }

        (*gradient_)[4][triangleId] = -1;
        if(tetraId != -1) {
          (*gradient_)[5][tetraId] = -1;
        }
      }
    }

    // at least three triangles in the lower star for one tetra
    if(ls[2].size() >= 3) {
      // store lower tetra
      const auto ncells = triangulation.getVertexStarNumber(a);
      ls[3].reserve(ncells);
      for(SimplexId i = 0; i < ncells; ++i) {
        SimplexId cellId;
        triangulation.getVertexStar(a, i, cellId);
        std::array<SimplexId, 3> lowVerts{-1, -1, -1};
        SimplexId v0{}, v1{}, v2{}, v3{};
        triangulation.getCellVertex(cellId, 0, v0);
        triangulation.getCellVertex(cellId, 1, v1);
        triangulation.getCellVertex(cellId, 2, v2);
        triangulation.getCellVertex(cellId, 3, v3);
        if(v0 == a) {
          lowVerts[0] = offsets[v1];
          lowVerts[1] = offsets[v2];
          lowVerts[2] = offsets[v3];
        } else if(v1 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v2];
          lowVerts[2] = offsets[v3];
        } else if(v2 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
          lowVerts[2] = offsets[v3];
        } else if(v3 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
          lowVerts[2] = offsets[v2];
        }
        if(offsets[a] > *std::max_element(
             lowVerts.begin(), lowVerts.end())) { // tetra in lowerStar

          // higher order vertex first
          std::sort(lowVerts.rbegin(), lowVerts.rend());

          uint8_t j{}, k{};
          // store triangles indices of current tetra
          std::array<uint8_t, 3> faces{};
          for(const auto &t : ls[2]) {
            // lowVerts & t.lowVerts are ordered, no need to check if
            // t.lowVerts[0] == lowVerts[2] or t.lowVerts[1] == lowVerts[0]
            if((t.lowVerts_[0] == lowVerts[0]
                && (t.lowVerts_[1] == lowVerts[1]
                    || t.lowVerts_[1] == lowVerts[2]))
               || (t.lowVerts_[0] == lowVerts[1]
                   && t.lowVerts_[1] == lowVerts[2])) {
              faces[k++] = j;
            }
            j++;
          }

          ls[3].emplace_back(CellExt{3, cellId, lowVerts, faces});
        }
      }

      if(updateMask != nullptr) {
        for(auto &vec : ls[3]) {
          int tetraId = vec.id_;
          int triangleId = (*gradient_)[5][tetraId];
          (*gradient_)[5][tetraId] = -1;
          if(triangleId != -1) {
            (*gradient_)[4][triangleId] = -1;
          }
        }
      }
    }
  }
}

template <typename triangulationType>
inline void
  DiscreteGradient::lowerStar(lowerStarType &ls,
                              const SimplexId a,
                              const SimplexId *const offsets,
                              const triangulationType &triangulation) const {

  // WARNING
  // If you modify this function, please make sure to also report your edit to
  // the other implementation of this function, `lowerStarWithMask`.

  // make sure that ls is cleared
  for(auto &vec : ls) {
    vec.clear();
  }

  // a belongs to its lower star
  CellExt const localCellExt{0, a};
  ls[0].emplace_back(localCellExt);

  // store lower edges
  const auto nedges = triangulation.getVertexEdgeNumber(a);
  ls[1].reserve(nedges);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    triangulation.getVertexEdge(a, i, edgeId);
    SimplexId vertexId;
    triangulation.getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == a) {
      triangulation.getEdgeVertex(edgeId, 1, vertexId);
    }
    if(offsets[vertexId] < offsets[a]) {
      ls[1].emplace_back(CellExt{1, edgeId, {offsets[vertexId], -1, -1}, {}});
    }
  }

  if(ls[1].size() < 2) {
    // at least two edges in the lower star for one triangle
    return;
  }

  const auto processTriangle
    = [&](const SimplexId triangleId, const SimplexId v0, const SimplexId v1,
          const SimplexId v2) {
        std::array<SimplexId, 3> lowVerts{-1, -1, -1};
        if(v0 == a) {
          lowVerts[0] = offsets[v1];
          lowVerts[1] = offsets[v2];
        } else if(v1 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v2];
        } else if(v2 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
        }
        // higher order vertex first
        if(lowVerts[0] < lowVerts[1]) {
          std::swap(lowVerts[0], lowVerts[1]);
        }
        if(offsets[a] > lowVerts[0]) { // triangle in lowerStar
          uint8_t j{}, k{};
          // store edges indices of current triangle
          std::array<uint8_t, 3> faces{};
          for(const auto &e : ls[1]) {
            if(e.lowVerts_[0] == lowVerts[0] || e.lowVerts_[0] == lowVerts[1]) {
              faces[k++] = j;
            }
            j++;
          }
          CellExt const localCellExt2{2, triangleId, lowVerts, faces};
          ls[2].emplace_back(localCellExt2);
        }
      };

  if(dimensionality_ == 2) {
    // store lower triangles

    // use optimised triangulation methods:
    // getVertexStar instead of getVertexTriangle
    // getCellVertex instead of getTriangleVertex
    const auto ncells = triangulation.getVertexStarNumber(a);
    ls[2].reserve(ncells);
    for(SimplexId i = 0; i < ncells; ++i) {
      SimplexId cellId;
      triangulation.getVertexStar(a, i, cellId);
      SimplexId v0{}, v1{}, v2{};
      triangulation.getCellVertex(cellId, 0, v0);
      triangulation.getCellVertex(cellId, 1, v1);
      triangulation.getCellVertex(cellId, 2, v2);
      processTriangle(cellId, v0, v1, v2);
    }
  } else if(dimensionality_ == 3) {
    // store lower triangles
    const auto ntri = triangulation.getVertexTriangleNumber(a);
    ls[2].reserve(ntri);
    for(SimplexId i = 0; i < ntri; i++) {
      SimplexId triangleId;
      triangulation.getVertexTriangle(a, i, triangleId);
      SimplexId v0{}, v1{}, v2{};
      triangulation.getTriangleVertex(triangleId, 0, v0);
      triangulation.getTriangleVertex(triangleId, 1, v1);
      triangulation.getTriangleVertex(triangleId, 2, v2);
      processTriangle(triangleId, v0, v1, v2);
    }

    // at least three triangles in the lower star for one tetra
    if(ls[2].size() >= 3) {
      // store lower tetra
      const auto ncells = triangulation.getVertexStarNumber(a);
      ls[3].reserve(ncells);
      for(SimplexId i = 0; i < ncells; ++i) {
        SimplexId cellId;
        triangulation.getVertexStar(a, i, cellId);
        std::array<SimplexId, 3> lowVerts{-1, -1, -1};
        SimplexId v0{}, v1{}, v2{}, v3{};
        triangulation.getCellVertex(cellId, 0, v0);
        triangulation.getCellVertex(cellId, 1, v1);
        triangulation.getCellVertex(cellId, 2, v2);
        triangulation.getCellVertex(cellId, 3, v3);
        if(v0 == a) {
          lowVerts[0] = offsets[v1];
          lowVerts[1] = offsets[v2];
          lowVerts[2] = offsets[v3];
        } else if(v1 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v2];
          lowVerts[2] = offsets[v3];
        } else if(v2 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
          lowVerts[2] = offsets[v3];
        } else if(v3 == a) {
          lowVerts[0] = offsets[v0];
          lowVerts[1] = offsets[v1];
          lowVerts[2] = offsets[v2];
        }
        if(offsets[a] > *std::max_element(
             lowVerts.begin(), lowVerts.end())) { // tetra in lowerStar

          // higher order vertex first
          std::sort(lowVerts.rbegin(), lowVerts.rend());

          uint8_t j{}, k{};
          // store triangles indices of current tetra
          std::array<uint8_t, 3> faces{};
          for(const auto &t : ls[2]) {
            // lowVerts & t.lowVerts are ordered, no need to check if
            // t.lowVerts[0] == lowVerts[2] or t.lowVerts[1] == lowVerts[0]
            if((t.lowVerts_[0] == lowVerts[0]
                && (t.lowVerts_[1] == lowVerts[1]
                    || t.lowVerts_[1] == lowVerts[2]))
               || (t.lowVerts_[0] == lowVerts[1]
                   && t.lowVerts_[1] == lowVerts[2])) {
              faces[k++] = j;
            }
            j++;
          }

          CellExt const localCellExt3{3, cellId, lowVerts, faces};
          ls[3].emplace_back(localCellExt3);
        }
      }
    }
  }
}

template <typename triangulationType>
inline void DiscreteGradient::pairCells(
  CellExt &alpha, CellExt &beta, const triangulationType &triangulation) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  char localBId{0}, localAId{0};
  SimplexId a{}, b{};

  if(beta.dim_ == 1) {

    for(SimplexId i = 0; i < 2; ++i) {
      triangulation.getEdgeVertex(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto nedges = triangulation.getVertexEdgeNumber(alpha.id_);
    for(SimplexId i = 0; i < nedges; ++i) {
      triangulation.getVertexEdge(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
        break;
      }
    }
  } else if(beta.dim_ == 2) {
    for(SimplexId i = 0; i < 3; ++i) {
      triangulation.getTriangleEdge(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto ntri = triangulation.getEdgeTriangleNumber(alpha.id_);
    for(SimplexId i = 0; i < ntri; ++i) {
      triangulation.getEdgeTriangle(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
        break;
      }
    }
  } else {
    for(SimplexId i = 0; i < 4; ++i) {
      triangulation.getCellTriangle(beta.id_, i, a);
      if(a == alpha.id_) {
        localAId = i;
        break;
      }
    }
    const auto ntetra = triangulation.getTriangleStarNumber(alpha.id_);
    for(SimplexId i = 0; i < ntetra; ++i) {
      triangulation.getTriangleStar(alpha.id_, i, b);
      if(b == beta.id_) {
        localBId = i;
        break;
      }
    }
  }
  (*gradient_)[2 * alpha.dim_][alpha.id_] = localBId;
  (*gradient_)[2 * alpha.dim_ + 1][beta.id_] = localAId;
#else
  TTK_FORCE_USE(triangulation);
  (*gradient_)[2 * alpha.dim_][alpha.id_] = beta.id_;
  (*gradient_)[2 * alpha.dim_ + 1][beta.id_] = alpha.id_;
#endif // TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  alpha.paired_ = true;
  beta.paired_ = true;
}

template <typename triangulationType>
int DiscreteGradient::processLowerStars(
  const SimplexId *const offsets, const triangulationType &triangulation) {

  // WARNING
  // If you modify this function, please make sure to also report your edit to
  // the other implementation of this function, `processLowerStarsWithMask`.

  /* Compute gradient */

  auto nverts = triangulation.getNumberOfVertices();

  // Comparison function for Cells inside priority queues
  const auto orderCells = [&](const CellExt &a, const CellExt &b) -> bool {
    return a.lowVerts_ > b.lowVerts_;
  };

  // Type alias for priority queues
  using pqType
    = std::priority_queue<std::reference_wrapper<CellExt>,
                          std::vector<std::reference_wrapper<CellExt>>,
                          decltype(orderCells)>;

  // To reduce allocations, priority queues and lowerStar objects are
  // cleaned & reused between iterations.

  // Priority queues are pushed at the beginning and popped at the
  // end. To pop the minimum, elements should be sorted in a
  // decreasing order.
  pqType pqZero{orderCells}, pqOne{orderCells};

  // store lower star structure
  lowerStarType Lx;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(Lx, pqZero, pqOne)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId x = 0; x < nverts; x++) {

    // clear priority queues (they should be empty at the end of the
    // previous iteration)
    while(!pqZero.empty()) {
      pqZero.pop();
    }
    while(!pqOne.empty()) {
      pqOne.pop();
    }

    // Insert into pqOne cofacets of cell c_alpha such as numUnpairedFaces == 1
    const auto insertCofacets = [&](const CellExt &ca, lowerStarType &ls) {
      if(ca.dim_ == 1) {
        for(auto &beta : ls[2]) {
          if(ls[1][beta.faces_[0]].id_ == ca.id_
             || ls[1][beta.faces_[1]].id_ == ca.id_) {
            // edge ca belongs to triangle beta
            if(numUnpairedFacesTriangle(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }

      } else if(ca.dim_ == 2) {
        for(auto &beta : ls[3]) {
          if(ls[2][beta.faces_[0]].id_ == ca.id_
             || ls[2][beta.faces_[1]].id_ == ca.id_
             || ls[2][beta.faces_[2]].id_ == ca.id_) {
            // triangle ca belongs to tetra beta
            if(numUnpairedFacesTetra(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }
      }
    };

    lowerStar(Lx, x, offsets, triangulation);
    // In case the vertex is a ghost, the gradient of the
    // simplices of its star is set to GHOST_GRADIENT
#ifdef TTK_ENABLE_MPI
    if(ttk::isRunningWithMPI()
       && triangulation.getVertexRank(x) != ttk::MPIrank_) {
      int sizeDim = Lx.size();
      for(int i = 0; i < sizeDim; i++) {
        int nCells = Lx[i].size();
        for(int j = 0; j < nCells; j++) {
          setCellToGhost(Lx[i][j].dim_, Lx[i][j].id_);
        }
      }
    } else
#endif // TTK_ENABLE_MPI

    {
      // Lx[1] empty => x is a local minimum
      if(!Lx[1].empty()) {
        // get delta: 1-cell (edge) with minimal G value (steeper gradient)
        size_t minId = 0;
        for(size_t i = 1; i < Lx[1].size(); ++i) {
          const auto &a = Lx[1][minId].lowVerts_[0];
          const auto &b = Lx[1][i].lowVerts_[0];
          if(a > b) {
            // edge[i] < edge[0]
            minId = i;
          }
        }

        auto &c_delta = Lx[1][minId];

        // store x (0-cell) -> delta (1-cell) V-path
        pairCells(Lx[0][0], c_delta, triangulation);

        // push every 1-cell in Lx that is not delta into pqZero
        for(auto &alpha : Lx[1]) {
          if(alpha.id_ != c_delta.id_) {
            pqZero.push(alpha);
          }
        }

        // push into pqOne every coface of delta in Lx (2-cells only,
        // 3-cells have not any facet paired yet) such that
        // numUnpairedFaces == 1
        insertCofacets(c_delta, Lx);

        while(!pqOne.empty() || !pqZero.empty()) {
          while(!pqOne.empty()) {
            auto &c_alpha = pqOne.top().get();
            pqOne.pop();
            auto unpairedFaces = numUnpairedFaces(c_alpha, Lx);
            if(unpairedFaces.first == 0) {
              pqZero.push(c_alpha);
            } else {
              auto &c_pair_alpha = Lx[c_alpha.dim_ - 1][unpairedFaces.second];

              // store (pair_alpha) -> (alpha) V-path
              pairCells(c_pair_alpha, c_alpha, triangulation);

              // add cofaces of c_alpha and c_pair_alpha to pqOne
              insertCofacets(c_alpha, Lx);
              insertCofacets(c_pair_alpha, Lx);
            }
          }

          // skip pair_alpha from pqZero:
          // cells in pqZero are not critical if already paired
          while(!pqZero.empty() && pqZero.top().get().paired_) {
            pqZero.pop();
          }

          if(!pqZero.empty()) {
            auto &c_gamma = pqZero.top().get();
            pqZero.pop();

            // gamma is a critical cell
            // mark gamma as paired
            c_gamma.paired_ = true;

            // add cofacets of c_gamma to pqOne
            insertCofacets(c_gamma, Lx);
          }
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::processLowerStarsWithMask(
  const SimplexId *const offsets,
  const triangulationType &triangulation,
  const std::vector<bool> *updateMask) {

  auto nverts = triangulation.getNumberOfVertices();

  /* Compute gradient */

  // Comparison function for Cells inside priority queues
  const auto orderCells = [&](const CellExt &a, const CellExt &b) -> bool {
    return a.lowVerts_ > b.lowVerts_;
  };

  // Type alias for priority queues
  using pqType
    = std::priority_queue<std::reference_wrapper<CellExt>,
                          std::vector<std::reference_wrapper<CellExt>>,
                          decltype(orderCells)>;

  // To reduce allocations, priority queues and lowerStar objects are
  // cleaned & reused between iterations.

  // Priority queues are pushed at the beginning and popped at the
  // end. To pop the minimum, elements should be sorted in a
  // decreasing order.
  pqType pqZero{orderCells}, pqOne{orderCells};

  // store lower star structure
  lowerStarType Lx;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
  firstprivate(Lx, pqZero, pqOne)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId x = 0; x < nverts; x++) {

    if((updateMask != nullptr) && ((*updateMask)[x] == false)) {
      continue;
    }

    // clear priority queues (they should be empty at the end of the
    // previous iteration)
    while(!pqZero.empty()) {
      pqZero.pop();
    }
    while(!pqOne.empty()) {
      pqOne.pop();
    }

    // Insert into pqOne cofacets of cell c_alpha such as numUnpairedFaces == 1
    const auto insertCofacets = [&](const CellExt &ca, lowerStarType &ls) {
      if(ca.dim_ == 1) {
        for(auto &beta : ls[2]) {
          if(ls[1][beta.faces_[0]].id_ == ca.id_
             || ls[1][beta.faces_[1]].id_ == ca.id_) {
            // edge ca belongs to triangle beta
            if(numUnpairedFacesTriangle(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }

      } else if(ca.dim_ == 2) {
        for(auto &beta : ls[3]) {
          if(ls[2][beta.faces_[0]].id_ == ca.id_
             || ls[2][beta.faces_[1]].id_ == ca.id_
             || ls[2][beta.faces_[2]].id_ == ca.id_) {
            // triangle ca belongs to tetra beta
            if(numUnpairedFacesTetra(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }
      }
    };

    lowerStarWithMask(Lx, x, offsets, triangulation, updateMask);
    // In case the vertex is a ghost, the gradient of the
    // simplices of its star is set to GHOST_GRADIENT
#ifdef TTK_ENABLE_MPI
    if(ttk::isRunningWithMPI()
       && triangulation.getVertexRank(x) != ttk::MPIrank_) {
      int sizeDim = Lx.size();
      for(int i = 0; i < sizeDim; i++) {
        int nCells = Lx[i].size();
        for(int j = 0; j < nCells; j++) {
          setCellToGhost(Lx[i][j].dim_, Lx[i][j].id_);
        }
      }
    } else
#endif // TTK_ENABLE_MPI

    {
      // Lx[1] empty => x is a local minimum
      if(!Lx[1].empty()) {
        // get delta: 1-cell (edge) with minimal G value (steeper gradient)
        size_t minId = 0;
        for(size_t i = 1; i < Lx[1].size(); ++i) {
          const auto &a = Lx[1][minId].lowVerts_[0];
          const auto &b = Lx[1][i].lowVerts_[0];
          if(a > b) {
            // edge[i] < edge[0]
            minId = i;
          }
        }

        auto &c_delta = Lx[1][minId];

        // store x (0-cell) -> delta (1-cell) V-path
        pairCells(Lx[0][0], c_delta, triangulation);

        // push every 1-cell in Lx that is not delta into pqZero
        for(auto &alpha : Lx[1]) {
          if(alpha.id_ != c_delta.id_) {
            pqZero.push(alpha);
          }
        }

        // push into pqOne every coface of delta in Lx (2-cells only,
        // 3-cells have not any facet paired yet) such that
        // numUnpairedFaces == 1
        insertCofacets(c_delta, Lx);

        while(!pqOne.empty() || !pqZero.empty()) {
          while(!pqOne.empty()) {
            auto &c_alpha = pqOne.top().get();
            pqOne.pop();
            auto unpairedFaces = numUnpairedFaces(c_alpha, Lx);
            if(unpairedFaces.first == 0) {
              pqZero.push(c_alpha);
            } else {
              auto &c_pair_alpha = Lx[c_alpha.dim_ - 1][unpairedFaces.second];

              // store (pair_alpha) -> (alpha) V-path
              pairCells(c_pair_alpha, c_alpha, triangulation);

              // add cofaces of c_alpha and c_pair_alpha to pqOne
              insertCofacets(c_alpha, Lx);
              insertCofacets(c_pair_alpha, Lx);
            }
          }

          // skip pair_alpha from pqZero:
          // cells in pqZero are not critical if already paired
          while(!pqZero.empty() && pqZero.top().get().paired_) {
            pqZero.pop();
          }

          if(!pqZero.empty()) {
            auto &c_gamma = pqZero.top().get();
            pqZero.pop();

            // gamma is a critical cell
            // mark gamma as paired
            c_gamma.paired_ = true;

            // add cofacets of c_gamma to pqOne
            insertCofacets(c_gamma, Lx);
          }
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::  processLowerStarsStochastic(
  const SimplexId *const offsets, const triangulationType &triangulation) {

  // WARNING
  // If you modify this function, please make sure to also report your edit to
  // the other implementation of this function, `processLowerStarsWithMask`.

  /* Compute gradient */

  auto nverts = triangulation.getNumberOfVertices();

  // Comparison function for Cells inside priority queues
  const auto orderCells = [&](const CellExt &a, const CellExt &b) -> bool {
    return a.lowVerts_ > b.lowVerts_;
  };

  // Type alias for priority queues
  using pqType
    = std::priority_queue<std::reference_wrapper<CellExt>,
                          std::vector<std::reference_wrapper<CellExt>>,
                          decltype(orderCells)>;

  // To reduce allocations, priority queues and lowerStar objects are
  // cleaned & reused between iterations.

  // Priority queues are pushed at the beginning and popped at the
  // end. To pop the minimum, elements should be sorted in a
  // decreasing order.
  pqType pqZero{orderCells}, pqOne{orderCells};

  // store lower star structure
  lowerStarType Lx;

//Compute edge length in directions dx, dy and dz
  float threshold = 10e-12;
  const auto nedges = triangulation.getVertexEdgeNumber(0);
  float x1, x2, x3;
  float y1, y2, y3;
  triangulation.getVertexPoint(0, x1, x2, x3);
  std::vector<float> edgeLengths(3);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    triangulation.getVertexEdge(0, i, edgeId);
    SimplexId vertexId;
    triangulation.getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == 0)
      triangulation.getEdgeVertex(edgeId, 1, vertexId);
    triangulation.getVertexPoint(vertexId, y1, y2, y3);
    if(std::abs(y1-x1)<threshold && std::abs(y2-x2) < threshold){
      edgeLengths[2]=std::abs(x3-y3);
    }
    if(std::abs(y2-x2)<threshold && std::abs(y3-x3) < threshold){
      edgeLengths[0]=std::abs(x1-y1);
    }
    if(std::abs(y1-x1)<threshold && std::abs(y3-x3) < threshold){
      edgeLengths[1]=std::abs(x2-y2);
    }
  }
  

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
firstprivate(Lx, pqZero, pqOne)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId x = 0; x < nverts; x++) {

    // clear priority queues (they should be empty at the end of the
    // previous iteration)
    while(!pqZero.empty()) {
      pqZero.pop();
    }
    while(!pqOne.empty()) {
      pqOne.pop();
    }

    // Insert into pqOne cofacets of cell c_alpha such as numUnpairedFaces == 1
    const auto insertCofacets = [&](const CellExt &ca, lowerStarType &ls) {
      if(ca.dim_ == 1) {
        for(auto &beta : ls[2]) {
          if(ls[1][beta.faces_[0]].id_ == ca.id_
             || ls[1][beta.faces_[1]].id_ == ca.id_) {
            // edge ca belongs to triangle beta
            if(numUnpairedFacesTriangle(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }

      } else if(ca.dim_ == 2) {
        for(auto &beta : ls[3]) {
          if(ls[2][beta.faces_[0]].id_ == ca.id_
             || ls[2][beta.faces_[1]].id_ == ca.id_
             || ls[2][beta.faces_[2]].id_ == ca.id_) {
            // triangle ca belongs to tetra beta
            if(numUnpairedFacesTetra(beta, ls).first == 1) {
              pqOne.push(beta);
            }
          }
        }
      }
    };

    lowerStar(Lx, x, offsets, triangulation);
    // In case the vertex is a ghost, the gradient of the
    // simplices of its star is set to GHOST_GRADIENT
#ifdef TTK_ENABLE_MPI
    if(ttk::isRunningWithMPI()
       && triangulation.getVertexRank(x) != ttk::MPIrank_) {
      int sizeDim = Lx.size();
      for(int i = 0; i < sizeDim; i++) {
        int nCells = Lx[i].size();
        for(int j = 0; j < nCells; j++) {
          setCellToGhost(Lx[i][j].dim_, Lx[i][j].id_);
        }
      }
    } else
#endif // TTK_ENABLE_MPI

    {
      // Lx[1] empty => x is a local minimum
      if(!Lx[1].empty()) {

        size_t minId = 0;
        std::array<float, 3> xCoords;
        triangulation.getVertexPoint(x, xCoords[0], xCoords[1], xCoords[2]);
        //build stencil
        std::vector<SimplexId> stencilIds;  //in order +dx, -dx, +dy, -dy, +dz, -dz 
        std::vector<float> stencilLengths = edgeLengths;
        buildStencil(x, xCoords, triangulation, stencilIds, stencilLengths);

        double grad[3];
        const double* scalars = static_cast<double const *>(this->inputScalarField_.first);
        grad[0]=(scalars[stencilIds[0]] - scalars[stencilIds[1]])/stencilLengths[0];
        grad[1]=(scalars[stencilIds[2]] - scalars[stencilIds[3]])/stencilLengths[1];
        grad[2]= triangulation.getDimensionality()== 3 ? (scalars[stencilIds[4]] - scalars[stencilIds[5]])/stencilLengths[2] : 0;

        std::vector<double> repartitionBounds;
        repartitionBounds.push_back(0);
        std::vector<int> indexInLowerStar;
        indexInLowerStar.push_back(0);
        double totalWeight=0;

        //std::cout<<"scalar products : "<<std::endl;
        for(size_t i = 0  ; i < Lx[1].size(); ++i){
          SimplexId vertexId;
          triangulation.getEdgeVertex(Lx[1][i].id_,0, vertexId);
          if(vertexId == x)triangulation.getEdgeVertex(Lx[1][i].id_ , 1, vertexId);
          if(std::find(stencilIds.begin(), stencilIds.end(), vertexId) == stencilIds.end())continue;
          std::array<float, 3> newCoords;
          triangulation.getVertexPoint(vertexId, newCoords[0], newCoords[1], newCoords[2]);

          float scalarProduct = -(newCoords[0]-xCoords[0])*grad[0] - (newCoords[1]-xCoords[1])*grad[1] - (newCoords[2]-xCoords[2])*grad[2];
          //std::cout<<scalarProduct<<"  ";
          if(scalarProduct > 0){
            double previousBound = repartitionBounds[repartitionBounds.size()-1];
            repartitionBounds.push_back(scalarProduct+previousBound);
            totalWeight+=scalarProduct;
            indexInLowerStar.push_back(i);
          }
        }
        //std::cout<<std::endl;

        for(size_t i = 1  ; i < repartitionBounds.size(); ++i) {
          repartitionBounds[i]/=totalWeight;
        }

        std::random_device rd;        
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        float random_number = dis(gen);
        size_t it=0;

        while(repartitionBounds[it] < random_number && repartitionBounds.size()>1){
          it++;
        }


        minId = indexInLowerStar[it];
        auto &c_delta = Lx[1][minId];

        //std::cout<<"center id : "<<x<<std::endl;
        //float x1{}, x2{}, x3{};
        //triangulation.getVertexPoint(x, x1, x2, x3);
        //std::cout<<"["<<x1<<", "<<x2<<", "<<x3<<"]"<<std::endl;
        //std::cout<<"lowerStar size : "<<Lx[1].size()<<std::endl;
        //std::cout<<"lowerStar points : "<<std::endl;
        //for (size_t z = 0 ; z < Lx[1].size(); z++){
        //  int id;
        //  triangulation.getEdgeVertex(Lx[1][z].id_,0, id);
        //  if(id == x)triangulation.getEdgeVertex(Lx[1][z].id_ , 1, id);
        //  float c1{}, c2{}, c3{};
        //  triangulation.getVertexPoint(id, c1, c2, c3);
        //  std::cout<<id<<" : ["<<c1<<", "<<c2<<", "<<c3<<"]"<<std::endl;
        //}
        //std::cout<<"grad(x) = ["<<grad[0]<<", "<<grad[1]<<", "<<grad[2]<<"]"<<std::endl;
        //std::cout<<"stencils : "<<std::endl;
        //for (size_t z = 0 ; z < stencilIds.size(); z++){
        // std::cout<<stencilIds[z]<<"  ";
        //}
        //std::cout<<std::endl;
        //std::cout<<"weights : ";
        //for (size_t z = 0 ; z < repartitionBounds.size(); z++){
        //  std::cout<<repartitionBounds[z]<<"  ";
        //}
        //std::cout<<std::endl;
        //std::cout<<"randomNumber : "<<random_number<<std::endl;
        //std::cout<<"arrow direction id : "<<Lx[1][minId].id_<<std::endl;
        //int arrowDirId;
        //triangulation.getEdgeVertex(Lx[1][minId].id_,0, arrowDirId);
        //if(arrowDirId == x)triangulation.getEdgeVertex(Lx[1][minId].id_ , 1, arrowDirId);
        //std::cout<<"vector : ["<<x<<", "<<arrowDirId<<"]"<<std::endl;
        //std::cout<<"====================================================="<<std::endl;
    //
    

        // store x (0-cell) -> delta (1-cell) V-path
        pairCells(Lx[0][0], c_delta, triangulation);

        // push every 1-cell in Lx that is not delta into pqZero
        for(auto &alpha : Lx[1]) {
          if(alpha.id_ != c_delta.id_) {
            pqZero.push(alpha);
          }
        }

        // push into pqOne every coface of delta in Lx (2-cells only,
        // 3-cells have not any facet paired yet) such that
        // numUnpairedFaces == 1
        insertCofacets(c_delta, Lx);

        while(!pqOne.empty() || !pqZero.empty()) {
          while(!pqOne.empty()) {
            auto &c_alpha = pqOne.top().get();
            pqOne.pop();
            auto unpairedFaces = numUnpairedFaces(c_alpha, Lx);
            if(unpairedFaces.first == 0) {
              pqZero.push(c_alpha);
            } else {
              auto &c_pair_alpha = Lx[c_alpha.dim_ - 1][unpairedFaces.second];

              // store (pair_alpha) -> (alpha) V-path
              pairCells(c_pair_alpha, c_alpha, triangulation);

              // add cofaces of c_alpha and c_pair_alpha to pqOne
              insertCofacets(c_alpha, Lx);
              insertCofacets(c_pair_alpha, Lx);
            }
          }

          // skip pair_alpha from pqZero:
          // cells in pqZero are not critical if already paired
          while(!pqZero.empty() && pqZero.top().get().paired_) {
            pqZero.pop();
          }

          if(!pqZero.empty()) {
            auto &c_gamma = pqZero.top().get();
            pqZero.pop();

            // gamma is a critical cell
            // mark gamma as paired
            c_gamma.paired_ = true;

            // add cofacets of c_gamma to pqOne
            insertCofacets(c_gamma, Lx);
          }
        }

      }
    }
  }
  return 0;
}

template<typename triangulationType>
void DiscreteGradient::buildStencil(const SimplexId &x,
                                    const std::array<float, 3> xCoords, 
                                    const triangulationType &triangulation,
                                    std::vector<SimplexId> &stencilIds, 
                                    std::vector<float> &stencilLength){
  float threshold = 10e-12;
  const auto nedges = triangulation.getVertexEdgeNumber(x);
  int dimension = triangulation.getDimensionality();
  stencilIds.resize(dimension*2, -1);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    triangulation.getVertexEdge(x, i, edgeId);
    SimplexId vertexId;
    triangulation.getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == x) {
      triangulation.getEdgeVertex(edgeId, 1, vertexId);
    }
    std::array<float, 3> currentCoords;
    triangulation.getVertexPoint(vertexId, currentCoords[0], currentCoords[1], currentCoords[2]);
    if(std::abs(currentCoords[0] - xCoords[0]) < threshold
        && std::abs(currentCoords[1] - xCoords[1]) < threshold && dimension==3){
        if(currentCoords[2] -xCoords[2] > threshold){
         stencilIds[4] = vertexId;
        }else if(currentCoords[2] -xCoords[2] < -threshold){
          stencilIds[5] = vertexId;
        }  
      }else if(std::abs(currentCoords[0] - xCoords[0]) < threshold 
                && std::abs(currentCoords[2] - xCoords[2]) < threshold){
        if(currentCoords[1] -xCoords[1] > threshold){
         stencilIds[2] = vertexId;
        }else if(currentCoords[1] -xCoords[1] < -threshold){
          stencilIds[3] = vertexId;
        }  
      }else if(std::abs(currentCoords[1] - xCoords[1]) < threshold 
                && std::abs(currentCoords[2] - xCoords[2]) < threshold){
        if(currentCoords[0] - xCoords[0] > threshold){
          stencilIds[0] = vertexId;
        }else if(currentCoords[0] - xCoords[0] < -threshold){
          stencilIds[1] = vertexId;
      }
    }
  }
  for (int i = 0 ; i < dimension; i++){
    if(stencilIds[2*i]!=-1 && stencilIds[2*i+1]!=-1){
      stencilLength[i]*=2;
    }
  }
  for (size_t i = 0 ; i < stencilIds.size(); i++){
    if(stencilIds[i]==-1)stencilIds[i]=x;
  }
          
}

template <typename triangulationType>
bool DiscreteGradient::isBoundary(
  const Cell &cell, const triangulationType &triangulation) const {

  if(cell.dim_ > this->dimensionality_ || cell.dim_ < 0) {
    return false;
  }

  const auto vert{this->getCellGreaterVertex(cell, triangulation)};
  return triangulation.isVertexOnBoundary(vert);
}

template <typename triangulationType>
SimplexId
  DiscreteGradient::getPairedCell(const Cell &cell,
                                  const triangulationType &triangulation,
                                  bool isReverse) const {

  // ensure that getPairedCell(Cell, boolean) calls are rejected
  static_assert(
    std::is_base_of<AbstractTriangulation, triangulationType>(),
    "triangulationType should be an AbstractTriangulation derivative");

#ifndef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  TTK_FORCE_USE(triangulation);
#endif // !TTK_ENABLE_DCG_OPTIMIZE_MEMORY

  if((cell.dim_ > this->dimensionality_ - 1 && !isReverse)
     || (cell.dim_ > this->dimensionality_ && isReverse) || cell.dim_ < 0) {
    return -1;
  }

  SimplexId id{-1};

  if(cell.dim_ == 0) {
    if(!isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[0][cell.id_]};
      if(locId != -1) {
        triangulation.getVertexEdge(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[0][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 1) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[1][cell.id_]};
      if(locId != -1) {
        triangulation.getEdgeVertex(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[1][cell.id_];
#endif
    } else {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[2][cell.id_]};
      if(locId != -1) {
        triangulation.getEdgeTriangle(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[2][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 2) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[3][cell.id_]};
      if(locId != -1) {
        triangulation.getTriangleEdge(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[3][cell.id_];
#endif
    } else {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[4][cell.id_]};
      if(locId != -1) {
        triangulation.getTriangleStar(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[4][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 3) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*gradient_)[5][cell.id_]};
      if(locId != -1) {
        triangulation.getCellTriangle(cell.id_, locId, id);
      }
#else
      id = (*gradient_)[5][cell.id_];
#endif
    }
  }

  return id;
}

template <typename triangulationType>
int DiscreteGradient::getDescendingPath(
  const Cell &cell,
  std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  if(cell.dim_ == 0) {
    // assume that cellId is a vertex
    SimplexId currentId = cell.id_;
    SimplexId connectedEdgeId;
    do {
      // add a vertex
      const Cell vertex(0, currentId);
      vpath.push_back(vertex);

      if(isCellCritical(vertex)) {
        break;
      }

      connectedEdgeId = getPairedCell(vertex, triangulation);
      if(connectedEdgeId == -1) {
        break;
      }

      // add an edge
      const Cell edge(1, connectedEdgeId);
      vpath.push_back(edge);

      if(isCellCritical(edge)) {
        break;
      }

      for(int i = 0; i < 2; ++i) {
        SimplexId vertexId;
        triangulation.getEdgeVertex(connectedEdgeId, i, vertexId);

        if(vertexId != currentId) {
          currentId = vertexId;
          break;
        }
      }

    } while(connectedEdgeId != -1);
  }

  return 0;
}

template <typename triangulationType>
bool DiscreteGradient::getDescendingPathThroughWall(
  const Cell &saddle2,
  const Cell &saddle1,
  const std::vector<bool> &isVisited,
  std::vector<Cell> *const vpath,
  const triangulationType &triangulation,
  const bool stopIfMultiConnected,
  const bool enableCycleDetector) const {

  // debug
  const SimplexId numberOfEdges = triangulation.getNumberOfEdges();
  std::vector<bool> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfEdges, false);
  }

  if(dimensionality_ == 3) {
    // add the 2-saddle to the path
    if(vpath != nullptr) {
      vpath->push_back(saddle2);
    }

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        triangulation.getTriangleEdge(saddle2.id_, i, edgeId);
        if(isVisited[edgeId]) {
          // saddle2 can be adjacent to saddle1 on the wall
          if(isSaddle1(Cell(1, edgeId))) {
            if(vpath != nullptr) {
              vpath->push_back(Cell(1, edgeId));
            }
            return false;
          }

          currentId = edgeId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }
    }

    int oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(not isCycle[currentId]) {
          isCycle[currentId] = true;
        } else {
          this->printErr("Cycle detected on the wall of 1-saddle "
                         + std::to_string(saddle1.id_));
          break;
        }
      }

      oldId = currentId;

      // add an edge
      const Cell edge(1, currentId);
      if(vpath != nullptr) {
        vpath->push_back(edge);
      }

      if(isCellCritical(edge)) {
        break;
      }

      const SimplexId connectedTriangleId = getPairedCell(edge, triangulation);

      // add a triangle
      const Cell triangle(2, connectedTriangleId);
      if(vpath != nullptr) {
        vpath->push_back(triangle);
      }

      if(isCellCritical(triangle)) {
        break;
      }

      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        triangulation.getTriangleEdge(connectedTriangleId, i, edgeId);

        if(isVisited[edgeId] and edgeId != oldId) {
          currentId = edgeId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return false;
}

template <typename triangulationType>
int DiscreteGradient::getAscendingPath(const Cell &cell,
                                       std::vector<Cell> &vpath,
                                       const triangulationType &triangulation,
                                       const bool enableCycleDetector) const {

  const SimplexId numberOfCells = triangulation.getNumberOfCells();
  std::vector<bool> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfCells, false);
  }

  if(dimensionality_ == 2) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {
        oldId = currentId;

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) {
          break;
        }

        const SimplexId connectedEdgeId
          = getPairedCell(triangle, triangulation, true);
        if(connectedEdgeId == -1) {
          break;
        }

        // add an edge
        const Cell edge(1, connectedEdgeId);
        vpath.push_back(edge);

        if(isCellCritical(edge)) {
          break;
        }

        const SimplexId starNumber
          = triangulation.getEdgeStarNumber(connectedEdgeId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          triangulation.getEdgeStar(connectedEdgeId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  } else if(dimensionality_ == 3) {
    if(cell.dim_ == 3) {
      // assume that cellId is a tetra
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {

        // debug
        if(enableCycleDetector) {
          if(not isCycle[currentId]) {
            isCycle[currentId] = true;
          } else {
            this->printErr("cycle detected in the path from tetra "
                           + std::to_string(cell.id_));
            break;
          }
        }

        oldId = currentId;

        // add a tetra
        const Cell tetra(3, currentId);
        vpath.push_back(tetra);

        if(isCellCritical(tetra)) {
          break;
        }

        const SimplexId connectedTriangleId
          = getPairedCell(tetra, triangulation, true);
        if(connectedTriangleId == -1) {
          break;
        }

        // add a triangle
        const Cell triangle(2, connectedTriangleId);
        vpath.push_back(triangle);

        if(isCellCritical(triangle)) {
          break;
        }

        const SimplexId starNumber
          = triangulation.getTriangleStarNumber(connectedTriangleId);
        for(SimplexId i = 0; i < starNumber; ++i) {
          SimplexId starId;
          triangulation.getTriangleStar(connectedTriangleId, i, starId);

          if(starId != currentId) {
            currentId = starId;
            break;
          }
        }

        // stop at convergence caused by boundary effect
      } while(currentId != oldId);
    }
  }

  return 0;
}

template <typename triangulationType>
bool DiscreteGradient::getAscendingPathThroughWall(
  const Cell &saddle1,
  const Cell &saddle2,
  const std::vector<bool> &isVisited,
  std::vector<Cell> *const vpath,
  const triangulationType &triangulation,
  const bool stopIfMultiConnected,
  const bool enableCycleDetector,
  bool *const cycleFound) const {

  // debug
  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  std::vector<bool> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfTriangles, false);
  }

  if(dimensionality_ == 3) {
    // add the 1-saddle to the path
    if(vpath != nullptr) {
      vpath->push_back(saddle1);
    }

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      const SimplexId triangleNumber
        = triangulation.getEdgeTriangleNumber(saddle1.id_);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        triangulation.getEdgeTriangle(saddle1.id_, i, triangleId);
        if(isVisited[triangleId]) {
          // saddle1 can be adjacent to saddle2 on the wall
          if(isSaddle2(Cell(2, triangleId))) {
            if(vpath != nullptr) {
              vpath->push_back(Cell(2, triangleId));
            }
            return false;
          }

          currentId = triangleId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }
    }

    if(currentId == -1) {
      return true;
    }

    SimplexId oldId;
    do {

      // debug
      if(enableCycleDetector) {
        if(not isCycle[currentId]) {
          isCycle[currentId] = true;
        } else {
          if(cycleFound)
            *cycleFound = true;
          else
            this->printErr("Cycle detected on the wall of 2-saddle "
                           + std::to_string(saddle2.id_));
          break;
        }
      }

      oldId = currentId;

      // add a triangle
      const Cell triangle(2, currentId);
      if(vpath != nullptr) {
        vpath->push_back(triangle);
      }

      if(isCellCritical(triangle)) {
        break;
      }

      const SimplexId connectedEdgeId
        = getPairedCell(triangle, triangulation, true);

      // add an edge
      const Cell edge(1, connectedEdgeId);
      if(vpath != nullptr) {
        vpath->push_back(edge);
      }

      if(isCellCritical(edge)) {
        break;
      }

      int nconnections = 0;
      const SimplexId triangleNumber
        = triangulation.getEdgeTriangleNumber(connectedEdgeId);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        triangulation.getEdgeTriangle(connectedEdgeId, i, triangleId);

        if(isVisited[triangleId] and triangleId != oldId) {
          currentId = triangleId;
          ++nconnections;
        }
      }
      if(stopIfMultiConnected && nconnections > 1) {
        return true;
      }

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }
  return false;
}

template <typename triangulationType>
bool DiscreteGradient::detectGradientCycle(
  const Cell &cell, const triangulationType &triangulation) const {
  if(dimensionality_ == 3) {
    if(cell.dim_ == 1) {
      const SimplexId originId = getPairedCell(cell, triangulation);
      if(originId == -1)
        return false;

      std::queue<SimplexId> bfs;
      bfs.push(originId);
      std::vector<bool> isVisited(triangulation.getNumberOfTriangles(), false);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId triangleId = bfs.front();
        bfs.pop();

        isVisited[triangleId] = true;

        for(int j = 0; j < 3; ++j) {
          SimplexId edgeId;
          triangulation.getTriangleEdge(triangleId, j, edgeId);

          const SimplexId pairedCellId
            = getPairedCell(Cell(1, edgeId), triangulation);

          if(triangleId == pairedCellId or pairedCellId == -1)
            continue;

          if(isVisited[pairedCellId])
            return true;
          else
            bfs.push(pairedCellId);
        }
      }
    }
  }
  return false;
}

template <typename triangulationType>
int DiscreteGradient::getDescendingWall(
  const Cell &cell,
  VisitedMask &mask,
  const triangulationType &triangulation,
  std::vector<Cell> *const wall,
  std::vector<SimplexId> *const saddles) const {

  if(saddles != nullptr) {
    saddles->clear();
  }

  if(dimensionality_ == 3) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      const SimplexId originId = cell.id_;

      std::queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId triangleId = bfs.front();
        bfs.pop();

        if(!mask.isVisited_[triangleId]) {
          mask.isVisited_[triangleId] = true;
          mask.visitedIds_.emplace_back(triangleId);

          // add the triangle
          if(wall != nullptr) {
            wall->push_back(Cell(2, triangleId));
          }

          for(int j = 0; j < 3; ++j) {
            SimplexId edgeId;
            triangulation.getTriangleEdge(triangleId, j, edgeId);

            if((saddles != nullptr) and isSaddle1(Cell(1, edgeId))) {
              saddles->emplace_back(edgeId);
            }

            const SimplexId pairedCellId
              = getPairedCell(Cell(1, edgeId), triangulation);

            if(pairedCellId != -1 and pairedCellId != triangleId) {
              bfs.push(pairedCellId);
            }
          }
        }
      }

      if(saddles != nullptr && saddles->size() > 1) {
        std::sort(saddles->begin(), saddles->end());
        const auto last = std::unique(saddles->begin(), saddles->end());
        saddles->erase(last, saddles->end());
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::getAscendingWall(
  const Cell &cell,
  VisitedMask &mask,
  const triangulationType &triangulation,
  std::vector<Cell> *const wall,
  std::vector<SimplexId> *const saddles) const {

  if(saddles != nullptr) {
    saddles->clear();
  }

  if(dimensionality_ == 3) {
    if(cell.dim_ == 1) {
      // assume that cellId is an edge
      const SimplexId originId = cell.id_;

      std::queue<SimplexId> bfs;
      bfs.push(originId);

      // BFS traversal
      while(!bfs.empty()) {
        const SimplexId edgeId = bfs.front();
        bfs.pop();

        if(!mask.isVisited_[edgeId]) {
          mask.isVisited_[edgeId] = true;
          mask.visitedIds_.emplace_back(edgeId);

          // add the edge
          if(wall != nullptr) {
            wall->push_back(Cell(1, edgeId));
          }

          const SimplexId triangleNumber
            = triangulation.getEdgeTriangleNumber(edgeId);
          for(SimplexId j = 0; j < triangleNumber; ++j) {
            SimplexId triangleId;
            triangulation.getEdgeTriangle(edgeId, j, triangleId);

            if((saddles != nullptr) and isSaddle2(Cell(2, triangleId))) {
              saddles->emplace_back(triangleId);
            }

            const SimplexId pairedCellId
              = getPairedCell(Cell(2, triangleId), triangulation, true);

            if(pairedCellId != -1 and pairedCellId != edgeId) {
              bfs.push(pairedCellId);
            }
          }
        }
      }

      if(saddles != nullptr && saddles->size() > 1) {
        std::sort(saddles->begin(), saddles->end());
        const auto last = std::unique(saddles->begin(), saddles->end());
        saddles->erase(last, saddles->end());
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::reverseAscendingPath(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  if(dimensionality_ == 2) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        triangulation.getCellEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          (*gradient_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeStarNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeStar(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*gradient_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*gradient_)[3][triangleId] = edgeId;
      (*gradient_)[2][edgeId] = triangleId;
#endif
    }
  } else if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId tetraId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 4; ++k) {
        SimplexId tmp;
        triangulation.getCellTriangle(tetraId, k, tmp);
        if(tmp == triangleId) {
          (*gradient_)[5][tetraId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getTriangleStarNumber(triangleId); ++k) {
        SimplexId tmp;
        triangulation.getTriangleStar(triangleId, k, tmp);
        if(tmp == tetraId) {
          (*gradient_)[4][triangleId] = k;
          break;
        }
      }
#else
      (*gradient_)[5][tetraId] = triangleId;
      (*gradient_)[4][triangleId] = tetraId;
#endif
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::reverseDescendingPath(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  // assume that the first cell is an edge
  for(size_t i = 0; i < vpath.size(); i += 2) {
    const SimplexId edgeId = vpath[i].id_;
    const SimplexId vertId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
    const auto nneighs = triangulation.getVertexEdgeNumber();
    for(int k = 0; k < nneighs; ++k) {
      SimplexId tmp;
      triangulation.getVertexEdge(vertId, k, tmp);
      if(tmp == edgeId) {
        (*gradient_)[0][vertId] = k;
        break;
      }
    }
    const auto nverts = triangulation.getEdgeStarNumber(edgeId);
    for(int k = 0; k < nverts; ++k) {
      SimplexId tmp;
      triangulation.getEdgeVertex(edgeId, k, tmp);
      if(tmp == vertId) {
        (*gradient_)[1][edgeId] = k;
        break;
      }
    }
#else
    TTK_FORCE_USE(triangulation);
    (*gradient_)[0][vertId] = edgeId;
    (*gradient_)[1][edgeId] = vertId;
#endif
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::reverseAscendingPathOnWall(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation,
  bool cancelReversal) const {

  if(dimensionality_ == 3) {
    // assume that the first cell is an edge
    if(cancelReversal) {
      if(vpath.empty())
        return 0;
      (*gradient_)[2][vpath[0].id_] = NULL_GRADIENT;
      // assume that the last cell is a triangle
      if(vpath.size() <= 1)
        return 0;
      (*gradient_)[3][vpath[vpath.size() - 1].id_] = NULL_GRADIENT;
    }
    const SimplexId numberOfCellsInPath = vpath.size();
    const SimplexId startIndex = (cancelReversal ? 2 : 0);
    for(SimplexId i = startIndex; i < numberOfCellsInPath; i += 2) {
      const SimplexId vpathEdgeIndex = i;
      const SimplexId vpathTriangleIndex = (cancelReversal ? i - 1 : i + 1);
      const SimplexId edgeId = vpath[vpathEdgeIndex].id_;
      const SimplexId triangleId = vpath[vpathTriangleIndex].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        triangulation.getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          (*gradient_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeTriangleNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*gradient_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*gradient_)[3][triangleId] = edgeId;
      (*gradient_)[2][edgeId] = triangleId;
#endif
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteGradient::reverseDescendingPathOnWall(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      const SimplexId edgeId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        triangulation.getTriangleEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          (*gradient_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeTriangleNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*gradient_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*gradient_)[2][edgeId] = triangleId;
      (*gradient_)[3][triangleId] = edgeId;
#endif
    }
  }

  return 0;
}

template <typename triangulationType>
ttk::SimplexId DiscreteGradient::getCellGreaterVertex(
  const Cell c, const triangulationType &triangulation) const {

  const auto offsets = this->inputOffsets_;

  auto cellDim = c.dim_;
  auto cellId = c.id_;

  SimplexId vertexId = -1;
  if(cellDim == 0) {
    vertexId = cellId;
  }

  else if(cellDim == 1) {
    SimplexId v0;
    SimplexId v1;
    triangulation.getEdgeVertex(cellId, 0, v0);
    triangulation.getEdgeVertex(cellId, 1, v1);

    if(offsets[v0] > offsets[v1]) {
      vertexId = v0;
    } else {
      vertexId = v1;
    }
  }

  else if(cellDim == 2) {
    SimplexId v0{}, v1{}, v2{};
    triangulation.getTriangleVertex(cellId, 0, v0);
    triangulation.getTriangleVertex(cellId, 1, v1);
    triangulation.getTriangleVertex(cellId, 2, v2);
    if(offsets[v0] > offsets[v1] && offsets[v0] > offsets[v2]) {
      vertexId = v0;
    } else if(offsets[v1] > offsets[v0] && offsets[v1] > offsets[v2]) {
      vertexId = v1;
    } else {
      vertexId = v2;
    }
  }

  else if(cellDim == 3) {
    SimplexId v0{}, v1{}, v2{}, v3{};
    triangulation.getCellVertex(cellId, 0, v0);
    triangulation.getCellVertex(cellId, 1, v1);
    triangulation.getCellVertex(cellId, 2, v2);
    triangulation.getCellVertex(cellId, 3, v3);
    if(offsets[v0] > offsets[v1] && offsets[v0] > offsets[v2]
       && offsets[v0] > offsets[v3]) {
      vertexId = v0;
    } else if(offsets[v1] > offsets[v0] && offsets[v1] > offsets[v2]
              && offsets[v1] > offsets[v3]) {
      vertexId = v1;
    } else if(offsets[v2] > offsets[v0] && offsets[v2] > offsets[v1]
              && offsets[v2] > offsets[v3]) {
      vertexId = v2;
    } else {
      vertexId = v3;
    }
  }
  return vertexId;
}

template <typename triangulationType>
ttk::SimplexId DiscreteGradient::getCellLowerVertex(
  const Cell c, const triangulationType &triangulation) const {

  const auto offsets = this->inputOffsets_;

  auto cellDim = c.dim_;
  auto cellId = c.id_;

  SimplexId vertexId = -1;
  if(cellDim == 0) {
    vertexId = cellId;
  }

  else if(cellDim == 1) {
    SimplexId v0;
    SimplexId v1;
    triangulation.getEdgeVertex(cellId, 0, v0);
    triangulation.getEdgeVertex(cellId, 1, v1);

    if(offsets[v0] < offsets[v1]) {
      vertexId = v0;
    } else {
      vertexId = v1;
    }
  }

  else if(cellDim == 2) {
    SimplexId v0{}, v1{}, v2{};
    triangulation.getTriangleVertex(cellId, 0, v0);
    triangulation.getTriangleVertex(cellId, 1, v1);
    triangulation.getTriangleVertex(cellId, 2, v2);
    if(offsets[v0] < offsets[v1] && offsets[v0] < offsets[v2]) {
      vertexId = v0;
    } else if(offsets[v1] < offsets[v0] && offsets[v1] < offsets[v2]) {
      vertexId = v1;
    } else {
      vertexId = v2;
    }
  }

  else if(cellDim == 3) {
    SimplexId v0{}, v1{}, v2{}, v3{};
    triangulation.getCellVertex(cellId, 0, v0);
    triangulation.getCellVertex(cellId, 1, v1);
    triangulation.getCellVertex(cellId, 2, v2);
    triangulation.getCellVertex(cellId, 3, v3);
    if(offsets[v0] < offsets[v1] && offsets[v0] < offsets[v2]
       && offsets[v0] < offsets[v3]) {
      vertexId = v0;
    } else if(offsets[v1] < offsets[v0] && offsets[v1] < offsets[v2]
              && offsets[v1] < offsets[v3]) {
      vertexId = v1;
    } else if(offsets[v2] < offsets[v0] && offsets[v2] < offsets[v1]
              && offsets[v2] < offsets[v3]) {
      vertexId = v2;
    } else {
      vertexId = v3;
    }
  }
  return vertexId;
}

template <typename triangulationType>
int DiscreteGradient::setGradientGlyphs(
  std::vector<std::array<float, 3>> &points,
  std::vector<char> &points_pairOrigins,
  std::vector<char> &cells_pairTypes,
  std::vector<SimplexId> &cellIds,
  std::vector<char> &cellDimensions,
  const triangulationType &triangulation) const {

  const auto nDims = this->getNumberOfDimensions();

  // number of glyphs per dimension
  std::vector<size_t> nGlyphsPerDim(nDims);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < nDims - 1; ++i) {
    const auto nCells = this->getNumberOfCells(i, triangulation);
    for(SimplexId j = 0; j < nCells; ++j) {
      if(this->getPairedCell(Cell{i, j}, triangulation) > -1) {
        nGlyphsPerDim[i]++;
      }
    }
  }

  // partial sum of number of gradient glyphs
  std::vector<size_t> offsets(nDims + 1);
  for(SimplexId i = 0; i < nDims; ++i) {
    offsets[i + 1] = offsets[i] + nGlyphsPerDim[i];
  }

  // total number of glyphs
  const auto nGlyphs = offsets.back();

  // resize arrays accordingly
  points.resize(2 * nGlyphs);
  points_pairOrigins.resize(2 * nGlyphs);
  cells_pairTypes.resize(nGlyphs);
  cellIds.resize(2 * nGlyphs);
  cellDimensions.resize(2 * nGlyphs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < nDims - 1; ++i) {
    const SimplexId nCells = getNumberOfCells(i, triangulation);
    size_t nProcessedGlyphs{offsets[i]};
    for(SimplexId j = 0; j < nCells; ++j) {
      const Cell c{i, j};
      const auto pcid = this->getPairedCell(c, triangulation);
      if(pcid > -1) {
        const Cell pc{i + 1, pcid};
        triangulation.getCellIncenter(
          c.id_, c.dim_, points[2 * nProcessedGlyphs].data());
        triangulation.getCellIncenter(
          pc.id_, pc.dim_, points[2 * nProcessedGlyphs + 1].data());
        points_pairOrigins[2 * nProcessedGlyphs] = 0;
        points_pairOrigins[2 * nProcessedGlyphs + 1] = 1;
        cells_pairTypes[nProcessedGlyphs] = i;
#ifdef TTK_ENABLE_MPI
        ttk::SimplexId globalId{-1};
        triangulation.getDistributedGlobalCellId(j, i, globalId);
        cellIds[2 * nProcessedGlyphs + 0] = globalId;
        triangulation.getDistributedGlobalCellId(pcid, i + 1, globalId);
        cellIds[2 * nProcessedGlyphs + 1] = globalId;
#else
        cellIds[2 * nProcessedGlyphs + 0] = j;
        cellIds[2 * nProcessedGlyphs + 1] = pcid;
#endif // TTK_ENABLE_MPI
        cellDimensions[2 * nProcessedGlyphs + 0] = i;
        cellDimensions[2 * nProcessedGlyphs + 1] = i + 1;
        nProcessedGlyphs++;
      }
    }
  }

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::simplifySaddleSaddleConnections1(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::vector<char> &isPL,
  const int iterationThreshold,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors,
  const triangulationType &triangulation) {
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<SimplexId> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle1, pl2dmt_saddle1,
                                 triangulation);

  std::vector<char> isRemovableSaddle2;
  std::vector<SimplexId> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle2, pl2dmt_saddle2,
                                 triangulation);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<SimplexId> saddle1Index;
  std::vector<SimplexId> saddle2Index;
  initializeSaddleSaddleConnections1<dataType>(
    isRemovableSaddle1, isRemovableSaddle2, allowBruteForce, vpaths,
    dmt_criticalPoints, saddle1Index, saddle2Index, triangulation);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>>
    S(cmp_f);
  orderSaddleSaddleConnections1<dataType>(vpaths, dmt_criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections1<dataType>(
    iterationThreshold, isPL, allowBoundary, allowBruteForce,
    returnSaddleConnectors, S, pl2dmt_saddle1, pl2dmt_saddle2,
    isRemovableSaddle1, isRemovableSaddle2, vpaths, dmt_criticalPoints,
    saddle1Index, saddle2Index, triangulation);

  this->printMsg("Saddle-Saddle pairs simplified #1", 1.0, t.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int DiscreteGradient::simplifySaddleSaddleConnections2(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::vector<char> &isPL,
  const int iterationThreshold,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors,
  const triangulationType &triangulation) {
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<SimplexId> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle1, pl2dmt_saddle1,
                                 triangulation);

  std::vector<char> isRemovableSaddle2;
  std::vector<SimplexId> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle2, pl2dmt_saddle2,
                                 triangulation);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<SimplexId> saddle1Index;
  std::vector<SimplexId> saddle2Index;
  initializeSaddleSaddleConnections2<dataType>(
    isRemovableSaddle1, isRemovableSaddle2, allowBruteForce, vpaths,
    dmt_criticalPoints, saddle1Index, saddle2Index, triangulation);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>>
    S(cmp_f);
  orderSaddleSaddleConnections2<dataType>(vpaths, dmt_criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections2<dataType>(
    iterationThreshold, isPL, allowBoundary, allowBruteForce,
    returnSaddleConnectors, S, pl2dmt_saddle1, pl2dmt_saddle2,
    isRemovableSaddle1, isRemovableSaddle2, vpaths, dmt_criticalPoints,
    saddle1Index, saddle2Index, triangulation);

  this->printMsg("Saddle-Saddle pairs simplified #2", 1.0, t.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::filterSaddleConnectors(
  const bool allowBoundary, const triangulationType &triangulation) {

  const bool allowBruteForce = false;
  const bool returnSaddleConnectors = true;

  // get the node type of a contour tree node (for compatibility with
  // ScalarFieldCriticalPoints)
  auto getNodeType = [&](const ftm::Node *node) {
    const int upDegree = node->getNumberOfUpSuperArcs();
    const int downDegree = node->getNumberOfDownSuperArcs();
    const int degree = upDegree + downDegree;

    // saddle point
    if(degree > 1) {
      if(upDegree == 2 and downDegree == 1) {
        return 2;
      } else if(upDegree == 1 and downDegree == 2) {
        return 1;
      }
    }
    // local extremum
    else {
      if(upDegree) {
        return 0;
      } else {
        return 3;
      }
    }

    return -1;
  };

  std::vector<std::pair<SimplexId, char>> cpset;

  const auto *const scalars
    = static_cast<const dataType *>(inputScalarField_.first);
  const auto *const offsets = inputOffsets_;

  ftm::FTMTree contourTree_{};
  contourTree_.setDebugLevel(debugLevel_);
  contourTree_.setVertexScalars(scalars);
  contourTree_.setTreeType(ftm::TreeType::Contour);
  contourTree_.setVertexSoSoffsets(offsets);
  contourTree_.setThreadNumber(threadNumber_);
  contourTree_.setSegmentation(false);
  contourTree_.build<dataType>(&triangulation);
  ftm::FTMTree_MT *tree = contourTree_.getTree(ftm::TreeType::Contour);

  const SimplexId numberOfNodes = tree->getNumberOfNodes();
  for(SimplexId nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
    const ftm::Node *node = tree->getNode(nodeId);
    const SimplexId vertexId = node->getVertexId();

    cpset.push_back(std::make_pair(vertexId, getNodeType(node)));
  }

  std::vector<char> isPL;
  getCriticalPointMap(cpset, isPL);

  simplifySaddleSaddleConnections1<dataType>(
    cpset, isPL, IterationThreshold, allowBoundary, allowBruteForce,
    returnSaddleConnectors, triangulation);
  simplifySaddleSaddleConnections2<dataType>(
    cpset, isPL, IterationThreshold, allowBoundary, allowBruteForce,
    returnSaddleConnectors, triangulation);

  return 0;
}

template <typename dataType, typename triangulationType>
int DiscreteGradient::getRemovableSaddles1(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const bool allowBoundary,
  std::vector<char> &isRemovableSaddle,
  std::vector<SimplexId> &pl2dmt_saddle,
  const triangulationType &triangulation) {

  const SimplexId numberOfEdges = triangulation.getNumberOfEdges();
  isRemovableSaddle.resize(numberOfEdges);

  std::vector<int>dmt1Saddle2PL_(numberOfEdges, -1);
  //dmt1Saddle2PL_.resize(numberOfEdges);
  //std::fill(dmt1Saddle2PL_.begin(), dmt1Saddle2PL_.end(), -1);

  // by default : 1-saddle is removable
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfEdges; ++i) {
    const Cell saddleCandidate(1, i);
    isRemovableSaddle[i] = isSaddle1(saddleCandidate);
  }

  // is [edgeId] in star of PL-1saddle?
  for(auto &criticalPoint : criticalPoints) {
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    if(criticalPointType == static_cast<char>(CriticalType::Saddle1)) {
      if(!allowBoundary and triangulation.isVertexOnBoundary(criticalPointId)) {
        continue;
      }

      SimplexId numberOfSaddles = 0;
      SimplexId saddleId = -1;
      const SimplexId edgeNumber
        = triangulation.getVertexEdgeNumber(criticalPointId);
      for(SimplexId i = 0; i < edgeNumber; ++i) {
        SimplexId edgeId;
        triangulation.getVertexEdge(criticalPointId, i, edgeId);
        const Cell saddleCandidate(1, edgeId);

        if(isSaddle1(saddleCandidate) and dmt1Saddle2PL_[edgeId] == -1) {
          saddleId = edgeId;
          ++numberOfSaddles;
        }
      }

      // only one DMT-1saddle in the star so this one is non-removable
      if(numberOfSaddles == 1) {
        if(dmt1Saddle2PL_[saddleId] == -1
           and pl2dmt_saddle[criticalPointId] == -1) {
          dmt1Saddle2PL_[saddleId] = criticalPointId;
          pl2dmt_saddle[criticalPointId] = saddleId;
          isRemovableSaddle[saddleId] = false;
        }
      }
    }
  }

  return 0;
}

template <typename dataType, typename triangulationType>
int DiscreteGradient::getRemovableSaddles2(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const bool allowBoundary,
  std::vector<char> &isRemovableSaddle,
  std::vector<SimplexId> &pl2dmt_saddle,
  const triangulationType &triangulation) {

  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  isRemovableSaddle.resize(numberOfTriangles);

  std::vector<int> dmt2Saddle2PL_(numberOfTriangles, -1);
  //dmt2Saddle2PL_.resize(numberOfTriangles);
  //std::fill(dmt2Saddle2PL_.begin(), dmt2Saddle2PL_.end(), -1);

  // by default : 2-saddle is removable
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfTriangles; ++i) {
    const Cell saddleCandidate(2, i);
    isRemovableSaddle[i] = isSaddle2(saddleCandidate);
  }

  // is [triangleId] in star of PL-2saddle?
  for(auto &criticalPoint : criticalPoints) {
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    if(criticalPointType == static_cast<char>(CriticalType::Saddle2)) {
      if(!allowBoundary and triangulation.isVertexOnBoundary(criticalPointId)) {
        continue;
      }

      SimplexId numberOfSaddles = 0;
      SimplexId saddleId = -1;
      const SimplexId triangleNumber
        = triangulation.getVertexTriangleNumber(criticalPointId);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        triangulation.getVertexTriangle(criticalPointId, i, triangleId);
        const Cell saddleCandidate(2, triangleId);

        if(isSaddle2(saddleCandidate) and dmt2Saddle2PL_[triangleId] == -1) {
          saddleId = triangleId;
          ++numberOfSaddles;
        }
      }

      // only one DMT-2saddle in the star so this one is non-removable
      if(numberOfSaddles == 1) {
        if(dmt2Saddle2PL_[saddleId] == -1
           and pl2dmt_saddle[criticalPointId] == -1) {
          dmt2Saddle2PL_[saddleId] = criticalPointId;
          pl2dmt_saddle[criticalPointId] = saddleId;
          isRemovableSaddle[saddleId] = false;
        }
      }
    }
  }

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::initializeSaddleSaddleConnections1(
  const std::vector<char> &isRemovableSaddle1,
  const std::vector<char> &isRemovableSaddle2,
  const bool allowBruteForce,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index,
  const triangulationType &triangulation) const {
  Timer t;

  const auto *const scalars
    = static_cast<const dataType *>(inputScalarField_.first);

  const int maximumDim = dimensionality_;
  const int saddle2Dim = maximumDim - 1;
  const int saddle1Dim = saddle2Dim - 1;

  // Part 1 : build initial structures
  // add the 2-saddles to CriticalPointList
  const SimplexId numberOfSaddle2Candidates
    = getNumberOfCells(saddle2Dim, triangulation);
  saddle2Index.resize(numberOfSaddle2Candidates, -1);
  for(SimplexId i = 0; i < numberOfSaddle2Candidates; ++i) {
    if(allowBruteForce or isRemovableSaddle2[i]) {
      const Cell saddle2Candidate(saddle2Dim, i);

      if(isSaddle2(saddle2Candidate)) {
        const SimplexId index = criticalPoints.size();
        saddle2Index[i] = index;
        criticalPoints.emplace_back(saddle2Candidate);
      }
    }
  }
  const SimplexId numberOf2Saddles = criticalPoints.size();

  // add the 1-saddles to CriticalPointList
  const SimplexId numberOfSaddle1Candidates
    = getNumberOfCells(saddle1Dim, triangulation);
  saddle1Index.resize(numberOfSaddle1Candidates, -1);
  for(SimplexId i = 0; i < numberOfSaddle1Candidates; ++i) {
    if(isRemovableSaddle1[i]) {
      const Cell saddle1Candidate(saddle1Dim, i);

      const SimplexId index = criticalPoints.size();
      saddle1Index[i] = index;
      criticalPoints.emplace_back(saddle1Candidate);
    }
  }

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
  std::vector<bool> isVisited(numberOfSaddle2Candidates, false);
  std::vector<SimplexId> visitedTriangles{};

  for(SimplexId i = 0; i < numberOf2Saddles; ++i) {
    const SimplexId destinationIndex = i;
    CriticalPoint &destination = criticalPoints[destinationIndex];
    const Cell &saddle2 = destination.cell_;
    VisitedMask mask{isVisited, visitedTriangles};

    std::vector<SimplexId> saddles1;
    getDescendingWall(saddle2, mask, triangulation, nullptr, &saddles1);

    for(auto &saddle1Id : saddles1) {
      if(!isRemovableSaddle1[saddle1Id]) {
        continue;
      }

      const Cell &saddle1 = Cell(1, saddle1Id);

      std::vector<Cell> path;
      const bool isMultiConnected = getAscendingPathThroughWall(
        saddle1, saddle2, isVisited, &path, triangulation, true);

      if(!isMultiConnected) {
        const SimplexId sourceIndex = saddle1Index[saddle1Id];
        CriticalPoint &source = criticalPoints[sourceIndex];

        // update source and destination
        const SimplexId sourceSlot = source.addSlot();
        const SimplexId destinationSlot = destination.addSlot();

        // update vpath
        const auto persistence
          = getPersistence(saddle2, saddle1, scalars, triangulation);

        vpaths.push_back(VPath(true, -1, sourceIndex, destinationIndex,
                               sourceSlot, destinationSlot, persistence));
      }
    }
  }

  // Part 3 : initialize the last structures
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    CriticalPoint &cp = criticalPoints[i];

    const SimplexId numberOfSlots = cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_ = 0;
  }

  const SimplexId numberOfVPaths = vpaths.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfVPaths; ++i) {
    const VPath &vpath = vpaths[i];

    if(vpath.isValid_) {
      const SimplexId sourceIndex = vpath.source_;
      const SimplexId destinationIndex = vpath.destination_;

      const SimplexId sourceSlot = vpath.sourceSlot_;
      const SimplexId destinationSlot = vpath.destinationSlot_;

      CriticalPoint &source = criticalPoints[sourceIndex];
      CriticalPoint &destination = criticalPoints[destinationIndex];

      source.vpaths_[sourceSlot] = i;
      destination.vpaths_[destinationSlot] = i;
    }
  }

  this->printMsg(
    " Initialization step #1", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::initializeSaddleSaddleConnections2(
  const std::vector<char> &isRemovableSaddle1,
  const std::vector<char> &isRemovableSaddle2,
  const bool allowBruteForce,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index,
  const triangulationType &triangulation) const {
  Timer t;

  const auto *const scalars
    = static_cast<const dataType *>(inputScalarField_.first);

  const int maximumDim = dimensionality_;
  const int saddle2Dim = maximumDim - 1;
  const int saddle1Dim = saddle2Dim - 1;

  // Part 1 : build initial structures
  // add the 1-saddles to CriticalPointList
  const SimplexId numberOfSaddle1Candidates
    = getNumberOfCells(saddle1Dim, triangulation);
  saddle1Index.resize(numberOfSaddle1Candidates, -1);
  for(SimplexId i = 0; i < numberOfSaddle1Candidates; ++i) {
    if(isRemovableSaddle1[i]) {
      const Cell saddle1Candidate(saddle1Dim, i);

      const SimplexId index = criticalPoints.size();
      saddle1Index[i] = index;
      criticalPoints.emplace_back(saddle1Candidate);
    }
  }
  const SimplexId numberOf1Saddles = criticalPoints.size();

  // add the 2-saddles to CriticalPointList
  const SimplexId numberOfSaddle2Candidates
    = getNumberOfCells(saddle2Dim, triangulation);
  saddle2Index.resize(numberOfSaddle2Candidates, -1);
  for(SimplexId i = 0; i < numberOfSaddle2Candidates; ++i) {
    if(allowBruteForce or isRemovableSaddle2[i]) {
      const Cell saddle2Candidate(saddle2Dim, i);

      if(isSaddle2(saddle2Candidate)) {
        const SimplexId index = criticalPoints.size();
        saddle2Index[i] = index;
        criticalPoints.emplace_back(saddle2Candidate);
      }
    }
  }

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
  std::vector<bool> isVisited(numberOfSaddle1Candidates, false);
  std::vector<SimplexId> visitedEdges{};
  for(SimplexId i = 0; i < numberOf1Saddles; ++i) {
    const SimplexId sourceIndex = i;
    CriticalPoint &source = criticalPoints[sourceIndex];
    const Cell &saddle1 = source.cell_;

    std::vector<SimplexId> saddles2;
    VisitedMask mask{isVisited, visitedEdges};
    getAscendingWall(saddle1, mask, triangulation, nullptr, &saddles2);

    for(auto &saddle2Id : saddles2) {
      if(!isRemovableSaddle2[saddle2Id]) {
        continue;
      }

      const Cell &saddle2 = Cell(2, saddle2Id);

      std::vector<Cell> path;
      const bool isMultiConnected = getDescendingPathThroughWall(
        saddle2, saddle1, isVisited, &path, triangulation, true);

      if(!isMultiConnected) {
        const SimplexId destinationIndex = saddle2Index[saddle2Id];
        CriticalPoint &destination = criticalPoints[destinationIndex];

        // update source and destination
        const SimplexId sourceSlot = source.addSlot();
        const SimplexId destinationSlot = destination.addSlot();

        // update vpath
        const auto persistence
          = getPersistence(saddle2, saddle1, scalars, triangulation);

        vpaths.push_back(VPath(true, -1, sourceIndex, destinationIndex,
                               sourceSlot, destinationSlot, persistence));
      }
    }
  }

  // Part 3 : initialize the last structures
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    CriticalPoint &cp = criticalPoints[i];

    const SimplexId numberOfSlots = cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_ = 0;
  }

  const SimplexId numberOfVPaths = vpaths.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfVPaths; ++i) {
    const VPath &vpath = vpaths[i];

    if(vpath.isValid_) {
      const SimplexId sourceIndex = vpath.source_;
      const SimplexId destinationIndex = vpath.destination_;

      const SimplexId sourceSlot = vpath.sourceSlot_;
      const SimplexId destinationSlot = vpath.destinationSlot_;

      CriticalPoint &source = criticalPoints[sourceIndex];
      CriticalPoint &destination = criticalPoints[destinationIndex];

      source.vpaths_[sourceSlot] = i;
      destination.vpaths_[destinationSlot] = i;
    }
  }

  this->printMsg(
    " Initialization step #2", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}


template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections1(
  const std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>> &S) {
  Timer t;

  const SimplexId numberOfVPaths = vpaths.size();
  for(SimplexId i = 0; i < numberOfVPaths; ++i) {
    const VPath &vpath = vpaths[i];

    if(vpath.isValid_) {
      const SimplexId saddleId = criticalPoints[vpath.destination_].cell_.id_;
      S.insert(std::make_tuple(vpath.persistence_, i, saddleId));
    }
  }

  this->printMsg(
    " Ordering of the vpaths #1", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}


template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections2(
  const std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>> &S) {
  Timer t;

  const SimplexId numberOfVPaths = vpaths.size();
  for(SimplexId i = 0; i < numberOfVPaths; ++i) {
    const VPath &vpath = vpaths[i];

    if(vpath.isValid_) {
      const SimplexId saddleId = criticalPoints[vpath.source_].cell_.id_;
      S.insert(std::make_tuple(vpath.persistence_, i, saddleId));
    }
  }

  this->printMsg(
    " Ordering of the vpaths #2", 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::processSaddleSaddleConnections1(
  const int iterationThreshold,
  const std::vector<char> &isPL,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors,
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>> &S,
  std::vector<SimplexId> &pl2dmt_saddle1,
  std::vector<SimplexId> &pl2dmt_saddle2,
  std::vector<char> &isRemovableSaddle1,
  std::vector<char> &isRemovableSaddle2,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index,
  const triangulationType &triangulation) {
  Timer t;

  const auto *const scalars
    = static_cast<const dataType *>(inputScalarField_.first);

  const SimplexId numberOfEdges = triangulation.getNumberOfEdges();
  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  const SimplexId optimizedSize = std::max(numberOfEdges, numberOfTriangles);
  std::vector<bool> isVisited(optimizedSize, false);
  std::vector<SimplexId> visitedCells{};

  int numberOfIterations{};
  while(!S.empty()) {
    if(iterationThreshold >= 0 and numberOfIterations >= iterationThreshold) {
      break;
    }

    auto ptr = S.begin();
    const SimplexId vpathId = std::get<1>(*ptr);
    S.erase(ptr);
    VPath &vpath = vpaths[vpathId];

    if(vpath.isValid_) {
      if(returnSaddleConnectors) {
        const dataType persistence = vpath.persistence_;
        if(persistence > SaddleConnectorsPersistenceThreshold) {
          break;
        }
      }

      const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
      const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;

      std::vector<SimplexId> saddles1;
      VisitedMask mask{isVisited, visitedCells};
      getDescendingWall(minSaddle2, mask, triangulation, nullptr, &saddles1);

      // check if at least one connection exists
      auto isFound = std::find_if(
      saddles1.begin(), saddles1.end(),
      [&](const auto &s) { return s == minSaddle1.id_; });      
      if(isFound == saddles1.end()) {
        ++numberOfIterations;
        continue;
      }

      // check if there is multiple connections
      std::vector<Cell> path;
      const bool isMultiConnected = getAscendingPathThroughWall(
        minSaddle1, minSaddle2, isVisited, &path, triangulation, true);

      if(isMultiConnected) {
        ++numberOfIterations;
        continue;
      }

      // filter by 1-saddle condition
      if(vpath.isValid_) {
        const Cell &dmt_saddle1 = criticalPoints[vpath.source_].cell_;
        const SimplexId dmt_saddle1Id = dmt_saddle1.id_;

        if(isSaddle1(dmt_saddle1)) {
          for(int i = 0; i < 2; ++i) {
            SimplexId vertexId;
            triangulation.getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle1)) {
              continue;
            }

            if(!allowBoundary and triangulation.isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle1[vertexId] == -1) {
              const SimplexId pl_saddle1Id = vertexId;

              SimplexId numberOfRemainingSaddles1 = 0;

              SimplexId savedId = -1;
              const SimplexId edgeNumber
                = triangulation.getVertexEdgeNumber(pl_saddle1Id);
              for(SimplexId j = 0; j < edgeNumber; ++j) {
                SimplexId edgeId;
                triangulation.getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId != dmt_saddle1Id and isSaddle1(Cell(1, edgeId))
                   and isRemovableSaddle1[edgeId]) {
                  ++numberOfRemainingSaddles1;
                  savedId = edgeId;
                }
              }

              if(numberOfRemainingSaddles1 == 0) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                pl2dmt_saddle1[vertexId] = dmt_saddle1Id;
                //dmt1Saddle2PL_[dmt_saddle1Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1 == 1) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                isRemovableSaddle1[savedId] = false;
                pl2dmt_saddle1[vertexId] = savedId;
                //dmt1Saddle2PL_[savedId] = vertexId;
                break;
              }
            } else if(pl2dmt_saddle1[vertexId] == dmt_saddle1Id) {
              vpath.invalidate();
              break;
            }
          }
        } else {
          vpath.invalidate();
        }
      }

      // filter by 2-saddle condition
      if(!allowBruteForce and vpath.isValid_) {
        const Cell &dmt_saddle2 = criticalPoints[vpath.destination_].cell_;
        const SimplexId dmt_saddle2Id = dmt_saddle2.id_;

        if(isSaddle2(dmt_saddle2)) {
          for(int i = 0; i < 3; ++i) {
            SimplexId vertexId;
            triangulation.getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle2)) {
              continue;
            }

            if(!allowBoundary and triangulation.isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle2[vertexId] == -1) {
              const SimplexId pl_saddle2Id = vertexId;

              SimplexId numberOfRemainingSaddles2 = 0;

              SimplexId savedId = -1;
              const SimplexId triangleNumber
                = triangulation.getVertexTriangleNumber(pl_saddle2Id);
              for(SimplexId j = 0; j < triangleNumber; ++j) {
                SimplexId triangleId;
                triangulation.getVertexTriangle(pl_saddle2Id, j, triangleId);

                if(triangleId != dmt_saddle2Id
                   and isSaddle2(Cell(2, triangleId))
                   and isRemovableSaddle2[triangleId]) {
                  ++numberOfRemainingSaddles2;
                  savedId = triangleId;
                }
              }

              if(numberOfRemainingSaddles2 == 0) {
                isRemovableSaddle2[dmt_saddle2Id] = false;
                pl2dmt_saddle2[vertexId] = dmt_saddle2Id;
                //dmt2Saddle2PL_[dmt_saddle2Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles2 == 1) {
                isRemovableSaddle2[dmt_saddle2Id] = false;
                isRemovableSaddle2[savedId] = false;
                pl2dmt_saddle2[vertexId] = savedId;
                //dmt2Saddle2PL_[savedId] = vertexId;
                break;
              }
            } else if(pl2dmt_saddle2[vertexId] == dmt_saddle2Id) {
              vpath.invalidate();
              break;
            }
          }
        } else {
          vpath.invalidate();
        }
      }

      if(vpath.isValid_) {
        reverseAscendingPathOnWall(path, triangulation);
      }
    }

    if(vpath.isValid_) {
      // add persistence pair to collection if necessary
      //if(CollectPersistencePairs and outputPersistencePairs_) {
      //  const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
      //  const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;
      //  outputPersistencePairs_->push_back({minSaddle1, minSaddle2});
      //}

      const SimplexId sourceId = vpath.source_;
      const SimplexId destinationId = vpath.destination_;

      // invalidate vpaths connected to destination
      std::vector<SimplexId> newSourceIds;
      CriticalPoint &destination = criticalPoints[destinationId];
      for(auto &destinationVPathId : destination.vpaths_) {
        VPath &destinationVPath = vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_ != sourceId) {
          // save critical point
          const SimplexId newSourceId = destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          // clear vpath
          destinationVPath.invalidate();
        }
      }

      // invalidate vpaths connected to source and save the critical points to
      // update
      std::vector<SimplexId> newDestinationIds;
      CriticalPoint &source = criticalPoints[sourceId];
      for(auto &sourceVPathId : source.vpaths_) {
        VPath &sourceVPath = vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_ != destinationId) {
          // save critical point
          const SimplexId newDestinationId = sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          CriticalPoint &newDestination = criticalPoints[newDestinationId];
          for(auto &newDestinationVPathId : newDestination.vpaths_) {
            VPath &newDestinationVPath = vpaths[newDestinationVPathId];
            if(newDestinationVPath.isValid_
               and newDestinationVPath.source_ != sourceId) {

              // clear vpath
              newDestinationVPath.invalidate();
            }
          }

          // clear vpath
          sourceVPath.invalidate();
        }
      }

      // finally invalidate current vpath and critical points
      vpath.invalidate();
      source.clear();
      destination.clear();

      // look at the gradient : reconnect locally the critical points
      for(auto &newDestinationId : newDestinationIds) {
        CriticalPoint &newDestination = criticalPoints[newDestinationId];
        const Cell &saddle2 = newDestination.cell_;

        std::vector<SimplexId> saddles1;
        VisitedMask mask{isVisited, visitedCells};
        getDescendingWall(saddle2, mask, triangulation, nullptr, &saddles1);

        for(auto &saddle1Id : saddles1) {
          const Cell saddle1(1, saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getAscendingPathThroughWall(
            saddle1, saddle2, isVisited, &path, triangulation, true);

          if(isMultiConnected) {
            continue;
          }

          SimplexId newSourceId = saddle1Index[saddle1Id];

          // connection to a new saddle1 (not present in the graph before)
          if(newSourceId == -1) {
            if(!isRemovableSaddle1[saddle1Id]) {
              continue;
            }

            const SimplexId newCriticalPointId = criticalPoints.size();
            saddle1Index[saddle1Id] = newCriticalPointId;
            criticalPoints.emplace_back(saddle1);

            newSourceId = newCriticalPointId;
          }
          CriticalPoint &newSource = criticalPoints[newSourceId];

          // update vpaths
          const SimplexId newVPathId = vpaths.size();
          const auto persistence
            = getPersistence(saddle2, saddle1, scalars, triangulation);

          vpaths.push_back(VPath(
            true, -1, newSourceId, newDestinationId, -1, -1, persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(
            std::make_tuple(persistence, newVPathId, newDestination.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(auto &newSourceId : newSourceIds) {
        CriticalPoint &newSource = criticalPoints[newSourceId];
        const Cell &saddle1 = newSource.cell_;

        std::vector<SimplexId> saddles2;
        VisitedMask mask{isVisited, visitedCells};
        getAscendingWall(saddle1, mask, triangulation, nullptr, &saddles2);

        for(auto &saddle2Id : saddles2) {
          const Cell saddle2(2, saddle2Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getDescendingPathThroughWall(
            saddle2, saddle1, isVisited, &path, triangulation, true);

          if(isMultiConnected) {
            continue;
          }

          const SimplexId newDestinationId = saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId == -1) {
            continue;
          }

          CriticalPoint &newDestination = criticalPoints[newDestinationId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists = false;
          for(auto &newDestinationVPathId : newDestination.vpaths_) {
            const VPath &newDestinationVPath = vpaths[newDestinationVPathId];

            if(newDestinationVPath.isValid_
               and newDestinationVPath.source_ == newSourceId) {
              alreadyExists = true;
              break;
            }
          }

          if(alreadyExists) {
            continue;
          }

          // update vpaths
          const SimplexId newVPathId = vpaths.size();
          const auto persistence
            = getPersistence(saddle2, saddle1, scalars, triangulation);

          vpaths.push_back(VPath(
            true, -1, newSourceId, newDestinationId, -1, -1, persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(
            std::make_tuple(persistence, newVPathId, newDestination.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  this->printMsg(" Processing of the vpaths #1", 1.0, t.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}


template <typename dataType, typename triangulationType>
int DiscreteGradient::  processSaddleSaddleConnections2(
  const int iterationThreshold,
  const std::vector<char> &isPL,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors,
  std::set<std::tuple<dataType, SimplexId, SimplexId>,
           SaddleSaddleVPathComparator<dataType>> &S,
  std::vector<SimplexId> &pl2dmt_saddle1,
  std::vector<SimplexId> &pl2dmt_saddle2,
  std::vector<char> &isRemovableSaddle1,
  std::vector<char> &isRemovableSaddle2,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index,
  const triangulationType &triangulation) {
  Timer t;

  this->printMsg("Saddle connector persistence threshold: "
                 + std::to_string(this->SaddleConnectorsPersistenceThreshold));

  const auto *const scalars
    = static_cast<const dataType *>(inputScalarField_.first);

  const SimplexId numberOfEdges = triangulation.getNumberOfEdges();
  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  const SimplexId optimizedSize = std::max(numberOfEdges, numberOfTriangles);
  std::vector<bool> isVisited(optimizedSize, false);
  std::vector<SimplexId> visitedIds{};

  int numberOfIterations{};
  while(!S.empty()) {
    if(iterationThreshold >= 0 and numberOfIterations >= iterationThreshold) {
      break;
    }

    auto ptr = S.begin();
    const SimplexId vpathId = std::get<1>(*ptr);
    S.erase(ptr);
    VPath &vpath = vpaths[vpathId];

    if(vpath.isValid_) {
      if(returnSaddleConnectors) {
        const dataType persistence = vpath.persistence_;
        if(persistence > SaddleConnectorsPersistenceThreshold) {
          break;
        }
      }

      const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
      const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;

      std::vector<SimplexId> saddles2;
      VisitedMask mask{isVisited, visitedIds};
      getAscendingWall(minSaddle1, mask, triangulation, nullptr, &saddles2);

      // check if at least one connection exists
      auto isFound = std::find_if(
      saddles2.begin(), saddles2.end(),
      [&](const auto &s) { return s == minSaddle2.id_; });      

      if(isFound == saddles2.end()) {
        ++numberOfIterations;
        continue;
      }

      // check if there is multiple connections
      std::vector<Cell> path;
      const bool isMultiConnected = getDescendingPathThroughWall(
        minSaddle2, minSaddle1, isVisited, &path, triangulation, true);

      if(isMultiConnected) {
        ++numberOfIterations;
        continue;
      }

      // filter by 1-saddle condition
      if(vpath.isValid_) {
        const Cell &dmt_saddle1 = criticalPoints[vpath.source_].cell_;
        const SimplexId dmt_saddle1Id = dmt_saddle1.id_;

        if(isSaddle1(dmt_saddle1)) {
          for(int i = 0; i < 2; ++i) {
            SimplexId vertexId;
            triangulation.getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle1)) {
              continue;
            }

            if(!allowBoundary and triangulation.isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle1[vertexId] == -1) {
              const SimplexId pl_saddle1Id = vertexId;

              SimplexId numberOfRemainingSaddles1 = 0;

              SimplexId savedId = -1;
              const SimplexId edgeNumber
                = triangulation.getVertexEdgeNumber(pl_saddle1Id);
              for(SimplexId j = 0; j < edgeNumber; ++j) {
                SimplexId edgeId;
                triangulation.getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId != dmt_saddle1Id and isSaddle1(Cell(1, edgeId))
                   and isRemovableSaddle1[edgeId]) {
                  ++numberOfRemainingSaddles1;
                  savedId = edgeId;
                }
              }

              if(numberOfRemainingSaddles1 == 0) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                pl2dmt_saddle1[vertexId] = dmt_saddle1Id;
                //dmt1Saddle2PL_[dmt_saddle1Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1 == 1) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                isRemovableSaddle1[savedId] = false;
                pl2dmt_saddle1[vertexId] = savedId;
                //dmt1Saddle2PL_[savedId] = vertexId;
                break;
              }
            } else if(pl2dmt_saddle1[vertexId] == dmt_saddle1Id) {
              vpath.invalidate();
              break;
            }
          }
        } else {
          vpath.invalidate();
        }
      }

      // filter by 2-saddle condition
      if(!allowBruteForce and vpath.isValid_) {
        const Cell &dmt_saddle2 = criticalPoints[vpath.destination_].cell_;
        const SimplexId dmt_saddle2Id = dmt_saddle2.id_;

        if(isSaddle2(dmt_saddle2)) {
          for(int i = 0; i < 3; ++i) {
            SimplexId vertexId;
            triangulation.getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle2)) {
              continue;
            }

            if(!allowBoundary and triangulation.isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle2[vertexId] == -1) {
              const SimplexId pl_saddle2Id = vertexId;

              SimplexId numberOfRemainingSaddles2 = 0;

              SimplexId savedId = -1;
              const SimplexId triangleNumber
                = triangulation.getVertexTriangleNumber(pl_saddle2Id);
              for(SimplexId j = 0; j < triangleNumber; ++j) {
                SimplexId triangleId;
                triangulation.getVertexTriangle(pl_saddle2Id, j, triangleId);

                if(triangleId != dmt_saddle2Id
                   and isSaddle2(Cell(2, triangleId))
                   and isRemovableSaddle2[triangleId]) {
                  ++numberOfRemainingSaddles2;
                  savedId = triangleId;
                }
              }

              if(!numberOfRemainingSaddles2) {
                isRemovableSaddle2[dmt_saddle2Id] = false;
                pl2dmt_saddle2[vertexId] = dmt_saddle2Id;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles2 == 1) {
                isRemovableSaddle2[dmt_saddle2Id] = false;
                isRemovableSaddle2[savedId] = false;
                pl2dmt_saddle2[vertexId] = savedId;
                break;
              }
            } else if(pl2dmt_saddle2[vertexId] == dmt_saddle2Id) {
              vpath.invalidate();
              break;
            }
          }
        } else {
          vpath.invalidate();
        }
      }

      if(vpath.isValid_) {
        reverseDescendingPathOnWall(path, triangulation);
      }
    }

    if(vpath.isValid_) {
      // add persistence pair to collection if necessary
      //if(CollectPersistencePairs and outputPersistencePairs_) {
      //  const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
      //  const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;
      //  outputPersistencePairs_->push_back({minSaddle1, minSaddle2});
      //}

      const SimplexId sourceId = vpath.source_;
      const SimplexId destinationId = vpath.destination_;

      // invalidate vpaths connected to source
      std::vector<SimplexId> newDestinationIds;
      CriticalPoint &source = criticalPoints[sourceId];
      for(auto &sourceVPathId : source.vpaths_) {
        VPath &sourceVPath = vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_ != destinationId) {
          // save critical point
          const SimplexId newDestinationId = sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          // clear vpath
          sourceVPath.invalidate();
        }
      }

      // invalidate vpaths connected to destination and save the critical
      // points to update
      std::vector<SimplexId> newSourceIds;
      CriticalPoint &destination = criticalPoints[destinationId];
      for(auto &destinationVPathId : destination.vpaths_) {
        VPath &destinationVPath = vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_ != sourceId) {
          // save critical point
          const SimplexId newSourceId = destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          CriticalPoint &newSource = criticalPoints[newSourceId];
          for(auto &newSourceVPathId : newSource.vpaths_) {
            VPath &newSourceVPath = vpaths[newSourceVPathId];
            if(newSourceVPath.isValid_
               and newSourceVPath.destination_ != destinationId) {

              // clear vpath
              newSourceVPath.invalidate();
            }
          }

          // clear vpath
          destinationVPath.invalidate();
        }
      }

      // finally invalidate current vpath and critical points
      vpath.invalidate();
      source.clear();
      destination.clear();

      // look at the gradient : reconnect locally the critical points
      for(auto &newSourceId : newSourceIds) {
        CriticalPoint &newSource = criticalPoints[newSourceId];
        const Cell &saddle1 = newSource.cell_;

        std::vector<SimplexId> saddles2;
        VisitedMask mask{isVisited, visitedIds};
        getAscendingWall(saddle1, mask, triangulation, nullptr, &saddles2);

        for(auto &saddle2Id : saddles2) {
          const Cell saddle2(2, saddle2Id);

          const bool isMultiConnected = getDescendingPathThroughWall(
            saddle2, saddle1, isVisited, nullptr, triangulation, true);

          if(isMultiConnected) {
            continue;
          }

          SimplexId newDestinationId = saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId == -1) {
            if(!isRemovableSaddle2[saddle2Id]) {
              continue;
            }

            const SimplexId newCriticalPointId = criticalPoints.size();
            saddle2Index[saddle2Id] = newCriticalPointId;
            criticalPoints.emplace_back(saddle2);

            newDestinationId = newCriticalPointId;
          }

          CriticalPoint &newDestination = criticalPoints[newDestinationId];

          // update vpaths
          const SimplexId newVPathId = vpaths.size();
          const auto persistence
            = getPersistence(saddle2, saddle1, scalars, triangulation);

          vpaths.push_back(VPath(
            true, -1, newSourceId, newDestinationId, -1, -1, persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(
            std::make_tuple(persistence, newVPathId, newSource.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(auto &newDestinationId : newDestinationIds) {
        CriticalPoint &newDestination = criticalPoints[newDestinationId];
        const Cell &saddle2 = newDestination.cell_;

        std::vector<SimplexId> saddles1;
        VisitedMask mask{isVisited, visitedIds};
        getDescendingWall(saddle2, mask, triangulation, nullptr, &saddles1);

        for(auto &saddle1Id : saddles1) {
          const Cell saddle1(1, saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getAscendingPathThroughWall(
            saddle1, saddle2, isVisited, &path, triangulation, true);

          if(isMultiConnected) {
            continue;
          }

          const SimplexId newSourceId = saddle1Index[saddle1Id];

          if(newSourceId == -1) {
            continue;
          }

          CriticalPoint &newSource = criticalPoints[newSourceId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists = false;
          for(auto &newSourceVPathId : newSource.vpaths_) {
            const VPath &newSourceVPath = vpaths[newSourceVPathId];

            if(newSourceVPath.isValid_
               and newSourceVPath.destination_ == newDestinationId) {
              alreadyExists = true;
              break;
            }
          }

          if(alreadyExists) {
            continue;
          }

          // update vpaths
          const SimplexId newVPathId = vpaths.size();
          const auto persistence
            = getPersistence(saddle2, saddle1, scalars, triangulation);

          vpaths.push_back(VPath(
            true, -1, newSourceId, newDestinationId, -1, -1, persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(
            std::make_tuple(persistence, newVPathId, newSource.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  this->printMsg(" Processing of the vpaths #2", 1.0, t.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

