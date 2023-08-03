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
         - scalars[getCellLowerVertex(down, triangulation)];
}

template <typename triangulationType>
int DiscreteGradient::buildGradient(const triangulationType &triangulation,
                                    bool bypassCache) {

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
    // compute gradient pairs
    this->processLowerStars(this->inputOffsets_, triangulation);

    this->printMsg(
      "Built discrete gradient", 1.0, tm.getElapsedTime(), this->threadNumber_);
  } else {
    this->printMsg("Fetched cached discrete gradient");
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
inline void
  DiscreteGradient::lowerStar(lowerStarType &ls,
                              const SimplexId a,
                              const SimplexId *const offsets,
                              const triangulationType &triangulation) const {

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
  std::vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfEdges, 0);
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
        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
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
  std::vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfCells, 0);
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
          if(isCycle[currentId] == 0) {
            isCycle[currentId] = 1;
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
  const bool enableCycleDetector) const {

  // debug
  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  std::vector<char> isCycle;
  if(enableCycleDetector) {
    isCycle.resize(numberOfTriangles, 0);
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
        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
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
  const triangulationType &triangulation) const {

  if(dimensionality_ == 3) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;
      const SimplexId triangleId = vpath[i + 1].id_;

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
