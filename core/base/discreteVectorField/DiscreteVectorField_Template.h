#pragma once

#include <DiscreteGradient.h>
#include <DiscreteVectorField.h>

using ttk::SimplexId;
using ttk::dcg::Cell;
using ttk::dcvf::CellOutExt;
using ttk::dcvf::DiscreteVectorField;

template <typename dataType, typename triangulationType>
float DiscreteVectorField::getPersistence(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {
  // Function will trace the provided vpath to calculate something similar to
  // 'persistence' value
  // Loop over vpath
  float result = 0;
  const SimplexId numberOfCellsInPath = vpath.size();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfCellsInPath > 0) {
    // Validating starting saddle assumption
    if(vpath[0].dim_ != 1 && this->dimensionality_ == 2) {
      this->printErr(
        "Persistence calculation must start with a critical saddle(dim=1)");
    }
  }
#endif
  for(SimplexId i = 0; i < numberOfCellsInPath - 1; i++) {
    Cell firstCell = vpath[i];
    vectorValue firstCenter;
    vectorValue firstVector;
    Cell secondCell = vpath[i + 1];
    vectorValue secondCenter;
    vectorValue secondVector;
    SimplexId v;
    vectorValue pointAdd;
    // Compute average point and vector value depending on dimension of cell
    switch(firstCell.dim_) {
      case 0:
        triangulation.getVertexPoint(
          firstCell.id_, firstCenter[0], firstCenter[1], firstCenter[2]);
        firstVector = getVectorValueAt<dataType>(firstCell.id_);
        break;
      case 1:
        triangulation.getEdgeVertex(firstCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, firstCenter[0], firstCenter[1], firstCenter[2]);
        firstVector = getVectorValueAt<dataType>(v);
        triangulation.getEdgeVertex(firstCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstCenter = firstCenter / 2;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        break;
      case 2:
        triangulation.getTriangleVertex(firstCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, firstCenter[0], firstCenter[1], firstCenter[2]);
        firstVector = getVectorValueAt<dataType>(v);
        triangulation.getTriangleVertex(firstCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        triangulation.getTriangleVertex(firstCell.id_, 2, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstCenter = firstCenter / 3;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        break;
      case 3:
        triangulation.getCellVertex(firstCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, firstCenter[0], firstCenter[1], firstCenter[2]);
        firstVector = getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(firstCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(firstCell.id_, 2, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(firstCell.id_, 3, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        firstCenter = firstCenter + pointAdd;
        firstCenter = firstCenter / 4;
        firstVector = firstVector + getVectorValueAt<dataType>(v);
        break;
    }
    switch(secondCell.dim_) {
      case 0:
        triangulation.getVertexPoint(
          secondCell.id_, secondCenter[0], secondCenter[1], secondCenter[2]);
        secondVector = getVectorValueAt<dataType>(secondCell.id_);
        break;
      case 1:
        triangulation.getEdgeVertex(secondCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, secondCenter[0], secondCenter[1], secondCenter[2]);
        secondVector = getVectorValueAt<dataType>(v);
        triangulation.getEdgeVertex(secondCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondCenter = secondCenter / 2;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        break;
      case 2:
        triangulation.getTriangleVertex(secondCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, secondCenter[0], secondCenter[1], secondCenter[2]);
        secondVector = getVectorValueAt<dataType>(v);
        triangulation.getTriangleVertex(secondCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        triangulation.getTriangleVertex(secondCell.id_, 2, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondCenter = secondCenter / 3;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        break;
      case 3:
        triangulation.getCellVertex(secondCell.id_, 0, v);
        triangulation.getVertexPoint(
          v, secondCenter[0], secondCenter[1], secondCenter[2]);
        secondVector = getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(secondCell.id_, 1, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(secondCell.id_, 2, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        triangulation.getCellVertex(secondCell.id_, 3, v);
        triangulation.getVertexPoint(v, pointAdd[0], pointAdd[1], pointAdd[2]);
        secondCenter = secondCenter + pointAdd;
        secondCenter = secondCenter / 4;
        secondVector = secondVector + getVectorValueAt<dataType>(v);
        break;
    }
    vectorValue changeDirection;
    // Handle discrete vector direction(toward higher dim)
    if(firstCell.dim_ < secondCell.dim_) {
      changeDirection = secondCenter - firstCenter;
    } else {
      changeDirection = firstCenter - secondCenter;
    }
    vectorValue avgVector = (firstVector + secondVector) / 2;
    if(i % 2 == 0) { // Handle alternation in Vpath calculation
      result -= Geometry::dotProduct(changeDirection.data(), avgVector.data());
    } else {
      result += Geometry::dotProduct(changeDirection.data(), avgVector.data());
    }

  } // End VPath loop

  return result;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::buildField(const triangulationType &triangulation) {

  // set member variables at each buildField() call
  this->dimensionality_ = triangulation.getCellVertexNumber(0) - 1;
  this->numberOfVertices_ = triangulation.getNumberOfVertices();

  this->vectors_ = &this->localVectors_;
  // allocate field memory
  this->initMemory(triangulation);

  Timer tm{};
  // compute discrete pairs
  this->processOutwardStars<dataType, triangulationType>(triangulation);

  this->printMsg(
    "Built discrete vectors", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::setCriticalPoints(
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
      isOnBoundary[o]
        = this->isBoundary<dataType, triangulationType>(cell, triangulation);
      PLVertexIdentifiers[o]
        = this->getCellGreaterVertex<dataType, triangulationType>(
          cell, triangulation);
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

template <typename dataType, typename triangulationType>
int DiscreteVectorField::setCriticalPoints(
  std::vector<std::array<float, 3>> &points,
  std::vector<char> &cellDimensions,
  std::vector<SimplexId> &cellIds,
  std::vector<char> &isOnBoundary,
  std::vector<SimplexId> &PLVertexIdentifiers,
  const triangulationType &triangulation) const {

  std::array<std::vector<SimplexId>, 4> criticalCellsByDim;
  getCriticalPoints(criticalCellsByDim, triangulation);
  setCriticalPoints<dataType, triangulationType>(
    criticalCellsByDim, points, cellDimensions, cellIds, isOnBoundary,
    PLVertexIdentifiers, triangulation);

  return 0;
}

template <typename triangulationType>
int DiscreteVectorField::getCriticalPoints(
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
      // std::cout << "Testing " << i << "," << j << std::endl;
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
SimplexId DiscreteVectorField::getNumberOfCells(
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

template <typename dataType, typename triangulationType>
inline bool DiscreteVectorField::compare(const triangulationType &triangulation,
                                         SimplexId vertexA,
                                         SimplexId vertexB,
                                         float &weightValue) const {
  weightValue = 0; // Silences unassignment warnings

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexA == -1 || vertexB == -1) {
    this->printErr("Passed Null value to compare()");
    return vertexA < vertexB;
  }
  SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  if(vertexA < 0 || vertexB < 0 || vertexA >= numberOfVertices
     || vertexB >= numberOfVertices) {
    std::cout << vertexA << "," << vertexB << std::endl;
    this->printErr("Passed invalid vertex values");
    return vertexA < vertexB;
  }
#endif
  // Compute the edge direction a->b (c(b)-c(a))
  vectorValue locationA;
  triangulation.getVertexPoint(
    vertexA, locationA[0], locationA[1], locationA[2]);

  vectorValue locationB;
  triangulation.getVertexPoint(
    vertexB, locationB[0], locationB[1], locationB[2]);
  vectorValue edgeDirection = locationB - locationA;

  // Averaged vector computation
  auto averageVector = (getVectorValueAt<dataType>(vertexA)
                        + getVectorValueAt<dataType>(vertexB))
                       / 2;
  weightValue
    = Geometry::dotProduct(edgeDirection.data(), averageVector.data());

  if(weightValue > 0)
    return true;
  if(weightValue < 0)
    return false;
  return vertexA < vertexB; // Tiebreaking scheme on id, lower ID 'wins'
}

template <typename dataType, typename triangulationType>
inline void DiscreteVectorField::OutwardStar(
  outwardStarType &os,
  const SimplexId a,
  const triangulationType &triangulation) const {

  // make sure that os is cleared
  for(auto &vec : os) {
    vec.clear();
  }

  // a belongs to its outward star
  os[0].emplace_back(CellOutExt{0, a});

  // store outward edges
  const auto nedges = triangulation.getVertexEdgeNumber(a);
  os[1].reserve(nedges);
  for(SimplexId i = 0; i < nedges; i++) {
    SimplexId edgeId;
    triangulation.getVertexEdge(a, i, edgeId);
    SimplexId vertexId;
    triangulation.getEdgeVertex(edgeId, 0, vertexId);
    if(vertexId == a) {
      triangulation.getEdgeVertex(edgeId, 1, vertexId);
    }
    float edgeWeight{0};
    if(compare<dataType, triangulationType>(
         triangulation, a, vertexId, edgeWeight)) {
      os[1].emplace_back(CellOutExt{1,
                                    edgeId,
                                    {vertexId, -1, -1},
                                    {-(edgeWeight), -INFINITY, -INFINITY},
                                    {}});
    }
  }

  if(os[1].size() < 2) {
    // at least two edges in the outward star for one triangle
    return;
  }

  const auto processTriangle
    = [&](const SimplexId triangleId, const SimplexId v0, const SimplexId v1,
          const SimplexId v2) {
        std::array<SimplexId, 3> lowVerts{-1, -1, -1};
        if(v0 == a) {
          lowVerts[0] = v1;
          lowVerts[1] = v2;
        } else if(v1 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v2;
        } else if(v2 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
        }

        float edgeWeight1, edgeWeight2;
        if(compare<dataType, triangulationType>(
             triangulation, a, lowVerts[0], edgeWeight1)
           && compare<dataType, triangulationType>(
             triangulation, a, lowVerts[1],
             edgeWeight2)) { // triangle in outwardStar
          uint8_t j{}, k{};
          // store edges indices(in outward star) of current triangle

          std::array<uint8_t, 3>
            faces{}; // Might need to find faces differently(simply replace if
                     // condition)
          for(const auto &e : os[1]) {
            if(e.lowVerts_[0] == lowVerts[0] || e.lowVerts_[0] == lowVerts[1]) {
              faces[k++] = j;
            }
            j++;
          }

          std::array<float, 3> lowVertWeights
            = {-(edgeWeight1), -(edgeWeight2), -INFINITY};
          if(edgeWeight2
             < edgeWeight1) { // Flip order if necessary (decreasing negated)
            lowVertWeights[0] = -(edgeWeight2);
            lowVertWeights[1] = -(edgeWeight1);
          }
          os[2].emplace_back(
            CellOutExt{2, triangleId, lowVerts, lowVertWeights, faces});
        }
      };

  if(dimensionality_ == 2) {
    // store outward triangles

    // use optimised triangulation methods:
    // getVertexStar instead of getVertexTriangle
    // getCellVertex instead of getTriangleVertex
    const auto ncells = triangulation.getVertexStarNumber(a);
    os[2].reserve(ncells);
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
    // store outward triangles
    const auto ntri = triangulation.getVertexTriangleNumber(a);
    os[2].reserve(ntri);
    for(SimplexId i = 0; i < ntri; i++) {
      SimplexId triangleId;
      triangulation.getVertexTriangle(a, i, triangleId);
      SimplexId v0{}, v1{}, v2{};
      triangulation.getTriangleVertex(triangleId, 0, v0);
      triangulation.getTriangleVertex(triangleId, 1, v1);
      triangulation.getTriangleVertex(triangleId, 2, v2);
      processTriangle(triangleId, v0, v1, v2);
    }

    // at least three triangles in the outward star for one tetra
    if(os[2].size() >= 3) {
      // store outward tetra
      const auto ncells = triangulation.getVertexStarNumber(a);
      os[3].reserve(ncells);
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
          lowVerts[0] = v1;
          lowVerts[1] = v2;
          lowVerts[2] = v3;
        } else if(v1 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v2;
          lowVerts[2] = v3;
        } else if(v2 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
          lowVerts[2] = v3;
        } else if(v3 == a) {
          lowVerts[0] = v0;
          lowVerts[1] = v1;
          lowVerts[2] = v2;
        }
        // Need to check each of the edges individually with compare function
        float eWeight0, eWeight1, eWeight2;
        if(compare<dataType, triangulationType>(
             triangulation, a, lowVerts[0], eWeight0)
           && compare<dataType, triangulationType>(
             triangulation, a, lowVerts[1], eWeight1)
           && compare<dataType, triangulationType>(
             triangulation, a, lowVerts[2], eWeight2)) { // tetra in outwardStar

          uint8_t j{}, k{};
          // store triangles indices of current tetra
          std::array<uint8_t, 3> faces{};
          for(const auto &t : os[2]) {
            // lowVerts & t.lowVerts are not ordered, so need to check if
            // t.lowVerts has both elements in lowVerts
            bool vert1Found
              = std::find(lowVerts.begin(), lowVerts.end(), t.lowVerts_[0])
                != lowVerts.end();
            bool vert2Found
              = std::find(lowVerts.begin(), lowVerts.end(), t.lowVerts_[1])
                != lowVerts.end();

            if(vert1Found && vert2Found) {
              faces[k++] = j;
            }
            j++;
          }
          std::array<float, 3> lowVertWeights
            = {-(eWeight0), -(eWeight1), -(eWeight2)};
          std::sort(lowVertWeights.rbegin(),
                    lowVertWeights
                      .rend()); // Sort in decreasing order e.g.,(-1, -2, -3)
          os[3].emplace_back(
            CellOutExt{3, cellId, lowVerts, lowVertWeights, faces});
        }
      }
    }
  }
}

template <typename triangulationType>
inline void DiscreteVectorField::pairCells(
  CellOutExt &alpha, CellOutExt &beta, const triangulationType &triangulation) {
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
  (*vectors_)[2 * alpha.dim_][alpha.id_] = localBId;
  (*vectors_)[2 * alpha.dim_ + 1][beta.id_] = localAId;
#else
  TTK_FORCE_USE(triangulation);
  (*vectors_)[2 * alpha.dim_][alpha.id_] = beta.id_;
  (*vectors_)[2 * alpha.dim_ + 1][beta.id_] = alpha.id_;
#endif // TTK_ENABLE_DCG_OPTIMIZE_MEMORY
  alpha.paired_ = true;
  beta.paired_ = true;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::processOutwardStars(
  const triangulationType &triangulation) {

  /* Compute discrete vector field */

  auto nverts = triangulation.getNumberOfVertices();

  // Comparison function for Cells inside priority queues
  const auto orderCells
    = [&](const CellOutExt &a, const CellOutExt &b) -> bool {
    return a.lowVertWeights_ > b.lowVertWeights_;
  };

  // Type alias for priority queues
  using pqType
    = std::priority_queue<std::reference_wrapper<CellOutExt>,
                          std::vector<std::reference_wrapper<CellOutExt>>,
                          decltype(orderCells)>;

  // To reduce allocations, priority queues and outwardStar objects are
  // cleaned & reused between iterations.

  // Priority queues are pushed at the beginning and popped at the
  // end. To pop the minimum, elements should be sorted in a
  // decreasing order.
  pqType pqZero{orderCells}, pqOne{orderCells};

  // store outward star structure
  outwardStarType Lx;

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
    const auto insertCofacets = [&](const CellOutExt &ca, outwardStarType &ls) {
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

    OutwardStar<dataType, triangulationType>(Lx, x, triangulation);
    // In case the vertex is a ghost, the gradient of the
    // simplices of its star is set to GHOST_CONNECTION
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
        // get delta: 1-cell (edge) with minimal negated weight value (steeper
        // 'descent')
        size_t minId = 0;
        for(size_t i = 1; i < Lx[1].size(); ++i) {
          const auto &a = Lx[1][minId].lowVertWeights_[0];
          const auto &b = Lx[1][i].lowVertWeights_[0];
          if(a > b) { // Check for the smallest negated weight
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

template <typename dataType, typename triangulationType>
bool DiscreteVectorField::isBoundary(
  const Cell &cell, const triangulationType &triangulation) const {

  if(cell.dim_ > this->dimensionality_ || cell.dim_ < 0 || cell.id_ < 0) {
    return false;
  }

  const auto vert{this->getCellGreaterVertex<dataType, triangulationType>(
    cell, triangulation)};
  if(vert == -1) {
    return false;
  }
  return triangulation.isVertexOnBoundary(vert);
}

template <typename triangulationType>
SimplexId
  DiscreteVectorField::getPairedCell(const Cell &cell,
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
      const auto locId{(*vectors_)[0][cell.id_]};
      if(locId != -1) {
        triangulation.getVertexEdge(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[0][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 1) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*vectors_)[1][cell.id_]};
      if(locId != -1) {
        triangulation.getEdgeVertex(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[1][cell.id_];
#endif
    } else {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*vectors_)[2][cell.id_]};
      if(locId != -1) {
        triangulation.getEdgeTriangle(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[2][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 2) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*vectors_)[3][cell.id_]};
      if(locId != -1) {
        triangulation.getTriangleEdge(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[3][cell.id_];
#endif
    } else {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*vectors_)[4][cell.id_]};
      if(locId != -1) {
        triangulation.getTriangleStar(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[4][cell.id_];
#endif
    }
  }

  else if(cell.dim_ == 3) {
    if(isReverse) {
#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      const auto locId{(*vectors_)[5][cell.id_]};
      if(locId != -1) {
        triangulation.getCellTriangle(cell.id_, locId, id);
      }
#else
      id = (*vectors_)[5][cell.id_];
#endif
    }
  }

  return id;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::getDescendingPath(
  const Cell &cell,
  std::vector<Cell> &vpath,
  const triangulationType &triangulation,
  const bool stopOnCycle) const {

  const SimplexId numberOfVerts = triangulation.getNumberOfVertices();
  std::vector<char> isCycle;
  isCycle.resize(numberOfVerts, 0);
  int numAlternates{0};

  if(cell.dim_ == 0) {
    // assume that cellId is a vertex
    SimplexId currentId = cell.id_;
    SimplexId connectedEdgeId;
    do {
      // Check if a cycle occured first
      if(isCycle[currentId] == 0) {
        isCycle[currentId] = 1;
      } else {
        if(stopOnCycle || dimensionality_ == 3) {
          break; // Don't trace beyond (not possible in 3D)
        }
        // Otherwise simulate the vpath connected with new saddle.
        if(!reverseFullOrbit) {
          // If desired, remove vpath around the cycle
          while(!(vpath.back() == Cell(0, currentId))) {
            vpath.pop_back();
          }
        }

        if(dimensionality_ == 2) {
          // Can trace the ascending paths (pick and add smallest)
          if(!reverseFullOrbit) {
            const Cell vertex(0, currentId);
            connectedEdgeId = getPairedCell(vertex, triangulation);
          } else {
            connectedEdgeId = vpath.back().id_;
            vpath.pop_back(); // Added back later
          }
          float bestWeight{std::numeric_limits<float>::max()};
          bool firstCellCritical{false};
          std::vector<Cell> vpathBest;
          int bestAlternates{0};
          SimplexId triangleNumber
            = triangulation.getEdgeTriangleNumber(connectedEdgeId);
          for(SimplexId i = 0; i < triangleNumber; ++i) {
            SimplexId triangle;
            std::vector<Cell> newVpath;
            newVpath.emplace_back(
              Cell(1, connectedEdgeId)); // Add saddle to vpath for persistence
                                         // calculation
            std::vector<char> ascIsCycle;
            const SimplexId numberOfCells = triangulation.getNumberOfCells();
            ascIsCycle.resize(numberOfCells, 0);
            triangulation.getEdgeTriangle(connectedEdgeId, i, triangle);
            Cell firstTriangle = Cell(2, triangle);
            int alternates
              = getAscendingPathRecursive<dataType, triangulationType>(
                firstTriangle, newVpath, triangulation, isCycle, ascIsCycle);
            float weight = getPersistence<dataType, triangulationType>(
              newVpath, triangulation);
            if(i == 0) {
              bestWeight = weight;
              vpathBest = newVpath;
              bestAlternates = alternates;
              firstCellCritical = isCellCritical(newVpath.back());
            } else if((!firstCellCritical && isCellCritical(newVpath.back()))
                      || (weight < bestWeight
                          && isCellCritical(newVpath.back()))) {
              bestWeight = weight;
              vpathBest = newVpath;
              bestAlternates = alternates;
              firstCellCritical = true;
            }
          }
          // Update the vpath with new path
          for(Cell newCell : vpathBest) {
            vpath.emplace_back(newCell);
          }
          numAlternates = 1 + bestAlternates;
          break;
        }
      }
      // add a vertex
      const Cell vertex(0, currentId);
      vpath.emplace_back(vertex);

      if(isCellCritical(vertex)) {
        break;
      }

      connectedEdgeId = getPairedCell(vertex, triangulation);
      if(connectedEdgeId == -1) {
        break;
      }

      // add an edge
      const Cell edge(1, connectedEdgeId);
      vpath.emplace_back(edge);

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

  return numAlternates;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::getDescendingPathRecursive(
  const Cell &cell,
  std::vector<Cell> &vpath,
  const triangulationType &triangulation,
  std::vector<char> &previousDescPaths,
  std::vector<char> &previousAscPaths) const {

  const SimplexId numberOfVerts = triangulation.getNumberOfVertices();
  std::vector<char> isCycle;
  isCycle.resize(numberOfVerts, 0);
  int numAlternates{0};

  if(cell.dim_ == 0) {
    // assume that cellId is a vertex
    SimplexId currentId = cell.id_;
    SimplexId connectedEdgeId;
    do {
      // Check if a cycle occured first
      if(isCycle[currentId] == 0) {
        if(previousDescPaths[currentId] == 1) {
          if(vpath.size() == 0) {
            // add a vertex
            const Cell vertex(0, currentId);
            vpath.emplace_back(vertex);
          }
          break; // Already seen in previous recursive loop
        }
        isCycle[currentId] = 1;
        previousDescPaths[currentId] = 1;
      } else {
        if(dimensionality_ == 3) {
          break; // Don't trace beyond (not possible in 3D)
        }
        // Otherwise simulate the vpath connected with new saddle
        if(!reverseFullOrbit) {
          // Remove vpath around the cycle(if desired)
          while(!(vpath.back() == Cell(0, currentId))) {
            vpath.pop_back();
          }
        }

        if(dimensionality_ == 2) {
          // Can trace the ascending paths (pick and add smallest)
          if(!reverseFullOrbit) {
            const Cell vertex(0, currentId);
            connectedEdgeId = getPairedCell(vertex, triangulation);
          } else {
            connectedEdgeId = vpath.back().id_;
            vpath.pop_back(); // Added back later
          }
          float bestWeight{std::numeric_limits<float>::max()};
          bool firstCellCritical{false};
          std::vector<Cell> vpathBest;
          int bestAlternates{0};
          SimplexId triangleNumber
            = triangulation.getEdgeTriangleNumber(connectedEdgeId);
          for(SimplexId i = 0; i < triangleNumber; ++i) {
            SimplexId triangle;
            std::vector<Cell> newVpath;
            newVpath.emplace_back(Cell(1, connectedEdgeId));
            triangulation.getEdgeTriangle(connectedEdgeId, i, triangle);
            Cell firstTriangle = Cell(2, triangle);
            int alternates
              = getAscendingPathRecursive<dataType, triangulationType>(
                firstTriangle, newVpath, triangulation, previousDescPaths,
                previousAscPaths);
            float weight = getPersistence<dataType, triangulationType>(
              newVpath, triangulation);
            if(i == 0) {
              bestWeight = weight;
              vpathBest = newVpath;
              bestAlternates = alternates;
              firstCellCritical = isCellCritical(newVpath.back());
            } else if((!firstCellCritical && isCellCritical(newVpath.back()))
                      || (weight < bestWeight
                          && isCellCritical(newVpath.back()))) {
              bestWeight = weight;
              vpathBest = newVpath;
              bestAlternates = alternates;
              firstCellCritical = true;
            }
          }
          // Update the vpath with new path
          for(Cell newCell : vpathBest) {
            vpath.emplace_back(newCell);
          }
          numAlternates = 1 + bestAlternates;
          break;
        }
      }
      // add a vertex
      const Cell vertex(0, currentId);
      vpath.emplace_back(vertex);

      if(isCellCritical(vertex)) {
        break;
      }

      connectedEdgeId = getPairedCell(vertex, triangulation);
      if(connectedEdgeId == -1) {
        break;
      }

      // add an edge
      const Cell edge(1, connectedEdgeId);
      vpath.emplace_back(edge);

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

  return numAlternates;
}

template <typename triangulationType>
bool DiscreteVectorField::getDescendingPathThroughWall(
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
      vpath->emplace_back(saddle2);
    }

    SimplexId currentId = -1;
    {
      int nconnections = 0;
      for(int i = 0; i < 3; ++i) {
        SimplexId edgeId;
        triangulation.getTriangleEdge(saddle2.id_, i, edgeId);
        if(isVisited[edgeId]) {
          // saddle2 can be adjacent to saddle1 on the wall
          if(isCellCritical(Cell(1, edgeId))) {
            if(vpath != nullptr) {
              vpath->emplace_back(Cell(1, edgeId));
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
        vpath->emplace_back(edge);
      }

      if(isCellCritical(edge)) {
        break;
      }

      const SimplexId connectedTriangleId = getPairedCell(edge, triangulation);

      // add a triangle
      const Cell triangle(2, connectedTriangleId);
      if(vpath != nullptr) {
        vpath->emplace_back(triangle);
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

template <typename dataType, typename triangulationType>
int DiscreteVectorField::getAscendingPath(
  const Cell &cell,
  std::vector<Cell> &vpath,
  const triangulationType &triangulation,
  const bool stopOnCycle) const {

  const SimplexId numberOfCells = triangulation.getNumberOfCells();
  std::vector<char> isCycle;
  isCycle.resize(numberOfCells, 0);
  int numAlternates{0};

  if(dimensionality_ == 2) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {
        oldId = currentId;
        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
          // Found cycle in path, don't add triangle already seen
          if(stopOnCycle) {
            break; // Don't trace beyond cycle
          }
          // Else, simulate finding the next path on the saddle
          if(!reverseFullOrbit) {
            // First, Remove vpath around the cycle
            while(!(vpath.back() == Cell(2, currentId))) {
              vpath.pop_back();
            }
          }

          // Can trace the descending paths (pick and add smallest)
          SimplexId connectedEdgeId;
          if(!reverseFullOrbit) {
            const Cell triangle(2, currentId);
            connectedEdgeId = getPairedCell(triangle, triangulation, true);
          } else {
            connectedEdgeId = vpath.back().id_;
            vpath.pop_back(); // Added back later
          }
          // vpath.emplace_back(Cell(1, connectedEdgeId));//Need to add the
          // saddle skipping over(handled later)
          float bestWeight{std::numeric_limits<float>::max()};
          int bestAlternates{0};
          bool firstCellCritical{false};
          std::vector<Cell> vpathBest;
          for(int i = 0; i < 2; ++i) {
            SimplexId vertex;
            std::vector<Cell> newVpath;
            newVpath.emplace_back(Cell(1, connectedEdgeId));
            const SimplexId numberOfVerts = triangulation.getNumberOfVertices();
            std::vector<char> descIsCycle;
            descIsCycle.resize(numberOfVerts, 0);
            triangulation.getEdgeVertex(connectedEdgeId, i, vertex);
            Cell firstVertex = Cell(0, vertex);
            int alternates
              = getDescendingPathRecursive<dataType, triangulationType>(
                firstVertex, newVpath, triangulation, descIsCycle, isCycle);
            float weight = getPersistence<dataType, triangulationType>(
              newVpath, triangulation);
            if(i == 0) {
              bestWeight = weight;
              bestAlternates = alternates;
              vpathBest = newVpath;
              firstCellCritical = isCellCritical(newVpath.back());
            } else if((!firstCellCritical && isCellCritical(newVpath.back()))
                      || (weight < bestWeight
                          && isCellCritical(newVpath.back()))) {
              bestWeight = weight;
              bestAlternates = alternates;
              vpathBest = newVpath;
              firstCellCritical = true;
            }
          }
          // Update the vpath with new path
          for(Cell newCell : vpathBest) {
            vpath.emplace_back(newCell);
          }
          numAlternates = 1 + bestAlternates;
          break;
        }

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.emplace_back(triangle);

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
        vpath.emplace_back(edge);

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

        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
          break; // Cycle detected
        }

        oldId = currentId;

        // add a tetra
        const Cell tetra(3, currentId);
        vpath.emplace_back(tetra);

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
        vpath.emplace_back(triangle);

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

  return numAlternates;
}

template <typename dataType, typename triangulationType>
int DiscreteVectorField::getAscendingPathRecursive(
  const Cell &cell,
  std::vector<Cell> &vpath,
  const triangulationType &triangulation,
  std::vector<char> &previousDescPaths,
  std::vector<char> &previousAscPaths) const {

  const SimplexId numberOfCells = triangulation.getNumberOfCells();
  std::vector<char> isCycle;
  isCycle.resize(numberOfCells, 0);
  int numAlternates{0};

  if(dimensionality_ == 2) {
    if(cell.dim_ == 2) {
      // assume that cellId is a triangle
      SimplexId currentId = cell.id_;
      SimplexId oldId;
      do {
        oldId = currentId;
        if(isCycle[currentId] == 0) {
          if(previousAscPaths[currentId] == 1) {
            if(vpath.size() == 0) {
              // add a triangle
              const Cell triangle(2, currentId);
              vpath.emplace_back(triangle);
            }
            break; // Already traced path previously
          }
          isCycle[currentId] = 1;
          previousAscPaths[currentId] = 1;
        } else {
          // Found cycle in path, don't add triangle already seen

          // Simulate finding the next path on the saddle
          // First, Remove vpath around the cycle(if we don't want to trace it)
          if(!reverseFullOrbit) {
            while(!(vpath.back() == Cell(2, currentId))) {
              vpath.pop_back();
            }
          }

          // Can trace the descending paths (pick and add smallest)
          SimplexId connectedEdgeId;
          if(!reverseFullOrbit) {
            const Cell triangle(2, currentId);
            connectedEdgeId = getPairedCell(triangle, triangulation, true);
          } else {
            connectedEdgeId = vpath.back().id_;
            vpath.pop_back();
          }
          float bestWeight{std::numeric_limits<float>::max()};
          int bestAlternates{0};
          bool firstCellCritical{false};
          std::vector<Cell> vpathBest;
          for(int i = 0; i < 2; ++i) {
            SimplexId vertex;
            std::vector<Cell> newVpath;
            newVpath.emplace_back(Cell(1, connectedEdgeId));
            triangulation.getEdgeVertex(connectedEdgeId, i, vertex);
            Cell firstVertex = Cell(0, vertex);
            int alternates
              = getDescendingPathRecursive<dataType, triangulationType>(
                firstVertex, newVpath, triangulation, previousDescPaths,
                previousAscPaths);
            float weight = getPersistence<dataType, triangulationType>(
              newVpath, triangulation);
            if(i == 0) {
              bestWeight = weight;
              bestAlternates = alternates;
              vpathBest = newVpath;
              firstCellCritical = isCellCritical(newVpath.back());
            } else if((!firstCellCritical && isCellCritical(newVpath.back()))
                      || (weight < bestWeight
                          && isCellCritical(newVpath.back()))) {
              bestWeight = weight;
              bestAlternates = alternates;
              vpathBest = newVpath;
              firstCellCritical = true;
            }
          }
          // Update the vpath with new path
          for(Cell newCell : vpathBest) {
            vpath.emplace_back(newCell);
          }
          numAlternates = 1 + bestAlternates;
          break;
        }

        // add a triangle
        const Cell triangle(2, currentId);
        vpath.emplace_back(triangle);

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
        vpath.emplace_back(edge);

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

        if(isCycle[currentId] == 0) {
          isCycle[currentId] = 1;
        } else {
          break; // Cycle detected
        }

        oldId = currentId;

        // add a tetra
        const Cell tetra(3, currentId);
        vpath.emplace_back(tetra);

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
        vpath.emplace_back(triangle);

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

  return numAlternates;
}

template <typename triangulationType>
void DiscreteVectorField::getAscendingPathThroughWall(
  const Cell &saddle1,
  const Cell &saddle2,
  const std::vector<bool> &isVisited,
  std::vector<Cell> *const vpath,
  const triangulationType &triangulation) const {

  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  std::vector<char> isCycle;
  isCycle.resize(numberOfTriangles, 0);
  std::vector<SimplexId> alternatePathStack; // Keep track of alternate paths

  if(dimensionality_ == 3) {
    // add the 1-saddle to the path
    if(vpath != nullptr) {
      vpath->emplace_back(saddle1);
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
          if(isCellCritical(Cell(2, triangleId))) {
            if(vpath != nullptr) {
              vpath->emplace_back(Cell(2, triangleId));
            }
            return;
          }
          if(nconnections < 1) {
            currentId = triangleId;
            ++nconnections;
          } else {
            alternatePathStack.emplace_back(triangleId);
          }
        }
      }
    }

    if(currentId == -1) { // Shouldn't happen
      this->printErr(
        "Current ID not updated for getAscendingPathThroughWall()");
      return;
    }

    SimplexId oldId;
    do {

      if(isCycle[currentId] == 0) {
        isCycle[currentId] = 1;
      } else {
        // Probably can't get here
        this->printErr("Cycle detected on the wall of 2-saddle "
                       + std::to_string(saddle2.id_));
        return;
      }

      oldId = currentId;

      // add a triangle
      const Cell triangle(2, currentId);
      if(vpath != nullptr) {
        vpath->emplace_back(triangle);
      }

      if(isCellCritical(triangle)) {
        break;
      }

      const SimplexId connectedEdgeId
        = getPairedCell(triangle, triangulation, true);

      // add an edge
      const Cell edge(1, connectedEdgeId);
      if(vpath != nullptr) {
        vpath->emplace_back(edge);
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

        if(isVisited[triangleId] and triangleId != oldId
           and isCycle[triangleId] == 0) {
          if(nconnections < 1) {
            currentId = triangleId;
            ++nconnections;
          } else {
            alternatePathStack.emplace_back(triangleId);
          }
        }
      }
      if(nconnections == 0) {
        // No path found, trace back using the stack
        if(!alternatePathStack.empty()) {
          currentId = alternatePathStack.back();
          alternatePathStack.pop_back();
        } else { // Should be impossible (assuming a direct path exists to
                 // saddle2)
          this->printErr(
            "No alternate paths detected for getAscendingPathThroughWall()");
          return;
        }
      }

      // stop at convergence caused by boundary effect
    } while(currentId != oldId);
  }

  return;
}

template <typename triangulationType>
int DiscreteVectorField::getDescendingWall(
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
            wall->emplace_back(Cell(2, triangleId));
          }

          for(int j = 0; j < 3; ++j) {
            SimplexId edgeId;
            triangulation.getTriangleEdge(triangleId, j, edgeId);

            if((saddles != nullptr) and isCellCritical(Cell(1, edgeId))) {
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
int DiscreteVectorField::getAscendingWall(
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
            wall->emplace_back(Cell(1, edgeId));
          }

          const SimplexId triangleNumber
            = triangulation.getEdgeTriangleNumber(edgeId);
          for(SimplexId j = 0; j < triangleNumber; ++j) {
            SimplexId triangleId;
            triangulation.getEdgeTriangle(edgeId, j, triangleId);

            if((saddles != nullptr) and isCellCritical(Cell(2, triangleId))) {
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
int DiscreteVectorField::reverseAscendingPath(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  if(dimensionality_ == 2) {
    // assume that the first cell is an edge
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId edgeId = vpath[i].id_;

      if((i + 1) == numberOfCellsInPath) {
        TTK_FORCE_USE(triangulation);
        (*vectors_)[2][edgeId] = NULL_CONNECTION;
        break;
      } // Handle cycle case
      const SimplexId triangleId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 3; ++k) {
        SimplexId tmp;
        triangulation.getCellEdge(triangleId, k, tmp);
        if(tmp == edgeId) {
          (*vectors_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeStarNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeStar(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*vectors_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*vectors_)[3][triangleId] = edgeId;
      (*vectors_)[2][edgeId] = triangleId;
#endif
    }
  } else if(dimensionality_ == 3) {
    // assume that the first cell is a triangle
    const SimplexId numberOfCellsInPath = vpath.size();
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      const SimplexId triangleId = vpath[i].id_;
      if((i + 1) == numberOfCellsInPath) {
        (*vectors_)[4][triangleId] = NULL_CONNECTION;
        break;
      } // Handle cycle case
      const SimplexId tetraId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
      for(int k = 0; k < 4; ++k) {
        SimplexId tmp;
        triangulation.getCellTriangle(tetraId, k, tmp);
        if(tmp == triangleId) {
          (*vectors_)[5][tetraId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getTriangleStarNumber(triangleId); ++k) {
        SimplexId tmp;
        triangulation.getTriangleStar(triangleId, k, tmp);
        if(tmp == tetraId) {
          (*vectors_)[4][triangleId] = k;
          break;
        }
      }
#else
      (*vectors_)[5][tetraId] = triangleId;
      (*vectors_)[4][triangleId] = tetraId;
#endif
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteVectorField::reverseDescendingPath(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {

  // assume that the first cell is an edge
  for(size_t i = 0; i < vpath.size(); i += 2) {
    const SimplexId edgeId = vpath[i].id_;
    if((i + 1) == vpath.size()) {
      TTK_FORCE_USE(triangulation);
      (*vectors_)[1][edgeId] = NULL_CONNECTION;
      continue;
    } // Handle cycle case
    const SimplexId vertId = vpath[i + 1].id_;

#ifdef TTK_ENABLE_DCG_OPTIMIZE_MEMORY
    const auto nneighs = triangulation.getVertexEdgeNumber();
    for(int k = 0; k < nneighs; ++k) {
      SimplexId tmp;
      triangulation.getVertexEdge(vertId, k, tmp);
      if(tmp == edgeId) {
        (*vectors_)[0][vertId] = k;
        break;
      }
    }
    const auto nverts = triangulation.getEdgeStarNumber(edgeId);
    for(int k = 0; k < nverts; ++k) {
      SimplexId tmp;
      triangulation.getEdgeVertex(edgeId, k, tmp);
      if(tmp == vertId) {
        (*vectors_)[1][edgeId] = k;
        break;
      }
    }
#else
    TTK_FORCE_USE(triangulation);
    (*vectors_)[0][vertId] = edgeId;
    (*vectors_)[1][edgeId] = vertId;
#endif
  }

  return 0;
}

template <typename triangulationType>
int DiscreteVectorField::reverseAlternatingPath(
  const std::vector<Cell> &vpath,
  const triangulationType &triangulation) const {
  if(dimensionality_ == 2) {
    const SimplexId numberOfCellsInPath = vpath.size();
    std::vector<Cell> currentVpath;
    int currentDim{-1};
    if(numberOfCellsInPath > 1) {
      currentDim = vpath[1].dim_;
    } else {
      this->printErr("Rotating vpath not large enough to call reverse");
    }
    // Assume the first cell is always an edge
    for(SimplexId i = 0; i < numberOfCellsInPath; i += 2) {
      currentVpath.emplace_back(vpath[i]);
      if(currentDim == vpath[i + 1].dim_) {
        currentVpath.emplace_back(vpath[i + 1]);
      } else {
        if(currentDim == 0) {
          reverseDescendingPath(currentVpath, triangulation);
        } else {
          reverseAscendingPath(currentVpath, triangulation);
        }
        currentVpath.clear();
        currentVpath.emplace_back(vpath[i]);
        currentDim = vpath[i + 1].dim_;
        currentVpath.emplace_back(vpath[i + 1]);
      }
    }
    if(currentDim == 0) {
      reverseDescendingPath(currentVpath, triangulation);
    } else {
      reverseAscendingPath(currentVpath, triangulation);
    }
  }
  return 0;
}

template <typename triangulationType>
int DiscreteVectorField::reverseAscendingPathOnWall(
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
          (*vectors_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeTriangleNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*vectors_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*vectors_)[3][triangleId] = edgeId;
      (*vectors_)[2][edgeId] = triangleId;
#endif
    }
  }

  return 0;
}

template <typename triangulationType>
int DiscreteVectorField::reverseDescendingPathOnWall(
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
          (*vectors_)[3][triangleId] = k;
          break;
        }
      }
      for(int k = 0; k < triangulation.getEdgeTriangleNumber(edgeId); ++k) {
        SimplexId tmp;
        triangulation.getEdgeTriangle(edgeId, k, tmp);
        if(tmp == triangleId) {
          (*vectors_)[2][edgeId] = k;
          break;
        }
      }
#else
      TTK_FORCE_USE(triangulation);
      (*vectors_)[2][edgeId] = triangleId;
      (*vectors_)[3][triangleId] = edgeId;
#endif
    }
  }

  return 0;
}

template <typename dataType, typename triangulationType>
ttk::SimplexId DiscreteVectorField::getCellGreaterVertex(
  const Cell c, const triangulationType &triangulation) const {

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
    float edgeWeight;
    if(compare<dataType, triangulationType>(
         triangulation, v0, v1, edgeWeight)) {
      vertexId = v0;
    } else {
      vertexId = v1;
    }
  }

  else if(cellDim == 2) {
    SimplexId v0, v1, v2;
    triangulation.getTriangleVertex(cellId, 0, v0);
    triangulation.getTriangleVertex(cellId, 1, v1);
    triangulation.getTriangleVertex(cellId, 2, v2);
    float edgeWeight;
    bool hasV0For1
      = compare<dataType, triangulationType>(triangulation, v0, v1, edgeWeight);
    bool hasV0For2
      = compare<dataType, triangulationType>(triangulation, v0, v2, edgeWeight);
    bool hasV1For2
      = compare<dataType, triangulationType>(triangulation, v1, v2, edgeWeight);
    if(hasV0For1 && hasV0For2) {
      vertexId = v0;
    } else if(!hasV0For1 && hasV1For2) {
      vertexId = v1;
    } else if(!hasV0For2 && !hasV1For2) {
      vertexId = v2;
    } else {
      vertexId = -1; // No clear greater vertex
    }
  }

  else if(cellDim == 3) {
    SimplexId v0, v1, v2, v3;
    triangulation.getCellVertex(cellId, 0, v0);
    triangulation.getCellVertex(cellId, 1, v1);
    triangulation.getCellVertex(cellId, 2, v2);
    triangulation.getCellVertex(cellId, 3, v3);
    float edgeWeight;
    bool hasV0For1
      = compare<dataType, triangulationType>(triangulation, v0, v1, edgeWeight);
    bool hasV0For2
      = compare<dataType, triangulationType>(triangulation, v0, v2, edgeWeight);
    bool hasV0For3
      = compare<dataType, triangulationType>(triangulation, v0, v3, edgeWeight);
    bool hasV1For2
      = compare<dataType, triangulationType>(triangulation, v1, v2, edgeWeight);
    bool hasV1For3
      = compare<dataType, triangulationType>(triangulation, v1, v3, edgeWeight);
    bool hasV2For3
      = compare<dataType, triangulationType>(triangulation, v2, v3, edgeWeight);
    if(hasV0For1 && hasV0For2 && hasV0For3) {
      vertexId = v0;
    } else if(!hasV0For1 && hasV1For2 && hasV1For3) {
      vertexId = v1;
    } else if(!hasV0For2 && !hasV1For2 && hasV2For3) {
      vertexId = v2;
    } else if(!hasV0For3 && !hasV1For3 && !hasV2For3) {
      vertexId = v3;
    } else {
      vertexId = -1; // No clear greater vertex
    }
  }
  return vertexId;
}

template <typename dataType, typename triangulationType>
ttk::SimplexId DiscreteVectorField::getCellLowerVertex(
  const Cell c, const triangulationType &triangulation) const {

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
    float edgeWeight;
    if(compare<dataType, triangulationType>(
         triangulation, v1, v0, edgeWeight)) {
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
    float edgeWeight;
    bool hasV0For1
      = compare<dataType, triangulationType>(triangulation, v0, v1, edgeWeight);
    bool hasV0For2
      = compare<dataType, triangulationType>(triangulation, v0, v2, edgeWeight);
    bool hasV1For2
      = compare<dataType, triangulationType>(triangulation, v1, v2, edgeWeight);
    if(!(hasV0For1) && !(hasV0For2)) {
      vertexId = v0;
    } else if(hasV0For1 && !hasV1For2) {
      vertexId = v1;
    } else if(hasV1For2 && hasV0For2) {
      vertexId = v2;
    } else {
      vertexId = -1; // No clear lower vertex
    }
  }

  else if(cellDim == 3) {
    SimplexId v0{}, v1{}, v2{}, v3{};
    triangulation.getCellVertex(cellId, 0, v0);
    triangulation.getCellVertex(cellId, 1, v1);
    triangulation.getCellVertex(cellId, 2, v2);
    triangulation.getCellVertex(cellId, 3, v3);
    float edgeWeight;
    bool hasV0For1
      = compare<dataType, triangulationType>(triangulation, v0, v1, edgeWeight);
    bool hasV0For2
      = compare<dataType, triangulationType>(triangulation, v0, v2, edgeWeight);
    bool hasV0For3
      = compare<dataType, triangulationType>(triangulation, v0, v3, edgeWeight);
    bool hasV1For2
      = compare<dataType, triangulationType>(triangulation, v1, v2, edgeWeight);
    bool hasV1For3
      = compare<dataType, triangulationType>(triangulation, v1, v3, edgeWeight);
    bool hasV2For3
      = compare<dataType, triangulationType>(triangulation, v2, v3, edgeWeight);
    if(!hasV0For1 && !hasV0For2 && !hasV0For3) {
      vertexId = v0;
    } else if(hasV0For1 && !hasV1For2 && !hasV1For3) {
      vertexId = v1;
    } else if(hasV0For2 && hasV1For2 && !hasV2For3) {
      vertexId = v2;
    } else if(hasV0For3 && hasV1For3 && hasV2For3) {
      vertexId = v3;
    } else {
      vertexId = -1; // No clear lower vertex
    }
  }
  return vertexId;
}

template <typename triangulationType>
int DiscreteVectorField::setVectorGlyphs(
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
