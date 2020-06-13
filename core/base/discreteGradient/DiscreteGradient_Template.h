/// \ingroup baseCode
/// \class ttk::DiscreteGradient
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
using ttk::dcg::CriticalPoint;
using ttk::dcg::DiscreteGradient;
using ttk::dcg::SaddleSaddleVPathComparator;
using ttk::dcg::VPath;

template <typename dataType>
dataType DiscreteGradient::scalarMax(const Cell &cell,
                                     const dataType *const scalars) const {
  dataType scalar{};

  if(dimensionality_ == 2) {
    switch(cell.dim_) {
      case 0:
        scalar = scalars[cell.id_];
        break;

      case 1:
        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId = -1;
          inputTriangulation_->getEdgeVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar < vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId vertexId = -1;
          inputTriangulation_->getCellVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar < vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cell.dim_) {
      case 0:
        scalar = scalars[cell.id_];
        break;

      case 1:
        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId = -1;
          inputTriangulation_->getEdgeVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar < vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId vertexId = -1;
          inputTriangulation_->getTriangleVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar < vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 3:
        for(int i = 0; i < 4; ++i) {
          SimplexId vertexId = -1;
          inputTriangulation_->getCellVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar < vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;
    }
  }

  return scalar;
}

template <typename dataType>
dataType DiscreteGradient::scalarMin(const Cell &cell,
                                     const dataType *const scalars) const {
  dataType scalar{};

  if(dimensionality_ == 2) {
    switch(cell.dim_) {
      case 0:
        scalar = scalars[cell.id_];
        break;

      case 1:
        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar > vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getCellVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar > vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cell.dim_) {
      case 0:
        scalar = scalars[cell.id_];
        break;

      case 1:
        for(int i = 0; i < 2; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar > vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 2:
        for(int i = 0; i < 3; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getTriangleVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar > vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;

      case 3:
        for(int i = 0; i < 4; ++i) {
          SimplexId vertexId;
          inputTriangulation_->getCellVertex(cell.id_, i, vertexId);
          const dataType vertexScalar = scalars[vertexId];

          if(!i or scalar > vertexScalar) {
            scalar = vertexScalar;
          }
        }
        break;
    }
  }

  return scalar;
}

template <typename dataType>
dataType DiscreteGradient::getPersistence(const Cell &up,
                                          const Cell &down,
                                          const dataType *const scalars) const {
  return scalarMax<dataType>(up, scalars) - scalarMin<dataType>(down, scalars);
}

template <typename dataType, typename idType>
int DiscreteGradient::buildGradient() {
  Timer t;

  const auto *const offsets = static_cast<const idType *>(inputOffsets_);
  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  const int numberOfDimensions = getNumberOfDimensions();

  // init number of cells by dimension
  std::vector<SimplexId> numberOfCells(numberOfDimensions);
  for(int i = 0; i < numberOfDimensions; ++i) {
    numberOfCells[i] = getNumberOfCells(i);
  }

  dmtMax2PL_.clear();
  dmt1Saddle2PL_.clear();
  dmt2Saddle2PL_.clear();
  gradient_.clear();
  gradient_.resize(dimensionality_);
  for(int i = 0; i < dimensionality_; ++i) {
    // init gradient memory
    gradient_[i].resize(numberOfDimensions);
    gradient_[i][i].resize(numberOfCells[i], -1);
    gradient_[i][i + 1].resize(numberOfCells[i + 1], -1);
  }

  sortVertices(numberOfCells[0], vertsOrder_, scalars, offsets);

  // compute gradient pairs
  processLowerStars();

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Data-set: " << numberOfVertices_ << " v., "
        << inputTriangulation_->getNumberOfEdges() << " e.";
    if(inputTriangulation_->getDimensionality() == 3) {
      msg << ", " << inputTriangulation_->getNumberOfTriangles() << " t., "
          << inputTriangulation_->getNumberOfCells() << " T." << std::endl;
    } else if(inputTriangulation_->getDimensionality() == 2) {
      msg << ", " << inputTriangulation_->getNumberOfCells() << " t."
          << std::endl;
    }

    msg << "[DiscreteGradient] Processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints(
  const std::vector<Cell> &criticalPoints,
  std::vector<size_t> &nCriticalPointsByDim) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputCriticalPoints_numberOfPoints_) {
    std::cerr << "[DiscreteGradient] critical points' pointer to "
                 "numberOfPoints is null."
              << std::endl;
    return -1;
  }
  if(!outputCriticalPoints_points_) {
    std::cerr
      << "[DiscreteGradient] critical points' pointer to points is null."
      << std::endl;
    return -1;
  }
  if(!inputScalarField_) {
    std::cerr << "[DiscreteGradient] critical points' pointer to the input "
                 "scalar field is null."
              << std::endl;
    return -1;
  }
#endif
  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);
  auto *outputCriticalPoints_points_cellScalars
    = static_cast<std::vector<dataType> *>(
      outputCriticalPoints_points_cellScalars_);

  const auto nCritPoints = criticalPoints.size();
  (*outputCriticalPoints_numberOfPoints_) = nCritPoints;

  const int numberOfDimensions = getNumberOfDimensions();
  nCriticalPointsByDim.resize(numberOfDimensions, 0);

  // sequential loop over critical points
  for(size_t i = 0; i < nCritPoints; ++i) {
    const Cell &cell = criticalPoints[i];
    nCriticalPointsByDim[cell.dim_]++;
  }

  outputCriticalPoints_points_->resize(3 * nCritPoints);
  if(outputCriticalPoints_points_cellDimensions_) {
    outputCriticalPoints_points_cellDimensions_->resize(nCritPoints);
  }
  if(outputCriticalPoints_points_cellIds_) {
    outputCriticalPoints_points_cellIds_->resize(nCritPoints);
  }
  if(outputCriticalPoints_points_cellScalars) {
    outputCriticalPoints_points_cellScalars->resize(nCritPoints);
  }
  if(outputCriticalPoints_points_isOnBoundary_) {
    outputCriticalPoints_points_isOnBoundary_->resize(nCritPoints);
  }
  if(outputCriticalPoints_points_PLVertexIdentifiers_) {
    outputCriticalPoints_points_PLVertexIdentifiers_->resize(nCritPoints);
  }

  // for all critical cells
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCritPoints; ++i) {
    const Cell &cell = criticalPoints[i];
    const int cellDim = cell.dim_;
    const SimplexId cellId = cell.id_;

    float incenter[3];
    inputTriangulation_->getCellIncenter(cell.id_, cell.dim_, incenter);

    const auto scalar = scalarMax<dataType>(cell, scalars);
    const char isOnBoundary = isBoundary(cell);

    (*outputCriticalPoints_points_)[3 * i] = incenter[0];
    (*outputCriticalPoints_points_)[3 * i + 1] = incenter[1];
    (*outputCriticalPoints_points_)[3 * i + 2] = incenter[2];

    if(outputCriticalPoints_points_cellDimensions_) {
      (*outputCriticalPoints_points_cellDimensions_)[i] = cellDim;
    }
    if(outputCriticalPoints_points_cellIds_) {
      (*outputCriticalPoints_points_cellIds_)[i] = cellId;
    }
    if(outputCriticalPoints_points_cellScalars) {
      (*outputCriticalPoints_points_cellScalars)[i] = scalar;
    }
    if(outputCriticalPoints_points_isOnBoundary_) {
      (*outputCriticalPoints_points_isOnBoundary_)[i] = isOnBoundary;
    }
    if(outputCriticalPoints_points_PLVertexIdentifiers_) {
      auto vertId = getCellGreaterVertex(cell);
      (*outputCriticalPoints_points_PLVertexIdentifiers_)[i] = vertId;
    }
  }

  {
    std::stringstream msg;
    for(int i = 0; i < numberOfDimensions; ++i) {
      msg << "[DiscreteGradient] " << nCriticalPointsByDim[i] << " " << i
          << "-cell(s)." << std::endl;
    }
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints() const {
  std::vector<Cell> criticalPoints;
  getCriticalPoints(criticalPoints);
  std::vector<size_t> nCriticalPointsByDim{};
  setCriticalPoints<dataType>(criticalPoints, nCriticalPointsByDim);

  return 0;
}

template <typename dataType>
int DiscreteGradient::getRemovableMaxima(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const bool allowBoundary,
  std::vector<char> &isRemovableMaximum,
  std::vector<SimplexId> &pl2dmt_maximum) {
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  const SimplexId numberOfCells = inputTriangulation_->getNumberOfCells();
  const int maximumDim = dimensionality_;

  // Detect DMT-max cells to remove
  isRemovableMaximum.resize(numberOfCells);

  dmtMax2PL_.resize(numberOfCells);
  std::fill(dmtMax2PL_.begin(), dmtMax2PL_.end(), -1);

  // by default : maximum is removable
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfCells; ++i) {
    const Cell maximumCandidate(maximumDim, i);
    isRemovableMaximum[i] = isMaximum(maximumCandidate);
  }

  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const std::pair<SimplexId, char> &criticalPoint = criticalPoints[i];
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    if(criticalPointType == static_cast<char>(CriticalType::Local_maximum)) {
      if(!allowBoundary
         and inputTriangulation_->isVertexOnBoundary(criticalPointId)) {
        continue;
      }

      SimplexId numberOfMaxima = 0;
      SimplexId maximumId = -1;
      const SimplexId starNumber
        = inputTriangulation_->getVertexStarNumber(criticalPointId);
      for(SimplexId j = 0; j < starNumber; ++j) {
        SimplexId starId;
        inputTriangulation_->getVertexStar(criticalPointId, j, starId);

        if(isMaximum(Cell(maximumDim, starId)) and dmtMax2PL_[starId] == -1) {
          maximumId = starId;
          ++numberOfMaxima;
        }
      }

      // a DMT-maximum in the star of only one PL-maximum cannot be removed
      // and is automatically associated to it.
      if(numberOfMaxima == 1) {
        if(dmtMax2PL_[maximumId] == -1
           and pl2dmt_maximum[criticalPointId] == -1) {
          dmtMax2PL_[maximumId] = criticalPointId;
          pl2dmt_maximum[criticalPointId] = maximumId;
          isRemovableMaximum[maximumId] = false;
        }
      }
    }
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::getRemovableSaddles1(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const bool allowBoundary,
  std::vector<char> &isRemovableSaddle,
  std::vector<SimplexId> &pl2dmt_saddle) {
  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  isRemovableSaddle.resize(numberOfEdges);

  dmt1Saddle2PL_.resize(numberOfEdges);
  std::fill(dmt1Saddle2PL_.begin(), dmt1Saddle2PL_.end(), -1);

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
      if(!allowBoundary
         and inputTriangulation_->isVertexOnBoundary(criticalPointId)) {
        continue;
      }

      SimplexId numberOfSaddles = 0;
      SimplexId saddleId = -1;
      const SimplexId edgeNumber
        = inputTriangulation_->getVertexEdgeNumber(criticalPointId);
      for(SimplexId i = 0; i < edgeNumber; ++i) {
        SimplexId edgeId;
        inputTriangulation_->getVertexEdge(criticalPointId, i, edgeId);
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

template <typename dataType>
int DiscreteGradient::getRemovableSaddles2(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const bool allowBoundary,
  std::vector<char> &isRemovableSaddle,
  std::vector<SimplexId> &pl2dmt_saddle) {
  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
  isRemovableSaddle.resize(numberOfTriangles);

  dmt2Saddle2PL_.resize(numberOfTriangles);
  std::fill(dmt2Saddle2PL_.begin(), dmt2Saddle2PL_.end(), -1);

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
      if(!allowBoundary
         and inputTriangulation_->isVertexOnBoundary(criticalPointId)) {
        continue;
      }

      SimplexId numberOfSaddles = 0;
      SimplexId saddleId = -1;
      const SimplexId triangleNumber
        = inputTriangulation_->getVertexTriangleNumber(criticalPointId);
      for(SimplexId i = 0; i < triangleNumber; ++i) {
        SimplexId triangleId;
        inputTriangulation_->getVertexTriangle(criticalPointId, i, triangleId);
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

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections1(
  const std::vector<char> &isRemovableSaddle1,
  const std::vector<char> &isRemovableSaddle2,
  const bool allowBruteForce,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index) const {
  Timer t;

  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  const int maximumDim = dimensionality_;
  const int saddle2Dim = maximumDim - 1;
  const int saddle1Dim = saddle2Dim - 1;

  // Part 1 : build initial structures
  // add the 2-saddles to CriticalPointList
  const SimplexId numberOfSaddle2Candidates = getNumberOfCells(saddle2Dim);
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
  const SimplexId numberOfSaddle1Candidates = getNumberOfCells(saddle1Dim);
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

    std::set<SimplexId> saddles1;
    getDescendingWall(saddle2, mask, nullptr, &saddles1);

    for(auto &saddle1Id : saddles1) {
      if(!isRemovableSaddle1[saddle1Id]) {
        continue;
      }

      const Cell &saddle1 = Cell(1, saddle1Id);

      std::vector<Cell> path;
      const bool isMultiConnected
        = getAscendingPathThroughWall(saddle1, saddle2, isVisited, &path, true);

      if(!isMultiConnected) {
        const SimplexId sourceIndex = saddle1Index[saddle1Id];
        CriticalPoint &source = criticalPoints[sourceIndex];

        // update source and destination
        const SimplexId sourceSlot = source.addSlot();
        const SimplexId destinationSlot = destination.addSlot();

        // update vpath
        const auto persistence
          = getPersistence<dataType>(saddle2, saddle1, scalars);

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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime()
        << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t"
        << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
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
  std::vector<SimplexId> &saddle2Index) {
  Timer t;

  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
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

      std::set<SimplexId> saddles1;
      VisitedMask mask{isVisited, visitedCells};
      getDescendingWall(minSaddle2, mask, nullptr, &saddles1);

      // check if at least one connection exists
      auto isFound = saddles1.find(minSaddle1.id_);
      if(isFound == saddles1.end()) {
        ++numberOfIterations;
        continue;
      }

      // check if there is multiple connections
      std::vector<Cell> path;
      const bool isMultiConnected = getAscendingPathThroughWall(
        minSaddle1, minSaddle2, isVisited, &path, true);

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
            inputTriangulation_->getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle1)) {
              continue;
            }

            if(!allowBoundary
               and inputTriangulation_->isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle1[vertexId] == -1) {
              const SimplexId pl_saddle1Id = vertexId;

              SimplexId numberOfRemainingSaddles1 = 0;

              SimplexId savedId = -1;
              const SimplexId edgeNumber
                = inputTriangulation_->getVertexEdgeNumber(pl_saddle1Id);
              for(SimplexId j = 0; j < edgeNumber; ++j) {
                SimplexId edgeId;
                inputTriangulation_->getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId != dmt_saddle1Id and isSaddle1(Cell(1, edgeId))
                   and isRemovableSaddle1[edgeId]) {
                  ++numberOfRemainingSaddles1;
                  savedId = edgeId;
                }
              }

              if(numberOfRemainingSaddles1 == 0) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                pl2dmt_saddle1[vertexId] = dmt_saddle1Id;
                dmt1Saddle2PL_[dmt_saddle1Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1 == 1) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                isRemovableSaddle1[savedId] = false;
                pl2dmt_saddle1[vertexId] = savedId;
                dmt1Saddle2PL_[savedId] = vertexId;
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
            inputTriangulation_->getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle2)) {
              continue;
            }

            if(!allowBoundary
               and inputTriangulation_->isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle2[vertexId] == -1) {
              const SimplexId pl_saddle2Id = vertexId;

              SimplexId numberOfRemainingSaddles2 = 0;

              SimplexId savedId = -1;
              const SimplexId triangleNumber
                = inputTriangulation_->getVertexTriangleNumber(pl_saddle2Id);
              for(SimplexId j = 0; j < triangleNumber; ++j) {
                SimplexId triangleId;
                inputTriangulation_->getVertexTriangle(
                  pl_saddle2Id, j, triangleId);

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
                dmt2Saddle2PL_[dmt_saddle2Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles2 == 1) {
                isRemovableSaddle2[dmt_saddle2Id] = false;
                isRemovableSaddle2[savedId] = false;
                pl2dmt_saddle2[vertexId] = savedId;
                dmt2Saddle2PL_[savedId] = vertexId;
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
        reverseAscendingPathOnWall(path);
      }
    }

    if(vpath.isValid_) {
      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_) {
        const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
        const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back({minSaddle1, minSaddle2});
      }

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

        std::set<SimplexId> saddles1;
        VisitedMask mask{isVisited, visitedCells};
        getDescendingWall(saddle2, mask, nullptr, &saddles1);

        for(auto &saddle1Id : saddles1) {
          const Cell saddle1(1, saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getAscendingPathThroughWall(
            saddle1, saddle2, isVisited, &path, true);

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
            = getPersistence<dataType>(saddle2, saddle1, scalars);

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

        std::set<SimplexId> saddles2;
        VisitedMask mask{isVisited, visitedCells};
        getAscendingWall(saddle1, mask, nullptr, &saddles2);

        for(auto &saddle2Id : saddles2) {
          const Cell saddle2(2, saddle2Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getDescendingPathThroughWall(
            saddle2, saddle1, isVisited, &path, true);

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
            = getPersistence<dataType>(saddle2, saddle1, scalars);

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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t"
        << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections1(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::vector<char> &isPL,
  const int iterationThreshold,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors) {
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<SimplexId> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(
    criticalPoints, allowBoundary, isRemovableSaddle1, pl2dmt_saddle1);

  std::vector<char> isRemovableSaddle2;
  std::vector<SimplexId> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(
    criticalPoints, allowBoundary, isRemovableSaddle2, pl2dmt_saddle2);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<SimplexId> saddle1Index;
  std::vector<SimplexId> saddle2Index;
  initializeSaddleSaddleConnections1<dataType>(
    isRemovableSaddle1, isRemovableSaddle2, allowBruteForce, vpaths,
    dmt_criticalPoints, saddle1Index, saddle2Index);

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
    saddle1Index, saddle2Index);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle-Saddle pairs simplified in "
        << t.getElapsedTime() << " s, " << threadNumber_ << " thread(s)."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections2(
  const std::vector<char> &isRemovableSaddle1,
  const std::vector<char> &isRemovableSaddle2,
  const bool allowBruteForce,
  std::vector<VPath> &vpaths,
  std::vector<CriticalPoint> &criticalPoints,
  std::vector<SimplexId> &saddle1Index,
  std::vector<SimplexId> &saddle2Index) const {
  Timer t;

  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  const int maximumDim = dimensionality_;
  const int saddle2Dim = maximumDim - 1;
  const int saddle1Dim = saddle2Dim - 1;

  // Part 1 : build initial structures
  // add the 1-saddles to CriticalPointList
  const SimplexId numberOfSaddle1Candidates = getNumberOfCells(saddle1Dim);
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
  const SimplexId numberOfSaddle2Candidates = getNumberOfCells(saddle2Dim);
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

    std::set<SimplexId> saddles2;
    VisitedMask mask{isVisited, visitedEdges};
    getAscendingWall(saddle1, mask, nullptr, &saddles2);

    for(auto &saddle2Id : saddles2) {
      if(!isRemovableSaddle2[saddle2Id]) {
        continue;
      }

      const Cell &saddle2 = Cell(2, saddle2Id);

      std::vector<Cell> path;
      const bool isMultiConnected = getDescendingPathThroughWall(
        saddle2, saddle1, isVisited, &path, true);

      if(!isMultiConnected) {
        const SimplexId destinationIndex = saddle2Index[saddle2Id];
        CriticalPoint &destination = criticalPoints[destinationIndex];

        // update source and destination
        const SimplexId sourceSlot = source.addSlot();
        const SimplexId destinationSlot = destination.addSlot();

        // update vpath
        const auto persistence
          = getPersistence<dataType>(saddle2, saddle1, scalars);

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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime()
        << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t"
        << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleSaddleConnections2(
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
  std::vector<SimplexId> &saddle2Index) {
  Timer t;

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle connector persistence threshold: "
        << SaddleConnectorsPersistenceThreshold << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  const SimplexId numberOfEdges = inputTriangulation_->getNumberOfEdges();
  const SimplexId numberOfTriangles
    = inputTriangulation_->getNumberOfTriangles();
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

      std::set<SimplexId> saddles2;
      VisitedMask mask{isVisited, visitedIds};
      getAscendingWall(minSaddle1, mask, nullptr, &saddles2);

      // check if at least one connection exists
      auto isFound = saddles2.find(minSaddle2.id_);
      if(isFound == saddles2.end()) {
        ++numberOfIterations;
        continue;
      }

      // check if there is multiple connections
      std::vector<Cell> path;
      const bool isMultiConnected = getDescendingPathThroughWall(
        minSaddle2, minSaddle1, isVisited, &path, true);

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
            inputTriangulation_->getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle1)) {
              continue;
            }

            if(!allowBoundary
               and inputTriangulation_->isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle1[vertexId] == -1) {
              const SimplexId pl_saddle1Id = vertexId;

              SimplexId numberOfRemainingSaddles1 = 0;

              SimplexId savedId = -1;
              const SimplexId edgeNumber
                = inputTriangulation_->getVertexEdgeNumber(pl_saddle1Id);
              for(SimplexId j = 0; j < edgeNumber; ++j) {
                SimplexId edgeId;
                inputTriangulation_->getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId != dmt_saddle1Id and isSaddle1(Cell(1, edgeId))
                   and isRemovableSaddle1[edgeId]) {
                  ++numberOfRemainingSaddles1;
                  savedId = edgeId;
                }
              }

              if(numberOfRemainingSaddles1 == 0) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                pl2dmt_saddle1[vertexId] = dmt_saddle1Id;
                dmt1Saddle2PL_[dmt_saddle1Id] = vertexId;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1 == 1) {
                isRemovableSaddle1[dmt_saddle1Id] = false;
                isRemovableSaddle1[savedId] = false;
                pl2dmt_saddle1[vertexId] = savedId;
                dmt1Saddle2PL_[savedId] = vertexId;
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
            inputTriangulation_->getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId] != static_cast<char>(CriticalType::Saddle2)) {
              continue;
            }

            if(!allowBoundary
               and inputTriangulation_->isVertexOnBoundary(vertexId)) {
              continue;
            }

            if(pl2dmt_saddle2[vertexId] == -1) {
              const SimplexId pl_saddle2Id = vertexId;

              SimplexId numberOfRemainingSaddles2 = 0;

              SimplexId savedId = -1;
              const SimplexId triangleNumber
                = inputTriangulation_->getVertexTriangleNumber(pl_saddle2Id);
              for(SimplexId j = 0; j < triangleNumber; ++j) {
                SimplexId triangleId;
                inputTriangulation_->getVertexTriangle(
                  pl_saddle2Id, j, triangleId);

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
        reverseDescendingPathOnWall(path);
      }
    }

    if(vpath.isValid_) {
      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_) {
        const Cell &minSaddle1 = criticalPoints[vpath.source_].cell_;
        const Cell &minSaddle2 = criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back({minSaddle1, minSaddle2});
      }

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

        std::set<SimplexId> saddles2;
        VisitedMask mask{isVisited, visitedIds};
        getAscendingWall(saddle1, mask, nullptr, &saddles2);

        for(auto &saddle2Id : saddles2) {
          const Cell saddle2(2, saddle2Id);

          const bool isMultiConnected = getDescendingPathThroughWall(
            saddle2, saddle1, isVisited, nullptr, true);

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
            = getPersistence<dataType>(saddle2, saddle1, scalars);

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

        std::set<SimplexId> saddles1;
        VisitedMask mask{isVisited, visitedIds};
        getDescendingWall(saddle2, mask, nullptr, &saddles1);

        for(auto &saddle1Id : saddles1) {
          const Cell saddle1(1, saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected = getAscendingPathThroughWall(
            saddle1, saddle2, isVisited, &path, true);

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
            = getPersistence<dataType>(saddle2, saddle1, scalars);

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

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t"
        << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections2(
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::vector<char> &isPL,
  const int iterationThreshold,
  const bool allowBoundary,
  const bool allowBruteForce,
  const bool returnSaddleConnectors) {
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<SimplexId> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(
    criticalPoints, allowBoundary, isRemovableSaddle1, pl2dmt_saddle1);

  std::vector<char> isRemovableSaddle2;
  std::vector<SimplexId> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(
    criticalPoints, allowBoundary, isRemovableSaddle2, pl2dmt_saddle2);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<SimplexId> saddle1Index;
  std::vector<SimplexId> saddle2Index;
  initializeSaddleSaddleConnections2<dataType>(
    isRemovableSaddle1, isRemovableSaddle2, allowBruteForce, vpaths,
    dmt_criticalPoints, saddle1Index, saddle2Index);

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
    saddle1Index, saddle2Index);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle-Saddle pairs simplified in "
        << t.getElapsedTime() << " s, " << threadNumber_ << " thread(s)."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType, typename idType>
int DiscreteGradient::filterSaddleConnectors(const bool allowBoundary) {
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

  const auto *const offsets = static_cast<const idType *>(inputOffsets_);
  const auto *const scalars = static_cast<const dataType *>(inputScalarField_);

  ftm::FTMTree contourTree;
  contourTree.setDebugLevel(debugLevel_);
  contourTree.setupTriangulation(inputTriangulation_, false);
  contourTree.setVertexScalars(scalars);
  contourTree.setTreeType(ftm::TreeType::Contour);
  contourTree.setVertexSoSoffsets(offsets);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.setSegmentation(false);
  contourTree.build<dataType, SimplexId>();
  ftm::FTMTree_MT *tree = contourTree.getTree(ftm::TreeType::Contour);

  const SimplexId numberOfNodes = tree->getNumberOfNodes();
  for(SimplexId nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
    const ftm::Node *node = tree->getNode(nodeId);
    const SimplexId vertexId = node->getVertexId();

    cpset.push_back(std::make_pair(vertexId, getNodeType(node)));
  }

  std::vector<char> isPL;
  getCriticalPointMap(cpset, isPL);

  simplifySaddleSaddleConnections1<dataType>(cpset, isPL, IterationThreshold,
                                             allowBoundary, allowBruteForce,
                                             returnSaddleConnectors);
  simplifySaddleSaddleConnections2<dataType>(cpset, isPL, IterationThreshold,
                                             allowBoundary, allowBruteForce,
                                             returnSaddleConnectors);

  return 0;
}

template <typename dataType, typename idType>
int DiscreteGradient::reverseGradient(bool detectCriticalPoints) {

  std::vector<std::pair<SimplexId, char>> criticalPoints{};

  if(detectCriticalPoints) {

    // get critical points as cells
    std::vector<Cell> criticalCells{};
    getCriticalPoints(criticalCells);

    criticalPoints.resize(criticalCells.size());

    // iterate over cells to get points (max vertex) and type
    std::transform(
      criticalCells.begin(), criticalCells.end(), criticalPoints.begin(),
      [&](const Cell &c) {
        const auto cellType = criticalTypeFromCellDimension(c.dim_);
        const auto vertexId = getCellGreaterVertex(c);
        return std::make_pair(vertexId, static_cast<char>(cellType));
      });

    // print number of critical cells
    {
      // foreach dimension
      const int numberOfDimensions = getNumberOfDimensions();
      std::vector<SimplexId> nDMTCriticalPoints(numberOfDimensions, 0);
      for(const auto &c : criticalCells) {
        ++nDMTCriticalPoints[c.dim_];
      }

      std::vector<SimplexId> nPLInteriorCriticalPoints(numberOfDimensions, 0);
      for(const auto &cp : criticalPoints) {
        if(!inputTriangulation_->isVertexOnBoundary(cp.first)) {
          ++nPLInteriorCriticalPoints[cp.second];
        }
      }

      {
        std::stringstream msg;
        for(int i = 0; i < numberOfDimensions; ++i) {
          msg << "[DiscreteGradient] " << nDMTCriticalPoints[i] << " " << i
              << "-cell(s)";
          msg << " and " << nPLInteriorCriticalPoints[i] << " interior PL."
              << std::endl;
        }

        dMsg(std::cout, msg.str(), infoMsg);
      }
    }
  }

  Timer t;

  const bool allowBoundary = true;
  const bool returnSaddleConnectors = false;
  bool allowBruteForce = false;

  std::vector<char> isPL;
  getCriticalPointMap(criticalPoints, isPL);

  dmt1Saddle2PL_.resize(inputTriangulation_->getNumberOfEdges());
  std::fill(dmt1Saddle2PL_.begin(), dmt1Saddle2PL_.end(), -1);

  if(dimensionality_ == 3) {
    simplifySaddleSaddleConnections1<dataType>(
      criticalPoints, isPL, IterationThreshold, allowBoundary, allowBruteForce,
      returnSaddleConnectors);
    simplifySaddleSaddleConnections2<dataType>(
      criticalPoints, isPL, IterationThreshold, allowBoundary, allowBruteForce,
      returnSaddleConnectors);
  }

  if(dimensionality_ == 3 and ReturnSaddleConnectors) {
    filterSaddleConnectors<dataType, idType>(allowBoundary);
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Gradient reversed in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}
