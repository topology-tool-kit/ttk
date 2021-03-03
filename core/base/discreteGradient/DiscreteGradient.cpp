#include <DiscreteGradient.h>

using namespace std;
using namespace ttk;
using namespace dcg;

int DiscreteGradient::getDimensionality() const {
  return dimensionality_;
}

int DiscreteGradient::getNumberOfDimensions() const {
  return dimensionality_ + 1;
}

void DiscreteGradient::initMemory(const AbstractTriangulation &triangulation) {

  Timer tm{};
  const int numberOfDimensions = this->getNumberOfDimensions();

  // init number of cells by dimension
  std::vector<SimplexId> numberOfCells(numberOfDimensions);
  for(int i = 0; i < numberOfDimensions; ++i) {
    numberOfCells[i] = this->getNumberOfCells(i, triangulation);
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

  std::vector<std::vector<std::string>> rows{
    {"#Vertices", std::to_string(numberOfCells[0])},
    {"#Edges", std::to_string(numberOfCells[1])},
    {"#Triangles", std::to_string(numberOfCells[2])}};

  if(dimensionality_ == 3) {
    rows.emplace_back(
      std::vector<std::string>{"#Tetras", std::to_string(numberOfCells[3])});
  }

  this->printMsg(rows);
  this->printMsg("Initialized discrete gradient memory", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFaces(const CellExt &c,
                                     const lowerStarType &ls) const {
  // c.dim_ cannot be <= 1
  if(c.dim_ == 2) {
    return numUnpairedFacesTriangle(c, ls);
  } else if(c.dim_ == 3) {
    return numUnpairedFacesTetra(c, ls);
  }

  return {0, -1};
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFacesTriangle(const CellExt &c,
                                             const lowerStarType &ls) const {
  // number of unpaired faces
  std::pair<size_t, SimplexId> res{0, -1};

  // loop over edge faces of triangle
  // (2 edges per triangle in lower star)
  for(size_t i = 0; i < 2; ++i) {
    if(!ls[1][c.faces_[i]].paired_) {
      res.first++;
      res.second = c.faces_[i];
    }
  }

  return res;
}

std::pair<size_t, SimplexId>
  DiscreteGradient::numUnpairedFacesTetra(const CellExt &c,
                                          const lowerStarType &ls) const {
  // number of unpaired faces
  std::pair<size_t, SimplexId> res{0, -1};

  // loop over triangle faces of tetra
  for(const auto f : c.faces_) {
    if(!ls[2][f].paired_) {
      res.first++;
      res.second = f;
    }
  }

  return res;
}

CriticalType
  DiscreteGradient::criticalTypeFromCellDimension(const int dim) const {
  if(dim == 0) {
    return CriticalType::Local_minimum;
  } else if(dim == 1) {
    return CriticalType::Saddle1;
  } else if(dim == 2 && dimensionality_ == 3) {
    return CriticalType::Saddle2;
  } else if(dim == dimensionality_) {
    return CriticalType::Local_maximum;
  } else {
    return CriticalType::Regular;
  }
}

bool DiscreteGradient::isMinimum(const Cell &cell) const {
  if(cell.dim_ == 0) {
    return (gradient_[0][0][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle1(const Cell &cell) const {
  if(cell.dim_ == 1) {
    return (gradient_[0][1][cell.id_] == -1
            and gradient_[1][1][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle2(const Cell &cell) const {
  if(dimensionality_ == 3 and cell.dim_ == 2) {
    return (gradient_[1][2][cell.id_] == -1
            and gradient_[2][2][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isMaximum(const Cell &cell) const {
  if(dimensionality_ == 2 and cell.dim_ == 2) {
    return (gradient_[1][2][cell.id_] == -1);
  }

  if(dimensionality_ == 3 and cell.dim_ == 3) {
    return (gradient_[2][3][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isCellCritical(const int cellDim,
                                      const SimplexId cellId) const {
  if(dimensionality_ == 2) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1);
        break;
    }
  } else if(dimensionality_ == 3) {
    switch(cellDim) {
      case 0:
        return (gradient_[0][0][cellId] == -1);
        break;

      case 1:
        return (gradient_[0][1][cellId] == -1
                and gradient_[1][1][cellId] == -1);
        break;

      case 2:
        return (gradient_[1][2][cellId] == -1
                and gradient_[2][2][cellId] == -1);
        break;

      case 3:
        return (gradient_[2][3][cellId] == -1);
        break;
    }
  }
  return false;
}

bool DiscreteGradient::isCellCritical(const Cell &cell) const {
  return isCellCritical(cell.dim_, cell.id_);
}

int DiscreteGradient::getCriticalPointMap(
  const vector<pair<SimplexId, char>> &criticalPoints, vector<char> &isPL) {
  isPL.resize(numberOfVertices_);
  std::fill(isPL.begin(), isPL.end(), 0);
  for(pair<SimplexId, char> criticalPoint : criticalPoints) {
    const SimplexId criticalPointId = criticalPoint.first;
    const char criticalPointType = criticalPoint.second;

    isPL[criticalPointId] = criticalPointType;
  }

  return 0;
}

int DiscreteGradient::setManifoldSize(
  const std::vector<Cell> &criticalPoints,
  const std::vector<size_t> &nCriticalPointsByDim,
  const std::vector<SimplexId> &maxSeeds,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold) {

  const auto nCritPoints = criticalPoints.size();
  const auto nDimensions = getNumberOfDimensions();

  outputCriticalPoints_points_manifoldSize_.resize(nCritPoints, 0);

  // pre-compute size of descending manifold cells
  std::map<SimplexId, size_t> descendingCellsSize{};
  for(SimplexId i = 0; i < numberOfVertices_; ++i) {
    descendingCellsSize[descendingManifold[i]]++;
  }

  // minima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCriticalPointsByDim[0]; ++i) {
    const Cell &cell = criticalPoints[i];
    const SimplexId seedId = descendingManifold[cell.id_];
    const SimplexId manifoldSize = descendingCellsSize[seedId];
    outputCriticalPoints_points_manifoldSize_[i] = manifoldSize;
  }

  // index of first maximum in critical points array
  size_t nFirstMaximum{};
  for(int i = 0; i < nDimensions - 1; ++i) {
    nFirstMaximum += nCriticalPointsByDim[i];
  }

  // pre-compute size of ascending manifold cells
  std::map<SimplexId, size_t> ascendingCellsSize{};
  for(SimplexId i = 0; i < numberOfVertices_; ++i) {
    ascendingCellsSize[ascendingManifold[i]]++;
  }

  // pre-compute maximum SimplexId -> index in maxSeeds
  std::map<SimplexId, size_t> seedsPos{};
  for(size_t i = 0; i < maxSeeds.size(); ++i) {
    seedsPos[maxSeeds[i]] = i;
  }

  // maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = nFirstMaximum; i < nCritPoints; ++i) {
    const Cell &cell = criticalPoints[i];
    if(seedsPos.find(cell.id_) != seedsPos.end()) {
      const auto seedId = seedsPos[cell.id_];
      const SimplexId manifoldSize = ascendingCellsSize[seedId];
      outputCriticalPoints_points_manifoldSize_[i] = manifoldSize;
    }
  }

  return 0;
}
