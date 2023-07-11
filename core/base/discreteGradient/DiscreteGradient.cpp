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

  // clear & init gradient memory
  for(int i = 0; i < dimensionality_; ++i) {
    (*gradient_)[2 * i].clear();
    (*gradient_)[2 * i].resize(numberOfCells[i], -1);
    (*gradient_)[2 * i + 1].clear();
    (*gradient_)[2 * i + 1].resize(numberOfCells[i + 1], -1);
  }

  std::vector<std::vector<std::string>> rows{
    {"#Vertices", std::to_string(numberOfCells[0])},
    {"#Edges", std::to_string(numberOfCells[1])},
  };

  if(dimensionality_ >= 2) {
    rows.emplace_back(
      std::vector<std::string>{"#Triangles", std::to_string(numberOfCells[2])});
  }

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
    return ((*gradient_)[0][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle1(const Cell &cell) const {
  if(cell.dim_ == 1) {
    return ((*gradient_)[1][cell.id_] == -1
            and (*gradient_)[2][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isSaddle2(const Cell &cell) const {
  if(dimensionality_ == 3 and cell.dim_ == 2) {
    return ((*gradient_)[3][cell.id_] == -1
            and (*gradient_)[4][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isMaximum(const Cell &cell) const {
  if(dimensionality_ == 1 and cell.dim_ == 1) {
    return ((*gradient_)[1][cell.id_] == -1);
  }

  if(dimensionality_ == 2 and cell.dim_ == 2) {
    return ((*gradient_)[3][cell.id_] == -1);
  }

  if(dimensionality_ == 3 and cell.dim_ == 3) {
    return ((*gradient_)[5][cell.id_] == -1);
  }

  return false;
}

bool DiscreteGradient::isCellCritical(const int cellDim,
                                      const SimplexId cellId) const {

  if(cellDim > this->dimensionality_) {
    return false;
  }

  if(cellDim == 0) {
    return ((*gradient_)[0][cellId] == NULL_GRADIENT);
  }

  if(cellDim == 1) {
    return (
      (*gradient_)[1][cellId] == NULL_GRADIENT
      && (dimensionality_ == 1 || (*gradient_)[2][cellId] == NULL_GRADIENT));
  }

  if(cellDim == 2) {
    return (
      (*gradient_)[3][cellId] == NULL_GRADIENT
      && (dimensionality_ == 2 || (*gradient_)[4][cellId] == NULL_GRADIENT));
  }

  if(cellDim == 3) {
    return ((*gradient_)[5][cellId] == NULL_GRADIENT);
  }

  return false;
}

bool DiscreteGradient::isCellCritical(const Cell &cell) const {
  return isCellCritical(cell.dim_, cell.id_);
}

int DiscreteGradient::setManifoldSize(
  const std::array<std::vector<SimplexId>, 4> &criticalCellsByDim,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  std::vector<SimplexId> &manifoldSize) const {

  const auto nCritPoints{
    criticalCellsByDim[0].size() + criticalCellsByDim[1].size()
    + criticalCellsByDim[2].size() + criticalCellsByDim[3].size()};

  const auto dim{this->dimensionality_};

  if(nCritPoints == 0
     || (criticalCellsByDim[0].empty() && criticalCellsByDim[dim].empty())) {
    // no critical points || no extrema
    return 0;
  }

  manifoldSize.resize(nCritPoints, 0);

  // descending manifold cells size
  if(!criticalCellsByDim[0].empty()) {
    for(SimplexId i = 0; i < numberOfVertices_; ++i) {
      if(descendingManifold[i] != -1) {
        manifoldSize[descendingManifold[i]]++;
      }
    }
  }

  if(!criticalCellsByDim[dim].empty()) {
    // index of first maximum in critical points array
    const auto nFirstMaximum{nCritPoints - criticalCellsByDim[dim].size()};

    // ascending manifold cells size
    for(SimplexId i = 0; i < numberOfVertices_; ++i) {
      if(ascendingManifold[i] != -1) {
        manifoldSize[ascendingManifold[i] + nFirstMaximum]++;
      }
    }
  }

  return 0;
}

#ifdef TTK_ENABLE_MPI
void DiscreteGradient::setCellToGhost(const int cellDim,
                                      const SimplexId cellId) {
  if(cellDim == 0) {
    (*gradient_)[0][cellId] = GHOST_GRADIENT;
  }

  if(cellDim == 1) {
    (*gradient_)[1][cellId] = GHOST_GRADIENT;
    (*gradient_)[2][cellId] = GHOST_GRADIENT;
  }

  if(cellDim == 2) {
    (*gradient_)[3][cellId] = GHOST_GRADIENT;
    (*gradient_)[4][cellId] = GHOST_GRADIENT;
  }

  if(cellDim == 3) {
    (*gradient_)[5][cellId] = GHOST_GRADIENT;
  }
}
#endif
