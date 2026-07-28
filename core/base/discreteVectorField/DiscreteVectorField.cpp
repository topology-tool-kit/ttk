#include <DiscreteVectorField.h>

using namespace std;
using namespace ttk;
using namespace dcvf;

int DiscreteVectorField::getDimensionality() const {
  return dimensionality_;
}

int DiscreteVectorField::getNumberOfDimensions() const {
  return dimensionality_ + 1;
}

void DiscreteVectorField::initMemory(
  const AbstractTriangulation &triangulation) {

  Timer tm{};
  const int numberOfDimensions = this->getNumberOfDimensions();

  // init number of cells by dimension
  std::vector<SimplexId> numberOfCells(numberOfDimensions);
  for(int i = 0; i < numberOfDimensions; ++i) {
    numberOfCells[i] = this->getNumberOfCells(i, triangulation);
  }

  // clear & init discrete vectors memory
  for(int i = 0; i < dimensionality_; ++i) {
    (*vectors_)[2 * i].clear();
    (*vectors_)[2 * i].resize(numberOfCells[i], -1);
    (*vectors_)[2 * i + 1].clear();
    (*vectors_)[2 * i + 1].resize(numberOfCells[i + 1], -1);
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
  this->printMsg("Initialized discrete vectors memory", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);
}

std::pair<size_t, SimplexId>
  DiscreteVectorField::numUnpairedFaces(const CellOutExt &c,
                                        const outwardStarType &ls) const {
  // c.dim_ cannot be <= 1
  if(c.dim_ == 2) {
    return numUnpairedFacesTriangle(c, ls);
  } else if(c.dim_ == 3) {
    return numUnpairedFacesTetra(c, ls);
  }

  return {0, -1};
}

std::pair<size_t, SimplexId> DiscreteVectorField::numUnpairedFacesTriangle(
  const CellOutExt &c, const outwardStarType &ls) const {
  // number of unpaired faces
  std::pair<size_t, SimplexId> res{0, -1};

  // loop over edge faces of triangle
  // (2 edges per triangle in outward star)
  for(size_t i = 0; i < 2; ++i) {
    if(!ls[1][c.faces_[i]].paired_) {
      res.first++;
      res.second = c.faces_[i];
    }
  }

  return res;
}

std::pair<size_t, SimplexId>
  DiscreteVectorField::numUnpairedFacesTetra(const CellOutExt &c,
                                             const outwardStarType &ls) const {
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

bool DiscreteVectorField::isCellCritical(const int cellDim,
                                         const SimplexId cellId) const {

  if(cellDim > this->dimensionality_) {
    return false;
  }
#ifndef TTK_ENABLE_KAMIKAZE
  if(cellId < 0) {
    this->printErr("Invalid cell ID given to isCellCritical");
    return false;
  }
#endif

  if(cellDim == 0) {
    return ((*vectors_)[0][cellId] == NULL_CONNECTION);
  }

  if(cellDim == 1) {
    return (
      (*vectors_)[1][cellId] == NULL_CONNECTION
      && (dimensionality_ == 1 || (*vectors_)[2][cellId] == NULL_CONNECTION));
  }

  if(cellDim == 2) {
    return (
      (*vectors_)[3][cellId] == NULL_CONNECTION
      && (dimensionality_ == 2 || (*vectors_)[4][cellId] == NULL_CONNECTION));
  }

  if(cellDim == 3) {
    return ((*vectors_)[5][cellId] == NULL_CONNECTION);
  }

  return false;
}

bool DiscreteVectorField::isCellCritical(const Cell &cell) const {
  return isCellCritical(cell.dim_, cell.id_);
}

int DiscreteVectorField::setManifoldSize(
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
    const SimplexId nMin = static_cast<SimplexId>(criticalCellsByDim[0].size());
    for(SimplexId i = 0; i < numberOfVertices_; ++i) {
      if(descendingManifold[i] != -1 && descendingManifold[i] < nMin) {
        manifoldSize[descendingManifold[i]]++;
      }
    }
  }

  if(!criticalCellsByDim[dim].empty()) {
    // index of first maximum in critical points array
    const auto nFirstMaximum{nCritPoints - criticalCellsByDim[dim].size()};
    const SimplexId nMax
      = static_cast<SimplexId>(criticalCellsByDim[dim].size());
    // ascending manifold cells size
    for(SimplexId i = 0; i < numberOfVertices_; ++i) {
      if(ascendingManifold[i] != -1 && ascendingManifold[i] < nMax) {
        manifoldSize[ascendingManifold[i] + nFirstMaximum]++;
      }
    }
  }

  return 0;
}

#ifdef TTK_ENABLE_MPI
void DiscreteVectorField::setCellToGhost(const int cellDim,
                                         const SimplexId cellId) {
  if(cellDim == 0) {
    (*vectors_)[0][cellId] = GHOST_CONNECTION;
  }

  if(cellDim == 1) {
    (*vectors_)[1][cellId] = GHOST_CONNECTION;
    (*vectors_)[2][cellId] = GHOST_CONNECTION;
  }

  if(cellDim == 2) {
    (*vectors_)[3][cellId] = GHOST_CONNECTION;
    (*vectors_)[4][cellId] = GHOST_CONNECTION;
  }

  if(cellDim == 3) {
    (*vectors_)[5][cellId] = GHOST_CONNECTION;
  }
}
#endif
