/// \ingroup base
/// \class ttk::MorseSmaleComplex3D::
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex3D processing package.
///
/// %MorseSmaleComplex3D is a TTK processing package that takes a scalar field
/// on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex3D.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX3D_H
#define _MORSESMALECOMPLEX3D_H

// base code includes
#include <AbstractMorseSmaleComplex.h>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 3D triangulation.
   */
  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex3D();
    ~MorseSmaleComplex3D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename dataType, typename idtype>
    int execute();

    /**
     * Compute the (saddle1, saddle2) pairs not detected by the
     * contour tree.
     */
    template <typename dataType, typename idType>
    int computePersistencePairs(
      std::vector<std::tuple<SimplexId, SimplexId, dataType>>
        &pl_saddleSaddlePairs);

    template <typename dataType>
    int setAugmentedCriticalPoints(const std::vector<dcg::Cell> &criticalPoints,
                                   SimplexId *ascendingManifold,
                                   SimplexId *descendingManifold) const;

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * gradient.
     */
    int getSaddleConnectors(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the maxima.
     */
    int getDescendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles) const;

    /**
     * Compute the geometrical embedding of the descending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename dataType>
    int setDescendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles) const;

    /**
     * Find all tetras in the star of edgeId
     *
     * (primal: star of edgeId -> dual: vertices of polygon)
     */
    int getDualPolygon(const SimplexId edgeId,
                       std::vector<SimplexId> &polygon) const;

    /**
     * Sort the polygon vertices to be clockwise
     */
    int sortDualPolygonVertices(std::vector<SimplexId> &polygon) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the minima.
     */
    int getAscendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles) const;

    /**
     * Compute the geometrical embedding of the ascending
     * 2-separatrices. This function needs the following
     * internal pointers to be set:
     * outputSeparatrices2_numberOfPoints_
     * outputSeparatrices2_points_
     * outputSeparatrices2_numberOfCells_
     * outputSeparatrices2_cells_
     * inputScalarField_
     */
    template <typename dataType>
    int setAscendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<ttk::Separatrix>> &separatrices,
      std::vector<std::vector<std::vector<ttk::dcg::Cell>>>
        &separatricesGeometry) const;
  };
} // namespace ttk

template <typename dataType>
int ttk::MorseSmaleComplex3D::setAscendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to "
                 "numberOfPoints is null."
              << std::endl;
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    std::cerr
      << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null."
      << std::endl;
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to "
                 "numberOfCells is null."
              << std::endl;
    return -1;
  }
  if(outputSeparatrices2_cells_ == nullptr) {
    std::cerr
      << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null."
      << std::endl;
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input "
                 "scalar field is null."
              << std::endl;
    return -1;
  }
#endif

  const auto scalars = static_cast<const dataType *>(inputScalarField_);
  auto separatrixFunctionMaxima = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionMaxima_);
  auto separatrixFunctionMinima = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionMinima_);
  auto separatrixFunctionDiffs = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionDiffs_);

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = (outputSeparatrices2_cells_separatrixIds_ != nullptr
       && !outputSeparatrices2_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
                            outputSeparatrices2_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices2_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices2_numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outputSeparatrices2_cells_->size()};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      ncells += separatricesGeometry[geomId].size();
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  struct PolygonCell {
    std::vector<SimplexId> tetras_{};
    dataType sepFuncMax_{}, sepFuncMin_{};
    SimplexId sourceId_{}, sepId_{};
    char onBoundary_{};
    bool valid_{false};
  };

  // store the polygonal cells tetras SimplexId
  std::vector<PolygonCell> polygonTetras(ncells - noldcells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto &sepSaddles = separatricesSaddles[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle1

    // compute separatrix function diff
    const dataType sepFuncMin = discreteGradient_.scalarMin(src, scalars);
    const auto maxId = *std::max_element(
      sepSaddles.begin(), sepSaddles.end(),
      [=](const SimplexId a, const SimplexId b) {
        return discreteGradient_.scalarMax(Cell{2, a}, scalars)
               < discreteGradient_.scalarMax(Cell{2, b}, scalars);
      });
    const dataType sepFuncMax
      = discreteGradient_.scalarMax(Cell{2, maxId}, scalars);

    // get boundary condition
    const char onBoundary = std::count_if(
      sepSaddles.begin(), sepSaddles.end(), [=](const SimplexId a) {
        return inputTriangulation_->isEdgeOnBoundary(a);
      });

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      // index of current cell in cell data arrays
      const auto k = geomCellsBegId[i] + j - noldcells;
      auto &polyCell = polygonTetras[k];

      // Transform to dual : edge -> polygon
      getDualPolygon(cell.id_, polyCell.tetras_);

      if(polyCell.tetras_.size() > 2) {
        sortDualPolygonVertices(polyCell.tetras_);
        polyCell.sepFuncMax_ = sepFuncMax;
        polyCell.sepFuncMin_ = sepFuncMin;
        polyCell.sourceId_ = src.id_;
        polyCell.sepId_ = sepId;
        polyCell.onBoundary_ = onBoundary;
        polyCell.valid_ = true;
      }
    }
  }

  decltype(polygonTetras) flatTetras{};
  flatTetras.reserve(polygonTetras.size());

  for(auto &&polyCell : polygonTetras) {
    if(polyCell.valid_) {
      flatTetras.emplace_back(std::move(polyCell));
    }
  }

  // count number of valid new cells and new points
  size_t nnewpoints{};
  std::vector<size_t> pointsPerCell(flatTetras.size() + 1);
  for(size_t i = 0; i < flatTetras.size(); ++i) {
    const auto &poly = flatTetras[i];
    nnewpoints += poly.tetras_.size();
    pointsPerCell[i + 1] = nnewpoints;
  }

  // reduce number of points to remove duplicates
  std::vector<SimplexId> cellVertsIds(nnewpoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < flatTetras.size(); ++i) {
    const auto &poly = flatTetras[i];
    for(size_t j = 0; j < poly.tetras_.size(); ++j) {
      cellVertsIds[pointsPerCell[i] + j] = poly.tetras_[j];
    }
  }

  PSORT(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<size_t> vertId2PointsId(inputTriangulation_->getNumberOfCells());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  ncells = noldcells + flatTetras.size();
  const auto nnewcellids = pointsPerCell.back() + flatTetras.size();

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);
  outputSeparatrices2_cells_->resize(firstCellId + nnewcellids);
  auto cells = &outputSeparatrices2_cells_->at(firstCellId);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(ncells);
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(ncells);
  if(separatrixFunctionDiffs != nullptr)
    separatrixFunctionDiffs->resize(ncells);
  if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
    outputSeparatrices2_cells_isOnBoundary_->resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    inputTriangulation_->getTetraIncenter(cellVertsIds[i], &points[3 * i]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < flatTetras.size(); ++i) {
    const auto &poly = flatTetras[i];
    const auto k = pointsPerCell[i] + i;
    cells[k + 0] = poly.tetras_.size();
    for(size_t j = 0; j < poly.tetras_.size(); ++j) {
      cells[k + 1 + j] = vertId2PointsId[poly.tetras_[j]];
    }
    const auto l = i + noldcells;
    if(outputSeparatrices2_cells_sourceIds_ != nullptr)
      (*outputSeparatrices2_cells_sourceIds_)[l] = poly.sourceId_;
    if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
      (*outputSeparatrices2_cells_separatrixIds_)[l] = poly.sepId_;
    if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
      (*outputSeparatrices2_cells_separatrixTypes_)[l] = 1;
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[l] = poly.sepFuncMax_;
    if(separatrixFunctionDiffs != nullptr)
      (*separatrixFunctionMinima)[l] = poly.sepFuncMin_;
    if(separatrixFunctionDiffs != nullptr)
      (*separatrixFunctionDiffs)[l] = poly.sepFuncMax_ - poly.sepFuncMin_;
    if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices2_cells_isOnBoundary_)[l] = poly.onBoundary_;
  }

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}

template <typename dataType>
int ttk::MorseSmaleComplex3D::setDescendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to "
                 "numberOfPoints is null."
              << std::endl;
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    std::cerr
      << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null."
      << std::endl;
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to "
                 "numberOfCells is null."
              << std::endl;
    return -1;
  }
  if(outputSeparatrices2_cells_ == nullptr) {
    std::cerr
      << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null."
      << std::endl;
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input "
                 "scalar field is null."
              << std::endl;
    return -1;
  }
#endif

  const auto scalars = static_cast<const dataType *>(inputScalarField_);
  auto separatrixFunctionMaxima = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionMaxima_);
  auto separatrixFunctionMinima = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionMinima_);
  auto separatrixFunctionDiffs = static_cast<std::vector<dataType> *>(
    outputSeparatrices2_cells_separatrixFunctionDiffs_);

  // max existing separatrix id + 1 or 0 if no previous separatrices
  const SimplexId separatrixId
    = (outputSeparatrices2_cells_separatrixIds_ != nullptr
       && !outputSeparatrices2_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
                            outputSeparatrices2_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices2_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices2_numberOfCells_)};
  // old number of separatrices cells
  const auto noldcells{ncells};
  // index of last vertex of last old cell + 1
  const auto firstCellId{outputSeparatrices2_cells_->size()};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      ncells += separatricesGeometry[geomId].size();
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  // resize arrays
  outputSeparatrices2_cells_->resize(
    firstCellId + 4 * (ncells - noldcells)); // triangles cells
  auto cells = &outputSeparatrices2_cells_->at(firstCellId);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(ncells);
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(ncells);
  if(separatrixFunctionDiffs != nullptr)
    separatrixFunctionDiffs->resize(ncells);
  if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
    outputSeparatrices2_cells_isOnBoundary_->resize(ncells);

  // store the cells/triangles vertices vertexId
  std::vector<SimplexId> cellVertsIds(3 * (ncells - noldcells));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto &sepSaddles = separatricesSaddles[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    const dcg::Cell &src = sep.source_; // saddle2
    const char sepType = 2;

    // compute separatrix function diff
    const dataType sepFuncMax = discreteGradient_.scalarMax(src, scalars);
    const auto minId = *std::min_element(
      sepSaddles.begin(), sepSaddles.end(),
      [=](const SimplexId a, const SimplexId b) {
        return discreteGradient_.scalarMin(Cell{1, a}, scalars)
               < discreteGradient_.scalarMin(Cell{1, b}, scalars);
      });
    const dataType sepFuncMin
      = discreteGradient_.scalarMin(Cell{1, minId}, scalars);
    const dataType sepFuncDiff = sepFuncMax - sepFuncMin;

    // get boundary condition
    const char onBoundary = std::count_if(
      sepSaddles.begin(), sepSaddles.end(), [=](const SimplexId a) {
        return inputTriangulation_->isEdgeOnBoundary(a);
      });

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];

      // first store the SimplexId of the cell/triangle vertices
      SimplexId v0{}, v1{}, v2{};
      inputTriangulation_->getTriangleVertex(cell.id_, 0, v0);
      inputTriangulation_->getTriangleVertex(cell.id_, 1, v1);
      inputTriangulation_->getTriangleVertex(cell.id_, 2, v2);

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j;
      // index of current cell among all new cells
      const auto m = l - noldcells;

      cells[4 * m + 0] = 3;
      cells[4 * m + 1] = v0;
      cells[4 * m + 2] = v1;
      cells[4 * m + 3] = v2;
      cellVertsIds[3 * m + 0] = v0;
      cellVertsIds[3 * m + 1] = v1;
      cellVertsIds[3 * m + 2] = v2;

      if(outputSeparatrices2_cells_sourceIds_ != nullptr)
        (*outputSeparatrices2_cells_sourceIds_)[l] = src.id_;
      if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
        (*outputSeparatrices2_cells_separatrixIds_)[l] = sepId;
      if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices2_cells_separatrixTypes_)[l] = sepType;
      if(separatrixFunctionMaxima != nullptr)
        (*separatrixFunctionMaxima)[l] = sepFuncMax;
      if(separatrixFunctionDiffs != nullptr)
        (*separatrixFunctionMinima)[l] = sepFuncMin;
      if(separatrixFunctionDiffs != nullptr)
        (*separatrixFunctionDiffs)[l] = sepFuncDiff;
      if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices2_cells_isOnBoundary_)[l] = onBoundary;
    }
  }

  // reduce the cell vertices ids
  // (cells are triangles sharing two vertices)
  PSORT(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<size_t> vertId2PointsId(
    inputTriangulation_->getNumberOfVertices());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    inputTriangulation_->getVertexPoint(
      cellVertsIds[i], points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < ncells - noldcells; ++i) {
    cells[4 * i + 1] = vertId2PointsId[cells[4 * i + 1]];
    cells[4 * i + 2] = vertId2PointsId[cells[4 * i + 2]];
    cells[4 * i + 3] = vertId2PointsId[cells[4 * i + 3]];
  }

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}

template <typename dataType, typename idType>
int ttk::MorseSmaleComplex3D::execute() {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField_) {
    std::cerr
      << "[MorseSmaleComplex3D] Error: input scalar field pointer is null."
      << std::endl;
    return -1;
  }

  if(!inputOffsets_) {
    std::cerr
      << "[MorseSmaleComplex3D] Error: input offset field pointer is null."
      << std::endl;
    return -1;
  }
#endif
  Timer t;

  // nullptr_t is implicitly convertible and comparable to any pointer type
  // or pointer-to-member type.
  SimplexId *ascendingManifold
    = static_cast<SimplexId *>(outputAscendingManifold_);
  SimplexId *descendingManifold
    = static_cast<SimplexId *>(outputDescendingManifold_);
  SimplexId *morseSmaleManifold
    = static_cast<SimplexId *>(outputMorseSmaleManifold_);

  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.setDebugLevel(debugLevel_);
  {
    Timer tmp;
    discreteGradient_.buildGradient<dataType, idType>();

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Discrete gradient overall computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ReturnSaddleConnectors) {
    discreteGradient_.reverseGradient<dataType, idType>();
  }

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  std::vector<std::vector<Separatrix>> separatrices1{};
  std::vector<std::vector<std::vector<dcg::Cell>>> separatricesGeometry1;

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getDescendingSeparatrices1(
      criticalPoints, separatrices1.back(), separatricesGeometry1.back());

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Descending 1-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getAscendingSeparatrices1(
      criticalPoints, separatrices1.back(), separatricesGeometry1.back());

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Ascending 1-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  // saddle-connectors
  if(ComputeSaddleConnectors) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getSaddleConnectors(
      criticalPoints, separatrices1.back(), separatricesGeometry1.back());

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Saddle connectors computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeDescendingSeparatrices1 || ComputeAscendingSeparatrices1
     || ComputeSaddleConnectors) {
    Timer tmp{};

    flattenSeparatricesVectors(separatrices1, separatricesGeometry1);
    setSeparatrices1<dataType>(separatrices1[0], separatricesGeometry1[0]);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] 1-separatrices set in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getDescendingSeparatrices2(
      criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    setDescendingSeparatrices2<dataType>(
      separatrices, separatricesGeometry, separatricesSaddles);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Descending 2-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getAscendingSeparatrices2(
      criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    setAscendingSeparatrices2<dataType>(
      separatrices, separatricesGeometry, separatricesSaddles);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Ascending 2-separatrices computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  std::vector<SimplexId> maxSeeds;
  {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ascendingManifold)
      setAscendingSegmentation(
        criticalPoints, maxSeeds, ascendingManifold, numberOfMaxima);

    if(descendingManifold)
      setDescendingSegmentation(
        criticalPoints, descendingManifold, numberOfMinima);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold,
                           descendingManifold, morseSmaleManifold);

    if(ascendingManifold or descendingManifold) {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Segmentation computed in "
          << tmp.getElapsedTime() << " s." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(outputCriticalPoints_numberOfPoints_ and outputSeparatrices1_points_) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints<dataType>(
      criticalPoints, nCriticalPointsByDim);

    if(ascendingManifold and descendingManifold) {
      discreteGradient_.setManifoldSize(criticalPoints, nCriticalPointsByDim,
                                        maxSeeds, ascendingManifold,
                                        descendingManifold);
    }
  }

  {
    const SimplexId numberOfVertices
      = inputTriangulation_->getNumberOfVertices();
    std::stringstream msg;
    msg << "[MorseSmaleComplex3D] Data-set (" << numberOfVertices
        << " points) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType, typename idType>
int ttk::MorseSmaleComplex3D::computePersistencePairs(
  std::vector<std::tuple<SimplexId, SimplexId, dataType>>
    &pl_saddleSaddlePairs) {

  const dataType *scalars = static_cast<const dataType *>(inputScalarField_);

  std::vector<std::array<dcg::Cell, 2>> dmt_pairs;
  {
    // simplify to be PL-conformant
    discreteGradient_.setDebugLevel(debugLevel_);
    discreteGradient_.setThreadNumber(threadNumber_);
    discreteGradient_.setCollectPersistencePairs(false);
    discreteGradient_.buildGradient<dataType, idType>();
    discreteGradient_.reverseGradient<dataType, idType>();

    // collect saddle-saddle connections
    discreteGradient_.setCollectPersistencePairs(true);
    discreteGradient_.setOutputPersistencePairs(&dmt_pairs);
    discreteGradient_.reverseGradient<dataType, idType>(false);
  }

  // transform DMT pairs into PL pairs
  for(const auto &pair : dmt_pairs) {
    const SimplexId v0 = discreteGradient_.getCellGreaterVertex(pair[0]);
    const SimplexId v1 = discreteGradient_.getCellGreaterVertex(pair[1]);
    const dataType persistence = scalars[v1] - scalars[v0];

    if(v0 != -1 and v1 != -1 and persistence >= 0) {
      if(!inputTriangulation_->isVertexOnBoundary(v0)
         or !inputTriangulation_->isVertexOnBoundary(v1)) {
        pl_saddleSaddlePairs.emplace_back(v0, v1, persistence);
      }
    }
  }
  return 0;
}

#endif // MORSESMALECOMPLEX3D_H
