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

#pragma once

// base code includes
#include <AbstractMorseSmaleComplex.h>

#include <set>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 3D triangulation.
   */
  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex3D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename dataType, typename triangulationType>
    int execute(const triangulationType &triangulation);

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the saddle-connectors by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getSaddleConnectors(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the maxima.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

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
    template <typename triangulationType>
    int setDescendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * Find all tetras in the star of edgeId
     *
     * (primal: star of edgeId -> dual: vertices of polygon)
     */
    template <typename triangulationType>
    int getDualPolygon(const SimplexId edgeId,
                       SimplexId *const polygon,
                       const size_t polSize,
                       const triangulationType &triangulation) const;

    /**
     * Sort the polygon vertices to be clockwise
     */
    template <typename triangulationType>
    int sortDualPolygonVertices(SimplexId *const polygon,
                                const size_t polSize,
                                const triangulationType &triangulation) const;

    /**
     * Compute the 2-separatrices by reading into the discrete
     * gradient from the minima.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices2(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

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
    template <typename triangulationType>
    int setAscendingSeparatrices2(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const std::vector<std::set<SimplexId>> &separatricesSaddles,
      const triangulationType &triangulation) const;

    /**
     * @brief Flatten the vectors of vectors into their first component
     */
    void flattenSeparatricesVectors(
      std::vector<std::vector<ttk::Separatrix>> &separatrices,
      std::vector<std::vector<std::vector<ttk::dcg::Cell>>>
        &separatricesGeometry) const;
  };
} // namespace ttk

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::setAscendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    this->printErr("2-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices2_cells_connectivity_ == nullptr
     || outputSeparatrices2_cells_offsets_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  const auto offsets = inputOffsets_;
  auto separatrixFunctionMaxima = outputS2_cells_separatrixFunctionMaximaId_;
  auto separatrixFunctionMinima = outputS2_cells_separatrixFunctionMinimaId_;

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
  const auto firstCellId{outputSeparatrices2_cells_connectivity_->size()};
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

  // store the separatrices info (one per separatrix)
  std::vector<SimplexId> sepSourceIds(validGeomIds.size());
  std::vector<SimplexId> sepIds(validGeomIds.size());
  std::vector<char> sepOnBoundary(validGeomIds.size());
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(separatrixId + validGeomIds.size());
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(separatrixId + validGeomIds.size());
  // store the polygonal cells tetras SimplexId
  std::vector<SimplexId> polygonNTetras(ncells - noldcells);
  std::vector<SimplexId> polygonEdgeIds(ncells - noldcells);
  std::vector<SimplexId> polygonSepInfosIds(ncells - noldcells);

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
    const auto sepFuncMin
      = discreteGradient_.getCellLowerVertex(src, triangulation);
    const auto maxId = *std::max_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
        return offsets[discreteGradient_.getCellGreaterVertex(
                 Cell{2, a}, triangulation)]
               < offsets[discreteGradient_.getCellGreaterVertex(
                 Cell{2, b}, triangulation)];
      });
    const auto sepFuncMax
      = discreteGradient_.getCellGreaterVertex(Cell{2, maxId}, triangulation);

    // get boundary condition
    const char onBoundary
      = std::count_if(sepSaddles.begin(), sepSaddles.end(),
                      [&triangulation](const SimplexId a) {
                        return triangulation.isEdgeOnBoundary(a);
                      });

    sepIds[i] = sepId;
    sepSourceIds[i] = src.id_;
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[sepId] = sepFuncMax;
    if(separatrixFunctionMinima != nullptr)
      (*separatrixFunctionMinima)[sepId] = sepFuncMin;
    sepOnBoundary[i] = onBoundary;

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      // index of current cell in cell data arrays
      const auto k = geomCellsBegId[i] + j - noldcells;

      polygonNTetras[k] = triangulation.getEdgeStarNumber(cell.id_);

      if(polygonNTetras[k] > 2) {
        polygonEdgeIds[k] = cell.id_;
        polygonSepInfosIds[k] = i;
      }
    }
  }

  // indices of valid polygon tetras
  std::vector<SimplexId> validTetraIds{};
  validTetraIds.reserve(polygonNTetras.size());

  for(size_t i = 0; i < polygonNTetras.size(); ++i) {
    if(polygonNTetras[i] > 2) {
      validTetraIds.emplace_back(i);
    }
  }

  // count number of valid new cells and new points
  size_t nnewpoints{};
  std::vector<SimplexId> pointsPerCell(validTetraIds.size() + 1);
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    nnewpoints += polygonNTetras[validTetraIds[i]];
    pointsPerCell[i + 1] = nnewpoints;
  }

  // resize connectivity array
  outputSeparatrices2_cells_connectivity_->resize(firstCellId + nnewpoints);
  auto cellsConn = &outputSeparatrices2_cells_connectivity_->at(firstCellId);
  // copy of cell connectivity array (for removing duplicates vertices)
  std::vector<SimplexId> cellVertsIds(nnewpoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    const auto k = validTetraIds[i];

    // get tetras in edge star
    getDualPolygon(polygonEdgeIds[k], &cellVertsIds[pointsPerCell[i]],
                   polygonNTetras[k], triangulation);
    // sort tetras (in-place)
    sortDualPolygonVertices(
      &cellVertsIds[pointsPerCell[i]], polygonNTetras[k], triangulation);

    for(SimplexId j = 0; j < polygonNTetras[k]; ++j) {
      cellsConn[pointsPerCell[i] + j] = cellVertsIds[pointsPerCell[i] + j];
    }
  }

  PSORT(this->threadNumber_)(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<SimplexId> vertId2PointsId(triangulation.getNumberOfCells());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  ncells = noldcells + validTetraIds.size();

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);
  outputSeparatrices2_cells_offsets_->resize(ncells + 1);
  outputSeparatrices2_cells_offsets_->at(0) = 0;
  auto cellsOff = &outputSeparatrices2_cells_offsets_->at(noldcells);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
    outputSeparatrices2_cells_isOnBoundary_->resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    triangulation.getTetraIncenter(cellVertsIds[i], &points[3 * i]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    const auto m = validTetraIds[i];
    const auto k = pointsPerCell[i];
    for(SimplexId j = 0; j < polygonNTetras[m]; ++j) {
      cellsConn[k + j] = vertId2PointsId[cellsConn[k + j]];
    }
    const auto l = i + noldcells;
    const auto n = polygonSepInfosIds[m];

    if(outputSeparatrices2_cells_sourceIds_ != nullptr)
      (*outputSeparatrices2_cells_sourceIds_)[l] = sepSourceIds[n];
    if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
      (*outputSeparatrices2_cells_separatrixIds_)[l] = sepIds[n];
    if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
      (*outputSeparatrices2_cells_separatrixTypes_)[l] = 1;
    if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices2_cells_isOnBoundary_)[l] = sepOnBoundary[n];
  }

  for(size_t i = 0; i < validTetraIds.size(); ++i) {
    // fill offsets sequentially (due to iteration dependencies)
    cellsOff[i + 1] = cellsOff[i] + polygonNTetras[validTetraIds[i]];
  }

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::setDescendingSeparatrices2(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    this->printErr("2-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices2_cells_connectivity_ == nullptr
     || outputSeparatrices2_cells_offsets_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  const auto offsets = inputOffsets_;
  auto separatrixFunctionMaxima = outputS2_cells_separatrixFunctionMaximaId_;
  auto separatrixFunctionMinima = outputS2_cells_separatrixFunctionMinimaId_;

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
  const auto firstCellId{outputSeparatrices2_cells_connectivity_->size()};

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
  outputSeparatrices2_cells_offsets_->resize(ncells + 1);
  outputSeparatrices2_cells_offsets_->at(0) = 0;
  outputSeparatrices2_cells_connectivity_->resize(
    firstCellId + 3 * (ncells - noldcells)); // triangles cells
  auto cellsOff = &outputSeparatrices2_cells_offsets_->at(noldcells);
  auto cellsConn = &outputSeparatrices2_cells_connectivity_->at(firstCellId);
  if(outputSeparatrices2_cells_sourceIds_ != nullptr)
    outputSeparatrices2_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
    outputSeparatrices2_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
    outputSeparatrices2_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(separatrixId + validGeomIds.size());
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(separatrixId + validGeomIds.size());
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
    const auto sepFuncMax
      = discreteGradient_.getCellGreaterVertex(src, triangulation);
    const auto minId = *std::min_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, offsets, this](const SimplexId a, const SimplexId b) {
        return offsets[discreteGradient_.getCellLowerVertex(
                 Cell{1, a}, triangulation)]
               < offsets[discreteGradient_.getCellLowerVertex(
                 Cell{1, b}, triangulation)];
      });
    const auto sepFuncMin
      = discreteGradient_.getCellLowerVertex(Cell{1, minId}, triangulation);
    if(separatrixFunctionMaxima != nullptr)
      (*separatrixFunctionMaxima)[sepId] = sepFuncMax;
    if(separatrixFunctionMinima != nullptr)
      (*separatrixFunctionMinima)[sepId] = sepFuncMin;

    // get boundary condition
    const char onBoundary
      = std::count_if(sepSaddles.begin(), sepSaddles.end(),
                      [&triangulation](const SimplexId a) {
                        return triangulation.isEdgeOnBoundary(a);
                      });

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];

      // first store the SimplexId of the cell/triangle vertices
      SimplexId v0{}, v1{}, v2{};
      triangulation.getTriangleVertex(cell.id_, 0, v0);
      triangulation.getTriangleVertex(cell.id_, 1, v1);
      triangulation.getTriangleVertex(cell.id_, 2, v2);

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j;
      // index of current cell among all new cells
      const auto m = l - noldcells;

      cellsConn[3 * m + 0] = v0;
      cellsConn[3 * m + 1] = v1;
      cellsConn[3 * m + 2] = v2;
      cellVertsIds[3 * m + 0] = v0;
      cellVertsIds[3 * m + 1] = v1;
      cellVertsIds[3 * m + 2] = v2;

      if(outputSeparatrices2_cells_sourceIds_ != nullptr)
        (*outputSeparatrices2_cells_sourceIds_)[l] = src.id_;
      if(outputSeparatrices2_cells_separatrixIds_ != nullptr)
        (*outputSeparatrices2_cells_separatrixIds_)[l] = sepId;
      if(outputSeparatrices2_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices2_cells_separatrixTypes_)[l] = sepType;
      if(outputSeparatrices2_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices2_cells_isOnBoundary_)[l] = onBoundary;
    }
  }

  // reduce the cell vertices ids
  // (cells are triangles sharing two vertices)
  PSORT(this->threadNumber_)(cellVertsIds.begin(), cellVertsIds.end());
  const auto last = std::unique(cellVertsIds.begin(), cellVertsIds.end());
  cellVertsIds.erase(last, cellVertsIds.end());

  // vertex Id to index in points array
  std::vector<size_t> vertId2PointsId(triangulation.getNumberOfVertices());

  const auto noldpoints{npoints};
  npoints += cellVertsIds.size();
  outputSeparatrices2_points_->resize(3 * npoints);
  auto points = &outputSeparatrices2_points_->at(3 * noldpoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < cellVertsIds.size(); ++i) {
    // vertex 3D coords
    triangulation.getVertexPoint(
      cellVertsIds[i], points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    // vertex index in cellVertsIds array (do not forget offset)
    vertId2PointsId[cellVertsIds[i]] = i + noldpoints;
  }

  const auto lastOffset = noldcells == 0 ? 0 : cellsOff[-1];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < ncells - noldcells; ++i) {
    cellsOff[i] = 3 * i + lastOffset;
    cellsConn[3 * i + 0] = vertId2PointsId[cellsConn[3 * i + 0]];
    cellsConn[3 * i + 1] = vertId2PointsId[cellsConn[3 * i + 1]];
    cellsConn[3 * i + 2] = vertId2PointsId[cellsConn[3 * i + 2]];
  }

  cellsOff[ncells - noldcells] = cellsOff[ncells - noldcells - 1] + 3;

  (*outputSeparatrices2_numberOfPoints_) = npoints;
  (*outputSeparatrices2_numberOfCells_) = ncells;

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplex3D::execute(const triangulationType &triangulation) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField_) {
    this->printErr("Input scalar field pointer is null.");
    return -1;
  }

  if(!inputOffsets_) {
    this->printErr("Input offset field pointer is null.");
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
    discreteGradient_.buildGradient<triangulationType>(triangulation);

    this->printMsg("Discrete gradient computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  if(ReturnSaddleConnectors) {
    discreteGradient_.reverseGradient<dataType>(triangulation);
  }

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints, triangulation);

  std::vector<std::vector<Separatrix>> separatrices1{};
  std::vector<std::vector<std::vector<dcg::Cell>>> separatricesGeometry1;

  // 1-separatrices
  if(ComputeDescendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getDescendingSeparatrices1(criticalPoints, separatrices1.back(),
                               separatricesGeometry1.back(), triangulation);

    this->printMsg("Descending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices1) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getAscendingSeparatrices1(criticalPoints, separatrices1.back(),
                              separatricesGeometry1.back(), triangulation);

    this->printMsg("Ascending 1-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  // saddle-connectors
  if(ComputeSaddleConnectors) {
    Timer tmp;
    separatrices1.emplace_back();
    separatricesGeometry1.emplace_back();

    getSaddleConnectors(criticalPoints, separatrices1.back(),
                        separatricesGeometry1.back(), triangulation);

    this->printMsg("Saddle connectors computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  if(ComputeDescendingSeparatrices1 || ComputeAscendingSeparatrices1
     || ComputeSaddleConnectors) {
    Timer tmp{};

    flattenSeparatricesVectors(separatrices1, separatricesGeometry1);
    setSeparatrices1(separatrices1[0], separatricesGeometry1[0], triangulation);

    this->printMsg(
      "1-separatrices set", 1.0, tmp.getElapsedTime(), this->threadNumber_);
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices,
                               separatricesGeometry, separatricesSaddles,
                               triangulation);
    setDescendingSeparatrices2(
      separatrices, separatricesGeometry, separatricesSaddles, triangulation);

    this->printMsg("Descending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  if(ComputeAscendingSeparatrices2) {
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<SimplexId>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices,
                              separatricesGeometry, separatricesSaddles,
                              triangulation);
    setAscendingSeparatrices2(
      separatrices, separatricesGeometry, separatricesSaddles, triangulation);

    this->printMsg("Ascending 2-separatrices computed", 1.0,
                   tmp.getElapsedTime(), this->threadNumber_);
  }

  std::vector<SimplexId> maxSeeds;
  {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ascendingManifold)
      setAscendingSegmentation(criticalPoints, maxSeeds, ascendingManifold,
                               numberOfMaxima, triangulation);

    if(descendingManifold)
      setDescendingSegmentation(
        criticalPoints, descendingManifold, numberOfMinima, triangulation);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold,
                           descendingManifold, morseSmaleManifold,
                           triangulation);

    if(ascendingManifold or descendingManifold) {
      this->printMsg("Segmentation computed", 1.0, tmp.getElapsedTime(),
                     this->threadNumber_);
    }
  }

  if(outputCriticalPoints_points_ != nullptr) {
    std::vector<size_t> nCriticalPointsByDim;
    discreteGradient_.setCriticalPoints(
      criticalPoints, nCriticalPointsByDim, *outputCriticalPoints_points_,
      *outputCriticalPoints_points_cellDimensions_,
      *outputCriticalPoints_points_cellIds_,
      *outputCriticalPoints_points_isOnBoundary_,
      *outputCriticalPoints_points_PLVertexIdentifiers_, triangulation);

    if(ascendingManifold and descendingManifold) {
      discreteGradient_.setManifoldSize(
        criticalPoints, nCriticalPointsByDim, maxSeeds, ascendingManifold,
        descendingManifold, *outputCriticalPoints_points_manifoldSize_);
    }
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::getAscendingSeparatrices1(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  std::vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 2)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori :
  // numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const SimplexId numberOfSeparatrices = 4 * numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle = criticalPoints[saddleIndex];

    // add ascending vpaths
    {
      const Cell &saddle2 = saddle;

      const SimplexId starNumber
        = triangulation.getTriangleStarNumber(saddle2.id_);
      for(SimplexId j = 0; j < starNumber; ++j) {
        const SimplexId shift = j;

        SimplexId tetraId;
        triangulation.getTriangleStar(saddle2.id_, j, tetraId);

        std::vector<Cell> vpath;
        vpath.push_back(saddle2);
        discreteGradient_.getAscendingPath(
          Cell(3, tetraId), vpath, triangulation);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 3 and discreteGradient_.isCellCritical(lastCell)) {
          const SimplexId separatrixIndex = 4 * i + shift;

          separatricesGeometry[separatrixIndex] = std::move(vpath);
          separatrices[separatrixIndex]
            = Separatrix(true, saddle, lastCell, false, separatrixIndex);
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::getSaddleConnectors(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  const auto nTriangles = triangulation.getNumberOfTriangles();
  // visited triangles (one vector per thread)
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  std::vector<std::vector<SimplexId>> visitedTriangles(this->threadNumber_);

  for(auto &vec : isVisited) {
    // resize threads outside of main loop
    vec.resize(nTriangles, false);
  }

  // list of 2-saddles
  std::vector<Cell> saddles2{};
  // copy cells instead of taking a reference?
  std::copy_if(criticalPoints.begin(), criticalPoints.end(),
               std::back_inserter(saddles2),
               [](const Cell &c) -> bool { return c.dim_ == 2; });

  using Vpath = std::vector<Cell>;
  using SepSads = std::pair<Cell, Cell>;

  std::vector<std::vector<SepSads>> sepsByThread(saddles2.size());
  std::vector<std::vector<Vpath>> sepsGeomByThread(saddles2.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < saddles2.size(); ++i) {
    const auto &s2{saddles2[i]};

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    std::set<SimplexId> saddles1{};
    dcg::VisitedMask mask{isVisited[tid], visitedTriangles[tid]};
    discreteGradient_.getDescendingWall(
      s2, mask, triangulation, nullptr, &saddles1);

    for(const auto saddle1Id : saddles1) {
      const Cell s1{1, saddle1Id};

      Vpath vpath;
      const bool isMultiConnected
        = discreteGradient_.getAscendingPathThroughWall(
          s1, s2, isVisited[tid], &vpath, triangulation);
      const auto &last = vpath.back();

      if(!isMultiConnected && last.dim_ == s2.dim_ && last.id_ == s2.id_) {
        sepsGeomByThread[i].emplace_back(std::move(vpath));
        sepsByThread[i].emplace_back(s1, s2);
      }
    }
  }

  // count total number of separatrices in sepsByThread
  std::vector<size_t> partialSepsId(sepsByThread.size() + 1, 0);

  for(size_t i = 0; i < sepsByThread.size(); ++i) {
    partialSepsId[i + 1] = partialSepsId[i] + sepsByThread[i].size();
  }

  // pre-allocate output vectors
  separatrices.resize(partialSepsId.back());
  separatricesGeometry.resize(partialSepsId.back());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sepsByThread.size(); ++i) {
    for(size_t j = 0; j < sepsByThread[i].size(); ++j) {
      const auto &sads = sepsByThread[i][j];
      const size_t k = partialSepsId[i] + j;
      separatrices[k] = Separatrix{
        true, sads.first, sads.second, false, static_cast<SimplexId>(k)};
      separatricesGeometry[k] = std::move(sepsGeomByThread[i][j]);
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::getAscendingSeparatrices2(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  std::vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 1)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const SimplexId numberOfEdges = triangulation.getNumberOfEdges();
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  for(auto &vec : isVisited) {
    vec.resize(numberOfEdges, false);
  }
  std::vector<std::vector<SimplexId>> visitedEdges(this->threadNumber_);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle1 = criticalPoints[saddleIndex];

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    std::vector<Cell> wall;
    dcg::VisitedMask mask{isVisited[tid], visitedEdges[tid]};
    discreteGradient_.getAscendingWall(
      saddle1, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle1, emptyCell, false, i);
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::getDescendingSeparatrices2(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  std::vector<std::set<SimplexId>> &separatricesSaddles,
  const triangulationType &triangulation) const {
  const Cell emptyCell;

  std::vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 2)
      saddleIndexes.push_back(i);
  }
  const SimplexId numberOfSaddles = saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfWalls =
  // numberOfSaddles
  const SimplexId numberOfSeparatrices = numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);
  separatricesSaddles.resize(numberOfSeparatrices);

  const SimplexId numberOfTriangles = triangulation.getNumberOfTriangles();
  std::vector<std::vector<bool>> isVisited(this->threadNumber_);
  for(auto &vec : isVisited) {
    vec.resize(numberOfTriangles, false);
  }
  std::vector<std::vector<SimplexId>> visitedTriangles(this->threadNumber_);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < numberOfSaddles; ++i) {
    const SimplexId saddleIndex = saddleIndexes[i];
    const Cell &saddle2 = criticalPoints[saddleIndex];

#ifdef TTK_ENABLE_OPENMP
    const size_t tid = omp_get_thread_num();
#else
    const size_t tid = 0;
#endif // TTK_ENABLE_OPENMP

    std::vector<Cell> wall;
    dcg::VisitedMask mask{isVisited[tid], visitedTriangles[tid]};
    discreteGradient_.getDescendingWall(
      saddle2, mask, triangulation, &wall, &separatricesSaddles[i]);

    separatricesGeometry[i] = std::move(wall);
    separatrices[i] = Separatrix(true, saddle2, emptyCell, false, i);
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::getDualPolygon(
  const SimplexId edgeId,
  SimplexId *const polygon,
  const size_t polSize,
  const triangulationType &triangulation) const {

  for(size_t i = 0; i < polSize; ++i) {
    SimplexId starId;
    triangulation.getEdgeStar(edgeId, i, starId);
    polygon[i] = starId;
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::sortDualPolygonVertices(
  SimplexId *const polygon,
  const size_t polSize,
  const triangulationType &triangulation) const {

  for(size_t i = 1; i < polSize; ++i) {

    // find polygon[i - 1] neighboring tetra in polygon[i..]
    bool isFound = false;
    size_t j = i;
    for(; j < polSize; ++j) {
      // check if current is the neighbor
      for(SimplexId k = 0;
          k < triangulation.getCellNeighborNumber(polygon[i - 1]); ++k) {
        SimplexId neighborId{};
        triangulation.getCellNeighbor(polygon[i - 1], k, neighborId);
        if(neighborId == polygon[j]) {
          isFound = true;
          break;
        }
      }
      if(isFound)
        break;
    }

    // place polygon[j] next to polygon[i - 1]
    if(isFound) {
      std::swap(polygon[j], polygon[i]);
    }
  }

  return 0;
}
