#include <MorseSmaleComplex3D.h>
#include <iterator>

using namespace ttk;

MorseSmaleComplex3D::MorseSmaleComplex3D() : AbstractMorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex3D");
}

MorseSmaleComplex3D::~MorseSmaleComplex3D() {
}

template <typename dataType, typename triangulationType>
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
  if(outputSeparatrices2_cells_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
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
    const dataType sepFuncMin
      = discreteGradient_.scalarMin(src, scalars, triangulation);
    const auto maxId = *std::max_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, scalars, this](const SimplexId a, const SimplexId b) {
        return discreteGradient_.scalarMax(Cell{2, a}, scalars, triangulation)
               < discreteGradient_.scalarMax(
                 Cell{2, b}, scalars, triangulation);
      });
    const dataType sepFuncMax
      = discreteGradient_.scalarMax(Cell{2, maxId}, scalars, triangulation);

    // get boundary condition
    const char onBoundary
      = std::count_if(sepSaddles.begin(), sepSaddles.end(),
                      [&triangulation](const SimplexId a) {
                        return triangulation.isEdgeOnBoundary(a);
                      });

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      // index of current cell in cell data arrays
      const auto k = geomCellsBegId[i] + j - noldcells;
      auto &polyCell = polygonTetras[k];

      // Transform to dual : edge -> polygon
      getDualPolygon(cell.id_, polyCell.tetras_, triangulation);

      if(polyCell.tetras_.size() > 2) {
        sortDualPolygonVertices(polyCell.tetras_, triangulation);
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
  std::vector<size_t> vertId2PointsId(triangulation.getNumberOfCells());

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
    triangulation.getTetraIncenter(cellVertsIds[i], &points[3 * i]);
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

template <typename dataType, typename triangulationType>
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
  if(outputSeparatrices2_cells_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr("2-separatrices pointer to the input scalar field is null.");
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
    const dataType sepFuncMax
      = discreteGradient_.scalarMax(src, scalars, triangulation);
    const auto minId = *std::min_element(
      sepSaddles.begin(), sepSaddles.end(),
      [&triangulation, scalars, this](const SimplexId a, const SimplexId b) {
        return discreteGradient_.scalarMin(Cell{1, a}, scalars, triangulation)
               < discreteGradient_.scalarMin(
                 Cell{1, b}, scalars, triangulation);
      });
    const dataType sepFuncMin
      = discreteGradient_.scalarMin(Cell{1, minId}, scalars, triangulation);
    const dataType sepFuncDiff = sepFuncMax - sepFuncMin;

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

template <typename dataType, typename idType, typename triangulationType>
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
    discreteGradient_.buildGradient<dataType, idType>(triangulation);

    this->printMsg("Discrete gradient computed", 1.0, tmp.getElapsedTime(),
                   this->threadNumber_);
  }

  if(ReturnSaddleConnectors) {
    discreteGradient_.reverseGradient<dataType, idType>(triangulation);
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
    setSeparatrices1<dataType>(
      separatrices1[0], separatricesGeometry1[0], triangulation);

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
    setDescendingSeparatrices2<dataType>(
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
    setAscendingSeparatrices2<dataType>(
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

  if(outputCriticalPoints_numberOfPoints_ and outputSeparatrices1_points_) {
    std::vector<size_t> nCriticalPointsByDim{};
    discreteGradient_.setCriticalPoints<dataType>(
      criticalPoints, nCriticalPointsByDim, triangulation);

    discreteGradient_.fetchOutputCriticalPoints(
      outputCriticalPoints_numberOfPoints_, outputCriticalPoints_points_,
      outputCriticalPoints_points_cellDimensions_,
      outputCriticalPoints_points_cellIds_,
      outputCriticalPoints_points_isOnBoundary_,
      outputCriticalPoints_points_PLVertexIdentifiers_);

    if(ascendingManifold and descendingManifold) {
      discreteGradient_.setManifoldSize(criticalPoints, nCriticalPointsByDim,
                                        maxSeeds, ascendingManifold,
                                        descendingManifold);
      discreteGradient_.fetchOutputManifoldSize(
        outputCriticalPoints_points_manifoldSize_);
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
  std::vector<SimplexId> &polygon,
  const triangulationType &triangulation) const {

  polygon.resize(triangulation.getEdgeStarNumber(edgeId));
  for(size_t i = 0; i < polygon.size(); ++i) {
    SimplexId starId;
    triangulation.getEdgeStar(edgeId, i, starId);
    polygon[i] = starId;
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplex3D::sortDualPolygonVertices(
  std::vector<SimplexId> &polygon,
  const triangulationType &triangulation) const {

  for(size_t i = 1; i < polygon.size(); ++i) {

    // find polygon[i - 1] neighboring tetra in polygon[i..]
    bool isFound = false;
    size_t j = i;
    for(; j < polygon.size(); ++j) {
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

void MorseSmaleComplex3D::flattenSeparatricesVectors(
  std::vector<std::vector<ttk::Separatrix>> &separatrices,
  std::vector<std::vector<std::vector<ttk::dcg::Cell>>> &separatricesGeometry)
  const {

  std::vector<size_t> partialSizes{0};
  for(const auto &sep : separatrices) {
    partialSizes.emplace_back(partialSizes.back() + sep.size());
  }
  separatrices[0].resize(partialSizes.back());
  separatricesGeometry[0].resize(partialSizes.back());

  for(size_t i = 1; i < separatrices.size(); ++i) {
    const auto offset = partialSizes[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < separatrices[i].size(); ++j) {
      // shift separatrices geometry_
      for(size_t k = 0; k < separatrices[i][j].geometry_.size(); ++k) {
        separatrices[i][j].geometry_[k] += offset;
      }
      // flatten separatrices1 and separatricesGeometry1
      separatrices[0][offset + j] = std::move(separatrices[i][j]);
      separatricesGeometry[0][offset + j]
        = std::move(separatricesGeometry[i][j]);
    }
  }
}

// explicit specializations to reduce compile time
#define MSC3D_SPECIALIZE(DATATYPE)                                        \
  template int ttk::MorseSmaleComplex3D::execute<DATATYPE, SimplexId,     \
                                                 ExplicitTriangulation>(  \
    const ExplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex3D::execute<DATATYPE, SimplexId,     \
                                                 ImplicitTriangulation>(  \
    const ImplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex3D::execute<                         \
    DATATYPE, SimplexId, PeriodicImplicitTriangulation>(                  \
    const PeriodicImplicitTriangulation &);                               \
  template int ttk::MorseSmaleComplex3D::execute<DATATYPE, LongSimplexId, \
                                                 ExplicitTriangulation>(  \
    const ExplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex3D::execute<DATATYPE, LongSimplexId, \
                                                 ImplicitTriangulation>(  \
    const ImplicitTriangulation &);                                       \
  template int ttk::MorseSmaleComplex3D::execute<                         \
    DATATYPE, LongSimplexId, PeriodicImplicitTriangulation>(              \
    const PeriodicImplicitTriangulation &)

MSC3D_SPECIALIZE(double);
MSC3D_SPECIALIZE(float);
MSC3D_SPECIALIZE(long long);
MSC3D_SPECIALIZE(unsigned long long);
MSC3D_SPECIALIZE(long);
MSC3D_SPECIALIZE(unsigned long);
MSC3D_SPECIALIZE(int);
MSC3D_SPECIALIZE(unsigned int);
MSC3D_SPECIALIZE(short);
MSC3D_SPECIALIZE(unsigned short);
MSC3D_SPECIALIZE(char);
MSC3D_SPECIALIZE(signed char);
MSC3D_SPECIALIZE(unsigned char);
