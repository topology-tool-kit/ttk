#include <MorseSmaleComplex.h>

ttk::MorseSmaleComplex::MorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
}

void ttk::MorseSmaleComplex::flattenSeparatricesVectors(
  std::vector<std::vector<Separatrix>> &separatrices,
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

void ttk::MorseSmaleComplex::clear() {
  // critical points
  criticalPoints_points_.clear();
  criticalPoints_points_cellDimensions_.clear();
  criticalPoints_points_cellIds_.clear();
  criticalPoints_points_isOnBoundary_.clear();
  criticalPoints_points_PLVertexIdentifiers_.clear();
  criticalPoints_points_manifoldSize_.clear();

  // 1-separatrices
  separatrices1_numberOfPoints_ = separatrices1_numberOfCells_ = 0;
  separatrices1_points_.clear();
  separatrices1_points_smoothingMask_.clear();
  separatrices1_points_cellDimensions_.clear();
  separatrices1_points_cellIds_.clear();
  separatrices1_cells_connectivity_.clear();
  separatrices1_cells_sourceIds_.clear();
  separatrices1_cells_destinationIds_.clear();
  separatrices1_cells_separatrixIds_.clear();
  separatrices1_cells_separatrixTypes_.clear();
  separatrices1_cells_isOnBoundary_.clear();
  separatrices1_cells_sepFuncMinId_.clear();
  separatrices1_cells_sepFuncMaxId_.clear();

  // 2-separatrices
  separatrices2_numberOfPoints_ = separatrices2_numberOfCells_ = 0;
  separatrices2_points_.clear();
  separatrices2_cells_offsets_.clear();
  separatrices2_cells_connectivity_.clear();
  separatrices2_cells_sourceIds_.clear();
  separatrices2_cells_separatrixIds_.clear();
  separatrices2_cells_separatrixTypes_.clear();
  separatrices2_cells_isOnBoundary_.clear();
  separatrices2_cells_sepFuncMinId_.clear();
  separatrices2_cells_sepFuncMaxId_.clear();
}
