/// \ingroup base
/// \class ttk::AbstractMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %abstractMorseSmaleComplex processing package.
///
/// %AbstractMorseSmaleComplex is a processing package managing
/// the inputs, outputs and dependencies of the %MorseSmaleComplex2D
/// and %MorseSmaleComplex3D packages.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <Triangulation.h>

#include <queue>

namespace ttk {

  /**
   * Utility class representing Ridge lines, Valley lines
   * and Saddle connectors.
   */
  struct Separatrix {
    // default :
    explicit Separatrix()
      : isValid_{}, source_{}, destination_{}, isReversed_{}, geometry_{} {
    }

    // initialization with one segment :
    explicit Separatrix(const bool isValid,
                        const dcg::Cell &saddle,
                        const dcg::Cell &extremum,
                        const bool isSegmentReversed,
                        const SimplexId segmentGeometry)
      : isValid_{isValid}, source_{saddle}, destination_{extremum} {
      isReversed_.push_back(isSegmentReversed);
      geometry_.push_back(segmentGeometry);
    }

    // initialization with multiple segments :
    explicit Separatrix(const bool isValid,
                        const dcg::Cell &saddle,
                        const dcg::Cell &extremum,
                        const std::vector<char> &isReversed,
                        const std::vector<SimplexId> &geometry)
      : isValid_{isValid}, source_{saddle}, destination_{extremum},
        isReversed_{isReversed}, geometry_{geometry} {
    }

    explicit Separatrix(const Separatrix &separatrix)
      : isValid_{separatrix.isValid_}, source_{separatrix.source_},
        destination_{separatrix.destination_},
        isReversed_{separatrix.isReversed_}, geometry_{separatrix.geometry_} {
    }

    explicit Separatrix(Separatrix &&separatrix) noexcept
      : isValid_{separatrix.isValid_}, source_{separatrix.source_},
        destination_{separatrix.destination_}, isReversed_{std::move(
                                                 separatrix.isReversed_)},
        geometry_{std::move(separatrix.geometry_)} {
    }

    Separatrix &operator=(Separatrix &&separatrix) noexcept {
      isValid_ = separatrix.isValid_;
      source_ = separatrix.source_;
      destination_ = separatrix.destination_;
      isReversed_ = std::move(separatrix.isReversed_);
      geometry_ = std::move(separatrix.geometry_);

      return *this;
    }

    /**
     * Flag indicating if this separatrix can be processed.
     */
    bool isValid_;

    /**
     * Source cell of the separatrix.
     */
    dcg::Cell source_;

    /**
     * Destination cell of the separatrix.
     */
    dcg::Cell destination_;

    /**
     * Container of flags, isReversed[i] indicates if the
     * element stored at id=geometry_[i] can be reversed.
     */
    std::vector<char> isReversed_;

    /**
     * Container of ids. Each id addresses a separate
     * container corresponding to a dense representation
     * of the geometry (i.e. separatricesGeometry).
     */
    std::vector<SimplexId> geometry_;
  };

  /**
   * Parent class containing convenience functions shared between
   * Morse-Smale Complex algorithms for 2D and 3D domains.
   */
  class AbstractMorseSmaleComplex : virtual public Debug {

  public:
    AbstractMorseSmaleComplex();
    virtual ~AbstractMorseSmaleComplex();

    /**
     * Set the threshold for the iterative gradient reversal process.
     * Disable thresholding with -1 (default).
     */
    int setIterationThreshold(const int iterationThreshold) {
      discreteGradient_.setIterationThreshold(iterationThreshold);
      return 0;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the ascending manifolds of the critical points
     * (disabled by default).
     */
    int setComputeAscendingSeparatrices1(const bool state) {
      ComputeAscendingSeparatrices1 = state;
      return 0;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the descending manifolds of the critical points
     * (disabled by default).
     */
    int setComputeDescendingSeparatrices1(const bool state) {
      ComputeDescendingSeparatrices1 = state;
      return 0;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the visible saddle-connectors (disabled by default).
     */
    int setComputeSaddleConnectors(const bool state) {
      ComputeSaddleConnectors = state;
      return 0;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the ascending 2-separatrices (disabled by default).
     */
    int setComputeAscendingSeparatrices2(const bool state) {
      ComputeAscendingSeparatrices2 = state;
      return 0;
    }

    /**
     * Enable/Disable computation of the geometrical embedding of
     * the descending 2-separatrices (disabled by default).
     */
    int setComputeDescendingSeparatrices2(const bool state) {
      ComputeDescendingSeparatrices2 = state;
      return 0;
    }

    /**
     * Enable/Disable post-processing gradient reversal of
     * the (saddle,...,saddle) vpaths under a given persistence
     * threshold (disabled by default).
     */
    int setReturnSaddleConnectors(const bool state) {
      ReturnSaddleConnectors = state;
      discreteGradient_.setReturnSaddleConnectors(state);
      return 0;
    }

    /**
     * Set the threshold value for post-processing of
     * (saddle,...,saddle) vpaths gradient reversal
     * (default value is 0.0).
     */
    int setSaddleConnectorsPersistenceThreshold(const double threshold) {
      SaddleConnectorsPersistenceThreshold = threshold;
      discreteGradient_.setSaddleConnectorsPersistenceThreshold(threshold);
      return 0;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      discreteGradient_.preconditionTriangulation(data);
      data->preconditionCellEdges();
      data->preconditionCellNeighbors();
    }

    /**
     * Set the input scalar field associated on the points of the data set.
     */
    inline int setInputScalarField(const void *const data) {
      inputScalarField_ = data;
      discreteGradient_.setInputScalarField(inputScalarField_);
      return 0;
    }

    /**
     * Set the input offset field associated on the points of the data set
     * (if none, identifiers are used instead).
     *
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p data buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    inline int setInputOffsets(const SimplexId *const data) {
      inputOffsets_ = data;
      discreteGradient_.setInputOffsets(inputOffsets_);
      return 0;
    }

    /**
     * Set the output critical points data pointers.
     */
    inline int setOutputCriticalPoints(
      SimplexId *const criticalPoints_numberOfPoints,
      std::vector<float> *const criticalPoints_points,
      std::vector<char> *const criticalPoints_points_cellDimensons,
      std::vector<SimplexId> *const criticalPoints_points_cellIds,
      void *const criticalPoints_points_cellScalars,
      std::vector<char> *const criticalPoints_points_isOnBoundary,
      std::vector<SimplexId> *const criticalPoints_points_PLVertexIdentifiers,
      std::vector<SimplexId> *criticalPoints_points_manifoldSize) {
      outputCriticalPoints_numberOfPoints_ = criticalPoints_numberOfPoints;
      outputCriticalPoints_points_ = criticalPoints_points;
      outputCriticalPoints_points_cellDimensions_
        = criticalPoints_points_cellDimensons;
      outputCriticalPoints_points_cellIds_ = criticalPoints_points_cellIds;
      outputCriticalPoints_points_cellScalars_
        = criticalPoints_points_cellScalars;
      outputCriticalPoints_points_isOnBoundary_
        = criticalPoints_points_isOnBoundary;
      outputCriticalPoints_points_PLVertexIdentifiers_
        = criticalPoints_points_PLVertexIdentifiers;
      outputCriticalPoints_points_manifoldSize_
        = criticalPoints_points_manifoldSize;

      discreteGradient_.setOutputCriticalPoints(
        criticalPoints_points_cellScalars);
      return 0;
    }

    /**
     * Set the data pointers to the output 1-separatrices.
     */
    inline int setOutputSeparatrices1(
      SimplexId *const separatrices1_numberOfPoints,
      std::vector<float> *const separatrices1_points,
      std::vector<char> *const separatrices1_points_smoothingMask,
      std::vector<char> *const separatrices1_points_cellDimensions,
      std::vector<SimplexId> *const separatrices1_points_cellIds,
      SimplexId *const separatrices1_numberOfCells,
      std::vector<long long> *const separatrices1_cells_connectivity,
      std::vector<SimplexId> *const separatrices1_cells_sourceIds,
      std::vector<SimplexId> *const separatrices1_cells_destinationIds,
      std::vector<SimplexId> *const separatrices1_cells_separatrixIds,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      std::vector<double> *const separatrices1_cells_separatrixFunctionMaxima,
      std::vector<double> *const separatrices1_cells_separatrixFunctionMinima,
      std::vector<double> *const separatrices1_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices1_cells_isOnBoundary) {
      outputSeparatrices1_numberOfPoints_ = separatrices1_numberOfPoints;
      outputSeparatrices1_points_ = separatrices1_points;
      outputSeparatrices1_points_smoothingMask_
        = separatrices1_points_smoothingMask;
      outputSeparatrices1_points_cellDimensions_
        = separatrices1_points_cellDimensions;
      outputSeparatrices1_points_cellIds_ = separatrices1_points_cellIds;
      outputSeparatrices1_numberOfCells_ = separatrices1_numberOfCells;
      outputSeparatrices1_cells_connectivity_
        = separatrices1_cells_connectivity;
      outputSeparatrices1_cells_sourceIds_ = separatrices1_cells_sourceIds;
      outputSeparatrices1_cells_destinationIds_
        = separatrices1_cells_destinationIds;
      outputSeparatrices1_cells_separatrixIds_
        = separatrices1_cells_separatrixIds;
      outputSeparatrices1_cells_separatrixTypes_
        = separatrices1_cells_separatrixTypes;
      outputSeparatrices1_cells_separatrixFunctionMaxima_
        = separatrices1_cells_separatrixFunctionMaxima;
      outputSeparatrices1_cells_separatrixFunctionMinima_
        = separatrices1_cells_separatrixFunctionMinima;
      outputSeparatrices1_cells_separatrixFunctionDiffs_
        = separatrices1_cells_separatrixFunctionDiffs;
      outputSeparatrices1_cells_isOnBoundary_
        = separatrices1_cells_isOnBoundary;
      return 0;
    }

    /**
     * Set the data pointers to the output 2-separatrices.
     */
    inline int setOutputSeparatrices2(
      SimplexId *const separatrices2_numberOfPoints,
      std::vector<float> *const separatrices2_points,
      SimplexId *const separatrices2_numberOfCells,
      std::vector<long long> *const separatrices2_cells_offsets,
      std::vector<long long> *const separatrices2_cells_connectivity,
      std::vector<SimplexId> *const separatrices2_cells_sourceIds,
      std::vector<SimplexId> *const separatrices2_cells_separatrixIds,
      std::vector<char> *const separatrices2_cells_separatrixTypes,
      std::vector<double> *const separatrices2_cells_separatrixFunctionMaxima,
      std::vector<double> *const separatrices2_cells_separatrixFunctionMinima,
      std::vector<double> *const separatrices2_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices2_cells_isOnBoundary) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_offsets_ = separatrices2_cells_offsets;
      outputSeparatrices2_cells_connectivity_
        = separatrices2_cells_connectivity;
      outputSeparatrices2_cells_sourceIds_ = separatrices2_cells_sourceIds;
      outputSeparatrices2_cells_separatrixIds_
        = separatrices2_cells_separatrixIds;
      outputSeparatrices2_cells_separatrixTypes_
        = separatrices2_cells_separatrixTypes;
      outputSeparatrices2_cells_separatrixFunctionMaxima_
        = separatrices2_cells_separatrixFunctionMaxima;
      outputSeparatrices2_cells_separatrixFunctionMinima_
        = separatrices2_cells_separatrixFunctionMinima;
      outputSeparatrices2_cells_separatrixFunctionDiffs_
        = separatrices2_cells_separatrixFunctionDiffs;
      outputSeparatrices2_cells_isOnBoundary_
        = separatrices2_cells_isOnBoundary;
      return 0;
    }

    /**
     * Set the data pointers to the output segmentation scalar fields.
     */
    inline int setOutputMorseComplexes(void *const ascendingManifold,
                                       void *const descendingManifold,
                                       void *const morseSmaleManifold) {
      outputAscendingManifold_ = ascendingManifold;
      outputDescendingManifold_ = descendingManifold;
      outputMorseSmaleManifold_ = morseSmaleManifold;
      return 0;
    }

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getDescendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the geometrical embedding of the 1-separatrices. This
     * function needs the following internal pointers to be set:
     * outputSeparatrices1_numberOfPoints_
     * outputSeparatrices1_points_
     * outputSeparatrices1_numberOfCells_
     * outputSeparatrices1_cells_
     * inputScalarField_
     */
    template <typename dataType, typename triangulationType>
    int setSeparatrices1(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;

    /**
     * Compute the ascending manifold of the maxima.
     */
    template <typename triangulationType>
    int setAscendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                 std::vector<SimplexId> &maxSeeds,
                                 SimplexId *const morseSmaleManifold,
                                 SimplexId &numberOfMaxima,
                                 const triangulationType &triangulation) const;

    /**
     * Compute the descending manifold of the minima.
     */
    template <typename triangulationType>
    int setDescendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                  SimplexId *const morseSmaleManifold,
                                  SimplexId &numberOfMinima,
                                  const triangulationType &triangulation) const;

    /**
     * Compute the final combinatorial Morse-Smale complex
     * segmentation.
     */
    template <typename triangulationType>
    int setFinalSegmentation(const SimplexId numberOfMaxima,
                             const SimplexId numberOfMinima,
                             const SimplexId *const ascendingManifold,
                             const SimplexId *const descendingManifold,
                             SimplexId *const morseSmaleManifold,
                             const triangulationType &triangulation) const;

  protected:
    bool ReverveSaddleMaximumConnection{false};
    bool ReverveSaddleSaddleConnection{false};
    bool ComputeAscendingSeparatrices1{true};
    bool ComputeDescendingSeparatrices1{true};
    bool ComputeSaddleConnectors{false};
    bool ComputeAscendingSeparatrices2{false};
    bool ComputeDescendingSeparatrices2{false};
    bool ReturnSaddleConnectors{false};
    double SaddleConnectorsPersistenceThreshold{};

    dcg::DiscreteGradient discreteGradient_{};

    const void *inputScalarField_{};
    const SimplexId *inputOffsets_{};

    SimplexId *outputCriticalPoints_numberOfPoints_{};
    std::vector<float> *outputCriticalPoints_points_{};
    std::vector<char> *outputCriticalPoints_points_cellDimensions_{};
    std::vector<SimplexId> *outputCriticalPoints_points_cellIds_{};
    void *outputCriticalPoints_points_cellScalars_{};
    std::vector<char> *outputCriticalPoints_points_isOnBoundary_{};
    std::vector<SimplexId> *outputCriticalPoints_points_PLVertexIdentifiers_{};
    std::vector<SimplexId> *outputCriticalPoints_points_manifoldSize_{};

    SimplexId *outputSeparatrices1_numberOfPoints_{};
    std::vector<float> *outputSeparatrices1_points_{};
    std::vector<char> *outputSeparatrices1_points_smoothingMask_{};
    std::vector<char> *outputSeparatrices1_points_cellDimensions_{};
    std::vector<SimplexId> *outputSeparatrices1_points_cellIds_{};
    SimplexId *outputSeparatrices1_numberOfCells_{};
    std::vector<long long> *outputSeparatrices1_cells_connectivity_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_sourceIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_destinationIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_separatrixIds_{};
    std::vector<char> *outputSeparatrices1_cells_separatrixTypes_{};
    std::vector<double> *outputSeparatrices1_cells_separatrixFunctionMaxima_{};
    std::vector<double> *outputSeparatrices1_cells_separatrixFunctionMinima_{};
    std::vector<double> *outputSeparatrices1_cells_separatrixFunctionDiffs_{};
    std::vector<char> *outputSeparatrices1_cells_isOnBoundary_{};

    SimplexId *outputSeparatrices2_numberOfPoints_{};
    std::vector<float> *outputSeparatrices2_points_{};
    SimplexId *outputSeparatrices2_numberOfCells_{};
    std::vector<long long> *outputSeparatrices2_cells_offsets_{};
    std::vector<long long> *outputSeparatrices2_cells_connectivity_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_sourceIds_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_separatrixIds_{};
    std::vector<char> *outputSeparatrices2_cells_separatrixTypes_{};
    std::vector<double> *outputSeparatrices2_cells_separatrixFunctionMaxima_{};
    std::vector<double> *outputSeparatrices2_cells_separatrixFunctionMinima_{};
    std::vector<double> *outputSeparatrices2_cells_separatrixFunctionDiffs_{};
    std::vector<char> *outputSeparatrices2_cells_isOnBoundary_{};

    void *outputAscendingManifold_{};
    void *outputDescendingManifold_{};
    void *outputMorseSmaleManifold_{};
  };
} // namespace ttk

template <typename triangulationType>
int ttk::AbstractMorseSmaleComplex::getDescendingSeparatrices1(
  const std::vector<Cell> &criticalPoints,
  std::vector<Separatrix> &separatrices,
  std::vector<std::vector<Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

  std::vector<SimplexId> saddleIndexes;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 1)
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

    // add descending vpaths
    {
      const Cell &saddle1 = saddle;

      for(int j = 0; j < 2; ++j) {
        const int shift = j + 2;

        SimplexId vertexId;
        triangulation.getEdgeVertex(saddle1.id_, j, vertexId);

        std::vector<Cell> vpath;
        vpath.push_back(saddle1);
        discreteGradient_.getDescendingPath(
          Cell(0, vertexId), vpath, triangulation);

        const Cell &lastCell = vpath.back();
        if(lastCell.dim_ == 0 and discreteGradient_.isCellCritical(lastCell)) {
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

template <typename dataType, typename triangulationType>
int ttk::AbstractMorseSmaleComplex::setSeparatrices1(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
  const triangulationType &triangulation) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices1_numberOfPoints_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices1_points_ == nullptr) {
    this->printErr("1-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices1_numberOfCells_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices1_cells_connectivity_ == nullptr) {
    this->printErr("1-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputScalarField_ == nullptr) {
    this->printErr(
      " 1-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  const auto scalars = static_cast<const dataType *>(inputScalarField_);
  auto separatrixFunctionMaxima
    = outputSeparatrices1_cells_separatrixFunctionMaxima_;
  auto separatrixFunctionMinima
    = outputSeparatrices1_cells_separatrixFunctionMinima_;
  auto separatrixFunctionDiffs
    = outputSeparatrices1_cells_separatrixFunctionDiffs_;

  // max existing separatrix id + 1 or 0
  const SimplexId separatrixId
    = (outputSeparatrices1_cells_separatrixIds_ != nullptr
       && !outputSeparatrices1_cells_separatrixIds_->empty())
        ? *std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
                            outputSeparatrices1_cells_separatrixIds_->end())
            + 1
        : 0;

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices1_numberOfPoints_)};
  // total number of separatrices cells
  auto ncells{static_cast<size_t>(*outputSeparatrices1_numberOfCells_)};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // corresponding separatrix index in separatrices array
  std::vector<SimplexId> geomIdSep{};
  // points beginning id for each separatrix geometry
  std::vector<size_t> geomPointsBegId{npoints};
  // cells beginning id for each separatrix geometry
  std::vector<size_t> geomCellsBegId{ncells};

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < separatrices.size(); ++i) {
    const auto &sep = separatrices[i];
    if(!sep.isValid_ || sep.geometry_.empty()) {
      continue;
    }
    for(const auto geomId : sep.geometry_) {
      const auto sepSize = separatricesGeometry[geomId].size();
      npoints += sepSize;
      ncells += sepSize - 1;
      geomPointsBegId.emplace_back(npoints);
      geomCellsBegId.emplace_back(ncells);
      validGeomIds.emplace_back(geomId);
      geomIdSep.emplace_back(i);
    }
  }

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;

  // resize arrays
  outputSeparatrices1_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize(2 * ncells);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  if(outputSeparatrices1_points_smoothingMask_ != nullptr)
    outputSeparatrices1_points_smoothingMask_->resize(npoints);
  if(outputSeparatrices1_points_cellDimensions_ != nullptr)
    outputSeparatrices1_points_cellDimensions_->resize(npoints);
  if(outputSeparatrices1_points_cellIds_ != nullptr)
    outputSeparatrices1_points_cellIds_->resize(npoints);
  if(outputSeparatrices1_cells_sourceIds_ != nullptr)
    outputSeparatrices1_cells_sourceIds_->resize(ncells);
  if(outputSeparatrices1_cells_destinationIds_ != nullptr)
    outputSeparatrices1_cells_destinationIds_->resize(ncells);
  if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
    outputSeparatrices1_cells_separatrixIds_->resize(ncells);
  if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
    outputSeparatrices1_cells_separatrixTypes_->resize(ncells);
  if(separatrixFunctionMaxima != nullptr)
    separatrixFunctionMaxima->resize(ncells);
  if(separatrixFunctionMinima != nullptr)
    separatrixFunctionMinima->resize(ncells);
  if(separatrixFunctionDiffs != nullptr)
    separatrixFunctionDiffs->resize(ncells);
  if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
    outputSeparatrices1_cells_isOnBoundary_->resize(ncells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < validGeomIds.size(); ++i) {
    const auto &sep = separatrices[geomIdSep[i]];
    const auto &sepGeom = separatricesGeometry[validGeomIds[i]];
    const auto sepId = separatrixId + i;
    // saddle (asc/desc sep) or saddle1 (saddle connector)
    const dcg::Cell &src = sep.source_;
    // extremum (asc/desc sep) or saddle2 (saddle connector)
    const dcg::Cell &dst = sep.destination_;

    // get separatrix type
    const auto saddleConnector
      = dimensionality == 3 && src.dim_ == 1 && dst.dim_ == 2;
    const char sepType
      = saddleConnector ? 1 : std::min(dst.dim_, dimensionality - 1);

    // compute separatrix function diff
    const auto sepFuncMax = static_cast<double>(std::max(
      scalars[discreteGradient_.getCellGreaterVertex(src, triangulation)],
      scalars[discreteGradient_.getCellGreaterVertex(dst, triangulation)]));
    const auto sepFuncMin = static_cast<double>(std::min(
      scalars[discreteGradient_.getCellLowerVertex(src, triangulation)],
      scalars[discreteGradient_.getCellLowerVertex(dst, triangulation)]));
    const auto sepFuncDiff = sepFuncMax - sepFuncMin;

    // get boundary condition
    const auto onBoundary
      = saddleConnector
          ? static_cast<char>(
            discreteGradient_.isBoundary(src, triangulation)
            && discreteGradient_.isBoundary(dst, triangulation))
          : static_cast<char>(discreteGradient_.isBoundary(src, triangulation))
              + static_cast<char>(
                discreteGradient_.isBoundary(dst, triangulation));

    for(size_t j = 0; j < sepGeom.size(); ++j) {
      const auto &cell = sepGeom[j];
      std::array<float, 3> pt{};
      triangulation.getCellIncenter(cell.id_, cell.dim_, pt.data());

      // index of current point in point data arrays
      const auto k = geomPointsBegId[i] + j;

      points[3 * k + 0] = pt[0];
      points[3 * k + 1] = pt[1];
      points[3 * k + 2] = pt[2];

      if(outputSeparatrices1_points_smoothingMask_ != nullptr)
        (*outputSeparatrices1_points_smoothingMask_)[k]
          = (j == 0 || j == sepGeom.size() - 1) ? 0 : 1;
      if(outputSeparatrices1_points_cellDimensions_ != nullptr)
        (*outputSeparatrices1_points_cellDimensions_)[k] = cell.dim_;
      if(outputSeparatrices1_points_cellIds_ != nullptr)
        (*outputSeparatrices1_points_cellIds_)[k] = cell.id_;

      // skip filling cell data for first geometry point
      if(j == 0)
        continue;

      // index of current cell in cell data arrays
      const auto l = geomCellsBegId[i] + j - 1;

      cellsConn[2 * l + 0] = k - 1;
      cellsConn[2 * l + 1] = k;

      if(outputSeparatrices1_cells_sourceIds_ != nullptr)
        (*outputSeparatrices1_cells_sourceIds_)[l] = src.id_;
      if(outputSeparatrices1_cells_destinationIds_ != nullptr)
        (*outputSeparatrices1_cells_destinationIds_)[l] = dst.id_;
      if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
        (*outputSeparatrices1_cells_separatrixIds_)[l] = sepId;
      if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices1_cells_separatrixTypes_)[l] = sepType;
      if(separatrixFunctionMaxima != nullptr)
        (*separatrixFunctionMaxima)[l] = sepFuncMax;
      if(separatrixFunctionDiffs != nullptr)
        (*separatrixFunctionMinima)[l] = sepFuncMin;
      if(separatrixFunctionDiffs != nullptr)
        (*separatrixFunctionDiffs)[l] = sepFuncDiff;
      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[l] = onBoundary;
    }
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = ncells;

  return 0;
}

template <typename triangulationType>
int ttk::AbstractMorseSmaleComplex::setAscendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  std::vector<SimplexId> &maxSeeds,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation) const {

  const SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  const SimplexId numberOfCells = triangulation.getNumberOfCells();
  std::vector<SimplexId> morseSmaleManifoldOnCells(numberOfCells, -1);
  const int cellDim = triangulation.getDimensionality();

  // get the seeds : maxima
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == cellDim)
      maxSeeds.push_back(criticalPoint.id_);
  }
  const SimplexId numberOfSeeds = maxSeeds.size();
  numberOfMaxima = numberOfSeeds;

  // Triangulation method pointers for 3D
  auto getCellFace = &triangulationType::getCellTriangle;
  auto getFaceStarNumber = &triangulationType::getTriangleStarNumber;
  auto getFaceStar = &triangulationType::getTriangleStar;
  if(cellDim == 2) {
    // Triangulation method pointers for 2D
    getCellFace = &triangulationType::getCellEdge;
    getFaceStarNumber = &triangulationType::getEdgeStarNumber;
    getFaceStar = &triangulationType::getEdgeStar;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(SimplexId i = 0; i < numberOfSeeds; ++i) {
    std::queue<SimplexId> bfs;

    // push the seed
    {
      const SimplexId seedId = maxSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()) {
      const SimplexId cofacetId = bfs.front();
      bfs.pop();

      if(morseSmaleManifoldOnCells[cofacetId] == -1) {
        morseSmaleManifoldOnCells[cofacetId] = i;

        for(int j = 0; j < (cellDim + 1); ++j) {
          SimplexId facetId = -1;
          (triangulation.*getCellFace)(cofacetId, j, facetId);

          SimplexId starNumber = (triangulation.*getFaceStarNumber)(facetId);
          for(SimplexId k = 0; k < starNumber; ++k) {
            SimplexId neighborId = -1;
            (triangulation.*getFaceStar)(facetId, k, neighborId);

            if(neighborId == cofacetId) {
              continue;
            }

            const SimplexId pairedCellId = discreteGradient_.getPairedCell(
              Cell(cellDim, neighborId), triangulation, true);

            if(pairedCellId == facetId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  // put segmentation infos from cells to points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < numberOfVertices; ++i) {
    SimplexId starId;
    triangulation.getVertexStar(i, 0, starId);
    morseSmaleManifold[i] = morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

template <typename triangulationType>
int ttk::AbstractMorseSmaleComplex::setDescendingSegmentation(
  const std::vector<Cell> &criticalPoints,
  SimplexId *const morseSmaleManifold,
  SimplexId &numberOfMinima,
  const triangulationType &triangulation) const {

  const SimplexId numberOfVertices = triangulation.getNumberOfVertices();
  std::fill(morseSmaleManifold, morseSmaleManifold + numberOfVertices, -1);

  // get the seeds : minima
  std::vector<SimplexId> seeds;
  const SimplexId numberOfCriticalPoints = criticalPoints.size();
  for(SimplexId i = 0; i < numberOfCriticalPoints; ++i) {
    const Cell &criticalPoint = criticalPoints[i];

    if(criticalPoint.dim_ == 0)
      seeds.push_back(criticalPoint.id_);
  }
  const SimplexId numberOfSeeds = seeds.size();
  numberOfMinima = numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(SimplexId i = 0; i < numberOfSeeds; ++i) {
    std::queue<SimplexId> bfs;

    // push the seed
    {
      const SimplexId seedId = seeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()) {
      const SimplexId vertexId = bfs.front();
      bfs.pop();

      if(morseSmaleManifold[vertexId] == -1) {
        morseSmaleManifold[vertexId] = i;

        const SimplexId edgeNumber
          = triangulation.getVertexEdgeNumber(vertexId);
        for(SimplexId j = 0; j < edgeNumber; ++j) {
          SimplexId edgeId;
          triangulation.getVertexEdge(vertexId, j, edgeId);

          for(int k = 0; k < 2; ++k) {
            SimplexId neighborId;
            triangulation.getEdgeVertex(edgeId, k, neighborId);

            if(neighborId == vertexId) {
              continue;
            }

            const SimplexId pairedCellId = discreteGradient_.getPairedCell(
              Cell(0, neighborId), triangulation);

            if(pairedCellId == edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

template <typename triangulationType>
int ttk::AbstractMorseSmaleComplex::setFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  const size_t nVerts = triangulation.getNumberOfVertices();

  // associate a unique "sparse region id" to each (ascending, descending) pair

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    const auto d = ascendingManifold[i];
    const auto a = descendingManifold[i];
    if(a == -1 || d == -1) {
      morseSmaleManifold[i] = -1;
    } else {
      morseSmaleManifold[i] = a * numberOfMaxima + d;
    }
  }

  // store the "sparse region ids" by copying the morseSmaleManifold output
  std::vector<SimplexId> sparseRegionIds(
    morseSmaleManifold, morseSmaleManifold + nVerts);

  // get unique "sparse region ids"
  PSORT(this->threadNumber_)(sparseRegionIds.begin(), sparseRegionIds.end());
  const auto last = std::unique(sparseRegionIds.begin(), sparseRegionIds.end());
  sparseRegionIds.erase(last, sparseRegionIds.end());

  // "sparse region id" -> "dense region id"
  std::map<SimplexId, size_t> sparseToDenseRegionId{};

  for(size_t i = 0; i < sparseRegionIds.size(); ++i) {
    sparseToDenseRegionId[sparseRegionIds[i]] = i;
  }

  // update region id on all vertices: "sparse id" -> "dense id"

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = sparseToDenseRegionId[morseSmaleManifold[i]];
  }

  return 0;
}
