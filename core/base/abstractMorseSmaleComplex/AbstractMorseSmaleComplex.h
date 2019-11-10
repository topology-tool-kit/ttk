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

#ifndef _ABSTRACTMORSESMALECOMPLEX_H
#define _ABSTRACTMORSESMALECOMPLEX_H

// base code includes
#include <DiscreteGradient.h>
#include <Triangulation.h>
#include <Wrapper.h>

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

    explicit Separatrix(Separatrix &&separatrix)
      : isValid_{separatrix.isValid_}, source_{std::move(separatrix.source_)},
        destination_{std::move(separatrix.destination_)},
        isReversed_{std::move(separatrix.isReversed_)},
        geometry_{std::move(separatrix.geometry_)} {
    }

    Separatrix &operator=(Separatrix &&separatrix) {
      isValid_ = separatrix.isValid_;
      source_ = std::move(separatrix.source_);
      destination_ = std::move(separatrix.destination_);
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
  class AbstractMorseSmaleComplex : public Debug {

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
     * Enable/Disable gradient reversal of (saddle,...,maximum) vpaths
     * (disabled by default).
     */
    int setReverveSaddleMaximumConnection(const bool state) {
      discreteGradient_.setReverseSaddleMaximumConnection(state);
      return 0;
    }

    /**
     * Enable/Disable gradient reversal of (saddle,...,saddle) vpaths
     * (disabled by default).
     */
    int setReverveSaddleSaddleConnection(const bool state) {
      discreteGradient_.setReverseSaddleSaddleConnection(state);
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
     * Enable/Disable compromise on memory allowing
     * computation the geometry of the 2-separatrices
     * in parallel.
     */
    int setPrioritizeSpeedOverMemory(const bool state) {
      PrioritizeSpeedOverMemory = state;
      return 0;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline int setupTriangulation(Triangulation *const data) {
      inputTriangulation_ = data;
      discreteGradient_.setupTriangulation(inputTriangulation_);

      inputTriangulation_->preprocessCellEdges();
      inputTriangulation_->preprocessCellNeighbors();
      return 0;
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
     */
    inline int setInputOffsets(const void *const data) {
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
        criticalPoints_numberOfPoints, criticalPoints_points,
        criticalPoints_points_cellDimensons, criticalPoints_points_cellIds,
        criticalPoints_points_cellScalars, criticalPoints_points_isOnBoundary,
        criticalPoints_points_PLVertexIdentifiers,
        criticalPoints_points_manifoldSize);
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
      std::vector<SimplexId> *const separatrices1_cells,
      std::vector<SimplexId> *const separatrices1_cells_sourceIds,
      std::vector<SimplexId> *const separatrices1_cells_destinationIds,
      std::vector<SimplexId> *const separatrices1_cells_separatrixIds,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      void *const separatrices1_cells_separatrixFunctionMaxima,
      void *const separatrices1_cells_separatrixFunctionMinima,
      void *const separatrices1_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices1_cells_isOnBoundary) {
      outputSeparatrices1_numberOfPoints_ = separatrices1_numberOfPoints;
      outputSeparatrices1_points_ = separatrices1_points;
      outputSeparatrices1_points_smoothingMask_
        = separatrices1_points_smoothingMask;
      outputSeparatrices1_points_cellDimensions_
        = separatrices1_points_cellDimensions;
      outputSeparatrices1_points_cellIds_ = separatrices1_points_cellIds;
      outputSeparatrices1_numberOfCells_ = separatrices1_numberOfCells;
      outputSeparatrices1_cells_ = separatrices1_cells;
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
      std::vector<SimplexId> *const separatrices2_cells,
      std::vector<SimplexId> *const separatrices2_cells_sourceIds,
      std::vector<SimplexId> *const separatrices2_cells_separatrixIds,
      std::vector<char> *const separatrices2_cells_separatrixTypes,
      void *const separatrices2_cells_separatrixFunctionMaxima,
      void *const separatrices2_cells_separatrixFunctionMinima,
      void *const separatrices2_cells_separatrixFunctionDiffs,
      std::vector<char> *const separatrices2_cells_isOnBoundary) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_ = separatrices2_cells;
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
    int getDescendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const;

    /**
     * Compute the geometrical embedding of the 1-separatrices. This
     * function needs the following internal pointers to be set:
     * outputSeparatrices1_numberOfPoints_
     * outputSeparatrices1_points_
     * outputSeparatrices1_numberOfCells_
     * outputSeparatrices1_cells_
     * inputScalarField_
     */
    template <typename dataType>
    int setSeparatrices1(
      const std::vector<Separatrix> &separatrices,
      const std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const;

    /**
     * Compute the ascending manifold of the maxima.
     */
    int setAscendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                 std::vector<SimplexId> &maxSeeds,
                                 SimplexId *const morseSmaleManifold,
                                 SimplexId &numberOfMaxima) const;

    /**
     * Compute the descending manifold of the minima.
     */
    int setDescendingSegmentation(const std::vector<dcg::Cell> &criticalPoints,
                                  SimplexId *const morseSmaleManifold,
                                  SimplexId &numberOfMinima) const;

    /**
     * Compute the final combinatorial Morse-Smale complex
     * segmentation.
     */
    int setFinalSegmentation(const SimplexId numberOfMaxima,
                             const SimplexId numberOfMinima,
                             const SimplexId *const ascendingManifold,
                             const SimplexId *const descendingManifold,
                             SimplexId *const morseSmaleManifold) const;

  protected:
    bool ReverveSaddleMaximumConnection;
    bool ReverveSaddleSaddleConnection;
    bool ComputeAscendingSeparatrices1;
    bool ComputeDescendingSeparatrices1;
    bool ComputeSaddleConnectors;
    bool ComputeAscendingSeparatrices2;
    bool ComputeDescendingSeparatrices2;
    bool ReturnSaddleConnectors;
    double SaddleConnectorsPersistenceThreshold;
    bool PrioritizeSpeedOverMemory;

    dcg::DiscreteGradient discreteGradient_;

    const void *inputScalarField_;
    Triangulation *inputTriangulation_;
    const void *inputOffsets_;

    SimplexId *outputCriticalPoints_numberOfPoints_;
    std::vector<float> *outputCriticalPoints_points_;
    std::vector<char> *outputCriticalPoints_points_cellDimensions_;
    std::vector<SimplexId> *outputCriticalPoints_points_cellIds_;
    void *outputCriticalPoints_points_cellScalars_;
    std::vector<char> *outputCriticalPoints_points_isOnBoundary_;
    std::vector<SimplexId> *outputCriticalPoints_points_PLVertexIdentifiers_;
    std::vector<SimplexId> *outputCriticalPoints_points_manifoldSize_;

    SimplexId *outputSeparatrices1_numberOfPoints_;
    std::vector<float> *outputSeparatrices1_points_;
    std::vector<char> *outputSeparatrices1_points_smoothingMask_;
    std::vector<char> *outputSeparatrices1_points_cellDimensions_;
    std::vector<SimplexId> *outputSeparatrices1_points_cellIds_;
    SimplexId *outputSeparatrices1_numberOfCells_;
    std::vector<SimplexId> *outputSeparatrices1_cells_;
    std::vector<SimplexId> *outputSeparatrices1_cells_sourceIds_;
    std::vector<SimplexId> *outputSeparatrices1_cells_destinationIds_;
    std::vector<SimplexId> *outputSeparatrices1_cells_separatrixIds_;
    std::vector<char> *outputSeparatrices1_cells_separatrixTypes_;
    void *outputSeparatrices1_cells_separatrixFunctionMaxima_;
    void *outputSeparatrices1_cells_separatrixFunctionMinima_;
    void *outputSeparatrices1_cells_separatrixFunctionDiffs_;
    std::vector<char> *outputSeparatrices1_cells_isOnBoundary_;

    SimplexId *outputSeparatrices2_numberOfPoints_;
    std::vector<float> *outputSeparatrices2_points_;
    SimplexId *outputSeparatrices2_numberOfCells_;
    std::vector<SimplexId> *outputSeparatrices2_cells_;
    std::vector<SimplexId> *outputSeparatrices2_cells_sourceIds_;
    std::vector<SimplexId> *outputSeparatrices2_cells_separatrixIds_;
    std::vector<char> *outputSeparatrices2_cells_separatrixTypes_;
    void *outputSeparatrices2_cells_separatrixFunctionMaxima_;
    void *outputSeparatrices2_cells_separatrixFunctionMinima_;
    void *outputSeparatrices2_cells_separatrixFunctionDiffs_;
    std::vector<char> *outputSeparatrices2_cells_isOnBoundary_;

    void *outputAscendingManifold_;
    void *outputDescendingManifold_;
    void *outputMorseSmaleManifold_;
  };
} // namespace ttk

template <typename dataType>
int ttk::AbstractMorseSmaleComplex::setSeparatrices1(
  const std::vector<Separatrix> &separatrices,
  const std::vector<std::vector<dcg::Cell>> &separatricesGeometry) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices1_numberOfPoints_) {
    std::cerr << "[AbstractMorseSmaleComplex] 1-separatrices pointer to "
                 "numberOfPoints is null."
              << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_points_) {
    std::cerr
      << "[AbstractMorseSmaleComplex] 1-separatrices pointer to points is null."
      << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_numberOfCells_) {
    std::cerr << "[AbstractMorseSmaleComplex] 1-separatrices pointer to "
                 "numberOfCells is null."
              << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_cells_) {
    std::cerr
      << "[AbstractMorseSmaleComplex] 1-separatrices pointer to cells is null."
      << std::endl;
    return -1;
  }
  if(!inputScalarField_) {
    std::cerr << "[AbstractMorseSmaleComplex] 1-separatrices pointer to the "
                 "input scalar field is null."
              << std::endl;
    return -1;
  }
#endif
  const dataType *const scalars
    = static_cast<const dataType *>(inputScalarField_);
  std::vector<dataType> *outputSeparatrices1_cells_separatrixFunctionMaxima
    = static_cast<std::vector<dataType> *>(
      outputSeparatrices1_cells_separatrixFunctionMaxima_);
  std::vector<dataType> *outputSeparatrices1_cells_separatrixFunctionMinima
    = static_cast<std::vector<dataType> *>(
      outputSeparatrices1_cells_separatrixFunctionMinima_);
  std::vector<dataType> *outputSeparatrices1_cells_separatrixFunctionDiffs
    = static_cast<std::vector<dataType> *>(
      outputSeparatrices1_cells_separatrixFunctionDiffs_);

  SimplexId pointId = (*outputSeparatrices1_numberOfPoints_);
  SimplexId cellId = (*outputSeparatrices1_numberOfCells_);
  SimplexId separatrixId = 0;
  if(outputSeparatrices1_cells_separatrixIds_
     and outputSeparatrices1_cells_separatrixIds_->size()) {
    separatrixId
      = *std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
                          outputSeparatrices1_cells_separatrixIds_->end())
        + 1;
  }

  const int dimensionality = inputTriangulation_->getCellVertexNumber(0) - 1;

  for(const Separatrix &separatrix : separatrices) {
    if(!separatrix.isValid_)
      continue;
    if(!separatrix.geometry_.size())
      continue;

    const dcg::Cell &saddle = separatrix.source_;
    const dcg::Cell &extremum = separatrix.destination_;

    // get separatrix type
    const char separatrixType = std::min(extremum.dim_, dimensionality - 1);

    // compute separatrix function diff
    const dataType separatrixFunctionMaximum
      = std::max(discreteGradient_.scalarMax<dataType>(saddle, scalars),
                 discreteGradient_.scalarMax<dataType>(extremum, scalars));
    const dataType separatrixFunctionMinimum
      = std::min(discreteGradient_.scalarMin<dataType>(saddle, scalars),
                 discreteGradient_.scalarMin<dataType>(extremum, scalars));
    const dataType separatrixFunctionDiff
      = separatrixFunctionMaximum - separatrixFunctionMinimum;

    // get boundary condition
    const char isOnBoundary = (char)discreteGradient_.isBoundary(saddle)
                              + (char)discreteGradient_.isBoundary(extremum);

    bool isFirst = true;
    for(const SimplexId geometryId : separatrix.geometry_) {
      SimplexId oldPointId = -1;
      for(auto cellIte = separatricesGeometry[geometryId].begin();
          cellIte != separatricesGeometry[geometryId].end(); ++cellIte) {
        const dcg::Cell &cell = *cellIte;

        float point[3];
        discreteGradient_.getCellIncenter(cell, point);

        outputSeparatrices1_points_->push_back(point[0]);
        outputSeparatrices1_points_->push_back(point[1]);
        outputSeparatrices1_points_->push_back(point[2]);

        if(outputSeparatrices1_points_smoothingMask_) {
          if(cellIte == separatricesGeometry[geometryId].begin()
             or cellIte == separatricesGeometry[geometryId].end() - 1)
            outputSeparatrices1_points_smoothingMask_->push_back(0);
          else
            outputSeparatrices1_points_smoothingMask_->push_back(1);
        }

        if(outputSeparatrices1_points_cellDimensions_)
          outputSeparatrices1_points_cellDimensions_->push_back(cell.dim_);

        if(outputSeparatrices1_points_cellIds_)
          outputSeparatrices1_points_cellIds_->push_back(cell.id_);

        if(oldPointId != -1) {
          outputSeparatrices1_cells_->push_back(2);
          outputSeparatrices1_cells_->push_back(oldPointId);
          outputSeparatrices1_cells_->push_back(pointId);

          if(outputSeparatrices1_cells_sourceIds_)
            outputSeparatrices1_cells_sourceIds_->push_back(saddle.id_);
          if(outputSeparatrices1_cells_destinationIds_)
            outputSeparatrices1_cells_destinationIds_->push_back(extremum.id_);
          if(outputSeparatrices1_cells_separatrixIds_)
            outputSeparatrices1_cells_separatrixIds_->push_back(separatrixId);
          if(outputSeparatrices1_cells_separatrixTypes_)
            outputSeparatrices1_cells_separatrixTypes_->push_back(
              separatrixType);
          if(outputSeparatrices1_cells_separatrixFunctionMaxima)
            outputSeparatrices1_cells_separatrixFunctionMaxima->push_back(
              separatrixFunctionMaximum);
          if(outputSeparatrices1_cells_separatrixFunctionMinima)
            outputSeparatrices1_cells_separatrixFunctionMinima->push_back(
              separatrixFunctionMinimum);
          if(outputSeparatrices1_cells_separatrixFunctionDiffs)
            outputSeparatrices1_cells_separatrixFunctionDiffs->push_back(
              separatrixFunctionDiff);
          if(outputSeparatrices1_cells_isOnBoundary_)
            outputSeparatrices1_cells_isOnBoundary_->push_back(isOnBoundary);

          ++cellId;
          isFirst = false;
        }

        oldPointId = pointId;
        ++pointId;
      }
    }

    if(!isFirst)
      ++separatrixId;
  }

  (*outputSeparatrices1_numberOfPoints_) = pointId;
  (*outputSeparatrices1_numberOfCells_) = cellId;

  return 0;
}

#endif // ABSTRACTMORSESMALECOMPLEX_H
