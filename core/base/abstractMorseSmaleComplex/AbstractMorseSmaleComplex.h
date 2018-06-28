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
#include<Wrapper.h>
#include<DiscreteGradient.h>
#include<Triangulation.h>

#include<queue>

namespace ttk{

  /**
   * Utility class representing Ridge lines, Valley lines
   * and Saddle connectors.
   */
  struct Separatrix{
    // default :
    explicit Separatrix():
      isValid_{},
      source_{},
      destination_{},
      isReversed_{},
      geometry_{}
    {}

    // initialization with one segment :
    explicit Separatrix(const bool isValid,
        const dcg::Cell& saddle,
        const dcg::Cell& extremum,
        const bool isSegmentReversed,
        const dcg::simplexId_t segmentGeometry):
      isValid_{isValid},
      source_{saddle},
      destination_{extremum}
    {
      isReversed_.push_back(isSegmentReversed);
      geometry_.push_back(segmentGeometry);
    }

    // initialization with multiple segments :
    explicit Separatrix(const bool isValid,
        const dcg::Cell& saddle,
        const dcg::Cell& extremum,
        const std::vector<char>& isReversed,
        const std::vector<dcg::simplexId_t>& geometry):
      isValid_{isValid},
      source_{saddle},
      destination_{extremum},
      isReversed_{isReversed},
      geometry_{geometry}
    {}

    explicit Separatrix(const Separatrix& separatrix):
      isValid_{separatrix.isValid_},
      source_{separatrix.source_},
      destination_{separatrix.destination_},
      isReversed_{separatrix.isReversed_},
      geometry_{separatrix.geometry_}
    {}

    explicit Separatrix(Separatrix&& separatrix):
      isValid_{separatrix.isValid_},
      source_{std::move(separatrix.source_)},
      destination_{std::move(separatrix.destination_)},
      isReversed_{std::move(separatrix.isReversed_)},
      geometry_{std::move(separatrix.geometry_)}
    {}

    Separatrix& operator=(Separatrix&& separatrix){
      isValid_=separatrix.isValid_;
      source_=std::move(separatrix.source_);
      destination_=std::move(separatrix.destination_);
      isReversed_=std::move(separatrix.isReversed_);
      geometry_=std::move(separatrix.geometry_);

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
    std::vector<dcg::simplexId_t> geometry_;
  };

  /**
   * Parent class containing convenience functions shared between
   * Morse-Smale Complex algorithms for 2D and 3D domains.
   */
  class AbstractMorseSmaleComplex : public Debug{

    public:

      AbstractMorseSmaleComplex();
      ~AbstractMorseSmaleComplex();

      /**
       * Set the threshold for the iterative gradient reversal process.
       * Disable thresholding with -1 (default).
       */
      int setIterationThreshold(const int iterationThreshold){
        discreteGradient_.setIterationThreshold(iterationThreshold);
        return 0;
      }

      /**
       * Enable/Disable gradient reversal of (saddle,...,maximum) vpaths
       * (enabled by default).
       */
      int setReverveSaddleMaximumConnection(const bool state){
        discreteGradient_.setReverseSaddleMaximumConnection(state);
        return 0;
      }

      /**
       * Enable/Disable gradient reversal of (saddle,...,saddle) vpaths
       * (enabled by default).
       */
      int setReverveSaddleSaddleConnection(const bool state){
        discreteGradient_.setReverseSaddleSaddleConnection(state);
        return 0;
      }

      /**
       * Enable/Disable computation of the geometrical embedding of
       * the ascending manifolds of the critical points
       * (enabled by default).
       */
      int setComputeAscendingSeparatrices1(const bool state){
        ComputeAscendingSeparatrices1=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the geometrical embedding of
       * the descending manifolds of the critical points
       * (enabled by default).
       */
      int setComputeDescendingSeparatrices1(const bool state){
        ComputeDescendingSeparatrices1=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the geometrical embedding of
       * the visible saddle-connectors (enabled by default).
       */
      int setComputeSaddleConnectors(const bool state){
        ComputeSaddleConnectors=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the geometrical embedding of
       * the ascending 2-separatrices (disabled by default).
       */
      int setComputeAscendingSeparatrices2(const bool state){
        ComputeAscendingSeparatrices2=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the geometrical embedding of
       * the descending 2-separatrices (disabled by default).
       */
      int setComputeDescendingSeparatrices2(const bool state){
        ComputeDescendingSeparatrices2=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the ascending manifold
       * of the maxima (enabled by default).
       */
      int setComputeAscendingSegmentation(const bool state){
        ComputeAscendingSegmentation=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the descending manifold
       * of the minima (enabled by default).
       */
      int setComputeDescendingSegmentation(const bool state){
        ComputeDescendingSegmentation=state;
        return 0;
      }

      /**
       * Enable/Disable computation of the final combinatorial
       * Morse-Smale Complex (segmentation as a field on vertices,
       * enabled by default).
       */
      int setComputeFinalSegmentation(const bool state){
        ComputeFinalSegmentation=state;
        return 0;
      }

      /**
       * Enable/Disable post-processing gradient reversal of
       * the (saddle,...,saddle) vpaths under a given persistence
       * threshold (disabled by default).
       */
      int setReturnSaddleConnectors(const bool state){
        ReturnSaddleConnectors=state;
        discreteGradient_.setReturnSaddleConnectors(state);
        return 0;
      }

      /**
       * Set the threshold value for post-processing of
       * (saddle,...,saddle) vpaths gradient reversal
       * (default value is 0.0).
       */
      int setSaddleConnectorsPersistenceThreshold(const double threshold){
        SaddleConnectorsPersistenceThreshold=threshold;
        discreteGradient_.setSaddleConnectorsPersistenceThreshold(threshold);
        return 0;
      }

      /**
       * Enable/Disable compromise on memory allowing
       * computation the geometry of the 2-separatrices
       * in parallel.
       */
      int setPrioritizeSpeedOverMemory(const bool state){
        PrioritizeSpeedOverMemory=state;
        return 0;
      }

      /**
       * Set the input triangulation and preprocess the needed
       * mesh traversal queries.
       */
      inline int setupTriangulation(Triangulation* const data){
        inputTriangulation_=data;
        discreteGradient_.setupTriangulation(inputTriangulation_);

        inputTriangulation_->preprocessCellEdges();
        inputTriangulation_->preprocessCellNeighbors();
        return 0;
      }

      /**
       * Set the input scalar field associated on the points of the data set.
       */
      inline int setInputScalarField(void* const data){
        inputScalarField_=data;
        discreteGradient_.setInputScalarField(inputScalarField_);
        return 0;
      }

      /**
       * Set the input offset field associated on the points of the data set
       * (if none, identifiers are used instead).
       */
      inline int setInputOffsets(void* const data){
        inputOffsets_=data;
        discreteGradient_.setInputOffsets(inputOffsets_);
        return 0;
      }

      /**
       * Set the output critical points data pointers.
       */
      inline int setOutputCriticalPoints(dcg::simplexId_t* const criticalPoints_numberOfPoints,
          std::vector<float>* const criticalPoints_points,
          std::vector<char>* const criticalPoints_points_cellDimensons,
          std::vector<dcg::simplexId_t>* const criticalPoints_points_cellIds,
          void* const criticalPoints_points_cellScalars,
          std::vector<char>* const criticalPoints_points_isOnBoundary,
          std::vector<dcg::simplexId_t>* const criticalPoints_points_PLVertexIdentifiers,
          std::vector<dcg::simplexId_t>* criticalPoints_points_manifoldSize){
        discreteGradient_.setOutputCriticalPoints(criticalPoints_numberOfPoints,
            criticalPoints_points,
            criticalPoints_points_cellDimensons,
            criticalPoints_points_cellIds,
            criticalPoints_points_cellScalars,
            criticalPoints_points_isOnBoundary,
            criticalPoints_points_PLVertexIdentifiers,
            criticalPoints_points_manifoldSize);
        return 0;
      }

      /**
       * Set the data pointers to the output 1-separatrices.
       */
      inline int setOutputSeparatrices1(dcg::simplexId_t* const separatrices1_numberOfPoints,
          std::vector<float>* const separatrices1_points,
          std::vector<char>* const separatrices1_points_smoothingMask,
          std::vector<char>* const separatrices1_points_cellDimensions,
          std::vector<dcg::simplexId_t>* const separatrices1_points_cellIds,
          dcg::simplexId_t* const separatrices1_numberOfCells,
          std::vector<dcg::simplexId_t>* const separatrices1_cells,
          std::vector<dcg::simplexId_t>* const separatrices1_cells_sourceIds,
          std::vector<dcg::simplexId_t>* const separatrices1_cells_destinationIds,
          std::vector<dcg::simplexId_t>* const separatrices1_cells_separatrixIds,
          std::vector<char>* const separatrices1_cells_separatrixTypes,
          void* const separatrices1_cells_separatrixFunctionMaxima,
          void* const separatrices1_cells_separatrixFunctionMinima,
          void* const separatrices1_cells_separatrixFunctionDiffs,
          std::vector<char>* const separatrices1_cells_isOnBoundary){
        outputSeparatrices1_numberOfPoints_=separatrices1_numberOfPoints;
        outputSeparatrices1_points_=separatrices1_points;
        outputSeparatrices1_points_smoothingMask_=separatrices1_points_smoothingMask;
        outputSeparatrices1_points_cellDimensions_=separatrices1_points_cellDimensions;
        outputSeparatrices1_points_cellIds_=separatrices1_points_cellIds;
        outputSeparatrices1_numberOfCells_=separatrices1_numberOfCells;
        outputSeparatrices1_cells_=separatrices1_cells;
        outputSeparatrices1_cells_sourceIds_=separatrices1_cells_sourceIds;
        outputSeparatrices1_cells_destinationIds_=separatrices1_cells_destinationIds;
        outputSeparatrices1_cells_separatrixIds_=separatrices1_cells_separatrixIds;
        outputSeparatrices1_cells_separatrixTypes_=separatrices1_cells_separatrixTypes;
        outputSeparatrices1_cells_separatrixFunctionMaxima_=separatrices1_cells_separatrixFunctionMaxima;
        outputSeparatrices1_cells_separatrixFunctionMinima_=separatrices1_cells_separatrixFunctionMinima;
        outputSeparatrices1_cells_separatrixFunctionDiffs_=separatrices1_cells_separatrixFunctionDiffs;
        outputSeparatrices1_cells_isOnBoundary_=separatrices1_cells_isOnBoundary;
        return 0;
      }

      /**
       * Set the data pointers to the output 2-separatrices.
       */
      inline int setOutputSeparatrices2(dcg::simplexId_t* const separatrices2_numberOfPoints,
          std::vector<float>* const separatrices2_points,
          dcg::simplexId_t* const separatrices2_numberOfCells,
          std::vector<dcg::simplexId_t>* const separatrices2_cells,
          std::vector<dcg::simplexId_t>* const separatrices2_cells_sourceIds,
          std::vector<dcg::simplexId_t>* const separatrices2_cells_separatrixIds,
          std::vector<char>* const separatrices2_cells_separatrixTypes,
          void* const separatrices2_cells_separatrixFunctionMaxima,
          void* const separatrices2_cells_separatrixFunctionMinima,
          void* const separatrices2_cells_separatrixFunctionDiffs,
          std::vector<char>* const separatrices2_cells_isOnBoundary){
        outputSeparatrices2_numberOfPoints_=separatrices2_numberOfPoints;
        outputSeparatrices2_points_=separatrices2_points;
        outputSeparatrices2_numberOfCells_=separatrices2_numberOfCells;
        outputSeparatrices2_cells_=separatrices2_cells;
        outputSeparatrices2_cells_sourceIds_=separatrices2_cells_sourceIds;
        outputSeparatrices2_cells_separatrixIds_=separatrices2_cells_separatrixIds;
        outputSeparatrices2_cells_separatrixTypes_=separatrices2_cells_separatrixTypes;
        outputSeparatrices2_cells_separatrixFunctionMaxima_=separatrices2_cells_separatrixFunctionMaxima;
        outputSeparatrices2_cells_separatrixFunctionMinima_=separatrices2_cells_separatrixFunctionMinima;
        outputSeparatrices2_cells_separatrixFunctionDiffs_=separatrices2_cells_separatrixFunctionDiffs;
        outputSeparatrices2_cells_isOnBoundary_=separatrices2_cells_isOnBoundary;
        return 0;
      }

      /**
       * Set the data pointers to the output segmentation scalar fields.
       */
      inline int setOutputMorseComplexes(void* const ascendingManifold,
          void* const descendingManifold,
          void* const morseSmaleManifold){
        outputAscendingManifold_=ascendingManifold;
        outputDescendingManifold_=descendingManifold;
        outputMorseSmaleManifold_=morseSmaleManifold;
        return 0;
      }

      /**
       * Compute the descending 1-separatrices by reading into the discrete
       * gradient.
       */
      int getDescendingSeparatrices1(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const;

      /**
       * Compute the geometrical embedding of the 1-separatrices.
       */
      template<typename dataType>
      int setSeparatrices1(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const;

      /**
       * Compute the ascending manifold of the maxima.
       */
      int setAscendingSegmentation(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<dcg::simplexId_t>& maxSeeds,
          dcg::simplexId_t* const morseSmaleManifold,
          dcg::simplexId_t& numberOfMaxima) const;

      /**
       * Compute the descending manifold of the minima.
       */
      int setDescendingSegmentation(const std::vector<dcg::Cell>& criticalPoints,
          dcg::simplexId_t* const morseSmaleManifold,
          dcg::simplexId_t& numberOfMinima) const;

      /**
       * Compute the final combinatorial Morse-Smale complex
       * segmentation.
       */
      int setFinalSegmentation(const dcg::simplexId_t numberOfMaxima,
          const dcg::simplexId_t numberOfMinima,
          const dcg::simplexId_t* const ascendingManifold,
          const dcg::simplexId_t* const descendingManifold,
          dcg::simplexId_t* const morseSmaleManifold) const;

    protected:

      bool ReverveSaddleMaximumConnection;
      bool ReverveSaddleSaddleConnection;
      bool ComputeAscendingSeparatrices1;
      bool ComputeDescendingSeparatrices1;
      bool ComputeSaddleConnectors;
      bool ComputeAscendingSeparatrices2;
      bool ComputeDescendingSeparatrices2;
      bool ComputeAscendingSegmentation;
      bool ComputeDescendingSegmentation;
      bool ComputeFinalSegmentation;
      bool ReturnSaddleConnectors;
      double SaddleConnectorsPersistenceThreshold;
      bool PrioritizeSpeedOverMemory;

      dcg::DiscreteGradient discreteGradient_;

      void* inputScalarField_;
      Triangulation* inputTriangulation_;
      void* inputOffsets_;

      dcg::simplexId_t* outputCriticalPoints_numberOfPoints_;
      std::vector<float>* outputCriticalPoints_points_;
      std::vector<char>* outputCriticalPoints_points_cellDimensions_;
      std::vector<dcg::simplexId_t>* outputCriticalPoints_points_cellIds_;
      void* outputCriticalPoints_points_cellScalars_;
      std::vector<char>* outputCriticalPoints_points_isOnBoundary_;
      std::vector<dcg::simplexId_t>* outputCriticalPoints_points_PLVertexIdentifiers_;
      std::vector<dcg::simplexId_t>* outputCriticalPoints_points_manifoldSize_;

      dcg::simplexId_t* outputSeparatrices1_numberOfPoints_;
      std::vector<float>* outputSeparatrices1_points_;
      std::vector<char>* outputSeparatrices1_points_smoothingMask_;
      std::vector<char>* outputSeparatrices1_points_cellDimensions_;
      std::vector<dcg::simplexId_t>* outputSeparatrices1_points_cellIds_;
      dcg::simplexId_t* outputSeparatrices1_numberOfCells_;
      std::vector<dcg::simplexId_t>* outputSeparatrices1_cells_;
      std::vector<dcg::simplexId_t>* outputSeparatrices1_cells_sourceIds_;
      std::vector<dcg::simplexId_t>* outputSeparatrices1_cells_destinationIds_;
      std::vector<dcg::simplexId_t>* outputSeparatrices1_cells_separatrixIds_;
      std::vector<char>* outputSeparatrices1_cells_separatrixTypes_;
      void* outputSeparatrices1_cells_separatrixFunctionMaxima_;
      void* outputSeparatrices1_cells_separatrixFunctionMinima_;
      void* outputSeparatrices1_cells_separatrixFunctionDiffs_;
      std::vector<char>* outputSeparatrices1_cells_isOnBoundary_;

      dcg::simplexId_t* outputSeparatrices2_numberOfPoints_;
      std::vector<float>* outputSeparatrices2_points_;
      dcg::simplexId_t* outputSeparatrices2_numberOfCells_;
      std::vector<dcg::simplexId_t>* outputSeparatrices2_cells_;
      std::vector<dcg::simplexId_t>* outputSeparatrices2_cells_sourceIds_;
      std::vector<dcg::simplexId_t>* outputSeparatrices2_cells_separatrixIds_;
      std::vector<char>* outputSeparatrices2_cells_separatrixTypes_;
      void* outputSeparatrices2_cells_separatrixFunctionMaxima_;
      void* outputSeparatrices2_cells_separatrixFunctionMinima_;
      void* outputSeparatrices2_cells_separatrixFunctionDiffs_;
      std::vector<char>* outputSeparatrices2_cells_isOnBoundary_;

      void* outputAscendingManifold_;
      void* outputDescendingManifold_;
      void* outputMorseSmaleManifold_;
  };
}

template<typename dataType>
int ttk::AbstractMorseSmaleComplex::setSeparatrices1(const 
std::vector<Separatrix>& separatrices,
    const std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionDiffs_);

  dcg::simplexId_t pointId=(*outputSeparatrices1_numberOfPoints_);
  dcg::simplexId_t cellId=(*outputSeparatrices1_numberOfCells_);
  dcg::simplexId_t separatrixId=0;
  if(outputSeparatrices1_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
        outputSeparatrices1_cells_separatrixIds_->end())+1;
  }

  const int dimensionality=inputTriangulation_->getCellVertexNumber(0)-1;

  for(const Separatrix& separatrix : separatrices){
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle=separatrix.source_;
    const dcg::Cell& extremum=separatrix.destination_;

    // get separatrix type
    const char separatrixType=std::min(extremum.dim_, dimensionality-1);

    // compute separatrix function diff
    const dataType separatrixFunctionMaximum=std::max(discreteGradient_.scalarMax<dataType>(saddle, scalars),
        discreteGradient_.scalarMax<dataType>(extremum, scalars));
    const dataType separatrixFunctionMinimum=std::min(discreteGradient_.scalarMin<dataType>(saddle, scalars),
        discreteGradient_.scalarMin<dataType>(extremum, scalars));
    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    // get boundary condition
    const char isOnBoundary=(char)discreteGradient_.isBoundary(saddle) + (char)discreteGradient_.isBoundary(extremum);

    bool isFirst=true;
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      dcg::simplexId_t oldPointId=-1;
      for(auto cellIte=separatricesGeometry[geometryId].begin(); cellIte!=separatricesGeometry[geometryId].end(); ++cellIte){
        const dcg::Cell& cell=*cellIte;

        float point[3];
        discreteGradient_.getCellIncenter(cell, point);

        outputSeparatrices1_points_->push_back(point[0]);
        outputSeparatrices1_points_->push_back(point[1]);
        outputSeparatrices1_points_->push_back(point[2]);

        if(cellIte==separatricesGeometry[geometryId].begin() or
            cellIte==separatricesGeometry[geometryId].end()-1)
          outputSeparatrices1_points_smoothingMask_->push_back(0);
        else
          outputSeparatrices1_points_smoothingMask_->push_back(1);
        outputSeparatrices1_points_cellDimensions_->push_back(cell.dim_);
        outputSeparatrices1_points_cellIds_->push_back(cell.id_);

        if(oldPointId!=-1){
          outputSeparatrices1_cells_->push_back(2);
          outputSeparatrices1_cells_->push_back(oldPointId);
          outputSeparatrices1_cells_->push_back(pointId);

          outputSeparatrices1_cells_sourceIds_->push_back(saddle.id_);
          outputSeparatrices1_cells_destinationIds_->push_back(extremum.id_);
          outputSeparatrices1_cells_separatrixIds_->push_back(separatrixId);
          outputSeparatrices1_cells_separatrixTypes_->push_back(separatrixType);
          outputSeparatrices1_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
          outputSeparatrices1_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
          outputSeparatrices1_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
          outputSeparatrices1_cells_isOnBoundary_->push_back(isOnBoundary);

          ++cellId;
          isFirst=false;
        }

        oldPointId=pointId;
        ++pointId;
      }
    }

    if(!isFirst)
      ++separatrixId;
  }

  (*outputSeparatrices1_numberOfPoints_)=pointId;
  (*outputSeparatrices1_numberOfCells_)=cellId;

  return 0;
}

#endif // ABSTRACTMORSESMALECOMPLEX_H
