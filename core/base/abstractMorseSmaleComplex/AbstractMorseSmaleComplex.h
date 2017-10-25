/// \ingroup base
/// \class ttk::AbstractMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %abstractMorseSmaleComplex processing package.
///
/// %AbstractMorseSmaleComplex is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkAbstractMorseSmaleComplex.cpp %for a usage example.

#ifndef _ABSTRACTMORSESMALECOMPLEX_H
#define _ABSTRACTMORSESMALECOMPLEX_H

// base code includes
#include<Wrapper.h>
#include<DiscreteGradient.h>
#include<Triangulation.h>

#include<queue>

namespace ttk{

  struct Separatrix{
    // default :
    Separatrix():
      isValid_{},
      source_{},
      destination_{},
      isReversed_{},
      geometry_{}
    {}

    // for just one segment :
    Separatrix(const bool isValid,
        const Cell& saddle,
        const Cell& extremum,
        const bool isSegmentReversed,
        const int segmentGeometry):
      isValid_{isValid},
      source_{saddle},
      destination_{extremum}
    {
      isReversed_.push_back(isSegmentReversed);
      geometry_.push_back(segmentGeometry);
    }

    // for multiple segments :
    Separatrix(const bool isValid,
        const Cell& saddle,
        const Cell& extremum,
        const vector<char>& isReversed,
        const vector<int>& geometry):
      isValid_{isValid},
      source_{saddle},
      destination_{extremum},
      isReversed_{isReversed},
      geometry_{geometry}
    {}

    Separatrix(const Separatrix& separatrix):
      isValid_{separatrix.isValid_},
      source_{separatrix.source_},
      destination_{separatrix.destination_},
      isReversed_{separatrix.isReversed_},
      geometry_{separatrix.geometry_}
    {}

    Separatrix(Separatrix&& separatrix):
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

    bool isValid_;
    Cell source_;
    Cell destination_;
    vector<char> isReversed_;
    vector<int> geometry_;
  };

  class AbstractMorseSmaleComplex : public Debug{

    public:

      AbstractMorseSmaleComplex();
      ~AbstractMorseSmaleComplex();

      int setIterationThreshold(const int iterationThreshold){
        discreteGradient_.setIterationThreshold(iterationThreshold);
        return 0;
      }

      int setReverveSaddleMaximumConnection(const bool state){
        discreteGradient_.setReverseSaddleMaximumConnection(state);
        return 0;
      }

      int setReverveSaddleSaddleConnection(const bool state){
        discreteGradient_.setReverseSaddleSaddleConnection(state);
        return 0;
      }

      int setComputeAscendingSeparatrices1(const bool state){
        ComputeAscendingSeparatrices1=state;
        return 0;
      }

      int setComputeDescendingSeparatrices1(const bool state){
        ComputeDescendingSeparatrices1=state;
        return 0;
      }

      int setComputeSaddleConnectors(const bool state){
        ComputeSaddleConnectors=state;
        return 0;
      }

      int setComputeAscendingSeparatrices2(const bool state){
        ComputeAscendingSeparatrices2=state;
        return 0;
      }

      int setComputeDescendingSeparatrices2(const bool state){
        ComputeDescendingSeparatrices2=state;
        return 0;
      }

      int setComputeAscendingSegmentation(const bool state){
        ComputeAscendingSegmentation=state;
        return 0;
      }

      int setComputeDescendingSegmentation(const bool state){
        ComputeDescendingSegmentation=state;
        return 0;
      }

      int setComputeFinalSegmentation(const bool state){
        ComputeFinalSegmentation=state;
        return 0;
      }

      int setReturnSaddleConnectors(const bool state){
        ReturnSaddleConnectors=state;
        return 0;
      }

      int setSaddleConnectorsPersistenceThreshold(const double threshold){
        SaddleConnectorsPersistenceThreshold=threshold;
        return 0;
      }

      inline int setupTriangulation(Triangulation* const data){
        inputTriangulation_=data;
        discreteGradient_.setupTriangulation(inputTriangulation_);

        inputTriangulation_->preprocessCellEdges();
        inputTriangulation_->preprocessCellNeighbors();
        return 0;
      }

      inline int setInputScalarField(void* const data){
        inputScalarField_=data;
        discreteGradient_.setInputScalarField(inputScalarField_);
        return 0;
      }

      inline int setInputOffsets(void* const data){
        inputOffsets_=data;
        discreteGradient_.setInputOffsets(inputOffsets_);
        return 0;
      }

      inline int setOutputCriticalPoints(int* const criticalPoints_numberOfPoints,
          vector<float>* const criticalPoints_points,
          vector<int>* const criticalPoints_points_cellDimensons,
          vector<int>* const criticalPoints_points_cellIds,
          void* const criticalPoints_points_cellScalars,
          vector<char>* const criticalPoints_points_isOnBoundary,
          vector<int>* const criticalPoints_points_PLVertexIdentifiers,
          vector<int>* criticalPoints_points_manifoldSize){
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

      inline int setOutputSeparatrices1(int* const separatrices1_numberOfPoints,
          vector<float>* const separatrices1_points,
          vector<char>* const separatrices1_points_smoothingMask,
          vector<int>* const separatrices1_points_cellDimensions,
          vector<int>* const separatrices1_points_cellIds,
          int* const separatrices1_numberOfCells,
          vector<int>* const separatrices1_cells,
          vector<int>* const separatrices1_cells_sourceIds,
          vector<int>* const separatrices1_cells_destinationIds,
          vector<int>* const separatrices1_cells_separatrixIds,
          vector<char>* const separatrices1_cells_separatrixTypes,
          void* const separatrices1_cells_separatrixFunctionMaxima,
          void* const separatrices1_cells_separatrixFunctionMinima,
          void* const separatrices1_cells_separatrixFunctionDiffs,
          vector<char>* const separatrices1_cells_isOnBoundary){
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

      inline int setOutputSeparatrices2(int* const separatrices2_numberOfPoints,
          vector<float>* const separatrices2_points,
          int* const separatrices2_numberOfCells,
          vector<int>* const separatrices2_cells,
          vector<int>* const separatrices2_cells_sourceIds,
          vector<int>* const separatrices2_cells_separatrixIds,
          vector<char>* const separatrices2_cells_separatrixTypes,
          void* const separatrices2_cells_separatrixFunctionMaxima,
          void* const separatrices2_cells_separatrixFunctionMinima,
          void* const separatrices2_cells_separatrixFunctionDiffs,
          vector<char>* const separatrices2_cells_isOnBoundary){
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

      inline int setOutputMorseComplexes(void* const ascendingManifold,
          void* const descendingManifold,
          void* const morseSmaleManifold){
        outputAscendingManifold_=ascendingManifold;
        outputDescendingManifold_=descendingManifold;
        outputMorseSmaleManifold_=morseSmaleManifold;
        return 0;
      }

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

      DiscreteGradient discreteGradient_;

      void* inputScalarField_;
      Triangulation* inputTriangulation_;
      void* inputOffsets_;

      int* outputCriticalPoints_numberOfPoints_;
      vector<float>* outputCriticalPoints_points_;
      vector<int>* outputCriticalPoints_points_cellDimensions_;
      vector<int>* outputCriticalPoints_points_cellIds_;
      void* outputCriticalPoints_points_cellScalars_;
      vector<char>* outputCriticalPoints_points_isOnBoundary_;
      vector<int>* outputCriticalPoints_points_PLVertexIdentifiers_;
      vector<int>* outputCriticalPoints_points_manifoldSize_;

      int* outputSeparatrices1_numberOfPoints_;
      vector<float>* outputSeparatrices1_points_;
      vector<char>* outputSeparatrices1_points_smoothingMask_;
      vector<int>* outputSeparatrices1_points_cellDimensions_;
      vector<int>* outputSeparatrices1_points_cellIds_;
      int* outputSeparatrices1_numberOfCells_;
      vector<int>* outputSeparatrices1_cells_;
      vector<int>* outputSeparatrices1_cells_sourceIds_;
      vector<int>* outputSeparatrices1_cells_destinationIds_;
      vector<int>* outputSeparatrices1_cells_separatrixIds_;
      vector<char>* outputSeparatrices1_cells_separatrixTypes_;
      void* outputSeparatrices1_cells_separatrixFunctionMaxima_;
      void* outputSeparatrices1_cells_separatrixFunctionMinima_;
      void* outputSeparatrices1_cells_separatrixFunctionDiffs_;
      vector<char>* outputSeparatrices1_cells_isOnBoundary_;

      int* outputSeparatrices2_numberOfPoints_;
      vector<float>* outputSeparatrices2_points_;
      int* outputSeparatrices2_numberOfCells_;
      vector<int>* outputSeparatrices2_cells_;
      vector<int>* outputSeparatrices2_cells_sourceIds_;
      vector<int>* outputSeparatrices2_cells_separatrixIds_;
      vector<char>* outputSeparatrices2_cells_separatrixTypes_;
      void* outputSeparatrices2_cells_separatrixFunctionMaxima_;
      void* outputSeparatrices2_cells_separatrixFunctionMinima_;
      void* outputSeparatrices2_cells_separatrixFunctionDiffs_;
      vector<char>* outputSeparatrices2_cells_isOnBoundary_;

      void* outputAscendingManifold_;
      void* outputDescendingManifold_;
      void* outputMorseSmaleManifold_;
  };
}

#endif // ABSTRACTMORSESMALECOMPLEX_H
