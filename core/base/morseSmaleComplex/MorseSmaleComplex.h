/// \ingroup base
/// \class ttk::MorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK processing package for the computation of Morse-Smale complexes.
///
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \sa ttk::Triangulation
/// \sa vtkMorseSmaleComplex.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX_H
#define _MORSESMALECOMPLEX_H

// base code includes
#include<AbstractMorseSmaleComplex.h>
#include<MorseSmaleComplex2D.h>
#include<MorseSmaleComplex3D.h>

namespace ttk{

  class MorseSmaleComplex : public AbstractMorseSmaleComplex{

    public:

      MorseSmaleComplex();
      ~MorseSmaleComplex();

      int setIterationThreshold(const int iterationThreshold){
        return abstractMorseSmaleComplex_->setIterationThreshold(iterationThreshold);
      }

      int setReverseSaddleMaximumConnection(const bool state){
        return abstractMorseSmaleComplex_->setReverveSaddleMaximumConnection(state);
      }

      int setReverseSaddleSaddleConnection(const bool state){
        return abstractMorseSmaleComplex_->setReverveSaddleSaddleConnection(state);
      }

      int setComputeAscendingSeparatrices1(const bool state){
        return abstractMorseSmaleComplex_->setComputeAscendingSeparatrices1(state);
      }

      int setComputeDescendingSeparatrices1(const bool state){
        return abstractMorseSmaleComplex_->setComputeDescendingSeparatrices1(state);
      }

      int setComputeSaddleConnectors(const bool state){
        return abstractMorseSmaleComplex_->setComputeSaddleConnectors(state);
      }

      int setComputeAscendingSeparatrices2(const bool state){
        return abstractMorseSmaleComplex_->setComputeAscendingSeparatrices2(state);
      }

      int setComputeDescendingSeparatrices2(const bool state){
        return abstractMorseSmaleComplex_->setComputeDescendingSeparatrices2(state);
      }

      int setComputeAscendingSegmentation(const bool state){
        return abstractMorseSmaleComplex_->setComputeAscendingSegmentation(state);
      }

      int setComputeDescendingSegmentation(const bool state){
        return abstractMorseSmaleComplex_->setComputeDescendingSegmentation(state);
      }

      int setComputeFinalSegmentation(const bool state){
        return abstractMorseSmaleComplex_->setComputeFinalSegmentation(state);
      }

      int setupTriangulation(Triangulation* const data){
        inputTriangulation_=data;
        const int dimensionality=inputTriangulation_->getCellVertexNumber(0)-1;

        switch(dimensionality){
          case 2:
            abstractMorseSmaleComplex_=&morseSmaleComplex2D_;
            break;

          case 3:
            abstractMorseSmaleComplex_=&morseSmaleComplex3D_;
            break;
        }

        abstractMorseSmaleComplex_->setupTriangulation(inputTriangulation_);

        return 0;
      }

      inline int setDebugLevel(const int& debugLevel){
        morseSmaleComplex2D_.setDebugLevel(debugLevel);
        morseSmaleComplex3D_.setDebugLevel(debugLevel);
        debugLevel_=debugLevel;

        return 0;
      }

      inline int setThreadNumber(const int& threadNumber){
        morseSmaleComplex2D_.setThreadNumber(threadNumber);
        morseSmaleComplex3D_.setThreadNumber(threadNumber);
        threadNumber_=threadNumber;

        return 0;
      }

      inline int setWrapper(const Wrapper* const wrapper){
        morseSmaleComplex2D_.setWrapper(wrapper);
        morseSmaleComplex3D_.setWrapper(wrapper);

        return 0;
      }

      inline int setInputScalarField(void* const data){
        morseSmaleComplex2D_.setInputScalarField(data);
        morseSmaleComplex3D_.setInputScalarField(data);
        return 0;
      }

      inline int setInputOffsets(void* const data){
        morseSmaleComplex2D_.setInputOffsets(data);
        morseSmaleComplex3D_.setInputOffsets(data);

        return 0;
      }

      inline int setOutputCriticalPoints(int* const criticalPoints_numberOfPoints,
          vector<float>* const criticalPoints_points,
          vector<int>* const criticalPoints_points_cellDimensons,
          vector<int>* const criticalPoints_points_cellIds,
          void* const criticalPoints_points_cellScalars,
          vector<char>* const criticalPoints_points_isOnBoundary,
          vector<int>* const criticalPoints_points_PLVertexIdentifiers,
          vector<int>* const criticalPoints_points_manifoldSize){
        morseSmaleComplex2D_.setOutputCriticalPoints(criticalPoints_numberOfPoints,
            criticalPoints_points,
            criticalPoints_points_cellDimensons,
            criticalPoints_points_cellIds,
            criticalPoints_points_cellScalars,
            criticalPoints_points_isOnBoundary,
            criticalPoints_points_PLVertexIdentifiers,
            criticalPoints_points_manifoldSize);
        morseSmaleComplex3D_.setOutputCriticalPoints(criticalPoints_numberOfPoints,
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
        morseSmaleComplex2D_.setOutputSeparatrices1(separatrices1_numberOfPoints,
            separatrices1_points,
            separatrices1_points_smoothingMask,
            separatrices1_points_cellDimensions,
            separatrices1_points_cellIds,
            separatrices1_numberOfCells,
            separatrices1_cells,
            separatrices1_cells_sourceIds,
            separatrices1_cells_destinationIds,
            separatrices1_cells_separatrixIds,
            separatrices1_cells_separatrixTypes,
            separatrices1_cells_separatrixFunctionMaxima,
            separatrices1_cells_separatrixFunctionMinima,
            separatrices1_cells_separatrixFunctionDiffs,
            separatrices1_cells_isOnBoundary);
        morseSmaleComplex3D_.setOutputSeparatrices1(separatrices1_numberOfPoints,
            separatrices1_points,
            separatrices1_points_smoothingMask,
            separatrices1_points_cellDimensions,
            separatrices1_points_cellIds,
            separatrices1_numberOfCells,
            separatrices1_cells,
            separatrices1_cells_sourceIds,
            separatrices1_cells_destinationIds,
            separatrices1_cells_separatrixIds,
            separatrices1_cells_separatrixTypes,
            separatrices1_cells_separatrixFunctionMaxima,
            separatrices1_cells_separatrixFunctionMinima,
            separatrices1_cells_separatrixFunctionDiffs,
            separatrices1_cells_isOnBoundary);

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
        morseSmaleComplex3D_.setOutputSeparatrices2(separatrices2_numberOfPoints,
            separatrices2_points,
            separatrices2_numberOfCells,
            separatrices2_cells,
            separatrices2_cells_sourceIds,
            separatrices2_cells_separatrixIds,
            separatrices2_cells_separatrixTypes,
            separatrices2_cells_separatrixFunctionMaxima,
            separatrices2_cells_separatrixFunctionMinima,
            separatrices2_cells_separatrixFunctionDiffs,
            separatrices2_cells_isOnBoundary);

        return 0;
      }

      inline int setOutputMorseComplexes(void* const ascendingManifold,
          void* const descendingManifold,
          void* const morseSmaleManifold){
        morseSmaleComplex2D_.setOutputMorseComplexes(ascendingManifold,
            descendingManifold,
            morseSmaleManifold);
        morseSmaleComplex3D_.setOutputMorseComplexes(ascendingManifold,
            descendingManifold,
            morseSmaleManifold);

        return 0;
      }

      template<typename dataType>
        int execute(){
          const int dimensionality=inputTriangulation_->getCellVertexNumber(0)-1;

          switch(dimensionality){
            case 2:
              morseSmaleComplex2D_.execute<dataType>();
              break;

            case 3:
              morseSmaleComplex3D_.execute<dataType>();
              break;
          }

          return 0;
        }

    protected:

      AbstractMorseSmaleComplex* abstractMorseSmaleComplex_;
      MorseSmaleComplex2D morseSmaleComplex2D_;
      MorseSmaleComplex3D morseSmaleComplex3D_;

  };
}

#endif // MORSESMALECOMPLEX_H
