#include<AbstractMorseSmaleComplex.h>

AbstractMorseSmaleComplex::AbstractMorseSmaleComplex():
  ReverveSaddleMaximumConnection{},
  ReverveSaddleSaddleConnection{},
  ComputeAscendingSeparatrices1{},
  ComputeDescendingSeparatrices1{},
  ComputeSaddleConnectors{},
  ComputeAscendingSeparatrices2{},
  ComputeDescendingSeparatrices2{},
  ComputeAscendingSegmentation{},
  ComputeDescendingSegmentation{},
  ComputeFinalSegmentation{},

  inputScalarField_{},
  inputTriangulation_{},
  inputOffsets_{},

  outputCriticalPoints_numberOfPoints_{},
  outputCriticalPoints_points_{},
  outputCriticalPoints_points_cellDimensions_{},
  outputCriticalPoints_points_cellIds_{},
  outputCriticalPoints_points_cellScalars_{},
  outputCriticalPoints_points_isOnBoundary_{},
  outputCriticalPoints_points_PLVertexIdentifiers_{},
  outputCriticalPoints_points_manifoldSize_{},

  outputSeparatrices1_numberOfPoints_{},
  outputSeparatrices1_points_{},
  outputSeparatrices1_points_smoothingMask_{},
  outputSeparatrices1_points_cellDimensions_{},
  outputSeparatrices1_points_cellIds_{},
  outputSeparatrices1_numberOfCells_{},
  outputSeparatrices1_cells_{},
  outputSeparatrices1_cells_sourceIds_{},
  outputSeparatrices1_cells_destinationIds_{},
  outputSeparatrices1_cells_separatrixIds_{},
  outputSeparatrices1_cells_separatrixTypes_{},
  outputSeparatrices1_cells_separatrixFunctionMaxima_{},
  outputSeparatrices1_cells_separatrixFunctionMinima_{},
  outputSeparatrices1_cells_separatrixFunctionDiffs_{},
  outputSeparatrices1_cells_isOnBoundary_{},

  outputSeparatrices2_numberOfPoints_{},
  outputSeparatrices2_points_{},
  outputSeparatrices2_numberOfCells_{},
  outputSeparatrices2_cells_{},
  outputSeparatrices2_cells_sourceIds_{},
  outputSeparatrices2_cells_separatrixIds_{},
  outputSeparatrices2_cells_separatrixTypes_{},
  outputSeparatrices2_cells_separatrixFunctionMaxima_{},
  outputSeparatrices2_cells_separatrixFunctionMinima_{},
  outputSeparatrices2_cells_separatrixFunctionDiffs_{},
  outputSeparatrices2_cells_isOnBoundary_{},

  outputAscendingManifold_{},
  outputDescendingManifold_{},
  outputMorseSmaleManifold_{}
{
  discreteGradient_.setReverseSaddleMaximumConnection(true);
  discreteGradient_.setReverseSaddleSaddleConnection(true);
  ComputeAscendingSeparatrices1 = true;
  ComputeDescendingSeparatrices1 = true;
  ComputeAscendingSegmentation = true;
  ComputeDescendingSegmentation = true;
  ComputeFinalSegmentation = true;
}

AbstractMorseSmaleComplex::~AbstractMorseSmaleComplex(){
}

