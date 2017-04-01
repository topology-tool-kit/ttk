/// \ingroup baseCode
/// \class ttk::MorseSmaleComplex2D
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex2D processing package.
///
/// %MorseSmaleComplex2D is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkMorseSmaleComplex2D.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX2D_H
#define _MORSESMALECOMPLEX2D_H

// base code includes
#include<AbstractMorseSmaleComplex.h>

namespace ttk{

  class MorseSmaleComplex2D : public AbstractMorseSmaleComplex{

    public:

      MorseSmaleComplex2D();
      ~MorseSmaleComplex2D();

      template<typename dataType>
        int execute();

      int getSeparatrices(const vector<Cell>& criticalPoints,
          vector<Separatrix>& separatrices,
          vector<vector<Cell>>& separatricesGeometry) const;

      template<typename dataType>
        int setSeparatrices(const vector<Separatrix>& separatrices,
            const vector<vector<Cell>>& separatricesGeometry) const;

      int setAscendingSegmentation(const vector<Cell>& criticalPoints,
          int* const morseSmaleManifold,
          int& numberOfMinima) const;

      int setDescendingSegmentation(const vector<Cell>& criticalPoints,
          int* const morseSmaleManifold,
          int& numberOfMinima) const;

      int setFinalSegmentation(const int numberOfMaxima,
          const int numberOfMinima,
          const int* const ascendingManifold,
          const int* const descendingManifold,
          int* const morseSmaleManifold) const;

  };
}

template<typename dataType>
int MorseSmaleComplex2D::setSeparatrices(const vector<Separatrix>& separatrices,
    const vector<vector<Cell>>& separatricesGeometry) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMaxima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMaxima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMinima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMinima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionDiffs=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionDiffs_);

  (*outputSeparatrices1_numberOfPoints_)=0;
  (*outputSeparatrices1_numberOfCells_)=0;

  const int dimensionality=inputTriangulation_->getCellVertexNumber(0)-1;

  int pointId{};
  int cellId{};
  int separatrixId{};
  for(const Separatrix& separatrix : separatrices){
    if(separatrix.isValid_){
      const Cell& saddle=separatrix.source_;
      const Cell& extremum=separatrix.destination_;

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

      for(const int geometryId : separatrix.geometry_){
        int oldPointId=-1;
        for(const Cell& cell : separatricesGeometry[geometryId]){
          float point[3];
          discreteGradient_.getCellIncenter(cell, point);

          outputSeparatrices1_points_->push_back(point[0]);
          outputSeparatrices1_points_->push_back(point[1]);
          outputSeparatrices1_points_->push_back(point[2]);

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
          }

          oldPointId=pointId;
          ++pointId;
        }
      }

      ++separatrixId;
    }
  }

  (*outputSeparatrices1_numberOfPoints_)=pointId;
  (*outputSeparatrices1_numberOfCells_)=cellId;

  return 0;
}

template<typename dataType>
int MorseSmaleComplex2D::execute(){
  Timer t;

  int* ascendingManifold=static_cast<int*>(outputAscendingManifold_);
  int* descendingManifold=static_cast<int*>(outputDescendingManifold_);
  int* morseSmaleManifold=static_cast<int*>(outputMorseSmaleManifold_);

  discreteGradient_.setDebugLevel(debugLevel_);
  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.buildGradient<dataType>();
  discreteGradient_.reverseGradient<dataType>();

  vector<Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);
  discreteGradient_.setCriticalPoints<dataType>(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1 or ComputeAscendingSeparatrices1){
    vector<Separatrix> separatrices;
    vector<vector<Cell>> separatricesGeometry;

    getSeparatrices(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices<dataType>(separatrices, separatricesGeometry);
  }

  {
    int numberOfMaxima{};
    int numberOfMinima{};

    if(ComputeAscendingSegmentation)
      setAscendingSegmentation(criticalPoints, ascendingManifold, numberOfMaxima);

    if(ComputeDescendingSegmentation)
      setDescendingSegmentation(criticalPoints, descendingManifold, numberOfMinima);

    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation and ComputeFinalSegmentation)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold, descendingManifold, morseSmaleManifold);
  }

  {
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    stringstream msg;
    msg << "[MorseSmaleComplex] Data-set (" << numberOfVertices
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // MORSESMALECOMPLEX2D_H
