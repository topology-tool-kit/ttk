/// \ingroup base
/// \class ttk::MorseSmaleComplex3D::
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex3D processing package.
///
/// %MorseSmaleComplex3D is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex3D.cpp %for a usage example.

#ifndef _MORSESMALECOMPLEX3D_H
#define _MORSESMALECOMPLEX3D_H

// base code includes
#include<AbstractMorseSmaleComplex.h>

namespace ttk{

  /**
   * Class specialized in building the Morse-Smale complex
   * of 3D triangulation.
   */
  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex{

    public:

      MorseSmaleComplex3D();
      ~MorseSmaleComplex3D();

      /**
       * Main function for computing the whole Morse-Smale complex.
       */
      template<typename dataType>
        int execute();

      /**
       * Compute the (saddle1, saddle2) pairs not detected by the
       * contour tree.
       */
      template<typename dataType>
        int computePersistencePairs(const std::vector<std::tuple<int, int, dataType>>& JTPairs,
            const std::vector<std::tuple<int, int, dataType>>& STPairs,
            std::vector<std::tuple<int,int,dataType>>& pl_saddleSaddlePairs);

      template <typename dataType>
        int setAugmentedCriticalPoints(const std::vector<Cell>& criticalPoints,
            int* ascendingManifold,
            int* descendingManifold) const;

      /**
       * Compute the descending 1-separatrices by reading into the discrete
       * gradient.
       */
      int getAscendingSeparatrices1(const std::vector<Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<Cell>>& separatricesGeometry) const;

      /**
       * Compute the saddle-connectors by reading into the discrete
       * gradient.
       */
      int getSaddleConnectors(const std::vector<Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<Cell>>& separatricesGeometry) const;

      /**
       * Compute the geometrical embedding of the saddle-connectors.
       */
      template<typename dataType>
      int setSaddleConnectors(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<Cell>>& separatricesGeometry) const;

      /**
       * Compute the 2-separatrices by reading into the discrete
       * gradient from the maxima.
       */
      int getDescendingSeparatrices2(const std::vector<Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<Cell>>& separatricesGeometry,
          std::vector<std::set<int>>& separatricesSaddles) const;

      /**
       * Compute the geometrical embedding of the descending
       * 2-separatrices.
       */
      template<typename dataType>
      int setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<Cell>>& separatricesGeometry,
          const std::vector<std::set<int>>& separatricesSaddles) const;
      template<typename dataType>
        int omp_setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
            const std::vector<std::vector<Cell>>& separatricesGeometry,
            const std::vector<std::set<int>>& separatricesSaddles) const;

      int getDualPolygon(const int edgeId, std::vector<int>& polygon) const;

      int sortDualPolygonVertices(std::vector<int>& polygon) const;

      /**
       * Compute the 2-separatrices by reading into the discrete
       * gradient from the minima.
       */
      int getAscendingSeparatrices2(const std::vector<Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<Cell>>& separatricesGeometry,
          std::vector<std::set<int>>& separatricesSaddles) const;

      /**
       * Compute the geometrical embedding of the ascending
       * 2-separatrices.
       */
      template<typename dataType>
      int setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<Cell>>& separatricesGeometry,
          const std::vector<std::set<int>>& separatricesSaddles) const;
  };
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::setSaddleConnectors(const 
std::vector<Separatrix>& separatrices,
    const std::vector<std::vector<Cell>>& separatricesGeometry) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices1_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices1_numberOfPoints_);
  int cellId=(*outputSeparatrices1_numberOfCells_);
  int separatrixId=0;
  if(outputSeparatrices1_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
        outputSeparatrices1_cells_separatrixIds_->end())+1;
  }

  for(const Separatrix& separatrix : separatrices){
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const Cell& saddle1=separatrix.source_;
    const Cell& saddle2=separatrix.destination_;

    // get separatrix type : saddle-connector
    const char separatrixType=1;

    // compute separatrix function diff
    const dataType separatrixFunctionMaximum=std::max(discreteGradient_.scalarMax<dataType>(saddle1, scalars),
        discreteGradient_.scalarMax<dataType>(saddle2, scalars));
    const dataType separatrixFunctionMinimum=std::min(discreteGradient_.scalarMin<dataType>(saddle1, scalars),
        discreteGradient_.scalarMin<dataType>(saddle2, scalars));
    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    // get boundary condition
    const char isOnBoundary=(discreteGradient_.isBoundary(saddle1) and discreteGradient_.isBoundary(saddle2));

    bool isFirst=true;
    for(const int geometryId : separatrix.geometry_){
      int oldPointId=-1;
      for(auto cellIte=separatricesGeometry[geometryId].begin(); cellIte!=separatricesGeometry[geometryId].end(); ++cellIte){
        const Cell& cell=*cellIte;
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

          outputSeparatrices1_cells_sourceIds_->push_back(saddle1.id_);
          outputSeparatrices1_cells_destinationIds_->push_back(saddle2.id_);
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

template<typename dataType>
int ttk::MorseSmaleComplex3D::setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
   const std::vector<std::vector<Cell>>& separatricesGeometry,
   const std::vector<std::set<int>>& separatricesSaddles) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices2_numberOfPoints_);
  int cellId=(*outputSeparatrices2_numberOfCells_);
  int separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  std::vector<int> isVisited(numberOfCells, -1);

  const int numberOfSeparatrices=separatrices.size();
  for(int i=0; i<numberOfSeparatrices; ++i){
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const Cell& saddle=separatrix.source_;
    const char separatrixType=1;
    const int saddleId=saddle.id_;

    const dataType separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(saddle, scalars);
    dataType separatrixFunctionMaximum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const int saddle2Id : separatricesSaddles[i]){
      if(inputTriangulation_->isTriangleOnBoundary(saddle2Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(Cell(2,saddle2Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMaximum=std::max(separatrixFunctionMaximum,
            discreteGradient_.scalarMax<dataType>(Cell(2,saddle2Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const int geometryId : separatrix.geometry_){
      for(const Cell& edge : separatricesGeometry[geometryId]){
        const int edgeId=edge.id_;

        // Transform to dual : edge -> polygon
        std::vector<int> polygon;
        getDualPolygon(edgeId, polygon);

        const int vertexNumber=polygon.size();
        if(vertexNumber>2){
          sortDualPolygonVertices(polygon);

          // add the polygon
          outputSeparatrices2_cells_->push_back(vertexNumber);

          float point[3];
          for(int i=0; i<vertexNumber; ++i){
            const int tetraId=polygon[i];
            discreteGradient_.getCellIncenter(Cell(3,tetraId), point);

            if(isVisited[tetraId]==-1){
              outputSeparatrices2_points_->push_back(point[0]);
              outputSeparatrices2_points_->push_back(point[1]);
              outputSeparatrices2_points_->push_back(point[2]);

              outputSeparatrices2_cells_->push_back(pointId);

              isVisited[tetraId]=pointId;
              ++pointId;
            }
            else
              outputSeparatrices2_cells_->push_back(isVisited[tetraId]);
          }

          outputSeparatrices2_cells_sourceIds_->push_back(saddleId);
          outputSeparatrices2_cells_separatrixIds_->push_back(separatrixId);
          outputSeparatrices2_cells_separatrixTypes_->push_back(separatrixType);
          outputSeparatrices2_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
          outputSeparatrices2_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
          outputSeparatrices2_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
          outputSeparatrices2_cells_isOnBoundary_->push_back(isOnBoundary);

          ++cellId;
          isFirst=false;
        }
      }
    }

    if(!isFirst)
      ++separatrixId;
  }

  (*outputSeparatrices2_numberOfPoints_)=pointId;
  (*outputSeparatrices2_numberOfCells_)=cellId;

  return 0;
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::omp_setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
    const std::vector<std::vector<Cell>>& separatricesGeometry,
    const std::vector<std::set<int>>& separatricesSaddles) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices2_numberOfPoints_);
  int cellId=(*outputSeparatrices2_numberOfCells_);
  int separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  std::vector<int> separatrixIds(threadNumber_, 0);
  std::vector<int> numberOfPoints(threadNumber_, 0);
  std::vector<std::vector<float>> separatrices2_points(threadNumber_);
  std::vector<int> numberOfCells(threadNumber_, 0);
  std::vector<std::vector<int>> separatrices2_cells(threadNumber_);
  std::vector<std::vector<int>> separatrices2_cells_sourceIds(threadNumber_);
  std::vector<std::vector<int>> separatrices2_cells_separatrixIds(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_separatrixTypes(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMaxima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMinima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionDiffs(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_isOnBoundary(threadNumber_);

  const int numberOfSeparatrices=separatrices.size();
#pragma omp parallel for num_threads(threadNumber_)
  for(int i=0; i<numberOfSeparatrices; ++i){
    const int threadId=omp_get_thread_num();
    const Separatrix& separatrix=separatrices[i];
    if(separatrix.isValid_){
      const Cell& saddle=separatrix.source_;
      const char separatrixType=2;
      const int saddleId=saddle.id_;

      const dataType separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(saddle, scalars);
      dataType separatrixFunctionMinimum{};

      // get separatrix infos
      char isOnBoundary{};
      bool isFirst=true;
      for(const int saddle1Id : separatricesSaddles[i]){
        if(inputTriangulation_->isEdgeOnBoundary(saddle1Id))
          ++isOnBoundary;

        if(isFirst){
          separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(Cell(1,saddle1Id), scalars);
          isFirst=false;
        }
        else{
          separatrixFunctionMinimum=std::min(separatrixFunctionMinimum,
              discreteGradient_.scalarMin<dataType>(Cell(1,saddle1Id), scalars));
        }
      }

      const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

      isFirst=true;
      for(const int geometryId : separatrix.geometry_){
        for(const Cell& cell : separatricesGeometry[geometryId]){
          const int triangleId=cell.id_;

          separatrices2_cells[threadId].push_back(3);
          for(int k=0; k<3; ++k){
            int vertexId;
            inputTriangulation_->getTriangleVertex(triangleId, k, vertexId);
            float point[3];
            inputTriangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);

            separatrices2_points[threadId].push_back(point[0]);
            separatrices2_points[threadId].push_back(point[1]);
            separatrices2_points[threadId].push_back(point[2]);
            separatrices2_cells[threadId].push_back(numberOfPoints[threadId]);
            ++numberOfPoints[threadId];
          }
          separatrices2_cells_sourceIds[threadId].push_back(saddleId);
          separatrices2_cells_separatrixIds[threadId].push_back(separatrixIds[threadId]);
          separatrices2_cells_separatrixTypes[threadId].push_back(separatrixType);
          separatrices2_cells_separatrixFunctionMaxima[threadId].push_back(separatrixFunctionMaximum);
          separatrices2_cells_separatrixFunctionMinima[threadId].push_back(separatrixFunctionMinimum);
          separatrices2_cells_separatrixFunctionDiffs[threadId].push_back(separatrixFunctionDiff);
          separatrices2_cells_isOnBoundary[threadId].push_back(isOnBoundary);
          ++numberOfCells[threadId];

          isFirst=false;
        }
      }

      if(!isFirst)
        ++separatrixIds[threadId];
    }
  }

  const int oldPointSize=outputSeparatrices2_points_->size();
  const int oldCellSize=outputSeparatrices2_cells_->size();
  const int oldFieldSize=outputSeparatrices2_cells_sourceIds_->size();
  {
    int npoints=0;
    int ncells=0;
    int nnpoints=pointId;
    int nncells=cellId;
    std::vector<int> offsetPoints(threadNumber_, 0);
    std::vector<int> offsetCells(threadNumber_, 0);
    std::vector<int> offsetNPoints(threadNumber_, 0);
    std::vector<int> offsetNCells(threadNumber_, 0);

    for(int i=0; i<threadNumber_; ++i){
      offsetPoints[i]=npoints;
      offsetCells[i]=ncells;
      offsetNPoints[i]=nnpoints;
      offsetNCells[i]=nncells;

      npoints+=separatrices2_points[i].size();
      ncells+=separatrices2_cells[i].size();
      nnpoints+=numberOfPoints[i];
      nncells+=numberOfCells[i];
    }

    outputSeparatrices2_points_->resize(oldPointSize+npoints);
    outputSeparatrices2_cells_->resize(oldCellSize+ncells);
    outputSeparatrices2_cells_sourceIds_->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_separatrixIds_->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_separatrixTypes_->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_separatrixFunctionMaxima->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_separatrixFunctionMinima->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_separatrixFunctionDiffs->resize(oldFieldSize+nncells);
    outputSeparatrices2_cells_isOnBoundary_->resize(oldFieldSize+nncells);
#pragma omp parallel for num_threads(threadNumber_)
    for(int i=0; i<threadNumber_; ++i){
      // reduce: points
      const int tmp_npoints=separatrices2_points[i].size();
      const int tmp_offsetPoints=offsetPoints[i];
      for(int j=0; j<tmp_npoints; ++j)
        (*outputSeparatrices2_points_)[oldPointSize+tmp_offsetPoints+j]=separatrices2_points[i][j];

      // reduce: cells
      const int tmp_ncells=separatrices2_cells[i].size();
      const int tmp_offsetCells=offsetCells[i];
      for(int j=0; j<tmp_ncells;){
        const int cellSize=separatrices2_cells[i][j];
        (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j]=cellSize;
        for(int k=0; k<cellSize; ++k){
          const int tmp_pointId=separatrices2_cells[i][j+k+1];
          (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j+k+1]=offsetNPoints[i]+tmp_pointId;
        }
        j+=(cellSize+1);
      }

      // reduce: fields
      for(int j=0; j<numberOfCells[i]; ++j){
        (*outputSeparatrices2_cells_sourceIds_)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_sourceIds[i][j];
        (*outputSeparatrices2_cells_separatrixIds_)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_separatrixIds[i][j];
        (*outputSeparatrices2_cells_separatrixTypes_)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_separatrixTypes[i][j];
        (*outputSeparatrices2_cells_separatrixFunctionMaxima)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_separatrixFunctionMaxima[i][j];
        (*outputSeparatrices2_cells_separatrixFunctionMinima)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_separatrixFunctionMinima[i][j];
        (*outputSeparatrices2_cells_separatrixFunctionDiffs)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_separatrixFunctionDiffs[i][j];
        (*outputSeparatrices2_cells_isOnBoundary_)[oldFieldSize+offsetNCells[i]+j]=
          separatrices2_cells_isOnBoundary[i][j];
      }
    }
  }

  {
    int totalNumberOfPoints=0;
    int totalNumberOfCells=0;
    for(int i=0; i<threadNumber_; ++i){
      totalNumberOfPoints+=numberOfPoints[i];
      totalNumberOfCells+=numberOfCells[i];
    }
    (*outputSeparatrices2_numberOfPoints_)+=totalNumberOfPoints;
    (*outputSeparatrices2_numberOfCells_)+=totalNumberOfCells;
  }

  return 0;
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
   const std::vector<std::vector<Cell>>& separatricesGeometry,
   const std::vector<std::set<int>>& separatricesSaddles) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices2_numberOfPoints_);
  int cellId=(*outputSeparatrices2_numberOfCells_);
  int separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::vector<int> isVisited(numberOfVertices, -1);

  const int numberOfSeparatrices=separatrices.size();
  for(int i=0; i<numberOfSeparatrices; ++i){
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const Cell& saddle=separatrix.source_;
    const char separatrixType=2;
    const int saddleId=saddle.id_;

    const dataType separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(saddle, scalars);
    dataType separatrixFunctionMinimum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const int saddle1Id : separatricesSaddles[i]){
      if(inputTriangulation_->isEdgeOnBoundary(saddle1Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(Cell(1,saddle1Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMinimum=std::min(separatrixFunctionMinimum,
            discreteGradient_.scalarMin<dataType>(Cell(1,saddle1Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const int geometryId : separatrix.geometry_){
      for(const Cell& cell : separatricesGeometry[geometryId]){
        const int triangleId=cell.id_;

        outputSeparatrices2_cells_->push_back(3);
        float point[3];
        for(int k=0; k<3; ++k){
          int vertexId;
          inputTriangulation_->getTriangleVertex(triangleId, k, vertexId);

          if(isVisited[vertexId]==-1){
            inputTriangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);

            outputSeparatrices2_points_->push_back(point[0]);
            outputSeparatrices2_points_->push_back(point[1]);
            outputSeparatrices2_points_->push_back(point[2]);

            outputSeparatrices2_cells_->push_back(pointId);

            isVisited[vertexId]=pointId;
            ++pointId;
          }
          else
            outputSeparatrices2_cells_->push_back(isVisited[vertexId]);
        }
        outputSeparatrices2_cells_sourceIds_->push_back(saddleId);
        outputSeparatrices2_cells_separatrixIds_->push_back(separatrixId);
        outputSeparatrices2_cells_separatrixTypes_->push_back(separatrixType);
        outputSeparatrices2_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
        outputSeparatrices2_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
        outputSeparatrices2_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
        outputSeparatrices2_cells_isOnBoundary_->push_back(isOnBoundary);

        ++cellId;
        isFirst=false;
      }
    }

    if(!isFirst)
      ++separatrixId;
  }

  (*outputSeparatrices2_numberOfPoints_)=pointId;
  (*outputSeparatrices2_numberOfCells_)=cellId;

  return 0;
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::execute(){
  Timer t;

  int* ascendingManifold=static_cast<int*>(outputAscendingManifold_);
  int* descendingManifold=static_cast<int*>(outputDescendingManifold_);
  int* morseSmaleManifold=static_cast<int*>(outputMorseSmaleManifold_);

  discreteGradient_.buildGradient<dataType>();
  discreteGradient_.buildGradient2<dataType>();
  discreteGradient_.buildGradient3<dataType>();
  discreteGradient_.reverseGradient<dataType>();

  std::vector<Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1){
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    getDescendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);
  }

  if(ComputeAscendingSeparatrices1){
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    getAscendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);
  }

  // saddle-connectors
  if(ComputeSaddleConnectors){
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    getSaddleConnectors(criticalPoints, separatrices, separatricesGeometry);
    setSaddleConnectors<dataType>(separatrices, separatricesGeometry);
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2){
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    std::vector<std::set<int>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    omp_setDescendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
  }

  if(ComputeAscendingSeparatrices2){
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<Cell>> separatricesGeometry;
    std::vector<std::set<int>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    setAscendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
  }

  std::vector<int> maxSeeds;
  {
    int numberOfMaxima{};
    int numberOfMinima{};

    if(ComputeAscendingSegmentation)
      setAscendingSegmentation(criticalPoints, maxSeeds, ascendingManifold, numberOfMaxima);

    if(ComputeDescendingSegmentation)
      setDescendingSegmentation(criticalPoints, descendingManifold, numberOfMinima);

    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation and ComputeFinalSegmentation)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold, descendingManifold, morseSmaleManifold);
  }

  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation)
    discreteGradient_.setAugmentedCriticalPoints<dataType>(criticalPoints,
        maxSeeds,
        ascendingManifold,
        descendingManifold);
  else
    discreteGradient_.setCriticalPoints<dataType>(criticalPoints);

  {
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    std::stringstream msg;
    msg << "[MorseSmaleComplex3D] Data-set (" << numberOfVertices
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::computePersistencePairs(const std::vector<std::tuple<int, int, dataType>>& JTPairs,
    const std::vector<std::tuple<int, int, dataType>>& STPairs,
    std::vector<std::tuple<int,int,dataType>>& pl_saddleSaddlePairs){
  dataType* scalars=static_cast<dataType*>(inputScalarField_);
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

  // get original list of critical points
  std::vector<std::pair<int,char>> pl_criticalPoints;
  {
    const int* const offsets=static_cast<int*>(inputOffsets_);
    std::vector<int> sosOffsets(numberOfVertices);
    for(int i=0; i<numberOfVertices; ++i)
      sosOffsets[i]=offsets[i];

    ScalarFieldCriticalPoints<dataType> scp;

    scp.setDebugLevel(debugLevel_);
    scp.setThreadNumber(threadNumber_);
    scp.setDomainDimension(inputTriangulation_->getDimensionality());
    scp.setScalarValues(inputScalarField_);
    scp.setVertexNumber(numberOfVertices);
    scp.setSosOffsets(&sosOffsets);
    scp.setupTriangulation(inputTriangulation_);
    scp.setOutput(&pl_criticalPoints);

    scp.execute();
  }

  // build accepting list
  std::vector<char> isAccepted(numberOfVertices, false);
  for(const auto& i : JTPairs){
    const int v0=std::get<0>(i);
    const int v1=std::get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }
  for(const auto& i : STPairs){
    const int v0=std::get<0>(i);
    const int v1=std::get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }

  // filter the critical points according to the filtering list and boundary condition
  std::vector<char> isSaddle1(numberOfVertices, false);
  std::vector<char> isSaddle2(numberOfVertices, false);
  std::vector<std::pair<int,char>> pl_filteredCriticalPoints;
  for(const auto& i : pl_criticalPoints){
    const int vertexId=i.first;
    const char type=i.second;
    if(isAccepted[vertexId]){
      pl_filteredCriticalPoints.push_back(i);

      switch(type){
        case 1:
          isSaddle1[vertexId]=true;
          break;

        case 2:
          isSaddle2[vertexId]=true;
          break;
      }
    }
  }

  std::vector<std::tuple<Cell,Cell>> dmt_pairs;
  {
    // simplify to be PL-conformant
    discreteGradient_.setDebugLevel(debugLevel_);
    discreteGradient_.setThreadNumber(threadNumber_);
    discreteGradient_.setReverseSaddleMaximumConnection(true);
    discreteGradient_.setReverseSaddleSaddleConnection(true);
    discreteGradient_.setCollectPersistencePairs(false);
    discreteGradient_.buildGradient<dataType>();
    discreteGradient_.buildGradient2<dataType>();
    discreteGradient_.buildGradient3<dataType>();
    discreteGradient_.reverseGradient<dataType>(pl_criticalPoints);

    // collect saddle-saddle connections
    discreteGradient_.setReverseSaddleMaximumConnection(true);
    discreteGradient_.setCollectPersistencePairs(true);
    discreteGradient_.setOutputPersistencePairs(&dmt_pairs);
    discreteGradient_.reverseGradient<dataType>(pl_filteredCriticalPoints);
  }

  // transform DMT pairs into PL pairs
  for(const auto& i : dmt_pairs){
    const Cell& saddle1=std::get<0>(i);
    const Cell& saddle2=std::get<1>(i);

    int v0=-1;
    for(int j=0; j<2; ++j){
      int vertexId;
      inputTriangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

      if(isSaddle1[vertexId]){
        v0=vertexId;
        break;
      }
    }
    if(v0==-1){
      dataType scalar{};
      for(int j=0; j<2; ++j){
        int vertexId;
        inputTriangulation_->getEdgeVertex(saddle1.id_,j,vertexId);
        const dataType vertexScalar=scalars[vertexId];

        if(!j or scalar>vertexScalar){
          v0=vertexId;
          scalar=vertexScalar;
        }
      }
    }

    int v1=-1;
    for(int j=0; j<3; ++j){
      int vertexId;
      inputTriangulation_->getTriangleVertex(saddle2.id_, j, vertexId);

      if(isSaddle2[vertexId]){
        v1=vertexId;
        break;
      }
    }
    if(v1==-1){
      dataType scalar{};
      for(int j=0; j<3; ++j){
        int vertexId;
        inputTriangulation_->getTriangleVertex(saddle2.id_,j,vertexId);
        const dataType vertexScalar=scalars[vertexId];

        if(!j or scalar<vertexScalar){
          v1=vertexId;
          scalar=vertexScalar;
        }
      }
    }

    const dataType persistence=scalars[v1]-scalars[v0];

    if(v0!=-1 and v1!=-1 and persistence>=0){
      if(!inputTriangulation_->isVertexOnBoundary(v0) or !inputTriangulation_->isVertexOnBoundary(v1)){
        pl_saddleSaddlePairs.push_back(std::make_tuple(v0,v1,persistence));
      }
    }
  }
  return 0;
}

#endif // MORSESMALECOMPLEX3D_H
