/// \ingroup base
/// \class ttk::MorseSmaleComplex3D
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

  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex{

    public:

      MorseSmaleComplex3D();
      ~MorseSmaleComplex3D();

      template<typename dataType>
        int execute();

      template<typename dataType>
        int computePersistencePairs(const vector<tuple<int, int, dataType>>& JTPairs,
            const vector<tuple<int, int, dataType>>& STPairs,
            vector<tuple<int,int,dataType>>& pl_saddleSaddlePairs);

      template <typename dataType>
        int setAugmentedCriticalPoints(const vector<Cell>& criticalPoints,
            int* ascendingManifold,
            int* descendingManifold) const;

      int getSeparatrices1(const vector<Cell>& criticalPoints,
          vector<Separatrix>& separatrices,
          vector<vector<Cell>>& separatricesGeometry) const;

      template<typename dataType>
      int setSeparatrices1(const vector<Separatrix>& separatrices,
          const vector<vector<Cell>>& separatricesGeometry) const;

      int getSaddleConnectors(const vector<Cell>& criticalPoints,
          vector<Separatrix>& separatrices,
          vector<vector<Cell>>& separatricesGeometry) const;

      template<typename dataType>
      int setSaddleConnectors(const vector<Separatrix>& separatrices,
          const vector<vector<Cell>>& separatricesGeometry) const;

      int getDescendingSeparatrices2(const vector<Cell>& criticalPoints,
          vector<Separatrix>& separatrices,
          vector<vector<Cell>>& separatricesGeometry,
          vector<set<int>>& separatricesSaddles) const;

      template<typename dataType>
      int setDescendingSeparatrices2(const vector<Separatrix>& separatrices,
          const vector<vector<Cell>>& separatricesGeometry,
          const vector<set<int>>& separatricesSaddles) const;

      int getDualPolygon(const int edgeId, vector<int>& polygon) const;

      int sortDualPolygonVertices(vector<int>& polygon) const;

      int getAscendingSeparatrices2(const vector<Cell>& criticalPoints,
          vector<Separatrix>& separatrices,
          vector<vector<Cell>>& separatricesGeometry,
          vector<set<int>>& separatricesSaddles) const;

      template<typename dataType>
      int setAscendingSeparatrices2(const vector<Separatrix>& separatrices,
          const vector<vector<Cell>>& separatricesGeometry,
          const vector<set<int>>& separatricesSaddles) const;

      int setAscendingSegmentation(const vector<Cell>& criticalPoints,
          vector<int>& maxSeeds,
          int* const morseSmaleManifold,
          int& numberOfMaxima);

      int setDescendingSegmentation(const vector<Cell>& criticalPoints,
          int* const morseSmaleManifold,
          int& numberOfMinima);

      int setFinalSegmentation(const int numberOfMaxima,
          const int numberOfMinima,
          const int* const ascendingManifold,
          const int* const descendingManifold,
          int* morseSmaleManifold) const;

  };
}

template<typename dataType>
int MorseSmaleComplex3D::setSeparatrices1(const vector<Separatrix>& separatrices,
    const vector<vector<Cell>>& separatricesGeometry) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMaxima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMaxima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMinima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMinima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionDiffs=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices1_numberOfPoints_);
  int cellId=(*outputSeparatrices1_numberOfCells_);
  int separatrixId=separatrices.size();

  const int dimensionality=inputTriangulation_->getCellVertexNumber(0)-1;

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
int MorseSmaleComplex3D::setSaddleConnectors(const vector<Separatrix>& separatrices,
    const vector<vector<Cell>>& separatricesGeometry) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMaxima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMaxima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionMinima=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionMinima_);
  vector<dataType>* outputSeparatrices1_cells_separatrixFunctionDiffs=
    static_cast<vector<dataType>*>(outputSeparatrices1_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices1_numberOfPoints_);
  int cellId=(*outputSeparatrices1_numberOfCells_);
  int separatrixId=separatrices.size();

  for(const Separatrix& separatrix : separatrices){
    if(separatrix.isValid_){
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
int MorseSmaleComplex3D::setAscendingSeparatrices2(const vector<Separatrix>& separatrices,
   const vector<vector<Cell>>& separatricesGeometry,
   const vector<set<int>>& separatricesSaddles) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices2_numberOfPoints_);
  int cellId=(*outputSeparatrices2_numberOfCells_);
  int separatrixId=separatrices.size();

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  vector<int> isVisited(numberOfCells, -1);

  int i=0;
  for(const Separatrix& separatrix : separatrices){
    if(separatrix.isValid_){
      const Cell& saddle=separatrix.source_;
      const char separatrixType=1;
      const int saddleId=saddle.id_;

      const dataType separatrixFunctionMinimum=discreteGradient_.scalarMax<dataType>(saddle, scalars);
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

      for(const int geometryId : separatrix.geometry_){
        for(const Cell& edge : separatricesGeometry[geometryId]){
          const int edgeId=edge.id_;

          // Transform to dual : edge -> polygon
          vector<int> polygon;
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
          }
        }
      }

      ++separatrixId;
    }
    ++i;
  }

  (*outputSeparatrices2_numberOfPoints_)=pointId;
  (*outputSeparatrices2_numberOfCells_)=cellId;

  return 0;
}

template<typename dataType>
int MorseSmaleComplex3D::setDescendingSeparatrices2(const vector<Separatrix>& separatrices,
   const vector<vector<Cell>>& separatricesGeometry,
   const vector<set<int>>& separatricesSaddles) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  int pointId=(*outputSeparatrices2_numberOfPoints_);
  int cellId=(*outputSeparatrices2_numberOfCells_);
  int separatrixId=separatrices.size();

  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  vector<int> isVisited(numberOfVertices, -1);

  int i=0;
  for(const Separatrix& separatrix : separatrices){
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
        }
      }

      ++separatrixId;
    }
    ++i;
  }

  (*outputSeparatrices2_numberOfPoints_)=pointId;
  (*outputSeparatrices2_numberOfCells_)=cellId;

  return 0;
}

template<typename dataType>
int MorseSmaleComplex3D::execute(){
  Timer t;

  int* ascendingManifold=static_cast<int*>(outputAscendingManifold_);
  int* descendingManifold=static_cast<int*>(outputDescendingManifold_);
  int* morseSmaleManifold=static_cast<int*>(outputMorseSmaleManifold_);

  discreteGradient_.setDebugLevel(debugLevel_);
  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.setCollectPersistencePairs(false);
  discreteGradient_.setReturnSaddleConnectors(ReturnSaddleConnectors);
  discreteGradient_.setSaddleConnectorsPersistenceThreshold(
      SaddleConnectorsPersistenceThreshold);
  discreteGradient_.buildGradient<dataType>();
  discreteGradient_.buildGradient2<dataType>();
  discreteGradient_.buildGradient3<dataType>();
  discreteGradient_.reverseGradient<dataType>();

  vector<Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1 or ComputeAscendingSeparatrices1){
    vector<Separatrix> separatrices;
    vector<vector<Cell>> separatricesGeometry;

    getSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);
  }

  // saddle-connectors
  if(ComputeSaddleConnectors){
    vector<Separatrix> separatrices;
    vector<vector<Cell>> separatricesGeometry;
    getSaddleConnectors(criticalPoints, separatrices, separatricesGeometry);
    setSaddleConnectors<dataType>(separatrices, separatricesGeometry);
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2){
    vector<Separatrix> separatrices;
    vector<vector<Cell>> separatricesGeometry;
    vector<set<int>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    setDescendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
  }

  if(ComputeAscendingSeparatrices2){
    vector<Separatrix> separatrices;
    vector<vector<Cell>> separatricesGeometry;
    vector<set<int>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
    setAscendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
  }

  vector<int> maxSeeds;
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

template<typename dataType>
int MorseSmaleComplex3D::computePersistencePairs(const vector<tuple<int, int, dataType>>& JTPairs,
    const vector<tuple<int, int, dataType>>& STPairs,
    vector<tuple<int,int,dataType>>& pl_saddleSaddlePairs){
  dataType* scalars=static_cast<dataType*>(inputScalarField_);
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

  // get original list of critical points
  vector<pair<int,char>> pl_criticalPoints;
  {
    const int* const offsets=static_cast<int*>(inputOffsets_);
    vector<int> sosOffsets(numberOfVertices);
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
  vector<char> isAccepted(numberOfVertices, false);
  for(const auto& i : JTPairs){
    const int v0=get<0>(i);
    const int v1=get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }
  for(const auto& i : STPairs){
    const int v0=get<0>(i);
    const int v1=get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }

  // filter the critical points according to the filtering list and boundary condition
  vector<char> isSaddle1(numberOfVertices, false);
  vector<char> isSaddle2(numberOfVertices, false);
  vector<pair<int,char>> pl_filteredCriticalPoints;
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

  vector<tuple<Cell,Cell>> dmt_pairs;
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
    const Cell& saddle1=get<0>(i);
    const Cell& saddle2=get<1>(i);

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
        pl_saddleSaddlePairs.push_back(make_tuple(v0,v1,persistence));
      }
    }
  }
  return 0;
}

#endif // MORSESMALECOMPLEX3D_H
