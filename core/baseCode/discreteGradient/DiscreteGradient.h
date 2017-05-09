/// \ingroup baseCode
/// \class ttk::DiscreteGradient
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2016.
///
/// \brief TTK %discreteGradient processing package.
///
/// %DiscreteGradient is a TTK processing package that handles discrete gradient
/// (in the sense of Discrete Morse Theory).
///
/// \sa ttk::Triangulation
/// \sa vtkDiscreteGradient.cpp %for a usage example.

#ifndef _DISCRETEGRADIENT_H
#define _DISCRETEGRADIENT_H

// base code includes
#include<Triangulation.h>
#include<Geometry.h>
#include<Wrapper.h>
#include<ScalarFieldCriticalPoints.h>

#include<queue>
#include<algorithm>
#include<set>
#include<array>

// ### DEBUG PURPOSE ONLY : these flags enable/disable different simplification implementations ###
// naive version based only on gradient modification (slowest version)
//#define NAIVE_VERSION

// hybrid version based on graph and look at the gradient for update
#define PROTO_VERSION

// ### DEBUG PURPOSE ONLY : these flags enable/disable different features ###
// debug feature : automatic cycle detection
//#define ALLOW_CYCLE_DETECTOR

// debug feature : automatic double-pairing detection
//#define ALLOW_DOUBLE_PAIRING_DETECTOR

// debug feature : automatic PL-compliant detection
//#define ALLOW_PL_COMPLIANT_DETECTOR

// debug feature : exit(-1) when a cycle is detected
//#define ALLOW_EXIT

// debug feature: cout some debug informations
//#define PRINT_DEBUG

// debug feature: cout error messages
#define PRINT_ERROR

// print some additional infos
//#define PRINT_INFOS

namespace ttk{
  using wallId_t=unsigned long long int;

  struct Cell{
    Cell():
      dim_{-1},
      id_{-1}
    {}

    Cell(const int dim,
        const int id):
      dim_{dim},
      id_{id}
    {}

    Cell(const Cell& cell):
      dim_{cell.dim_},
      id_{cell.id_}
    {}

    int dim_;
    int id_;
  };

  struct Segment{
    Segment():
      orientation_{},
      isValid_{}
    {}

    Segment(const bool orientation,
        const vector<Cell>& cells,
        const bool isValid):
      orientation_{orientation},
      cells_{cells},
      isValid_{isValid}
    {}

    Segment(const bool orientation,
        vector<Cell>&& cells,
        const bool isValid):
      orientation_{orientation},
      cells_{cells},
      isValid_{isValid}
    {}

    Segment(const Segment& segment):
      orientation_{segment.orientation_},
      cells_{segment.cells_},
      isValid_{segment.isValid_}
    {}

    Segment(Segment&& segment):
      orientation_{segment.orientation_},
      cells_{segment.cells_},
      isValid_{segment.isValid_}
    {}

    int invalidate(){
      isValid_=false;
      clear();

      return 0;
    }

    int clear(){
      cells_.clear();

      return 0;
    }

    bool orientation_;
    vector<Cell> cells_;
    bool isValid_;
  };

  struct VPath{
    VPath():
      isValid_{},
      source_{-1},
      destination_{-1},
      sourceSlot_{-1},
      destinationSlot_{-1}
    {}

    VPath(const bool isValid,
        const int segmentId,
        const int source,
        const int destination,
        const int sourceSlot,
        const int destinationSlot,
        const double persistence):
      isValid_{isValid},
      states_{1},
      segments_{segmentId},
      source_{source},
      destination_{destination},
      sourceSlot_{sourceSlot},
      destinationSlot_{destinationSlot},
      persistence_{persistence}
    {}

    VPath(const bool isValid,
        const vector<char>& states,
        const vector<int>& segments,
        const int source,
        const int destination,
        const int sourceSlot,
        const int destinationSlot,
        const double persistence):
      isValid_{isValid},
      states_{states},
      segments_{segments},
      source_{source},
      destination_{destination},
      sourceSlot_{sourceSlot},
      destinationSlot_{destinationSlot},
      persistence_{persistence}
    {}

    VPath(const bool isValid,
        vector<char>&& states,
        vector<int>&& segments,
        const int source,
        const int destination,
        const int sourceSlot,
        const int destinationSlot,
        const double persistence):
      isValid_{isValid},
      states_{states},
      segments_{segments},
      source_{source},
      destination_{destination},
      sourceSlot_{sourceSlot},
      destinationSlot_{destinationSlot},
      persistence_{persistence}
    {}

    VPath(const VPath& vpath):
      isValid_{vpath.isValid_},
      states_{vpath.states_},
      segments_{vpath.segments_},
      source_{vpath.source_},
      destination_{vpath.destination_},
      sourceSlot_{vpath.sourceSlot_},
      destinationSlot_{vpath.destinationSlot_},
      persistence_{vpath.persistence_}
    {}

    VPath(VPath&& vpath):
      isValid_{vpath.isValid_},
      states_{vpath.states_},
      segments_{vpath.segments_},
      source_{vpath.source_},
      destination_{vpath.destination_},
      sourceSlot_{vpath.sourceSlot_},
      destinationSlot_{vpath.destinationSlot_},
      persistence_{vpath.persistence_}
    {}

    int invalidate(){
      isValid_=false;
      source_=-1;
      destination_=-1;
      persistence_=-1;
      clear();

      return 0;
    }

    int clear(){
      states_.clear();
      segments_.clear();

      return 0;
    }

    bool isValid_;
    vector<char> states_;
    vector<int> segments_;
    int source_;
    int destination_;
    int sourceSlot_;
    int destinationSlot_;
    double persistence_;
  };

  struct CriticalPoint{
    CriticalPoint():
      numberOfSlots_{}
    {}

    CriticalPoint(const Cell& cell):
      cell_{cell},
      numberOfSlots_{}
    {}

    CriticalPoint(const Cell& cell,
        const vector<int>& vpaths):
      cell_{cell},
      vpaths_{vpaths},
      numberOfSlots_{}
    {}

    CriticalPoint(const Cell& cell,
        vector<int>&& vpaths):
      cell_{cell},
      vpaths_{vpaths},
      numberOfSlots_{}
    {}

    CriticalPoint(const CriticalPoint& criticalPoint):
      cell_{criticalPoint.cell_},
      vpaths_{criticalPoint.vpaths_},
      numberOfSlots_{criticalPoint.numberOfSlots_}
    {}

    CriticalPoint(CriticalPoint&& criticalPoint):
      cell_{criticalPoint.cell_},
      vpaths_{criticalPoint.vpaths_},
      numberOfSlots_{criticalPoint.numberOfSlots_}
    {}

    int omp_addSlot(){
      int numberOfSlots=0;

#ifdef withOpenMP
# pragma omp atomic capture
#endif
      numberOfSlots=(numberOfSlots_++);

      return numberOfSlots;
    }

    int addSlot(){
      return (numberOfSlots_++);
    }

    int clear(){
      vpaths_.clear();

      return 0;
    }

    Cell cell_;
    vector<int> vpaths_;
    int numberOfSlots_;
  };

  template <typename dataType>
    struct SaddleMaximumVPathComparator{
      bool operator()(const pair<dataType,int>& v1, const pair<dataType,int>& v2) const{
        const dataType persistence1=v1.first;
        const dataType persistence2=v2.first;

        const int vpathId1=v1.second;
        const int vpathId2=v2.second;

        if(persistence1!=persistence2)
          return (persistence1<persistence2);

        return (vpathId1<vpathId2);
      };
    };

  template <typename dataType>
    struct SaddleSaddleVPathComparator{
      bool operator()(const tuple<dataType,int,int>& v1, const tuple<dataType,int,int>& v2) const{
        const dataType persistence1=get<0>(v1);
        const dataType persistence2=get<0>(v2);

        const int vpathId1=get<1>(v1);
        const int vpathId2=get<1>(v2);

        const int saddleId1=get<2>(v1);
        const int saddleId2=get<2>(v2);

        if(persistence1!=persistence2)
          return (persistence1<persistence2);

        if(saddleId1!=saddleId2)
          return (saddleId1<saddleId2);

        return (vpathId1<vpathId2);
      };
    };

  class DiscreteGradient : public Debug{

    public:

      DiscreteGradient();
      ~DiscreteGradient();

      int setIterationThreshold(const int iterationThreshold){
        IterationThreshold=iterationThreshold;
        return 0;
      }

      int setReverseSaddleMaximumConnection(const bool state){
        ReverseSaddleMaximumConnection=state;
        return 0;
      }

      int setReverseSaddleSaddleConnection(const bool state){
        ReverseSaddleSaddleConnection=state;
        return 0;
      }

      int setAllowReversingWithNonRemovable(const bool state){
        AllowReversingWithNonRemovable=state;
        return 0;
      }

      int setCollectPersistencePairs(const bool state){
        CollectPersistencePairs=state;
        return 0;
      }

      int setOutputPersistencePairs(vector<tuple<Cell,Cell>>* const data){
        outputPersistencePairs_=data;
        return 0;
      }

      template <typename dataType>
        dataType scalarMax(const Cell& cell, const dataType* const scalars) const;

      template <typename dataType>
        dataType scalarMin(const Cell& cell, const dataType* const scalars) const;

      template <typename dataType>
        dataType getPersistence(const Cell& up, const Cell& down, const dataType* const scalars) const;

      template <typename dataType>
        bool isHigherThan(const int vertexA,
            const int vertexB,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        bool isLowerThan(const int vertexA,
            const int vertexB,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        int cellMax(const int cellDim,
            const int cellA,
            const int cellB,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        int cellMin(const int cellDim,
            const int cellA,
            const int cellB,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        int g0(const int cellDim,
            const int cellId,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        int g0_second(const int cellDim,
            const int cellId,
            const dataType* const scalars,
            const int* const offsets) const;

      template <typename dataType>
        int assignGradient(const int alphaDim,
            const dataType* const scalars,
            const int* const offsets,
            vector<vector<int>>& gradient) const;

      template <typename dataType>
        int assignGradient2(const int alphaDim,
            const dataType* const scalars,
            const int* const offsets,
            vector<vector<int>>& gradient) const;

      template <typename dataType>
        int buildGradient();

      template <typename dataType>
        int buildGradient2();

      template <typename dataType>
        int getRemovableMaxima(const vector<pair<int,char>>& criticalPoints,
            vector<char>& isRemovable) const;

      template <typename dataType>
        int getRemovableSaddles1(const vector<pair<int,char>>& criticalPoints,
            vector<char>& isRemovable) const;

      template <typename dataType>
        int getRemovableSaddles2(const vector<pair<int,char>>& criticalPoints,
            vector<char>& isRemovable) const;

      template <typename dataType>
        int initializeSaddleMaximumConnections(const vector<char>& isRemovableMaximum,
            vector<char>& isRemovableSaddle,
            const bool allowBoundary,
            vector<Segment>& segments,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints) const;

      template <typename dataType>
        int orderSaddleMaximumConnections(const vector<VPath>& vpaths,
            set<pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>& S);

      template <typename dataType>
        int computeCoefficients(const bool isDense,
            vector<char>& denseCoefficients,
            vector<Segment>& segments,
            const CriticalPoint& source,
            VPath& newVPath,
            const vector<VPath>& vpaths) const;

      template <typename dataType>
        int processSaddleMaximumConnections(const int iterationThreshold,
            set<pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>& S,
            vector<Segment>& segments,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints) const;

      template <typename dataType>
        int reverseSaddleMaximumConnections(const vector<Segment>& segments);

      template <typename dataType>
        int simplifySaddleMaximumConnections(const vector<char>& isRemovableMaximum,
            vector<char>& isRemovableSaddle,
            const int iterationThreshold,
            const bool allowBoundary);

      template <typename dataType>
        int naive_reverseSaddleSaddleConnection(vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2);

      template <typename dataType>
        int initializeSaddleSaddleConnections1(const vector<char>& isRemovableSaddle1,
            const vector<char>& isRemovableSaddle2,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            vector<int>& saddle1Index,
            vector<int>& saddle2Index) const;

      template <typename dataType>
        int orderSaddleSaddleConnections1(const vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S);

      template <typename dataType>
        int processSaddleSaddleConnections1(const int iterationThreshold,
            set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S,
            vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            vector<int>& saddle1Index,
            vector<int>& saddle2Index);

      template <typename dataType>
        int simplifySaddleSaddleConnections1(vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2,
            const int iterationThreshold);

      template <typename dataType>
        int naive_reverseSaddleSaddleConnection2(vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2);

      template <typename dataType>
        int initializeSaddleSaddleConnections2(const vector<char>& isRemovableSaddle1,
            const vector<char>& isRemovableSaddle2,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            vector<int>& saddle1Index,
            vector<int>& saddle2Index) const;

      template <typename dataType>
        int orderSaddleSaddleConnections2(const vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S);

      template <typename dataType>
        int processSaddleSaddleConnections2(const int iterationThreshold,
            set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S,
            vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2,
            vector<VPath>& vpaths,
            vector<CriticalPoint>& criticalPoints,
            vector<int>& saddle1Index,
            vector<int>& saddle2Index);

      template <typename dataType>
        int simplifySaddleSaddleConnections2(vector<char>& isRemovableSaddle1,
            vector<char>& isRemovableSaddle2,
            const int iterationThreshold);

      template<typename dataType>
        int reverseGradient(const vector<pair<int,char>>& criticalPoints);

      template <typename dataType>
        int reverseGradient();

      inline int setInputScalarField(void* const data){
        inputScalarField_=data;
        return 0;
      }

      inline int setupTriangulation(Triangulation* const data){
        inputTriangulation_=data;
        if(inputTriangulation_){
          dimensionality_=inputTriangulation_->getCellVertexNumber(0)-1;

          inputTriangulation_->preprocessBoundaryVertices();
          inputTriangulation_->preprocessBoundaryEdges();
          inputTriangulation_->preprocessVertexNeighbors();
          inputTriangulation_->preprocessVertexEdges();
          inputTriangulation_->preprocessVertexStars();
          inputTriangulation_->preprocessEdges();
          inputTriangulation_->preprocessEdgeStars();
          if(dimensionality_==2){
            inputTriangulation_->preprocessCellEdges();
          }
          else if(dimensionality_==3){
            inputTriangulation_->preprocessBoundaryTriangles();
            inputTriangulation_->preprocessVertexTriangles();
            inputTriangulation_->preprocessEdgeTriangles();
            inputTriangulation_->preprocessTriangles();
            inputTriangulation_->preprocessTriangleEdges();
            inputTriangulation_->preprocessTriangleStars();
            inputTriangulation_->preprocessCellTriangles();
          }
        }
        return 0;
      }

      inline int setInputOffsets(void* const data){
        inputOffsets_=data;
        return 0;
      }

      inline int setOutputCriticalPoints(int* const criticalPoints_numberOfPoints,
          vector<float>* const criticalPoints_points,
          vector<int>* const criticalPoints_points_cellDimensons,
          vector<int>* const criticalPoints_points_cellIds,
          void* const criticalPoints_points_cellScalars,
          vector<char>* const criticalPoints_points_isOnBoundary){
        outputCriticalPoints_numberOfPoints_=criticalPoints_numberOfPoints;
        outputCriticalPoints_points_=criticalPoints_points;
        outputCriticalPoints_points_cellDimensions_=criticalPoints_points_cellDimensons;
        outputCriticalPoints_points_cellIds_=criticalPoints_points_cellIds;
        outputCriticalPoints_points_cellScalars_=criticalPoints_points_cellScalars;
        outputCriticalPoints_points_isOnBoundary_=criticalPoints_points_isOnBoundary;
        return 0;
      }

      inline int setOutputGradientGlyphs(int* const gradientGlyphs_numberOfPoints,
          vector<float>* const gradientGlyphs_points,
          vector<int>* const gradientGlyphs_points_pairOrigins,
          int* const gradientGlyphs_numberOfCells,
          vector<int>* const gradientGlyphs_cells,
          vector<int>* const gradientGlyphs_cells_pairTypes){
        outputGradientGlyphs_numberOfPoints_=gradientGlyphs_numberOfPoints;
        outputGradientGlyphs_points_=gradientGlyphs_points;
        outputGradientGlyphs_points_pairOrigins_=gradientGlyphs_points_pairOrigins;
        outputGradientGlyphs_numberOfCells_=gradientGlyphs_numberOfCells;
        outputGradientGlyphs_cells_=gradientGlyphs_cells;
        outputGradientGlyphs_cells_pairTypes_=gradientGlyphs_cells_pairTypes;
        return 0;
      }

      int getDimensionality() const;

      int getNumberOfDimensions() const;

      int getNumberOfCells(const int dimension) const;

      bool isMinimum(const Cell& cell) const;

      bool isSaddle1(const Cell& cell) const;

      bool isSaddle2(const Cell& cell) const;

      bool isMaximum(const Cell& cell) const;

      bool isCellCritical(const Cell& cell) const;

      bool isBoundary(const Cell& cell) const;

      int getPairedCell(const Cell& cell, bool isReverse=false) const;

      int getCriticalPoints(vector<Cell>& criticalPoints) const;

      int getAscendingPath(const Cell& cell,
          vector<Cell>& vpath,
          const bool enableCycleDetector=false) const;

      int getDescendingPath(const Cell& cell, vector<Cell>& vpath) const;

      int getDescendingPathThroughWall(const wallId_t wallId,
          const Cell& saddle2,
          const Cell& saddle1,
          const vector<wallId_t>& isVisited,
          vector<Cell>* const vpath,
          const bool enableCycleDetector=false) const;

      bool getAscendingPathThroughWall(const wallId_t wallId,
          const Cell& saddle1,
          const Cell& saddle2,
          const vector<wallId_t>& isVisited,
          vector<Cell>* const vpath,
          const bool enableCycleDetector=false) const;

      int getDescendingWall(const wallId_t wallId,
          const Cell& cell,
          vector<wallId_t>& isVisited,
          vector<Cell>* const wall=nullptr,
          set<int>* const saddles=nullptr) const;

      int getAscendingWall(const wallId_t wallId,
          const Cell& cell,
          vector<wallId_t>& isVisited,
          vector<Cell>* const wall=nullptr,
          set<int>* const saddles=nullptr) const;

      int reverseAscendingPath(const vector<Cell>& vpath);

      int reverseAscendingPathOnWall(const vector<Cell>& vpath);

      int reverseDescendingPathOnWall(const vector<Cell>& vpath);

      int getEdgeIncenter(int edgeId, float incenter[3]) const;

      int getTriangleIncenter(int triangleId, float incenter[3]) const;

      int getTetraIncenter(int tetraId, float incenter[3]) const;

      int getCellIncenter(const Cell& cell, float incenter[3]) const;

      template <typename dataType>
        int setCriticalPoints(const vector<Cell>& criticalPoints) const;

      template <typename dataType>
        int setCriticalPoints() const;

      int setGradientGlyphs() const;

    protected:
      bool MustBeRemovable;
      bool AllowBoundary;
      bool ForceBoundary;

      int IterationThreshold;
      bool ReverseSaddleMaximumConnection;
      bool ReverseSaddleSaddleConnection;
      bool AllowReversingWithNonRemovable;
      bool CollectPersistencePairs;

      int dimensionality_;
      vector<vector<vector<int>>> gradient_;
      vector<pair<int,char>> criticalPoints_;

      void* inputScalarField_;
      void* inputOffsets_;
      Triangulation* inputTriangulation_;

      int* outputCriticalPoints_numberOfPoints_;
      vector<float>* outputCriticalPoints_points_;
      vector<int>* outputCriticalPoints_points_cellDimensions_;
      vector<int>* outputCriticalPoints_points_cellIds_;
      void* outputCriticalPoints_points_cellScalars_;
      vector<char>* outputCriticalPoints_points_isOnBoundary_;

      int* outputGradientGlyphs_numberOfPoints_;
      vector<float>* outputGradientGlyphs_points_;
      vector<int>* outputGradientGlyphs_points_pairOrigins_;
      int* outputGradientGlyphs_numberOfCells_;
      vector<int>* outputGradientGlyphs_cells_;
      vector<int>* outputGradientGlyphs_cells_pairTypes_;

      vector<tuple<Cell,Cell>>* outputPersistencePairs_;
  };
}

template <typename dataType>
dataType DiscreteGradient::scalarMax(const Cell& cell, const dataType* const scalars) const{
  dataType scalar{};

  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          int vertexId;
          inputTriangulation_->getCellVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          int vertexId;
          inputTriangulation_->getTriangleVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 3:
        for(int i=0; i<4; ++i){
          int vertexId;
          inputTriangulation_->getCellVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;
    }
  }

  return scalar;
}

template <typename dataType>
dataType DiscreteGradient::scalarMin(const Cell& cell, const dataType* const scalars) const{
  dataType scalar{};

  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          int vertexId;
          inputTriangulation_->getCellVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          int vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          int vertexId;
          inputTriangulation_->getTriangleVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 3:
        for(int i=0; i<4; ++i){
          int vertexId;
          inputTriangulation_->getCellVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;
    }
  }

  return scalar;
}

template <typename dataType>
dataType DiscreteGradient::getPersistence(const Cell& up, const Cell& down, const dataType* const scalars) const{
  return scalarMax<dataType>(up,scalars)-scalarMin<dataType>(down,scalars);
}

template <typename dataType>
bool DiscreteGradient::isHigherThan(const int vertexA,
    const int vertexB,
    const dataType* const scalars,
    const int* const offsets) const{
  if(scalars[vertexA] != scalars[vertexB]) return scalars[vertexA]>scalars[vertexB];
  else return offsets[vertexA]>offsets[vertexB];
}

template <typename dataType>
bool DiscreteGradient::isLowerThan(const int vertexA,
    const int vertexB,
    const dataType* const scalars,
    const int* const offsets) const{
  if(scalars[vertexA] != scalars[vertexB]) return scalars[vertexA]<scalars[vertexB];
  else return offsets[vertexA]<offsets[vertexB];
}

template <typename dataType>
int DiscreteGradient::cellMax(const int cellDim,
    const int cellA,
    const int cellB,
    const dataType* const scalars,
    const int* const offsets) const{
  const int vertexNumber=cellDim+1;

  if(dimensionality_==2){
    array<int,3> vsetA;
    array<int,3> vsetB;
    const auto sosGreaterThan=[&scalars,&offsets](const int a, const int b){
      if(scalars[a] != scalars[b]) return scalars[a]>scalars[b];
      else return offsets[a]>offsets[b];
    };

    switch(cellDim){
      case 0:
        return (isHigherThan<dataType>(cellA, cellB, scalars, offsets))? cellA : cellB;

      case 1:
        for(int k=0; k<2; ++k){
          inputTriangulation_->getEdgeVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getEdgeVertex(cellB, k, vsetB[k]);
        }
        break;

      case 2:
        for(int k=0; k<3; ++k){
          inputTriangulation_->getCellVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getCellVertex(cellB, k, vsetB[k]);
        }
        break;

      default: return -1;
    }

    sort(vsetA.begin(),vsetA.begin()+vertexNumber,sosGreaterThan);
    sort(vsetB.begin(),vsetB.begin()+vertexNumber,sosGreaterThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isHigherThan<dataType>(vsetA[k], vsetB[k], scalars,offsets))? cellA : cellB;
    }
  }
  else if(dimensionality_==3){
    array<int,4> vsetA;
    array<int,4> vsetB;
    const auto sosGreaterThan=[&scalars,&offsets](const int a, const int b){
      if(scalars[a] != scalars[b]) return scalars[a]>scalars[b];
      else return offsets[a]>offsets[b];
    };

    switch(cellDim){
      case 0:
        return (isHigherThan<dataType>(cellA, cellB, scalars, offsets)? cellA : cellB);

      case 1:
        for(int k=0; k<2; ++k){
          inputTriangulation_->getEdgeVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getEdgeVertex(cellB, k, vsetB[k]);
        }
        break;

      case 2:
        for(int k=0; k<3; ++k){
          inputTriangulation_->getTriangleVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getTriangleVertex(cellB, k, vsetB[k]);
        }
        break;

      case 3:
        for(int k=0; k<4; ++k){
          inputTriangulation_->getCellVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getCellVertex(cellB, k, vsetB[k]);
        }
        break;

      default: return -1;
    }

    sort(vsetA.begin(),vsetA.begin()+vertexNumber,sosGreaterThan);
    sort(vsetB.begin(),vsetB.begin()+vertexNumber,sosGreaterThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isHigherThan<dataType>(vsetA[k], vsetB[k], scalars,offsets)? cellA : cellB);
    }
  }

  return -1;
}

template <typename dataType>
int DiscreteGradient::cellMin(const int cellDim,
    const int cellA,
    const int cellB,
    const dataType* const scalars,
    const int* const offsets) const{
  const int vertexNumber=cellDim+1;

  if(dimensionality_==2){
    array<int,3> vsetA;
    array<int,3> vsetB;
    const auto sosLowerThan=[&scalars,&offsets](const int a, const int b){
      if(scalars[a] != scalars[b]) return scalars[a]<scalars[b];
      else return offsets[a]<offsets[b];
    };

    switch(cellDim){
      case 0:
        return (isLowerThan<dataType>(cellA, cellB, scalars, offsets))? cellA : cellB;

      case 1:
        for(int k=0; k<2; ++k){
          inputTriangulation_->getEdgeVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getEdgeVertex(cellB, k, vsetB[k]);
        }
        break;

      case 2:
        for(int k=0; k<3; ++k){
          inputTriangulation_->getCellVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getCellVertex(cellB, k, vsetB[k]);
        }
        break;

      default: return -1;
    }

    sort(vsetA.begin(),vsetA.begin()+vertexNumber,sosLowerThan);
    sort(vsetB.begin(),vsetB.begin()+vertexNumber,sosLowerThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isLowerThan<dataType>(vsetA[k], vsetB[k], scalars,offsets))? cellA : cellB;
    }
  }
  else if(dimensionality_==3){
    array<int,4> vsetA;
    array<int,4> vsetB;
    const auto sosLowerThan=[&scalars,&offsets](const int a, const int b){
      if(scalars[a] != scalars[b]) return scalars[a]<scalars[b];
      else return offsets[a]<offsets[b];
    };

    switch(cellDim){
      case 0:
        return (isLowerThan<dataType>(cellA, cellB, scalars, offsets))? cellA : cellB;

      case 1:
        for(int k=0; k<2; ++k){
          inputTriangulation_->getEdgeVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getEdgeVertex(cellB, k, vsetB[k]);
        }
        break;

      case 2:
        for(int k=0; k<3; ++k){
          inputTriangulation_->getTriangleVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getTriangleVertex(cellB, k, vsetB[k]);
        }
        break;

      case 3:
        for(int k=0; k<4; ++k){
          inputTriangulation_->getCellVertex(cellA, k, vsetA[k]);
          inputTriangulation_->getCellVertex(cellB, k, vsetB[k]);
        }
        break;

      default: return -1;
    }

    sort(vsetA.begin(),vsetA.begin()+vertexNumber,sosLowerThan);
    sort(vsetB.begin(),vsetB.begin()+vertexNumber,sosLowerThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isLowerThan<dataType>(vsetA[k], vsetB[k], scalars,offsets))? cellA : cellB;
    }
  }

  return -1;
}

template <typename dataType>
int DiscreteGradient::g0(const int cellDim,
    const int cellId,
    const dataType* const scalars,
    const int* const offsets) const{
  int facet0;
  int facet1;
  int facetMax{-1};

  if(dimensionality_==2){
    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facet0);
        inputTriangulation_->getEdgeVertex(cellId, 1, facet1);
        facetMax=cellMax<dataType>(0, facet0, facet1, scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getCellEdge(cellId,0,facet0);
        inputTriangulation_->getCellEdge(cellId,1,facet1);
        facetMax=cellMax<dataType>(1, facet0, facet1, scalars,offsets);

        inputTriangulation_->getCellEdge(cellId,2,facet0);
        facetMax=cellMax<dataType>(1, facet0, facetMax, scalars,offsets);
        break;

      default: return -1;
    }
  }
  else if(dimensionality_==3){
    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facet0);
        inputTriangulation_->getEdgeVertex(cellId, 1, facet1);
        facetMax=cellMax<dataType>(0, facet0, facet1, scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getTriangleEdge(cellId,0,facet0);
        inputTriangulation_->getTriangleEdge(cellId,1,facet1);
        facetMax=cellMax<dataType>(1, facet0, facet1, scalars,offsets);

        inputTriangulation_->getTriangleEdge(cellId,2,facet0);
        facetMax=cellMax<dataType>(1, facet0, facetMax, scalars,offsets);
        break;

      case 3:
        inputTriangulation_->getCellTriangle(cellId,0,facet0);
        inputTriangulation_->getCellTriangle(cellId,1,facet1);
        facetMax=cellMax<dataType>(2, facet0, facet1, scalars,offsets);

        inputTriangulation_->getCellTriangle(cellId,2,facet0);
        facetMax=cellMax<dataType>(2, facet0, facetMax, scalars,offsets);

        inputTriangulation_->getCellTriangle(cellId,3,facet0);
        facetMax=cellMax<dataType>(2, facet0, facetMax, scalars,offsets);
        break;

      default: return -1;
    }
  }

  return facetMax;
}

template <typename dataType>
int DiscreteGradient::g0_second(const int cellDim,
    const int cellId,
    const dataType* const scalars,
    const int* const offsets) const{
  int facetMax{-1};
  int facetMaxSecond{-1};

  if(dimensionality_==2){
    int facets[3];

    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facets[0]);
        inputTriangulation_->getEdgeVertex(cellId, 1, facets[1]);
        facetMaxSecond=cellMin<dataType>(0, facets[0], facets[1], scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getCellEdge(cellId,0,facets[0]);
        inputTriangulation_->getCellEdge(cellId,1,facets[1]);
        facetMax=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);

        inputTriangulation_->getCellEdge(cellId,2,facets[2]);
        facetMax=cellMax<dataType>(1, facets[2], facetMax, scalars,offsets);

        if(facetMax==facets[0])
          facetMaxSecond=cellMax<dataType>(1, facets[1], facets[2], scalars,offsets);
        else if(facetMax==facets[1])
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[2], scalars,offsets);
        else
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);
        break;
    }
  }
  else if(dimensionality_==3){
    int facets[4];

    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facets[0]);
        inputTriangulation_->getEdgeVertex(cellId, 1, facets[1]);
        facetMaxSecond=cellMin<dataType>(0, facets[0], facets[1], scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getTriangleEdge(cellId,0,facets[0]);
        inputTriangulation_->getTriangleEdge(cellId,1,facets[1]);
        facetMax=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);

        inputTriangulation_->getTriangleEdge(cellId,2,facets[2]);
        facetMax=cellMax<dataType>(1, facets[2], facetMax, scalars,offsets);

        if(facetMax==facets[0])
          facetMaxSecond=cellMax<dataType>(1, facets[1], facets[2], scalars,offsets);
        else if(facetMax==facets[1])
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[2], scalars,offsets);
        else
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);
        break;

      case 3:
        inputTriangulation_->getCellTriangle(cellId,0,facets[0]);
        inputTriangulation_->getCellTriangle(cellId,1,facets[1]);
        inputTriangulation_->getCellTriangle(cellId,2,facets[2]);
        inputTriangulation_->getCellTriangle(cellId,3,facets[3]);

        if(facets[0]==cellMax<dataType>(2,facets[0],facets[1],scalars,offsets)){
          facetMax=facets[0];
          facetMaxSecond=facets[1];
        }
        else{
          facetMax=facets[1];
          facetMaxSecond=facets[0];
        }

        for(int i=2; i<4; i++){
          if(facets[i]==cellMax<dataType>(2,facets[i],facetMax,scalars,offsets)){
            facetMaxSecond=facetMax;
            facetMax=facets[i];
          }
          else if(facets[i]==cellMax<dataType>(2,facets[i],facetMaxSecond,scalars,offsets)){
            facetMaxSecond=facets[i];
          }
        }
        break;
    }
  }

  return facetMaxSecond;
}

template <typename dataType>
int DiscreteGradient::assignGradient(const int alphaDim,
    const dataType* const scalars,
    const int* const offsets,
    vector<vector<int>>& gradient) const{
  const int betaDim=alphaDim+1;
  const int alphaNumber=gradient[alphaDim].size();

  if(dimensionality_==2){
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int alpha=0; alpha<alphaNumber; ++alpha){
      int betaNumber{};
      switch(alphaDim){
        case 0: betaNumber=inputTriangulation_->getVertexEdgeNumber(alpha); break;
        case 1: betaNumber=inputTriangulation_->getEdgeStarNumber(alpha); break;
      }
      int gamma{-1};
      for(int k=0; k<betaNumber; ++k){
        int beta;
        switch(alphaDim){
          case 0: inputTriangulation_->getVertexEdge(alpha,k,beta); break;
          case 1: inputTriangulation_->getEdgeStar(alpha,k,beta); break;
        }
        // take beta such that alpha is the highest facet of beta
        if(alpha==g0<dataType>(betaDim,beta,scalars,offsets)){
          if(gamma==-1)
            gamma=beta;
          else
            gamma=cellMin<dataType>(betaDim,beta,gamma,scalars,offsets);
        }
      }
      if(gamma!=-1){
        gradient[alphaDim][alpha]=gamma;
        gradient[betaDim][gamma]=alpha;
      }
    }
  }
  else if(dimensionality_==3){
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int alpha=0; alpha<alphaNumber; ++alpha){
      int betaNumber{};
      switch(alphaDim){
        case 0: betaNumber=inputTriangulation_->getVertexEdgeNumber(alpha); break;
        case 1: betaNumber=inputTriangulation_->getEdgeTriangleNumber(alpha); break;
        case 2: betaNumber=inputTriangulation_->getTriangleStarNumber(alpha); break;
      }
      int gamma{-1};
      for(int k=0; k<betaNumber; ++k){
        int beta;
        switch(alphaDim){
          case 0: inputTriangulation_->getVertexEdge(alpha,k,beta); break;
          case 1: inputTriangulation_->getEdgeTriangle(alpha,k,beta); break;
          case 2: inputTriangulation_->getTriangleStar(alpha,k,beta); break;
        }
        // take beta such that alpha is the highest facet of beta
        if(alpha==g0<dataType>(betaDim,beta,scalars,offsets)){
          if(gamma==-1)
            gamma=beta;
          else
            gamma=cellMin<dataType>(betaDim,beta,gamma,scalars,offsets);
        }
      }
      if(gamma!=-1){
        gradient[alphaDim][alpha]=gamma;
        gradient[betaDim][gamma]=alpha;
      }
    }
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::assignGradient2(const int alphaDim,
    const dataType* const scalars,
    const int* const offsets,
    vector<vector<int>>& gradient) const{
  if(alphaDim>0){
    const int betaDim=alphaDim+1;
    const int alphaNumber=gradient[alphaDim].size();

    if(dimensionality_==2){
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int alpha=0; alpha<alphaNumber; ++alpha){
        // alpha must be unpaired
        if(gradient[alphaDim][alpha]==-1){
          int betaNumber{};
          switch(alphaDim){
            case 1: betaNumber=inputTriangulation_->getEdgeStarNumber(alpha); break;
          }
          int gamma{-1};
          for(int k=0; k<betaNumber; ++k){
            int beta;
            switch(alphaDim){
              case 1: inputTriangulation_->getEdgeStar(alpha,k,beta); break;
            }
            // take beta such that alpha is the second highest facet of beta
            if(alpha==g0_second<dataType>(betaDim,beta,scalars,offsets)){
              if(gamma==-1)
                gamma=beta;
              else
                gamma=cellMin<dataType>(betaDim,beta,gamma,scalars,offsets);
            }
          }

          if(gamma!=-1 and gradient[betaDim][gamma]==-1){
            gradient[alphaDim][alpha]=gamma;
            gradient[betaDim][gamma]=alpha;
          }
        }
      }
    }
    else if(dimensionality_==3){
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int alpha=0; alpha<alphaNumber; ++alpha){
        // alpha must be unpaired
        if(gradient[alphaDim][alpha]==-1){
          int betaNumber{};
          switch(alphaDim){
            case 1: betaNumber=inputTriangulation_->getEdgeTriangleNumber(alpha); break;
            case 2: betaNumber=inputTriangulation_->getTriangleStarNumber(alpha); break;
          }
          int gamma{-1};
          for(int k=0; k<betaNumber; ++k){
            int beta;
            switch(alphaDim){
              case 1: inputTriangulation_->getEdgeTriangle(alpha,k,beta); break;
              case 2: inputTriangulation_->getTriangleStar(alpha,k,beta); break;
            }
            // take beta such that alpha is the second highest facet of beta
            if(alpha==g0_second<dataType>(betaDim,beta,scalars,offsets)){
              if(gamma==-1)
                gamma=beta;
              else
                gamma=cellMin<dataType>(betaDim,beta,gamma,scalars,offsets);
            }
          }

          if(gamma!=-1 and gradient[betaDim][gamma]==-1){
            gradient[alphaDim][alpha]=gamma;
            gradient[betaDim][gamma]=alpha;
          }
        }
      }
    }
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::buildGradient(){
  Timer t;

  const int* const offsets=static_cast<int*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfDimensions=getNumberOfDimensions();

  // init number of cells by dimension
  vector<int> numberOfCells(numberOfDimensions);
  for(int i=0; i<numberOfDimensions; ++i)
    numberOfCells[i]=getNumberOfCells(i);

  gradient_.clear();
  gradient_.resize(dimensionality_);
  for(int i=0; i<dimensionality_; ++i){
    // init gradient memory
    gradient_[i].resize(numberOfDimensions);
    gradient_[i][i].resize(numberOfCells[i], -1);
    gradient_[i][i+1].resize(numberOfCells[i+1], -1);

    // compute gradient pairs
    assignGradient<dataType>(i, scalars, offsets, gradient_[i]);
  }

  {
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

    stringstream msg;
    msg << "[DiscreteGradient] Data-set (" << numberOfVertices
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::buildGradient2(){
  Timer t;

  const int* const offsets=static_cast<int*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  for(int i=1; i<dimensionality_; ++i)
    assignGradient2<dataType>(i, scalars, offsets, gradient_[i]);

  {
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

    stringstream msg;
    msg << "[DiscreteGradient] Data-set (" << numberOfVertices
      << " points) processed twice in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints(const vector<Cell>& criticalPoints) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  vector<dataType>* outputCriticalPoints_points_cellScalars=
    static_cast<vector<dataType>*>(outputCriticalPoints_points_cellScalars_);

  (*outputCriticalPoints_numberOfPoints_)=0;

  const int numberOfDimensions=getNumberOfDimensions();
  vector<int> numberOfCriticalPointsByDimension(numberOfDimensions,0);

  // for all critical cells
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& cell=criticalPoints[i];
    const int cellDim=cell.dim_;
    const int cellId=cell.id_;
    numberOfCriticalPointsByDimension[cellDim]++;

    float incenter[3];
    getCellIncenter(cell, incenter);

    const dataType scalar=scalarMax<dataType>(cell, scalars);
    const char isOnBoundary=isBoundary(cell);

    outputCriticalPoints_points_->push_back(incenter[0]);
    outputCriticalPoints_points_->push_back(incenter[1]);
    outputCriticalPoints_points_->push_back(incenter[2]);

    outputCriticalPoints_points_cellDimensions_->push_back(cellDim);
    outputCriticalPoints_points_cellIds_->push_back(cellId);
    outputCriticalPoints_points_cellScalars->push_back(scalar);
    outputCriticalPoints_points_isOnBoundary_->push_back(isOnBoundary);

    (*outputCriticalPoints_numberOfPoints_)++;
  }

  {
    stringstream msg;
    for(int i=0; i<numberOfDimensions; ++i)
      msg << "[DiscreteGradient] " << numberOfCriticalPointsByDimension[i] << " " << i << "-cell(s)." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints() const{
  vector<Cell> criticalPoints;
  getCriticalPoints(criticalPoints);

  setCriticalPoints<dataType>(criticalPoints);

  return 0;
}

template <typename dataType>
int DiscreteGradient::getRemovableMaxima(const vector<pair<int,char>>& criticalPoints,
    vector<char>& isRemovable) const{
  const int* const offsets=static_cast<int*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfCriticalPoints=criticalPoints.size();
  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  const int maximumDim=dimensionality_;

  // Detect DMT-max cells to remove
  isRemovable.resize(numberOfCells);
  std::fill(isRemovable.begin(), isRemovable.end(), false);
  {
    vector<char> isAuthorized(numberOfCells,false);
    for(int i=0; i<numberOfCriticalPoints; ++i){
      const pair<int,char>& criticalPoint=criticalPoints[i];
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==maximumDim){
        const int starNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);

        // find maxStarId
        int maxStarId=-1;
        for(int j=0; j<starNumber; ++j){
          int starId;
          inputTriangulation_->getVertexStar(criticalPointId, j, starId);

          const Cell star(maximumDim,starId);
          if(isMaximum(star)){
            if(maxStarId==-1)
              maxStarId=starId;
            else
              maxStarId=cellMax<dataType>(maximumDim, maxStarId, starId, scalars, offsets);
          }
        }

        if(maxStarId==-1){
          if(!inputTriangulation_->isVertexOnBoundary(criticalPointId)){
            stringstream msg;
            msg << "[DiscreteGradient] No DMT-maxima connected to non-boundary PL-maximum vertexId=" << criticalPointId << endl;
            dMsg(cout, msg.str(), advancedInfoMsg);
          }

          continue;
        }
        else
          isAuthorized[maxStarId]=true;
      }
    }

    for(int i=0; i<numberOfCells; ++i){
      const Cell cell(maximumDim,i);
      if(isMaximum(cell) and !isAuthorized[i])
        isRemovable[i]=true;
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::getRemovableSaddles1(const vector<pair<int,char>>& criticalPoints,
    vector<char>& isRemovable) const{
  const int* const offsets=static_cast<int*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfCriticalPoints=criticalPoints.size();
  const int numberOfSaddles=inputTriangulation_->getNumberOfEdges();
  const char saddleDim=1;

  // Detect 1-saddles to remove
  isRemovable.resize(numberOfSaddles);
  std::fill(isRemovable.begin(), isRemovable.end(), false);
  {
    vector<char> isAuthorized(numberOfSaddles,false);
    for(int i=0; i<numberOfCriticalPoints; ++i){
      const pair<int,char>& criticalPoint=criticalPoints[i];
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==saddleDim){
        const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(criticalPointId);

        // find maxEdgeId
        int maxEdgeId=-1;
        for(int j=0; j<edgeNumber; ++j){
          int edgeId;
          inputTriangulation_->getVertexEdge(criticalPointId, j, edgeId);

          const Cell edge(saddleDim,edgeId);
          if(isSaddle1(edge)){
            if(maxEdgeId==-1)
              maxEdgeId=edgeId;
            else
              maxEdgeId=cellMax<dataType>(saddleDim, maxEdgeId, edgeId, scalars, offsets);
          }
        }

        if(maxEdgeId==-1){
          if(!inputTriangulation_->isVertexOnBoundary(criticalPointId)){
            stringstream msg;
            msg << "[DiscreteGradient] No DMT-1saddle connected to non-boundary PL-1saddle vertexId=" << criticalPointId << endl;
            dMsg(cout, msg.str(), advancedInfoMsg);
          }

          continue;
        }
        else
          isAuthorized[maxEdgeId]=true;
      }
    }

    if(CollectPersistencePairs){
      for(int i=0; i<numberOfSaddles; ++i){
        const Cell cell(saddleDim,i);
        if(isSaddle1(cell) and isAuthorized[i])
          isRemovable[i]=true;
      }
    }
    else{
      for(int i=0; i<numberOfSaddles; ++i){
        const Cell cell(saddleDim,i);
        if(isSaddle1(cell) and !isAuthorized[i])
          isRemovable[i]=true;
      }
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::getRemovableSaddles2(const vector<pair<int,char>>& criticalPoints,
    vector<char>& isRemovable) const{
  const int* const offsets=static_cast<int*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

  const int numberOfCriticalPoints=criticalPoints.size();
  const int numberOfSaddleCandidates=inputTriangulation_->getNumberOfTriangles();
  const char saddleDim=2;

  // Detect 1-saddles to remove
  isRemovable.resize(numberOfSaddleCandidates);
  std::fill(isRemovable.begin(), isRemovable.end(), false);
  {
    vector<char> isPLSaddle2(numberOfVertices, false);
    for(int i=0; i<numberOfCriticalPoints; ++i){
      const pair<int,char>& criticalPoint=criticalPoints[i];
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==saddleDim)
        isPLSaddle2[criticalPointId]=true;
    }

    vector<int> maxTriangle(numberOfVertices, -1);
    for(int i=0; i<numberOfSaddleCandidates; ++i){
      const Cell triangle(2,i);
      if(isSaddle2(triangle)){
        for(int j=0; j<3; ++j){
          int vertexId;
          inputTriangulation_->getTriangleVertex(i, j, vertexId);

          const int maxTriangleId=maxTriangle[vertexId];

          if(isPLSaddle2[vertexId]){
            if(maxTriangleId==-1)
              maxTriangle[vertexId]=i;
            else
              maxTriangle[vertexId]=cellMax<dataType>(saddleDim, maxTriangleId, i, scalars, offsets);
          }
        }
      }
    }

    vector<char> isAuthorized(numberOfSaddleCandidates, false);
    for(int i=0; i<numberOfVertices; ++i){
      const int maxTriangleId=maxTriangle[i];

      if(maxTriangleId!=-1)
        isAuthorized[maxTriangleId]=true;
    }

    if(CollectPersistencePairs){
      for(int i=0; i<numberOfSaddleCandidates; ++i){
        const Cell cell(saddleDim,i);
        if(isSaddle2(cell) and isAuthorized[i])
          isRemovable[i]=true;
      }
    }
    else{
      for(int i=0; i<numberOfSaddleCandidates; ++i){
        const Cell cell(saddleDim,i);
        if(isSaddle2(cell) and !isAuthorized[i])
          isRemovable[i]=true;
      }
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleMaximumConnections(const vector<char>& isRemovableMaximum,
    vector<char>& isRemovableSaddle,
    const bool allowBoundary,
    vector<Segment>& segments,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints) const{

  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int maximumDim=dimensionality_;
  const int saddleDim=maximumDim-1;

  // Part 1 : build initial structures
  // add the saddles to CriticalPointList and count them
  const int numberOfSaddleCandidates=getNumberOfCells(saddleDim);

  for(int i=0; i<numberOfSaddleCandidates; ++i){
    if(isRemovableSaddle[i]){
      const Cell saddleCandidate(saddleDim, i);

      if(!allowBoundary and isBoundary(saddleCandidate))
        continue;

      if(isCellCritical(saddleCandidate))
        criticalPoints.push_back(CriticalPoint(saddleCandidate));
    }
  }
  const int numberOfSaddles=criticalPoints.size();

  // add the maxima to CriticalPointList and build MaxIndex
  const int numberOfMaximumCandidates=getNumberOfCells(maximumDim);

  // experimental
  vector<char> mask(numberOfMaximumCandidates, false);

  vector<int> maximumIndex(numberOfMaximumCandidates,-1);
  for(int i=0; i<numberOfMaximumCandidates; ++i){
    if(isRemovableMaximum[i]){
      const int index=criticalPoints.size();
      maximumIndex[i]=index;

      const Cell maximum(maximumDim, i);
      criticalPoints.push_back(CriticalPoint(maximum));
      mask[i]=true;
    }
  }

  const int numberOfVPaths=2*numberOfSaddles;
  vpaths.resize(numberOfVPaths);
  segments.resize(numberOfVPaths);

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int sourceIndex=i;
    CriticalPoint& source=criticalPoints[sourceIndex];

    const Cell& saddle=source.cell_;
    const int saddleId=saddle.id_;

    int starNumber{};
    if(maximumDim==2)
      starNumber=inputTriangulation_->getEdgeStarNumber(saddleId);
    else if(maximumDim==3)
      starNumber=inputTriangulation_->getTriangleStarNumber(saddleId);

    vector<vector<Cell>> paths(starNumber);
    for(int j=0; j<starNumber; ++j){
      int starId;
      if(maximumDim==2)
        inputTriangulation_->getEdgeStar(saddleId, j, starId);
      else if(maximumDim==3)
        inputTriangulation_->getTriangleStar(saddleId, j, starId);

      const Cell star(maximumDim, starId);

      paths[j].push_back(saddle);
      getAscendingPath(star, paths[j]);
    }

    // detect initial double-connection
    if(starNumber>1){
      bool isDoubleConnected=false;
      const Cell& lastCell0=paths[0].back();
      for(int j=1; j<starNumber; ++j){
        const Cell& lastCell=paths[j].back();

        if(lastCell0.id_==lastCell.id_){
          isDoubleConnected=true;
          break;
        }
      }
      if(isDoubleConnected)
        continue;
    }

    for(int j=0; j<starNumber; ++j){
      const int shift=j;

      // apriori: there is at least 1 one cell
      const Cell& lastCell=paths[j].back();
      if(isMaximum(lastCell) and isRemovableMaximum[lastCell.id_]){
        const Cell& maximum=lastCell;

        const int destinationIndex=maximumIndex[maximum.id_];
        CriticalPoint& destination=criticalPoints[destinationIndex];

        mask[maximum.id_]=false;

        // update source and destination
        const int sourceSlot=source.omp_addSlot();
        const int destinationSlot=destination.omp_addSlot();

        // update vpath
        const int vpathIndex=2*sourceIndex+shift;
        VPath& vpath=vpaths[vpathIndex];
        vpath.source_=sourceIndex;
        vpath.destination_=destinationIndex;
        vpath.sourceSlot_=sourceSlot;
        vpath.destinationSlot_=destinationSlot;
        vpath.states_.push_back(1);
        vpath.segments_.push_back(vpathIndex);
        vpath.persistence_=getPersistence<dataType>(maximum, saddle, scalars);
        vpath.isValid_=true;

        // update segment
        Segment& segment=segments[vpathIndex];
        segment.orientation_=true;
        segment.cells_=std::move(paths[j]);
        segment.isValid_=true;
      }
    }
  }

  cout << "MSC max not covered : " << std::count(mask.begin(), mask.end(), true) << endl;

  // experimental
  if(dimensionality_==3){
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    vector<set<int>> vtriangles(numberOfVertices);

    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    for(int i=0; i<numberOfTriangles; ++i){
      for(int j=0; j<3; ++j){
        int vertexId;
        inputTriangulation_->getTriangleVertex(i, j, vertexId);
        vtriangles[vertexId].insert(i);
      }
    }

    vector<int> saddles2;
    for(pair<int,char> criticalPoint : criticalPoints_){
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==2)
        saddles2.push_back(criticalPointId);
    }

    int n=0;
    for(int saddle2Id : saddles2){
      bool isFound=false;
      for(int triangleId : vtriangles[saddle2Id]){
        const Cell saddle(2,triangleId);
        const int saddleId=triangleId;

        if(!isRemovableSaddle[triangleId] and isCellCritical(saddle)){
          const int starNumber=inputTriangulation_->getTriangleStarNumber(saddleId);

          vector<vector<Cell>> paths(starNumber);
          for(int j=0; j<starNumber; ++j){
            int starId;
            inputTriangulation_->getTriangleStar(saddleId, j, starId);

            const Cell star(maximumDim, starId);

            paths[j].push_back(saddle);
            getAscendingPath(star, paths[j]);
          }

          // detect initial double-connection
          if(starNumber>1){
            bool isDoubleConnected=false;
            const Cell& lastCell0=paths[0].back();
            for(int j=1; j<starNumber; ++j){
              const Cell& lastCell=paths[j].back();

              if(lastCell0.id_==lastCell.id_){
                isDoubleConnected=true;
                break;
              }
            }
            if(isDoubleConnected)
              continue;
          }

          for(int j=0; j<starNumber; ++j){
            const Cell& lastCell=paths[j].back();
            // apriori: triangleId is connected by the gradient to a removable max
            if(isMaximum(lastCell) and isRemovableMaximum[lastCell.id_] and mask[lastCell.id_]){
              for(int tId : vtriangles[saddle2Id]){
                // apriori: triangleId has a removable neighbor
                if(tId!=triangleId and isRemovableSaddle[tId]){
                  isFound=true;
                  break;
                }
              }
            }
          }

        }
      }
      if(isFound) ++n;
    }
    cout << "number of PL 2-saddle connected to a non-removable MSC 2-saddle connected to a removable not-covered MSC maximum : " << n << endl;
  }
  if(dimensionality_==3){
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    vector<set<int>> vtriangles(numberOfVertices);

    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    for(int i=0; i<numberOfTriangles; ++i){
      for(int j=0; j<3; ++j){
        int vertexId;
        inputTriangulation_->getTriangleVertex(i, j, vertexId);
        vtriangles[vertexId].insert(i);
      }
    }

    vector<int> saddles2;
    for(pair<int,char> criticalPoint : criticalPoints_){
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==2)
        saddles2.push_back(criticalPointId);
    }

    int n=0;
    for(int saddle2Id : saddles2){
      bool isFound=false;
      for(int triangleId : vtriangles[saddle2Id]){
        const Cell saddle(2,triangleId);
        const int saddleId=triangleId;

        if(!isRemovableSaddle[triangleId] and isCellCritical(saddle)){
          const int starNumber=inputTriangulation_->getTriangleStarNumber(saddleId);

          vector<vector<Cell>> paths(starNumber);
          for(int j=0; j<starNumber; ++j){
            int starId;
            inputTriangulation_->getTriangleStar(saddleId, j, starId);

            const Cell star(maximumDim, starId);

            paths[j].push_back(saddle);
            getAscendingPath(star, paths[j]);
          }

          // detect initial double-connection
          if(starNumber>1){
            bool isDoubleConnected=false;
            const Cell& lastCell0=paths[0].back();
            for(int j=1; j<starNumber; ++j){
              const Cell& lastCell=paths[j].back();

              if(lastCell0.id_==lastCell.id_){
                isDoubleConnected=true;
                break;
              }
            }
            if(isDoubleConnected)
              continue;
          }

          for(int j=0; j<starNumber; ++j){
            const Cell& lastCell=paths[j].back();
            // apriori: triangleId is connected by the gradient to a removable max
            if(isMaximum(lastCell) and isRemovableMaximum[lastCell.id_] and mask[lastCell.id_]){
              for(int tId : vtriangles[saddle2Id]){
                // apriori: triangleId has a removable neighbor
                if(tId!=triangleId and isRemovableSaddle[tId]){
                  isFound=true;
                  break;
                }
              }
            }
            else
              isFound=false;
          }

        }
      }
      if(isFound) ++n;
    }
    cout << "number of PL 2-saddle connected to a non-removable MSC 2-saddle connected only to removable not-covered MSC maximum : " << n << endl;
  }
  if(dimensionality_==3){
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    vector<set<int>> vtriangles(numberOfVertices);

    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    for(int i=0; i<numberOfTriangles; ++i){
      for(int j=0; j<3; ++j){
        int vertexId;
        inputTriangulation_->getTriangleVertex(i, j, vertexId);
        vtriangles[vertexId].insert(i);
      }
    }

    vector<int> saddles2;
    for(pair<int,char> criticalPoint : criticalPoints_){
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(criticalPointType==2)
        saddles2.push_back(criticalPointId);
    }

    int n=0;
    for(int saddle2Id : saddles2){
      bool isFound=false;
      for(int triangleId : vtriangles[saddle2Id]){
        const Cell saddle(2,triangleId);
        const int saddleId=triangleId;

        if(isRemovableSaddle[triangleId]){
          const int starNumber=inputTriangulation_->getTriangleStarNumber(saddleId);

          vector<vector<Cell>> paths(starNumber);
          for(int j=0; j<starNumber; ++j){
            int starId;
            inputTriangulation_->getTriangleStar(saddleId, j, starId);

            const Cell star(maximumDim, starId);

            paths[j].push_back(saddle);
            getAscendingPath(star, paths[j]);
          }

          // detect initial double-connection
          if(starNumber>1){
            bool isDoubleConnected=false;
            const Cell& lastCell0=paths[0].back();
            for(int j=1; j<starNumber; ++j){
              const Cell& lastCell=paths[j].back();

              if(lastCell0.id_==lastCell.id_){
                isDoubleConnected=true;
                break;
              }
            }
            if(isDoubleConnected)
              continue;
          }

          for(int j=0; j<starNumber; ++j){
            const Cell& lastCell=paths[j].back();
            // apriori: triangleId is connected by the gradient to a removable max
            if(isMaximum(lastCell) and !isRemovableMaximum[lastCell.id_]){
              for(int tId : vtriangles[saddle2Id]){
                // apriori: triangleId has a removable neighbor
                if(tId!=triangleId and isRemovableSaddle[tId]){
                  isFound=true;
                  break;
                }
              }
            }
            else
              isFound=false;
          }
        }
      }
      if(isFound) ++n;
    }
    cout << "number of PL 2-saddle connected to at least one removable MSC 2-saddle connected only to non-removable MSC maximum : " << n << endl;
  }
  //

  // Part 3 : initialize the last structures
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    CriticalPoint& cp=criticalPoints[i];

    const int numberOfSlots=cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_=0;
  }

#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const int sourceIndex=vpath.source_;
      const int destinationIndex=vpath.destination_;

      const int sourceSlot=vpath.sourceSlot_;
      const int destinationSlot=vpath.destinationSlot_;

      CriticalPoint& source=criticalPoints[sourceIndex];
      CriticalPoint& destination=criticalPoints[destinationIndex];

      source.vpaths_[sourceSlot]=i;
      destination.vpaths_[destinationSlot]=i;
    }
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleMaximumConnections(const vector<VPath>& vpaths,
    set<pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>& S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_)
      S.insert(make_pair(vpath.persistence_,i));
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::computeCoefficients(const bool isDense,
    vector<char>& denseCoefficients,
    vector<Segment>& segments,
    const CriticalPoint& source,
    VPath& newVPath,
    const vector<VPath>& vpaths) const{
  if(isDense){
    const int numberOfSegments=segments.size();
    // apriori : the following will make only one allocation, the size is fixed
    denseCoefficients.resize(numberOfSegments);

    std::fill(denseCoefficients.begin(), denseCoefficients.end(), 0);

    // 1) initialize accumulator
    const int numberOfNewVPathSegments=newVPath.segments_.size();
    for(int i=0; i<numberOfNewVPathSegments; ++i){
      const int segmentId=newVPath.segments_[i];
      const char segmentState=newVPath.states_[i];

      denseCoefficients[segmentId]=segmentState;
    }

    // 2) add source.vpaths.segments to accumulator
    const int numberOfSourceVPaths=source.vpaths_.size();
    for(int i=0; i<numberOfSourceVPaths; ++i){
      const int sourceVPathId=source.vpaths_[i];
      const VPath& sourceVPath=vpaths[sourceVPathId];

      if(sourceVPath.isValid_){
        const int numberOfSourceVPathSegments=sourceVPath.segments_.size();
        for(int j=0; j<numberOfSourceVPathSegments; ++j){
          const int segmentId=sourceVPath.segments_[j];
          const char segmentState=sourceVPath.states_[j];

          denseCoefficients[segmentId]+=segmentState;
        }
      }
    }

    // 3) update newVPath to the result of accumulation
    newVPath.states_.clear();
    newVPath.segments_.clear();
    for(int i=0; i<numberOfSegments; ++i){
      const int segmentId=i;
      const char segmentState=denseCoefficients[segmentId];

      if(segmentState!=0){
        newVPath.states_.push_back(segmentState);
        newVPath.segments_.push_back(segmentId);
      }
    }
  }
  else{
    vector<pair<int,char>> sparseCoefficients;

    // 1) initialize accumulator
    const int numberOfNewVPathSegments=newVPath.segments_.size();
    for(int i=0; i<numberOfNewVPathSegments; ++i){
      const int segmentId=newVPath.segments_[i];
      const char segmentState=newVPath.states_[i];

      sparseCoefficients.push_back(make_pair(segmentId,segmentState));
    }

    // 2) add source.vpaths.segments to accumulator
    const int numberOfSourceVPaths=source.vpaths_.size();
    for(int i=0; i<numberOfSourceVPaths; ++i){
      const int sourceVPathId=source.vpaths_[i];
      const VPath& sourceVPath=vpaths[sourceVPathId];

      if(sourceVPath.isValid_){
        const int numberOfSourceVPathSegments=sourceVPath.segments_.size();
        for(int j=0; j<numberOfSourceVPathSegments; ++j){
          const int segmentId=sourceVPath.segments_[j];
          const char segmentState=sourceVPath.states_[j];

          bool isIn=false;
          const int sparseCoefficientsSize=sparseCoefficients.size();
          for(int k=0; k<sparseCoefficientsSize; ++k){
            const int savedSegmentId=sparseCoefficients[k].first;
            const char savedSegmentState=sparseCoefficients[k].second;

            if(segmentId==savedSegmentId){
              sparseCoefficients[k].second=segmentState+savedSegmentState;
              isIn=true;
            }
          }
          if(!isIn)
            sparseCoefficients.push_back(make_pair(segmentId,segmentState));
        }
      }
    }

    // 3) update newVPath to the result of accumulation
    newVPath.states_.clear();
    newVPath.segments_.clear();
    const int sparseCoefficientsSize=sparseCoefficients.size();
    for(int i=0; i<sparseCoefficientsSize; ++i){
      const int segmentId=sparseCoefficients[i].first;
      const char segmentState=sparseCoefficients[i].second;

      // apriori : sparseCoefficients store coefficient zero; we must remove them
      if(segmentState!=0){
        newVPath.states_.push_back(segmentState);
        newVPath.segments_.push_back(segmentId);
      }
    }
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleMaximumConnections(const int iterationThreshold,
    set<pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>& S,
    vector<Segment>& segments,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  int numberOfIterations{};
  vector<char> denseCoefficients;
  while(S.size()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=ptr->second;
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    if(vpath.isValid_){
#ifdef PRINT_DEBUG
      {
        const int sourceId=vpath.source_;
        const int destinationId=vpath.destination_;
        cout << "reverse saddle2.id=" << criticalPoints[sourceId].cell_.id_;
        cout << " max.id=" << criticalPoints[destinationId].cell_.id_ << endl;
      }
#endif

      // all segments of the selected vpath are reversed
      const int numberOfVPathSegments=vpath.segments_.size();
      for(int i=0; i<numberOfVPathSegments; ++i){
        const int segmentId=vpath.segments_[i];
        Segment& segment=segments[segmentId];

        segment.orientation_=!segment.orientation_;
        vpath.states_[i]*=-1;
      }

      // search new destination for newVPath
      int newDestinationId=-1;
      const int sourceId=vpath.source_;
      const int destinationId=vpath.destination_;
      CriticalPoint& source=criticalPoints[sourceId];
      CriticalPoint& destination=criticalPoints[destinationId];
      const int numberOfSourceVPaths=source.vpaths_.size();
      const int numberOfDestinationVPaths=destination.vpaths_.size();
      for(int i=0; i<numberOfSourceVPaths; ++i){
        const int sourceVPathId=source.vpaths_[i];
        const VPath& sourceVPath=vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_!=destinationId){
          newDestinationId=sourceVPath.destination_;
          break;
        }
      }

      // no valid destination so continue
      const bool hasInvalidDestination=(newDestinationId==-1);
      if(hasInvalidDestination)
        vpath.invalidate();

      // update destination.vpaths
      for(int i=0; i<numberOfDestinationVPaths; ++i){
        // newVPath = destination.vpath
        const int newVPathId=destination.vpaths_[i];
        VPath& newVPath=vpaths[newVPathId];

        if(newVPathId==vpathId) continue;

        if(!newVPath.isValid_) continue;

        if(hasInvalidDestination){
          newVPath.invalidate();
          continue;
        }

        // check for double-connections in newVPath
        const int newSourceId=newVPath.source_;
        CriticalPoint& newSource=criticalPoints[newSourceId];
        bool isDoubleConnected=false;
        const int numberOfNewSourceVPaths=newSource.vpaths_.size();
        for(int j=0; j<numberOfNewSourceVPaths; ++j){
          const int newSourceVPathId=newSource.vpaths_[j];
          VPath& newSourceVPath=vpaths[newSourceVPathId];

          if(newSourceVPath.isValid_ and newSourceVPath.destination_==newDestinationId){
            isDoubleConnected=true;
            newSourceVPath.invalidate();
            break;
          }
        }

        // invalid newVPath
        if(isDoubleConnected){
          newVPath.invalidate();
          continue;
        }

        // compute final coefficients of newVPath with sparse representation
        computeCoefficients<dataType>(false, denseCoefficients, segments, source, newVPath, vpaths);

        // update the destination of newVPath
        newVPath.destination_=newDestinationId;

        // add newVPath to newDestination connectivity
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        newDestination.vpaths_.push_back(newVPathId);

        // erase newVPath
        S.erase(make_pair(newVPath.persistence_,newVPathId));

        // update persistence
        newVPath.persistence_=getPersistence<dataType>(newDestination.cell_,newSource.cell_, scalars);

        // repush newVPath to confirm update
        S.insert(make_pair(newVPath.persistence_,newVPathId));
      }

      // invalid source.vpaths
      for(int i=0; i<numberOfSourceVPaths; ++i){
        const int sourceVPathId=source.vpaths_[i];
        VPath& sourceVPath=vpaths[sourceVPathId];

        sourceVPath.invalidate();
      }

      // erase connectivity of source and destination
      source.vpaths_.clear();
      destination.vpaths_.clear();
    }

    ++numberOfIterations;
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::reverseSaddleMaximumConnections(const vector<Segment>& segments){
  Timer t;

  const int numberOfSegments=segments.size();

  for(int i=0; i<numberOfSegments; ++i){
    const Segment& segment=segments[i];
    if(segment.isValid_ and segment.orientation_==false)
      reverseAscendingPath(segment.cells_);
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Gradient reversal step :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleMaximumConnections(const vector<char>& isRemovableMaximum,
    vector<char>& isRemovableSaddle,
    const int iterationThreshold,
    const bool allowBoundary){
  Timer t;

  // Part 1 : initialization
  vector<Segment> segments;
  vector<VPath> vpaths;
  vector<CriticalPoint> criticalPoints;
  initializeSaddleMaximumConnections<dataType>(isRemovableMaximum, isRemovableSaddle, allowBoundary, segments, vpaths, criticalPoints);

  // Part 2 : push the vpaths and order by persistence
  SaddleMaximumVPathComparator<dataType> cmp_f;
  set<pair<dataType,int>, SaddleMaximumVPathComparator<dataType>> S(cmp_f);
  orderSaddleMaximumConnections<dataType>(vpaths, S);

  // Part 3 : process the vpaths
  processSaddleMaximumConnections<dataType>(iterationThreshold, S, segments, vpaths, criticalPoints);

  // Part 4 : gradient reversal
  reverseSaddleMaximumConnections<dataType>(segments);

  {
    stringstream msg;
    msg << "[DiscreteGradient] Saddle-Maximum pairs ";
    if(allowBoundary)
      msg <<  "on boundary ";
    msg << "simplified in "
      << t.getElapsedTime() << " s, "<< threadNumber_ << " thread(s)." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::naive_reverseSaddleSaddleConnection(vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2){
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  wallId_t descendingWallId=1;
  vector<wallId_t> isVisited(numberOfTriangles, 0);

  bool isFirst=true;
  Cell minSaddle1;
  Cell minSaddle2;
  vector<Cell> minVPath;
  dataType minPersistence=-1;

  for(int i=0; i<numberOfTriangles; ++i){
    // process removable non-boundary 2-saddle
    if(isRemovableSaddle2[i] and !inputTriangulation_->isTriangleOnBoundary(i)){
      const Cell saddle2(2,i);

      set<int> saddles1;
      const wallId_t savedDescendingWallId=descendingWallId;
      getDescendingWall(descendingWallId, saddle2, isVisited, nullptr, &saddles1);
      ++descendingWallId;

      if(saddles1.empty()){
#ifdef PRINT_INFOS
        cout << "[DiscreteGradient] No 1-saddle on the wall of 2-saddle.id_=" << saddle2.id_ << endl;
#endif
        continue;
      }

      bool isPairable=false;
      const int numberOfSaddles1=saddles1.size();
      vector<vector<Cell>> vpaths(numberOfSaddles1);
      vector<char> isMultiConnected(numberOfSaddles1, false);
      int j=0;
      for(const int saddle1Id : saddles1){
        const Cell saddle1(1,saddle1Id);

        isMultiConnected[j]=getAscendingPathThroughWall(savedDescendingWallId, saddle1, saddle2, isVisited, &(vpaths[j]));

        if(!isMultiConnected[j]){
          const Cell& lastCell=vpaths[j].back();
          if(lastCell.dim_!=saddle2.dim_ or lastCell.id_!=saddle2.id_){
            ++j;
            continue;
          }
        }

        if(!isMultiConnected[j]){
          // pairing with a removable 1-saddle if possible
          if(isRemovableSaddle1[saddle1.id_]){
            isPairable=true;

            const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
            if(isFirst or persistence<minPersistence){
              // update minimizer
              isFirst=false;
              minSaddle1=saddle1;
              minSaddle2=saddle2;
              minVPath=vpaths[j];
              minPersistence=persistence;
            }
          }
        }
        ++j;
      }

      if(!isPairable and AllowReversingWithNonRemovable){
        // pairing with a non-removable 1-saddle anyway
        j=0;
        for(const int saddle1Id : saddles1){
          const Cell saddle1(1,saddle1Id);

          if(!isMultiConnected[j]){
            isPairable=true;

            const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
            if(isFirst or persistence<minPersistence){
              // update minimizer
              isFirst=false;
              minSaddle1=saddle1;
              minSaddle2=saddle2;
              minVPath=vpaths[j];
              minPersistence=persistence;
            }
          }
          ++j;
        }
      }

      if(!isPairable){
#ifdef PRINT_INFOS
        cout << "[DiscreteGradient] Non-pairable 2-saddle.id_=" << saddle2.id_ << endl;
#endif
      }
    }
  }

  if(!isFirst){
#ifdef PRINT_DEBUG
    for(const Cell& cell : minVPath)
      cout << cell.id_ << " ";
    cout << " persistence=" << minPersistence << endl;
#endif

    if(CollectPersistencePairs and outputPersistencePairs_)
      outputPersistencePairs_->push_back(make_tuple(minSaddle1, minSaddle2));

    const Cell& lastCell=minVPath.back();
    if(lastCell.dim_==minSaddle2.dim_ and lastCell.id_==minSaddle2.id_)
      reverseAscendingPathOnWall(minVPath);
    else{
#ifdef PRINT_ERROR
      cout << "[DiscreteGradient] Error : lastCell in ascending path is not saddle2.id=" << minSaddle2.id_ << endl;
#endif
    }
  }

  if(isFirst) return -1;

  return 0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections1(const vector<char>& isRemovableSaddle1,
    const vector<char>& isRemovableSaddle2,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    vector<int>& saddle1Index,
    vector<int>& saddle2Index) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int maximumDim=dimensionality_;
  const int saddle2Dim=maximumDim-1;
  const int saddle1Dim=saddle2Dim-1;

  // Part 1 : build initial structures
  // add the 2-saddles to CriticalPointList
  const int numberOfSaddle2Candidates=getNumberOfCells(saddle2Dim);
  saddle2Index.resize(numberOfSaddle2Candidates, -1);
  for(int i=0; i<numberOfSaddle2Candidates; ++i){
    //if(MustBeRemovable and !isRemovableSaddle2[i]) continue;
    if(!isRemovableSaddle2[i]) continue;

    const Cell saddle2Candidate(saddle2Dim, i);

    if(!AllowBoundary and isBoundary(saddle2Candidate)) continue;

    if(isSaddle2(saddle2Candidate)){
      const int index=criticalPoints.size();
      saddle2Index[i]=index;
      criticalPoints.push_back(CriticalPoint(saddle2Candidate));
    }
  }
  const int numberOf2Saddles=criticalPoints.size();

  // add the 1-saddles to CriticalPointList
  const int numberOfSaddle1Candidates=getNumberOfCells(saddle1Dim);
  saddle1Index.resize(numberOfSaddle1Candidates, -1);
  for(int i=0; i<numberOfSaddle1Candidates; ++i){
    if(MustBeRemovable and !isRemovableSaddle1[i]) continue;

    const Cell saddle1Candidate(saddle1Dim, i);

    if(ForceBoundary and !isBoundary(saddle1Candidate)) continue;

    if(isSaddle1(saddle1Candidate)){
      const int index=criticalPoints.size();
      saddle1Index[i]=index;
      criticalPoints.push_back(CriticalPoint(saddle1Candidate));
    }
  }

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
  wallId_t descendingWallId=1;
  vector<wallId_t> isVisited(numberOfSaddle2Candidates, 0);
  for(int i=0; i<numberOf2Saddles; ++i){
    const int destinationIndex=i;
    CriticalPoint& destination=criticalPoints[destinationIndex];
    const Cell& saddle2=destination.cell_;

    set<int> saddles1;
    const wallId_t savedDescendingWallId=descendingWallId;
    getDescendingWall(descendingWallId, saddle2, isVisited, nullptr, &saddles1);
    ++descendingWallId;

    for(const int saddle1Id : saddles1){
      if(MustBeRemovable and !isRemovableSaddle1[saddle1Id]) continue;

      const Cell& saddle1=Cell(1,saddle1Id);

      if(ForceBoundary and !isBoundary(saddle1)) continue;

      vector<Cell> path;
      const bool isMultiConnected=getAscendingPathThroughWall(savedDescendingWallId, saddle1, saddle2, isVisited, &path);

      if(!isMultiConnected){
        const int sourceIndex=saddle1Index[saddle1Id];
        CriticalPoint& source=criticalPoints[sourceIndex];

        // update source and destination
        const int sourceSlot=source.addSlot();
        const int destinationSlot=destination.addSlot();

        // update vpath
        const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
        vpaths.push_back(VPath(true,-1,sourceIndex,destinationIndex,sourceSlot,destinationSlot,persistence));
      }
    }
  }

  // Part 3 : initialize the last structures
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    CriticalPoint& cp=criticalPoints[i];

    const int numberOfSlots=cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_=0;
  }

  const int numberOfVPaths=vpaths.size();
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const int sourceIndex=vpath.source_;
      const int destinationIndex=vpath.destination_;

      const int sourceSlot=vpath.sourceSlot_;
      const int destinationSlot=vpath.destinationSlot_;

      CriticalPoint& source=criticalPoints[sourceIndex];
      CriticalPoint& destination=criticalPoints[destinationIndex];

      source.vpaths_[sourceSlot]=i;
      destination.vpaths_[destinationSlot]=i;
    }
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections1(const vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const int saddleId=criticalPoints[vpath.destination_].cell_.id_;
      S.insert(make_tuple(vpath.persistence_,i,saddleId));
    }
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleSaddleConnections1(const int iterationThreshold,
    set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S,
    vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    vector<int>& saddle1Index,
    vector<int>& saddle2Index){
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  const int optimizedSize=std::max(numberOfEdges, numberOfTriangles);
  wallId_t wallId=1;
  vector<wallId_t> isVisited(optimizedSize, 0);

  int numberOfIterations{};
  while(!S.empty()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=get<1>(*ptr);
#ifdef PRINT_DEBUG
    const dataType vpathPersistence=get<0>(*ptr);
#endif
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    if(vpath.isValid_){
      {
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;

        set<int> saddles1;
        const wallId_t savedWallId=wallId;
        getDescendingWall(wallId, minSaddle2, isVisited, nullptr, &saddles1);
        ++wallId;

        // check if at least one connection exists
        auto isFound=saddles1.find(minSaddle1.id_);
        if(isFound==saddles1.end()){
          ++numberOfIterations;
          continue;
        }

        // check if there is multiple connections
        vector<Cell> path;
        const bool isMultiConnected=getAscendingPathThroughWall(savedWallId, minSaddle1, minSaddle2, isVisited, &path);
        if(isMultiConnected){
          ++numberOfIterations;
          continue;
        }

#ifdef PRINT_DEBUG
        for(const Cell& cell : path)
          cout << cell.id_ << " ";
        cout << " persistence=" << vpathPersistence << endl;
#endif

        reverseAscendingPathOnWall(path);
      }

      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_){
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back(make_tuple(minSaddle1, minSaddle2));
      }

      const int sourceId=vpath.source_;
      const int destinationId=vpath.destination_;

      // invalidate vpaths connected to destination
      vector<int> newSourceIds;
      CriticalPoint& destination=criticalPoints[destinationId];
      for(const int destinationVPathId : destination.vpaths_){
        VPath& destinationVPath=vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_!=sourceId){
          // save critical point
          const int newSourceId=destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          // clear vpath
          destinationVPath.invalidate();
        }
      }

      // invalidate vpaths connected to source and save the critical points to update
      vector<int> newDestinationIds;
      CriticalPoint& source=criticalPoints[sourceId];
      for(const int sourceVPathId : source.vpaths_){
        VPath& sourceVPath=vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_!=destinationId){
          // save critical point
          const int newDestinationId=sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          CriticalPoint& newDestination=criticalPoints[newDestinationId];
          for(const int newDestinationVPathId : newDestination.vpaths_){
            VPath& newDestinationVPath=vpaths[newDestinationVPathId];
            if(newDestinationVPath.isValid_ and newDestinationVPath.source_!=sourceId){

              // clear vpath
              newDestinationVPath.invalidate();
            }
          }

          // clear vpath
          sourceVPath.invalidate();
        }
      }

      // finally invalidate current vpath and critical points
      vpath.invalidate();
      source.clear();
      destination.clear();

      // look at the gradient : reconnect locally the critical points
      for(const int newDestinationId : newDestinationIds){
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        const Cell& saddle2=newDestination.cell_;

        set<int> saddles1;
        const wallId_t savedWallId=wallId;
        getDescendingWall(wallId, saddle2, isVisited, nullptr, &saddles1);
        ++wallId;

        for(const int saddle1Id : saddles1){
          const Cell saddle1(1,saddle1Id);

          vector<Cell> path;
          const bool isMultiConnected=getAscendingPathThroughWall(savedWallId, saddle1, saddle2, isVisited, &path);
          if(isMultiConnected)
            continue;

          int newSourceId=saddle1Index[saddle1Id];

          // connection to a new saddle1 (not present in the graph before)
          if(newSourceId==-1){
            if(MustBeRemovable and !isRemovableSaddle1[saddle1Id]) continue;

            const int newCriticalPointId=criticalPoints.size();
            saddle1Index[saddle1Id]=newCriticalPointId;
            criticalPoints.push_back(CriticalPoint(saddle1));

            newSourceId=newCriticalPointId;
          }
          CriticalPoint& newSource=criticalPoints[newSourceId];

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(make_tuple(persistence,newVPathId,newDestination.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(const int newSourceId : newSourceIds){
        CriticalPoint& newSource=criticalPoints[newSourceId];
        const Cell& saddle1=newSource.cell_;

        set<int> saddles2;
        const wallId_t savedWallId=wallId;
        getAscendingWall(wallId, saddle1, isVisited, nullptr, &saddles2);
        ++wallId;

        for(const int saddle2Id : saddles2){
          const Cell saddle2(2,saddle2Id);

          vector<Cell> path;
          const bool isMultiConnected=getDescendingPathThroughWall(savedWallId, saddle2, saddle1, isVisited, &path);
          if(isMultiConnected)
            continue;

          const int newDestinationId=saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId==-1)
            continue;

          CriticalPoint& newDestination=criticalPoints[newDestinationId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists=false;
          for(const int newDestinationVPathId : newDestination.vpaths_){
            const VPath& newDestinationVPath=vpaths[newDestinationVPathId];

            if(newDestinationVPath.isValid_ and newDestinationVPath.source_==newSourceId){
              alreadyExists=true;
              break;
            }
          }

          if(alreadyExists)
            continue;

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(make_tuple(persistence,newVPathId,newDestination.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections1(vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2,
    const int iterationThreshold){
  // Part 1 : initialization
  vector<VPath> vpaths;
  vector<CriticalPoint> criticalPoints;
  vector<int> saddle1Index;
  vector<int> saddle2Index;
  initializeSaddleSaddleConnections1<dataType>(isRemovableSaddle1,
      isRemovableSaddle2,
      vpaths,
      criticalPoints,
      saddle1Index,
      saddle2Index);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  set<tuple<dataType,int,int>, SaddleSaddleVPathComparator<dataType>> S(cmp_f);
  orderSaddleSaddleConnections1<dataType>(vpaths, criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections1<dataType>(iterationThreshold,
      S,
      isRemovableSaddle1,
      isRemovableSaddle2,
      vpaths,
      criticalPoints,
      saddle1Index,
      saddle2Index);

  return 0;
}

template <typename dataType>
int DiscreteGradient::naive_reverseSaddleSaddleConnection2(vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2){
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
  wallId_t ascendingWallId=1;
  vector<wallId_t> isVisited(numberOfEdges, 0);

  bool isFirst=true;
  Cell minSaddle1;
  Cell minSaddle2;
  vector<Cell> minVPath;
  dataType minPersistence=-1;

  for(int i=0; i<numberOfEdges; ++i){
    // process removable non-boundary 2-saddle
    if(isRemovableSaddle1[i] and !inputTriangulation_->isEdgeOnBoundary(i)){
      const Cell saddle1(1,i);

      set<int> saddles2;
      const wallId_t savedAscendingWallId=ascendingWallId;
      getAscendingWall(ascendingWallId, saddle1, isVisited, nullptr, &saddles2);
      ++ascendingWallId;

      if(saddles2.empty()){
#ifdef PRINT_INFOS
        cout << "[DiscreteGradient] No 2-saddle on the wall of 1-saddle.id_=" << saddle1.id_ << endl;
#endif
        continue;
      }

      bool isPairable=false;
      const int numberOfSaddles2=saddles2.size();
      vector<vector<Cell>> vpaths(numberOfSaddles2);
      vector<char> isMultiConnected(numberOfSaddles2, false);
      int j=0;
      for(const int saddle2Id : saddles2){
        const Cell saddle2(2,saddle2Id);

        isMultiConnected[j]=getDescendingPathThroughWall(savedAscendingWallId, saddle2, saddle1, isVisited, &(vpaths[j]));

        if(!isMultiConnected[j]){
          const Cell& lastCell=vpaths[j].back();
          if(lastCell.dim_!=saddle1.dim_ or lastCell.id_!=saddle1.id_){
            ++j;
            continue;
          }
        }

        if(!isMultiConnected[j]){
          // pairing with a removable 2-saddle if possible
          if(isRemovableSaddle2[saddle2.id_]){
            isPairable=true;

            const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
            if(isFirst or persistence<minPersistence){
              // update minimizer
              isFirst=false;
              minSaddle1=saddle1;
              minSaddle2=saddle2;
              minVPath=vpaths[j];
              minPersistence=persistence;
            }
          }
        }
        ++j;
      }

      if(!isPairable and AllowReversingWithNonRemovable){
        // pairing with a non-removable 1-saddle anyway
        j=0;
        for(const int saddle2Id : saddles2){
          const Cell saddle2(2,saddle2Id);

          if(!isMultiConnected[j]){
            isPairable=true;

            const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
            if(isFirst or persistence<minPersistence){
              // update minimizer
              isFirst=false;
              minSaddle1=saddle1;
              minSaddle2=saddle2;
              minVPath=vpaths[j];
              minPersistence=persistence;
            }
          }
          ++j;
        }
      }

      if(!isPairable){
#ifdef PRINT_INFOS
        cout << "[DiscreteGradient] Non-pairable 1-saddle.id_=" << saddle1.id_ << endl;
#endif
      }
    }
  }

  if(!isFirst){
#ifdef PRINT_DEBUG
    for(const Cell& cell : minVPath)
      cout << cell.id_ << " ";
    cout << " persistence=" << minPersistence << endl;
#endif

    if(CollectPersistencePairs and outputPersistencePairs_)
      outputPersistencePairs_->push_back(make_tuple(minSaddle1, minSaddle2));

    const Cell& lastCell=minVPath.back();
    if(lastCell.dim_==minSaddle1.dim_ and lastCell.id_==minSaddle1.id_)
      reverseDescendingPathOnWall(minVPath);
#ifdef PRINT_ERROR
    else{
      cout << "[DiscreteGradient] Error : lastCell in descending path is not saddle1.id=" << minSaddle1.id_ << endl;
    }
#endif
  }

  // could not reverse any pair
  if(isFirst) return -1;

  return 0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections2(const vector<char>& isRemovableSaddle1,
    const vector<char>& isRemovableSaddle2,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    vector<int>& saddle1Index,
    vector<int>& saddle2Index) const{
    Timer t;

    const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

    const int maximumDim=dimensionality_;
    const int saddle2Dim=maximumDim-1;
    const int saddle1Dim=saddle2Dim-1;

    // Part 1 : build initial structures
    // add the 1-saddles to CriticalPointList
    const int numberOfSaddle1Candidates=getNumberOfCells(saddle1Dim);
    saddle1Index.resize(numberOfSaddle1Candidates, -1);
    for(int i=0; i<numberOfSaddle1Candidates; ++i){
      //if(MustBeRemovable and !isRemovableSaddle1[i]) continue;
      if(!isRemovableSaddle1[i]) continue;

      const Cell saddle1Candidate(saddle1Dim, i);

      if(!AllowBoundary and isBoundary(saddle1Candidate)) continue;

      if(isSaddle1(saddle1Candidate)){
        const int index=criticalPoints.size();
        saddle1Index[i]=index;
        criticalPoints.push_back(CriticalPoint(saddle1Candidate));
      }
    }
    const int numberOf1Saddles=criticalPoints.size();

    // add the 2-saddles to CriticalPointList
    const int numberOfSaddle2Candidates=getNumberOfCells(saddle2Dim);
    saddle2Index.resize(numberOfSaddle2Candidates, -1);
    for(int i=0; i<numberOfSaddle2Candidates; ++i){
      if(MustBeRemovable and !isRemovableSaddle2[i]) continue;

      const Cell saddle2Candidate(saddle2Dim, i);

      if(ForceBoundary and !isBoundary(saddle2Candidate)) continue;

      if(isSaddle2(saddle2Candidate)){
        const int index=criticalPoints.size();
        saddle2Index[i]=index;
        criticalPoints.push_back(CriticalPoint(saddle2Candidate));
      }
    }

    // Part 2 : update the structures
    // apriori: by default construction, the vpaths and segments are not valid
    wallId_t ascendingWallId=1;
    vector<wallId_t> isVisited(numberOfSaddle1Candidates, 0);
    for(int i=0; i<numberOf1Saddles; ++i){
      const int sourceIndex=i;
      CriticalPoint& source=criticalPoints[sourceIndex];
      const Cell& saddle1=source.cell_;

      set<int> saddles2;
      const wallId_t savedAscendingWallId=ascendingWallId;
      getAscendingWall(ascendingWallId, saddle1, isVisited, nullptr, &saddles2);
      ++ascendingWallId;

      for(const int saddle2Id : saddles2){
        if(MustBeRemovable and !isRemovableSaddle2[saddle2Id]) continue;

        const Cell& saddle2=Cell(2,saddle2Id);

        if(ForceBoundary and !isBoundary(saddle2)) continue;

        vector<Cell> path;
        const bool isMultiConnected=getDescendingPathThroughWall(savedAscendingWallId, saddle2, saddle1, isVisited, &path);

        if(!isMultiConnected){
          const int destinationIndex=saddle2Index[saddle2Id];
          CriticalPoint& destination=criticalPoints[destinationIndex];

          // update source and destination
          const int sourceSlot=source.addSlot();
          const int destinationSlot=destination.addSlot();

          // update vpath
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
          vpaths.push_back(VPath(true,-1,sourceIndex,destinationIndex,sourceSlot,destinationSlot,persistence));
        }
      }
    }

    // Part 3 : initialize the last structures
    const int numberOfCriticalPoints=criticalPoints.size();
    for(int i=0; i<numberOfCriticalPoints; ++i){
      CriticalPoint& cp=criticalPoints[i];

      const int numberOfSlots=cp.numberOfSlots_;
      cp.vpaths_.resize(numberOfSlots);
      cp.numberOfSlots_=0;
    }

    const int numberOfVPaths=vpaths.size();
#ifdef withOpenMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i=0; i<numberOfVPaths; ++i){
      const VPath& vpath=vpaths[i];

      if(vpath.isValid_){
        const int sourceIndex=vpath.source_;
        const int destinationIndex=vpath.destination_;

        const int sourceSlot=vpath.sourceSlot_;
        const int destinationSlot=vpath.destinationSlot_;

        CriticalPoint& source=criticalPoints[sourceIndex];
        CriticalPoint& destination=criticalPoints[destinationIndex];

        source.vpaths_[sourceSlot]=i;
        destination.vpaths_[destinationSlot]=i;
      }
    }

    {
      stringstream msg;
      msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  return 0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections2(const vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const int saddleId=criticalPoints[vpath.source_].cell_.id_;
      S.insert(make_tuple(vpath.persistence_,i,saddleId));
    }
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleSaddleConnections2(const int iterationThreshold,
    set<tuple<dataType,int,int>,SaddleSaddleVPathComparator<dataType>>& S,
    vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2,
    vector<VPath>& vpaths,
    vector<CriticalPoint>& criticalPoints,
    vector<int>& saddle1Index,
    vector<int>& saddle2Index){
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
  const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  const int optimizedSize=std::max(numberOfEdges, numberOfTriangles);
  wallId_t wallId=1;
  vector<wallId_t> isVisited(optimizedSize, 0);

  int numberOfIterations{};
  while(!S.empty()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=get<1>(*ptr);
#ifdef PRINT_DEBUG
    const dataType vpathPersistence=get<0>(*ptr);
#endif
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    if(vpath.isValid_){
      {
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;

        set<int> saddles2;
        const wallId_t savedWallId=wallId;
        getAscendingWall(wallId, minSaddle1, isVisited, nullptr, &saddles2);
        ++wallId;

        // check if at least one connection exists
        auto isFound=saddles2.find(minSaddle2.id_);
        if(isFound==saddles2.end()){
          ++numberOfIterations;
          continue;
        }

        // check if there is multiple connections
        vector<Cell> path;
        const bool isMultiConnected=getDescendingPathThroughWall(savedWallId, minSaddle2, minSaddle1, isVisited, &path);
        if(isMultiConnected){
          ++numberOfIterations;
          continue;
        }

#ifdef PRINT_DEBUG
        for(const Cell& cell : path)
          cout << cell.id_ << " ";
        cout << " persistence=" << vpathPersistence << endl;
#endif

        reverseDescendingPathOnWall(path);
      }

      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_){
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back(make_tuple(minSaddle1, minSaddle2));
      }

      const int sourceId=vpath.source_;
      const int destinationId=vpath.destination_;

      // invalidate vpaths connected to source
      vector<int> newDestinationIds;
      CriticalPoint& source=criticalPoints[sourceId];
      for(const int sourceVPathId : source.vpaths_){
        VPath& sourceVPath=vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_!=destinationId){
          // save critical point
          const int newDestinationId=sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          // clear vpath
          sourceVPath.invalidate();
        }
      }

      // invalidate vpaths connected to destination and save the critical points to update
      vector<int> newSourceIds;
      CriticalPoint& destination=criticalPoints[destinationId];
      for(const int destinationVPathId : destination.vpaths_){
        VPath& destinationVPath=vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_!=sourceId){
          // save critical point
          const int newSourceId=destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          CriticalPoint& newSource=criticalPoints[newSourceId];
          for(const int newSourceVPathId : newSource.vpaths_){
            VPath& newSourceVPath=vpaths[newSourceVPathId];
            if(newSourceVPath.isValid_ and newSourceVPath.destination_!=destinationId){

              // clear vpath
              newSourceVPath.invalidate();
            }
          }

          // clear vpath
          destinationVPath.invalidate();
        }
      }

      // finally invalidate current vpath and critical points
      vpath.invalidate();
      source.clear();
      destination.clear();

      // look at the gradient : reconnect locally the critical points
      for(const int newSourceId : newSourceIds){
        CriticalPoint& newSource=criticalPoints[newSourceId];
        const Cell& saddle1=newSource.cell_;

        set<int> saddles2;
        const wallId_t savedWallId=wallId;
        getAscendingWall(wallId, saddle1, isVisited, nullptr, &saddles2);
        ++wallId;

        for(const int saddle2Id : saddles2){
          const Cell saddle2(2,saddle2Id);

          const bool isMultiConnected=getDescendingPathThroughWall(savedWallId, saddle2, saddle1, isVisited, nullptr);
          if(isMultiConnected)
            continue;

          int newDestinationId=saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId==-1){
            if(MustBeRemovable and !isRemovableSaddle2[saddle2Id]) continue;

            const int newCriticalPointId=criticalPoints.size();
            saddle2Index[saddle2Id]=newCriticalPointId;
            criticalPoints.push_back(CriticalPoint(saddle2));

            newDestinationId=newCriticalPointId;
          }

          CriticalPoint& newDestination=criticalPoints[newDestinationId];

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(make_tuple(persistence,newVPathId,newSource.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(const int newDestinationId : newDestinationIds){
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        const Cell& saddle2=newDestination.cell_;

        set<int> saddles1;
        const wallId_t savedWallId=wallId;
        getDescendingWall(wallId, saddle2, isVisited, nullptr, &saddles1);
        ++wallId;

        for(const int saddle1Id : saddles1){
          const Cell saddle1(1,saddle1Id);

          vector<Cell> path;
          const bool isMultiConnected=getAscendingPathThroughWall(savedWallId, saddle1, saddle2, isVisited, &path);
          if(isMultiConnected)
            continue;

          const int newSourceId=saddle1Index[saddle1Id];

          if(newSourceId==-1)
            continue;

          CriticalPoint& newSource=criticalPoints[newSourceId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists=false;
          for(const int newSourceVPathId : newSource.vpaths_){
            const VPath& newSourceVPath=vpaths[newSourceVPathId];

            if(newSourceVPath.isValid_ and newSourceVPath.destination_==newDestinationId){
              alreadyExists=true;
              break;
            }
          }

          if(alreadyExists)
            continue;

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1, scalars);
          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(make_tuple(persistence,newVPathId,newSource.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  {
    stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections2(vector<char>& isRemovableSaddle1,
    vector<char>& isRemovableSaddle2,
    const int iterationThreshold){

  // Part 1 : initialization
  vector<VPath> vpaths;
  vector<CriticalPoint> criticalPoints;
  vector<int> saddle1Index;
  vector<int> saddle2Index;
  initializeSaddleSaddleConnections2<dataType>(isRemovableSaddle1,
      isRemovableSaddle2,
      vpaths,
      criticalPoints,
      saddle1Index,
      saddle2Index);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  set<tuple<dataType,int,int>, SaddleSaddleVPathComparator<dataType>> S(cmp_f);
  orderSaddleSaddleConnections2<dataType>(vpaths, criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections2<dataType>(iterationThreshold,
      S,
      isRemovableSaddle1,
      isRemovableSaddle2,
      vpaths,
      criticalPoints,
      saddle1Index,
      saddle2Index);

  return 0;
}

template<typename dataType>
int DiscreteGradient::reverseGradient(const vector<pair<int,char>>& criticalPoints){
  if(ReverseSaddleMaximumConnection){
    vector<char> isRemovableMaximum;
    getRemovableMaxima<dataType>(criticalPoints, isRemovableMaximum);

    vector<char> isRemovableSaddle;
    if(dimensionality_==2){
      const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
      isRemovableSaddle.resize(numberOfEdges, true);
    }
    else if(dimensionality_==3)
      getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle);

    const int needFirstPass=std::count(isRemovableMaximum.begin(), isRemovableMaximum.end(), true);
    if(!needFirstPass){
      stringstream msg;
      msg << "[DiscreteGradient] No DMT-maxima to remove, first pass of gradient reversal is disabled." << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }

    if(needFirstPass){
      simplifySaddleMaximumConnections<dataType>(isRemovableMaximum, isRemovableSaddle, IterationThreshold, true);

      int n=0;
      bool needSecondPass=false;
      const int maximumDim=dimensionality_;
      const int numberOfMaximumCandidates=getNumberOfCells(maximumDim);
      for(int i=0; i<numberOfMaximumCandidates; ++i){
        const Cell maximum(maximumDim, i);
        if(isCellCritical(maximum) and isRemovableMaximum[i]){
          needSecondPass=true;
          ++n;
        }
      }

      if(needSecondPass){
        cout << "MSC max not removed: " << n<< endl;
      }

      /*
      if(needSecondPass)
        simplifySaddleMaximumConnections<dataType>(isRemovableMaximum, IterationThreshold, true);
      else{
        stringstream msg;
        msg << "[DiscreteGradient] No DMT-maxima to remove, second pass of gradient reversal is disabled." << endl;
        dMsg(cout, msg.str(), advancedInfoMsg);
      }
      */
    }
    else{
      stringstream msg;
      msg << "[DiscreteGradient] No DMT-maxima to remove, first pass of gradient reversal is disabled." << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
  }

#ifdef ALLOW_PL_COMPLIANT_DETECTOR
  if(dimensionality_==3){
    for(pair<int,char> criticalPoint : criticalPoints){
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(inputTriangulation_->isVertexOnBoundary(criticalPointId))
        continue;

      switch(criticalPointType){
        case 0:
          if(!isMinimum(Cell(0,criticalPointId))){
            cout << "[DiscreteGradient] Non PL-compliant PL-minimum vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
            exit(-1);
#endif
          }
          break;

        case 1:
          {
            bool isFound=false;
            const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(criticalPointId);
            for(int i=0; i<edgeNumber; ++i){
              int edgeId;
              inputTriangulation_->getVertexEdge(criticalPointId, i, edgeId);

              if(isSaddle1(Cell(1,edgeId))){
                isFound=true;
                break;
              }
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-saddle1 vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;

        case 2:
          {
            bool isFound=false;
            const int tetraNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);
            for(int i=0; i<tetraNumber; ++i){
              int tetraId;
              inputTriangulation_->getVertexStar(criticalPointId, i, tetraId);

              for(int j=0; j<4; ++j){
                int triangleId;
                inputTriangulation_->getCellTriangle(tetraId, j, triangleId);

                bool isValid=false;
                for(int k=0; k<3; ++k){
                  int vertexId;
                  inputTriangulation_->getTriangleVertex(triangleId, k, vertexId);

                  if(vertexId==criticalPointId){
                    isValid=true;
                    break;
                  }
                }

                if(isValid and isSaddle2(Cell(2,triangleId))){
                  isFound=true;
                  break;
                }
              }
              if(isFound)
                break;
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-saddle2 vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;

        case 3:
          {
            bool isFound=false;
            const int tetraNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);
            for(int i=0; i<tetraNumber; ++i){
              int tetraId;
              inputTriangulation_->getVertexStar(criticalPointId, i, tetraId);

              if(isMaximum(Cell(3,tetraId))){
                isFound=true;
                break;
              }
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-maximum vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;
      }
    }
    cout << "[DiscreteGradient] DMT cells are PL-compliant." << endl;
  }
#endif

#ifdef ALLOW_CYCLE_DETECTOR
  if(dimensionality_==3){
    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    for(int i=0; i<numberOfTriangles; ++i){
      if(isSaddle2(Cell(2,i))){
        const int starNumber=inputTriangulation_->getTriangleStarNumber(i);
        for(int j=0; j<starNumber; ++j){
          int starId;
          inputTriangulation_->getTriangleStar(i, j, starId);

          vector<Cell> path;
          getAscendingPath(Cell(3,starId), path, true);
        }
      }
    }
  }
  cout << "[DiscreteGradient] No cycle detected." << endl;
#endif

  // experimental
  if(dimensionality_==3 and ReverseSaddleSaddleConnection){
    vector<char> isRemovableSaddle1;
    vector<char> isRemovableSaddle2;

#ifdef PRINT_DEBUG
    {
      const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
      const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();

      getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
      getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);

      int nbsaddle1{};
      for(int i=0; i<numberOfEdges; ++i){
        if(isRemovableSaddle1[i] and !inputTriangulation_->isEdgeOnBoundary(i))
          ++nbsaddle1;
      }

      int nbsaddle2{};
      for(int i=0; i<numberOfTriangles; ++i){
        if(isRemovableSaddle2[i] and !inputTriangulation_->isTriangleOnBoundary(i))
          ++nbsaddle2;
      }

      int bsaddle1{};
      for(int i=0; i<numberOfEdges; ++i){
        if(isRemovableSaddle1[i] and inputTriangulation_->isEdgeOnBoundary(i))
          ++bsaddle1;
      }

      int bsaddle2{};
      for(int i=0; i<numberOfTriangles; ++i){
        if(isRemovableSaddle2[i] and inputTriangulation_->isTriangleOnBoundary(i))
          ++bsaddle2;
      }
      cout << "\nnumber of non-boundary 1-saddles to remove : "<< nbsaddle1 << endl;
      cout << "number of non-boundary 2-saddles to remove : "<< nbsaddle2 << endl;
      cout << "\nnumber of boundary 1-saddles to remove : "<< bsaddle1 << endl;
      cout << "number of boundary 2-saddles to remove : "<< bsaddle2 << endl;
    }
#endif

    Timer t;
#ifdef NAIVE_VERSION
    {
      int ret=0;
      do{
        getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
        getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);
        ret=naive_reverseSaddleSaddleConnection<dataType>(isRemovableSaddle1, isRemovableSaddle2);
      }while(!ret);

      ret=0;
      do{
        getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
        getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);
        ret=naive_reverseSaddleSaddleConnection2<dataType>(isRemovableSaddle1, isRemovableSaddle2);
      }while(!ret);
    }
#endif

#ifdef PROTO_VERSION
    getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
    getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);

    // walls computed from 2-saddles
    simplifySaddleSaddleConnections1<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);

    getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
    getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);

    // walls computed from 1-saddles
    simplifySaddleSaddleConnections2<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);

    // experimental
//#define REVERSE_LITTLE_PAIRS_INSIDE
#ifdef REVERSE_LITTLE_PAIRS_INSIDE
    {
      getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
      getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);

      // 2-saddle
      const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
      vector<char> isPLSaddle2(numberOfVertices, false);
      for(const pair<int, char>& criticalPoint : criticalPoints){
        const int criticalPointId=criticalPoint.first;
        const char criticalPointType=criticalPoint.second;

        if(criticalPointType==2)
          isPLSaddle2[criticalPointId]=true;
      }

      // add 2-saddles of the stars
      vector<vector<int>> saddles2(numberOfVertices);
      const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
      for(int i=0; i<numberOfTriangles; ++i){
        const Cell saddle2Candidate(2,i);

        if(!isSaddle2(saddle2Candidate) or isBoundary(saddle2Candidate)) continue;

        for(int j=0; j<3; ++j){
          int vertexId;
          inputTriangulation_->getTriangleVertex(i, j, vertexId);

          if(isPLSaddle2[vertexId])
            saddles2[vertexId].push_back(i);
        }
      }

      // add 2-saddles of the link
      /*
      for(int i=0; i<numberOfTriangles; ++i){
        const Cell saddle2Candidate(2,i);

        if(!isSaddle2(saddle2Candidate) or isBoundary(saddle2Candidate)) continue;

        vector<int> pl_saddles2;
        const int starNumber=inputTriangulation_->getTriangleStarNumber(i);
        for(int j=0; j<starNumber; ++j){
          int starId;
          inputTriangulation_->getTriangleStar(i, j, starId);

          for(int k=0; k<4; ++k){
            int vertexId;
            inputTriangulation_->getCellVertex(starId, k, vertexId);

            if(isPLSaddle2[vertexId])
              pl_saddles2.push_back(vertexId);
          }
        }

        if(pl_saddles2.size()==1){
          const int vertexId=pl_saddles2[0];
          saddles2[vertexId].push_back(i);
        }
      }
      */

      // invert
      for(int i=0; i<numberOfVertices; ++i){
        if(isPLSaddle2[i] and saddles2[i].size()>1){
          for(const int saddle2Id : saddles2[i]){
            isRemovableSaddle2[saddle2Id]=!isRemovableSaddle2[saddle2Id];
          }
        }
      }

      simplifySaddleSaddleConnections1<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);
    }
#endif

//#define LAST_PASS
#ifdef LAST_PASS
    {
      MustBeRemovable=false;
      ForceBoundary=true;
      simplifySaddleSaddleConnections1<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);
      simplifySaddleSaddleConnections2<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);
      ForceBoundary=false;
      MustBeRemovable=true;
    }
#endif

//#define REVERSE_LITTLE_PAIRS_BOUNDARY
#ifdef REVERSE_LITTLE_PAIRS_BOUNDARY
    {
      AllowBoundary=true;
      // 1-saddle
      // add 1-saddles of the stars
      const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
      vector<char> isPLSaddle1(numberOfVertices, false);
      vector<vector<int>> saddles1(numberOfVertices);
      for(const pair<int, char>& criticalPoint : criticalPoints){
        const int criticalPointId=criticalPoint.first;
        const char criticalPointType=criticalPoint.second;

        if(!inputTriangulation_->isVertexOnBoundary(criticalPointId))
          continue;

        if(criticalPointType==1){
          isPLSaddle1[criticalPointId]=true;

          const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(criticalPointId);
          for(int j=0; j<edgeNumber; ++j){
            int edgeId;
            inputTriangulation_->getVertexEdge(criticalPointId, j, edgeId);

            if(isSaddle1(Cell(1,edgeId)) and inputTriangulation_->isEdgeOnBoundary(edgeId))
              saddles1[criticalPointId].push_back(edgeId);
          }
        }
      }

      // invert
      for(int i=0; i<numberOfVertices; ++i){
        if(isPLSaddle1[i] and saddles1[i].size()>1){
          for(const int saddle1Id : saddles1[i]){
            isRemovableSaddle1[saddle1Id]=!isRemovableSaddle1[saddle1Id];
          }
        }
      }

      simplifySaddleSaddleConnections2<dataType>(isRemovableSaddle1, isRemovableSaddle2, IterationThreshold);
      AllowBoundary=false;
    }
#endif

#ifdef PRINT_PB1
    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    vector<wallId_t> isVisited(numberOfTriangles, -1);
    set<int> saddles1;
    const Cell start(2,13124973);
    getDescendingWall(0, start, isVisited, nullptr, &saddles1);

    cout << "saddles=" << endl;
    for(int k : saddles1){
      cout << k << " isOnBoundary=" << inputTriangulation_->isEdgeOnBoundary(k);
      bool isMultiConnected=getAscendingPathThroughWall(0, Cell(1,k), start, isVisited, nullptr);
      cout << " isMultiConnected=" << isMultiConnected;
      cout << " isRemovable=" << (int)isRemovableSaddle1[k] << endl;
    }
    cout << endl;
#endif

#endif
    {
      stringstream msg;
      msg << "[DiscreteGradient] Saddle-Saddle pairs ";
      if(CollectPersistencePairs)
        msg << "collected in ";
      else
        msg << "simplified in ";
      msg << t.getElapsedTime() << " s, "<< threadNumber_ << " thread(s)." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }

#ifdef ALLOW_CYCLE_DETECTOR
    // cycle detector 1)
    {
      const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();

      vector<Cell> cp;
      getCriticalPoints(cp);

      vector<wallId_t> isVisited(numberOfTriangles, 0);

      const int numberOfCriticalPoints=cp.size();
      for(int i=0; i<numberOfCriticalPoints; ++i){
        const Cell& criticalPoint=cp[i];

        if(criticalPoint.dim_==2){
          const Cell saddle2=criticalPoint;

          set<int> saddles1;
          getDescendingWall(saddle2.id_, saddle2, isVisited, nullptr, &saddles1);

          for(const int saddle1Id : saddles1){
            const Cell& saddle1=Cell(1,saddle1Id);

            getAscendingPathThroughWall(saddle2.id_, saddle1, saddle2, isVisited, nullptr, true);
          }
        }
      }
    }

    // cycle detector 2)
    {
      const int numberOfEdges=inputTriangulation_->getNumberOfEdges();

      vector<Cell> cp;
      getCriticalPoints(cp);

      vector<wallId_t> isVisited(numberOfEdges, 0);

      const int numberOfCriticalPoints=cp.size();
      for(int i=0; i<numberOfCriticalPoints; ++i){
        const Cell& criticalPoint=cp[i];

        if(criticalPoint.dim_==1){
          const Cell saddle1=criticalPoint;

          set<int> saddles2;
          getAscendingWall(saddle1.id_, saddle1, isVisited, nullptr, &saddles2);

          for(const int saddle2Id : saddles2){
            const Cell& saddle2=Cell(2,saddle2Id);

            getDescendingPathThroughWall(saddle1.id_, saddle2, saddle1, isVisited, nullptr, true);
          }
        }
      }
    }

    cout << "[DiscreteGradient] No cycle detected." << endl;
#endif

#ifdef ALLOW_PL_COMPLIANT_DETECTOR
  if(dimensionality_==3){
    for(pair<int,char> criticalPoint : criticalPoints){
      const int criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(inputTriangulation_->isVertexOnBoundary(criticalPointId))
        continue;

      switch(criticalPointType){
        case 0:
          if(!isMinimum(Cell(0,criticalPointId))){
            cout << "[DiscreteGradient] Non PL-compliant PL-minimum vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
            exit(-1);
#endif
          }
          break;

        case 1:
          {
            bool isFound=false;
            const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(criticalPointId);
            for(int i=0; i<edgeNumber; ++i){
              int edgeId;
              inputTriangulation_->getVertexEdge(criticalPointId, i, edgeId);

              if(isSaddle1(Cell(1,edgeId))){
                isFound=true;
                break;
              }
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-saddle1 vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;

        case 2:
          {
            bool isFound=false;
            const int tetraNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);
            for(int i=0; i<tetraNumber; ++i){
              int tetraId;
              inputTriangulation_->getVertexStar(criticalPointId, i, tetraId);

              for(int j=0; j<4; ++j){
                int triangleId;
                inputTriangulation_->getCellTriangle(tetraId, j, triangleId);

                bool isValid=false;
                for(int k=0; k<3; ++k){
                  int vertexId;
                  inputTriangulation_->getTriangleVertex(triangleId, k, vertexId);

                  if(vertexId==criticalPointId){
                    isValid=true;
                    break;
                  }
                }

                if(isValid and isSaddle2(Cell(2,triangleId))){
                  isFound=true;
                  break;
                }
              }
              if(isFound)
                break;
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-saddle2 vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;

        case 3:
          {
            bool isFound=false;
            const int tetraNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);
            for(int i=0; i<tetraNumber; ++i){
              int tetraId;
              inputTriangulation_->getVertexStar(criticalPointId, i, tetraId);

              if(isMaximum(Cell(3,tetraId))){
                isFound=true;
                break;
              }
            }
            if(!isFound){
              cout << "[DiscreteGradient] Non PL-compliant PL-maximum vertexId=" << criticalPointId << endl;
#ifdef ALLOW_EXIT
              exit(-1);
#endif
            }
          }
          break;
      }
    }
    cout << "[DiscreteGradient] DMT cells are PL-compliant." << endl;
  }
#endif

#ifdef PRINT_DEBUG
    {
      const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
      const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();

      getRemovableSaddles1<dataType>(criticalPoints, isRemovableSaddle1);
      getRemovableSaddles2<dataType>(criticalPoints, isRemovableSaddle2);

      int nbsaddle1{};
      for(int i=0; i<numberOfEdges; ++i){
        if(isRemovableSaddle1[i] and !inputTriangulation_->isEdgeOnBoundary(i))
          ++nbsaddle1;
      }

      int nbsaddle2{};
      for(int i=0; i<numberOfTriangles; ++i){
        if(isRemovableSaddle2[i] and !inputTriangulation_->isTriangleOnBoundary(i))
          ++nbsaddle2;
      }

      int bsaddle1{};
      for(int i=0; i<numberOfEdges; ++i){
        if(isRemovableSaddle1[i] and inputTriangulation_->isEdgeOnBoundary(i))
          ++bsaddle1;
      }

      int bsaddle2{};
      for(int i=0; i<numberOfTriangles; ++i){
        if(isRemovableSaddle2[i] and inputTriangulation_->isTriangleOnBoundary(i))
          ++bsaddle2;
      }
      cout << "\nnumber of non-boundary 1-saddles to remove : "<< nbsaddle1 << endl;
      cout << "number of non-boundary 2-saddles to remove : "<< nbsaddle2 << endl;
      cout << "\nnumber of boundary 1-saddles to remove : "<< bsaddle1 << endl;
      cout << "number of boundary 2-saddles to remove : "<< bsaddle2 << endl;
    }
#endif
  }

#ifdef ALLOW_DOUBLE_PAIRING_DETECTOR
  if(dimensionality_==3){
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
    const int numberOfEdges=inputTriangulation_->getNumberOfEdges();
    const int numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
    const int numberOfCells=inputTriangulation_->getNumberOfCells();

    vector<int> edgeDetector(numberOfEdges, -1);
    for(int i=0; i<numberOfVertices; ++i){
      const int pairedCellId=gradient_[0][0][i];
      if(pairedCellId!=-1){
        if(edgeDetector[pairedCellId]==-1)
          edgeDetector[pairedCellId]=i;
        else{
          cout << "Double-paired edgeId=" << pairedCellId << endl;
#ifdef ALLOW_EXIT
          exit(-1);
#endif
        }
      }
    }
    for(int i=0; i<numberOfTriangles; ++i){
      const int pairedCellId=gradient_[1][2][i];
      if(pairedCellId!=-1){
        if(edgeDetector[pairedCellId]==-1)
          edgeDetector[pairedCellId]=i;
        else{
          cout << "Double-paired edgeId=" << pairedCellId << endl;
#ifdef ALLOW_EXIT
          exit(-1);
#endif
        }
      }
    }

    vector<int> triangleDetector(numberOfTriangles, -1);
    for(int i=0; i<numberOfEdges; ++i){
      const int pairedCellId=gradient_[1][1][i];
      if(pairedCellId!=-1){
        if(triangleDetector[pairedCellId]==-1)
          triangleDetector[pairedCellId]=i;
        else{
          cout << "Double-paired triangleId=" << pairedCellId << endl;
#ifdef ALLOW_EXIT
          exit(-1);
#endif
        }
      }
    }
    for(int i=0; i<numberOfCells; ++i){
      const int pairedCellId=gradient_[2][3][i];
      if(pairedCellId!=-1){
        if(triangleDetector[pairedCellId]==-1)
          triangleDetector[pairedCellId]=i;
        else{
          cout << "Double-paired triangleId=" << pairedCellId << endl;
#ifdef ALLOW_EXIT
          exit(-1);
#endif
        }
      }
    }

    cout << "[DiscreteGradient] No double-paired cell detected." << endl;
  }
#endif

  return 0;
}

template<typename dataType>
int DiscreteGradient::reverseGradient(){
  {
    // foreach dimension
    const int numberOfDimensions=getNumberOfDimensions();
    vector<int> numberOfCriticalPointsByDimension(numberOfDimensions,0);
    for(int i=0; i<numberOfDimensions; ++i){

      // foreach cell of that dimension
      const int numberOfCells=getNumberOfCells(i);
      for(int j=0; j<numberOfCells; ++j){
        const Cell cell(i,j);

        if(isCellCritical(cell))
          ++numberOfCriticalPointsByDimension[i];
      }
    }

    {
      stringstream msg;
      for(int i=0; i<numberOfDimensions; ++i)
        msg << "[DiscreteGradient] " << numberOfCriticalPointsByDimension[i] << " " << i << "-cell(s)." << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  }

  //vector<pair<int,char>> criticalPoints;
  criticalPoints_.clear();

  // get the PL critical points
  if(ReverseSaddleMaximumConnection or ReverseSaddleSaddleConnection){
    const int numberOfVertices=inputTriangulation_->getNumberOfVertices();

    const int* const offsets=static_cast<int*>(inputOffsets_);
    vector<int> sosOffsets(numberOfVertices);
    for(int i=0; i<numberOfVertices; ++i)
      sosOffsets[i]=offsets[i];

    ScalarFieldCriticalPoints<dataType> scp;

    scp.setDebugLevel(debugLevel_);
    scp.setThreadNumber(threadNumber_);
    scp.setDomainDimension(dimensionality_);
    scp.setScalarValues(inputScalarField_);
    scp.setVertexNumber(numberOfVertices);
    scp.setSosOffsets(&sosOffsets);
    scp.setupTriangulation(inputTriangulation_);
    scp.setOutput(&criticalPoints_);

    scp.execute();
  }

  reverseGradient<dataType>(criticalPoints_);

  return 0;
}

#endif // DISCRETEGRADIENT_H
