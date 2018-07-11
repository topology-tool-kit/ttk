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
      template<typename dataType, typename idtype>
        int execute();

      /**
       * Compute the (saddle1, saddle2) pairs not detected by the
       * contour tree.
       */
      template<typename dataType, typename idType>
        int computePersistencePairs(const std::vector<std::tuple<dcg::simplexId_t, dcg::simplexId_t, dataType>>& JTPairs,
            const std::vector<std::tuple<dcg::simplexId_t, dcg::simplexId_t, dataType>>& STPairs,
            std::vector<std::tuple<dcg::simplexId_t, dcg::simplexId_t, dataType>>& pl_saddleSaddlePairs);

      template <typename dataType>
        int setAugmentedCriticalPoints(const std::vector<dcg::Cell>& criticalPoints,
            dcg::simplexId_t* ascendingManifold,
            dcg::simplexId_t* descendingManifold) const;

      /**
       * Compute the descending 1-separatrices by reading into the discrete
       * gradient.
       */
      int getAscendingSeparatrices1(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const;

      /**
       * Compute the saddle-connectors by reading into the discrete
       * gradient.
       */
      int getSaddleConnectors(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const;

      /**
       * Compute the geometrical embedding of the saddle-connectors. This
       * function needs the following internal pointers to be set:
       * outputSeparatrices1_numberOfPoints_
       * outputSeparatrices1_points_
       * outputSeparatrices1_numberOfCells_
       * outputSeparatrices1_cells_
       * inputScalarField_
       */
      template<typename dataType>
      int setSaddleConnectors(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const;

      /**
       * Compute the 2-separatrices by reading into the discrete
       * gradient from the maxima.
       */
      int getDescendingSeparatrices2(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
          std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;

      /**
       * Compute the geometrical embedding of the descending
       * 2-separatrices. This function needs the following
       * internal pointers to be set:
       * outputSeparatrices2_numberOfPoints_
       * outputSeparatrices2_points_
       * outputSeparatrices2_numberOfCells_
       * outputSeparatrices2_cells_
       * inputScalarField_
       */
      template<typename dataType>
      int setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
          const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;
#ifdef TTK_ENABLE_OPENMP      
      template<typename dataType>
        int omp_setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
            const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
            const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;
#endif

      int getDualPolygon(const dcg::simplexId_t edgeId, std::vector<dcg::simplexId_t>& polygon) const;

      int sortDualPolygonVertices(std::vector<dcg::simplexId_t>& polygon) const;

      /**
       * Compute the 2-separatrices by reading into the discrete
       * gradient from the minima.
       */
      int getAscendingSeparatrices2(const std::vector<dcg::Cell>& criticalPoints,
          std::vector<Separatrix>& separatrices,
          std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
          std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;

      /**
       * Compute the geometrical embedding of the ascending
       * 2-separatrices. This function needs the following
       * internal pointers to be set:
       * outputSeparatrices2_numberOfPoints_
       * outputSeparatrices2_points_
       * outputSeparatrices2_numberOfCells_
       * outputSeparatrices2_cells_
       * inputScalarField_
       */
      template<typename dataType>
      int setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
          const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
          const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;

#ifdef TTK_ENABLE_OPENMP      
      template<typename dataType>
        int omp_setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
            const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
            const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const;
#endif      
  };
}

template<typename dataType>
int ttk::MorseSmaleComplex3D::setSaddleConnectors(const 
std::vector<Separatrix>& separatrices,
    const std::vector<std::vector<dcg::Cell>>& separatricesGeometry) const{
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices1_numberOfPoints_){
    std::cerr << "[MorseSmaleComplex3D] 1-separatrices pointer to numberOfPoints is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_points_){
    std::cerr << "[MorseSmaleComplex3D] 1-separatrices pointer to points is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_numberOfCells_){
    std::cerr << "[MorseSmaleComplex3D] 1-separatrices pointer to numberOfCells is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices1_cells_){
    std::cerr << "[MorseSmaleComplex3D] 1-separatrices pointer to cells is null." << std::endl;
    return -1;
  }
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] 1-separatrices pointer to the input scalar field is null." << std::endl;
    return -1;
  }
#endif
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
  if(outputSeparatrices1_cells_separatrixIds_ and
      outputSeparatrices1_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices1_cells_separatrixIds_->begin(),
        outputSeparatrices1_cells_separatrixIds_->end())+1;
  }

  for(const Separatrix& separatrix : separatrices){
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle1=separatrix.source_;
    const dcg::Cell& saddle2=separatrix.destination_;

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
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      dcg::simplexId_t oldPointId=-1;
      for(auto cellIte=separatricesGeometry[geometryId].begin(); cellIte!=separatricesGeometry[geometryId].end(); ++cellIte){
        const dcg::Cell& cell=*cellIte;
        float point[3];
        discreteGradient_.getCellIncenter(cell, point);

        outputSeparatrices1_points_->push_back(point[0]);
        outputSeparatrices1_points_->push_back(point[1]);
        outputSeparatrices1_points_->push_back(point[2]);

        if(outputSeparatrices1_points_smoothingMask_){
          if(cellIte==separatricesGeometry[geometryId].begin() or
              cellIte==separatricesGeometry[geometryId].end()-1)
            outputSeparatrices1_points_smoothingMask_->push_back(0);
          else
            outputSeparatrices1_points_smoothingMask_->push_back(1);
        }
        if(outputSeparatrices1_points_cellDimensions_)
          outputSeparatrices1_points_cellDimensions_->push_back(cell.dim_);
        if(outputSeparatrices1_points_cellIds_)
          outputSeparatrices1_points_cellIds_->push_back(cell.id_);

        if(oldPointId!=-1){
          outputSeparatrices1_cells_->push_back(2);
          outputSeparatrices1_cells_->push_back(oldPointId);
          outputSeparatrices1_cells_->push_back(pointId);

          if(outputSeparatrices1_cells_sourceIds_)
            outputSeparatrices1_cells_sourceIds_->push_back(saddle1.id_);
          if(outputSeparatrices1_cells_destinationIds_)
            outputSeparatrices1_cells_destinationIds_->push_back(saddle2.id_);
          if(outputSeparatrices1_cells_separatrixIds_)
            outputSeparatrices1_cells_separatrixIds_->push_back(separatrixId);
          if(outputSeparatrices1_cells_separatrixTypes_)
            outputSeparatrices1_cells_separatrixTypes_->push_back(separatrixType);
          if(outputSeparatrices1_cells_separatrixFunctionMaxima)
            outputSeparatrices1_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
          if(outputSeparatrices1_cells_separatrixFunctionMinima)
            outputSeparatrices1_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
          if(outputSeparatrices1_cells_separatrixFunctionDiffs)
            outputSeparatrices1_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
          if(outputSeparatrices1_cells_isOnBoundary_)
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

#ifdef TTK_ENABLE_OPENMP      
template<typename dataType>
int ttk::MorseSmaleComplex3D::omp_setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
   const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
   const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const{
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices2_numberOfPoints_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfPoints is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_points_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_numberOfCells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfCells is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_cells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null." << std::endl;
    return -1;
  }
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input scalar field is null." << std::endl;
    return -1;
  }
#endif

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  dcg::simplexId_t pointId=(*outputSeparatrices2_numberOfPoints_);
  dcg::simplexId_t separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_ and
      outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  std::vector<dcg::simplexId_t> separatrixIds(threadNumber_, 0);
  std::vector<dcg::simplexId_t> numberOfPoints(threadNumber_, 0);
  std::vector<std::vector<float>> separatrices2_points(threadNumber_);
  std::vector<dcg::simplexId_t> numberOfCells(threadNumber_, 0);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells(threadNumber_);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells_sourceIds(threadNumber_);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells_separatrixIds(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_separatrixTypes(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMaxima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMinima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionDiffs(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_isOnBoundary(threadNumber_);

  const dcg::simplexId_t numberOfSeparatrices=separatrices.size();
#pragma omp parallel for num_threads(threadNumber_)
  for(dcg::simplexId_t i=0; i<numberOfSeparatrices; ++i){
    const ThreadId threadId=omp_get_thread_num();
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle=separatrix.source_;
    const char separatrixType=1;
    const dcg::simplexId_t saddleId=saddle.id_;

    const dataType separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(saddle, scalars);
    dataType separatrixFunctionMaximum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const dcg::simplexId_t saddle2Id : separatricesSaddles[i]){
      if(inputTriangulation_->isTriangleOnBoundary(saddle2Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(dcg::Cell(2,saddle2Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMaximum=std::max(separatrixFunctionMaximum,
            discreteGradient_.scalarMax<dataType>(dcg::Cell(2,saddle2Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      for(const dcg::Cell& edge : separatricesGeometry[geometryId]){
        const dcg::simplexId_t edgeId=edge.id_;

        // Transform to dual : edge -> polygon
        std::vector<dcg::simplexId_t> polygon;
        getDualPolygon(edgeId, polygon);

        const dcg::simplexId_t vertexNumber=polygon.size();
        if(vertexNumber>2){
          sortDualPolygonVertices(polygon);

          // add the polygon
          separatrices2_cells[threadId].push_back(vertexNumber);
          for(dcg::simplexId_t k=0; k<vertexNumber; ++k){
            const dcg::simplexId_t tetraId=polygon[k];
            float point[3];
            discreteGradient_.getCellIncenter(dcg::Cell(3,tetraId), point);

            separatrices2_points[threadId].push_back(point[0]);
            separatrices2_points[threadId].push_back(point[1]);
            separatrices2_points[threadId].push_back(point[2]);
            separatrices2_cells[threadId].push_back(numberOfPoints[threadId]);
            ++numberOfPoints[threadId];
          }

          if(outputSeparatrices2_cells_sourceIds_)
            separatrices2_cells_sourceIds[threadId].push_back(saddleId);
          if(outputSeparatrices2_cells_separatrixIds_)
            separatrices2_cells_separatrixIds[threadId].push_back(separatrixIds[threadId]);
          if(outputSeparatrices2_cells_separatrixTypes_)
            separatrices2_cells_separatrixTypes[threadId].push_back(separatrixType);
          if(outputSeparatrices2_cells_separatrixFunctionMaxima)
            separatrices2_cells_separatrixFunctionMaxima[threadId].push_back(separatrixFunctionMaximum);
          if(outputSeparatrices2_cells_separatrixFunctionMinima)
            separatrices2_cells_separatrixFunctionMinima[threadId].push_back(separatrixFunctionMinimum);
          if(outputSeparatrices2_cells_separatrixFunctionDiffs)
            separatrices2_cells_separatrixFunctionDiffs[threadId].push_back(separatrixFunctionDiff);
          if(outputSeparatrices2_cells_isOnBoundary_)
            separatrices2_cells_isOnBoundary[threadId].push_back(isOnBoundary);
          ++numberOfCells[threadId];

          isFirst=false;
        }
      }
    }

    if(!isFirst)
      ++separatrixIds[threadId];
  }

  const dcg::simplexId_t oldPointSize=outputSeparatrices2_points_->size();
  const dcg::simplexId_t oldCellSize=outputSeparatrices2_cells_->size();
  const dcg::simplexId_t oldFieldSize=outputSeparatrices2_cells_sourceIds_->size();
  dcg::simplexId_t totalNumberOfPoints=0;
  dcg::simplexId_t totalNumberOfCells=0;
  {
    dcg::simplexId_t npoints=0;
    dcg::simplexId_t ncells=0;
    dcg::simplexId_t nnpoints=pointId;
    dcg::simplexId_t nncells=0;
    dcg::simplexId_t tmp_separatrixId=separatrixId;
    std::vector<dcg::simplexId_t> offsetPoints(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetCells(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetNPoints(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetNCells(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetSeparatrixIds(threadNumber_, 0);

    for(ThreadId i=0; i<threadNumber_; ++i){
      offsetPoints[i]=npoints;
      offsetCells[i]=ncells;
      offsetNPoints[i]=nnpoints;
      offsetNCells[i]=nncells;
      offsetSeparatrixIds[i]=tmp_separatrixId;

      npoints+=separatrices2_points[i].size();
      ncells+=separatrices2_cells[i].size();
      nnpoints+=numberOfPoints[i];
      nncells+=numberOfCells[i];
      tmp_separatrixId+=separatrixIds[i];

      totalNumberOfPoints+=numberOfPoints[i];
      totalNumberOfCells+=numberOfCells[i];
    }

    outputSeparatrices2_points_->resize(oldPointSize+npoints);
    outputSeparatrices2_cells_->resize(oldCellSize+ncells);
    if(outputSeparatrices2_cells_sourceIds_)
      outputSeparatrices2_cells_sourceIds_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixIds_)
      outputSeparatrices2_cells_separatrixIds_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixTypes_)
      outputSeparatrices2_cells_separatrixTypes_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionMaxima)
      outputSeparatrices2_cells_separatrixFunctionMaxima->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionMinima)
      outputSeparatrices2_cells_separatrixFunctionMinima->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionDiffs)
      outputSeparatrices2_cells_separatrixFunctionDiffs->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_isOnBoundary_)
      outputSeparatrices2_cells_isOnBoundary_->resize(oldFieldSize+totalNumberOfCells);
#pragma omp parallel for num_threads(threadNumber_)
    for(ThreadId i=0; i<threadNumber_; ++i){
      // reduce: points
      const dcg::simplexId_t tmp_npoints=separatrices2_points[i].size();
      const dcg::simplexId_t tmp_offsetPoints=offsetPoints[i];
      for(dcg::simplexId_t j=0; j<tmp_npoints; ++j)
        (*outputSeparatrices2_points_)[oldPointSize+tmp_offsetPoints+j]=separatrices2_points[i][j];

      // reduce: cells
      const dcg::simplexId_t tmp_ncells=separatrices2_cells[i].size();
      const dcg::simplexId_t tmp_offsetCells=offsetCells[i];
      for(dcg::simplexId_t j=0; j<tmp_ncells;){
        const dcg::simplexId_t cellSize=separatrices2_cells[i][j];
        (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j]=cellSize;
        for(dcg::simplexId_t k=0; k<cellSize; ++k){
          const dcg::simplexId_t tmp_pointId=separatrices2_cells[i][j+k+1];
          (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j+k+1]=offsetNPoints[i]+tmp_pointId;
        }
        j+=(cellSize+1);
      }

      // reduce: fields
      for(dcg::simplexId_t j=0; j<numberOfCells[i]; ++j){
        if(outputSeparatrices2_cells_sourceIds_)
          (*outputSeparatrices2_cells_sourceIds_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_sourceIds[i][j];
        if(outputSeparatrices2_cells_separatrixIds_)
          (*outputSeparatrices2_cells_separatrixIds_)[oldFieldSize+offsetNCells[i]+j]=
            offsetSeparatrixIds[i]+separatrices2_cells_separatrixIds[i][j];
        if(outputSeparatrices2_cells_separatrixTypes_)
          (*outputSeparatrices2_cells_separatrixTypes_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixTypes[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionMaxima)
          (*outputSeparatrices2_cells_separatrixFunctionMaxima)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionMaxima[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionMinima)
          (*outputSeparatrices2_cells_separatrixFunctionMinima)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionMinima[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionDiffs)
          (*outputSeparatrices2_cells_separatrixFunctionDiffs)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionDiffs[i][j];
        if(outputSeparatrices2_cells_isOnBoundary_)
          (*outputSeparatrices2_cells_isOnBoundary_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_isOnBoundary[i][j];
      }
    }
  }

  (*outputSeparatrices2_numberOfPoints_)+=totalNumberOfPoints;
  (*outputSeparatrices2_numberOfCells_)+=totalNumberOfCells;

  return 0;
}
#endif

template<typename dataType>
int ttk::MorseSmaleComplex3D::setAscendingSeparatrices2(const std::vector<Separatrix>& separatrices,
   const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
   const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const{
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices2_numberOfPoints_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfPoints is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_points_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_numberOfCells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfCells is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_cells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null." << std::endl;
    return -1;
  }
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input scalar field is null." << std::endl;
    return -1;
  }
#endif

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  dcg::simplexId_t pointId=(*outputSeparatrices2_numberOfPoints_);
  dcg::simplexId_t cellId=(*outputSeparatrices2_numberOfCells_);
  dcg::simplexId_t separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_ and
      outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  const dcg::simplexId_t numberOfCells=inputTriangulation_->getNumberOfCells();
  std::vector<dcg::simplexId_t> isVisited(numberOfCells, -1);

  const dcg::simplexId_t numberOfSeparatrices=separatrices.size();
  for(dcg::simplexId_t i=0; i<numberOfSeparatrices; ++i){
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle=separatrix.source_;
    const char separatrixType=1;
    const dcg::simplexId_t saddleId=saddle.id_;

    const dataType separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(saddle, scalars);
    dataType separatrixFunctionMaximum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const dcg::simplexId_t saddle2Id : separatricesSaddles[i]){
      if(inputTriangulation_->isTriangleOnBoundary(saddle2Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(dcg::Cell(2,saddle2Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMaximum=std::max(separatrixFunctionMaximum,
            discreteGradient_.scalarMax<dataType>(dcg::Cell(2,saddle2Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      for(const dcg::Cell& edge : separatricesGeometry[geometryId]){
        const dcg::simplexId_t edgeId=edge.id_;

        // Transform to dual : edge -> polygon
        std::vector<dcg::simplexId_t> polygon;
        getDualPolygon(edgeId, polygon);

        const dcg::simplexId_t vertexNumber=polygon.size();
        if(vertexNumber>2){
          sortDualPolygonVertices(polygon);

          // add the polygon
          outputSeparatrices2_cells_->push_back(vertexNumber);

          float point[3];
          for(dcg::simplexId_t i=0; i<vertexNumber; ++i){
            const dcg::simplexId_t tetraId=polygon[i];
            discreteGradient_.getCellIncenter(dcg::Cell(3,tetraId), point);

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

          if(outputSeparatrices2_cells_sourceIds_)
            outputSeparatrices2_cells_sourceIds_->push_back(saddleId);
          if(outputSeparatrices2_cells_separatrixIds_)
            outputSeparatrices2_cells_separatrixIds_->push_back(separatrixId);
          if(outputSeparatrices2_cells_separatrixTypes_)
            outputSeparatrices2_cells_separatrixTypes_->push_back(separatrixType);
          if(outputSeparatrices2_cells_separatrixFunctionMaxima)
            outputSeparatrices2_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
          if(outputSeparatrices2_cells_separatrixFunctionMinima)
            outputSeparatrices2_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
          if(outputSeparatrices2_cells_separatrixFunctionDiffs)
            outputSeparatrices2_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
          if(outputSeparatrices2_cells_isOnBoundary_)
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

#ifdef TTK_ENABLE_OPENMP      
template<typename dataType>
int ttk::MorseSmaleComplex3D::omp_setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
    const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
    const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const{
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices2_numberOfPoints_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfPoints is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_points_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_numberOfCells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfCells is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_cells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null." << std::endl;
    return -1;
  }
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input scalar field is null." << std::endl;
    return -1;
  }
#endif

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  dcg::simplexId_t pointId=(*outputSeparatrices2_numberOfPoints_);
  dcg::simplexId_t separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_ and
      outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  std::vector<dcg::simplexId_t> separatrixIds(threadNumber_, 0);
  std::vector<dcg::simplexId_t> numberOfPoints(threadNumber_, 0);
  std::vector<std::vector<float>> separatrices2_points(threadNumber_);
  std::vector<dcg::simplexId_t> numberOfCells(threadNumber_, 0);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells(threadNumber_);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells_sourceIds(threadNumber_);
  std::vector<std::vector<dcg::simplexId_t>> separatrices2_cells_separatrixIds(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_separatrixTypes(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMaxima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionMinima(threadNumber_);
  std::vector<std::vector<dataType>> separatrices2_cells_separatrixFunctionDiffs(threadNumber_);
  std::vector<std::vector<char>> separatrices2_cells_isOnBoundary(threadNumber_);

  const dcg::simplexId_t numberOfSeparatrices=separatrices.size();
#pragma omp parallel for num_threads(threadNumber_)
  for(dcg::simplexId_t i=0; i<numberOfSeparatrices; ++i){
    const ThreadId threadId=omp_get_thread_num();
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle=separatrix.source_;
    const char separatrixType=2;
    const dcg::simplexId_t saddleId=saddle.id_;

    const dataType separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(saddle, scalars);
    dataType separatrixFunctionMinimum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const dcg::simplexId_t saddle1Id : separatricesSaddles[i]){
      if(inputTriangulation_->isEdgeOnBoundary(saddle1Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(dcg::Cell(1,saddle1Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMinimum=std::min(separatrixFunctionMinimum,
            discreteGradient_.scalarMin<dataType>(dcg::Cell(1,saddle1Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      for(const dcg::Cell& cell : separatricesGeometry[geometryId]){
        const dcg::simplexId_t triangleId=cell.id_;

        separatrices2_cells[threadId].push_back(3);
        for(int k=0; k<3; ++k){
          dcg::simplexId_t vertexId;
          inputTriangulation_->getTriangleVertex(triangleId, k, vertexId);
          float point[3];
          inputTriangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);

          separatrices2_points[threadId].push_back(point[0]);
          separatrices2_points[threadId].push_back(point[1]);
          separatrices2_points[threadId].push_back(point[2]);
          separatrices2_cells[threadId].push_back(numberOfPoints[threadId]);
          ++numberOfPoints[threadId];
        }

        if(outputSeparatrices2_cells_sourceIds_)
          separatrices2_cells_sourceIds[threadId].push_back(saddleId);
        if(outputSeparatrices2_cells_separatrixIds_)
          separatrices2_cells_separatrixIds[threadId].push_back(separatrixIds[threadId]);
        if(outputSeparatrices2_cells_separatrixTypes_)
          separatrices2_cells_separatrixTypes[threadId].push_back(separatrixType);
        if(outputSeparatrices2_cells_separatrixFunctionMaxima)
          separatrices2_cells_separatrixFunctionMaxima[threadId].push_back(separatrixFunctionMaximum);
        if(outputSeparatrices2_cells_separatrixFunctionMinima)
          separatrices2_cells_separatrixFunctionMinima[threadId].push_back(separatrixFunctionMinimum);
        if(outputSeparatrices2_cells_separatrixFunctionDiffs)
          separatrices2_cells_separatrixFunctionDiffs[threadId].push_back(separatrixFunctionDiff);
        if(outputSeparatrices2_cells_isOnBoundary_)
          separatrices2_cells_isOnBoundary[threadId].push_back(isOnBoundary);
        ++numberOfCells[threadId];

        isFirst=false;
      }
    }

    if(!isFirst)
      ++separatrixIds[threadId];
  }

  const dcg::simplexId_t oldPointSize=outputSeparatrices2_points_->size();
  const dcg::simplexId_t oldCellSize=outputSeparatrices2_cells_->size();
  const dcg::simplexId_t oldFieldSize=outputSeparatrices2_cells_sourceIds_->size();
  dcg::simplexId_t totalNumberOfPoints=0;
  dcg::simplexId_t totalNumberOfCells=0;
  {
    dcg::simplexId_t npoints=0;
    dcg::simplexId_t ncells=0;
    dcg::simplexId_t nnpoints=pointId;
    dcg::simplexId_t nncells=0;
    dcg::simplexId_t tmp_separatrixId=separatrixId;
    std::vector<dcg::simplexId_t> offsetPoints(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetCells(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetNPoints(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetNCells(threadNumber_, 0);
    std::vector<dcg::simplexId_t> offsetSeparatrixIds(threadNumber_, 0);

    for(ThreadId i=0; i<threadNumber_; ++i){
      offsetPoints[i]=npoints;
      offsetCells[i]=ncells;
      offsetNPoints[i]=nnpoints;
      offsetNCells[i]=nncells;
      offsetSeparatrixIds[i]=tmp_separatrixId;

      npoints+=separatrices2_points[i].size();
      ncells+=separatrices2_cells[i].size();
      nnpoints+=numberOfPoints[i];
      nncells+=numberOfCells[i];
      tmp_separatrixId+=separatrixIds[i];

      totalNumberOfPoints+=numberOfPoints[i];
      totalNumberOfCells+=numberOfCells[i];
    }

    outputSeparatrices2_points_->resize(oldPointSize+npoints);
    outputSeparatrices2_cells_->resize(oldCellSize+ncells);
    if(outputSeparatrices2_cells_sourceIds_)
      outputSeparatrices2_cells_sourceIds_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixIds_)
      outputSeparatrices2_cells_separatrixIds_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixTypes_)
      outputSeparatrices2_cells_separatrixTypes_->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionMaxima)
      outputSeparatrices2_cells_separatrixFunctionMaxima->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionMinima)
      outputSeparatrices2_cells_separatrixFunctionMinima->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_separatrixFunctionDiffs)
      outputSeparatrices2_cells_separatrixFunctionDiffs->resize(oldFieldSize+totalNumberOfCells);
    if(outputSeparatrices2_cells_isOnBoundary_)
      outputSeparatrices2_cells_isOnBoundary_->resize(oldFieldSize+totalNumberOfCells);
#pragma omp parallel for num_threads(threadNumber_)
    for(ThreadId i=0; i<threadNumber_; ++i){
      // reduce: points
      const dcg::simplexId_t tmp_npoints=separatrices2_points[i].size();
      const dcg::simplexId_t tmp_offsetPoints=offsetPoints[i];
      for(dcg::simplexId_t j=0; j<tmp_npoints; ++j)
        (*outputSeparatrices2_points_)[oldPointSize+tmp_offsetPoints+j]=separatrices2_points[i][j];

      // reduce: cells
      const dcg::simplexId_t tmp_ncells=separatrices2_cells[i].size();
      const dcg::simplexId_t tmp_offsetCells=offsetCells[i];
      for(dcg::simplexId_t j=0; j<tmp_ncells;){
        const dcg::simplexId_t cellSize=separatrices2_cells[i][j];
        (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j]=cellSize;
        for(dcg::simplexId_t k=0; k<cellSize; ++k){
          const dcg::simplexId_t tmp_pointId=separatrices2_cells[i][j+k+1];
          (*outputSeparatrices2_cells_)[oldCellSize+tmp_offsetCells+j+k+1]=offsetNPoints[i]+tmp_pointId;
        }
        j+=(cellSize+1);
      }

      // reduce: fields
      for(dcg::simplexId_t j=0; j<numberOfCells[i]; ++j){
        if(outputSeparatrices2_cells_sourceIds_)
          (*outputSeparatrices2_cells_sourceIds_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_sourceIds[i][j];
        if(outputSeparatrices2_cells_separatrixIds_)
          (*outputSeparatrices2_cells_separatrixIds_)[oldFieldSize+offsetNCells[i]+j]=
            offsetSeparatrixIds[i]+separatrices2_cells_separatrixIds[i][j];
        if(outputSeparatrices2_cells_separatrixTypes_)
          (*outputSeparatrices2_cells_separatrixTypes_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixTypes[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionMaxima)
          (*outputSeparatrices2_cells_separatrixFunctionMaxima)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionMaxima[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionMinima)
          (*outputSeparatrices2_cells_separatrixFunctionMinima)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionMinima[i][j];
        if(outputSeparatrices2_cells_separatrixFunctionDiffs)
          (*outputSeparatrices2_cells_separatrixFunctionDiffs)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_separatrixFunctionDiffs[i][j];
        if(outputSeparatrices2_cells_isOnBoundary_)
          (*outputSeparatrices2_cells_isOnBoundary_)[oldFieldSize+offsetNCells[i]+j]=
            separatrices2_cells_isOnBoundary[i][j];
      }
    }
  }

  (*outputSeparatrices2_numberOfPoints_)+=totalNumberOfPoints;
  (*outputSeparatrices2_numberOfCells_)+=totalNumberOfCells;

  return 0;
}
#endif

template<typename dataType>
int ttk::MorseSmaleComplex3D::setDescendingSeparatrices2(const std::vector<Separatrix>& separatrices,
   const std::vector<std::vector<dcg::Cell>>& separatricesGeometry,
   const std::vector<std::set<dcg::simplexId_t>>& separatricesSaddles) const{
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputSeparatrices2_numberOfPoints_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfPoints is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_points_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to points is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_numberOfCells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to numberOfCells is null." << std::endl;
    return -1;
  }
  if(!outputSeparatrices2_cells_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to cells is null." << std::endl;
    return -1;
  }
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] 2-separatrices pointer to the input scalar field is null." << std::endl;
    return -1;
  }
#endif

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMaxima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMaxima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionMinima=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionMinima_);
  std::vector<dataType>* outputSeparatrices2_cells_separatrixFunctionDiffs=
    static_cast<std::vector<dataType>*>(outputSeparatrices2_cells_separatrixFunctionDiffs_);

  dcg::simplexId_t pointId=(*outputSeparatrices2_numberOfPoints_);
  dcg::simplexId_t cellId=(*outputSeparatrices2_numberOfCells_);
  dcg::simplexId_t separatrixId=0;
  if(outputSeparatrices2_cells_separatrixIds_ and
      outputSeparatrices2_cells_separatrixIds_->size()){
    separatrixId=*std::max_element(outputSeparatrices2_cells_separatrixIds_->begin(),
        outputSeparatrices2_cells_separatrixIds_->end())+1;
  }

  const dcg::simplexId_t numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::vector<dcg::simplexId_t> isVisited(numberOfVertices, -1);

  const dcg::simplexId_t numberOfSeparatrices=separatrices.size();
  for(dcg::simplexId_t i=0; i<numberOfSeparatrices; ++i){
    const Separatrix& separatrix=separatrices[i];
    if(!separatrix.isValid_) continue;
    if(!separatrix.geometry_.size()) continue;

    const dcg::Cell& saddle=separatrix.source_;
    const char separatrixType=2;
    const dcg::simplexId_t saddleId=saddle.id_;

    const dataType separatrixFunctionMaximum=discreteGradient_.scalarMax<dataType>(saddle, scalars);
    dataType separatrixFunctionMinimum{};

    // get separatrix infos
    char isOnBoundary{};
    bool isFirst=true;
    for(const dcg::simplexId_t saddle1Id : separatricesSaddles[i]){
      if(inputTriangulation_->isEdgeOnBoundary(saddle1Id))
        ++isOnBoundary;

      if(isFirst){
        separatrixFunctionMinimum=discreteGradient_.scalarMin<dataType>(dcg::Cell(1,saddle1Id), scalars);
        isFirst=false;
      }
      else{
        separatrixFunctionMinimum=std::min(separatrixFunctionMinimum,
            discreteGradient_.scalarMin<dataType>(dcg::Cell(1,saddle1Id), scalars));
      }
    }

    const dataType separatrixFunctionDiff=separatrixFunctionMaximum-separatrixFunctionMinimum;

    isFirst=true;
    for(const dcg::simplexId_t geometryId : separatrix.geometry_){
      for(const dcg::Cell& cell : separatricesGeometry[geometryId]){
        const dcg::simplexId_t triangleId=cell.id_;

        outputSeparatrices2_cells_->push_back(3);
        float point[3];
        for(int k=0; k<3; ++k){
          dcg::simplexId_t vertexId;
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
        if(outputSeparatrices2_cells_sourceIds_)
          outputSeparatrices2_cells_sourceIds_->push_back(saddleId);
        if(outputSeparatrices2_cells_separatrixIds_)
          outputSeparatrices2_cells_separatrixIds_->push_back(separatrixId);
        if(outputSeparatrices2_cells_separatrixTypes_)
          outputSeparatrices2_cells_separatrixTypes_->push_back(separatrixType);
        if(outputSeparatrices2_cells_separatrixFunctionMaxima)
          outputSeparatrices2_cells_separatrixFunctionMaxima->push_back(separatrixFunctionMaximum);
        if(outputSeparatrices2_cells_separatrixFunctionMinima)
          outputSeparatrices2_cells_separatrixFunctionMinima->push_back(separatrixFunctionMinimum);
        if(outputSeparatrices2_cells_separatrixFunctionDiffs)
          outputSeparatrices2_cells_separatrixFunctionDiffs->push_back(separatrixFunctionDiff);
        if(outputSeparatrices2_cells_isOnBoundary_)
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

template<typename dataType, typename idType>
int ttk::MorseSmaleComplex3D::execute(){
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField_){
    std::cerr << "[MorseSmaleComplex3D] Error: input scalar field pointer is null." << std::endl;
    return -1;
  }

  if(!inputOffsets_){
    std::cerr << "[MorseSmaleComplex3D] Error: input offset field pointer is null." << std::endl;
    return -1;
  }
#endif
  Timer t;

  // nullptr_t is implicitly convertible and comparable to any pointer type
  // or pointer-to-member type.
  dcg::simplexId_t* ascendingManifold=static_cast<dcg::simplexId_t*>(outputAscendingManifold_);
  dcg::simplexId_t* descendingManifold=static_cast<dcg::simplexId_t*>(outputDescendingManifold_);
  dcg::simplexId_t* morseSmaleManifold=static_cast<dcg::simplexId_t*>(outputMorseSmaleManifold_);

  discreteGradient_.setThreadNumber(threadNumber_);
  discreteGradient_.setDebugLevel(debugLevel_);
  {
    Timer tmp;
    discreteGradient_.buildGradient<dataType,idType>();
    discreteGradient_.buildGradient2<dataType,idType>();
    discreteGradient_.buildGradient3<dataType,idType>();

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Discrete gradient overall computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }
  discreteGradient_.reverseGradient<dataType,idType>();

  std::vector<dcg::Cell> criticalPoints;
  discreteGradient_.getCriticalPoints(criticalPoints);

  // 1-separatrices
  if(ComputeDescendingSeparatrices1){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getDescendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Descending 1-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices1){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getAscendingSeparatrices1(criticalPoints, separatrices, separatricesGeometry);
    setSeparatrices1<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Ascending 1-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  // saddle-connectors
  if(ComputeSaddleConnectors){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    getSaddleConnectors(criticalPoints, separatrices, separatricesGeometry);
    setSaddleConnectors<dataType>(separatrices, separatricesGeometry);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Saddle connectors computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  // 2-separatrices
  if(ComputeDescendingSeparatrices2){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<dcg::simplexId_t>> separatricesSaddles;
    getDescendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
#ifdef TTK_ENABLE_OPENMP
    if(PrioritizeSpeedOverMemory)
      omp_setDescendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
    else
#endif
    setDescendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Descending 2-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(ComputeAscendingSeparatrices2){
    Timer tmp;
    std::vector<Separatrix> separatrices;
    std::vector<std::vector<dcg::Cell>> separatricesGeometry;
    std::vector<std::set<dcg::simplexId_t>> separatricesSaddles;
    getAscendingSeparatrices2(criticalPoints, separatrices, separatricesGeometry, separatricesSaddles);
#ifdef TTK_ENABLE_OPENMP
    if(PrioritizeSpeedOverMemory)
      omp_setAscendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);
    else
#endif
    setAscendingSeparatrices2<dataType>(separatrices, separatricesGeometry, separatricesSaddles);

    {
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Ascending 2-separatrices computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  std::vector<dcg::simplexId_t> maxSeeds;
  {
    Timer tmp;

    dcg::simplexId_t numberOfMaxima{};
    dcg::simplexId_t numberOfMinima{};

    if(ascendingManifold)
      setAscendingSegmentation(criticalPoints, maxSeeds, ascendingManifold, numberOfMaxima);

    if(descendingManifold)
      setDescendingSegmentation(criticalPoints, descendingManifold, numberOfMinima);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      setFinalSegmentation(numberOfMaxima, numberOfMinima, ascendingManifold, descendingManifold, morseSmaleManifold);

    if(ascendingManifold or descendingManifold){
      std::stringstream msg;
      msg << "[MorseSmaleComplex3D] Segmentation computed in "
        << tmp.getElapsedTime() << " s."
        << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  if(outputCriticalPoints_numberOfPoints_ and outputSeparatrices1_points_){
    if(ascendingManifold and descendingManifold)
      discreteGradient_.setAugmentedCriticalPoints<dataType,idType>(criticalPoints,
          maxSeeds,
          ascendingManifold,
          descendingManifold);
    else
      discreteGradient_.setCriticalPoints<dataType,idType>(criticalPoints);
  }

  {
    const dcg::simplexId_t numberOfVertices=inputTriangulation_->getNumberOfVertices();
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

template<typename dataType, typename idType>
int ttk::MorseSmaleComplex3D::computePersistencePairs(const std::vector<std::tuple<dcg::simplexId_t, dcg::simplexId_t, dataType>>& JTPairs,
    const std::vector<std::tuple<dcg::simplexId_t, dcg::simplexId_t, dataType>>& STPairs,
    std::vector<std::tuple<dcg::simplexId_t,dcg::simplexId_t,dataType>>& pl_saddleSaddlePairs){
  dataType* scalars=static_cast<dataType*>(inputScalarField_);
  const dcg::simplexId_t numberOfVertices=inputTriangulation_->getNumberOfVertices();

  // get original list of critical points
  std::vector<std::pair<dcg::simplexId_t,char>> pl_criticalPoints;
  {
    const dcg::simplexId_t* const offsets=static_cast<dcg::simplexId_t*>(inputOffsets_);
    std::vector<dcg::simplexId_t> sosOffsets(numberOfVertices);
    for(dcg::simplexId_t i=0; i<numberOfVertices; ++i)
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
    const dcg::simplexId_t v0=std::get<0>(i);
    const dcg::simplexId_t v1=std::get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }
  for(const auto& i : STPairs){
    const dcg::simplexId_t v0=std::get<0>(i);
    const dcg::simplexId_t v1=std::get<1>(i);
    isAccepted[v0]=true;
    isAccepted[v1]=true;
  }

  // filter the critical points according to the filtering list and boundary condition
  std::vector<char> isSaddle1(numberOfVertices, false);
  std::vector<char> isSaddle2(numberOfVertices, false);
  std::vector<std::pair<dcg::simplexId_t,char>> pl_filteredCriticalPoints;
  for(const auto& i : pl_criticalPoints){
    const dcg::simplexId_t vertexId=i.first;
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

  std::vector<std::tuple<dcg::Cell,dcg::Cell>> dmt_pairs;
  {
    // simplify to be PL-conformant
    discreteGradient_.setDebugLevel(debugLevel_);
    discreteGradient_.setThreadNumber(threadNumber_);
    discreteGradient_.setReverseSaddleMaximumConnection(true);
    discreteGradient_.setReverseSaddleSaddleConnection(true);
    discreteGradient_.setCollectPersistencePairs(false);
    discreteGradient_.buildGradient<dataType,idType>();
    discreteGradient_.buildGradient2<dataType,idType>();
    discreteGradient_.buildGradient3<dataType,idType>();
    discreteGradient_.reverseGradient<dataType,idType>(pl_criticalPoints);

    // collect saddle-saddle connections
    discreteGradient_.setReverseSaddleMaximumConnection(true);
    discreteGradient_.setCollectPersistencePairs(true);
    discreteGradient_.setOutputPersistencePairs(&dmt_pairs);
    discreteGradient_.reverseGradient<dataType,idType>(pl_filteredCriticalPoints);
  }

  // transform DMT pairs into PL pairs
  for(const auto& i : dmt_pairs){
    const dcg::Cell& saddle1=std::get<0>(i);
    const dcg::Cell& saddle2=std::get<1>(i);

    dcg::simplexId_t v0=-1;
    for(dcg::simplexId_t j=0; j<2; ++j){
      dcg::simplexId_t vertexId;
      inputTriangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

      if(isSaddle1[vertexId]){
        v0=vertexId;
        break;
      }
    }
    if(v0==-1){
      dataType scalar{};
      for(int j=0; j<2; ++j){
        dcg::simplexId_t vertexId;
        inputTriangulation_->getEdgeVertex(saddle1.id_,j,vertexId);
        const dataType vertexScalar=scalars[vertexId];

        if(!j or scalar>vertexScalar){
          v0=vertexId;
          scalar=vertexScalar;
        }
      }
    }

    dcg::simplexId_t v1=-1;
    for(int j=0; j<3; ++j){
      dcg::simplexId_t vertexId;
      inputTriangulation_->getTriangleVertex(saddle2.id_, j, vertexId);

      if(isSaddle2[vertexId]){
        v1=vertexId;
        break;
      }
    }
    if(v1==-1){
      dataType scalar{};
      for(int j=0; j<3; ++j){
        dcg::simplexId_t vertexId;
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
