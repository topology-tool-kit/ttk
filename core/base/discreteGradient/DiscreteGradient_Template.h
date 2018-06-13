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

#ifndef _DISCRETEGRADIENT_TPL_H
#define _DISCRETEGRADIENT_TPL_H

#include<DiscreteGradient.h>

template <typename dataType>
dataType DiscreteGradient::scalarMax(const Cell& cell, const dataType*
const scalars) const{
  dataType scalar{};

  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          simplexId_t vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          simplexId_t vertexId;
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
          simplexId_t vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          simplexId_t vertexId;
          inputTriangulation_->getTriangleVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar<vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 3:
        for(int i=0; i<4; ++i){
          simplexId_t vertexId;
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
dataType DiscreteGradient::scalarMin(const Cell& cell, const dataType*
const scalars) const{
  dataType scalar{};

  if(dimensionality_==2){
    switch(cell.dim_){
      case 0:
        scalar=scalars[cell.id_];
        break;

      case 1:
        for(int i=0; i<2; ++i){
          simplexId_t vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          simplexId_t vertexId;
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
          simplexId_t vertexId;
          inputTriangulation_->getEdgeVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 2:
        for(int i=0; i<3; ++i){
          simplexId_t vertexId;
          inputTriangulation_->getTriangleVertex(cell.id_,i,vertexId);
          const dataType vertexScalar=scalars[vertexId];

          if(!i or scalar>vertexScalar)
            scalar=vertexScalar;
        }
        break;

      case 3:
        for(int i=0; i<4; ++i){
          simplexId_t vertexId;
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
dataType DiscreteGradient::getPersistence(const Cell& up, const Cell& down,
                                          const dataType* const scalars) const{
  return scalarMax<dataType>(up,scalars)-scalarMin<dataType>(down,scalars);
}

template <typename dataType>
bool DiscreteGradient::isHigherThan(const simplexId_t vertexA,
                                    const simplexId_t vertexB,
                                    const dataType* const scalars,
                                    const simplexId_t* const offsets) const{
  if(scalars[vertexA] != scalars[vertexB]) return
      scalars[vertexA]>scalars[vertexB];
  else return offsets[vertexA]>offsets[vertexB];
}

template <typename dataType>
bool DiscreteGradient::isLowerThan(const simplexId_t vertexA,
                                   const simplexId_t vertexB,
                                   const dataType* const scalars,
                                   const simplexId_t* const offsets) const{
  if(scalars[vertexA] != scalars[vertexB]) return
      scalars[vertexA]<scalars[vertexB];
  else return offsets[vertexA]<offsets[vertexB];
}

template <typename dataType>
int DiscreteGradient::cellMax(const int cellDim,
                              const simplexId_t cellA,
                              const simplexId_t cellB,
                              const dataType* const scalars,
                              const simplexId_t* const offsets) const{
  const int vertexNumber=cellDim+1;
  simplexId_t vsetA[4];
  simplexId_t vsetB[4];

  const auto sosGreaterThan=[&scalars,&offsets](const simplexId_t a, const simplexId_t b){
    if(scalars[a] != scalars[b]) return scalars[a]>scalars[b];
    else return offsets[a]>offsets[b];
  };

  if(dimensionality_==2){
    switch(cellDim){
      case 0:
        return (isHigherThan<dataType>(cellA, cellB, scalars, offsets))? cellA
                                                                       :
               cellB;

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

    std::sort(vsetA, vsetA+vertexNumber, sosGreaterThan);
    std::sort(vsetB, vsetB+vertexNumber, sosGreaterThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isHigherThan<dataType>(vsetA[k], vsetB[k],
                                          scalars,offsets))?
                  cellA : cellB;
    }
  }
  else if(dimensionality_==3){
    switch(cellDim){
      case 0:
        return (isHigherThan<dataType>(cellA, cellB, scalars, offsets)? cellA :
                cellB);

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

    std::sort(vsetA, vsetA+vertexNumber, sosGreaterThan);
    std::sort(vsetB, vsetB+vertexNumber, sosGreaterThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isHigherThan<dataType>(vsetA[k], vsetB[k], scalars,offsets)?
                   cellA : cellB);
    }
  }

  return -1;
}

template <typename dataType>
int DiscreteGradient::cellMin(const int cellDim,
                              const simplexId_t cellA,
                              const simplexId_t cellB,
                              const dataType* const scalars,
                              const simplexId_t* const offsets) const{
  const int vertexNumber=cellDim+1;
  simplexId_t vsetA[4];
  simplexId_t vsetB[4];

  const auto sosLowerThan=[&scalars,&offsets](const simplexId_t a, const simplexId_t b){
    if(scalars[a] != scalars[b]) return scalars[a]<scalars[b];
    else return offsets[a]<offsets[b];
  };

  if(dimensionality_==2){
    switch(cellDim){
      case 0:
        return (isLowerThan<dataType>(cellA, cellB, scalars, offsets))? cellA :
               cellB;

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

    std::sort(vsetA, vsetA+vertexNumber, sosLowerThan);
    std::sort(vsetB, vsetB+vertexNumber, sosLowerThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isLowerThan<dataType>(vsetA[k], vsetB[k], scalars,offsets))?
                  cellA : cellB;
    }
  }
  else if(dimensionality_==3){
    switch(cellDim){
      case 0:
        return (isLowerThan<dataType>(cellA, cellB, scalars, offsets))? cellA :
               cellB;

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

    std::sort(vsetA, vsetA+vertexNumber, sosLowerThan);
    std::sort(vsetB, vsetB+vertexNumber, sosLowerThan);
    for(int k=0; k<vertexNumber; ++k){
      if(vsetA[k]==vsetB[k]) continue;
      else return (isLowerThan<dataType>(vsetA[k], vsetB[k], scalars,offsets))?
                  cellA : cellB;
    }
  }

  return -1;
}

template <typename dataType>
int DiscreteGradient::g0(const int cellDim,
                         const simplexId_t cellId,
                         const dataType* const scalars,
                         const simplexId_t* const offsets) const{
  simplexId_t facet0;
  simplexId_t facet1;
  simplexId_t facetMax{-1};

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
                                const simplexId_t cellId,
                                const dataType* const scalars,
                                const simplexId_t* const offsets) const{
  simplexId_t facetMax{-1};
  simplexId_t facetMaxSecond{-1};
  simplexId_t facets[4];

  if(dimensionality_==2){
    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facets[0]);
        inputTriangulation_->getEdgeVertex(cellId, 1, facets[1]);
        facetMaxSecond=cellMin<dataType>(0, facets[0], facets[1],
                                         scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getCellEdge(cellId,0,facets[0]);
        inputTriangulation_->getCellEdge(cellId,1,facets[1]);
        facetMax=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);

        inputTriangulation_->getCellEdge(cellId,2,facets[2]);
        facetMax=cellMax<dataType>(1, facets[2], facetMax, scalars,offsets);

        if(facetMax==facets[0])
          facetMaxSecond=cellMax<dataType>(1, facets[1], facets[2],
                                           scalars,offsets);
        else if(facetMax==facets[1])
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[2],
                                           scalars,offsets);
        else
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[1],
                                           scalars,offsets);
        break;
    }
  }
  else if(dimensionality_==3){
    switch(cellDim){
      case 1:
        inputTriangulation_->getEdgeVertex(cellId, 0, facets[0]);
        inputTriangulation_->getEdgeVertex(cellId, 1, facets[1]);
        facetMaxSecond=cellMin<dataType>(0, facets[0], facets[1],
                                         scalars,offsets);
        break;

      case 2:
        inputTriangulation_->getTriangleEdge(cellId,0,facets[0]);
        inputTriangulation_->getTriangleEdge(cellId,1,facets[1]);
        facetMax=cellMax<dataType>(1, facets[0], facets[1], scalars,offsets);

        inputTriangulation_->getTriangleEdge(cellId,2,facets[2]);
        facetMax=cellMax<dataType>(1, facets[2], facetMax, scalars,offsets);

        if(facetMax==facets[0])
          facetMaxSecond=cellMax<dataType>(1, facets[1], facets[2],
                                           scalars,offsets);
        else if(facetMax==facets[1])
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[2],
                                           scalars,offsets);
        else
          facetMaxSecond=cellMax<dataType>(1, facets[0], facets[1],
                                           scalars,offsets);
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
          else
          if(facets[i]==cellMax<dataType>(2,facets[i],facetMaxSecond,scalars,offsets)){
            facetMaxSecond=facets[i];
          }
        }
        break;
    }
  }

  return facetMaxSecond;
}

template <typename dataType>
int DiscreteGradient::g0_third(const int cellDim,
                               const simplexId_t cellId,
                               const dataType* const scalars,
                               const simplexId_t* const offsets) const{
  simplexId_t facetMaxThird{-1};

  if(dimensionality_==3){
    simplexId_t facetMax{-1};
    simplexId_t facetMaxSecond{-1};
    simplexId_t facetMin{-1};
    simplexId_t facets[4];

    switch(cellDim){
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

        facetMin=cellMin<dataType>(2,facets[0],facets[1],scalars,offsets);

        for(int i=2; i<4; i++){
          if(facets[i]==cellMax<dataType>(2,facets[i],facetMax,scalars,offsets)){
            facetMaxSecond=facetMax;
            facetMax=facets[i];
          }
          else
          if(facets[i]==cellMax<dataType>(2,facets[i],facetMaxSecond,scalars,offsets)){
            facetMaxSecond=facets[i];
          }

          facetMin=cellMin<dataType>(2,facets[i],facetMin,scalars,offsets);
        }

        for(int i=0; i<4; i++){
          if(facets[i]!=facetMax and facets[i]!=facetMaxSecond and
             facets[i]!=facetMin){
            facetMaxThird=facets[i];
            break;
          }
        }
        break;
    }
  }

  return facetMaxThird;
}

template <typename dataType>
int DiscreteGradient::assignGradient(const int alphaDim,
                                     const dataType* const scalars,
                                     const simplexId_t* const offsets,
                                     std::vector<std::vector<simplexId_t>>& gradient) const{
  const int betaDim=alphaDim+1;
  const simplexId_t alphaNumber=gradient[alphaDim].size();

  if(dimensionality_==2){
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
    for(simplexId_t alpha=0; alpha<alphaNumber; ++alpha){
      simplexId_t betaNumber{};
      switch(alphaDim){
        case 0: betaNumber=inputTriangulation_->getVertexEdgeNumber(alpha);
          break;
        case 1: betaNumber=inputTriangulation_->getEdgeStarNumber(alpha); break;
      }
      simplexId_t gamma{-1};
      for(simplexId_t k=0; k<betaNumber; ++k){
        simplexId_t beta;
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
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
    for(simplexId_t alpha=0; alpha<alphaNumber; ++alpha){
      simplexId_t betaNumber{};
      switch(alphaDim){
        case 0: betaNumber=inputTriangulation_->getVertexEdgeNumber(alpha);
          break;
        case 1: betaNumber=inputTriangulation_->getEdgeTriangleNumber(alpha);
          break;
        case 2: betaNumber=inputTriangulation_->getTriangleStarNumber(alpha);
          break;
      }
      simplexId_t gamma{-1};
      for(simplexId_t k=0; k<betaNumber; ++k){
        simplexId_t beta;
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
                                      const simplexId_t* const offsets,
                                      std::vector<std::vector<simplexId_t>>& gradient) const{
  if(alphaDim>0){
    const int betaDim=alphaDim+1;
    const simplexId_t alphaNumber=gradient[alphaDim].size();

    if(dimensionality_==2){
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
      for(simplexId_t alpha=0; alpha<alphaNumber; ++alpha){
        // alpha must be unpaired
        if(gradient[alphaDim][alpha]==-1){
          simplexId_t betaNumber{};
          switch(alphaDim){
            case 1: betaNumber=inputTriangulation_->getEdgeStarNumber(alpha);
              break;
          }
          simplexId_t gamma{-1};
          for(simplexId_t k=0; k<betaNumber; ++k){
            simplexId_t beta;
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
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
      for(simplexId_t alpha=0; alpha<alphaNumber; ++alpha){
        // alpha must be unpaired
        if(gradient[alphaDim][alpha]==-1){
          simplexId_t betaNumber{};
          switch(alphaDim){
            case 1:
              betaNumber=inputTriangulation_->getEdgeTriangleNumber(alpha); break;
            case 2:
              betaNumber=inputTriangulation_->getTriangleStarNumber(alpha); break;
          }
          simplexId_t gamma{-1};
          for(simplexId_t k=0; k<betaNumber; ++k){
            simplexId_t beta;
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
int DiscreteGradient::assignGradient3(const int alphaDim,
                                      const dataType* const scalars,
                                      const simplexId_t* const offsets,
                                      std::vector<std::vector<simplexId_t>>& gradient) const{
  if(alphaDim>0){
    const int betaDim=alphaDim+1;
    const simplexId_t alphaNumber=gradient[alphaDim].size();

    if(dimensionality_==3){
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
      for(simplexId_t alpha=0; alpha<alphaNumber; ++alpha){
        // alpha must be unpaired
        if(gradient[alphaDim][alpha]==-1){
          simplexId_t betaNumber{};
          switch(alphaDim){
            case 2:
              betaNumber=inputTriangulation_->getTriangleStarNumber(alpha); break;
          }
          simplexId_t gamma{-1};
          for(simplexId_t k=0; k<betaNumber; ++k){
            simplexId_t beta;
            switch(alphaDim){
              case 2: inputTriangulation_->getTriangleStar(alpha,k,beta); break;
            }
            // take beta such that alpha is the second highest facet of beta
            if(alpha==g0_third<dataType>(betaDim,beta,scalars,offsets)){
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

  const simplexId_t* const offsets=static_cast<simplexId_t*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int numberOfDimensions=getNumberOfDimensions();

  // init number of cells by dimension
  std::vector<simplexId_t> numberOfCells(numberOfDimensions);
  for(int i=0; i<numberOfDimensions; ++i)
    numberOfCells[i]=getNumberOfCells(i);

  dmtMax2PL_.clear();
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
    std::stringstream msg;
    msg << "[DiscreteGradient] Data-set (" << numberOfVertices_
        << " points) processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::buildGradient2(){
  Timer t;

  const simplexId_t* const offsets=static_cast<simplexId_t*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  for(int i=1; i<dimensionality_; ++i)
    assignGradient2<dataType>(i, scalars, offsets, gradient_[i]);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Data-set (" << numberOfVertices_
        << " points) post-processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::buildGradient3(){
  Timer t;

  const simplexId_t* const offsets=static_cast<simplexId_t*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  for(int i=2; i<dimensionality_; ++i)
    assignGradient3<dataType>(i, scalars, offsets, gradient_[i]);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Data-set (" << numberOfVertices_
        << " points) post-processed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints(const std::vector<Cell>&
criticalPoints) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputCriticalPoints_points_cellScalars=
    static_cast<std::vector<dataType>*>(outputCriticalPoints_points_cellScalars_);

  (*outputCriticalPoints_numberOfPoints_)=0;

  const int numberOfDimensions=getNumberOfDimensions();
  std::vector<int> numberOfCriticalPointsByDimension(numberOfDimensions,0);

  // for all critical cells
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& cell=criticalPoints[i];
    const int cellDim=cell.dim_;
    const simplexId_t cellId=cell.id_;
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
    if(dmtMax2PL_.size()){
      if(cellDim==0)
        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(cellId);
      else if(cellDim==dimensionality_)

        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(dmtMax2PL_[cellId]);
      else
        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(-1);
    }
    else
      outputCriticalPoints_points_PLVertexIdentifiers_->push_back(-1);

    (*outputCriticalPoints_numberOfPoints_)++;
  }

  {
    std::stringstream msg;
    for(int i=0; i<numberOfDimensions; ++i)
      msg << "[DiscreteGradient] "
          << numberOfCriticalPointsByDimension[i] << " " << i << "-cell(s)."
          << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::setCriticalPoints() const{
  std::vector<Cell> criticalPoints;
  getCriticalPoints(criticalPoints);

  setCriticalPoints<dataType>(criticalPoints);

  return 0;
}

template <typename dataType>
int DiscreteGradient::setAugmentedCriticalPoints(const std::vector<Cell>& criticalPoints,
                                                 std::vector<simplexId_t>& maxSeeds,
                                                 simplexId_t* ascendingManifold,
                                                 simplexId_t* descendingManifold) const{
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);
  std::vector<dataType>* outputCriticalPoints_points_cellScalars=
    static_cast<std::vector<dataType>*>(outputCriticalPoints_points_cellScalars_);
  (*outputCriticalPoints_numberOfPoints_)=0;

  const int numberOfDimensions=getNumberOfDimensions();
  std::vector<int> numberOfCriticalPointsByDimension(numberOfDimensions,0);

  // for all critical cells
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& cell=criticalPoints[i];
    const int cellDim=cell.dim_;
    const simplexId_t cellId=cell.id_;
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
    if(dmtMax2PL_.size()){
      if(cellDim==0)
        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(cellId);
      else if(cellDim==dimensionality_)

        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(dmtMax2PL_[cellId]);
      else
        outputCriticalPoints_points_PLVertexIdentifiers_->push_back(-1);
    }
    else
      outputCriticalPoints_points_PLVertexIdentifiers_->push_back(-1);

    simplexId_t manifoldSize=0;
    if(cellDim==0){
      const simplexId_t seedId=descendingManifold[cellId];
      manifoldSize=std::count(descendingManifold,
                              descendingManifold+numberOfVertices_, seedId);
    }
    else if(cellDim==dimensionality_){
      auto ite=std::find(maxSeeds.begin(), maxSeeds.end(), cellId);
      if(ite!=maxSeeds.end()){
        const simplexId_t seedId=std::distance(maxSeeds.begin(), ite);
        manifoldSize=std::count(ascendingManifold,
                                ascendingManifold+numberOfVertices_, seedId);
      }
    }
    outputCriticalPoints_points_manifoldSize_->push_back(manifoldSize);

    (*outputCriticalPoints_numberOfPoints_)++;
  }

  {
    std::stringstream msg;
    for(int i=0; i<numberOfDimensions; ++i)
      msg << "[DiscreteGradient] " << numberOfCriticalPointsByDimension[i]
          << " " << i << "-cell(s)." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::getRemovableMaxima(const std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                         const bool allowBoundary,
                                         std::vector<char>& isRemovableMaximum,
                                         std::vector<simplexId_t>& pl2dmt_maximum){
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  const simplexId_t numberOfCells=inputTriangulation_->getNumberOfCells();
  const int maximumDim=dimensionality_;

  // Detect DMT-max cells to remove
  isRemovableMaximum.resize(numberOfCells);

  dmtMax2PL_.resize(numberOfCells);
  std::fill(dmtMax2PL_.begin(), dmtMax2PL_.end(), -1);

  // by default : maximum is removable
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(simplexId_t i=0; i<numberOfCells; ++i){
    const Cell maximumCandidate(maximumDim, i);
    isRemovableMaximum[i]=isMaximum(maximumCandidate);
  }

  for(int i=0; i<numberOfCriticalPoints; ++i){
    const std::pair<simplexId_t,char>& criticalPoint=criticalPoints[i];
    const simplexId_t criticalPointId=criticalPoint.first;
    const char criticalPointType=criticalPoint.second;

    if(criticalPointType==maximumDim){
      if(!allowBoundary and
         inputTriangulation_->isVertexOnBoundary(criticalPointId)) continue;

      int numberOfMaxima=0;
      simplexId_t maximumId=-1;
      const simplexId_t
        starNumber=inputTriangulation_->getVertexStarNumber(criticalPointId);
      for(simplexId_t j=0; j<starNumber; ++j){
        simplexId_t starId;
        inputTriangulation_->getVertexStar(criticalPointId, j, starId);

        if(isMaximum(Cell(maximumDim, starId)) and dmtMax2PL_[starId]==-1){
          maximumId=starId;
          ++numberOfMaxima;
        }
      }

      // a DMT-maximum in the star of only one PL-maximum cannot be removed
      // and is automatically associated to it.
      if(numberOfMaxima==1){
        if(dmtMax2PL_[maximumId]==-1 and pl2dmt_maximum[criticalPointId]==-1){
          dmtMax2PL_[maximumId]=criticalPointId;
          pl2dmt_maximum[criticalPointId]=maximumId;
          isRemovableMaximum[maximumId]=false;
        }
      }
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::getRemovableSaddles1(const std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                           const bool allowBoundary,
                                           std::vector<char>& isRemovableSaddle,
                                           std::vector<simplexId_t>& pl2dmt_saddle) const{
  const simplexId_t numberOfEdges=inputTriangulation_->getNumberOfEdges();
  isRemovableSaddle.resize(numberOfEdges);

  std::vector<char> dmt2PL(numberOfEdges, false);

  // by default : 1-saddle is removable
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(simplexId_t i=0; i<numberOfEdges; ++i){
    const Cell saddleCandidate(1, i);
    isRemovableSaddle[i]=isSaddle1(saddleCandidate);
  }

  // is [edgeId] in star of PL-1saddle?
  for(auto& criticalPoint : criticalPoints){
    const simplexId_t criticalPointId=criticalPoint.first;
    const char criticalPointType=criticalPoint.second;

    if(criticalPointType==1){
      if(!allowBoundary and
         inputTriangulation_->isVertexOnBoundary(criticalPointId)) continue;

      int numberOfSaddles=0;
      simplexId_t saddleId=-1;
      const simplexId_t
        edgeNumber=inputTriangulation_->getVertexEdgeNumber(criticalPointId);
      for(simplexId_t i=0; i<edgeNumber; ++i){
        simplexId_t edgeId;
        inputTriangulation_->getVertexEdge(criticalPointId, i, edgeId);
        const Cell saddleCandidate(1, edgeId);

        if(isSaddle1(saddleCandidate) and !dmt2PL[edgeId]){
          saddleId=edgeId;
          ++numberOfSaddles;
        }
      }

      // only one DMT-1saddle in the star so this one is non-removable
      if(numberOfSaddles==1){
        if(!dmt2PL[saddleId] and pl2dmt_saddle[criticalPointId]==-1){
          dmt2PL[saddleId]=true;
          pl2dmt_saddle[criticalPointId]=saddleId;
          isRemovableSaddle[saddleId]=false;
        }
      }
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::getRemovableSaddles2(const
                                           std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                           const bool allowBoundary,
                                           std::vector<char>& isRemovableSaddle,
                                           std::vector<simplexId_t>& pl2dmt_saddle) const{
  const simplexId_t numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  isRemovableSaddle.resize(numberOfTriangles);

  std::vector<char> dmt2PL(numberOfTriangles, false);

  // by default : 2-saddle is removable
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(simplexId_t i=0; i<numberOfTriangles; ++i){
    const Cell saddleCandidate(2, i);
    isRemovableSaddle[i]=isSaddle2(saddleCandidate);
  }

  // is [triangleId] in star of PL-2saddle?
  for(auto& criticalPoint : criticalPoints){
    const int criticalPointId=criticalPoint.first;
    const char criticalPointType=criticalPoint.second;

    if(criticalPointType==2){
      if(!allowBoundary and
         inputTriangulation_->isVertexOnBoundary(criticalPointId)) continue;

      simplexId_t numberOfSaddles=0;
      simplexId_t saddleId=-1;
      const simplexId_t
        triangleNumber=inputTriangulation_->getVertexTriangleNumber(criticalPointId);
      for(simplexId_t i=0; i<triangleNumber; ++i){
        simplexId_t triangleId;
        inputTriangulation_->getVertexTriangle(criticalPointId, i, triangleId);
        const Cell saddleCandidate(2, triangleId);

        if(isSaddle2(saddleCandidate) and !dmt2PL[triangleId]){
          saddleId=triangleId;
          ++numberOfSaddles;
        }
      }

      // only one DMT-2saddle in the star so this one is non-removable
      if(numberOfSaddles==1){
        if(dmt2PL[saddleId]==false and pl2dmt_saddle[criticalPointId]==-1){
          dmt2PL[saddleId]=true;
          pl2dmt_saddle[criticalPointId]=saddleId;
          isRemovableSaddle[saddleId]=false;
        }
      }
    }
  }

  return  0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleMaximumConnections(const
                                                    std::vector<VPath>& vpaths,
                                                    std::set<std::pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>&
                                                    S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_)
      S.insert(std::make_pair(vpath.persistence_,i));
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::computeCoefficients(const bool isDense,
                                          std::vector<char>& denseCoefficients,
                                          std::vector<Segment>& segments,
                                          const CriticalPoint& source,
                                          VPath& newVPath,
                                          const std::vector<VPath>& vpaths) const{
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
    std::vector<std::pair<int,char>> sparseCoefficients;

    // 1) initialize accumulator
    const int numberOfNewVPathSegments=newVPath.segments_.size();
    for(int i=0; i<numberOfNewVPathSegments; ++i){
      const int segmentId=newVPath.segments_[i];
      const char segmentState=newVPath.states_[i];

      sparseCoefficients.push_back(std::make_pair(segmentId,segmentState));
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

            sparseCoefficients.push_back(std::make_pair(segmentId,segmentState));
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

      // apriori : sparseCoefficients store coefficient zero; we must remove
      // them
      if(segmentState!=0){
        newVPath.states_.push_back(segmentState);
        newVPath.segments_.push_back(segmentId);
      }
    }
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::reverseSaddleMaximumConnections(const
                                                      std::vector<Segment>& segments){
  Timer t;

  const int numberOfSegments=segments.size();

  for(int i=0; i<numberOfSegments; ++i){
    const Segment& segment=segments[i];
    if(segment.isValid_ and segment.orientation_==false)
      reverseAscendingPath(segment.cells_);
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Gradient reversal step :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int
DiscreteGradient::initializeSaddleMaximumConnections(std::vector<char>&
isRemovableMaximum,
                                                     std::vector<char>& isRemovableSaddle,
                                                     const bool allowBruteForce,
                                                     std::vector<Segment>& segments,
                                                     std::vector<VPath>& vpaths,
                                                     std::vector<CriticalPoint>& criticalPoints) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int maximumDim=dimensionality_;
  const int saddleDim=maximumDim-1;

  // Part 1 : build initial structures
  // add the saddles to CriticalPointList and count them
  const int numberOfSaddleCandidates=getNumberOfCells(saddleDim);
  for(int i=0; i<numberOfSaddleCandidates; ++i){
    if(allowBruteForce or isRemovableSaddle[i]){
      const Cell saddleCandidate(saddleDim, i);

      if(isCellCritical(saddleCandidate))
        criticalPoints.push_back(CriticalPoint(saddleCandidate));
    }
  }
  const int numberOfSaddles=criticalPoints.size();

  // add the maxima to CriticalPointList and build MaxIndex
  const int numberOfMaximumCandidates=getNumberOfCells(maximumDim);
  std::vector<int> maximumIndex(numberOfMaximumCandidates,-1);
  for(int i=0; i<numberOfMaximumCandidates; ++i){
    if(isRemovableMaximum[i]){
      const Cell maximumCandidate(maximumDim, i);

      const int index=criticalPoints.size();
      maximumIndex[i]=index;

      criticalPoints.push_back(CriticalPoint(maximumCandidate));
    }
  }

  const int numberOfVPaths=2*numberOfSaddles;
  vpaths.resize(numberOfVPaths);
  segments.resize(numberOfVPaths);

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int sourceIndex=i;
    CriticalPoint& source=criticalPoints[sourceIndex];

    const Cell& saddle=source.cell_;
    const simplexId_t saddleId=saddle.id_;

    simplexId_t starNumber{};
    if(maximumDim==2)
      starNumber=inputTriangulation_->getEdgeStarNumber(saddleId);
    else if(maximumDim==3)
      starNumber=inputTriangulation_->getTriangleStarNumber(saddleId);

    std::vector<std::vector<Cell>> paths(starNumber);
    for(simplexId_t j=0; j<starNumber; ++j){
      simplexId_t starId;
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
      for(simplexId_t j=1; j<starNumber; ++j){
        const Cell& lastCell=paths[j].back();

        if(lastCell0.id_==lastCell.id_){
          isDoubleConnected=true;
          break;
        }
      }
      if(isDoubleConnected)
        continue;
    }

    for(simplexId_t j=0; j<starNumber; ++j){
      const simplexId_t shift=j;

      // apriori: there is at least 1 one cell
      const Cell& lastCell=paths[j].back();
      if(isMaximum(lastCell) and isRemovableMaximum[lastCell.id_]){
        const Cell& maximum=lastCell;

        const int destinationIndex=maximumIndex[maximum.id_];
        CriticalPoint& destination=criticalPoints[destinationIndex];

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

  // Part 3 : initialize the last structures
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    CriticalPoint& cp=criticalPoints[i];

    const int numberOfSlots=cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_=0;
  }

#ifdef TTK_ENABLE_OPENMP
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
    std::stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime()
        << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleMaximumConnections(const int
                                                      iterationThreshold,
                                                      const std::vector<char>& isPL,
                                                      const bool allowBoundary,
                                                      const bool allowBruteForce,
                                                      std::set<std::pair<dataType,int>,SaddleMaximumVPathComparator<dataType>>& S,
                                                      std::vector<simplexId_t>& pl2dmt_saddle,
                                                      std::vector<simplexId_t>& pl2dmt_maximum,
                                                      std::vector<Segment>& segments,
                                                      std::vector<VPath>& vpaths,
                                                      std::vector<CriticalPoint>& criticalPoints) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  simplexId_t numberOfSaddleCandidates=0;
  if(dimensionality_==2)
    numberOfSaddleCandidates=inputTriangulation_->getNumberOfEdges();
  else if(dimensionality_==3)
    numberOfSaddleCandidates=inputTriangulation_->getNumberOfTriangles();
  const simplexId_t numberOfMaximumCandidates=inputTriangulation_->getNumberOfCells();

  const int maximumDim=dimensionality_;
  const int saddleDim=maximumDim-1;

  std::vector<char> isRemovedSaddle(numberOfSaddleCandidates, false);
  std::vector<char> isRemovedMaximum(numberOfMaximumCandidates, false);

  int numberOfIterations{};
  std::vector<char> denseCoefficients;
  while(S.size()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=ptr->second;
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    // filter by saddle condition
    int toRemoveSaddle=0;
    if(!allowBruteForce and vpath.isValid_){
      const int sourceId=vpath.source_;
      const simplexId_t dmt_saddleId=criticalPoints[sourceId].cell_.id_;

      if(!isRemovedSaddle[dmt_saddleId]){
        for(int i=0; i<(saddleDim+1); ++i){
          simplexId_t vertexId=-1;
          if(dimensionality_==2)
            inputTriangulation_->getEdgeVertex(dmt_saddleId, i, vertexId);
          else if(dimensionality_==3)
            inputTriangulation_->getTriangleVertex(dmt_saddleId, i, vertexId);

          if(isPL[vertexId]!=saddleDim) continue;

          if(!allowBoundary and
             inputTriangulation_->isVertexOnBoundary(vertexId)){
            toRemoveSaddle=1;
            continue;
          }

          if(pl2dmt_saddle[vertexId]==-1){
            const simplexId_t pl_saddleId=vertexId;

            int numberOfRemainingSaddles=0;

            simplexId_t saddleCandidateNumber=0;
            if(dimensionality_==2)

              saddleCandidateNumber=inputTriangulation_->getVertexEdgeNumber(pl_saddleId);
            else if(dimensionality_==3)

              saddleCandidateNumber=inputTriangulation_->getVertexTriangleNumber(pl_saddleId);

            for(simplexId_t j=0; j<saddleCandidateNumber; ++j){
              simplexId_t saddleCandidateId=-1;
              if(dimensionality_==2)
                inputTriangulation_->getVertexEdge(pl_saddleId, j,
                                                   saddleCandidateId);
              else if(dimensionality_==3)
                inputTriangulation_->getVertexTriangle(pl_saddleId, j,
                                                       saddleCandidateId);

              if(saddleCandidateId!=dmt_saddleId and
                 isCellCritical(Cell(saddleDim,saddleCandidateId)) and
                 !isRemovedSaddle[saddleCandidateId])
                ++numberOfRemainingSaddles;
            }

            if(!numberOfRemainingSaddles){
              pl2dmt_saddle[vertexId]=dmt_saddleId;
              vpath.invalidate();
              toRemoveSaddle=-1;
              break;
            }
          }
          else if(pl2dmt_saddle[vertexId]==dmt_saddleId){
            vpath.invalidate();
            toRemoveSaddle=-1;
            break;
          }
          else{
            toRemoveSaddle=1;
            break;
          }
        }

        if(vpath.isValid_){
          toRemoveSaddle=1;
        }
      }
      else{
        vpath.invalidate();
        toRemoveSaddle=-1;
      }
    }

    // filter by maximum condition
    int toRemoveMaximum=0;
    if(vpath.isValid_){
      const int destinationId=vpath.destination_;
      const simplexId_t dmt_maxId=criticalPoints[destinationId].cell_.id_;

      if(!isRemovedMaximum[dmt_maxId]){
        for(int i=0; i<(maximumDim+1); ++i){
          simplexId_t vertexId;
          inputTriangulation_->getCellVertex(dmt_maxId, i, vertexId);

          if(isPL[vertexId]!=maximumDim) continue;

          if(!allowBoundary and
             inputTriangulation_->isVertexOnBoundary(vertexId)){
            toRemoveMaximum=1;
            continue;
          }

          if(pl2dmt_maximum[vertexId]==-1){
            const simplexId_t pl_maxId=vertexId;

            simplexId_t numberOfRemainingMaxima=0;
            const simplexId_t
              starNumber=inputTriangulation_->getVertexStarNumber(pl_maxId);
            for(simplexId_t j=0; j<starNumber; ++j){
              simplexId_t starId;
              inputTriangulation_->getVertexStar(pl_maxId, j, starId);
              if(starId!=dmt_maxId and isMaximum(Cell(maximumDim,starId)) and
                 !isRemovedMaximum[starId])
                ++numberOfRemainingMaxima;
            }

            if(!numberOfRemainingMaxima){
              pl2dmt_maximum[vertexId]=dmt_maxId;
              vpath.invalidate();
              toRemoveMaximum=-1;
              break;
            }
          }
          else if(pl2dmt_maximum[vertexId]==dmt_maxId){
            vpath.invalidate();
            toRemoveMaximum=-1;
            break;
          }
          else{
            toRemoveMaximum=1;
            break;
          }
        }

        if(vpath.isValid_){
          toRemoveMaximum=1;
        }
      }
      else{
        vpath.invalidate();
        toRemoveMaximum=-1;
      }
    }

    // sync removed-state
    if(vpath.isValid_){
      if((toRemoveSaddle>=0 and toRemoveMaximum>0) or (toRemoveSaddle>0 and
                                                       toRemoveMaximum>=0)){
        const int sourceId=vpath.source_;
        const simplexId_t dmt_saddleId=criticalPoints[sourceId].cell_.id_;

        const int destinationId=vpath.destination_;
        const simplexId_t dmt_maxId=criticalPoints[destinationId].cell_.id_;

        isRemovedSaddle[dmt_saddleId]=true;
        isRemovedMaximum[dmt_maxId]=true;
      }
    }

    if(vpath.isValid_){
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

          if(newSourceVPath.isValid_ and
             newSourceVPath.destination_==newDestinationId){
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
        computeCoefficients<dataType>(false, denseCoefficients, segments,
                                      source, newVPath, vpaths);

        // update the destination of newVPath
        newVPath.destination_=newDestinationId;

        // add newVPath to newDestination connectivity
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        newDestination.vpaths_.push_back(newVPathId);

        // erase newVPath
        S.erase(std::make_pair(newVPath.persistence_,newVPathId));

        // update persistence

        newVPath.persistence_=getPersistence<dataType>(newDestination.cell_,newSource.
          cell_, scalars);

        // repush newVPath to confirm update
        S.insert(std::make_pair(newVPath.persistence_,newVPathId));
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
    std::stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleMaximumConnections(const
                                                       std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                                       const std::vector<char>& isPL,
                                                       const int iterationThreshold,
                                                       const bool allowBoundary,
                                                       const bool allowBruteForce){
  Timer t;

  // Part 0 : select the cells to keep or remove (gradient
  // is not modified).
  std::vector<char> isRemovableMaximum;
  std::vector<simplexId_t> pl2dmt_maximum(numberOfVertices_, -1);
  getRemovableMaxima<dataType>(criticalPoints, allowBoundary,
                               isRemovableMaximum, pl2dmt_maximum);

  std::vector<char> isRemovableSaddle;
  std::vector<simplexId_t> pl2dmt_saddle(numberOfVertices_, -1);
  if(!allowBruteForce){
    if(dimensionality_==2)
      getRemovableSaddles1<dataType>(criticalPoints, allowBoundary,
                                     isRemovableSaddle, pl2dmt_saddle);
    else if(dimensionality_==3)
      getRemovableSaddles2<dataType>(criticalPoints, allowBoundary,
                                     isRemovableSaddle, pl2dmt_saddle);
  }

  // Part 1 : build a virtual but complete MSC structure as
  // initialization (gradient is not modified).
  std::vector<Segment> segments;
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  initializeSaddleMaximumConnections<dataType>(isRemovableMaximum,
                                               isRemovableSaddle,
                                               allowBruteForce,
                                               segments,
                                               vpaths,
                                               dmt_criticalPoints);

  // Part 2 : push the vpaths into a set to order them by persistence
  // value - lower to higher (gradient is not modified).
  SaddleMaximumVPathComparator<dataType> cmp_f;
  std::set<std::pair<dataType,int>, SaddleMaximumVPathComparator<dataType>>
    S(cmp_f);
  orderSaddleMaximumConnections<dataType>(vpaths, S);

  // Part 3 : iteratively process the vpaths, virtually reverse the
  // selected vpath and update the structure accordingly (gradient is
  // not modified).
  processSaddleMaximumConnections<dataType>(iterationThreshold,
                                            isPL,
                                            allowBoundary,
                                            allowBruteForce,
                                            S,
                                            pl2dmt_saddle,
                                            pl2dmt_maximum,
                                            segments,
                                            vpaths,
                                            dmt_criticalPoints);

  // Part 4 : from the last version of the virtual MSC: use the
  // informations stored in the virtual structure to actually
  // reverse the vpaths in the gradient (gradient is modified).
  reverseSaddleMaximumConnections<dataType>(segments);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle-Maximum pairs simplified in "
        << t.getElapsedTime() << " s, "<< threadNumber_ << " thread(s)." <<
        std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections1(const
                                                         std::vector<char>& isRemovableSaddle1,
                                                         const std::vector<char>& isRemovableSaddle2,
                                                         const bool allowBruteForce,
                                                         std::vector<VPath>& vpaths,
                                                         std::vector<CriticalPoint>& criticalPoints,
                                                         std::vector<simplexId_t>& saddle1Index,
                                                         std::vector<simplexId_t>& saddle2Index) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int maximumDim=dimensionality_;
  const int saddle2Dim=maximumDim-1;
  const int saddle1Dim=saddle2Dim-1;

  // Part 1 : build initial structures
  // add the 2-saddles to CriticalPointList
  const simplexId_t numberOfSaddle2Candidates=getNumberOfCells(saddle2Dim);
  saddle2Index.resize(numberOfSaddle2Candidates, -1);
  for(simplexId_t i=0; i<numberOfSaddle2Candidates; ++i){
    if(allowBruteForce or isRemovableSaddle2[i]){
      const Cell saddle2Candidate(saddle2Dim, i);

      if(isSaddle2(saddle2Candidate)){
        const simplexId_t index=criticalPoints.size();
        saddle2Index[i]=index;
        criticalPoints.push_back(CriticalPoint(saddle2Candidate));
      }
    }
  }
  const simplexId_t numberOf2Saddles=criticalPoints.size();

  // add the 1-saddles to CriticalPointList
  const simplexId_t numberOfSaddle1Candidates=getNumberOfCells(saddle1Dim);
  saddle1Index.resize(numberOfSaddle1Candidates, -1);
  for(simplexId_t i=0; i<numberOfSaddle1Candidates; ++i){
    if(isRemovableSaddle1[i]){
      const Cell saddle1Candidate(saddle1Dim, i);

      const simplexId_t index=criticalPoints.size();
      saddle1Index[i]=index;
      criticalPoints.push_back(CriticalPoint(saddle1Candidate));
    }
  }

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
  wallId_t descendingWallId=1;
  std::vector<wallId_t> isVisited(numberOfSaddle2Candidates, 0);
  for(simplexId_t i=0; i<numberOf2Saddles; ++i){
    const simplexId_t destinationIndex=i;
    CriticalPoint& destination=criticalPoints[destinationIndex];
    const Cell& saddle2=destination.cell_;

    std::set<simplexId_t> saddles1;
    const wallId_t savedDescendingWallId=descendingWallId;
    getDescendingWall(descendingWallId, saddle2, isVisited, nullptr, &saddles1);
    ++descendingWallId;

    for(auto& saddle1Id : saddles1){
      if(!isRemovableSaddle1[saddle1Id]) continue;

      const Cell& saddle1=Cell(1,saddle1Id);

      std::vector<Cell> path;
      const bool
        isMultiConnected=getAscendingPathThroughWall(savedDescendingWallId, saddle1,
                                                     saddle2, isVisited, &path);

      if(!isMultiConnected){
        const simplexId_t sourceIndex=saddle1Index[saddle1Id];
        CriticalPoint& source=criticalPoints[sourceIndex];

        // update source and destination
        const int sourceSlot=source.addSlot();
        const int destinationSlot=destination.addSlot();

        // update vpath
        const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                            scalars);

        vpaths.push_back(VPath(true,-1,sourceIndex,destinationIndex,sourceSlot,
                               destinationSlot,persistence));
      }
    }
  }

  // Part 3 : initialize the last structures
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    CriticalPoint& cp=criticalPoints[i];

    const int numberOfSlots=cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_=0;
  }

  const int numberOfVPaths=vpaths.size();
#ifdef TTK_ENABLE_OPENMP
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
    std::stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" << t.getElapsedTime()
        << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections1(const
                                                    std::vector<VPath>& vpaths,
                                                    std::vector<CriticalPoint>& criticalPoints,
                                                    std::set<std::tuple<dataType,int,simplexId_t>,
                                                      SaddleSaddleVPathComparator<dataType>> &S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const simplexId_t saddleId=criticalPoints[vpath.destination_].cell_.id_;
      S.insert(std::make_tuple(vpath.persistence_,i,saddleId));
    }
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleSaddleConnections1(const int
                                                      iterationThreshold,
                                                      const std::vector<char>& isPL,
                                                      const bool allowBoundary,
                                                      const bool allowBruteForce,
                                                      const bool returnSaddleConnectors,
                                                      std::set<std::tuple<dataType,int,simplexId_t>,
                                                        SaddleSaddleVPathComparator<dataType>>& S,
                                                      std::vector<simplexId_t>& pl2dmt_saddle1,
                                                      std::vector<simplexId_t>& pl2dmt_saddle2,
                                                      std::vector<char>& isRemovableSaddle1,
                                                      std::vector<char>& isRemovableSaddle2,
                                                      std::vector<VPath>& vpaths,
                                                      std::vector<CriticalPoint>& criticalPoints,
                                                      std::vector<simplexId_t>& saddle1Index,
                                                      std::vector<simplexId_t>& saddle2Index){
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const simplexId_t numberOfEdges=inputTriangulation_->getNumberOfEdges();
  const simplexId_t numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  const simplexId_t optimizedSize=std::max(numberOfEdges, numberOfTriangles);
  wallId_t wallId=1;
  std::vector<wallId_t> isVisited(optimizedSize, 0);

  int numberOfIterations{};
  while(!S.empty()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=std::get<1>(*ptr);
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    if(vpath.isValid_){
      if(returnSaddleConnectors){
        const dataType persistence=vpath.persistence_;
        if(persistence>SaddleConnectorsPersistenceThreshold) break;
      }

      const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
      const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;

      std::set<simplexId_t> saddles1;
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
      std::vector<Cell> path;
      const bool isMultiConnected=getAscendingPathThroughWall(savedWallId,
                                                              minSaddle1, minSaddle2, isVisited, &path);
      if(isMultiConnected){
        ++numberOfIterations;
        continue;
      }

      // filter by 1-saddle condition
      if(vpath.isValid_){
        const Cell& dmt_saddle1=criticalPoints[vpath.source_].cell_;
        const simplexId_t dmt_saddle1Id=dmt_saddle1.id_;

        if(isSaddle1(dmt_saddle1)){
          for(int i=0; i<2; ++i){
            simplexId_t vertexId;
            inputTriangulation_->getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId]!=1) continue;

            if(!allowBoundary and
               inputTriangulation_->isVertexOnBoundary(vertexId)) continue;

            if(pl2dmt_saddle1[vertexId]==-1){
              const simplexId_t pl_saddle1Id=vertexId;

              simplexId_t numberOfRemainingSaddles1=0;

              simplexId_t savedId=-1;
              const simplexId_t
                edgeNumber=inputTriangulation_->getVertexEdgeNumber(pl_saddle1Id);
              for(simplexId_t j=0; j<edgeNumber; ++j){
                simplexId_t edgeId;
                inputTriangulation_->getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId!=dmt_saddle1Id and isSaddle1(Cell(1,edgeId)) and
                   isRemovableSaddle1[edgeId]){
                  ++numberOfRemainingSaddles1;
                  savedId=edgeId;
                }
              }

              if(numberOfRemainingSaddles1==0){
                isRemovableSaddle1[dmt_saddle1Id]=false;
                pl2dmt_saddle1[vertexId]=dmt_saddle1Id;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1==1){
                isRemovableSaddle1[dmt_saddle1Id]=false;
                isRemovableSaddle1[savedId]=false;
                pl2dmt_saddle1[vertexId]=savedId;
                break;
              }
            }
            else if(pl2dmt_saddle1[vertexId]==dmt_saddle1Id){
              vpath.invalidate();
              break;
            }
          }
        }
        else
          vpath.invalidate();
      }

      // filter by 2-saddle condition
      if(!allowBruteForce and vpath.isValid_){
        const Cell& dmt_saddle2=criticalPoints[vpath.destination_].cell_;
        const simplexId_t dmt_saddle2Id=dmt_saddle2.id_;

        if(isSaddle2(dmt_saddle2)){
          for(int i=0; i<3; ++i){
            simplexId_t vertexId;
            inputTriangulation_->getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId]!=2) continue;

            if(!allowBoundary and
               inputTriangulation_->isVertexOnBoundary(vertexId)) continue;

            if(pl2dmt_saddle2[vertexId]==-1){
              const simplexId_t pl_saddle2Id=vertexId;

              simplexId_t numberOfRemainingSaddles2=0;

              simplexId_t savedId=-1;
              const simplexId_t
                triangleNumber=inputTriangulation_->getVertexTriangleNumber(pl_saddle2Id);
              for(simplexId_t j=0; j<triangleNumber; ++j){
                simplexId_t triangleId;
                inputTriangulation_->getVertexTriangle(pl_saddle2Id, j,
                                                       triangleId);

                if(triangleId!=dmt_saddle2Id and isSaddle2(Cell(2,triangleId))
                   and isRemovableSaddle2[triangleId]){
                  ++numberOfRemainingSaddles2;
                  savedId=triangleId;
                }
              }

              if(numberOfRemainingSaddles2==0){
                isRemovableSaddle2[dmt_saddle2Id]=false;
                pl2dmt_saddle2[vertexId]=dmt_saddle2Id;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles2==1){
                isRemovableSaddle2[dmt_saddle2Id]=false;
                isRemovableSaddle2[savedId]=false;
                pl2dmt_saddle2[vertexId]=savedId;
                break;
              }
            }
            else if(pl2dmt_saddle2[vertexId]==dmt_saddle2Id){
              vpath.invalidate();
              break;
            }
          }
        }
        else
          vpath.invalidate();
      }

      if(vpath.isValid_)
        reverseAscendingPathOnWall(path);
    }

    if(vpath.isValid_){
      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_){
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back(std::make_tuple(minSaddle1, minSaddle2));
      }

      const int sourceId=vpath.source_;
      const int destinationId=vpath.destination_;

      // invalidate vpaths connected to destination
      std::vector<int> newSourceIds;
      CriticalPoint& destination=criticalPoints[destinationId];
      for(auto& destinationVPathId : destination.vpaths_){
        VPath& destinationVPath=vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_!=sourceId){
          // save critical point
          const int newSourceId=destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          // clear vpath
          destinationVPath.invalidate();
        }
      }

      // invalidate vpaths connected to source and save the critical points to
      // update
      std::vector<int> newDestinationIds;
      CriticalPoint& source=criticalPoints[sourceId];
      for(auto& sourceVPathId : source.vpaths_){
        VPath& sourceVPath=vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_!=destinationId){
          // save critical point
          const int newDestinationId=sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          CriticalPoint& newDestination=criticalPoints[newDestinationId];
          for(auto& newDestinationVPathId : newDestination.vpaths_){
            VPath& newDestinationVPath=vpaths[newDestinationVPathId];
            if(newDestinationVPath.isValid_ and
               newDestinationVPath.source_!=sourceId){

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
      for(auto& newDestinationId : newDestinationIds){
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        const Cell& saddle2=newDestination.cell_;

        std::set<simplexId_t> saddles1;
        const wallId_t savedWallId=wallId;
        getDescendingWall(wallId, saddle2, isVisited, nullptr, &saddles1);
        ++wallId;

        for(auto& saddle1Id : saddles1){
          const Cell saddle1(1,saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected=getAscendingPathThroughWall(savedWallId,
                                                                  saddle1, saddle2, isVisited, &path);
          if(isMultiConnected)
            continue;

          simplexId_t newSourceId=saddle1Index[saddle1Id];

          // connection to a new saddle1 (not present in the graph before)
          if(newSourceId==-1){
            if(!isRemovableSaddle1[saddle1Id]) continue;

            const simplexId_t newCriticalPointId=criticalPoints.size();
            saddle1Index[saddle1Id]=newCriticalPointId;
            criticalPoints.push_back(CriticalPoint(saddle1));

            newSourceId=newCriticalPointId;
          }
          CriticalPoint& newSource=criticalPoints[newSourceId];

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                              scalars);

          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(std::make_tuple(persistence,newVPathId,newDestination.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(auto& newSourceId : newSourceIds){
        CriticalPoint& newSource=criticalPoints[newSourceId];
        const Cell& saddle1=newSource.cell_;

        std::set<simplexId_t> saddles2;
        const wallId_t savedWallId=wallId;
        getAscendingWall(wallId, saddle1, isVisited, nullptr, &saddles2);
        ++wallId;

        for(auto& saddle2Id : saddles2){
          const Cell saddle2(2,saddle2Id);

          std::vector<Cell> path;
          const bool isMultiConnected=getDescendingPathThroughWall(savedWallId,
                                                                   saddle2, saddle1, isVisited, &path);
          if(isMultiConnected)
            continue;

          const simplexId_t newDestinationId=saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId==-1)
            continue;

          CriticalPoint& newDestination=criticalPoints[newDestinationId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists=false;
          for(auto& newDestinationVPathId : newDestination.vpaths_){
            const VPath& newDestinationVPath=vpaths[newDestinationVPathId];

            if(newDestinationVPath.isValid_ and
               newDestinationVPath.source_==newSourceId){
              alreadyExists=true;
              break;
            }
          }

          if(alreadyExists)
            continue;

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                              scalars);

          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(std::make_tuple(persistence,newVPathId,newDestination.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections1(const
                                                       std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                                       const std::vector<char>& isPL,
                                                       const int iterationThreshold,
                                                       const bool allowBoundary,
                                                       const bool allowBruteForce,
                                                       const bool returnSaddleConnectors){
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<simplexId_t> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle1, pl2dmt_saddle1);

  std::vector<char> isRemovableSaddle2;
  std::vector<simplexId_t> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle2, pl2dmt_saddle2);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<simplexId_t> saddle1Index;
  std::vector<simplexId_t> saddle2Index;
  initializeSaddleSaddleConnections1<dataType>(isRemovableSaddle1,
                                               isRemovableSaddle2,
                                               allowBruteForce,
                                               vpaths,
                                               dmt_criticalPoints,
                                               saddle1Index,
                                               saddle2Index);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  std::set<std::tuple<dataType,int,simplexId_t>, SaddleSaddleVPathComparator<dataType>>
    S(cmp_f);
  orderSaddleSaddleConnections1<dataType>(vpaths, dmt_criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections1<dataType>(iterationThreshold,
                                            isPL,
                                            allowBoundary,
                                            allowBruteForce,
                                            returnSaddleConnectors,
                                            S,
                                            pl2dmt_saddle1,
                                            pl2dmt_saddle2,
                                            isRemovableSaddle1,
                                            isRemovableSaddle2,
                                            vpaths,
                                            dmt_criticalPoints,
                                            saddle1Index,
                                            saddle2Index);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle-Saddle pairs simplified in "
        << t.getElapsedTime() << " s, "<< threadNumber_ << " thread(s)." <<
        std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::initializeSaddleSaddleConnections2(const
                                                         std::vector<char>& isRemovableSaddle1,
                                                         const std::vector<char>& isRemovableSaddle2,
                                                         const bool allowBruteForce,
                                                         std::vector<VPath>& vpaths,
                                                         std::vector<CriticalPoint>& criticalPoints,
                                                         std::vector<simplexId_t>& saddle1Index,
                                                         std::vector<simplexId_t>& saddle2Index) const{
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const int maximumDim=dimensionality_;
  const int saddle2Dim=maximumDim-1;
  const int saddle1Dim=saddle2Dim-1;

  // Part 1 : build initial structures
  // add the 1-saddles to CriticalPointList
  const simplexId_t numberOfSaddle1Candidates=getNumberOfCells(saddle1Dim);
  saddle1Index.resize(numberOfSaddle1Candidates, -1);
  for(simplexId_t i=0; i<numberOfSaddle1Candidates; ++i){
    if(isRemovableSaddle1[i]){
      const Cell saddle1Candidate(saddle1Dim, i);

      const simplexId_t index=criticalPoints.size();
      saddle1Index[i]=index;
      criticalPoints.push_back(CriticalPoint(saddle1Candidate));
    }
  }
  const simplexId_t numberOf1Saddles=criticalPoints.size();

  // add the 2-saddles to CriticalPointList
  const simplexId_t numberOfSaddle2Candidates=getNumberOfCells(saddle2Dim);
  saddle2Index.resize(numberOfSaddle2Candidates, -1);
  for(simplexId_t i=0; i<numberOfSaddle2Candidates; ++i){
    if(allowBruteForce or isRemovableSaddle2[i]){
      const Cell saddle2Candidate(saddle2Dim, i);

      if(isSaddle2(saddle2Candidate)){
        const simplexId_t index=criticalPoints.size();
        saddle2Index[i]=index;
        criticalPoints.push_back(CriticalPoint(saddle2Candidate));
      }
    }
  }

  // Part 2 : update the structures
  // apriori: by default construction, the vpaths and segments are not valid
  wallId_t ascendingWallId=1;
  std::vector<wallId_t> isVisited(numberOfSaddle1Candidates, 0);
  for(simplexId_t i=0; i<numberOf1Saddles; ++i){
    const simplexId_t sourceIndex=i;
    CriticalPoint& source=criticalPoints[sourceIndex];
    const Cell& saddle1=source.cell_;

    std::set<simplexId_t> saddles2;
    const wallId_t savedAscendingWallId=ascendingWallId;
    getAscendingWall(ascendingWallId, saddle1, isVisited, nullptr, &saddles2);
    ++ascendingWallId;

    for(auto& saddle2Id : saddles2){
      if(!isRemovableSaddle2[saddle2Id]) continue;

      const Cell& saddle2=Cell(2,saddle2Id);

      std::vector<Cell> path;
      const bool
        isMultiConnected=getDescendingPathThroughWall(savedAscendingWallId, saddle2,
                                                      saddle1, isVisited, &path);

      if(!isMultiConnected){
        const simplexId_t destinationIndex=saddle2Index[saddle2Id];
        CriticalPoint& destination=criticalPoints[destinationIndex];

        // update source and destination
        const int sourceSlot=source.addSlot();
        const int destinationSlot=destination.addSlot();

        // update vpath
        const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                            scalars);

        vpaths.push_back(VPath(true,-1,sourceIndex,destinationIndex,sourceSlot,
                               destinationSlot,persistence));
      }
    }
  }

  // Part 3 : initialize the last structures
  const simplexId_t numberOfCriticalPoints=criticalPoints.size();
  for(simplexId_t i=0; i<numberOfCriticalPoints; ++i){
    CriticalPoint& cp=criticalPoints[i];

    const int numberOfSlots=cp.numberOfSlots_;
    cp.vpaths_.resize(numberOfSlots);
    cp.numberOfSlots_=0;
  }

  const int numberOfVPaths=vpaths.size();
#ifdef TTK_ENABLE_OPENMP
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
    std::stringstream msg;
    msg << "[DiscreteGradient]  Initialization step :\t" <<
        t.getElapsedTime()
        << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }
  return 0;
}

template <typename dataType>
int DiscreteGradient::orderSaddleSaddleConnections2(const
                                                    std::vector<VPath>& vpaths,
                                                    std::vector<CriticalPoint>& criticalPoints,
                                                    std::set<std::tuple<dataType,int,simplexId_t>,
                                                      SaddleSaddleVPathComparator<dataType>> &S){
  Timer t;

  const int numberOfVPaths=vpaths.size();
  for(int i=0; i<numberOfVPaths; ++i){
    const VPath& vpath=vpaths[i];

    if(vpath.isValid_){
      const simplexId_t saddleId=criticalPoints[vpath.source_].cell_.id_;
      S.insert(std::make_tuple(vpath.persistence_,i,saddleId));
    }
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Ordering of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::processSaddleSaddleConnections2(const int
                                                      iterationThreshold,
                                                      const std::vector<char>& isPL,
                                                      const bool allowBoundary,
                                                      const bool allowBruteForce,
                                                      const bool returnSaddleConnectors,
                                                      std::set<std::tuple<dataType,int,simplexId_t>,
                                                        SaddleSaddleVPathComparator<dataType>>& S,
                                                      std::vector<simplexId_t>& pl2dmt_saddle1,
                                                      std::vector<simplexId_t>& pl2dmt_saddle2,
                                                      std::vector<char>& isRemovableSaddle1,
                                                      std::vector<char>& isRemovableSaddle2,
                                                      std::vector<VPath>& vpaths,
                                                      std::vector<CriticalPoint>& criticalPoints,
                                                      std::vector<simplexId_t>& saddle1Index,
                                                      std::vector<simplexId_t>& saddle2Index){
  Timer t;

  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  const simplexId_t numberOfEdges=inputTriangulation_->getNumberOfEdges();
  const simplexId_t numberOfTriangles=inputTriangulation_->getNumberOfTriangles();
  const simplexId_t optimizedSize=std::max(numberOfEdges, numberOfTriangles);
  wallId_t wallId=1;
  std::vector<wallId_t> isVisited(optimizedSize, 0);

  int numberOfIterations{};
  while(!S.empty()){
    if(iterationThreshold>=0 and numberOfIterations>=iterationThreshold) break;

    auto ptr=S.begin();
    const int vpathId=std::get<1>(*ptr);
    S.erase(ptr);
    VPath& vpath=vpaths[vpathId];

    if(vpath.isValid_){
      if(returnSaddleConnectors){
        const dataType persistence=vpath.persistence_;
        if(persistence>SaddleConnectorsPersistenceThreshold) break;
      }

      const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
      const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;

      std::set<simplexId_t> saddles2;
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
      std::vector<Cell> path;
      const bool isMultiConnected=getDescendingPathThroughWall(savedWallId,
                                                               minSaddle2, minSaddle1, isVisited, &path);
      if(isMultiConnected){
        ++numberOfIterations;
        continue;
      }

      // filter by 1-saddle condition
      if(vpath.isValid_){
        const Cell& dmt_saddle1=criticalPoints[vpath.source_].cell_;
        const simplexId_t dmt_saddle1Id=dmt_saddle1.id_;

        if(isSaddle1(dmt_saddle1)){
          for(int i=0; i<2; ++i){
            simplexId_t vertexId;
            inputTriangulation_->getEdgeVertex(dmt_saddle1Id, i, vertexId);

            if(isPL[vertexId]!=1) continue;

            if(!allowBoundary and
               inputTriangulation_->isVertexOnBoundary(vertexId)) continue;

            if(pl2dmt_saddle1[vertexId]==-1){
              const simplexId_t pl_saddle1Id=vertexId;

              simplexId_t numberOfRemainingSaddles1=0;

              simplexId_t savedId=-1;
              const simplexId_t
                edgeNumber=inputTriangulation_->getVertexEdgeNumber(pl_saddle1Id);
              for(simplexId_t j=0; j<edgeNumber; ++j){
                simplexId_t edgeId;
                inputTriangulation_->getVertexEdge(pl_saddle1Id, j, edgeId);

                if(edgeId!=dmt_saddle1Id and isSaddle1(Cell(1,edgeId)) and
                   isRemovableSaddle1[edgeId]){
                  ++numberOfRemainingSaddles1;
                  savedId=edgeId;
                }
              }

              if(numberOfRemainingSaddles1==0){
                isRemovableSaddle1[dmt_saddle1Id]=false;
                pl2dmt_saddle1[vertexId]=dmt_saddle1Id;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles1==1){
                isRemovableSaddle1[dmt_saddle1Id]=false;
                isRemovableSaddle1[savedId]=false;
                pl2dmt_saddle1[vertexId]=savedId;
                break;
              }
            }
            else if(pl2dmt_saddle1[vertexId]==dmt_saddle1Id){
              vpath.invalidate();
              break;
            }
          }
        }
        else
          vpath.invalidate();
      }

      // filter by 2-saddle condition
      if(!allowBruteForce and vpath.isValid_){
        const Cell& dmt_saddle2=criticalPoints[vpath.destination_].cell_;
        const simplexId_t dmt_saddle2Id=dmt_saddle2.id_;

        if(isSaddle2(dmt_saddle2)){
          for(int i=0; i<3; ++i){
            simplexId_t vertexId;
            inputTriangulation_->getTriangleVertex(dmt_saddle2Id, i, vertexId);

            if(isPL[vertexId]!=2) continue;

            if(!allowBoundary and
               inputTriangulation_->isVertexOnBoundary(vertexId)) continue;

            if(pl2dmt_saddle2[vertexId]==-1){
              const simplexId_t pl_saddle2Id=vertexId;

              simplexId_t numberOfRemainingSaddles2=0;

              simplexId_t savedId=-1;
              const simplexId_t
                triangleNumber=inputTriangulation_->getVertexTriangleNumber(pl_saddle2Id);
              for(simplexId_t j=0; j<triangleNumber; ++j){
                simplexId_t triangleId;
                inputTriangulation_->getVertexTriangle(pl_saddle2Id, j,
                                                       triangleId);

                if(triangleId!=dmt_saddle2Id and isSaddle2(Cell(2,triangleId))
                   and isRemovableSaddle2[triangleId]){
                  ++numberOfRemainingSaddles2;
                  savedId=triangleId;
                }
              }

              if(!numberOfRemainingSaddles2){
                isRemovableSaddle2[dmt_saddle2Id]=false;
                pl2dmt_saddle2[vertexId]=dmt_saddle2Id;
                vpath.invalidate();
                break;
              }
              if(numberOfRemainingSaddles2==1){
                isRemovableSaddle2[dmt_saddle2Id]=false;
                isRemovableSaddle2[savedId]=false;
                pl2dmt_saddle2[vertexId]=savedId;
                break;
              }
            }
            else if(pl2dmt_saddle2[vertexId]==dmt_saddle2Id){
              vpath.invalidate();
              break;
            }
          }
        }
        else
          vpath.invalidate();
      }

      if(vpath.isValid_)
        reverseDescendingPathOnWall(path);
    }

    if(vpath.isValid_){
      // add persistence pair to collection if necessary
      if(CollectPersistencePairs and outputPersistencePairs_){
        const Cell& minSaddle1=criticalPoints[vpath.source_].cell_;
        const Cell& minSaddle2=criticalPoints[vpath.destination_].cell_;
        outputPersistencePairs_->push_back(std::make_tuple(minSaddle1, minSaddle2));
      }

      const int sourceId=vpath.source_;
      const int destinationId=vpath.destination_;

      // invalidate vpaths connected to source
      std::vector<int> newDestinationIds;
      CriticalPoint& source=criticalPoints[sourceId];
      for(auto& sourceVPathId : source.vpaths_){
        VPath& sourceVPath=vpaths[sourceVPathId];

        if(sourceVPath.isValid_ and sourceVPath.destination_!=destinationId){
          // save critical point
          const int newDestinationId=sourceVPath.destination_;
          newDestinationIds.push_back(newDestinationId);

          // clear vpath
          sourceVPath.invalidate();
        }
      }

      // invalidate vpaths connected to destination and save the critical
      // points to update
      std::vector<int> newSourceIds;
      CriticalPoint& destination=criticalPoints[destinationId];
      for(auto& destinationVPathId : destination.vpaths_){
        VPath& destinationVPath=vpaths[destinationVPathId];

        if(destinationVPath.isValid_ and destinationVPath.source_!=sourceId){
          // save critical point
          const int newSourceId=destinationVPath.source_;
          newSourceIds.push_back(newSourceId);

          CriticalPoint& newSource=criticalPoints[newSourceId];
          for(auto& newSourceVPathId : newSource.vpaths_){
            VPath& newSourceVPath=vpaths[newSourceVPathId];
            if(newSourceVPath.isValid_ and
               newSourceVPath.destination_!=destinationId){

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
      for(auto& newSourceId : newSourceIds){
        CriticalPoint& newSource=criticalPoints[newSourceId];
        const Cell& saddle1=newSource.cell_;

        std::set<simplexId_t> saddles2;
        const wallId_t savedWallId=wallId;
        getAscendingWall(wallId, saddle1, isVisited, nullptr, &saddles2);
        ++wallId;

        for(auto& saddle2Id : saddles2){
          const Cell saddle2(2,saddle2Id);

          const bool isMultiConnected=getDescendingPathThroughWall(savedWallId,
                                                                   saddle2, saddle1, isVisited, nullptr);
          if(isMultiConnected)
            continue;

          int newDestinationId=saddle2Index[saddle2Id];

          // connection to a new saddle2 (not present in the graph before)
          if(newDestinationId==-1){
            if(!isRemovableSaddle2[saddle2Id]) continue;

            const simplexId_t newCriticalPointId=criticalPoints.size();
            saddle2Index[saddle2Id]=newCriticalPointId;
            criticalPoints.push_back(CriticalPoint(saddle2));

            newDestinationId=newCriticalPointId;
          }

          CriticalPoint& newDestination=criticalPoints[newDestinationId];

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                              scalars);

          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(std::make_tuple(persistence,newVPathId,newSource.cell_.id_));
        }
      }

      // look at the gradient : get the links not predicted by the graph
      for(auto& newDestinationId : newDestinationIds){
        CriticalPoint& newDestination=criticalPoints[newDestinationId];
        const Cell& saddle2=newDestination.cell_;

        std::set<simplexId_t> saddles1;
        const wallId_t savedWallId=wallId;
        getDescendingWall(wallId, saddle2, isVisited, nullptr, &saddles1);
        ++wallId;

        for(auto& saddle1Id : saddles1){
          const Cell saddle1(1,saddle1Id);

          std::vector<Cell> path;
          const bool isMultiConnected=getAscendingPathThroughWall(savedWallId,
                                                                  saddle1, saddle2, isVisited, &path);
          if(isMultiConnected)
            continue;

          const simplexId_t newSourceId=saddle1Index[saddle1Id];

          if(newSourceId==-1)
            continue;

          CriticalPoint& newSource=criticalPoints[newSourceId];

          // check existence of the possibly newVPath in the graph
          bool alreadyExists=false;
          for(auto& newSourceVPathId : newSource.vpaths_){
            const VPath& newSourceVPath=vpaths[newSourceVPathId];

            if(newSourceVPath.isValid_ and
               newSourceVPath.destination_==newDestinationId){
              alreadyExists=true;
              break;
            }
          }

          if(alreadyExists)
            continue;

          // update vpaths
          const int newVPathId=vpaths.size();
          const dataType persistence=getPersistence<dataType>(saddle2, saddle1,
                                                              scalars);

          vpaths.push_back(VPath(true,-1,newSourceId,newDestinationId,-1,-1,persistence));

          // update criticalPoints
          newDestination.vpaths_.push_back(newVPathId);
          newSource.vpaths_.push_back(newVPathId);

          // update set
          S.insert(std::make_tuple(persistence,newVPathId,newSource.cell_.id_));
        }
      }
    }

    ++numberOfIterations;
  }

  {
    std::stringstream msg;
    msg << "[DiscreteGradient]  Processing of the vpaths :\t" <<
        t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template <typename dataType>
int DiscreteGradient::simplifySaddleSaddleConnections2(const
                                                       std::vector<std::pair<simplexId_t,char>>& criticalPoints,
                                                       const std::vector<char>& isPL,
                                                       const int iterationThreshold,
                                                       const bool allowBoundary,
                                                       const bool allowBruteForce,
                                                       const bool returnSaddleConnectors){
  Timer t;

  // Part 0 : get removable cells
  std::vector<char> isRemovableSaddle1;
  std::vector<simplexId_t> pl2dmt_saddle1(numberOfVertices_, -1);
  getRemovableSaddles1<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle1, pl2dmt_saddle1);

  std::vector<char> isRemovableSaddle2;
  std::vector<simplexId_t> pl2dmt_saddle2(numberOfVertices_, -1);
  getRemovableSaddles2<dataType>(criticalPoints, allowBoundary,
                                 isRemovableSaddle2, pl2dmt_saddle2);

  // Part 1 : initialization
  std::vector<VPath> vpaths;
  std::vector<CriticalPoint> dmt_criticalPoints;
  std::vector<simplexId_t> saddle1Index;
  std::vector<simplexId_t> saddle2Index;
  initializeSaddleSaddleConnections2<dataType>(isRemovableSaddle1,
                                               isRemovableSaddle2,
                                               allowBruteForce,
                                               vpaths,
                                               dmt_criticalPoints,
                                               saddle1Index,
                                               saddle2Index);

  // Part 2 : push the vpaths and order by persistence
  SaddleSaddleVPathComparator<dataType> cmp_f;
  std::set<std::tuple<dataType,int,simplexId_t>, SaddleSaddleVPathComparator<dataType>>
    S(cmp_f);
  orderSaddleSaddleConnections2<dataType>(vpaths, dmt_criticalPoints, S);

  // Part 3 : process the vpaths
  processSaddleSaddleConnections2<dataType>(iterationThreshold,
                                            isPL,
                                            allowBoundary,
                                            allowBruteForce,
                                            returnSaddleConnectors,
                                            S,
                                            pl2dmt_saddle1,
                                            pl2dmt_saddle2,
                                            isRemovableSaddle1,
                                            isRemovableSaddle2,
                                            vpaths,
                                            dmt_criticalPoints,
                                            saddle1Index,
                                            saddle2Index);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Saddle-Saddle pairs simplified in "
        << t.getElapsedTime() << " s, "<< threadNumber_ << " thread(s)." <<
        std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template<typename dataType>
int DiscreteGradient::filterSaddleConnectors(const bool allowBoundary){
  const bool allowBruteForce=false;
  const bool returnSaddleConnectors=true;

  // get the node type of a contour tree node (for compatibility with
  // ScalarFieldCriticalPoints)
  auto getNodeType=[&](const ftm::FTMTree_MT* tree, const ftm::Node* node){
    const int upDegree   = node->getNumberOfUpSuperArcs();
    const int downDegree = node->getNumberOfDownSuperArcs();
    const int degree = upDegree + downDegree;

    // saddle point
    if (degree > 1) {
      if (upDegree == 2 and downDegree == 1)
        return 2;
      else if (upDegree == 1 and downDegree == 2)
        return 1;
    }
      // local extremum
    else {
      if (upDegree)
        return 0;
      else
        return 3;
    }

    return -1;
  };

  std::vector<std::pair<simplexId_t,char>> cpset;

  simplexId_t* const offsets=static_cast<simplexId_t*>(inputOffsets_);
  const dataType* const scalars=static_cast<dataType*>(inputScalarField_);

  ftm::FTMTree contourTree;
  contourTree.setupTriangulation(inputTriangulation_, false);
  contourTree.setVertexScalars(scalars);
  contourTree.setTreeType(ftm::TreeType::Contour);
  contourTree.setVertexSoSoffsets(offsets);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.setSegmentation(false);
  contourTree.build<dataType>();
  ftm::FTMTree_MT* tree=contourTree.getTree(ftm::TreeType::Contour);

  const ftm::idVertex numberOfNodes=tree->getNumberOfNodes();
  for (ftm::idVertex nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
    const ftm::Node* node = tree->getNode(nodeId);
    const ftm::idVertex vertexId = node->getVertexId();

    cpset.push_back(std::make_pair(vertexId, getNodeType(tree,node)));
  }

  std::vector<char> isPL;
  getCriticalPointMap(cpset, isPL);

  simplifySaddleSaddleConnections1<dataType>(cpset, isPL,
                                             IterationThreshold, allowBoundary, allowBruteForce,
                                             returnSaddleConnectors);
  simplifySaddleSaddleConnections2<dataType>(cpset, isPL,
                                             IterationThreshold, allowBoundary, allowBruteForce,
                                             returnSaddleConnectors);

  return 0;
}

template<typename dataType>
int DiscreteGradient::reverseGradient(const
                                      std::vector<std::pair<simplexId_t,char>>& criticalPoints){
  Timer t;

  const bool allowBoundary=true;
  const bool returnSaddleConnectors=false;
  bool allowBruteForce=false;

  std::vector<char> isPL;
  getCriticalPointMap(criticalPoints, isPL);

  if(ReverseSaddleMaximumConnection)
    simplifySaddleMaximumConnections<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce);

  if(dimensionality_==3 and ReverseSaddleSaddleConnection){
    simplifySaddleSaddleConnections1<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce,
                                               returnSaddleConnectors);
    simplifySaddleSaddleConnections2<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce,
                                               returnSaddleConnectors);
  }

  allowBruteForce=true;

  if(ReverseSaddleMaximumConnection)
    simplifySaddleMaximumConnections<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce);

  if(dimensionality_==3 and ReverseSaddleSaddleConnection){
    simplifySaddleSaddleConnections1<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce,
                                               returnSaddleConnectors);
    simplifySaddleSaddleConnections2<dataType>(criticalPoints, isPL,
                                               IterationThreshold, allowBoundary, allowBruteForce,
                                               returnSaddleConnectors);
  }

  if(dimensionality_==3 and ReverseSaddleMaximumConnection and
     ReverseSaddleSaddleConnection and ReturnSaddleConnectors)
    filterSaddleConnectors<dataType>(allowBoundary);

  {
    std::stringstream msg;
    msg << "[DiscreteGradient] Gradient reversed in "
        << t.getElapsedTime() << " s. (" << threadNumber_
        << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return 0;
}

template<typename dataType>
int DiscreteGradient::reverseGradient(){
  std::vector<std::pair<simplexId_t,char>> criticalPoints;

  // get the PL critical points
  if(ReverseSaddleMaximumConnection or ReverseSaddleSaddleConnection){
    const simplexId_t* const offsets=static_cast<simplexId_t*>(inputOffsets_);
    std::vector<simplexId_t> sosOffsets(numberOfVertices_);
    for(simplexId_t i=0; i<numberOfVertices_; ++i)
      sosOffsets[i]=offsets[i];

    ScalarFieldCriticalPoints<dataType> scp;

    scp.setDebugLevel(debugLevel_);
    scp.setThreadNumber(threadNumber_);
    scp.setDomainDimension(dimensionality_);
    scp.setScalarValues(inputScalarField_);
    scp.setVertexNumber(numberOfVertices_);
    scp.setSosOffsets(&sosOffsets);
    scp.setupTriangulation(inputTriangulation_);
    scp.setOutput(&criticalPoints);

    scp.execute();
  }

  // print number of critical cells
  {
    // foreach dimension
    const int numberOfDimensions=getNumberOfDimensions();
    std::vector<simplexId_t> numberOfDMTCriticalPointsByDimension(numberOfDimensions,0);
    for(int i=0; i<numberOfDimensions; ++i){

      // foreach cell of that dimension
      const simplexId_t numberOfCells=getNumberOfCells(i);
      for(simplexId_t j=0; j<numberOfCells; ++j){
        const Cell cell(i,j);

        if(isCellCritical(cell))
          ++numberOfDMTCriticalPointsByDimension[i];
      }
    }

    std::vector<simplexId_t> numberOfPLInteriorCriticalPoints(numberOfDimensions,0);
    for(auto& criticalPoint : criticalPoints){
      const simplexId_t criticalPointId=criticalPoint.first;
      const char criticalPointType=criticalPoint.second;

      if(!inputTriangulation_->isVertexOnBoundary(criticalPointId) and
         criticalPointType!=-1)
        ++numberOfPLInteriorCriticalPoints[criticalPointType];
    }

    {
      std::stringstream msg;
      for(int i=0; i<numberOfDimensions; ++i){
        msg << "[DiscreteGradient] " << numberOfDMTCriticalPointsByDimension[i]
            << " " << i << "-cell(s)";
        msg << " and " << numberOfPLInteriorCriticalPoints[i] << " interior PL."
            << std::endl;
      }

      dMsg(std::cout, msg.str(), infoMsg);
    }
  }

  reverseGradient<dataType>(criticalPoints);

  return 0;
}

#endif // DISCRETEGRADIENT_TPL_H
