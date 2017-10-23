/// \ingroup base
/// \class ttk::TopologicalSimplification
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK processing package for the topological simplification of scalar 
/// data.
///
/// Given an input scalar field and a list of critical points to remove, this 
/// class minimally edits the scalar field such that the listed critical points
/// disappear. This procedure is useful to speedup subsequent topological data 
/// analysis when outlier critical points can be easily identified. It is 
/// also useful for data simplification.
///
/// \b Related \b publication \n
/// "Generalized Topological Simplification of Scalar Fields on Surfaces" \n
/// Julien Tierny, Valerio Pascucci \n
/// Proc. of IEEE VIS 2012.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2012.
///
/// \sa vtkTopologicalSimplification.cpp %for a usage example.

#ifndef _TOPOLOGICALSIMPLIFICATION_H
#define _TOPOLOGICALSIMPLIFICATION_H

// base code includes
#include                  <Wrapper.h>

#include<Triangulation.h>
#include<set>
#include<tuple>
#include<type_traits>

namespace ttk{
  
  struct SweepCmp{
    private :
      bool isIncreasingOrder_;

    public:
      SweepCmp():
        isIncreasingOrder_{}
      {}

      SweepCmp(bool isIncreasingOrder):
        isIncreasingOrder_{isIncreasingOrder}
      {}

      int setIsIncreasingOrder(bool isIncreasingOrder){
        isIncreasingOrder_=isIncreasingOrder;
        return 0;
      }

      template <typename dataType>
        bool operator() (const tuple<dataType,int,int> &v0,
            const tuple<dataType,int,int> &v1) const{
          if(isIncreasingOrder_){
            return (get<0>(v0) < get<0>(v1) or
                (get<0>(v0) == get<0>(v1) and get<1>(v0) < get<1>(v1)));
          }
          else{
            return (get<0>(v0) > get<0>(v1) or
                (get<0>(v0) == get<0>(v1) and get<1>(v0) > get<1>(v1)));
          }
        };
  };

  class TopologicalSimplification : public Debug{

    public:

      TopologicalSimplification();

      ~TopologicalSimplification();

      template <typename dataType>
        bool isLowerThan(int a, int b, dataType* scalars, int* offsets) const;

      template <typename dataType>
        bool isHigherThan(int a, int b, dataType* scalars, int* offsets) const;

      template <typename dataType>
        int getCriticalType(int vertexId, dataType* scalars, int* offsets) const;

      template <typename dataType>
        int getCriticalPoints(dataType* scalars,
            int* offsets,
            vector<int>& minList,
            vector<int>& maxList) const;

      template <typename dataType>
        int getCriticalPoints(dataType* scalars,
            int* offsets,
            vector<int>& minList,
            vector<int>& maxList,
            vector<bool>& blackList) const;

      template <typename dataType>
        int addPerturbation(dataType* scalars, int* offsets) const;

      template <typename dataType>
        int execute() const;

      inline int setupTriangulation(Triangulation* triangulation){
        triangulation_=triangulation;
        if(triangulation_){
          vertexNumber_ = triangulation_->getNumberOfVertices();
          triangulation_->preprocessVertexNeighbors();
        }
        return 0;
      }

      inline int setVertexNumber(int vertexNumber){
        vertexNumber_=vertexNumber;
        return 0;
      }

      inline int setConstraintNumber(int constraintNumber){
        constraintNumber_=constraintNumber;
        return 0;
      }

      inline int setInputScalarFieldPointer(void *data){
        inputScalarFieldPointer_=data;
        return 0;
      }

      inline int setVertexIdentifierScalarFieldPointer(void* data){
        vertexIdentifierScalarFieldPointer_=data;
        return 0;
      }

      inline int setInputOffsetScalarFieldPointer(void* data){
        inputOffsetScalarFieldPointer_=data;
        return 0;
      }

      inline int setConsiderIdentifierAsBlackList(bool onOff){
        considerIdentifierAsBlackList_=onOff;
        return 0;
      }

      inline int setAddPerturbation(bool onOff){
        addPerturbation_=onOff;
        return 0;
      }

      inline int setOutputScalarFieldPointer(void *data){
        outputScalarFieldPointer_=data;
        return 0;
      }

      inline int setOutputOffsetScalarFieldPointer(void *data){
        outputOffsetScalarFieldPointer_=data;
        return 0;
      }
    protected:

      Triangulation* triangulation_;
      int vertexNumber_;
      int constraintNumber_;
      void* inputScalarFieldPointer_;
      void* vertexIdentifierScalarFieldPointer_;
      void* inputOffsetScalarFieldPointer_;
      bool considerIdentifierAsBlackList_;
      bool addPerturbation_;
      void* outputScalarFieldPointer_;
      void* outputOffsetScalarFieldPointer_;
  };
}

// if the package is a pure template typename, uncomment the following line
// #include                  <TopologicalSimplification.cpp>

template <typename dataType>
bool TopologicalSimplification::isLowerThan(int a, int b, dataType* scalars, int* offsets) const{
  return (scalars[a]<scalars[b] or
      (scalars[a]==scalars[b] and offsets[a]<offsets[b]));
}

template <typename dataType>
bool TopologicalSimplification::isHigherThan(int a, int b, dataType* scalars, int* offsets) const{
  return (scalars[a]>scalars[b] or
      (scalars[a]==scalars[b] and offsets[a]>offsets[b]));
}

template <typename dataType>
int TopologicalSimplification::getCriticalType(int vertex, dataType* scalars, int* offsets) const{
  bool isMinima{true};
  bool isMaxima{true};
  int neighborNumber=triangulation_->getVertexNeighborNumber(vertex);
  for(int i=0; i<neighborNumber; ++i){
    int neighbor;
    triangulation_->getVertexNeighbor(vertex,i,neighbor);
    
    if(isLowerThan<dataType>(neighbor,vertex,scalars,offsets)) isMinima=false;
    if(isHigherThan<dataType>(neighbor,vertex,scalars,offsets)) isMaxima=false;
    if(!isMinima and !isMaxima){
      return 0;
    }
  }

  if(isMinima) return -1;
  if(isMaxima) return 1;

  return 0;
}

template <typename dataType>
int TopologicalSimplification::getCriticalPoints(dataType* scalars,
    int* offsets,
    vector<int>& minima,
    vector<int>& maxima) const{
	vector<int> type(vertexNumber_, 0);
  
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
  for(int k=0; k<vertexNumber_; ++k)
    type[k]=getCriticalType<dataType>(k,scalars,offsets);

  for(int k=0; k<vertexNumber_; ++k){
    if(type[k]<0) minima.push_back(k);
    else if(type[k]>0) maxima.push_back(k);
  }
  return 0;
}

template <typename dataType>
int TopologicalSimplification::getCriticalPoints(dataType* scalars,
    int* offsets,
    vector<int>& minima,
    vector<int>& maxima,
    vector<bool>& extrema) const{
  vector<int> type(vertexNumber_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int k=0; k<vertexNumber_; ++k){
    if(considerIdentifierAsBlackList_ xor extrema[k]){
      type[k]=getCriticalType<dataType>(k,scalars,offsets);
    }
  }

  for(int k=0; k<vertexNumber_; ++k){
    if(type[k]<0) minima.push_back(k);
    else if(type[k]>0) maxima.push_back(k);
  }
  return 0;
}

template <typename dataType>
int TopologicalSimplification::addPerturbation(dataType* scalars, int* offsets) const{
  dataType epsilon{};

  if(is_same<dataType,double>::value) epsilon=pow10(1-DBL_DIG);
  else if(is_same<dataType,float>::value) epsilon=pow10(1-FLT_DIG);
  else return -1;

  vector<tuple<dataType,int,int>> perturbation(vertexNumber_);
  for(int i=0; i<vertexNumber_; ++i){
    get<0>(perturbation[i])=scalars[i];
    get<1>(perturbation[i])=offsets[i];
    get<2>(perturbation[i])=i;
  }

  SweepCmp cmp(true);
  sort(perturbation.begin(), perturbation.end(), cmp);

  for(int i=0; i<vertexNumber_; ++i){
    if(i){
      if(get<0>(perturbation[i]) <= get<0>(perturbation[i-1]))
        get<0>(perturbation[i])=get<0>(perturbation[i-1]) + epsilon;
    }
    scalars[get<2>(perturbation[i])]=get<0>(perturbation[i]);
  }

  return 0;
}

template <typename dataType>
int TopologicalSimplification::execute() const{
  
  // get input data
  dataType* inputScalars=static_cast<dataType*>(inputScalarFieldPointer_);
  dataType* scalars=static_cast<dataType*>(outputScalarFieldPointer_);
  int* identifiers=static_cast<int*>(vertexIdentifierScalarFieldPointer_);
  int* inputOffsets=static_cast<int*>(inputOffsetScalarFieldPointer_);
  int* offsets=static_cast<int*>(outputOffsetScalarFieldPointer_);

  Timer t;

  // pre-processing
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int k=0; k<vertexNumber_; ++k){
    scalars[k]=inputScalars[k];
    if(isnan(scalars[k]))
      scalars[k] = 0;
    offsets[k]=inputOffsets[k];
  }

  // get the user extremum list
  vector<bool> extrema(vertexNumber_, false);
  for(int k=0; k<constraintNumber_; ++k){
    const int identifierId=identifiers[k];

#ifndef TTK_ENABLE_KAMIKAZE
    if(identifierId>=0 and identifierId<vertexNumber_)
#endif
      extrema[identifierId]=true;
  }

  vector<int> authorizedMinima;
  vector<int> authorizedMaxima;
  vector<bool> authorizedExtrema(vertexNumber_, false);
  
  getCriticalPoints<dataType>(
    scalars, offsets,
    authorizedMinima,
    authorizedMaxima,
    extrema);
  
  {
    stringstream msg;
    msg << "[TopologicalSimplification] Maintaining "
      << constraintNumber_ 
      << " constraints ("
      << authorizedMinima.size() << " minima and "
      << authorizedMaxima.size() << " maxima)." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  // declare the tuple-comparison functor
  SweepCmp cmp;

  // processing
  int iteration{};
  for(int i=0; i<vertexNumber_; ++i){
   
    {
      stringstream msg;
      msg << "[TopologicalSimplification] Starting simplifying iteration #"
        << i << "..." << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
    
    for(int j=0; j<2; ++j){
      
      bool isIncreasingOrder=!j;

      cmp.setIsIncreasingOrder(isIncreasingOrder);
      set<tuple<dataType,int,int>, decltype(cmp)> sweepFront(cmp);
      vector<bool> visitedVertices(vertexNumber_, false);
      vector<int> adjustmentSequence(vertexNumber_);

      // add the seeds
      if(isIncreasingOrder){
        for(int k : authorizedMinima){
          authorizedExtrema[k]=true;
          sweepFront.emplace(scalars[k],offsets[k],k);
          visitedVertices[k]=true;
        }
      }
      else{
        for(int k : authorizedMaxima){
          authorizedExtrema[k]=true;
          sweepFront.emplace(scalars[k],offsets[k],k);
          visitedVertices[k]=true;
        }
      }
      
      // growth by neighborhood of the seeds
      int adjustmentPos = 0;
      do{
        auto front=sweepFront.begin();
        if(front==sweepFront.end()) return -1;

        int vertexId=get<2>(*front);
        sweepFront.erase(front);

        int neighborNumber=triangulation_->getVertexNeighborNumber(vertexId);
        for(int k=0; k<neighborNumber; ++k){
          int neighbor;
          triangulation_->getVertexNeighbor(vertexId,k,neighbor);
          if(!visitedVertices[neighbor]){
            sweepFront.emplace(scalars[neighbor],offsets[neighbor],neighbor);
            visitedVertices[neighbor]=true;
          }
        }
        adjustmentSequence[adjustmentPos]=vertexId;
        ++adjustmentPos;
      }while(!sweepFront.empty());

      // save offsets and rearrange scalars
      int offset = (isIncreasingOrder ? 0 : vertexNumber_ + 1);
      
      for(int k=0; k<vertexNumber_; ++k){
        
        if(isIncreasingOrder){
          if(k and scalars[adjustmentSequence[k]] <= scalars[adjustmentSequence[k-1]])
            scalars[adjustmentSequence[k]]=scalars[adjustmentSequence[k-1]];
          ++offset;
        }
        else{
          if(k and scalars[adjustmentSequence[k]] >= scalars[adjustmentSequence[k-1]])
            scalars[adjustmentSequence[k]]=scalars[adjustmentSequence[k-1]];
          --offset;
        }
        offsets[adjustmentSequence[k]]=offset;
      }
    }

    // test convergence
    bool needForMoreIterations{false};
    vector<int> minima;
    vector<int> maxima;
    getCriticalPoints<dataType>(scalars,offsets,minima,maxima);
    
    if(maxima.size() > authorizedMaxima.size()) needForMoreIterations=true;
    if(minima.size() > authorizedMinima.size()) needForMoreIterations=true;
    
    {
      stringstream msg;
      msg << "[TopologicalSimplification] Current status: "
        << minima.size() << " minima, " 
        << maxima.size() << " maxima." << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
    }
    
    if(!needForMoreIterations){
      for(int k : minima){
        if(!authorizedExtrema[k]){
          needForMoreIterations=true;
          break;
        }
      }
    }
    if(!needForMoreIterations){
      for(int k : maxima){
        if(!authorizedExtrema[k]){
          needForMoreIterations=true;
          break;
        }
      }
    }

    // optional adding of perturbation
    if(addPerturbation_) addPerturbation<dataType>(scalars, offsets);

    ++iteration;
    if(!needForMoreIterations) break;
  }

  {
    stringstream msg;
    msg << "[TopologicalSimplification] Scalar field simplified"
      << " in " << t.getElapsedTime() << " s. (" << threadNumber_
      << " threads(s), " << iteration << " ite.)."
      << endl;
    dMsg(cout,msg.str(),timeMsg);
  }

  return 0;
}
#endif // TOPOLOGICALSIMPLIFICATION_H
