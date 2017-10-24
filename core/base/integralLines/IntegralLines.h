/// \ingroup base
/// \class ttk::IntegralLines
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK processing package for the computation of edge-based integral 
/// lines of the gradient of an input scalar field defined on a PL manifold. 
/// 
/// Given a list of sources, the package produces forward or backward integral
/// lines along the edges of the input triangulation.
///
/// \sa vtkIntegralLines.cpp %for a usage example.

#ifndef _DISCRETESTREAMLINE_H
#define _DISCRETESTREAMLINE_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>

// std includes
#include<unordered_set>

namespace ttk{
  enum Direction{
    Forward=0,
    Backward
  };

  class IntegralLines : public Debug{

    public:

      IntegralLines();
      ~IntegralLines();

      template<typename dataType>
        inline float getDistance(const int& a, const int& b) const{
          float p0[3];
          triangulation_->getVertexPoint(a,p0[0],p0[1],p0[2]);
          float p1[3];
          triangulation_->getVertexPoint(b,p1[0],p1[1],p1[2]);

          return Geometry::distance(p0,p1,3);
        }

      template<typename dataType>
        inline float getGradient(const int& a, const int& b, dataType* scalars) const{
          return fabs(scalars[b]-scalars[a])/getDistance<dataType>(a,b);
        }

      template<typename dataType>
        int execute() const;

      template<typename dataType, class Compare>
        int execute(Compare cmp) const;

      inline int setVertexNumber(const int &vertexNumber){
        vertexNumber_=vertexNumber;
        return 0;
      }

      inline int setSeedNumber(const int &seedNumber){
        seedNumber_=seedNumber;
        return 0;
      }

      inline int setDirection(int direction){
        direction_=direction;
        return 0;
      }

      inline int setupTriangulation(Triangulation* triangulation){
        triangulation_=triangulation;
        if(triangulation_){
          triangulation_->preprocessVertexNeighbors();
        }
        return 0;
      }

      inline int setInputScalarField(void *data){
        inputScalarField_=data;
        return 0;
      }

      inline int setInputOffsets(void* data){
        inputOffsets_=data;
        return 0;
      }

      inline int setVertexIdentifierScalarField(void *data){
        vertexIdentifierScalarField_=data;
        return 0;
      }

      inline int setOutputTrajectories(vector<vector<int>>* trajectories){
        outputTrajectories_=trajectories;
        return 0;
      }

    protected:

      int vertexNumber_;
      int seedNumber_;
      int direction_;
      Triangulation* triangulation_;
      void* inputScalarField_;
      void* inputOffsets_;
      void* vertexIdentifierScalarField_;
      vector<vector<int>>* outputTrajectories_;
  };
}

template<typename dataType>
int IntegralLines::execute() const{
  int* offsets=static_cast<int*>(inputOffsets_);
  int* identifiers=static_cast<int*>(vertexIdentifierScalarField_);
  dataType* scalars=static_cast<dataType*>(inputScalarField_);
  vector<vector<int>>* trajectories=outputTrajectories_;

  Timer t;

  // get the seeds
  unordered_set<int> isSeed;
  for(int k=0; k<seedNumber_; ++k)
    isSeed.insert(identifiers[k]);
  vector<int> seeds;
  for(auto k : isSeed)
    seeds.push_back(k);
  isSeed.clear();

  trajectories->resize(seeds.size());
  for(unsigned int i=0; i<seeds.size(); ++i){
    int v{seeds[i]};
    (*trajectories)[i].push_back(v);

    bool isMax{};
    while(!isMax){
      int vnext{-1};
      float fnext=numeric_limits<float>::min();
      int neighborNumber=triangulation_->getVertexNeighborNumber(v);
      bool isLocalMax=true;
      bool isLocalMin=true;
      for(int k=0; k<neighborNumber; ++k){
        int n;
        triangulation_->getVertexNeighbor(v,k,n);

        if(scalars[n]<=scalars[v]) isLocalMax=false;
        if(scalars[n]>=scalars[v]) isLocalMin=false;

        if((direction_==static_cast<int>(Direction::Forward)) xor (scalars[n]<scalars[v])){
          const float f=getGradient<dataType>(v,n,scalars);
          if(f>fnext){
            vnext=n;
            fnext=f;
          }
        }
      }

      if(vnext==-1 and !isLocalMax and !isLocalMin){
        int onext=-1;
        for(int k=0; k<neighborNumber; ++k){
          int n;
          triangulation_->getVertexNeighbor(v,k,n);

          if(scalars[n]==scalars[v]){
            const int o=offsets[n];
            if((direction_==static_cast<int>(Direction::Forward)) xor (o<offsets[v])){
              if(o>onext){
                vnext=n;
                onext=o;
              }
            }
          }
        }
      }

      if(vnext==-1) isMax=true;
      else{
        v=vnext;
        (*trajectories)[i].push_back(v);
      }
    }
  }

  {
    stringstream msg;
    msg << "[IntegralLines] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

template<typename dataType, class Compare>
int IntegralLines::execute(Compare cmp) const{
  int* offsets=static_cast<int*>(inputOffsets_);
  int* identifiers=static_cast<int*>(vertexIdentifierScalarField_);
  dataType* scalars=static_cast<dataType*>(inputScalarField_);
  vector<vector<int>>* trajectories=outputTrajectories_;

  Timer t;

  // get the seeds
  unordered_set<int> isSeed;
  for(int k=0; k<seedNumber_; ++k)
    isSeed.insert(identifiers[k]);
  vector<int> seeds;
  for(auto k : isSeed)
    seeds.push_back(k);
  isSeed.clear();

  trajectories->resize(seeds.size());
  for(unsigned int i=0; i<seeds.size(); ++i){
    int v{seeds[i]};
    (*trajectories)[i].push_back(v);

    bool isMax{};
    while(!isMax){
      int vnext{-1};
      float fnext=numeric_limits<float>::min();
      int neighborNumber=triangulation_->getVertexNeighborNumber(v);
      bool isLocalMax=true;
      bool isLocalMin=true;
      for(int k=0; k<neighborNumber; ++k){
        int n;
        triangulation_->getVertexNeighbor(v,k,n);

        if(scalars[n]<=scalars[v]) isLocalMax=false;
        if(scalars[n]>=scalars[v]) isLocalMin=false;

        if((direction_==static_cast<int>(Direction::Forward)) xor (scalars[n]<scalars[v])){
          const float f=getGradient<dataType>(v,n,scalars);
          if(f>fnext){
            vnext=n;
            fnext=f;
          }
        }
      }

      if(vnext==-1 and !isLocalMax and !isLocalMin){
        int onext=-1;
        for(int k=0; k<neighborNumber; ++k){
          int n;
          triangulation_->getVertexNeighbor(v,k,n);

          if(scalars[n]==scalars[v]){
            const int o=offsets[n];
            if((direction_==static_cast<int>(Direction::Forward)) xor (o<offsets[v])){
              if(o>onext){
                vnext=n;
                onext=o;
              }
            }
          }
        }
      }

      if(vnext==-1) isMax=true;
      else{
        v=vnext;
        (*trajectories)[i].push_back(v);

        if(cmp(v)) isMax=true;
      }
    }
  }

  {
    stringstream msg;
    msg << "[IntegralLines] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}


#endif // DISCRETESTREAMLINE_H
