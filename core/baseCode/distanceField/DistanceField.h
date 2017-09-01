/// \ingroup baseCode
/// \class ttk::DistanceField
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date March 2016
///
/// \brief TTK processing package for distance field computation on PL 
/// manifolds.
///
/// This package takes a list of sources (a set of points with their global 
/// identifiers attached to them) and produces a distance field to the closest
/// source.
///
/// \b Related \b publication \n
/// "A note on two problems in connexion with graphs" \n
/// Edsger W. Dijkstra \n
/// Numerische Mathematik, 1959.
/// 
/// \sa vtkDistanceField.cpp %for a usage example.

#ifndef _DISTANCEFIELD_H
#define _DISTANCEFIELD_H

// base code includes
#include<Wrapper.h>
#include<Geometry.h>
#include<Triangulation.h>

// std includes
#include<limits>
#include<set>

namespace ttk{

  class DistanceField : public Debug{

    public:

      DistanceField();
      ~DistanceField();

      template <typename dataType>
        dataType getDistance(const int a, const int b) const;

      template <typename dataType>
        int execute() const;

      inline int setVertexNumber(int vertexNumber){
        vertexNumber_=vertexNumber;
        return 0;
      }

      inline int setSourceNumber(int sourceNumber){
        sourceNumber_=sourceNumber;
        return 0;
      }

      inline int setupTriangulation(Triangulation* triangulation){
        triangulation_=triangulation;
        if(triangulation_){
          triangulation_->preprocessVertexNeighbors();
        }
        return 0;
      }

      inline int setVertexIdentifierScalarFieldPointer(void* data){
        vertexIdentifierScalarFieldPointer_=data;
        return 0;
      }

      inline int setOutputScalarFieldPointer(void* data){
        outputScalarFieldPointer_=data;
        return 0;
      }

      inline int setOutputIdentifiers(void* data){
        outputIdentifiers_=data;
        return 0;
      }

      inline int setOutputSegmentation(void* data){
        outputSegmentation_=data;
        return 0;
      }

    protected:
      int vertexNumber_;
      int sourceNumber_;
      Triangulation* triangulation_;
      void* vertexIdentifierScalarFieldPointer_;
      void* outputScalarFieldPointer_;
      void* outputIdentifiers_;
      void* outputSegmentation_;
  };
}

template <typename dataType>
dataType DistanceField::getDistance(const int a, const int b) const{
  float p0[3];
  triangulation_->getVertexPoint(a,p0[0],p0[1],p0[2]);
  float p1[3];
  triangulation_->getVertexPoint(b,p1[0],p1[1],p1[2]);
  return Geometry::distance(p0,p1,3);
}

template <typename dataType>
int DistanceField::execute() const{
  int* identifiers=static_cast<int*>(vertexIdentifierScalarFieldPointer_);
  dataType* dist=static_cast<dataType*>(outputScalarFieldPointer_);
  int* origin=static_cast<int*>(outputIdentifiers_);
  int* seg=static_cast<int*>(outputSegmentation_);

  Timer t;

  fill(dist,dist+vertexNumber_,numeric_limits<dataType>::max());
  fill(origin,origin+vertexNumber_,-1);

  // get the sources
  set<int> isSource;
  for(int k=0; k<sourceNumber_; ++k)
    isSource.insert(identifiers[k]);
  vector<int> sources;
  for(auto s : isSource)
    sources.push_back(s);
  isSource.clear();

  // comparison lambda
  auto cmp=[](const pair<dataType,int>& a, const pair<dataType,int>& b){
    if(a.first != b.first) return a.first<b.first;
    else return a.second<b.second;
  };

  // prepare output
  vector<vector<dataType>> scalars(sources.size());
  for(auto& k : scalars)
    k.resize(vertexNumber_, numeric_limits<dataType>::max());

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i< (int) sources.size(); ++i){
    vector<bool> visited(vertexNumber_,false);
    set<pair<dataType,int>,decltype(cmp)> S(cmp);

    {
      const int s=sources[i];
      scalars[i][s]=0;
      visited[s]=true;

      const int neighborNumber=triangulation_->getVertexNeighborNumber(s);
      for(int k=0; k<neighborNumber; ++k){
        int neighbor;
        triangulation_->getVertexNeighbor(s,k,neighbor);
        if(!visited[neighbor]){
          scalars[i][neighbor]=getDistance<dataType>(s,neighbor);
          S.emplace(scalars[i][neighbor],neighbor);
        }
      }
    }

    while(!S.empty()){
      auto it=S.begin();
      const int vertex=it->second;

      if(!visited[vertex]){
        const dataType vertexScalar=scalars[i][vertex];
        const int neighborNumber=triangulation_->getVertexNeighborNumber(vertex);
        for(int k=0; k<neighborNumber; ++k){
          int neighbor;
          triangulation_->getVertexNeighbor(vertex,k,neighbor);

          const dataType delta=getDistance<dataType>(vertex,neighbor);
          if(vertexScalar+delta<scalars[i][neighbor]){
            scalars[i][neighbor]=vertexScalar+delta;
            S.emplace(scalars[i][neighbor],neighbor);
          }
        }
        visited[vertex]=true;
      }
      S.erase(it);
    }
  }

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int k=0; k<vertexNumber_; ++k){
    for(unsigned int i=0; i<sources.size(); ++i){
      if(i==0 or dist[k]>scalars[i][k]){
        dist[k]=scalars[i][k];
        origin[k]=sources[i];
        seg[k]=i;
      }
    }
  }

  {
    stringstream msg;
    msg << "[DistanceField] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // DISTANCEFIELD_H
