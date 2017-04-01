/// \ingroup baseCode
/// \class ttk::PersistenceDiagram
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2016.
///
/// \brief TTK processing package for the computation of persistence diagrams.
///
/// This package computes the persistence diagram of the extremum-saddle pairs 
/// of an input scalar field. The X-coordinate of each pair corresponds to its
/// birth, while its smallest and highest Y-coordinates correspond to its birth
/// and death respectively. 
///
/// In practice, each extremity of a persistence pair is represented by its 
/// vertexId and critical type. Based on that, the persistence of the pair 
/// and its 2D embedding can easily be obtained.
///
/// Persistence diagrams are useful and stable concise representations of the 
/// topological features of a data-set. It is useful to fine-tune persistence 
/// thresholds for topological simplification or for fast similarity 
/// estimations for instance. 
///
/// \b Related \b publication \n
/// "Computational Topology: An Introduction" \n
/// Herbert Edelsbrunner and John Harer \n
/// American Mathematical Society, 2010
///
/// \sa vtkPersistenceDiagram.cpp %for a usage example.

#ifndef _PERSISTENCEDIAGRAM_H
#define _PERSISTENCEDIAGRAM_H

// base code includes
#include<Wrapper.h>
#include<Triangulation.h>
#include<DiscreteGradient.h>
#include<ContourForests.h>

namespace ttk{

  class PersistenceDiagram : public Debug{

    public:

      PersistenceDiagram();
      ~PersistenceDiagram();

      inline int setComputeSaddleConnectors(bool state){
        ComputeSaddleConnectors=state;
        return 0;
      }

      NodeType getNodeType(MergeTree* tree,
          TreeType treeType,
          const int vertexId) const;

      template <typename scalarType>
        int sortPersistenceDiagram(
            vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>>& diagram,
            scalarType* scalars) const;

      template <typename scalarType>
        int computeCTPersistenceDiagram(ContourForests& tree,
            const vector<tuple<idVertex, idVertex, scalarType, bool>>& pairs,
            vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType, idVertex>>& diagram,
            scalarType* scalars) const;

      template <class scalarType>
        int execute() const;

      inline int setDMTPairs(vector<tuple<Cell,Cell>>* data){
        dmt_pairs=data;
        return 0;
      }

      inline int setupTriangulation(Triangulation* data){
        triangulation_ = data;
        if(triangulation_){
          ContourForests contourTree;
          contourTree.setupTriangulation(triangulation_);

          triangulation_->preprocessBoundaryVertices();
        }
        return 0;
      }

      inline int setInputScalars(void* data){
        inputScalars_ = data;
        return 0;
      }

      inline int setInputOffsets(void* data){
        inputOffsets_ = data;
        return 0;
      }

      inline int setOutputCTDiagram(void* data){
        CTDiagram_=data;
        return 0;
      }

    protected:

      vector<tuple<Cell,Cell>>* dmt_pairs;

      bool ComputeSaddleConnectors;

      Triangulation* triangulation_;
      void* inputScalars_;
      void* inputOffsets_;
      void* CTDiagram_;
  };
}

template <typename scalarType>
int PersistenceDiagram::sortPersistenceDiagram(
    vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>>& diagram,
    scalarType* scalars) const{
  auto cmp=[scalars](const tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>& a,
      const tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>& b){
    return scalars[get<0>(a)] < scalars[get<0>(b)];
  };

  std::sort(diagram.begin(), diagram.end(), cmp);

  return 0;
}

template <typename scalarType>
int PersistenceDiagram::computeCTPersistenceDiagram(ContourForests& tree,
    const vector<tuple<idVertex, idVertex, scalarType, bool>>& pairs,
    vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>>& diagram,
    scalarType* scalars) const{
  const idVertex numberOfPairs=pairs.size();
  diagram.resize(numberOfPairs);
  for(idVertex i=0; i<numberOfPairs; ++i){
    const idVertex v0=get<0>(pairs[i]);
    const idVertex v1=get<1>(pairs[i]);
    const scalarType persistenceValue=get<2>(pairs[i]);
    const bool type=get<3>(pairs[i]);

    get<4>(diagram[i])=persistenceValue;
    if(type==true){
      get<0>(diagram[i])=v0;
      get<1>(diagram[i])=getNodeType(tree.getJoinTree(),TreeType::Join,v0);
      get<2>(diagram[i])=v1;
      get<3>(diagram[i])=getNodeType(tree.getJoinTree(),TreeType::Join,v1);
      get<5>(diagram[i])=0;
    }
    else{
      get<0>(diagram[i])=v1;
      get<1>(diagram[i])=getNodeType(tree.getSplitTree(),TreeType::Split,v1);
      get<2>(diagram[i])=v0;
      get<3>(diagram[i])=getNodeType(tree.getSplitTree(),TreeType::Split,v0);
      get<5>(diagram[i])=2;
    }
  }

  return 0;
}

template <typename scalarType>
int PersistenceDiagram::execute() const{
  // get data
  vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>>& CTDiagram=
    *static_cast<vector<tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex>>*>(CTDiagram_);
  scalarType* scalars=static_cast<scalarType*>(inputScalars_);
  int* offsets=static_cast<int*>(inputOffsets_);

  const idVertex numberOfVertices=triangulation_->getNumberOfVertices();
  // convert offsets into a valid format for contour forests
  vector<idVertex> voffsets(numberOfVertices);
  std::copy(offsets,offsets+numberOfVertices,voffsets.begin());

  // get contour tree
  ContourForests contourTree;
  contourTree.setupTriangulation(triangulation_, false);
  contourTree.setVertexScalars(inputScalars_);
  contourTree.setTreeType(TreeType::Join);
  contourTree.setVertexSoSoffsets(voffsets);
  contourTree.setLessPartition(true);
  // for now, only one thread is supported
  //contourTree_->setThreadNumber(threadNumber_);
  contourTree.setThreadNumber(1);
  contourTree.build<scalarType>();

  // get persistence pairs
  vector<tuple<idVertex,idVertex,scalarType>> JTPairs;
  vector<tuple<idVertex,idVertex,scalarType>> STPairs;
#ifdef withOpenMP
#pragma omp parallel sections
#endif
  {
#ifdef withOpenMP
#pragma omp section
#endif
    contourTree.getJoinTree()->computePersistencePairs<scalarType>(JTPairs);
#ifdef withOpenMP
#pragma omp section
#endif
    contourTree.getSplitTree()->computePersistencePairs<scalarType>(STPairs);
  }

  // merge pairs
  vector<tuple<idVertex, idVertex, scalarType, bool>> 
    CTPairs(JTPairs.size()+STPairs.size());
  const idVertex JTSize=JTPairs.size();
  for(idVertex i=0; i<JTSize; ++i){
    const auto& x=JTPairs[i];
    CTPairs[i]=make_tuple(get<0>(x),get<1>(x),get<2>(x),true);
  }
  const idVertex STSize=STPairs.size();
  for(idVertex i=0; i<STSize; ++i){
    const auto& x=STPairs[i];
    CTPairs[JTSize+i]=make_tuple(get<0>(x),get<1>(x),get<2>(x),false);
  }

  // remove the last pair which is present two times (global extrema pair)
  {
    auto cmp=[](const tuple<idVertex,idVertex,scalarType,bool>& a,
        const tuple<idVertex,idVertex,scalarType,bool>& b){
      return get<2>(a) < get<2>(b);
    };

    std::sort(CTPairs.begin(), CTPairs.end(), cmp);
    CTPairs.erase(CTPairs.end()-1);
  }

  // get the saddle-saddle pairs
  vector<tuple<int,int,scalarType>> pl_saddleSaddlePairs;
  const int dimensionality=triangulation_->getDimensionality();
  if(dimensionality==3 and ComputeSaddleConnectors){
    // get original list of critical points
    vector<pair<int,char>> pl_criticalPoints;
    {
      const int* const offsets=static_cast<int*>(inputOffsets_);
      vector<int> sosOffsets(numberOfVertices);
      for(int i=0; i<numberOfVertices; ++i)
        sosOffsets[i]=offsets[i];

      ScalarFieldCriticalPoints<scalarType> scp;

      scp.setDebugLevel(debugLevel_);
      scp.setThreadNumber(threadNumber_);
      scp.setDomainDimension(dimensionality);
      scp.setScalarValues(inputScalars_);
      scp.setVertexNumber(numberOfVertices);
      scp.setSosOffsets(&sosOffsets);
      scp.setupTriangulation(triangulation_);
      scp.setOutput(&pl_criticalPoints);

      scp.execute();
    }

    // build rejecting list
    vector<char> isRejected(numberOfVertices, false);
    for(const auto& i : CTPairs){
      const int v0=get<0>(i);
      const int v1=get<1>(i);
      isRejected[v0]=true;
      isRejected[v1]=true;
    }

    // filter the critical points according to the filtering list and boundary condition
    vector<char> isSaddle1(numberOfVertices, false);
    vector<char> isSaddle2(numberOfVertices, false);
    vector<pair<int,char>> pl_filteredCriticalPoints;
    for(const auto& i : pl_criticalPoints){
      const int vertexId=i.first;
      const char type=i.second;
      if(!isRejected[vertexId]){
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
      DiscreteGradient discreteGradient;
      discreteGradient.setDebugLevel(debugLevel_);
      discreteGradient.setThreadNumber(threadNumber_);
      discreteGradient.setupTriangulation(triangulation_);
      discreteGradient.setIterationThreshold(-1);
      discreteGradient.setReverseSaddleMaximumConnection(true);
      discreteGradient.setReverseSaddleSaddleConnection(true);
      discreteGradient.setAllowReversingWithNonRemovable(false);
      discreteGradient.setCollectPersistencePairs(false);
      discreteGradient.setInputScalarField(inputScalars_);
      discreteGradient.setInputOffsets(inputOffsets_);
      discreteGradient.buildGradient<scalarType>();
      discreteGradient.reverseGradient<scalarType>(pl_criticalPoints);

      // collect saddle-saddle connections
      discreteGradient.setReverseSaddleMaximumConnection(false);
      discreteGradient.setCollectPersistencePairs(true);
      discreteGradient.setOutputPersistencePairs(&dmt_pairs);
      discreteGradient.reverseGradient<scalarType>(pl_filteredCriticalPoints);
    }

    // transform DMT pairs into PL pairs
    for(const auto& i : dmt_pairs){
      const Cell& saddle1=get<0>(i);
      const Cell& saddle2=get<1>(i);

      int v0=-1;
      for(int j=0; j<2; ++j){
        int vertexId;
        triangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

        if(isSaddle1[vertexId]){
          v0=vertexId;
          break;
        }
      }
      if(v0==-1){
        scalarType scalar{};
        for(int j=0; j<2; ++j){
          int vertexId;
          triangulation_->getEdgeVertex(saddle1.id_,j,vertexId);
          const scalarType vertexScalar=scalars[vertexId];

          if(!j or scalar>vertexScalar){
            v0=vertexId;
            scalar=vertexScalar;
          }
        }
      }

      int v1=-1;
      for(int j=0; j<3; ++j){
        int vertexId;
        triangulation_->getTriangleVertex(saddle2.id_, j, vertexId);

        if(isSaddle2[vertexId]){
          v1=vertexId;
          break;
        }
      }
      if(v1==-1){
        scalarType scalar{};
        for(int j=0; j<3; ++j){
          int vertexId;
          triangulation_->getTriangleVertex(saddle2.id_,j,vertexId);
          const scalarType vertexScalar=scalars[vertexId];

          if(!j or scalar<vertexScalar){
            v1=vertexId;
            scalar=vertexScalar;
          }
        }
      }

      const scalarType persistence=scalars[v1]-scalars[v0];

      if(v0!=-1 and v1!=-1 and persistence>=0)
        pl_saddleSaddlePairs.push_back(make_tuple(v0,v1,persistence));
    }
  }

  // get persistence diagrams
  computeCTPersistenceDiagram<scalarType>(
      contourTree,CTPairs, CTDiagram, scalars);

  // add saddle-saddle pairs to the diagram if needed
  if(dimensionality==3 and ComputeSaddleConnectors){
    for(const auto& i : pl_saddleSaddlePairs){
      const idVertex v0=get<0>(i);
      const idVertex v1=get<1>(i);
      const scalarType persistenceValue=get<2>(i);

      tuple<idVertex,NodeType,idVertex,NodeType,scalarType,idVertex> t;

      get<0>(t)=v0;
      get<1>(t)=NodeType::Saddle1;
      get<2>(t)=v1;
      get<3>(t)=NodeType::Saddle2;
      get<4>(t)=persistenceValue;
      get<5>(t)=1;

      CTDiagram.push_back(t);
    }
  }

  // finally sort the diagram
  sortPersistenceDiagram(CTDiagram, scalars);

  return 0;
}

#endif // PERSISTENCEDIAGRAM_H
