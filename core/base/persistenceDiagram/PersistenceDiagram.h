/// \ingroup base
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
#include<FTMTreePP.h>
#include<MorseSmaleComplex3D.h>

namespace ttk{

  class PersistenceDiagram : public Debug{

    public:

      PersistenceDiagram();
      ~PersistenceDiagram();

      inline int setComputeSaddleConnectors(bool state){
        ComputeSaddleConnectors=state;
        return 0;
      }

      ftm::NodeType getNodeType(ftm::FTMTree_MT* tree,
                                ftm::TreeType    treeType,
                                const int        vertexId) const;

      template <typename scalarType>
      int sortPersistenceDiagram(vector<tuple<ftm::idVertex,
                                              ftm::NodeType,
                                              ftm::idVertex,
                                              ftm::NodeType,
                                              scalarType,
                                              ftm::idVertex>>& diagram,
                                 scalarType*                   scalars) const;

      template <typename scalarType>
      int computeCTPersistenceDiagram(
          ftm::FTMTreePP& tree,
          const vector<tuple<ftm::idVertex, ftm::idVertex, scalarType, bool>>& pairs,
          vector<tuple<ftm::idVertex,
                       ftm::NodeType,
                       ftm::idVertex,
                       ftm::NodeType,
                       scalarType,
                       ftm::idVertex>>& diagram,
          scalarType*                   scalars) const;

      template <class scalarType>
        int execute() const;

      inline int setDMTPairs(vector<tuple<Cell,Cell>>* data){
        dmt_pairs=data;
        return 0;
      }

      inline int setupTriangulation(Triangulation* data){
        triangulation_ = data;
        if(triangulation_){
          ftm::FTMTreePP contourTree;
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
    vector<tuple<ftm::idVertex,ftm::NodeType,ftm::idVertex,ftm::NodeType,scalarType,ftm::idVertex>>& diagram,
    scalarType* scalars) const{
  auto cmp=[scalars](const tuple<ftm::idVertex,ftm::NodeType,ftm::idVertex,ftm::NodeType,scalarType,ftm::idVertex>& a,
      const tuple<ftm::idVertex,ftm::NodeType,ftm::idVertex,ftm::NodeType,scalarType,ftm::idVertex>& b){
    return scalars[get<0>(a)] < scalars[get<0>(b)];
  };

  std::sort(diagram.begin(), diagram.end(), cmp);

  return 0;
}

template <typename scalarType>
int PersistenceDiagram::computeCTPersistenceDiagram(
    ftm::FTMTreePP& tree,
    const vector<tuple<ftm::idVertex, ftm::idVertex, scalarType, bool>>& pairs,
    vector<tuple<ftm::idVertex,
                 ftm::NodeType,
                 ftm::idVertex,
                 ftm::NodeType,
                 scalarType,
                 ftm::idVertex>>& diagram,
    scalarType*                   scalars) const
{
   const ftm::idVertex numberOfPairs = pairs.size();
   diagram.resize(numberOfPairs);
   for (ftm::idVertex i = 0; i < numberOfPairs; ++i) {
      const ftm::idVertex v0               = get<0>(pairs[i]);
      const ftm::idVertex v1               = get<1>(pairs[i]);
      const scalarType    persistenceValue = get<2>(pairs[i]);
      const bool          type             = get<3>(pairs[i]);

      get<4>(diagram[i]) = persistenceValue;
      if (type == true) {
         get<0>(diagram[i]) = v0;
         get<1>(diagram[i]) = getNodeType(tree.getJoinTree(), ftm::TreeType::Join, v0);
         get<2>(diagram[i]) = v1;
         get<3>(diagram[i]) = getNodeType(tree.getJoinTree(), ftm::TreeType::Join, v1);
         get<5>(diagram[i]) = 0;
      } else {
         get<0>(diagram[i]) = v1;
         get<1>(diagram[i]) = getNodeType(tree.getSplitTree(), ftm::TreeType::Split, v1);
         get<2>(diagram[i]) = v0;
         get<3>(diagram[i]) = getNodeType(tree.getSplitTree(), ftm::TreeType::Split, v0);
         get<5>(diagram[i]) = 2;
      }
   }

   return 0;
}

template <typename scalarType>
int PersistenceDiagram::execute() const{
  // get data
  vector<tuple<ftm::idVertex,ftm::NodeType,ftm::idVertex,ftm::NodeType,scalarType,ftm::idVertex>>& CTDiagram=
    *static_cast<vector<tuple<ftm::idVertex,ftm::NodeType,ftm::idVertex,ftm::NodeType,scalarType,ftm::idVertex>>*>(CTDiagram_);
  scalarType* scalars=static_cast<scalarType*>(inputScalars_);
  int* offsets=static_cast<int*>(inputOffsets_);

  const ftm::idVertex numberOfVertices=triangulation_->getNumberOfVertices();
  // convert offsets into a valid format for contour forests
  vector<ftm::idVertex> voffsets(numberOfVertices);
  std::copy(offsets,offsets+numberOfVertices,voffsets.begin());

  // get contour tree
  ftm::FTMTreePP contourTree;
  contourTree.setupTriangulation(triangulation_, false);
  contourTree.setVertexScalars(inputScalars_);
  contourTree.setTreeType(ftm::TreeType::Join_Split);
  contourTree.setVertexSoSoffsets(voffsets);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.setSegmentation(false);
  contourTree.build<scalarType>();

  // get persistence pairs
  vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>> JTPairs;
  vector<tuple<ftm::idVertex,ftm::idVertex,scalarType>> STPairs;
  contourTree.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree.computePersistencePairs<scalarType>(STPairs, false);

  // merge pairs
  vector<tuple<ftm::idVertex, ftm::idVertex, scalarType, bool>> CTPairs(JTPairs.size() +
                                                                        STPairs.size());
  const ftm::idVertex JTSize = JTPairs.size();
  for (ftm::idVertex i = 0; i < JTSize; ++i) {
     const auto& x = JTPairs[i];
     CTPairs[i]    = make_tuple(get<0>(x), get<1>(x), get<2>(x), true);
  }
  const ftm::idVertex STSize = STPairs.size();
  for (ftm::idVertex i = 0; i < STSize; ++i) {
     const auto& x       = STPairs[i];
     CTPairs[JTSize + i] = make_tuple(get<0>(x), get<1>(x), get<2>(x), false);
  }

  // remove the last pair which is present two times (global extrema pair)
  {
     auto cmp = [](const tuple<ftm::idVertex, ftm::idVertex, scalarType, bool>& a,
                   const tuple<ftm::idVertex, ftm::idVertex, scalarType, bool>& b) {
        return get<2>(a) < get<2>(b);
     };

     std::sort(CTPairs.begin(), CTPairs.end(), cmp);
     CTPairs.erase(CTPairs.end() - 1);
  }

  // get the saddle-saddle pairs
  vector<tuple<int,int,scalarType>> pl_saddleSaddlePairs;
  const int dimensionality=triangulation_->getDimensionality();
  if(dimensionality==3 and ComputeSaddleConnectors){
    MorseSmaleComplex3D morseSmaleComplex;
    morseSmaleComplex.setDebugLevel(debugLevel_);
    morseSmaleComplex.setThreadNumber(threadNumber_);
    morseSmaleComplex.setupTriangulation(triangulation_);
    morseSmaleComplex.setInputScalarField(inputScalars_);
    morseSmaleComplex.setInputOffsets(inputOffsets_);
    morseSmaleComplex.computePersistencePairs<scalarType>(JTPairs, STPairs, pl_saddleSaddlePairs);
  }

  // get persistence diagrams
  computeCTPersistenceDiagram<scalarType>(
      contourTree,CTPairs, CTDiagram, scalars);

  // add saddle-saddle pairs to the diagram if needed
  if(dimensionality==3 and ComputeSaddleConnectors){
    for(const auto& i : pl_saddleSaddlePairs){
      const ftm::idVertex v0=get<0>(i);
      const ftm::idVertex v1=get<1>(i);
      const scalarType persistenceValue=get<2>(i);

      tuple<ftm::idVertex, ftm::NodeType, ftm::idVertex, ftm::NodeType, scalarType, ftm::idVertex> t;

      get<0>(t)=v0;
      get<1>(t)=ftm::NodeType::Saddle1;
      get<2>(t)=v1;
      get<3>(t)=ftm::NodeType::Saddle2;
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
