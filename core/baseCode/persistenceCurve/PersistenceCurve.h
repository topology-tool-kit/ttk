/// \ingroup baseCode
/// \class ttk::PersistenceCurve
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2016.
///
/// \brief TTK processing package for the computation of persistence curves.
///
/// This package computes the list of extremum-saddle pairs and computes the 
/// number of pairs as a function of persistence (i.e. the number of pairs 
/// whose persistence is higher than a threshold).
///
/// These curves provide useful visual clues in order to fine-tune persistence
/// simplification thresholds.
///
/// \sa vtkPersistenceCurve.cpp %for a usage example.

#ifndef _PERSISTENCECURVE_H
#define _PERSISTENCECURVE_H

// base code includes
#include<Wrapper.h>
#include<Triangulation.h>
#include<FTMTreePP.h>
#include<MorseSmaleComplex3D.h>

namespace ttk{

  class PersistenceCurve : public Debug{

    public:

      PersistenceCurve();
      ~PersistenceCurve();

      inline int setComputeSaddleConnectors(bool state){
        ComputeSaddleConnectors=state;
        return 0;
      }

      template <typename scalarType>
        int computePersistencePlot(const vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>>& pairs,
            vector<pair<scalarType, ftm::idVertex>> &plot) const;

      template <class scalarType>
        int execute() const;

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

      inline int setOutputJTPlot(void* data){
        JTPlot_=data;
        return 0;
      }

      inline int setOutputSTPlot(void* data){
        STPlot_=data;
        return 0;
      }

      inline int setOutputMSCPlot(void* data){
        MSCPlot_=data;
        return 0;
      }

      inline int setOutputCTPlot(void* data){
        CTPlot_=data;
        return 0;
      }

    protected:

      bool ComputeSaddleConnectors;

      Triangulation* triangulation_;
      void* inputScalars_;
      void* inputOffsets_;
      void* JTPlot_;
      void* MSCPlot_;
      void* STPlot_;
      void* CTPlot_;
  };
}

template <typename scalarType>
int PersistenceCurve::computePersistencePlot(const vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>>& pairs,
    vector<pair<scalarType, ftm::idVertex>> &plot) const{

  ftm::idVertex nbElmnt = pairs.size();
  plot.resize(nbElmnt);

  // build curve
  const scalarType epsilon=static_cast<scalarType>(pow(10, -REAL_SIGNIFICANT_DIGITS));
  for (ftm::idVertex i = 0; i < nbElmnt; ++i) {
    plot[i].first  = std::max(get<2>(pairs[i]), epsilon);
    plot[i].second = pairs.size() - i;
  }

  return 0;
}

template <typename scalarType>
int PersistenceCurve::execute() const{
  // get data
  vector<pair<scalarType, ftm::idVertex>>& JTPlot  = *static_cast<vector<pair<scalarType, ftm::idVertex>>*>(JTPlot_);
  vector<pair<scalarType, ftm::idVertex>>& STPlot  = *static_cast<vector<pair<scalarType, ftm::idVertex>>*>(STPlot_);
  vector<pair<scalarType, ftm::idVertex>>& MSCPlot = *static_cast<vector<pair<scalarType, ftm::idVertex>>*>(MSCPlot_);
  vector<pair<scalarType, ftm::idVertex>>& CTPlot  = *static_cast<vector<pair<scalarType, ftm::idVertex>>*>(CTPlot_);
  int* offsets                                     = static_cast<int*>(inputOffsets_);

  const ftm::idVertex numberOfVertices=triangulation_->getNumberOfVertices();
  // convert offsets into a valid format for contour tree
  vector<ftm::idVertex> voffsets(numberOfVertices);
  std::copy(offsets,offsets+numberOfVertices,voffsets.begin());

  // get contour tree
  ftm::FTMTreePP contourTree;
  contourTree.setupTriangulation(triangulation_, false);
  contourTree.setVertexScalars(inputScalars_);
  contourTree.setTreeType(ftm::TreeType::Join_Split);
  contourTree.setVertexSoSoffsets(voffsets.data());
  contourTree.setSegmentation(false);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.build<scalarType>();

  // get persistence pairs
  vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>> JTPairs;
  vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>> STPairs;
  contourTree.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree.computePersistencePairs<scalarType>(STPairs, false);

  // merge pairs
  vector<tuple<ftm::idVertex, ftm::idVertex, scalarType>> CTPairs(JTPairs.size() + STPairs.size());
  std::copy(JTPairs.begin(), JTPairs.end(), CTPairs.begin());
  std::copy(STPairs.begin(), STPairs.end(), CTPairs.begin() + JTPairs.size());
  {
    auto cmp=[](const tuple<ftm::idVertex,ftm::idVertex,scalarType>& a,
        const tuple<ftm::idVertex,ftm::idVertex,scalarType>& b){
      return get<2>(a) < get<2>(b);
    };
    std::sort(CTPairs.begin(), CTPairs.end(), cmp);
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

    // sort the saddle-saddle pairs by persistence value and compute curve
    {
       auto cmp = [](const tuple<ftm::idVertex, ftm::idVertex, scalarType>& a,
                     const tuple<ftm::idVertex, ftm::idVertex, scalarType>& b) {
          return get<2>(a) < get<2>(b);
       };
       std::sort(pl_saddleSaddlePairs.begin(), pl_saddleSaddlePairs.end(), cmp);
    }

    computePersistencePlot<scalarType>(pl_saddleSaddlePairs, MSCPlot);
  }

  // get persistence curves
#ifdef withOpenMP
#pragma omp parallel sections
#endif
  {
#ifdef withOpenMP
#pragma omp section
#endif
    if(JTPlot_)
      computePersistencePlot<scalarType>(JTPairs, JTPlot);
#ifdef withOpenMP
#pragma omp section
#endif
    if(STPlot_)
      computePersistencePlot<scalarType>(STPairs, STPlot);
#ifdef withOpenMP
#pragma omp section
#endif
    if(CTPlot_)
      computePersistencePlot<scalarType>(CTPairs, CTPlot);
  }

  return 0;
}

#endif // PERSISTENCECURVE_H
