/// \ingroup base
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
/// \sa ttkPersistenceCurve.cpp %for a usage example.

#ifndef _PERSISTENCECURVE_H
#define _PERSISTENCECURVE_H

// base code includes
#include <FTMTreePP.h>
#include <MorseSmaleComplex3D.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk {

  /**
   * Compute the persistence curve of a function on a triangulation.
   * TTK assumes that the input dataset is made of only one connected component.
   */
  class PersistenceCurve : public Debug {

  public:
    PersistenceCurve();
    ~PersistenceCurve();

    inline int setComputeSaddleConnectors(bool state) {
      ComputeSaddleConnectors = state;
      return 0;
    }

    template <typename scalarType>
    int computePersistencePlot(
      const std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
      std::vector<std::pair<scalarType, SimplexId>> &plot) const;

    template <typename scalarType, typename idType>
    int execute() const;

    inline int setupTriangulation(Triangulation *data) {
      triangulation_ = data;
      if(triangulation_) {
        ftm::FTMTreePP contourTree;
        contourTree.setDebugLevel(debugLevel_);
        contourTree.setupTriangulation(triangulation_);

        triangulation_->preprocessBoundaryVertices();
      }
      return 0;
    }

    inline int setInputScalars(void *data) {
      inputScalars_ = data;
      return 0;
    }

    inline int setInputOffsets(void *data) {
      inputOffsets_ = data;
      return 0;
    }

    inline int setOutputJTPlot(void *data) {
      JTPlot_ = data;
      return 0;
    }

    inline int setOutputSTPlot(void *data) {
      STPlot_ = data;
      return 0;
    }

    inline int setOutputMSCPlot(void *data) {
      MSCPlot_ = data;
      return 0;
    }

    inline int setOutputCTPlot(void *data) {
      CTPlot_ = data;
      return 0;
    }

  protected:
    bool ComputeSaddleConnectors;

    Triangulation *triangulation_;
    void *inputScalars_;
    void *inputOffsets_;
    void *JTPlot_;
    void *MSCPlot_;
    void *STPlot_;
    void *CTPlot_;
  };
} // namespace ttk

template <typename scalarType>
int ttk::PersistenceCurve::computePersistencePlot(
  const std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
  std::vector<std::pair<scalarType, SimplexId>> &plot) const {

  SimplexId nbElmnt = pairs.size();
  plot.resize(nbElmnt);

  // build curve
  const scalarType epsilon
    = static_cast<scalarType>(pow(10, -REAL_SIGNIFICANT_DIGITS));
  for(SimplexId i = 0; i < nbElmnt; ++i) {
    plot[i].first = std::max(std::get<2>(pairs[i]), epsilon);
    plot[i].second = pairs.size() - i;
  }

  return 0;
}

template <typename scalarType, typename idType>
int ttk::PersistenceCurve::execute() const {
  // get data
  std::vector<std::pair<scalarType, SimplexId>> &JTPlot
    = *static_cast<std::vector<std::pair<scalarType, SimplexId>> *>(JTPlot_);
  std::vector<std::pair<scalarType, SimplexId>> &STPlot
    = *static_cast<std::vector<std::pair<scalarType, SimplexId>> *>(STPlot_);
  std::vector<std::pair<scalarType, SimplexId>> &MSCPlot
    = *static_cast<std::vector<std::pair<scalarType, SimplexId>> *>(MSCPlot_);
  std::vector<std::pair<scalarType, SimplexId>> &CTPlot
    = *static_cast<std::vector<std::pair<scalarType, SimplexId>> *>(CTPlot_);
  SimplexId *offsets = static_cast<SimplexId *>(inputOffsets_);

  const SimplexId numberOfVertices = triangulation_->getNumberOfVertices();
  // convert offsets into a valid format for contour tree
  std::vector<SimplexId> voffsets(numberOfVertices);
  std::copy(offsets, offsets + numberOfVertices, voffsets.begin());

  // get contour tree
  ftm::FTMTreePP contourTree;
  contourTree.setDebugLevel(debugLevel_);
  contourTree.setupTriangulation(triangulation_, false);
  contourTree.setVertexScalars(inputScalars_);
  contourTree.setTreeType(ftm::TreeType::Join_Split);
  contourTree.setVertexSoSoffsets(voffsets.data());
  contourTree.setSegmentation(false);
  contourTree.setThreadNumber(threadNumber_);
  contourTree.build<scalarType, idType>();

  // get persistence pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> JTPairs;
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> STPairs;
  contourTree.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree.computePersistencePairs<scalarType>(STPairs, false);

  // merge pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> CTPairs(
    JTPairs.size() + STPairs.size());
  std::copy(JTPairs.begin(), JTPairs.end(), CTPairs.begin());
  std::copy(STPairs.begin(), STPairs.end(), CTPairs.begin() + JTPairs.size());
  {
    auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType> &a,
                  const std::tuple<SimplexId, SimplexId, scalarType> &b) {
      return std::get<2>(a) < std::get<2>(b);
    };
    std::sort(CTPairs.begin(), CTPairs.end(), cmp);
  }

  // get the saddle-saddle pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>>
    pl_saddleSaddlePairs;
  const int dimensionality = triangulation_->getDimensionality();
  if(dimensionality == 3 and ComputeSaddleConnectors) {
    MorseSmaleComplex3D morseSmaleComplex;
    morseSmaleComplex.setDebugLevel(debugLevel_);
    morseSmaleComplex.setThreadNumber(threadNumber_);
    morseSmaleComplex.setupTriangulation(triangulation_);
    morseSmaleComplex.setInputScalarField(inputScalars_);
    morseSmaleComplex.setInputOffsets(inputOffsets_);
    morseSmaleComplex.computePersistencePairs<scalarType, idType>(
      JTPairs, STPairs, pl_saddleSaddlePairs);

    // sort the saddle-saddle pairs by persistence value and compute curve
    {
      auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType> &a,
                    const std::tuple<SimplexId, SimplexId, scalarType> &b) {
        return std::get<2>(a) < std::get<2>(b);
      };
      std::sort(pl_saddleSaddlePairs.begin(), pl_saddleSaddlePairs.end(), cmp);
    }

    computePersistencePlot<scalarType>(pl_saddleSaddlePairs, MSCPlot);
  }

  // get persistence curves
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(JTPlot_)
      computePersistencePlot<scalarType>(JTPairs, JTPlot);
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(STPlot_)
      computePersistencePlot<scalarType>(STPairs, STPlot);
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(CTPlot_)
      computePersistencePlot<scalarType>(CTPairs, CTPlot);
  }

  return 0;
}

#endif // PERSISTENCECURVE_H
