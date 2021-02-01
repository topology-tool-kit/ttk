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

#pragma once

// base code includes
#include <DiscreteGradient.h>
#include <FTMTreePP.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * Compute the persistence curve of a function on a triangulation.
   * TTK assumes that the input dataset is made of only one connected component.
   */
  class PersistenceCurve : virtual public Debug {

  public:
    PersistenceCurve();

    inline void setComputeSaddleConnectors(bool state) {
      ComputeSaddleConnectors = state;
    }

    template <typename scalarType>
    int computePersistencePlot(
      const std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
      std::vector<std::pair<scalarType, SimplexId>> &plot) const;

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p inputOffsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    template <typename scalarType,
              class triangulationType = ttk::AbstractTriangulation>
    int execute(const scalarType *inputScalars,
                const SimplexId *inputOffsets,
                const triangulationType *triangulation);

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation) {
        triangulation->preconditionBoundaryVertices();
        contourTree_.setDebugLevel(debugLevel_);
        contourTree_.preconditionTriangulation(triangulation);
        if(this->ComputeSaddleConnectors) {
          dcg_.setDebugLevel(debugLevel_);
          dcg_.setThreadNumber(threadNumber_);
          dcg_.preconditionTriangulation(triangulation);
        }
      }
    }

    inline void setOutputJTPlot(void *const data) {
      JTPlot_ = data;
    }
    inline void setOutputSTPlot(void *const data) {
      STPlot_ = data;
    }
    inline void setOutputCTPlot(void *const data) {
      CTPlot_ = data;
    }
    inline void setOutputMSCPlot(void *const data) {
      MSCPlot_ = data;
    }

  protected:
    void *JTPlot_{};
    void *STPlot_{};
    void *CTPlot_{};
    void *MSCPlot_{};
    bool ComputeSaddleConnectors{false};
    ftm::FTMTreePP contourTree_{};
    dcg::DiscreteGradient dcg_{};
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
    = static_cast<scalarType>(Geometry::powIntTen(-REAL_SIGNIFICANT_DIGITS));
  for(SimplexId i = 0; i < nbElmnt; ++i) {
    plot[i].first = std::max(std::get<2>(pairs[i]), epsilon);
    plot[i].second = pairs.size() - i;
  }

  return 0;
}

template <typename scalarType, class triangulationType>
int ttk::PersistenceCurve::execute(const scalarType *inputScalars,
                                   const SimplexId *inputOffsets,
                                   const triangulationType *triangulation) {

  printMsg(ttk::debug::Separator::L1);

  Timer timer;

  using plotType = std::vector<std::pair<scalarType, SimplexId>>;

  auto JTPlot = static_cast<plotType *>(JTPlot_);
  auto STPlot = static_cast<plotType *>(STPlot_);
  auto CTPlot = static_cast<plotType *>(CTPlot_);
  auto MSCPlot = static_cast<plotType *>(MSCPlot_);

  contourTree_.setVertexScalars(inputScalars);
  contourTree_.setTreeType(ftm::TreeType::Join_Split);
  contourTree_.setVertexSoSoffsets(inputOffsets);
  contourTree_.setSegmentation(false);
  contourTree_.setThreadNumber(threadNumber_);
  contourTree_.build<scalarType>(triangulation);

  // get persistence pairs
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> JTPairs;
  std::vector<std::tuple<SimplexId, SimplexId, scalarType>> STPairs;
  contourTree_.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree_.computePersistencePairs<scalarType>(STPairs, false);

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
  const int dimensionality = triangulation->getDimensionality();
  if(dimensionality == 3 and ComputeSaddleConnectors and MSCPlot_ != nullptr) {
    std::vector<std::tuple<SimplexId, SimplexId, scalarType>>
      pl_saddleSaddlePairs;
    dcg_.setInputScalarField(inputScalars);
    dcg_.setInputOffsets(inputOffsets);
    dcg_.computeSaddleSaddlePersistencePairs<scalarType>(
      pl_saddleSaddlePairs, *triangulation);

    // sort the saddle-saddle pairs by persistence value and compute curve
    {
      auto cmp = [](const std::tuple<SimplexId, SimplexId, scalarType> &a,
                    const std::tuple<SimplexId, SimplexId, scalarType> &b) {
        return std::get<2>(a) < std::get<2>(b);
      };
      std::sort(pl_saddleSaddlePairs.begin(), pl_saddleSaddlePairs.end(), cmp);
    }

    computePersistencePlot<scalarType>(pl_saddleSaddlePairs, *MSCPlot);
  }

  // get persistence curves
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(JTPlot_ != nullptr) {
      computePersistencePlot<scalarType>(JTPairs, *JTPlot);
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(STPlot_ != nullptr) {
      computePersistencePlot<scalarType>(STPairs, *STPlot);
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
    if(CTPlot_ != nullptr) {
      computePersistencePlot<scalarType>(CTPairs, *CTPlot);
    }
  }

  printMsg(
    "Base execution completed", 1, timer.getElapsedTime(), threadNumber_);
  printMsg(ttk::debug::Separator::L1);

  return 0;
}
