#pragma once

#include <MergeTreePrincipalGeodesics.h>

namespace ttk {
  template <class dataType>
  dataType MergeTreePrincipalGeodesics::computeVarianceFromDistances(
    std::vector<dataType> &distances) {
    dataType variance = 0.0;
    for(auto distance : distances)
      variance += distance * distance;
    return variance / distances.size();
  }

  template <class dataType>
  double MergeTreePrincipalGeodesics::computeExplainedVariance(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<ftm::MergeTree<dataType>> &trees,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v2,
    std::vector<double> &ts,
    bool globalVariance) {

    std::vector<ftm::MergeTree<dataType>> allInterpolated(trees.size());
    ftm::MergeTree<dataType> barycenterInterpolated;
    if(not globalVariance) {
      for(unsigned int i = 0; i < trees.size(); ++i)
        getInterpolation<dataType>(
          barycenter, v, v2, ts[i], allInterpolated[i]);
      computeOneBarycenter<dataType>(allInterpolated, barycenterInterpolated);
    }

    std::vector<dataType> distances(trees.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees.size(); ++i) {
      if(globalVariance) {
        computeOneDistance(barycenter, trees[i], distances[i], true);
      } else {
        computeOneDistance(
          barycenterInterpolated, allInterpolated[i], distances[i], true);
      }
    }
    return computeVarianceFromDistances(distances);
  }

  template <class dataType>
  double MergeTreePrincipalGeodesics::computeGlobalVariance(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<ftm::MergeTree<dataType>> &trees) {
    std::vector<double> ts;
    std::vector<std::vector<double>> v, v2;
    return computeExplainedVariance<dataType>(
      barycenter, trees, v, v2, ts, true);
  }

  template <class dataType>
  double MergeTreePrincipalGeodesics::computeSurfaceExplainedVariance(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<ftm::MergeTree<dataType>> &trees,
    std::vector<std::vector<std::vector<double>>> &vS,
    std::vector<std::vector<std::vector<double>>> &v2s,
    std::vector<std::vector<double>> &ts) {
    std::vector<ftm::MergeTree<dataType>> allInterpolated(trees.size());
    ftm::MergeTree<dataType> barycenterInterpolated;
    for(unsigned int i = 0; i < trees.size(); ++i) {
      std::vector<double> treeTs(ts.size());
      for(unsigned int j = 0; j < ts.size(); ++j)
        treeTs[j] = ts[j][i];
      getMultiInterpolation<dataType>(
        barycenter, vS, v2s, treeTs, allInterpolated[i]);
    }
    computeOneBarycenter<dataType>(allInterpolated, barycenterInterpolated);

    std::vector<dataType> distances(trees.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) \
  num_threads(this->threadNumber_) if(parallelize_)
#endif
    for(unsigned int i = 0; i < trees.size(); ++i) {
      computeOneDistance(
        barycenterInterpolated, allInterpolated[i], distances[i], true);
    }
    return computeVarianceFromDistances(distances);
  }

  template <class dataType>
  void MergeTreePrincipalGeodesics::computeProjectionDistances(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v2,
    std::vector<double> &ts,
    std::vector<double> &distances,
    bool useDoubleInput,
    bool isFirstInput) {
    ftm::MergeTree<dataType> extremity1, extremity2;
    getInterpolation(barycenter, v, v2, 0.0, extremity1);
    getInterpolation(barycenter, v, v2, 1.0, extremity2);
    dataType distance;
    computeOneDistance<dataType>(
      extremity1, extremity2, distance, false, useDoubleInput, isFirstInput);

    double tBarycenter = 0.0;
    for(auto t : ts)
      tBarycenter += t / ts.size();

    distances.resize(ts.size());
    for(unsigned int i = 0; i < ts.size(); ++i)
      distances[i] = (ts[i] - tBarycenter) * distance;
  }

  template <class dataType>
  double MergeTreePrincipalGeodesics::computeExplainedVarianceT(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v2,
    std::vector<double> &ts) {
    std::vector<double> distances;
    computeProjectionDistances(barycenter, v, v2, ts, distances);
    return computeVarianceFromDistances(distances);
  }

  template <class dataType>
  double MergeTreePrincipalGeodesics::computeExplainedVarianceT(
    ftm::MergeTree<dataType> &barycenter,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v2,
    ftm::MergeTree<dataType> &barycenter2,
    std::vector<std::vector<double>> &trees2V,
    std::vector<std::vector<double>> &trees2V2,
    std::vector<double> &ts) {
    std::vector<double> distances, distances2;
    computeProjectionDistances(barycenter, v, v2, ts, distances, true);
    computeProjectionDistances(
      barycenter2, trees2V, trees2V2, ts, distances2, true, false);
    for(unsigned int i = 0; i < distances.size(); ++i)
      distances[i] = mixDistances(distances[i], distances2[i]);
    return computeVarianceFromDistances(distances);
  }
} // namespace ttk
