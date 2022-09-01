#pragma once

namespace ttk {
  template <class dataType>
  void MergeTreePrincipalGeodesics::getTransportationMatrix(
    ftm::MergeTree<dataType> &tree1,
    ftm::MergeTree<dataType> &tree2,
    std::vector<int> &tree1Corr,
    std::vector<int> &tree2Corr,
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
    std::vector<std::vector<double>> &transportationMatrix) {
    unsigned int tree1CorrRealSize = tree1.tree.getRealNumberOfNodes(),
                 tree2CorrRealSize = tree2.tree.getRealNumberOfNodes();
    unsigned int totalSize = tree1CorrRealSize + tree2CorrRealSize;
    transportationMatrix = std::vector<std::vector<double>>(
      totalSize, std::vector<double>(totalSize, 0.0));
    std::vector<ftm::idNode> tree1MatchingVector, tree2MatchingVector;
    getMatchingVector(tree1, tree2, matching, tree1MatchingVector);
    getInverseMatchingVector(tree1, tree2, matching, tree2MatchingVector);
    for(unsigned int i = 0; i < tree1MatchingVector.size(); ++i) {
      if(tree1.tree.isNodeAlone(i))
        continue;
      int index = ((int)tree1MatchingVector[i] != -1
                     ? tree2Corr[tree1MatchingVector[i]]
                     : tree1Corr[i] + tree2CorrRealSize);
      transportationMatrix[tree1Corr[i]][index] = 1.0 / totalSize;
    }
    for(unsigned int i = 0; i < tree2MatchingVector.size(); ++i) {
      if(tree2.tree.isNodeAlone(i))
        continue;
      if((int)tree2MatchingVector[i] == -1) {
        int index = tree2Corr[i] + tree1CorrRealSize;
        transportationMatrix[index][tree2Corr[i]] = 1.0 / totalSize;
      }
    }
  }

  template <class dataType>
  void MergeTreePrincipalGeodesics::getTreeMatrix(
    ftm::MergeTree<dataType> &tree,
    std::vector<std::vector<double>> &treeMatrix,
    std::vector<int> &treeCorr) {
    auto noNodes = tree.tree.getNumberOfNodes();
    auto realNoNodes = tree.tree.getRealNumberOfNodes();
    treeMatrix = std::vector<std::vector<double>>(
      realNoNodes, std::vector<double>(2, 0));
    treeCorr = std::vector<int>(noNodes, -1);
    ftm::FTMTree_MT *ftmTree = &(tree.tree);
    int cpt = 0;
    for(unsigned int i = 0; i < noNodes; ++i) {
      if(tree.tree.isNodeAlone(i))
        continue;
      auto birthDeath = getParametrizedBirthDeath<dataType>(ftmTree, i);
      treeMatrix[cpt][0] = std::get<0>(birthDeath);
      treeMatrix[cpt][1] = std::get<1>(birthDeath);
      treeCorr[i] = cpt;
      ++cpt;
    }
  }

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
    bool computeGlobalVariance) {
    std::vector<ftm::MergeTree<dataType>> allInterpolated(trees.size());
    ftm::MergeTree<dataType> barycenterInterpolated;
    if(not computeGlobalVariance) {
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
      if(computeGlobalVariance) {
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

    distances = std::vector<double>(ts.size());
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

  template <class dataType>
  std::tuple<dataType, dataType>
    MergeTreePrincipalGeodesics::getParametrizedBirthDeath(
      ftm::FTMTree_MT *tree, ftm::idNode node) {
    return normalizedWasserstein_
             ? getNormalizedBirthDeath<dataType>(tree, node)
             : tree->getBirthDeath<dataType>(node);
  }
} // namespace ttk
