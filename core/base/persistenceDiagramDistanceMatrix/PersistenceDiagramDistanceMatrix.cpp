#include <PersistenceDiagramDistanceMatrix.h>

std::vector<int> PersistenceDiagramDistanceMatrix::execute(
  std::vector<std::vector<DiagramTuple>> &intermediateDiagrams,
  std::vector<std::vector<std::vector<MatchingTuple>>> &all_matchings) {

  Timer tm;
  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Clustering " << numberOfInputs_
        << " diagrams in " << n_clusters_ << " cluster(s)." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  std::vector<std::vector<DiagramTuple>> data_min(numberOfInputs_);
  std::vector<std::vector<DiagramTuple>> data_sad(numberOfInputs_);
  std::vector<std::vector<DiagramTuple>> data_max(numberOfInputs_);

  std::vector<std::vector<int>> data_min_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_sad_idx(numberOfInputs_);
  std::vector<std::vector<int>> data_max_idx(numberOfInputs_);

  std::vector<int> inv_clustering(numberOfInputs_);

  bool do_min = false;
  bool do_sad = false;
  bool do_max = false;

  // Create diagrams for min, saddle and max persistence pairs
  for(int i = 0; i < numberOfInputs_; i++) {
    std::vector<DiagramTuple> &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      DiagramTuple t = CTDiagram[j];

      BNodeType nt1 = std::get<1>(t);
      BNodeType nt2 = std::get<3>(t);

      double dt = std::get<4>(t);
      // if (abs<double>(dt) < zeroThresh) continue;
      if(dt > 0) {
        if(nt1 == BLocalMin && nt2 == BLocalMax) {
          data_max[i].push_back(t);
          data_max_idx[i].push_back(j);
          do_max = true;
        } else {
          if(nt1 == BLocalMax || nt2 == BLocalMax) {
            data_max[i].push_back(t);
            data_max_idx[i].push_back(j);
            do_max = true;
          }
          if(nt1 == BLocalMin || nt2 == BLocalMin) {
            data_min[i].push_back(t);
            data_min_idx[i].push_back(j);
            do_min = true;
          }
          if((nt1 == BSaddle1 && nt2 == BSaddle2)
             || (nt1 == BSaddle2 && nt2 == BSaddle1)) {
            data_sad[i].push_back(t);
            data_sad_idx[i].push_back(j);
            do_sad = true;
          }
        }
      }
    }
  }

  {
    std::stringstream msg;
    switch(pairTypeClustering_) {
      case(0):
        msg << "[PersistenceDiagramDistanceMatrix] Only MIN-SAD Pairs";
        do_max = false;
        do_sad = false;
        break;
      case(1):
        msg << "[PersistenceDiagramDistanceMatrix] Only SAD-SAD Pairs";
        do_max = false;
        do_min = false;
        break;
      case(2):
        msg << "[PersistenceDiagramDistanceMatrix] Only SAD-MAX Pairs";
        do_min = false;
        do_sad = false;
        break;
      default:
        msg << "[PersistenceDiagramDistanceMatrix] All critical pairs: "
               "global clustering";
        break;
    }
    msg << std::endl;
    dMsg(std::cout, msg.str(), advancedInfoMsg);
  }

  vector<vector<vector<vector<MatchingTuple>>>>
    all_matchings_per_type_and_cluster;
  PDDistMat KMeans{};
  KMeans.setWasserstein(wasserstein_);
  KMeans.setThreadNumber(threadNumber_);
  KMeans.setNumberOfInputs(numberOfInputs_);
  KMeans.setUseProgressive(use_progressive_);
  KMeans.setAccelerated(use_accelerated_);
  KMeans.setUseKDTree(true);
  KMeans.setTimeLimit(time_limit_);
  KMeans.setGeometricalFactor(alpha_);
  KMeans.setLambda(lambda_);
  KMeans.setDeterministic(deterministic_);
  KMeans.setForceUseOfAlgorithm(forceUseOfAlgorithm_);
  KMeans.setDebugLevel(debugLevel_);
  KMeans.setDeltaLim(deltaLim_);
  KMeans.setUseDeltaLim(useDeltaLim_);
  KMeans.setDistanceWritingOptions(distanceWritingOptions_);
  KMeans.setKMeanspp(use_kmeanspp_);
  KMeans.setK(n_clusters_);
  KMeans.setDiagrams(&data_min, &data_sad, &data_max);
  KMeans.setDos(do_min, do_sad, do_max);
  KMeans.setOutputDistanceMatrix(outputDistanceMatrix_);
  KMeans.setUseFullDiagrams(useFullDiagrams_);
  KMeans.setPerClusterDistanceMatrix(perClusterDistanceMatrix_);
  inv_clustering = KMeans.execute();
  vector<vector<int>> centroids_sizes = KMeans.get_centroids_sizes();

  centroidsDistMat_ = KMeans.getCentroidsDistanceMatrix();
  diagramsDistMat_ = KMeans.getDiagramsDistanceMatrix();
  distanceToCentroid_ = KMeans.getDistanceToCentroid();

  /// Reconstruct matchings
  //
  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(numberOfInputs_);

  for(int j = 0; j < numberOfInputs_; ++j) {
    unsigned int c = inv_clustering[j];
    if(c + 1 > cluster_size.size()) {
      cluster_size.resize(c + 1);
      cluster_size[c] = 1;
      idxInCluster[j] = 0;
    } else {
      cluster_size[c]++;
      idxInCluster[j] = cluster_size[c] - 1;
    }
    if(debugLevel_ > 20) {
      cout << "id in cluster " << idxInCluster[j] << endl;
    }
  }

  all_matchings.resize(n_clusters_);
  for(int c = 0; c < n_clusters_; c++) {
    all_matchings[c].resize(numberOfInputs_);
  }

  {
    std::stringstream msg;
    msg << "[PersistenceDiagramDistanceMatrix] Processed in "
        << tm.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))."
        << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }
  return inv_clustering;
}
