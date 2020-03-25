#include <PersistenceDiagramDistanceMatrix.h>

void PersistenceDiagramDistanceMatrix::execute(
  std::vector<std::vector<DiagramTuple>> &intermediateDiagrams) {

  PDDistMat KMeans{};
  KMeans.setWasserstein(wasserstein_);
  KMeans.setThreadNumber(threadNumber_);
  KMeans.setNumberOfInputs(numberOfInputs_);
  KMeans.setUseProgressive(use_progressive_);
  KMeans.setAccelerated(use_accelerated_);
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
  KMeans.setOutputDistanceMatrix(outputDistanceMatrix_);
  KMeans.setUseFullDiagrams(useFullDiagrams_);
  KMeans.setPairTypeClustering(pairTypeClustering_);
  KMeans.execute(intermediateDiagrams);

  diagramsDistMat_ = KMeans.getDiagramsDistanceMatrix();
}
