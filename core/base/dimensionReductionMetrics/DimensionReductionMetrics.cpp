#include <DimensionReductionMetrics.h>

#include <random>

#include <PersistenceDiagramWarmRestartAuction.h>
#include <ripser.h>

ttk::DimensionReductionMetrics::DimensionReductionMetrics() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("DimensionReductionMetrics");
}

void ttk::DimensionReductionMetrics::execute(
  std::vector<std::vector<double>> const &input,
  std::vector<std::vector<double>> const &latent) {
  Timer tm{};

  if(input.size() != latent.size()) {
    printErr("Input and representation have different sizes");
    return;
  }
  if(input.size() == 0) {
    printErr("Empty input / representation");
    return;
  }

  n_ = input.size();
  dimHigh_ = input[0].size();
  dimLow_ = latent[0].size();
  printMsg("#pts: " + std::to_string(n_) + ", highDim: "
           + std::to_string(dimHigh_) + ", lowDim: " + std::to_string(dimLow_));

  inputCompressedDistanceMatrix_.clear();
  latentCompressedDistanceMatrix_.clear();

  for(unsigned i = 1; i < n_; ++i) {
    for(unsigned j = 0; j < i; ++j) {
      double sHigh = 0., sLow = 0.;
      for(unsigned d = 0; d < dimHigh_; ++d)
        sHigh += (input[i][d] - input[j][d]) * (input[i][d] - input[j][d]);
      inputCompressedDistanceMatrix_.push_back(sqrt(sHigh));
      for(unsigned d = 0; d < dimLow_; ++d)
        sLow += (latent[i][d] - latent[j][d]) * (latent[i][d] - latent[j][d]);
      latentCompressedDistanceMatrix_.push_back(sqrt(sLow));
    }
  }

  printMsg("Computing Wasserstein metrics", 0, tm.getElapsedTime());
  computeTopologicalMetrics();
  printMsg("Computing triplet accuracy", 0, tm.getElapsedTime());
  computeTripletAccuracy();
  printMsg("Computing distance-based measures", 0, tm.getElapsedTime());
  computePairwiseDistanceBasedMetrics();
  printMsg("Computing rank-based measures", 0, tm.getElapsedTime());
  computeRankBasedMetrics();
  printMsg("Complete", 1, tm.getElapsedTime());
}

void ttk::DimensionReductionMetrics::computeTopologicalMetrics() {
  rpd::MultidimensionalDiagram inputPD, latentPD;
  ripser::ripser({inputCompressedDistanceMatrix_}, inputPD, rpd::inf, 1, true);
  ripser::ripser(
    {latentCompressedDistanceMatrix_}, latentPD, rpd::inf, 1, true);
  inputPD[0].pop_back();
  latentPD[0].pop_back();

  PersistenceDiagramWarmRestartAuction auction0(inputPD[0]);
  auction0.setNewBidder(latentPD[0]);
  auction0.setWasserstein(Wasserstein);
  w0_ = auction0.runAuction();

  PersistenceDiagramWarmRestartAuction auction1(inputPD[1]);
  auction1.setNewBidder(latentPD[1]);
  auction1.setWasserstein(Wasserstein);
  w1_ = auction1.runAuction();
}

bool ttk::DimensionReductionMetrics::tripletOrderPreserved(int i,
                                                           int j,
                                                           int k) const {
  std::vector<std::pair<double, char>> inputTriangle
    = {{inputDM(i, j), 0}, {inputDM(i, k), 1}, {inputDM(j, k), 2}};
  std::vector<std::pair<double, char>> latentTriangle
    = {{latentDM(i, j), 0}, {latentDM(i, k), 1}, {latentDM(j, k), 2}};
  std::sort(inputTriangle.begin(), inputTriangle.end());
  std::sort(latentTriangle.begin(), latentTriangle.end());
  return inputTriangle[0].second == latentTriangle[0].second
         && inputTriangle[1].second == latentTriangle[1].second
         && inputTriangle[2].second == latentTriangle[2].second;
}

void ttk::DimensionReductionMetrics::computeTripletAccuracy() {
  unsigned stableTriplets = 0;
  if(SampleSize <= -1) {
    for(unsigned i = 0; i < n_; ++i) {
      for(unsigned j = i + 1; j < n_; ++j) {
        for(unsigned k = j + 1; k < n_; ++k) {
          if(tripletOrderPreserved(i, j, k))
            ++stableTriplets;
        }
      }
    }
    ta_ = double(stableTriplets) * 6 / (n_ * (n_ - 1) * (n_ - 2));
  } else {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> uniform(0, n_ - 1);
    for(int s = 0; s < SampleSize; ++s) {
      if(tripletOrderPreserved(uniform(gen), uniform(gen), uniform(gen)))
        ++stableTriplets;
    }
    ta_ = double(stableTriplets) / SampleSize;
  }
}

void ttk::DimensionReductionMetrics::computePairwiseDistanceBasedMetrics() {
  double sumX = 0., sumY = 0., sumXY = 0.;
  double sqSumX = 0., sqSumY = 0.;
  double sqSumDiff = 0.;
  const unsigned Ndis = n_ * (n_ - 1) / 2;
  for(unsigned i = 0; i < Ndis; ++i) {
    const double x = inputCompressedDistanceMatrix_[i];
    const double y = latentCompressedDistanceMatrix_[i];
    sumX += x;
    sqSumX += x * x;
    sumY += y;
    sqSumY += y * y;
    sumXY += x * y;
    sqSumDiff += (x - y) * (x - y);
  }
  lc_ = (Ndis * sumXY - sumX * sumY)
        / sqrt((Ndis * sqSumX - sumX * sumX) * (Ndis * sqSumY - sumY * sumY));
  rmse_ = sqrt(sqSumDiff / Ndis);
}

void ttk::DimensionReductionMetrics::computeRankBasedMetrics() {
  NeighborhoodSize = std::min(NeighborhoodSize, n_ - 1);

  int trustworthinessSum = 0;
  int continuitySum = 0;
  int LCMCSum = 0;
  double inputMRRESum = 0;
  double latentMRRESum = 0;
  const int normalizingTC
    = (NeighborhoodSize < n_ / 2)
        ? n_ * NeighborhoodSize * (2 * n_ - 3 * NeighborhoodSize - 1)
        : n_ * (n_ - NeighborhoodSize) * (n_ - NeighborhoodSize - 1);
  double normalizingMRRE = 0.;
  for(unsigned k = 1; k < NeighborhoodSize; ++k)
    normalizingMRRE += n_ * std::abs(double(n_) - 2 * k + 1) / k;

  for(unsigned i = 0; i < n_; ++i) {
    std::vector<std::pair<double, unsigned>> inputNeighborhood;
    std::vector<std::pair<double, unsigned>> latentNeighborhood;
    for(unsigned j = 0; j < n_; ++j) {
      inputNeighborhood.emplace_back(inputDM(i, j), j);
      latentNeighborhood.emplace_back(latentDM(i, j), j);
    }
    std::sort(inputNeighborhood.begin(), inputNeighborhood.end());
    std::sort(latentNeighborhood.begin(), latentNeighborhood.end());

    std::vector<unsigned> inputRanks(n_);
    std::vector<unsigned> latentRanks(n_);
    for(unsigned s = 0; s < n_; ++s) {
      inputRanks[inputNeighborhood[s].second] = s;
      latentRanks[latentNeighborhood[s].second] = s;
    }

    for(unsigned j = 0; j < n_; ++j) {
      if(latentRanks[j] <= NeighborhoodSize && inputRanks[j] > NeighborhoodSize)
        trustworthinessSum += inputRanks[j] - NeighborhoodSize;
      else if(inputRanks[j] <= NeighborhoodSize
              && latentRanks[j] > NeighborhoodSize)
        continuitySum += latentRanks[j] - NeighborhoodSize;
      else if(inputRanks[j] <= NeighborhoodSize
              && latentRanks[j] <= NeighborhoodSize && i != j)
        LCMCSum += 1;
    }

    for(unsigned s = 1; s <= NeighborhoodSize; ++s) {
      const unsigned inputNeighbor = inputNeighborhood[s].second;
      inputMRRESum += double(std::abs(int(s) - int(latentRanks[inputNeighbor])))
                      / latentRanks[inputNeighbor];
      const unsigned latentNeighbor = latentNeighborhood[s].second;
      latentMRRESum
        += double(std::abs(int(s) - int(inputRanks[latentNeighbor])))
           / inputRanks[latentNeighbor];
    }
  }

  trust_ = 1. - 2 * double(trustworthinessSum) / normalizingTC;
  cont_ = 1. - 2 * double(continuitySum) / normalizingTC;
  lcmc_ = (double(LCMCSum) / (n_ * NeighborhoodSize)
           - double(NeighborhoodSize) / (n_ - 1))
          / (1 - double(NeighborhoodSize) / (n_ - 1));
  mrreh_ = inputMRRESum / normalizingMRRE;
  mrrel_ = latentMRRESum / normalizingMRRE;
}
