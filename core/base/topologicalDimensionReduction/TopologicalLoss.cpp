#include "TopologicalLoss.h"

#ifdef TTK_ENABLE_TORCH

using namespace torch::indexing;

ttk::TopologicalLoss::TopologicalLoss(
  const torch::Tensor &input,
  const std::vector<std::vector<double>> &points,
  REGUL regul)
  : input_(input), points_(points), regul_(regul),
    device(input.device().type()) {
  precomputeInputPersistence();
}

torch::Tensor ttk::TopologicalLoss::computeLoss(const torch::Tensor &latent) {
  latent_ = latent;
  latentDimension = latent.size(1);

  if(regul_ == REGUL::TOPOAE)
    return diffTopoAELoss();
  else if(regul_ == REGUL::TOPOAE_DIM1)
    return diffTopoAELossDim1();
  else if(regul_ == REGUL::CASCADE)
    return diffCascadeAELoss();
  else if(regul_ == REGUL::ASYMMETRIC_CASCADE)
    return diffAsymmetricCascadeAELoss();
  else if(regul_ == REGUL::W_DIM1)
    return diffW1() + diffTopoAELoss();

  return torch::zeros({1}, device);
}

void ttk::TopologicalLoss::precomputeInputPersistence() {
  if(regul_ == REGUL::TOPOAE || regul_ == REGUL::W_DIM1) {
    rpd::EdgeSets3 inputCritical;
    ripser::ripser(points_, inputCritical, rpd::inf, 0, false);
    inputCriticalPairIndices = {pairsToTorch(inputCritical[0])};
  }
  if(regul_ == REGUL::TOPOAE_DIM1) {
    rpd::EdgeSets3 inputCritical;
    ripser::ripser(points_, inputCritical, rpd::inf, 1, false);
    for(int i = 0; i <= 2; ++i)
      inputCriticalPairIndices[i] = pairsToTorch(inputCritical[i]);
  } else if(regul_ == REGUL::W_DIM1) {
    ripser::ripser(points_, inputPD, rpd::inf, 1, false);
    auction = std::make_unique<
      PersistenceDiagramWarmRestartAuction<rpd::PersistencePair>>(inputPD[1]);
  } else if(regul_ == REGUL::CASCADE || regul_ == REGUL::ASYMMETRIC_CASCADE) {
    // first compute the PD with Ripser
    rpd::PairCellsWithOracle::callOracle(points_, inputPD);
    // use it to quickly compute the cascade
    rpd::PairCellsWithOracle pc(points_, inputPD, false, false);
    pc.run();

    rpd::EdgeSets4 inputCriticalAndCascade; // [0] is MST, [1] is RNG-MST, [2]
                                            // is MML, [3] is strict cascade
    pc.getCascades(inputCriticalAndCascade);
    for(int i = 0; i <= 3; ++i)
      inputCriticalPairIndices[i] = pairsToTorch(inputCriticalAndCascade[i]);
  }
}

void ttk::TopologicalLoss::computeLatent0Persistence(
  rpd::EdgeSet &latent0PD) const {
#ifdef TTK_ENABLE_CGAL
  if(latentDimension == 2)
    rpd::FastRipsPersistenceDiagram2(
      latent_.cpu().data_ptr<float>(), latent_.size(0))
      .compute0Persistence(latent0PD);
  else {
    rpd::EdgeSets3 latentPD;
    ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                   latent_.size(1), latentPD, rpd::inf, 0, false);
    latent0PD = latentPD[0];
  }
#else
  rpd::EdgeSets3 latentPD;
  ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                 latent_.size(1), latentPD, rpd::inf, 0, false);
  latent0PD = latentPD[0];
#endif
}

template <typename PersistenceType>
void ttk::TopologicalLoss::computeLatent0And1Persistence(
  PersistenceType &latentPD) const {
#ifdef TTK_ENABLE_CGAL
  if(latentDimension == 2)
    rpd::FastRipsPersistenceDiagram2(
      latent_.cpu().data_ptr<float>(), latent_.size(0))
      .computeRips0And1Persistence(latentPD, false, false);
  else
    ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                   latent_.size(1), latentPD, rpd::inf, 1, false);
#else
  ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                 latent_.size(1), latentPD, rpd::inf, 1, false);
#endif
}
template void ttk::TopologicalLoss::computeLatent0And1Persistence(
  rpd::MultidimensionalDiagram &latentPD) const;
template void ttk::TopologicalLoss::computeLatent0And1Persistence(
  rpd::EdgeSets3 &latentPD) const;

void ttk::TopologicalLoss::computeLatentCascades(
  rpd::EdgeSets4 &latentCriticalAndCascades) const {
#ifdef TTK_ENABLE_CGAL
  if(latentDimension == 2)
    rpd::FastRipsPersistenceDiagram2(
      latent_.cpu().data_ptr<float>(), latent_.size(0))
      .computeRips0And1Persistence(latentCriticalAndCascades, false, false);
  else {
    rpd::MultidimensionalDiagram latentPD;
    ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                   latent_.size(1), latentPD, rpd::inf, 1, false);
    rpd::PairCellsWithOracle pc(latent_.cpu().data_ptr<float>(),
                                latent_.size(0), latent_.size(1), latentPD,
                                false);
    pc.run();
    pc.getCascades(latentCriticalAndCascades);
  }
#else
  rpd::MultidimensionalDiagram latentPD;
  ripser::ripser(latent_.cpu().data_ptr<float>(), latent_.size(0),
                 latent_.size(1), latentPD, rpd::inf, 1, false);
  rpd::PairCellsWithOracle pc(latent_.cpu().data_ptr<float>(), latent_.size(0),
                              latent_.size(1), latentPD, false);
  pc.run();
  pc.getCascades(latentCriticalAndCascades);
#endif
}

/*** TopoAE-like distances ***/

torch::Tensor ttk::TopologicalLoss::diffTopoAELoss() const {
  rpd::EdgeSet latent0Critical;
  computeLatent0Persistence(latent0Critical);

  return diffEdgeSetMSE(inputCriticalPairIndices[0])
         + diffEdgeSetMSE(pairsToTorch(latent0Critical));
}

torch::Tensor ttk::TopologicalLoss::diffTopoAELossDim1() const {
  rpd::EdgeSets3 latentCritical;
  computeLatent0And1Persistence(latentCritical);
  return diffRNGMML(latentCritical);
}

torch::Tensor ttk::TopologicalLoss::diffCascadeAELoss() const {
  rpd::EdgeSets4 latentCriticalAndCascade;
  computeLatentCascades(latentCriticalAndCascade);

  return diffRNGMML(latentCriticalAndCascade)
         + diffEdgeSetMSE(inputCriticalPairIndices[3])
         + diffEdgeSetMSE(pairsToTorch(latentCriticalAndCascade[rpd::CASC1]));
}

torch::Tensor ttk::TopologicalLoss::diffAsymmetricCascadeAELoss() const {
  rpd::EdgeSets3 latentCritical;
  computeLatent0And1Persistence(latentCritical);

  return diffRNGMML(latentCritical)
         + diffEdgeSetMSE(inputCriticalPairIndices[3]);
}

/*** Wasserstein distances ***/

void ttk::TopologicalLoss::performAuction(
  const rpd::Diagram &latentPD,
  std::vector<unsigned> &directMatchingLatent,
  std::vector<unsigned> &directMatchingInput,
  std::vector<unsigned> &diagonalMatchingLatent) const {
  auction->setNewBidder(latentPD);
  std::vector<MatchingType> matchings;
  auction->runAuction(matchings);

  for(const MatchingType &m : matchings) {
    if(std::get<0>(m) < 0)
      break;
    else {
      if(std::get<1>(m) >= 0) {
        directMatchingLatent.push_back(std::get<0>(m));
        directMatchingInput.push_back(std::get<1>(m));
      } else
        diagonalMatchingLatent.push_back(std::get<0>(m));
    }
  }
}

torch::Tensor ttk::TopologicalLoss::diffW1() const {
  std::vector<rpd::Diagram> latentPD(0);
  computeLatent0And1Persistence(latentPD);

  if(latentPD[1].empty())
    return torch::zeros(1, device);

  std::vector<unsigned> directMatchingLatent(0);
  std::vector<unsigned> directMatchingInput(0);
  std::vector<unsigned> diagonalMatchingLatent(0);
  performAuction(latentPD[1], directMatchingLatent, directMatchingInput,
                 diagonalMatchingLatent);

  torch::Tensor directMatchedInputPD
    = torch::zeros({int(directMatchingInput.size()), 2});
  float *inputData = directMatchedInputPD.data_ptr<float>();
  for(unsigned i = 0; i < directMatchingInput.size(); ++i) {
    inputData[2 * i]
      = float(inputPD[1][directMatchingInput[i]].first.second); // birth
    inputData[2 * i + 1]
      = float(inputPD[1][directMatchingInput[i]].second.second); // death
  }
  directMatchedInputPD = directMatchedInputPD.to(device);

  const torch::Tensor directMatchedLatentPD
    = diffPD(latent_, latentPD[1], directMatchingLatent);
  const torch::Tensor diagonalMatchedLatentPD
    = diffPD(latent_, latentPD[1], diagonalMatchingLatent);
  const torch::Tensor diagProj
    = sqrt(2) / 2
      * (diagonalMatchedLatentPD.index({Slice(), 1})
         - diagonalMatchedLatentPD.index({Slice(), 0}));

  const torch::Tensor costs = torch::cat(
    {(directMatchedInputPD - directMatchedLatentPD).norm(2, 1), diagProj});
  return costs.norm(2);
}

/*** tensor tools ***/

torch::Tensor
  ttk::TopologicalLoss::pairsToTorch(const rpd::EdgeSet &edges) const {
  const torch::Tensor pairIndices
    = torch::zeros({(int)edges.size(), 2}, torch::kInt);
  int *data = pairIndices.data_ptr<int>();
  for(unsigned i = 0; i < edges.size(); ++i) {
    data[2 * i] = edges[i].first;
    data[2 * i + 1] = edges[i].second;
  }
  return pairIndices.to(device);
}

torch::Tensor
  ttk::TopologicalLoss::diffDistances(const torch::Tensor &data,
                                      const torch::Tensor &indices) {
  const torch::Tensor points = data.index({indices});
  return (points.index({Slice(), 0}) - points.index({Slice(), 1})).norm(2, 1);
}

torch::Tensor
  ttk::TopologicalLoss::diffEdgeSetMSE(const torch::Tensor &indices) const {
  if(indices.size(0) > 0) // compute only if the edge set is non-empty
    return torch::mse_loss(diffDistances(latent_, indices),
                           diffDistances(input_, indices), reduction_);
  else
    return torch::zeros({1}, device);
}

torch::Tensor
  ttk::TopologicalLoss::diffPD(const torch::Tensor &points,
                               const rpd::Diagram &PD,
                               const std::vector<unsigned int> &indices) const {
  torch::Tensor edgesIndices
    = torch::zeros({int(indices.size()), 2, 2}, torch::kInt);
  int *data = edgesIndices.data_ptr<int>();
  for(unsigned i = 0; i < indices.size(); ++i) {
    const unsigned index = indices[i];
    data[4 * i] = PD[index].first.first[0];
    data[4 * i + 1] = PD[index].first.first[1];
    data[4 * i + 2] = PD[index].second.first[0];
    data[4 * i + 3] = PD[index].second.first[1];
  }
  edgesIndices = edgesIndices.to(device);
  const torch::Tensor PDPoints = points.index({edgesIndices});
  return (PDPoints.index({Slice(), Slice(), 0})
          - PDPoints.index({Slice(), Slice(), 1}))
    .norm(2, 2);
}

#endif
