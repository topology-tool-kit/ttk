#include <MergeTreeAutoencoderDecoding.h>
#include <MergeTreeAutoencoderUtils.h>

ttk::MergeTreeAutoencoderDecoding::MergeTreeAutoencoderDecoding() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MergeTreeAutoencoderDecoding");
}

void ttk::MergeTreeAutoencoderDecoding::execute(
  std::vector<ttk::ftm::MergeTree<float>> &originsTrees,
  std::vector<ttk::ftm::MergeTree<float>> &originsPrimeTrees,
  std::vector<unsigned int *> &allRevNodeCorr,
  std::vector<unsigned int *> &allRevNodeCorrPrime,
  std::vector<unsigned int> &allRevNodeCorrSize,
  std::vector<unsigned int> &allRevNodeCorrPrimeSize) {
#ifndef TTK_ENABLE_TORCH
  TTK_FORCE_USE(originsTrees);
  TTK_FORCE_USE(originsPrimeTrees);
  TTK_FORCE_USE(allRevNodeCorr);
  TTK_FORCE_USE(allRevNodeCorrPrime);
  TTK_FORCE_USE(allRevNodeCorrSize);
  TTK_FORCE_USE(allRevNodeCorrPrimeSize);
  printErr("This module requires Torch.");
#else
  // --- Preprocessing
  if(not isPersistenceDiagram_) {
    for(unsigned int i = 0; i < originsPrimeTrees.size(); ++i) {
      bool const useMinMax = true;
      bool const cleanTree = false;
      bool const pt = 0.0;
      std::vector<int> nodeCorr;
      preprocessingPipeline<float>(originsTrees[i], 0.0, 100.0, 100.0,
                                   branchDecomposition_, useMinMax, cleanTree,
                                   pt, nodeCorr, false);
      preprocessingPipeline<float>(originsPrimeTrees[i], 0.0, 100.0, 100.0,
                                   branchDecomposition_, useMinMax, cleanTree,
                                   pt, nodeCorr, false);
    }
  }
  mergeTreesToTorchTrees(originsTrees, origins_, normalizedWasserstein_,
                         allRevNodeCorr, allRevNodeCorrSize);
  mergeTreesToTorchTrees(originsPrimeTrees, originsPrime_,
                         normalizedWasserstein_, allRevNodeCorrPrime,
                         allRevNodeCorrPrimeSize);

  // --- Execute
  if(allAlphas_[0].size() != originsPrime_.size()) {
    customAlphas_.resize(allAlphas_.size());
    for(unsigned int i = 0; i < customAlphas_.size(); ++i)
      customAlphas_[i] = std::vector<float>(
        allAlphas_[i][0].data_ptr<float>(),
        allAlphas_[i][0].data_ptr<float>() + allAlphas_[i][0].numel());
    allAlphas_.clear();
    createCustomRecs(origins_, originsPrime_);
  } else {
    recs_.resize(allAlphas_.size());
    for(unsigned int i = 0; i < recs_.size(); ++i) {
      recs_[i].resize(allAlphas_[i].size());
      for(unsigned int l = 0; l < allAlphas_[i].size(); ++l) {
        torch::Tensor act
          = (activate_ ? activation(allAlphas_[i][l]) : allAlphas_[i][l]);
        getMultiInterpolation(
          originsPrime_[l], vSPrimeTensor_[l], act, recs_[i][l]);
      }
    }
  }

  // --- Postprocessing
  for(unsigned int l = 0; l < origins_.size(); ++l) {
    postprocessingPipeline<float>(&(origins_[l].mTree.tree));
    postprocessingPipeline<float>(&(originsPrime_[l].mTree.tree));
  }
  if(!recs_.empty()) {
    for(unsigned int j = 0; j < recs_[0].size(); ++j) {
      for(unsigned int i = 0; i < recs_.size(); ++i) {
        wae::fixTreePrecisionScalars(recs_[i][j].mTree);
        postprocessingPipeline<float>(&(recs_[i][j].mTree.tree));
      }
    }
  }
#endif
}
