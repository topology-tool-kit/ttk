/// \ingroup base
/// \class ttk::MergeTreeAutoencoder
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// This module defines the %MergeTreeAutoencoder class that computes an
/// Auto-Encoding of merge trees or persistence diagrams.
///
/// \b Related \b publication: \n
/// "Wasserstein Auto-Encoders of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2023
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeWAE/">Merge
///   Tree Wasserstein Auto-Encoder example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramWAE/">Persistence
///   Diagram Wasserstein Auto-Encoder example</a> \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <MergeTreeNeuralLayer.h>
#include <MergeTreeNeuralNetwork.h>
#include <MergeTreeTorchUtils.h>

#ifdef TTK_ENABLE_TORCH
#include <torch/torch.h>
#endif

namespace ttk {

  /**
   * The MergeTreeAutoencoder class provides methods to compute an Auto-Encoding
   * of merge trees or persistence diagrams.
   */
  class MergeTreeAutoencoder : virtual public Debug,
                               public MergeTreeNeuralNetwork {

  protected:
    // Model hyper-parameters;
    int encoderNoLayers_ = 1;
    bool scaleLayerAfterLatent_ = false;
    unsigned int inputNumberOfAxes_ = 16;
    double inputOriginPrimeSizePercent_ = 15;
    double latentSpaceOriginPrimeSizePercent_ = 10;
    double reconstructionLossWeight_ = 1;
    double trackingLossWeight_ = 0;
    double metricLossWeight_ = 0;
    double clusteringLossWeight_ = 0;
    float clusteringLossTemp_ = 10;
    bool customLossDynamicWeight_ = false;
    bool customLossSpace_ = false;
    bool customLossActivate_ = false;
    bool normalizeMetricLoss_ = false;
    bool trackingLossDecoding_ = false;
    double trackingLossInitRandomness_ = 0.0;

    // Old hyper-parameters
    bool fullSymmetricAE_ = false;

#ifdef TTK_ENABLE_TORCH
    // Model optimized parameters
    std::vector<torch::Tensor> bestLatentCentroids_, latentCentroids_;

    std::vector<torch::Tensor> vSTensorCopy_, vSPrimeTensorCopy_;

    std::vector<mtu::TorchMergeTree<float>> customRecs_;
#endif

    // Filled by the algorithm
    double baseRecLoss_, baseRecLoss2_;
    std::vector<std::vector<float>> distanceMatrix_, customAlphas_;

  public:
    MergeTreeAutoencoder();

#ifdef TTK_ENABLE_TORCH
    //  -----------------------------------------------------------------------
    //  --- Init
    //  -----------------------------------------------------------------------
    void initClusteringLossParameters();

    bool initResetOutputBasis(unsigned int l,
                              unsigned int layerNoAxes,
                              double layerOriginPrimeSizePercent,
                              std::vector<mtu::TorchMergeTree<float>> &trees,
                              std::vector<mtu::TorchMergeTree<float>> &trees2,
                              std::vector<bool> &isTrain) override;

    void initOutputBasisSpecialCase(
      unsigned int l,
      unsigned int layerNoAxes,
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &trees2);

    float initParameters(std::vector<mtu::TorchMergeTree<float>> &trees,
                         std::vector<mtu::TorchMergeTree<float>> &trees2,
                         std::vector<bool> &isTrain,
                         bool computeError = false) override;

    //  -----------------------------------------------------------------------
    //  --- Backward
    //  -----------------------------------------------------------------------
    bool backwardStep(
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &outs,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<mtu::TorchMergeTree<float>> &outs2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      std::vector<std::vector<torch::Tensor>> &alphas,
      torch::optim::Optimizer &optimizer,
      std::vector<unsigned int> &indexes,
      std::vector<bool> &isTrain,
      std::vector<torch::Tensor> &torchCustomLoss) override;

    //  -----------------------------------------------------------------------
    //  --- Convergence
    //  -----------------------------------------------------------------------
    float computeOneLoss(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &out,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &out2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2,
      std::vector<torch::Tensor> &alphas,
      unsigned int treeIndex) override;

    //  -----------------------------------------------------------------------
    //  --- Main Functions
    //  -----------------------------------------------------------------------
    void
      customInit(std::vector<mtu::TorchMergeTree<float>> &torchTrees,
                 std::vector<mtu::TorchMergeTree<float>> &torchTrees2) override;

    void addCustomParameters(std::vector<torch::Tensor> &parameters) override;

    void computeCustomLosses(
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
      std::vector<std::vector<torch::Tensor>> &bestAlphas,
      std::vector<unsigned int> &indexes,
      std::vector<bool> &isTrain,
      unsigned int iteration,
      std::vector<std::vector<float>> &gapCustomLosses,
      std::vector<std::vector<float>> &iterationCustomLosses,
      std::vector<torch::Tensor> &torchCustomLoss) override;

    float computeIterationTotalLoss(
      float iterationLoss,
      std::vector<std::vector<float>> &iterationCustomLosses,
      std::vector<float> &iterationCustomLoss) override;

    void printCustomLosses(std::vector<float> &customLoss,
                           std::stringstream &prefix,
                           const debug::Priority &priority
                           = debug::Priority::INFO) override;

    void
      printGapLoss(float loss,
                   std::vector<std::vector<float>> &gapCustomLosses) override;

    //  -----------------------------------------------------------------------
    //  --- Custom Losses
    //  -----------------------------------------------------------------------
    double getCustomLossDynamicWeight(double recLoss, double &baseLoss);

    void computeMetricLoss(
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
      std::vector<std::vector<torch::Tensor>> alphas,
      std::vector<std::vector<float>> &baseDistanceMatrix,
      std::vector<unsigned int> &indexes,
      torch::Tensor &metricLoss);

    void computeClusteringLoss(std::vector<std::vector<torch::Tensor>> &alphas,
                               std::vector<unsigned int> &indexes,
                               torch::Tensor &clusteringLoss,
                               torch::Tensor &asgn);

    void computeTrackingLoss(torch::Tensor &trackingLoss);

    //  ---------------------------------------------------------------------------
    //  --- End Functions
    //  ---------------------------------------------------------------------------
    void createCustomRecs();

    //  -----------------------------------------------------------------------
    //  --- Utils
    //  -----------------------------------------------------------------------
    unsigned int getLatentLayerIndex();

    void copyCustomParams(bool get) override;

    //  ---------------------------------------------------------------------------
    //  --- Main Functions
    //  ---------------------------------------------------------------------------
    void
      executeEndFunction(std::vector<ftm::MergeTree<float>> &trees,
                         std::vector<ftm::MergeTree<float>> &trees2) override;
#endif
  }; // MergeTreeAutoencoder class

} // namespace ttk
