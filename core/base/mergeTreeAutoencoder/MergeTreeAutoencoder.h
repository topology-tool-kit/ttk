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
#include <MergeTreeAxesAlgorithmBase.h>
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
                               public MergeTreeAxesAlgorithmBase {

  protected:
    bool doCompute_;
    bool hasComputedOnce_ = false;

    // Model hyper-parameters;
    int encoderNoLayers_ = 1;
    bool scaleLayerAfterLatent_ = false;
    unsigned int inputNumberOfAxes_ = 16;
    double inputOriginPrimeSizePercent_ = 15;
    double latentSpaceOriginPrimeSizePercent_ = 10;
    unsigned int minIteration_ = 0;
    unsigned int maxIteration_ = 0;
    unsigned int iterationGap_ = 100;
    double batchSize_ = 1;
    int optimizer_ = 0;
    double gradientStepSize_ = 0.1;
    double beta1_ = 0.9;
    double beta2_ = 0.999;
    double reconstructionLossWeight_ = 1;
    double trackingLossWeight_ = 0;
    double metricLossWeight_ = 0;
    double clusteringLossWeight_ = 0;
    float clusteringLossTemp_ = 10;
    bool customLossDynamicWeight_ = false;
    bool customLossSpace_ = false;
    bool customLossActivate_ = false;
    bool normalizeMetricLoss_ = false;
    unsigned int noInit_ = 4;
    bool euclideanVectorsInit_ = false;
    bool initOriginPrimeStructByCopy_ = true;
    bool trackingLossDecoding_ = false;
    double trackingLossInitRandomness_ = 0.0;
    bool activate_ = true;
    unsigned int activationFunction_ = 1;
    bool activateOutputInit_ = false;

    bool createOutput_ = true;

    // Old hyper-parameters
    bool fullSymmetricAE_ = false;

#ifdef TTK_ENABLE_TORCH
    // Model optimized parameters
    std::vector<torch::Tensor> vSTensor_, vSPrimeTensor_, vS2Tensor_,
      vS2PrimeTensor_, latentCentroids_;
    std::vector<mtu::TorchMergeTree<float>> origins_, originsPrime_, origins2_,
      origins2Prime_;

    std::vector<mtu::TorchMergeTree<float>> originsCopy_, originsPrimeCopy_;

    // Filled by the algorithm
    std::vector<std::vector<torch::Tensor>> allAlphas_, allScaledAlphas_,
      allActAlphas_, allActScaledAlphas_;
    std::vector<std::vector<mtu::TorchMergeTree<float>>> recs_, recs2_;
    std::vector<mtu::TorchMergeTree<float>> customRecs_;
#endif

    // Filled by the algorithm
    unsigned noLayers_;
    double baseRecLoss_, baseRecLoss2_;
    float bestLoss_;
    std::vector<unsigned int> clusterAsgn_;
    std::vector<std::vector<float>> distanceMatrix_, customAlphas_;
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      baryMatchings_L0_, baryMatchings2_L0_;
    std::vector<double> inputToBaryDistances_L0_;

    // Tracking matchings
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      originsMatchings_, reconstMatchings_, customMatchings_;
    std::vector<
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>>
      dataMatchings_;
    std::vector<std::vector<double>> branchesCorrelationMatrix_,
      persCorrelationMatrix_;

    // Testing
    double t_allVectorCopy_time_ = 0.0;
    std::vector<unsigned int> originsNoZeroGrad_, originsPrimeNoZeroGrad_,
      vSNoZeroGrad_, vSPrimeNoZeroGrad_, origins2NoZeroGrad_,
      origins2PrimeNoZeroGrad_, vS2NoZeroGrad_, vS2PrimeNoZeroGrad_;
    bool outputInit_ = true;
#ifdef TTK_ENABLE_TORCH
    std::vector<mtu::TorchMergeTree<float>> initOrigins_, initOriginsPrime_,
      initRecs_;
#endif
    float bigValuesThreshold_ = 0;

  public:
    MergeTreeAutoencoder();

#ifdef TTK_ENABLE_TORCH
    //  -----------------------------------------------------------------------
    //  --- Init
    //  -----------------------------------------------------------------------
    void initOutputBasisTreeStructure(mtu::TorchMergeTree<float> &originPrime,
                                      bool isJT,
                                      mtu::TorchMergeTree<float> &baseOrigin);

    void initOutputBasis(unsigned int l, unsigned int dim, unsigned int dim2);

    void initOutputBasisVectors(unsigned int l,
                                torch::Tensor &w,
                                torch::Tensor &w2);

    void initOutputBasisVectors(unsigned int l,
                                unsigned int dim,
                                unsigned int dim2);

    void initInputBasisOrigin(
      std::vector<ftm::MergeTree<float>> &treesToUse,
      std::vector<ftm::MergeTree<float>> &trees2ToUse,
      double barycenterSizeLimitPercent,
      unsigned int barycenterMaxNoPairs,
      unsigned int barycenterMaxNoPairs2,
      mtu::TorchMergeTree<float> &origin,
      mtu::TorchMergeTree<float> &origin2,
      std::vector<double> &inputToBaryDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2);

    void initInputBasisVectors(
      std::vector<mtu::TorchMergeTree<float>> &tmTreesToUse,
      std::vector<mtu::TorchMergeTree<float>> &tmTrees2ToUse,
      std::vector<ftm::MergeTree<float>> &treesToUse,
      std::vector<ftm::MergeTree<float>> &trees2ToUse,
      mtu::TorchMergeTree<float> &origin,
      mtu::TorchMergeTree<float> &origin2,
      unsigned int noVectors,
      std::vector<std::vector<torch::Tensor>> &allAlphasInit,
      unsigned int l,
      std::vector<double> &inputToBaryDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2,
      torch::Tensor &vSTensor,
      torch::Tensor &vS2Tensor);

    void initClusteringLossParameters();

    float initParameters(std::vector<mtu::TorchMergeTree<float>> &trees,
                         std::vector<mtu::TorchMergeTree<float>> &trees2,
                         bool computeReconstructionError = false);

    void initStep(std::vector<mtu::TorchMergeTree<float>> &trees,
                  std::vector<mtu::TorchMergeTree<float>> &trees2);

    //  -----------------------------------------------------------------------
    //  --- Interpolation
    //  -----------------------------------------------------------------------
    void interpolationDiagonalProjection(
      mtu::TorchMergeTree<float> &interpolationTensor);

    void
      interpolationNestingProjection(mtu::TorchMergeTree<float> &interpolation);

    void interpolationProjection(mtu::TorchMergeTree<float> &interpolation);

    void getMultiInterpolation(mtu::TorchMergeTree<float> &origin,
                               torch::Tensor &vS,
                               torch::Tensor &alphas,
                               mtu::TorchMergeTree<float> &interpolation);

    //  -----------------------------------------------------------------------
    //  --- Forward
    //  -----------------------------------------------------------------------
    void getAlphasOptimizationTensors(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &origin,
      torch::Tensor &vSTensor,
      mtu::TorchMergeTree<float> &interpolated,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      torch::Tensor &reorderedTreeTensor,
      torch::Tensor &deltaOrigin,
      torch::Tensor &deltaA,
      torch::Tensor &originTensor_f,
      torch::Tensor &vSTensor_f);

    void computeAlphas(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &origin,
      torch::Tensor &vSTensor,
      mtu::TorchMergeTree<float> &interpolated,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &origin2,
      torch::Tensor &vS2Tensor,
      mtu::TorchMergeTree<float> &interpolated2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2,
      torch::Tensor &alphasOut);

    float assignmentOneData(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &origin,
      torch::Tensor &vSTensor,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &origin2,
      torch::Tensor &vS2Tensor,
      unsigned int k,
      torch::Tensor &alphasInit,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching2,
      torch::Tensor &bestAlphas,
      bool isCalled = false);

    float assignmentOneData(mtu::TorchMergeTree<float> &tree,
                            mtu::TorchMergeTree<float> &origin,
                            torch::Tensor &vSTensor,
                            mtu::TorchMergeTree<float> &tree2,
                            mtu::TorchMergeTree<float> &origin2,
                            torch::Tensor &vS2Tensor,
                            unsigned int k,
                            torch::Tensor &alphasInit,
                            torch::Tensor &bestAlphas,
                            bool isCalled = false);

    torch::Tensor activation(torch::Tensor &in);

    void outputBasisReconstruction(mtu::TorchMergeTree<float> &originPrime,
                                   torch::Tensor &vSPrimeTensor,
                                   mtu::TorchMergeTree<float> &origin2Prime,
                                   torch::Tensor &vS2PrimeTensor,
                                   torch::Tensor &alphas,
                                   mtu::TorchMergeTree<float> &out,
                                   mtu::TorchMergeTree<float> &out2,
                                   bool activate = true);

    bool forwardOneLayer(mtu::TorchMergeTree<float> &tree,
                         mtu::TorchMergeTree<float> &origin,
                         torch::Tensor &vSTensor,
                         mtu::TorchMergeTree<float> &originPrime,
                         torch::Tensor &vSPrimeTensor,
                         mtu::TorchMergeTree<float> &tree2,
                         mtu::TorchMergeTree<float> &origin2,
                         torch::Tensor &vS2Tensor,
                         mtu::TorchMergeTree<float> &origin2Prime,
                         torch::Tensor &vS2PrimeTensor,
                         unsigned int k,
                         torch::Tensor &alphasInit,
                         mtu::TorchMergeTree<float> &out,
                         mtu::TorchMergeTree<float> &out2,
                         torch::Tensor &bestAlphas,
                         float &bestDistance);

    bool forwardOneLayer(mtu::TorchMergeTree<float> &tree,
                         mtu::TorchMergeTree<float> &origin,
                         torch::Tensor &vSTensor,
                         mtu::TorchMergeTree<float> &originPrime,
                         torch::Tensor &vSPrimeTensor,
                         mtu::TorchMergeTree<float> &tree2,
                         mtu::TorchMergeTree<float> &origin2,
                         torch::Tensor &vS2Tensor,
                         mtu::TorchMergeTree<float> &origin2Prime,
                         torch::Tensor &vS2PrimeTensor,
                         unsigned int k,
                         torch::Tensor &alphasInit,
                         mtu::TorchMergeTree<float> &out,
                         mtu::TorchMergeTree<float> &out2,
                         torch::Tensor &bestAlphas);

    bool forwardOneData(mtu::TorchMergeTree<float> &tree,
                        mtu::TorchMergeTree<float> &tree2,
                        unsigned int treeIndex,
                        unsigned int k,
                        std::vector<torch::Tensor> &alphasInit,
                        mtu::TorchMergeTree<float> &out,
                        mtu::TorchMergeTree<float> &out2,
                        std::vector<torch::Tensor> &dataAlphas,
                        std::vector<mtu::TorchMergeTree<float>> &outs,
                        std::vector<mtu::TorchMergeTree<float>> &outs2);

    bool forwardStep(
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<unsigned int> &indexes,
      unsigned int k,
      std::vector<std::vector<torch::Tensor>> &allAlphasInit,
      bool computeReconstructionError,
      std::vector<mtu::TorchMergeTree<float>> &outs,
      std::vector<mtu::TorchMergeTree<float>> &outs2,
      std::vector<std::vector<torch::Tensor>> &bestAlphas,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      float &loss);

    bool forwardStep(std::vector<mtu::TorchMergeTree<float>> &trees,
                     std::vector<mtu::TorchMergeTree<float>> &trees2,
                     std::vector<unsigned int> &indexes,
                     unsigned int k,
                     std::vector<std::vector<torch::Tensor>> &allAlphasInit,
                     std::vector<mtu::TorchMergeTree<float>> &outs,
                     std::vector<mtu::TorchMergeTree<float>> &outs2,
                     std::vector<std::vector<torch::Tensor>> &bestAlphas);

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
      torch::optim::Optimizer &optimizer,
      std::vector<unsigned int> &indexes,
      torch::Tensor &metricLoss,
      torch::Tensor &clusteringLoss,
      torch::Tensor &trackingLoss);

    //  -----------------------------------------------------------------------
    //  --- Projection
    //  -----------------------------------------------------------------------
    void projectionStep();

    //  -----------------------------------------------------------------------
    //  --- Convergence
    //  -----------------------------------------------------------------------
    float computeOneLoss(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &out,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &out2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2);

    float computeLoss(
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &outs,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<mtu::TorchMergeTree<float>> &outs2,
      std::vector<unsigned int> &indexes,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2);

    bool isBestLoss(float loss, float &minLoss, unsigned int &cptBlocked);

    bool convergenceStep(float loss,
                         float &oldLoss,
                         float &minLoss,
                         unsigned int &cptBlocked);

    //  -----------------------------------------------------------------------
    //  --- Main Functions
    //  -----------------------------------------------------------------------
    void fit(std::vector<ftm::MergeTree<float>> &trees,
             std::vector<ftm::MergeTree<float>> &trees2);

    //  -----------------------------------------------------------------------
    //  --- Custom Losses
    //  -----------------------------------------------------------------------
    double getCustomLossDynamicWeight(double recLoss, double &baseLoss);

    void getDistanceMatrix(std::vector<mtu::TorchMergeTree<float>> &tmts,
                           std::vector<std::vector<float>> &distanceMatrix,
                           bool useDoubleInput = false,
                           bool isFirstInput = true);

    void getDistanceMatrix(std::vector<mtu::TorchMergeTree<float>> &tmts,
                           std::vector<mtu::TorchMergeTree<float>> &tmts2,
                           std::vector<std::vector<float>> &distanceMatrix);

    void getDifferentiableDistanceFromMatchings(
      mtu::TorchMergeTree<float> &tree1,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &tree1_2,
      mtu::TorchMergeTree<float> &tree2_2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings2,
      torch::Tensor &tensorDist,
      bool doSqrt);

    void getDifferentiableDistance(mtu::TorchMergeTree<float> &tree1,
                                   mtu::TorchMergeTree<float> &tree2,
                                   mtu::TorchMergeTree<float> &tree1_2,
                                   mtu::TorchMergeTree<float> &tree2_2,
                                   torch::Tensor &tensorDist,
                                   bool isCalled,
                                   bool doSqrt);

    void getDifferentiableDistance(mtu::TorchMergeTree<float> &tree1,
                                   mtu::TorchMergeTree<float> &tree2,
                                   torch::Tensor &tensorDist,
                                   bool isCalled,
                                   bool doSqrt);

    void getDifferentiableDistanceMatrix(
      std::vector<mtu::TorchMergeTree<float> *> &trees,
      std::vector<mtu::TorchMergeTree<float> *> &trees2,
      std::vector<std::vector<torch::Tensor>> &outDistMat);

    void getAlphasTensor(std::vector<std::vector<torch::Tensor>> &alphas,
                         std::vector<unsigned int> &indexes,
                         unsigned int layerIndex,
                         torch::Tensor &alphasOut);

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
    void
      createCustomRecs(std::vector<mtu::TorchMergeTree<float>> &origins,
                       std::vector<mtu::TorchMergeTree<float>> &originsPrime);

    void computeTrackingInformation();

    void
      createScaledAlphas(std::vector<std::vector<torch::Tensor>> &alphas,
                         std::vector<torch::Tensor> &vSTensor,
                         std::vector<std::vector<torch::Tensor>> &scaledAlphas);

    void createScaledAlphas();

    void createActivatedAlphas();

    //  -----------------------------------------------------------------------
    //  --- Utils
    //  -----------------------------------------------------------------------
    void copyParams(std::vector<mtu::TorchMergeTree<float>> &srcOrigins,
                    std::vector<mtu::TorchMergeTree<float>> &srcOriginsPrime,
                    std::vector<torch::Tensor> &srcVS,
                    std::vector<torch::Tensor> &srcVSPrime,
                    std::vector<mtu::TorchMergeTree<float>> &srcOrigins2,
                    std::vector<mtu::TorchMergeTree<float>> &srcOrigins2Prime,
                    std::vector<torch::Tensor> &srcVS2,
                    std::vector<torch::Tensor> &srcVS2Prime,
                    std::vector<std::vector<torch::Tensor>> &srcAlphas,
                    std::vector<mtu::TorchMergeTree<float>> &dstOrigins,
                    std::vector<mtu::TorchMergeTree<float>> &dstOriginsPrime,
                    std::vector<torch::Tensor> &dstVS,
                    std::vector<torch::Tensor> &dstVSPrime,
                    std::vector<mtu::TorchMergeTree<float>> &dstOrigins2,
                    std::vector<mtu::TorchMergeTree<float>> &dstOrigins2Prime,
                    std::vector<torch::Tensor> &dstVS2,
                    std::vector<torch::Tensor> &dstVS2Prime,
                    std::vector<std::vector<torch::Tensor>> &dstAlphas);

    void copyParams(std::vector<std::vector<mtu::TorchMergeTree<float>>> &src,
                    std::vector<std::vector<mtu::TorchMergeTree<float>>> &dst);

    unsigned int getLatentLayerIndex();

    //  -----------------------------------------------------------------------
    //  --- Testing
    //  -----------------------------------------------------------------------
    bool isTreeHasBigValues(ftm::MergeTree<float> &mTree,
                            float threshold = 10000);
#endif

    //  ---------------------------------------------------------------------------
    //  --- Main Functions
    //  ---------------------------------------------------------------------------
    void execute(std::vector<ftm::MergeTree<float>> &trees,
                 std::vector<ftm::MergeTree<float>> &trees2);
  }; // MergeTreeAutoencoder class

} // namespace ttk
