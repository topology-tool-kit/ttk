/// \ingroup base
/// \class ttk::MergeTreeNeuralNetwork
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// This module defines the %MergeTreeNeuralNetwork class providing functions to
/// define a neural network able to process merge trees or persistence diagrams.
///
/// This is an abstract class, to implement a derived class you need to define
/// the following functions:
///
/// - "initParameters" : initializes the network, like the different layers.
/// A strategy to initialize a sequence of layers consist in initializing a
/// first layer with the input topological representations, then pass them
/// through this layer to initialize the second one and so on.
/// A simple loop (whose number of iterations corresponds to the number of
/// layers) can do this using the "initInputBasis" and the "initOutputBasis"
/// function, then the "initGetReconstructed" function to pass the
/// representations to the layer that just have been initialized.
///
/// - "initResetOutputBasis" : please refer to the documentation of this
/// function.
///
/// - "customInit" : called just before the "initStep" function (that call the
/// "initParameters" function), is is intended to do custom operations depending
/// on the architecture and the optimization you want to define (such as
/// computing the distance matrix for the metric loss in the autoencoder case).
/// This function can be empty.
///
/// - "backwardStep" : optimizes the parameters of the network. A loss using
/// differentiable torch operations should be computed using the output of some
/// layers of the network (usually the output of the last layer but it can also
/// be any other layers). You can either use the torch coordinates of the
/// representations in a layer or their torch tensors to compute the loss. Then
/// use the torch::Tensor "backward" function to compute the gradients, then the
/// torch::optim::Optimizer "step" function to update the model parameters,
/// after this, the torch::optim::Optimizer "zero_grad" function should be
/// called to reset the gradient. If you have correctly created the
/// MergeTreeNeuralLayer objects (refer to the corresponding class
/// documentation), basically by calling the "requires_grad" function (with true
/// as parameter) for each layer after initializing its parameters, then
/// everything would be automatically handled to backpropagate the gradient of
/// the loss through the layers.
///
/// - "addCustomParameters" : adds custom parameters to the parameter list that
/// will be given to the optimizer, depending on the architecture and the
/// optimization you want to define (such as the centroids for the cluster loss
/// in the autoencoder case). This function can be empty.
///
/// - "computeOneLoss" : computes the loss for one input topological
/// representation, the loss computed here does not need to be differentiable
/// because it will only be used to print it in the console and to check
/// convergence of the method (i.e. it is not called in the "backwardStep"
/// function).
///
/// - "computeCustomLosses" : computes custom losses for all input topological
/// representations depending on the architecture and the optimization you want
/// to define (such as the clustering and the metric loss in the autoendoer
/// case). Like "computeOneLoss", the losses do not need to be differentiable
/// because they will only be used to print them in the console and to check
/// convergence of the method. This function can be empty.
///
/// - "computeIterationTotalLoss"
///
/// - "printCustomLosses" : prints the custom loss depending on the architecture
/// and the optimization you want to define (such as the clustering and the
/// metric loss in the autoendoer case). This function can be empty.
///
/// - "printGapLoss" : prints the "gap" loss, the aggregated loss over
/// iterationGap_ iterations.
///
/// - "copyCustomParams" : copy the custom parameters (for instance to save them
/// when a better loss is reached during the optimization) depending on the
/// architecture and the optimization you want to define (such as the centroids
/// for the cluster loss in the autoencoder case). This function can be empty.
///
/// - "executeEndFunction" : does specific operations at the end of the
/// optimization, such as calling the "computeTrackingInformation" and the
/// "computeCorrelationMatrix" functions.
///
/// \b Related \b publication: \n
/// "Wasserstein Auto-Encoders of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont, Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2023
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <MergeTreeNeuralBase.h>
#include <MergeTreeNeuralLayer.h>
#include <MergeTreeTorchUtils.h>

#ifdef TTK_ENABLE_TORCH
#include <torch/torch.h>
#endif

namespace ttk {

  /**
   * The MergeTreeNeuralNetwork class provides methods to define a neural
   * network able to process merge trees or persistence diagrams.
   */
  class MergeTreeNeuralNetwork : virtual public Debug,
                                 public MergeTreeNeuralBase {

  protected:
    // Minimum number of iterations to run
    unsigned int minIteration_ = 0;
    // Maximum number of iterations to run
    unsigned int maxIteration_ = 0;
    // Number of iterations between each print
    unsigned int iterationGap_ = 100;
    // Batch size between 0 and 1
    double batchSize_ = 1;
    // Optimizer
    // 0 : Adam
    // 1 : Stochastic Gradient Descent
    // 2 : RMS Prop
    int optimizer_ = 0;
    // Gradient Step/Learning rate
    double gradientStepSize_ = 0.1;
    // Adam parameters
    double beta1_ = 0.9;
    double beta2_ = 0.999;
    // Number of initializations to do (the better will be kept)
    unsigned int noInit_ = 4;
    // If activation functions should be used during the initialization
    bool activateOutputInit_ = false;
    // Limit in the size of the origin in output basis as a percentage of the
    // input total number of nodes
    double originPrimeSizePercent_ = 15;
    // Proportion between the train set and the validation/test set
    double trainTestSplit_ = 1.0;
    // If the input data should be shuffled before splitted
    bool shuffleBeforeSplit_ = true;

    bool createOutput_ = true;

#ifdef TTK_ENABLE_TORCH
    // Model optimized parameters
    std::vector<MergeTreeNeuralLayer> layers_;

    // Filled by the algorithm
    std::vector<std::vector<torch::Tensor>> allAlphas_, allScaledAlphas_,
      allActAlphas_, allActScaledAlphas_;
    std::vector<std::vector<mtu::TorchMergeTree<float>>> recs_, recs2_;

    std::vector<mtu::TorchMergeTree<float>> originsCopy_, originsPrimeCopy_;
#endif

    // Tracking matchings
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      originsMatchings_, reconstMatchings_, customMatchings_;
    std::vector<
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>>
      dataMatchings_;

    // Filled by the algorithm
    unsigned noLayers_;
    float bestLoss_;
    std::vector<unsigned int> clusterAsgn_;
    std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
      baryMatchings_L0_, baryMatchings2_L0_;
    std::vector<double> inputToBaryDistances_L0_;
    std::vector<std::vector<double>> branchesCorrelationMatrix_,
      persCorrelationMatrix_;

    // Testing
    double t_allVectorCopy_time_ = 0.0;
    std::vector<unsigned int> originsNoZeroGrad_, originsPrimeNoZeroGrad_,
      vSNoZeroGrad_, vSPrimeNoZeroGrad_, origins2NoZeroGrad_,
      origins2PrimeNoZeroGrad_, vS2NoZeroGrad_, vS2PrimeNoZeroGrad_;

  public:
    MergeTreeNeuralNetwork();

#ifdef TTK_ENABLE_TORCH
    //  -----------------------------------------------------------------------
    //  --- Init
    //  -----------------------------------------------------------------------
    /**
     * @brief Initialize an input basis.
     *
     * @param[in] l index of the layer to initialize.
     * @param[in] layerNoAxes number of axes in the basis.
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[out] allAlphasInit best estimation of the coordinates of each
     * input tree in the basis.
     */
    void initInputBasis(unsigned int l,
                        unsigned int layerNoAxes,
                        std::vector<mtu::TorchMergeTree<float>> &trees,
                        std::vector<mtu::TorchMergeTree<float>> &trees2,
                        std::vector<bool> &isTrain,
                        std::vector<std::vector<torch::Tensor>> &allAlphasInit);

    /**
     * @brief Initialize an output basis.
     *
     * @param[in] l index of the layer to initialize.
     * @param[in] layerOriginPrimeSizePercent the maximum number of nodes
     * allowed for the origin as a percentage of the total number of nodes
     * in the input trees (0 for no effect).
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     */
    void initOutputBasis(unsigned int l,
                         double layerOriginPrimeSizePercent,
                         std::vector<mtu::TorchMergeTree<float>> &trees,
                         std::vector<mtu::TorchMergeTree<float>> &trees2,
                         std::vector<bool> &isTrain);

    /**
     * @brief It is possible for the output basis of a layer to be badly
     * initialized such a merge tree has no nodes after being passed through it.
     * This function is intended to initialize a new output basis to avoid this
     * problem. It is called in the initGetReconstructed function during the
     * initialization procedure.
     * A simple call to initOutputBasis could be enough.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] l index of the layer to initialize.
     * @param[in] layerNoAxes number of axes in the basis.
     * @param[in] layerOriginPrimeSizePercent the maximum number of nodes
     * allowed for the origin as a percentage of the total number of nodes
     * in the input trees (0 for no effect).
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     *
     * @return true if it is not possible to initialize a new output basis,
     * false otherwise.
     */
    virtual bool
      initResetOutputBasis(unsigned int l,
                           unsigned int layerNoAxes,
                           double layerOriginPrimeSizePercent,
                           std::vector<mtu::TorchMergeTree<float>> &trees,
                           std::vector<mtu::TorchMergeTree<float>> &trees2,
                           std::vector<bool> &isTrain)
      = 0;

    /**
     * @brief Pass trees through a layer that just have been initialized.
     *
     * @param[in] l index of the layer being initialized.
     * @param[in] layerNoAxes number of axes in the basis.
     * @param[in] layerOriginPrimeSizePercent the maximum number of nodes
     * allowed for the origin as a percentage of the total number of nodes
     * in the input trees (0 for no effect).
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[out] recs the input trees after being passed through the layer.
     * @param[out] recs2 same as recs but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[out] allAlphasInit best estimation of the coordinates of each
     * input tree in the basis.
     *
     * @return true if one output merge tree has no nodes.
     */
    bool initGetReconstructed(
      unsigned int l,
      unsigned int layerNoAxes,
      double layerOriginPrimeSizePercent,
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<bool> &isTrain,
      std::vector<mtu::TorchMergeTree<float>> &recs,
      std::vector<mtu::TorchMergeTree<float>> &recs2,
      std::vector<std::vector<torch::Tensor>> &allAlphasInit);

    /**
     * @brief Initialize the parameters of the network.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[in] computeError boolean stating whether or not the initialization
     * error should be computed.
     *
     * @return return an initialization error, a value stating how bad the
     * initialization is (the higher the worst).
     */
    virtual float
      initParameters(std::vector<mtu::TorchMergeTree<float>> &trees,
                     std::vector<mtu::TorchMergeTree<float>> &trees2,
                     std::vector<bool> &isTrain,
                     bool computeError = false)
      = 0;

    /**
     * @brief Initialize the parameters of the network a specific number of
     * times and keep the best one.
     *
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     */
    void initStep(std::vector<mtu::TorchMergeTree<float>> &trees,
                  std::vector<mtu::TorchMergeTree<float>> &trees2,
                  std::vector<bool> &isTrain);

    /**
     * @brief Pass all the parameters from this class to a MergeTreeNeuralLayer
     * object.
     * Should be called just after creating a layer and before using it.
     *
     * @param[in] layer layer to process.
     */
    void passLayerParameters(MergeTreeNeuralLayer &layer);

    //  -----------------------------------------------------------------------
    //  --- Forward
    //  -----------------------------------------------------------------------
    /**
     * @brief Pass one tree through all layers of the network.
     *
     * @param[in] tree input tree.
     * @param[in] trees2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] treeIndex index of the input tree in the ensemble.
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] alphasInit the initial coordinates (for each layer) to use
     * when estimating the coordinates in the input basis of the input merge
     * tree.
     * @param[out] out the final output merge tree.
     * @param[out] out2 same as out but for second input, if any (i.e. when join
     * trees and split trees are given).
     * @param[out] dataAlphas the best estimated coordinates of the input tree
     * at each layer.
     * @param[out] outs the output merge tree of each layer (except the last
     * one).
     * @param[out] outs2 same as outs but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] train true if the input merge tree is in the training set
     * (false if validation/testing set).
     *
     * @return true if the output merge tree of a layer has no nodes.
     */
    bool forwardOneData(mtu::TorchMergeTree<float> &tree,
                        mtu::TorchMergeTree<float> &tree2,
                        unsigned int treeIndex,
                        unsigned int k,
                        std::vector<torch::Tensor> &alphasInit,
                        mtu::TorchMergeTree<float> &out,
                        mtu::TorchMergeTree<float> &out2,
                        std::vector<torch::Tensor> &dataAlphas,
                        std::vector<mtu::TorchMergeTree<float>> &outs,
                        std::vector<mtu::TorchMergeTree<float>> &outs2,
                        bool train = false);

    /**
     * @brief Pass all trees through all layers of the network.
     *
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] indexes batch indexes of the input trees to process.
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] allAlphasInit the initial coordinates for each input tree for
     * each layer to use when estimating the coordinates in the input basis.
     * @param[in] computeError true if the loss for each processed input tree
     * should be computed, false otherwise.
     * @param[out] outs the final output of each merge tree.
     * @param[out] outs2 same as out but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[out] bestAlphas the best estimated coordinates for each input tree
     * at each layer.
     * @param[out] layersOuts the output for each input merge tree of each layer
     * (except the last one).
     * @param[out] layersOuts2 same as layersOuts but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[out] matchings stores a matching for each processed input tree if
     * the loss involves an assignment problem.
     * @param[out] matchings2 same as matchings but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[out] loss the loss of the trees in the training set.
     * @param[out] testLoss the loss of the trees not in the training set.
     *
     * @return true if an output merge tree of a layer has no nodes.
     */
    bool forwardStep(
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<unsigned int> &indexes,
      std::vector<bool> &isTrain,
      unsigned int k,
      std::vector<std::vector<torch::Tensor>> &allAlphasInit,
      bool computeError,
      std::vector<mtu::TorchMergeTree<float>> &outs,
      std::vector<mtu::TorchMergeTree<float>> &outs2,
      std::vector<std::vector<torch::Tensor>> &bestAlphas,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &matchings2,
      float &loss,
      float &testLoss);

    /**
     * @brief Pass all trees through all layers of the network.
     *
     * @param[in] trees input trees.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] indexes batch indexes of the input trees to process.
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] allAlphasInit the initial coordinates for each input tree for
     * each layer to use when estimating the coordinates in the input basis.
     * @param[in] computeError true if the loss for each processed input tree
     * should be computed, false otherwise.
     * @param[out] outs the final output of each merge tree.
     * @param[out] outs2 same as out but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[out] bestAlphas the best estimated coordinates for each input tree
     * at each layer.
     * @param[out] layersOuts the output for each input merge tree of each layer
     * (except the last one).
     * @param[out] layersOuts2 same as layersOuts but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[out] matchings stores a matching for each processed input tree if
     * the loss involves an assignment problem.
     * @param[out] matchings2 same as matchings but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[out] loss the loss of the trees in the training set.
     *
     * @return true if an output merge tree of a layer has no nodes.
     */
    bool forwardStep(
      std::vector<mtu::TorchMergeTree<float>> &trees,
      std::vector<mtu::TorchMergeTree<float>> &trees2,
      std::vector<unsigned int> &indexes,
      unsigned int k,
      std::vector<std::vector<torch::Tensor>> &allAlphasInit,
      bool computeError,
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

    //  -----------------------------------------------------------------------
    //  --- Backward
    //  -----------------------------------------------------------------------
    /**
     * @brief Updates the parameters of the network to minimize the error.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] trees input trees.
     * @param[in] outs the final output of each merge tree.
     * @param[in] matchings stores a matching for each processed input tree if
     * the loss involves an assignment problem.
     * @param[in] trees2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] outs2 same as out but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] matchings2 same as matchings but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] alphas the best estimated coordinates for each input tree
     * at each layer.
     * @param[in] optimizer optimizer to use to modify the parameters.
     * @param[in] indexes batch indexes of the input trees to process.
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[in] torchCustomLoss custom losses that can be added to the
     * optimization (such as the loss ensuring the preservation of the clusters
     * or the distances in the autoencoder case).
     *
     * @return not used (false).
     */
    virtual bool backwardStep(
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
      std::vector<torch::Tensor> &torchCustomLoss)
      = 0;

    //  -----------------------------------------------------------------------
    //  --- Projection
    //  -----------------------------------------------------------------------
    /**
     * @brief Projection that ensures that the origins of the input and output
     * bases of each layer respect the elder rule.
     */
    void projectionStep();

    //  -----------------------------------------------------------------------
    //  --- Convergence
    //  -----------------------------------------------------------------------
    /**
     * @brief Computes the loss for one input tree.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] tree input tree.
     * @param[in] out the final output of the input merge tree.
     * @param[in] tree2 same as trees but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] out2 same as out but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] matchings stores a matching involving the input tree if
     * the loss involves an assignment problem.
     * @param[in] matchings2 same as matchings but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] alphas the best estimated coordinates for the input tree
     * at each layer.
     * @param[in] treeIndex index of the input tree in the ensemble.
     *
     * @return loss value.
     */
    virtual float computeOneLoss(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &out,
      mtu::TorchMergeTree<float> &tree2,
      mtu::TorchMergeTree<float> &out2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching2,
      std::vector<torch::Tensor> &alphas,
      unsigned int treeIndex)
      = 0;

    /**
     * @brief Tests if the current loss is the lowest one.
     *
     * @param[in] loss the current loss.
     * @param[in,out] minLoss the minimum loss that was reached, will be updated
     * if it is higher than the current loss.
     * @param[in,out] cptBlocked value that will be reset to 0 if the current
     * loss is better than the best one.
     *
     * @return true if the current loss is the lowest one.
     */
    bool isBestLoss(float loss, float &minLoss, unsigned int &cptBlocked);

    /**
     * @brief Tests if the optimization is done.
     *
     * @param[in] loss the current loss.
     * @param[in,out] oldLoss the loss of the previous iteration, will be set to
     * the current loss at the end of this function.
     * @param[in] minLoss the minimum loss that was reached.
     * @param[in,out] cptBlocked number of iterations during the minimum loss
     * was not updated.
     *
     * @return true if the optimization is done.
     */
    bool convergenceStep(float loss,
                         float &oldLoss,
                         float &minLoss,
                         unsigned int &cptBlocked);

    //  -----------------------------------------------------------------------
    //  --- Main Functions
    //  -----------------------------------------------------------------------
    /**
     * @brief Custom operations that need to be done before starting the
     * optimization.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] torchTrees the input trees.
     * @param[in] torchTrees2 same as torchTrees but for second input, if any
     * (i.e. when join trees and split trees are given).
     */
    virtual void
      customInit(std::vector<mtu::TorchMergeTree<float>> &torchTrees,
                 std::vector<mtu::TorchMergeTree<float>> &torchTrees2)
      = 0;

    /**
     * @brief This function adds parameters to the parameter list that will be
     * given to the optimizer. The vector should NOT be reinitialized in this
     * function, only use emplace_back.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in,out] parameters list of parameters to optimize that will be
     * given to the optimizer, should not be reset.
     */
    virtual void addCustomParameters(std::vector<torch::Tensor> &parameters)
      = 0;

    /**
     * @brief Compute the custom losses (such as the metric or cluster loss in
     * the autoencoder case).
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] layersOuts the output for each input merge tree of each layer
     * (except the last one).
     * @param[in] layersOuts2 same as layersOuts but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] bestAlphas the best estimated coordinates for each input tree
     * at each layer.
     * @param[in] indexes batch indexes of the input trees to process.
     * @param[in] isTrain vector stating for each input tree if it is in the
     * training set (true) or not (false).
     * @param[in] iteration iteration number.
     * @param[out] gapCustomLosses vector needing to be resized in order to have
     * a size corresponding to the number of different custom losses and should
     * be appended each custom loss computed here.
     * @param[out] iterationCustomLosses vector needing to be resized in order
     * to have a size corresponding to the number of different custom losses and
     * should be appended each custom loss computed here.
     * @param[out] torchCustomLoss vector needing to be resized in order to have
     * a size corresponding to the number of different custom losses and should
     * be appended each custom loss computed here.
     */
    virtual void computeCustomLosses(
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts,
      std::vector<std::vector<mtu::TorchMergeTree<float>>> &layersOuts2,
      std::vector<std::vector<torch::Tensor>> &bestAlphas,
      std::vector<unsigned int> &indexes,
      std::vector<bool> &isTrain,
      unsigned int iteration,
      std::vector<std::vector<float>> &gapCustomLosses,
      std::vector<std::vector<float>> &iterationCustomLosses,
      std::vector<torch::Tensor> &torchCustomLoss)
      = 0;

    /**
     * @brief Computes the total loss of the current iteration.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] iterationLoss loss of the current iteration.
     * @param[in] iterationCustomLosses all the custom losses of this iteration,
     * one for each different custom loss for each batch.
     * @param[out] iterationCustomLoss a vector containing the aggregated value
     * of each custom loss (usually the mean over each batch).
     *
     * @return the total loss of the current iteration.
     */
    virtual float computeIterationTotalLoss(
      float iterationLoss,
      std::vector<std::vector<float>> &iterationCustomLosses,
      std::vector<float> &iterationCustomLoss)
      = 0;

    /**
     * @brief Print the custom losses.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] customLoss a vector containing the aggregated value
     * of each custom loss (usually the mean over each batch).
     * @param[in] prefix the prefix of each message.
     * @param[in] priority the priority of the TTK message.
     */
    virtual void printCustomLosses(std::vector<float> &customLoss,
                                   std::stringstream &prefix,
                                   const debug::Priority &priority
                                   = debug::Priority::INFO)
      = 0;

    /**
     * @brief Print the gap loss (the aggregated loss over iterationGap_
     * iterations).
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] loss the loss to print.
     * @param[in] gapCustomLosses the custom losses to print, a vector
     * containing a vector for each custom loss.
     */
    virtual void printGapLoss(float loss,
                              std::vector<std::vector<float>> &gapCustomLosses)
      = 0;

    /**
     * @brief Initiliazes the network and trains it.
     *
     * @param[in] trees the input trees.
     * @param[in] trees2 same as torchTrees but for second input, if any (i.e.
     * when join trees and split trees are given).
     */
    void fit(std::vector<ftm::MergeTree<float>> &trees,
             std::vector<ftm::MergeTree<float>> &trees2);

    //  ---------------------------------------------------------------------------
    //  --- End Functions
    //  ---------------------------------------------------------------------------
    /**
     * @brief Computes the tracking information, i.e, the matchings between the
     * origin of the output basis of two consecutive layers (from the first one
     * to the one specified by the endLayer parameter), the matchings
     * between the input representations and the origin of the input basis of
     * the first layer, and the matching between the latter and the origin of
     * the output basis of the first layer.
     *
     * @param[in] endLayer layer number to stop the computation.
     */
    void computeTrackingInformation(unsigned int endLayer);

    /**
     * @brief Computes the correlation matrix between the pairs of the input
     * trees and the basis at the layer specified in parameter. The
     * "computeTrackingInformation" should have been called before.
     *
     * @param[in] trees the input trees.
     * @param[in] layer the layer at which the correlation should be computed.
     */
    void computeCorrelationMatrix(std::vector<ftm::MergeTree<float>> &trees,
                                  unsigned int layer);

    /**
     * @brief Scales the coordinates given in input by the norm of the basis at
     * the corresponding layer.
     *
     * @param[in] alphas coordinates for each input topological representation
     * for each layer.
     * @param[out] scaledAlphas scaled coordinates.
     */
    void
      createScaledAlphas(std::vector<std::vector<torch::Tensor>> &alphas,
                         std::vector<std::vector<torch::Tensor>> &scaledAlphas);

    /**
     * @brief Scales the coordinates of the input topological representations by
     * the norm of the basis at the corresponding layer.
     */
    void createScaledAlphas();

    /**
     * @brief Scales the activated coordinates (the coordinates passed through
     * the activation function) of the input topological representations by the
     * norm of the basis at the corresponding layer.
     */
    void createActivatedAlphas();

    //  -----------------------------------------------------------------------
    //  --- Utils
    //  -----------------------------------------------------------------------
    void copyParams(std::vector<mtu::TorchMergeTree<float>> &origins,
                    std::vector<mtu::TorchMergeTree<float>> &originsPrime,
                    std::vector<torch::Tensor> &vS,
                    std::vector<torch::Tensor> &vSPrime,
                    std::vector<mtu::TorchMergeTree<float>> &origins2,
                    std::vector<mtu::TorchMergeTree<float>> &origins2Prime,
                    std::vector<torch::Tensor> &vS2,
                    std::vector<torch::Tensor> &vS2Prime,
                    std::vector<std::vector<torch::Tensor>> &srcAlphas,
                    std::vector<std::vector<torch::Tensor>> &dstAlphas,
                    bool get);

    void copyParams(std::vector<std::vector<mtu::TorchMergeTree<float>>> &src,
                    std::vector<std::vector<mtu::TorchMergeTree<float>>> &dst);

    /**
     * @brief Set/Get custom parameters.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] get if the internal parameters should set from a copy (true)
     * or copied (false).
     */
    virtual void copyCustomParams(bool get) = 0;

    /**
     * @brief Construct a matrix as a torch tensor object containing all the
     * coordinates of the input topological representations in a specific layer.
     *
     * @param[in] alphas all the coordinates of the input topological
     * representations.
     * @param[in] indexes indexes of the topological representations to get.
     * @param[in] toGet boolean vector stating if a topological representation
     * should be processed or not.
     * @param[in] layerIndex index of the layer to process.
     * @param[in] alphasOut output torch tensor.
     */
    void getAlphasTensor(std::vector<std::vector<torch::Tensor>> &alphas,
                         std::vector<unsigned int> &indexes,
                         std::vector<bool> &toGet,
                         unsigned int layerIndex,
                         torch::Tensor &alphasOut);

    /**
     * @brief Construct a matrix as a torch tensor object containing all the
     * coordinates of the input topological representations in a specific layer.
     *
     * @param[in] alphas all the coordinates of the input topological
     * representations.
     * @param[in] indexes indexes of the topological representations to get.
     * @param[in] layerIndex index of the layer to process.
     * @param[in] alphasOut output torch tensor.
     */
    void getAlphasTensor(std::vector<std::vector<torch::Tensor>> &alphas,
                         std::vector<unsigned int> &indexes,
                         unsigned int layerIndex,
                         torch::Tensor &alphasOut);

    /**
     * @brief Construct a matrix as a torch tensor object containing all the
     * coordinates of the input topological representations in a specific layer.
     *
     * @param[in] alphas all the coordinates of the input topological
     * representations.
     * @param[in] layerIndex index of the layer to process.
     * @param[in] alphasOut output torch tensor.
     */
    void getAlphasTensor(std::vector<std::vector<torch::Tensor>> &alphas,
                         unsigned int layerIndex,
                         torch::Tensor &alphasOut);

    //  -----------------------------------------------------------------------
    //  --- Testing
    //  -----------------------------------------------------------------------
    void checkZeroGrad(unsigned int l, bool checkOutputBasis = true);

    bool isTreeHasBigValues(const ftm::MergeTree<float> &mTree,
                            float threshold = 10000);

    //  ---------------------------------------------------------------------------
    //  --- Main Functions
    //  ---------------------------------------------------------------------------
    /**
     * @brief Specific operations that can be done at the end of the
     * optimization. Like calling the "computeTrackingInformation" and the
     * "computeCorrelationMatrix" functions.
     *
     * This is a pure virtual function to define in derived classes.
     *
     * @param[in] trees the input trees.
     * @param[in] trees2 same as torchTrees but for second input, if any (i.e.
     * when join trees and split trees are given).
     */
    virtual void executeEndFunction(std::vector<ftm::MergeTree<float>> &trees,
                                    std::vector<ftm::MergeTree<float>> &trees2)
      = 0;
#endif

    void execute(std::vector<ftm::MergeTree<float>> &trees,
                 std::vector<ftm::MergeTree<float>> &trees2);
  }; // MergeTreeNeuralNetwork class

} // namespace ttk
