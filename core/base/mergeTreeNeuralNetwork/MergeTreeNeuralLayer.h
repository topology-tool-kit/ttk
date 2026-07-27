/// \ingroup base
/// \class ttk::MergeTreeNeuralLayer
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// This module defines the %MergeTreeNeuralLayer class that provide methods to
/// define and use a wasserstein layer able to process merge trees or
/// persistence diagrams.
///
/// To initialize the layer you can use the following functions:
/// "initInputBasisOrigin" and "initInputBasisVectors" for the input basis, then
/// "initOutputBasis" for the output basis. Please refer to the documentation of
/// these functions for how to use them.
///
/// Then, you should call the "requires_grad" function, with a parameter sets to
/// true, to enable torch to compute gradient for this layer for it to be
/// optimized.
///
/// The layer can then be used with the "forward" function to pass a topological
/// representation as input and get the topological representation at the output
/// of the layer (the transformed input).
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
#include <MergeTreeTorchUtils.h>

#ifdef TTK_ENABLE_TORCH
#include <torch/torch.h>
#endif

namespace ttk {

  /**
   * The MergeTreeNeuralLayer class provides methods to define and use a
   * wasserstein layer able to process merge trees or persistence diagrams.
   */
  class MergeTreeNeuralLayer : virtual public Debug,
                               public MergeTreeNeuralBase {

#ifdef TTK_ENABLE_TORCH
    // Layer parameters
    torch::Tensor vSTensor_, vSPrimeTensor_, vS2Tensor_, vS2PrimeTensor_;
    mtu::TorchMergeTree<float> origin_, originPrime_, origin2_, origin2Prime_;
#endif

  public:
    MergeTreeNeuralLayer();

#ifdef TTK_ENABLE_TORCH
    //  -----------------------------------------------------------------------
    //  --- Getter/Setter
    //  -----------------------------------------------------------------------
    const mtu::TorchMergeTree<float> &getOrigin() const;

    const mtu::TorchMergeTree<float> &getOriginPrime() const;

    const mtu::TorchMergeTree<float> &getOrigin2() const;

    const mtu::TorchMergeTree<float> &getOrigin2Prime() const;

    const torch::Tensor &getVSTensor() const;

    const torch::Tensor &getVSPrimeTensor() const;

    const torch::Tensor &getVS2Tensor() const;

    const torch::Tensor &getVS2PrimeTensor() const;

    void setOrigin(const mtu::TorchMergeTree<float> &tmt);

    void setOriginPrime(const mtu::TorchMergeTree<float> &tmt);

    void setOrigin2(const mtu::TorchMergeTree<float> &tmt);

    void setOrigin2Prime(const mtu::TorchMergeTree<float> &tmt);

    void setVSTensor(const torch::Tensor &vS);

    void setVSPrimeTensor(const torch::Tensor &vS);

    void setVS2Tensor(const torch::Tensor &vS);

    void setVS2PrimeTensor(const torch::Tensor &vS);

    //  -----------------------------------------------------------------------
    //  --- Init
    //  -----------------------------------------------------------------------
    /**
     * @brief Initialize the tree structure of the origin in an output basis
     * whose scalars have already been initialized.
     *
     * @param[out] originPrime origin merge tree of an output basis.
     * @param[in] isJT if the tree is a join tree.
     * @param[in] baseOrigin a merge tree whose tree structure can be used to
     * initialize the tree structure of originPrime (typically, it can be the
     * origin of the input basis).
     */
    void initOutputBasisTreeStructure(mtu::TorchMergeTree<float> &originPrime,
                                      bool isJT,
                                      mtu::TorchMergeTree<float> &baseOrigin);

    /**
     * @brief Initialize the output basis.
     *
     * @param[in] dim the number of nodes in the origin of the output basis
     * (corresponds to twice the number of persistence pairs).
     * @param[in] dim2 same as dim but for second input, if any (i.e. when join
     * trees and split trees are given).
     * @param[in] baseTensor the scalars of a merge tree that can be used to
     * initialize the scalars of the origin in the outuput basis (typically, it
     * can be the scalars of the origin of the previous output basis or the
     * origin of the input basis of this layer if this is the first one).
     */
    void initOutputBasis(const unsigned int dim,
                         const unsigned int dim2,
                         const torch::Tensor &baseTensor);

    /**
     * @brief Initialize the axes of the output basis.
     *
     * @param[in] w a matrix that will be used to compute the axes of the output
     * basis, if B is a matrix corresponding to the axes of the input basis then
     * the axes of the output basis will be initialized as wB.
     * @param[in] w2 same as w but for second input, if any (i.e. when join
     * trees and split trees are given).
     */
    void initOutputBasisVectors(torch::Tensor &w, torch::Tensor &w2);

    /**
     * @brief Initialize the axes of the output basis.
     *
     * @param[in] dim the number of nodes in the origin of the output basis
     * (corresponds to twice the number of persistence pairs).
     * @param[in] dim2 same as dim but for second input, if any (i.e. when join
     * trees and split trees are given).
     */
    void initOutputBasisVectors(unsigned int dim, unsigned int dim2);

    /**
     * @brief Initialize the origin of the input basis as the barycenter of an
     * ensemble of trees.
     *
     * @param[in] treesToUse the trees to use for the initialization (typically,
     * the input trees of this layer).
     * @param[in] trees2ToUse same as treesToUse but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] barycenterSizeLimitPercent the maximum number of nodes allowed
     * for the barycenter as a percentage of the total number of nodes in the
     * input trees (0 for no effect).
     * @param[in] barycenterMaxNoPairs the maximum number of nodes in the
     * barycenter (0 for no effect).
     * @param[in] barycenterMaxNoPairs2 same as barycenterMaxNoPairs but for
     * second input, if any (i.e. when join trees and split trees are given).
     * @param[out] inputToBaryDistances the distances of the input trees to the
     * origin of the basis.
     * @param[out] baryMatchings the matchings between the input trees and the
     * origin of the basis.
     * @param[out] baryMatchings2 same as baryMatchings but for second input, if
     * any (i.e. when join trees and split trees are given).
     */
    void initInputBasisOrigin(
      std::vector<ftm::MergeTree<float>> &treesToUse,
      std::vector<ftm::MergeTree<float>> &trees2ToUse,
      double barycenterSizeLimitPercent,
      unsigned int barycenterMaxNoPairs,
      unsigned int barycenterMaxNoPairs2,
      std::vector<double> &inputToBaryDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2);

    /**
     * @brief Initialize the axes of the input basis.
     *
     * @param[in] tmTrees the trees to use for the initialization as
     * TorchMergeTree objects (typically, the input trees of this layer).
     * @param[in] tmTrees2 same as tmTrees but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] trees the trees to use for the initialization as MergeTree
     * objects (typically, the input trees of this layer).
     * @param[in] trees2 same as tmTrees but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] noVectors number of axes in the basis.
     * @param[out] allAlphasInit the coordinates of each input tree in the
     * basis.
     * @param[in] inputToBaryDistances the distances of the input trees to the
     * origin of the basis.
     * @param[in] baryMatchings the matchings between the input trees and the
     * origin of the basis.
     * @param[in] baryMatchings2 same as baryMatchings but for second input, if
     * any (i.e. when join trees and split trees are given).
     * @param[in] origin the origin of the basis.
     * @param[in] origin2 same as origin but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[out] vSTensor the tensor representing the basis.
     * @param[out] vS2Tensor same as vS2Tensor but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] useInputBasis this boolean allows this function to also be
     * used for initializing the output basis (by setting this parameter to
     * false).
     */
    void initInputBasisVectors(
      std::vector<mtu::TorchMergeTree<float>> &tmTrees,
      std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
      std::vector<ftm::MergeTree<float>> &trees,
      std::vector<ftm::MergeTree<float>> &trees2,
      unsigned int noVectors,
      std::vector<torch::Tensor> &allAlphasInit,
      std::vector<double> &inputToBaryDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2,
      mtu::TorchMergeTree<float> &origin,
      mtu::TorchMergeTree<float> &origin2,
      torch::Tensor &vSTensor,
      torch::Tensor &vS2Tensor,
      bool useInputBasis = true);

    /**
     * @brief Overloaded function that initialize the axes of the input basis of
     * this instantiated MergeTreeNeuralLayer object.
     *
     * @param[in] tmTrees the trees to use for the initialization as
     * TorchMergeTree objects (typically, the input trees of this layer).
     * @param[in] tmTrees2 same as tmTrees but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] trees the trees to use for the initialization as MergeTree
     * objects (typically, the input trees of this layer).
     * @param[in] trees2 same as tmTrees but for second input, if any
     * (i.e. when join trees and split trees are given).
     * @param[in] noVectors number of axes in the basis.
     * @param[out] allAlphasInit the coordinates of each input tree in the
     * basis.
     * @param[in] inputToBaryDistances the distances of the input trees to the
     * origin of the basis.
     * @param[in] baryMatchings the matchings between the input trees and the
     * origin of the basis.
     * @param[out] baryMatchings2 same as baryMatchings but for second input, if
     * any (i.e. when join trees and split trees are given).
     * @param[in] useInputBasis this boolean allows this function to also be
     * used for initializing the output basis (by setting this parameter to
     * false).
     */
    void initInputBasisVectors(
      std::vector<mtu::TorchMergeTree<float>> &tmTrees,
      std::vector<mtu::TorchMergeTree<float>> &tmTrees2,
      std::vector<ftm::MergeTree<float>> &trees,
      std::vector<ftm::MergeTree<float>> &trees2,
      unsigned int noVectors,
      std::vector<torch::Tensor> &allAlphasInit,
      std::vector<double> &inputToBaryDistances,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings,
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode, double>>>
        &baryMatchings2,
      bool useInputBasis = true);

    void requires_grad(const bool requireGrad);

    void cuda();

    //  -----------------------------------------------------------------------
    //  --- Interpolation
    //  -----------------------------------------------------------------------
    /**
     * @brief Projection ensuring that no pairs are below diagonal.
     * Warning: this function only updates the Tensor object and not the scalars
     * of the MergeTree object, for this please call interpolationProjection.
     *
     * @param[in,out] interpolationTensor merge tree to process.
     */
    void interpolationDiagonalProjection(
      mtu::TorchMergeTree<float> &interpolationTensor);

    /**
     * @brief Projection ensuring the nesting condition.
     * Warning: this function only updates the Tensor object and not the scalars
     * of the MergeTree object, for this please call interpolationProjection.
     *
     * @param[in,out] interpolation merge tree to process.
     */
    void
      interpolationNestingProjection(mtu::TorchMergeTree<float> &interpolation);

    /**
     * @brief Projection ensuring the elder rule (no pairs below diagonal and
     * nesting condition), updates the Tensor object AND the scalars of the
     * MergeTree object.
     *
     * @param[in,out] interpolation merge tree to process.
     */
    void interpolationProjection(mtu::TorchMergeTree<float> &interpolation);

    /**
     * @brief Creates the merge tree at coordinates alphas of the basis with
     * origin as origin and vS as axes, followed by a projection ensuring the
     * elder rule.
     *
     * @param[in] origin origin of the basis.
     * @param[in] vS axes of the basis.
     * @param[in] alphas coordinates on the basis to evaluate.
     * @param[out] interpolation output merge tree.
     */
    void getMultiInterpolation(const mtu::TorchMergeTree<float> &origin,
                               const torch::Tensor &vS,
                               torch::Tensor &alphas,
                               mtu::TorchMergeTree<float> &interpolation);

    //  -----------------------------------------------------------------------
    //  --- Forward
    //  -----------------------------------------------------------------------
    /**
     * @brief Computes the necessary tensors used in the computation of the best
     * coordinates in the input basis of the input merge tree at constant
     * assignment.
     * According Appendix B of the reference "Wasserstein Auto-Encoders of Merge
     * Trees (and Persistence Diagrams), reorderedTreeTensor corresponds to
     * Beta_1', deltaOrigin to Beta_3', deltaA to B_2', originTensor_f to O' and
     * vSTensor_f to (B(O'))'.
     *
     * @param[in] tree input merge tree.
     * @param[in] origin origin of the basis.
     * @param[in] vSTensor axes of the basis.
     * @param[in] interpolated current estimation of the input tree in the
     * basis.
     * @param[in] matching matching between the input tree and the current
     * estimation in the basis.
     * @param[out] reorderedTreeTensor reordered tensor of the input tree given
     * the matching to its estimation on the basis (with zero on non-matched
     * pairs).
     * @param[out] deltaOrigin tensor of the projected pairs on the
     * diagonal of the origin in the input tree.
     * @param[out] deltaA tensor corresponding to the linear combination of the
     * axes of the basis given the coordinates for the projected pairs on the
     * diagonal of the estimated tree in the input tree.
     * @param[out] originTensor_f tensor of the origin.
     * @param[out] vSTensor_f tensor of the basis.
     */
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

    /**
     * @brief Computes the best coordinates in the input basis of the input
     * merge tree at constant assignment.
     *
     * @param[in] tree input merge tree.
     * @param[in] origin origin of the basis.
     * @param[in] vSTensor axes of the basis.
     * @param[in] interpolated current estimation of the input tree in the
     * basis.
     * @param[in] matching matching between the input tree and the current
     * estimation in the basis.
     * @param[in] tree2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] origin2 same as origin but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] vSTensor2 same as vSTensor but for second input, if any (i.e.
     * when join trees and split trees are given).
     * @param[in] interpolated2 same as interpolation but for second input, if
     * any (i.e. when join trees and split trees are given).
     * @param[in] matching2 same as matching but for second input, if any (i.e.
     * when join trees and split trees are given).
     * @param[out] alphasOut best coordinates of the input tree in the basis at
     * constant assignment.
     */
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

    /**
     * @brief Estimates the coordinates in the input basis of the input merge
     * tree.
     *
     * @param[in] tree input merge tree.
     * @param[in] tree2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] alphasInit the initial coordinates to use when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] bestMatching matching between the input tree and the tree at
     * the best estimation of the coordinates in the basis.
     * @param[in] bestMatching2 same as bestMatching but for second input, if
     * any (i.e. when join trees and split trees are given).
     * @param[out] bestAlphas the best estimation of the coordinates in the
     * input basis of the input merge tree.
     * @param[in] isCalled true if this function is called from a parallalized
     * context, i.e. if a team of threads has already been created and that
     * therefore it is not needed to create one.
     * @param[in] useInputBasis this boolean allows this function to also be
     * used for the output basis (by setting this parameter to false).
     *
     * @return the distance between the input merge tree and its best estimated
     * projection in the input basis.
     */
    float assignmentOneData(
      mtu::TorchMergeTree<float> &tree,
      mtu::TorchMergeTree<float> &tree2,
      unsigned int k,
      torch::Tensor &alphasInit,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &bestMatching2,
      torch::Tensor &bestAlphas,
      bool isCalled = false,
      bool useInputBasis = true);

    /**
     * @brief Estimates the coordinates in the input basis of the input merge
     * tree.
     *
     * @param[in] tree input merge tree.
     * @param[in] tree2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] alphasInit the initial coordinates to use when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[out] bestAlphas the best estimation of the coordinates in the
     * input basis of the input merge tree.
     * @param[in] isCalled true if this function is called from a parallalized
     * context, i.e. if a team of threads has already been created and that
     * therefore it is not needed to create one.
     * @param[in] useInputBasis this boolean allows this function to also be
     * used for the output basis (by setting this parameter to false).
     *
     * @return the distance between the input merge tree and its best estimated
     * projection in the input basis.
     */
    float assignmentOneData(mtu::TorchMergeTree<float> &tree,
                            mtu::TorchMergeTree<float> &tree2,
                            unsigned int k,
                            torch::Tensor &alphasInit,
                            torch::Tensor &bestAlphas,
                            bool isCalled = false,
                            bool useInputBasis = true);

    /**
     * @brief Reconstruct an ouput merge tree given coordinates.
     *
     * @param[in] alphas coordinates to use in the output basis.
     * @param[out] out the output merge tree.
     * @param[out] out2 same as out but for second input, if any (i.e. when join
     * trees and split trees are given).
     * @param[in] activate true if activation function should be used, false
     * otherwise.
     * @param[in] train true if the input merge tree is in the training set
     * (false if validation/testing set).
     */
    void outputBasisReconstruction(torch::Tensor &alphas,
                                   mtu::TorchMergeTree<float> &out,
                                   mtu::TorchMergeTree<float> &out2,
                                   bool activate = true,
                                   bool train = false);

    /**
     * @brief Pass a merge tree through the layer and get the output.
     *
     * @param[in] tree input merge tree.
     * @param[in] tree2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] alphasInit the initial coordinates to use when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[out] out the output merge tree.
     * @param[out] out2 same as out but for second input, if any (i.e. when join
     * trees and split trees are given).
     * @param[out] bestAlphas the best estimation of the coordinates in the
     * input basis of the input merge tree.
     * @param[out] bestDistance the distance between the input merge tree and
     * its best estimated projection in the input basis.
     * @param[in] train true if the input merge tree is in the training set
     * (false if validation/testing set).
     *
     * @return true if the output merge tree has no nodes.
     */
    bool forward(mtu::TorchMergeTree<float> &tree,
                 mtu::TorchMergeTree<float> &tree2,
                 unsigned int k,
                 torch::Tensor &alphasInit,
                 mtu::TorchMergeTree<float> &out,
                 mtu::TorchMergeTree<float> &out2,
                 torch::Tensor &bestAlphas,
                 float &bestDistance,
                 bool train = false);

    /**
     * @brief Pass a merge tree through the layer and get the output.
     *
     * @param[in] tree input merge tree.
     * @param[in] tree2 same as tree but for second input, if any (i.e. when
     * join trees and split trees are given).
     * @param[in] k number of projection steps to do when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[in] alphasInit the initial coordinates to use when estimating the
     * coordinates in the input basis of the input merge tree.
     * @param[out] out the output merge tree.
     * @param[out] out2 same as out but for second input, if any (i.e. when join
     * trees and split trees are given).
     * @param[out] bestAlphas the best estimation of the coordinates in the
     * input basis of the input merge tree.
     * @param[in] train true if the input merge tree is in the training set
     * (false if validation/testing set).
     *
     * @return true if the output merge tree has no nodes.
     */
    bool forward(mtu::TorchMergeTree<float> &tree,
                 mtu::TorchMergeTree<float> &tree2,
                 unsigned int k,
                 torch::Tensor &alphasInit,
                 mtu::TorchMergeTree<float> &out,
                 mtu::TorchMergeTree<float> &out2,
                 torch::Tensor &bestAlphas,
                 bool train = false);

    //  -----------------------------------------------------------------------
    //  --- Projection
    //  -----------------------------------------------------------------------
    /**
     * @brief Projection that ensures that the origins of the input and output
     * bases respect the elder rule.
     */
    void projectionStep();

    //  -----------------------------------------------------------------------
    //  --- Utils
    //  -----------------------------------------------------------------------
    void copyParams(mtu::TorchMergeTree<float> &origin,
                    mtu::TorchMergeTree<float> &originPrime,
                    torch::Tensor &vS,
                    torch::Tensor &vSPrime,
                    mtu::TorchMergeTree<float> &origin2,
                    mtu::TorchMergeTree<float> &origin2Prime,
                    torch::Tensor &vS2,
                    torch::Tensor &vS2Prime,
                    bool get);

    /**
     * @brief Fix the scalars of a merge tree to ensure that the nesting
     * condition is respected.
     *
     * @param[in] scalarsVector scalars array to process.
     * @param[in] node node to adjust.
     * @param[in] refNode reference node.
     */
    void adjustNestingScalars(std::vector<float> &scalarsVector,
                              ftm::idNode node,
                              ftm::idNode refNode);

    /**
     * @brief Create a balanced BDT structure (for output basis initialization).
     *
     * @param[in] parents vector containing the possible parents for each node.
     * @param[in] children vector containing the possible children for each
     * node.
     * @param[in] scalarsVector vector containing the scalars value.
     * @param[out] childrenFinal output vector containing the children of each
     * node, representing the tree structure.
     */
    void
      createBalancedBDT(std::vector<std::vector<ftm::idNode>> &parents,
                        std::vector<std::vector<ftm::idNode>> &children,
                        std::vector<float> &scalarsVector,
                        std::vector<std::vector<ftm::idNode>> &childrenFinal);

    //  -----------------------------------------------------------------------
    //  --- Testing
    //  -----------------------------------------------------------------------
    bool isTreeHasBigValues(ftm::MergeTree<float> &mTree,
                            float threshold = 10000);
#endif
  }; // MergeTreeNeuralLayer class

} // namespace ttk
