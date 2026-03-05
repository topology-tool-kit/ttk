/// \ingroup base
/// \class ttk::MergeTreeNeuralBase
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2024.
///
/// This module defines the %MergeTreeNeuralNetwork abstract class providing
/// functions to define a neural network processing merge trees or persistence
/// diagrams from end to end.
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
#include <MergeTreeAxesAlgorithmBase.h>
#include <MergeTreeTorchUtils.h>

#ifdef TTK_ENABLE_TORCH
#include <torch/torch.h>
#endif

namespace ttk {

  /**
   * This module defines the %MergeTreeNeuralNetwork abstract class providing
   * functions to define a neural network processing merge trees or persistence
   * diagrams from end to end.
   */
  class MergeTreeNeuralBase : virtual public Debug,
                              public MergeTreeAxesAlgorithmBase {

  protected:
    // ======== Model hyper-parameters
    // Dropout to use when training.
    double dropout_ = 0.0;
    // If the vectors should be initialized using euclidean distance (between
    // vectors representing the topological abstractions ordered given the
    // assignment to the barycenter), faster but less accurate than using
    // Wasserstein distance.
    bool euclideanVectorsInit_ = false;
    // If the vectors should be initialized randomly.
    bool randomAxesInit_ = false;
    // When computing the origin of the input basis, if the barycenter algorihm
    // should be initialized randomly (instead to the topological representation
    // minimizing the distance to the set), faster but less accurate.
    bool initBarycenterRandom_ = false;
    // When computing the origin of the input basis, if the barycenter algorithm
    // should run for only one iteration, faster but less accurate.
    bool initBarycenterOneIter_ = false;
    // If the structure of the origin of the output basis should be initialized
    // by copying the structure of the input basis.
    bool initOriginPrimeStructByCopy_ = true;
    // If the scalar values of the origin of the output basis should be
    // initialized by copying the values of the input basis.
    bool initOriginPrimeValuesByCopy_ = true;
    // Value between 0 and 1 allowing to add some randomness to the values of
    // the origin of the output basis when initOriginPrimeValuesByCopy_ is set
    // to true.
    double initOriginPrimeValuesByCopyRandomness_ = 0.0;
    // If activation functions should be used.
    bool activate_ = true;
    // Choice of the activation function
    // 0 : ReLU
    // 1 : Leaky ReLU
    unsigned int activationFunction_ = 1;

    bool useGpu_ = false;

    // ======== Testing
    float bigValuesThreshold_ = 0;

  public:
    MergeTreeNeuralBase();

#ifdef TTK_ENABLE_TORCH
    //  -----------------------------------------------------------------------
    //  --- Setter
    //  -----------------------------------------------------------------------
    void setDropout(const double dropout);

    void setEuclideanVectorsInit(const bool euclideanVectorsInit);

    void setRandomAxesInit(const bool randomAxesInit);

    void setInitBarycenterRandom(const bool initBarycenterRandom);

    void setInitBarycenterOneIter(const bool initBarycenterOneIter);

    void setInitOriginPrimeStructByCopy(const bool initOriginPrimeStructByCopy);

    void setInitOriginPrimeValuesByCopy(const bool initOriginPrimeValuesByCopy);

    void setInitOriginPrimeValuesByCopyRandomness(
      const double initOriginPrimeValuesByCopyRandomness);

    void setActivate(const bool activate);

    void setActivationFunction(const unsigned int activationFunction);

    void setUseGpu(const bool useGpu);

    void setBigValuesThreshold(const float bigValuesThreshold);

    //  -----------------------------------------------------------------------
    //  --- Utils
    //  -----------------------------------------------------------------------
    torch::Tensor activation(torch::Tensor &in);

    /**
     * @brief Fix the scalars of a merge tree to ensure that the nesting
     * condition is respected.
     *
     * @param[in] mTree Merge tree to process.
     */
    void fixTreePrecisionScalars(ftm::MergeTree<float> &mTree);

    /**
     * @brief Util function to print pairs of a merge tree.
     *
     * @param[in] mTree merge tree to process.
     * @param[in] useBD if the merge tree is in branch decomposition mode or
     * not.
     */
    void printPairs(const ftm::MergeTree<float> &mTree, bool useBD = true);

    //  -----------------------------------------------------------------------
    //  --- Distance
    //  -----------------------------------------------------------------------
    void getDistanceMatrix(const std::vector<mtu::TorchMergeTree<float>> &tmts,
                           std::vector<std::vector<float>> &distanceMatrix,
                           bool useDoubleInput = false,
                           bool isFirstInput = true);

    void getDistanceMatrix(const std::vector<mtu::TorchMergeTree<float>> &tmts,
                           const std::vector<mtu::TorchMergeTree<float>> &tmts2,
                           std::vector<std::vector<float>> &distanceMatrix);

    void getDifferentiableDistanceFromMatchings(
      const mtu::TorchMergeTree<float> &tree1,
      const mtu::TorchMergeTree<float> &tree2,
      const mtu::TorchMergeTree<float> &tree1_2,
      const mtu::TorchMergeTree<float> &tree2_2,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings,
      std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matchings2,
      torch::Tensor &tensorDist,
      bool doSqrt);

    void getDifferentiableDistance(const mtu::TorchMergeTree<float> &tree1,
                                   const mtu::TorchMergeTree<float> &tree2,
                                   const mtu::TorchMergeTree<float> &tree1_2,
                                   const mtu::TorchMergeTree<float> &tree2_2,
                                   torch::Tensor &tensorDist,
                                   bool isCalled,
                                   bool doSqrt);

    void getDifferentiableDistance(const mtu::TorchMergeTree<float> &tree1,
                                   const mtu::TorchMergeTree<float> &tree2,
                                   torch::Tensor &tensorDist,
                                   bool isCalled,
                                   bool doSqrt);

    void getDifferentiableDistanceMatrix(
      const std::vector<mtu::TorchMergeTree<float> *> &trees,
      const std::vector<mtu::TorchMergeTree<float> *> &trees2,
      std::vector<std::vector<torch::Tensor>> &outDistMat);
#endif
  }; // MergeTreeNeuralBase class

} // namespace ttk
