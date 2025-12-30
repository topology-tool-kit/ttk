/// \ingroup base
/// \class ttk::TopologicalDimensionReduction
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date 2024.
///
/// \brief TTK base class that embeds points into 2D, under topological
/// constraints
///
/// This module defines the %TopologicalDimensionReduction
/// class that serves as a backend for the DimensionReduction module. It embeds
/// high-dimensional point clouds into 2D, under topological constraints, with
/// an autoencoder-based approach.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/topoAEppTeaser/">Topological
///   Autoencoders++ Teaser example</a> \n
///
/// \b Related \b publications \n
///
/// "Topological Autoencoders" \n
/// Michael Moor, Max Horn, Bastian Rieck, Karsten Borgwardt, \n
/// Proceedings of the 37th International Conference on Machine Learning,
/// 2020. \n
///
/// "Optimizing persistent homology-based functions" \n
/// Mathieu Carriere, Frederic Chazal, Marc Glisse, Yuichi Ike,
/// Hariprasad Kannan, Yuhei Umeda, \n
/// Proceedings of the 38th International Conference on Machine Learning,
/// 2021. \n
///
/// "Topological Autoencoders++: Fast and Accurate Cycle-Aware Dimensionality
/// Reduction" \n
/// Mattéo Clémot, Julie Digne, Julien Tierny, \n
/// IEEE Transactions on Visualization and Computer Graphics.
/// Accepted, to be presented at IEEE VIS 2026.
///
/// \sa DimensionReduction.cpp %for a usage example.

#pragma once

// ttk includes
#include <Debug.h>
#include <DimensionReductionModel.h>
#include <TopologicalLoss.h>

namespace ttk {

  /**
   * The TopologicalDimensionReduction class provides a
   * backend for dimension reduction using autoencoders, with possible
   * constraints on the preservation of the topology of the input high
   * dimensional point cloud when projecting in low dimension
   */
  class TopologicalDimensionReduction : virtual public Debug {

  public:
    enum class OPTIMIZER : std::uint8_t {
      /** Adaptive Moment Estimation */
      ADAM = 0,
      /** Stochastic Gradient Descent */
      SGD = 1,
      /** Limited-memory Broyden–Fletcher–Goldfarb–Shanno */
      LBFGS = 2,
    };

    enum class MODEL : std::uint8_t {
      /** AutoEncoder architecture */
      AUTOENCODER = 0,
      /** AutoDecoder architecture */
      AUTODECODER = 1,
      /** Direct optimization*/
      DIRECT = 2,
    };

    using REGUL = TopologicalLoss::REGUL;

#ifdef TTK_ENABLE_TORCH

    TopologicalDimensionReduction(bool useCUDA,
                                  bool deterministic,
                                  int seed,
                                  int numberOfComponents,
                                  int epochs,
                                  double learningRate,
                                  OPTIMIZER optimizer,
                                  REGUL method,
                                  MODEL modelType,
                                  const std::string &architecture,
                                  const std::string &activation,
                                  int batchSize,
                                  bool batchNormalization,
                                  double regCoefficient,
                                  bool inputIsImages,
                                  bool preOptimize,
                                  int preOptimizeEpochs);

    /**
     * @brief Computes the projection with an AutoEncoder
     *
     * @param[out] outputEmbedding the final coordinates of the points
     *
     * @param[in] inputMatrix the high-dimension coordinates of the points
     *
     * @param[in] n the number of input points
     *
     * @return 0 in case of success.
     */
    int execute(std::vector<std::vector<double>> &outputEmbedding,
                const std::vector<double> &inputMatrix,
                size_t n);

  protected:
    const int NumberOfComponents;
    const int Epochs;
    const double LearningRate;
    const OPTIMIZER Optimizer;
    const REGUL Method;
    const MODEL ModelType;
    const bool InputIsImages;
    const std::string Architecture;
    const std::string Activation;
    const int BatchSize;
    const bool BatchNormalization;
    const double RegCoefficient;
    const bool PreOptimize;
    const int PreOptimizeEpochs;

  private:
    torch::DeviceType device{torch::kCPU};
    std::unique_ptr<DimensionReductionModel> model{nullptr};
    std::unique_ptr<torch::optim::Optimizer> torchOptimizer{nullptr};
    std::unique_ptr<TopologicalLoss> topologicalLossContainer{nullptr};

    int initializeModel(int inputSize, int inputDimension);
    void initializeOptimizer();

    void preOptimize(const torch::Tensor &input,
                     const torch::Tensor &target) const;

    void optimize(const torch::Tensor &input) const;
    void optimizeSimple(const torch::Tensor &input) const;

    inline void printLoss(int epoch, int maxEpoch, double loss) const;

#endif

  }; // TopologicalDimensionReduction class

} // namespace ttk
