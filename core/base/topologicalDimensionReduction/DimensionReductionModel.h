/// \ingroup base
/// \class ttk::DimensionReductionModel
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date July 2024.
///
/// \brief TTK base class for containing reduction dimension models
///
/// This file defines the %DimensionReductionModel abstract base class
/// for representing a dimension reduction model, and the following inherited
/// classes:
/// AutoEncoder
/// AutoDecoder
/// DirectOptimization
/// ConvolutionalAutoEncoder
///
/// \sa TopologicalDimensionReduction.cpp %for a usage example.

#pragma once

#ifdef TTK_ENABLE_TORCH

#include <torch/torch.h>

namespace ttk {

  /**
   * Abstract base class for dimension reduction models
   */
  class DimensionReductionModel : public torch::nn::Module {
  public:
    virtual torch::Tensor forward(torch::Tensor const &x) = 0;
    virtual torch::Tensor encode(torch::Tensor const &x) = 0;
    virtual torch::Tensor decode(torch::Tensor const &x) = 0;
  }; // DimensionReductionModel class

  /**
   * The AutoEncoder class provides a Torch-based autoencoder class for
   * autoencoder-based dimension reduction
   */
  class AutoEncoder : public DimensionReductionModel {
  public:
    AutoEncoder(int inputDim,
                int latentDim,
                const std::string &layersDescription,
                const std::string &activation = "ReLU",
                bool useBN = true);

    inline torch::Tensor forward(torch::Tensor const &x) override {
      return decoder->forward(encoder->forward(x));
    }

    inline torch::Tensor encode(torch::Tensor const &x) override {
      return encoder->forward(x);
    }

    inline torch::Tensor decode(torch::Tensor const &x) override {
      return decoder->forward(x);
    }

    static bool isStringValid(const std::string &s);

  private:
    torch::nn::Sequential encoder;
    torch::nn::Sequential decoder;
  }; // AutoEncoder class

  /**
   * The AutoDecoder class provides a Torch-based autodecoder class for
   * autodecoder-based dimension reduction
   */
  class AutoDecoder : public DimensionReductionModel {
  public:
    AutoDecoder(int inputDim,
                int inputSize,
                int latentDim,
                const std::string &layersDescription,
                const std::string &activation = "ReLU",
                bool useBN = true);

    inline torch::Tensor forward(torch::Tensor const & /*x*/) override {
      return decoder->forward(latent);
    }

    inline torch::Tensor encode(torch::Tensor const & /*x*/) override {
      return latent;
    }

    inline torch::Tensor decode(torch::Tensor const & /*x*/) override {
      return decoder->forward(latent);
    }

  private:
    torch::Tensor latent;
    torch::nn::Sequential decoder;
  }; // AutoDecoder class

  /**
   * The DirectOptimization class provides a Torch-based dummy class for
   * dimension reduction based on the direct optimization of the point cloud
   */
  class DirectOptimization : public DimensionReductionModel {
  public:
    DirectOptimization(int inputSize, int latentDim);

    inline torch::Tensor forward(torch::Tensor const &x) override {
      return x;
    }

    inline torch::Tensor encode(torch::Tensor const &x) override {
      input = x;
      return latent;
    }

    inline torch::Tensor decode(torch::Tensor const & /*x*/) override {
      return input;
    }

  private:
    torch::Tensor latent;
    torch::Tensor input;
  }; // DirectOptimization class

  /**
   * The ConvolutionalAutoEncoder class provides a Torch-based convolutional
   * autoencoder class for autoencoder-based dimension reduction of image
   * datasets
   */
  class ConvolutionalAutoEncoder : public DimensionReductionModel {
  public:
    ConvolutionalAutoEncoder(int imageSide,
                             int latentDim,
                             const std::string &layersDescription,
                             bool useBN);

    inline torch::Tensor forward(torch::Tensor const &x) override {
      return decoder->forward(encoder->forward(x));
    }

    inline torch::Tensor encode(torch::Tensor const &x) override {
      return encoder->forward(x);
    }

    inline torch::Tensor decode(torch::Tensor const &x) override {
      return decoder->forward(x);
    }

    static bool isStringValid(const std::string &s);

  private:
    torch::nn::Sequential encoder;
    torch::nn::Sequential decoder;
  }; // ConvolutionalAutoEncoder class

} // namespace ttk

#endif
