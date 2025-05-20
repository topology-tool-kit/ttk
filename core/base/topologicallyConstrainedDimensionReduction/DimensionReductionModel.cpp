#include "DimensionReductionModel.h"

#ifdef TTK_ENABLE_TORCH

ttk::AutoEncoder::AutoEncoder(int inputDim, int latentDim, const std::string &layersDescription,
                              const std::string &activation, bool useBN) {

  std::istringstream iss (layersDescription);
  const std::vector<std::string> hiddenDimsParsed(std::istream_iterator<std::string>{iss},
                                                  std::istream_iterator<std::string>());
  std::vector<unsigned> dims (1, inputDim);
  for (const std::string &s : hiddenDimsParsed)
    dims.push_back(std::stoi(s));
  dims.push_back(latentDim);

  const int n = dims.size() - 1;
  encoder = torch::nn::Sequential(torch::nn::Linear(dims[0],dims[1]));
  decoder = torch::nn::Sequential(torch::nn::Linear(dims[n],dims[n-1]));
  for (unsigned i=1; i<dims.size()-1; ++i) {
    if (activation == "ReLU") {
      encoder->push_back(torch::nn::ReLU());
      decoder->push_back(torch::nn::ReLU());
    } else if (activation == "Tanh") {
      encoder->push_back(torch::nn::Tanh());
      decoder->push_back(torch::nn::Tanh());
    }
    if (useBN) {
      encoder->push_back(torch::nn::BatchNorm1d(dims[i]));
      decoder->push_back(torch::nn::BatchNorm1d(dims[n-i]));
    }
    encoder->push_back(torch::nn::Linear(dims[i], dims[i+1]));
    decoder->push_back(torch::nn::Linear(dims[n-i], dims[n-(i+1)]));
  }
  register_module("encoder", encoder);
  register_module("decoder", decoder);
}

bool ttk::AutoEncoder::isStringValid(const std::string& s) {
  return std::regex_match(s, std::regex("([0-9]+( )*)*"));
}

ttk::AutoDecoder::AutoDecoder(int inputDim, int inputSize, int latentDim, const std::string &layersDescription,
                              const std::string &activation, bool useBN) {

  std::istringstream iss (layersDescription);
  const std::vector<std::string> hiddenDimsParsed(std::istream_iterator<std::string>{iss},
                                                  std::istream_iterator<std::string>());
  std::vector<unsigned> dims (1, inputDim);
  for (const std::string &s : hiddenDimsParsed)
    dims.push_back(std::stoi(s));
  dims.push_back(latentDim);

  latent = torch::rand({inputSize, latentDim});

  const int n = dims.size() - 1;
  decoder = torch::nn::Sequential(torch::nn::Linear(dims[n],dims[n-1]));
  for (unsigned i=1; i<dims.size()-1; ++i) {
    if (activation == "ReLU")
      decoder->push_back(torch::nn::ReLU());
    else if (activation == "Tanh")
      decoder->push_back(torch::nn::Tanh());
    if (useBN)
      decoder->push_back(torch::nn::BatchNorm1d(dims[n-i]));
    decoder->push_back(torch::nn::Linear(dims[n-i], dims[n-(i+1)]));
  }
  register_parameter("latent", latent, true);
  register_module("decoder", decoder);
}

ttk::DirectOptimization::DirectOptimization(int inputSize, int latentDim) {
  latent = torch::rand({inputSize, latentDim});
  register_parameter("latent", latent, true);
}

ttk::ConvolutionalAutoEncoder::ConvolutionalAutoEncoder(int imageSide, int latentDim, const std::string &layersDescription, bool useBN) {
  std::istringstream iss (layersDescription);
  const std::vector<std::string> hiddenDimsParsed(std::istream_iterator<std::string>{iss},
                                                  std::istream_iterator<std::string>());
  std::vector<int> denseLayersSizes = {-1};
  std::vector<int> convolutionalLayersChannels = {1};
  std::vector<int> convolutionalLayersStrides;
  for (const std::string &s : hiddenDimsParsed) {
    if (std::regex_match(s, std::regex("c[0-9]+/[0-9]+"))) {
      convolutionalLayersChannels.push_back(std::stoi(s.substr(1, s.find('/'))));
      convolutionalLayersStrides.push_back(std::stoi(s.substr(s.find('/')+1, s.length())));
    }
    else
      denseLayersSizes.push_back(std::stoi(s));
  }
  denseLayersSizes.push_back(latentDim);
  int convolutionalOutputImageSide = imageSide;

  /*** convolutional encoder ***/
  /** unflattening for convolutional layers **/
  encoder->push_back(torch::nn::Unflatten(torch::nn::UnflattenOptions(1, {1, imageSide, imageSide})));
  /** convolutional layers **/
  for (unsigned c = 0; c<convolutionalLayersChannels.size()-1; ++c) {
    encoder->push_back(torch::nn::Conv2d(torch::nn::Conv2dOptions(convolutionalLayersChannels[c],convolutionalLayersChannels[c+1],3).padding(1).stride(convolutionalLayersStrides[c])));
    encoder->push_back(torch::nn::ReLU());
    convolutionalOutputImageSide /= convolutionalLayersStrides[c];
  }
  /** flattening for dense layers **/
  encoder->push_back(torch::nn::Flatten());

  /*** dense encoder / decoder ***/
  denseLayersSizes[0] = convolutionalLayersChannels[convolutionalLayersChannels.size()-1] * convolutionalOutputImageSide * convolutionalOutputImageSide;
  const int n = denseLayersSizes.size() - 1;
  encoder->push_back(torch::nn::Linear(denseLayersSizes[0],denseLayersSizes[1]));
  decoder->push_back(torch::nn::Linear(denseLayersSizes[n],denseLayersSizes[n-1]));
  for (unsigned i=1; i<denseLayersSizes.size()-1; ++i) {
    encoder->push_back(torch::nn::ReLU());
    decoder->push_back(torch::nn::ReLU());
    if (useBN) {
      encoder->push_back(torch::nn::BatchNorm1d(denseLayersSizes[i]));
      decoder->push_back(torch::nn::BatchNorm1d(denseLayersSizes[n-i]));
    }
    encoder->push_back(torch::nn::Linear(denseLayersSizes[i], denseLayersSizes[i+1]));
    decoder->push_back(torch::nn::Linear(denseLayersSizes[n-i], denseLayersSizes[n-(i+1)]));
  }

  /*** convolutional decoder ***/
  /** unflattening for convolutional layers **/
  decoder->push_back(torch::nn::Unflatten(torch::nn::UnflattenOptions(1, {convolutionalLayersChannels[convolutionalLayersChannels.size()-1], convolutionalOutputImageSide, convolutionalOutputImageSide})));
  /** convolutional layers **/
  for (unsigned c = convolutionalLayersChannels.size()-1; c>0; --c) {
    decoder->push_back(torch::nn::ReLU());
    decoder->push_back(torch::nn::ConvTranspose2d(torch::nn::ConvTranspose2dOptions(convolutionalLayersChannels[c],convolutionalLayersChannels[c-1],3).padding(1).stride(convolutionalLayersStrides[c-1]).output_padding(1)));
  }
  /** flattening for output **/
  decoder->push_back(torch::nn::Flatten());

  register_module("encoder", encoder);
  register_module("decoder", decoder);
}

bool ttk::ConvolutionalAutoEncoder::isStringValid(const std::string& s) {
  return std::regex_match(s, std::regex("(c[0-9]+/[0-9]+( )*)*([0-9]+( )*)*"));
}

#endif