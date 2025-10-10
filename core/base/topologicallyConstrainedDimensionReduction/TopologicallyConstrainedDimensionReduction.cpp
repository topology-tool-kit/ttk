#include <TopologicallyConstrainedDimensionReduction.h>

#ifdef TTK_ENABLE_TORCH

using namespace torch::indexing;

ttk::TCDR::TCDR(bool useCUDA,
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
                bool inputIsImages)
  : NumberOfComponents(numberOfComponents), Epochs(epochs),
    LearningRate(learningRate), Optimizer(optimizer), Method(method),
    ModelType(modelType), InputIsImages(inputIsImages),
    Architecture(architecture), Activation(activation), BatchSize(batchSize),
    BatchNormalization(batchNormalization), RegCoefficient(regCoefficient) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TCDR");

  if(torch::cuda::is_available() && useCUDA && !deterministic)
    device = torch::kCUDA;

  if(deterministic) {
    at::manual_seed(seed);
    at::globalContext().setDeterministicFillUninitializedMemory(true);
    at::globalContext().setDeterministicAlgorithms(true, true);
  } else {
    at::globalContext().setDeterministicFillUninitializedMemory(false);
    at::globalContext().setDeterministicAlgorithms(false, false);
  }
}

int ttk::TCDR::execute(std::vector<std::vector<double>> &outputEmbedding,
                       const std::vector<double> &inputMatrix,
                       size_t n) {
  Timer tm{};
  printMsg("Initialization", 0., tm.getElapsedTime());

  const int inputSize = n;
  const int inputDimension = inputMatrix.size() / n;
  if(!InputIsImages)
    this->printMsg("input dimension: " + std::to_string(inputDimension), 0.0,
                   tm.getElapsedTime());
  else
    this->printMsg("input dimension: " + std::to_string(inputDimension) + " = "
                     + std::to_string((int)sqrt(inputDimension)) + " x "
                     + std::to_string((int)sqrt(inputDimension)) + " images",
                   .0, tm.getElapsedTime());
  this->printMsg("output dimension: " + std::to_string(NumberOfComponents), 0.0,
                 tm.getElapsedTime());
  this->printMsg(
    "#elements: " + std::to_string(inputSize), 0.0, tm.getElapsedTime());

#ifndef TTK_ENABLE_CGAL
  printWrn("TTK not compiled with CGAL enabled: this backend could be slow.");
#endif

  if((!InputIsImages && !AutoEncoder::isStringValid(Architecture))
     || (InputIsImages
         && !ConvolutionalAutoEncoder::isStringValid(Architecture))) {
    printErr("Invalid string for layers description.");
    return -1;
  }
  if(ModelType == MODEL::AUTOENCODER) {
    if(!InputIsImages)
      model = std::make_shared<AutoEncoder>(inputDimension, NumberOfComponents,
                                            Architecture, Activation,
                                            BatchNormalization);
    else
      model = std::make_shared<ConvolutionalAutoEncoder>(
        sqrt(inputDimension), NumberOfComponents, Architecture,
        BatchNormalization);
  } else if(ModelType == MODEL::AUTODECODER)
    model = std::make_shared<AutoDecoder>(inputDimension, inputSize,
                                          NumberOfComponents, Architecture,
                                          Activation, false);
  else if(ModelType == MODEL::DIRECT)
    model = std::make_shared<DirectOptimization>(inputSize, NumberOfComponents);
  model->to(device);

  if(Optimizer == OPTIMIZER::ADAM)
    torchOptimizer = std::make_shared<torch::optim::Adam>(
      model->parameters(), /*lr=*/LearningRate);
  else if(Optimizer == OPTIMIZER::SGD)
    torchOptimizer = std::make_shared<torch::optim::SGD>(
      model->parameters(), /*lr=*/LearningRate);
  else if(Optimizer == OPTIMIZER::LBFGS)
    torchOptimizer = std::make_shared<torch::optim::LBFGS>(
      model->parameters(), /*lr=*/LearningRate);

  const torch::Tensor input
    = torch::from_blob(const_cast<double *>(inputMatrix.data()),
                       {inputSize, inputDimension}, torch::kFloat64)
        .to(torch::kFloat32)
        .to(device);

  rpd::PointCloud points(inputSize, std::vector<double>(inputDimension));
  for(int i = 0; i < inputSize; ++i) {
    for(int j = 0; j < inputDimension; ++j)
      points[i][j] = inputMatrix[inputDimension * i + j];
  }

  if(Method == REGUL::NO_REGUL) {
    printMsg("Starting optimization", 0., tm.getElapsedTime());
    optimizeSimple(input);
  } else {
    printMsg("Computing input persistence", 0., tm.getElapsedTime());
    topologicalLossContainer
      = std::make_shared<TopologicalLoss>(input, points, Method);
    printMsg("Starting optimization", 0., tm.getElapsedTime());
    optimize(input);
  }

  const torch::Tensor latent = model->encode(input).cpu();
  for(int i = 0; i < inputSize; ++i) {
    for(int j = 0; j < NumberOfComponents; ++j)
      outputEmbedding[j][i] = latent[i][j].item<double>();
  }

  printMsg("Complete", 1., tm.getElapsedTime());
  return 0;
}

void ttk::TCDR::optimizeSimple(const torch::Tensor &input) const {
  int epoch = 0;

  auto closure = [&] {
    TensorIndex indices = Slice();
    if(BatchSize > 0)
      indices = torch::randint(input.size(0), {BatchSize}, torch::kInt).to(device);

    // step initialization
    torchOptimizer->zero_grad();
    const torch::Tensor prediction = model->forward(input.index(indices));

    // loss and optimizer step
    const torch::Tensor loss
      = torch::mse_loss(prediction, input.index(indices));
    loss.backward();

    // IO
    printLoss(epoch, loss.item<double>());

    return loss;
  };

  for(; epoch < Epochs; ++epoch)
    torchOptimizer->step(closure);
}

void ttk::TCDR::optimize(const torch::Tensor &input) const {
  int epoch = 0;

  auto closure = [&] {
    // step initialization
    torchOptimizer->zero_grad();
    const torch::Tensor latent = model->encode(input);
    const torch::Tensor prediction = model->decode(latent);

    // loss and optimizer step
    const torch::Tensor topologicalLoss
      = RegCoefficient * topologicalLossContainer->computeLoss(latent);
    const torch::Tensor reconstructionLoss = torch::mse_loss(prediction, input);
    const torch::Tensor loss = reconstructionLoss + topologicalLoss;
    loss.backward();

    // IO
    printLoss(epoch, loss.item<double>());

    return loss;
  };

  for(; epoch < Epochs; ++epoch)
    torchOptimizer->step(closure);
}

void ttk::TCDR::preOptimize(const torch::Tensor &input, const torch::Tensor &target) const {
  int epoch = 0;

  auto closure = [&] {
    // step initialization
    torchOptimizer->zero_grad();
    const torch::Tensor latent = model->encode(input);
    const torch::Tensor prediction = model->decode(latent);

    // loss and optimizer step
    const torch::Tensor loss = torch::mse_loss(latent, target) +
                               torch::mse_loss(prediction, input);
    loss.backward();

    // IO
    printLoss(epoch, loss.item<double>());

    return loss;
  };

  for(; epoch < Epochs; ++epoch)
    torchOptimizer->step(closure);
}

void ttk::TCDR::printLoss(int epoch, double loss) const {
  if(epoch % std::max(1, Epochs / 10) == 0)
    printMsg(
      "Loss at epoch " + std::to_string(epoch) + " : " + std::to_string(loss),
      double(epoch) / Epochs, -1, -1, debug::LineMode::REPLACE);
}

#endif