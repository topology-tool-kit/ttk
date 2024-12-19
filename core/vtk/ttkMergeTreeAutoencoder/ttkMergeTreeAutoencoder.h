/// \ingroup vtk
/// \class ttkMergeTreeAutoencoder
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeAutoencoder module.
///
/// This VTK filter uses the ttk::MergeTreeAutoencoder module to compute an
/// auto-encoder of merge trees or persistence diagrams.
///
/// \param Input vtkMultiBlockDataSet Input trees
/// \param Input (optional) vtkMultiBlockDataSet Input trees. If input are merge
/// trees, then this input can be used to process join and split trees together.
/// Pass as input either join or split trees in the first input and the other
/// type of trees in the second input. If input are persistence diagrams, then
/// this has no effect to use this input.
/// \param Input (optional) vtkTable Info (such as clustering assigment)
/// \param Output vtkMultiBlockDataSet Origins
/// \param Output vtkMultiBlockDataSet Bases Axes
/// \param Output vtkMultiBlockDataSet Coefficients
/// \param Output vtkMultiBlockDataSet Processed Input Trees
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication: \n
/// "Wasserstein Auto-Encoders of Merge Trees (and Persistence Diagrams)" \n
/// Mathieu Pont,  Julien Tierny.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2023
///
/// \sa ttk::MergeTreeAutoencoder
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeWAE/">Merge
///   tree Wasserstein Auto-Encoder example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramWAE/">Persistence
///   Diagram Wasserstein Auto-Encoder example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeAutoencoderModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <FTMTree.h>
#include <MergeTreeAutoencoder.h>

class TTKMERGETREEAUTOENCODER_EXPORT ttkMergeTreeAutoencoder
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeAutoencoder // and we inherit from the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Input options
  double oldEpsilonTree1;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> intermediateDTrees;
  std::vector<vtkUnstructuredGrid *> treesNodes, treesNodes2;
  std::vector<vtkUnstructuredGrid *> treesArcs, treesArcs2;
  std::vector<vtkDataSet *> treesSegmentation, treesSegmentation2;

  void setDataVisualization(int ttkNotUsed(numInputs),
                            int ttkNotUsed(numInputs2)) {
  }

  void resetDataVisualization() {
    setDataVisualization(0, 0);
    treesSegmentation.clear();
    treesSegmentation2.clear();
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
  using vtkAlgorithm::SetInputArrayToProcess;
  void SetInputArrayToProcess(const char *name) {
    vtkAlgorithm::SetInputArrayToProcess(0, 2, 0, 6, name);
  }

  void SetDoCompute(bool doCompute) {
    doCompute_ = doCompute;
    Modified();
    resetDataVisualization();
  }
  bool GetDoCompute() {
    return doCompute_;
  }

  void SetNormalizedWasserstein(bool nW) {
    normalizedWasserstein_ = nW;
    Modified();
    resetDataVisualization();
  }
  bool GetNormalizedWasserstein() {
    return normalizedWasserstein_;
  }

  void SetNumberOfEncoderLayers(unsigned int numberOfEncoderLayers) {
    encoderNoLayers_ = numberOfEncoderLayers;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfEncoderLayers() {
    return encoderNoLayers_;
  }

  void SetScaleLayerAfterLatent(unsigned int scaleLayerAfterLatent) {
    scaleLayerAfterLatent_ = scaleLayerAfterLatent;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetScaleLayerAfterLatent() {
    return scaleLayerAfterLatent_;
  }

  void SetInputNumberOfAxes(unsigned int numberOfAxes) {
    inputNumberOfAxes_ = numberOfAxes;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetInputNumberOfAxes() {
    return inputNumberOfAxes_;
  }

  void SetInputOriginPrimeSizePercent(double originSize) {
    inputOriginPrimeSizePercent_ = originSize;
    Modified();
    resetDataVisualization();
  }
  double GetInputOriginPrimeSizePercent() {
    return inputOriginPrimeSizePercent_;
  }

  //  Latent space number of axes
  void SetNumberOfAxes(unsigned int numberOfAxes) {
    numberOfAxes_ = numberOfAxes;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfAxes() {
    return numberOfAxes_;
  }

  void SetLatentSpaceOriginPrimeSizePercent(double originSize) {
    latentSpaceOriginPrimeSizePercent_ = originSize;
    Modified();
    resetDataVisualization();
  }
  double GetLatentSpaceOriginPrimeSizePercent() {
    return latentSpaceOriginPrimeSizePercent_;
  }

  void SetNumberOfProjectionSteps(unsigned int noSteps) {
    k_ = noSteps;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfProjectionSteps() {
    return k_;
  }

  void SetBarycenterSizeLimitPercent(double percent) {
    barycenterSizeLimitPercent_ = percent;
    Modified();
    resetDataVisualization();
  }
  double GetBarycenterSizeLimitPercent() {
    return barycenterSizeLimitPercent_;
  }

  void SetMinIteration(unsigned int minIteration) {
    minIteration_ = minIteration;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetMinIteration() {
    return minIteration_;
  }

  void SetMaxIteration(unsigned int maxIteration) {
    maxIteration_ = maxIteration;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetMaxIteration() {
    return maxIteration_;
  }

  void SetIterationGap(unsigned int iterationGap) {
    iterationGap_ = iterationGap;
    Modified();
    resetDataVisualization();
  }
  double GetIterationGap() {
    return iterationGap_;
  }

  void SetBatchSize(double bs) {
    batchSize_ = bs;
    Modified();
    resetDataVisualization();
  }
  double GetBatchSize() {
    return batchSize_;
  }

  void SetOptimizer(int optimizer) {
    optimizer_ = optimizer;
    Modified();
    resetDataVisualization();
  }
  int GetOptimizer() {
    return optimizer_;
  }

  void SetGradientStepSize(double lr) {
    gradientStepSize_ = lr;
    Modified();
    resetDataVisualization();
  }
  double GetGradientStepSize() {
    return gradientStepSize_;
  }

  void SetBeta1(double beta) {
    beta1_ = beta;
    Modified();
    resetDataVisualization();
  }
  double GetBeta1() {
    return beta1_;
  }

  void SetBeta2(double beta) {
    beta2_ = beta;
    Modified();
    resetDataVisualization();
  }
  double GetBeta2() {
    return beta2_;
  }

  void SetReconstructionLossWeight(double reconstructionLossWeight) {
    reconstructionLossWeight_ = reconstructionLossWeight;
    Modified();
    resetDataVisualization();
  }
  double GetReconstructionLossWeight() {
    return reconstructionLossWeight_;
  }

  void SetTrackingLossWeight(double trackingLossWeight) {
    trackingLossWeight_ = trackingLossWeight;
    Modified();
    resetDataVisualization();
  }
  double GetTrackingLossWeight() {
    return trackingLossWeight_;
  }

  void SetMetricLossWeight(double metricLossWeight) {
    metricLossWeight_ = metricLossWeight;
    Modified();
    resetDataVisualization();
  }
  double GetMetricLossWeight() {
    return metricLossWeight_;
  }

  void SetClusteringLossWeight(double clusteringLossWeight) {
    clusteringLossWeight_ = clusteringLossWeight;
    Modified();
    resetDataVisualization();
  }
  double GetClusteringLossWeight() {
    return clusteringLossWeight_;
  }

  void SetCustomLossSpace(bool customLossSpace) {
    customLossSpace_ = customLossSpace;
    Modified();
    resetDataVisualization();
  }
  bool GetCustomLossSpace() {
    return customLossSpace_;
  }

  void SetCustomLossActivate(bool customLossActivate) {
    customLossActivate_ = customLossActivate;
    Modified();
    resetDataVisualization();
  }
  bool GetCustomLossActivate() {
    return customLossActivate_;
  }

  void SetNormalizeMetricLoss(bool normalizeMetricLoss) {
    normalizeMetricLoss_ = normalizeMetricLoss;
    Modified();
    resetDataVisualization();
  }
  bool GetNormalizeMetricLoss() {
    return normalizeMetricLoss_;
  }

  void SetClusteringLossTemperature(double clusteringLossTemperature) {
    clusteringLossTemp_ = clusteringLossTemperature;
    Modified();
    resetDataVisualization();
  }
  double GetClusteringLossTemperature() {
    return clusteringLossTemp_;
  }

  void SetCustomLossDynamicWeight(bool customLossDynamicWeight) {
    customLossDynamicWeight_ = customLossDynamicWeight;
    Modified();
    resetDataVisualization();
  }
  bool GetCustomLossDynamicWeight() {
    return customLossDynamicWeight_;
  }

  void SetNumberOfInit(unsigned int noInit) {
    noInit_ = noInit;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetNumberOfInit() {
    return noInit_;
  }

  void SetEuclideanVectorsInit(bool euclideanVectorsInit) {
    euclideanVectorsInit_ = euclideanVectorsInit;
    Modified();
    resetDataVisualization();
  }
  bool GetEuclideanVectorsInit() {
    return euclideanVectorsInit_;
  }

  void SetInitOriginPrimeStructByCopy(bool initOriginPrimeStructByCopy) {
    initOriginPrimeStructByCopy_ = initOriginPrimeStructByCopy;
    Modified();
    resetDataVisualization();
  }
  bool GetInitOriginPrimeStructByCopy() {
    return initOriginPrimeStructByCopy_;
  }

  void SetTrackingLossDecoding(double trackingLossDecoding) {
    trackingLossDecoding_ = trackingLossDecoding;
    Modified();
    resetDataVisualization();
  }
  double GetTrackingLossDecoding() {
    return trackingLossDecoding_;
  }

  void SetTrackingLossInitRandomness(double trackingLossInitRandomness) {
    trackingLossInitRandomness_ = trackingLossInitRandomness;
    Modified();
    resetDataVisualization();
  }
  double GetTrackingLossInitRandomness() {
    return trackingLossInitRandomness_;
  }

  void SetDeterministic(bool deterministic) {
    deterministic_ = deterministic;
    Modified();
    resetDataVisualization();
  }
  bool GetDeterministic() {
    return deterministic_;
  }

  void SetActivate(bool activate) {
    activate_ = activate;
    Modified();
    resetDataVisualization();
  }
  bool GetActivate() {
    return activate_;
  }

  void SetActivationFunction(unsigned int activationFunction) {
    activationFunction_ = activationFunction;
    Modified();
    resetDataVisualization();
  }
  unsigned int GetActivationFunction() {
    return activationFunction_;
  }

  void SetFullSymmetricAE(bool fullSymmetricAE) {
    fullSymmetricAE_ = fullSymmetricAE;
    Modified();
    resetDataVisualization();
  }
  bool GetFullSymmetricAE() {
    return fullSymmetricAE_;
  }

  void SetActivateOutputInit(bool activateOutputInit) {
    activateOutputInit_ = activateOutputInit;
    Modified();
    resetDataVisualization();
  }
  bool GetActivateOutputInit() {
    return activateOutputInit_;
  }

  void SetJoinSplitMixtureCoefficient(double joinSplitMixtureCoefficient) {
    mixtureCoefficient_ = joinSplitMixtureCoefficient;
    Modified();
    resetDataVisualization();
  }
  double GetJoinSplitMixtureCoefficient() {
    return mixtureCoefficient_;
  }

  void SetEpsilon1UseFarthestSaddle(bool epsilon1UseFarthestSaddle) {
    epsilon1UseFarthestSaddle_ = epsilon1UseFarthestSaddle;
    Modified();
    resetDataVisualization();
  }
  bool GetEpsilon1UseFarthestSaddle() {
    return epsilon1UseFarthestSaddle_;
  }

  void SetEpsilonTree1(double epsilonTree1) {
    epsilonTree1_ = epsilonTree1;
    oldEpsilonTree1 = epsilonTree1_;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilonTree1() {
    return epsilonTree1_;
  }

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    epsilon2Tree1_ = epsilon2Tree1;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilon2Tree1() {
    return epsilon2Tree1_;
  }

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    epsilon3Tree1_ = epsilon3Tree1;
    Modified();
    resetDataVisualization();
  }
  double GetEpsilon3Tree1() {
    return epsilon3Tree1_;
  }

  void SetPersistenceThreshold(double persistenceThreshold) {
    persistenceThreshold_ = persistenceThreshold;
    Modified();
    resetDataVisualization();
  }
  double GetPersistenceThreshold() {
    return persistenceThreshold_;
  }

  void SetDeleteMultiPersPairs(bool delMultiPersPairs) {
    deleteMultiPersPairs_ = delMultiPersPairs;
    Modified();
    resetDataVisualization();
  }
  bool GetDeleteMultiPersPairs() {
    return deleteMultiPersPairs_;
  }

  // Testing options
  void SetNodePerTask(int nodePerTask) {
    nodePerTask_ = nodePerTask;
    Modified();
    resetDataVisualization();
  }
  int GetNodePerTask() {
    return nodePerTask_;
  }

  // Output options
  void SetCreateOutput(bool createOutput) {
    createOutput_ = createOutput;
    Modified();
    resetDataVisualization();
  }
  bool GetCreateOutput() {
    return createOutput_;
  }

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeAutoencoder *New();
  vtkTypeMacro(ttkMergeTreeAutoencoder, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreeAutoencoder();

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int run(vtkInformationVector *outputVector,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
          std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);

  int runCompute(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);

  int runOutput(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);
};
