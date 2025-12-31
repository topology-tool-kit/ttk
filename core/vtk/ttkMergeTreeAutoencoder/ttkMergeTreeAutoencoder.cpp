#include <MergeTreeAxesAlgorithmUtils.h>
#include <ttkMergeTreeAutoencoder.h>
#include <ttkMergeTreeNeuralNetworkUtils.h>
#include <ttkMergeTreeUtils.h>
#include <ttkMergeTreeVisualization.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnsignedIntArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeAutoencoder);

/**
 * Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMergeTreeAutoencoder::ttkMergeTreeAutoencoder() {
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(4);
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeAutoencoder::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else if(port == 2) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else
    return 0;
  return 1;
}

/**
 * Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMergeTreeAutoencoder::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0 or port == 1 or port == 2 or port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else {
    return 0;
  }
  return 1;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkMergeTreeAutoencoder::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {
#ifndef TTK_ENABLE_TORCH
  TTK_FORCE_USE(inputVector);
  TTK_FORCE_USE(outputVector);
  printErr("This filter requires Torch.");
  return 0;
#else
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);
  auto table = vtkTable::GetData(inputVector[2], 0);

  // ------------------------------------------------------------------------------------
  // --- Load blocks
  // ------------------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTrees, inputTrees2;
  ttk::ftm::loadBlocks(inputTrees, blocks);
  ttk::ftm::loadBlocks(inputTrees2, blocks2);

  // Load table
  clusterAsgn_.clear();
  vtkAbstractArray *clusterAsgn;
  if(table) {
    clusterAsgn = this->GetInputArrayToProcess(0, inputVector);
    if(clusterAsgn) {
      clusterAsgn_.resize(clusterAsgn->GetNumberOfValues());
      for(unsigned int i = 0; i < clusterAsgn_.size(); ++i)
        clusterAsgn_[i] = clusterAsgn->GetVariantValue(i).ToInt();
    }
  }
  if((not table or not clusterAsgn) and clusteringLossWeight_ != 0) {
    printErr(
      "You must provide a table column in info input to use clustering loss");
    return 0;
  }
  if(clusteringLossWeight_ != 0) {
    std::stringstream ss;
    for(auto &e : clusterAsgn_)
      ss << e << " ";
    printMsg(ss.str());
  }

  // ------------------------------------------------------------------------------------
  // If we have already computed once but the input has changed
  if((treesNodes.size() != 0 and inputTrees[0]->GetBlock(0) != treesNodes[0])
     or (treesNodes2.size() != inputTrees2.size()))
    resetDataVisualization();

  // Parameters
  branchDecomposition_ = true;
  if(not normalizedWasserstein_) {
    oldEpsilonTree1 = epsilonTree1_;
    epsilonTree1_ = 100;
  } else
    epsilonTree1_ = oldEpsilonTree1;
  if(normalizedWasserstein_)
    printMsg("Computation with normalized Wasserstein.");
  else
    printMsg("Computation without normalized Wasserstein.");

  return run(outputVector, inputTrees, inputTrees2);
#endif
}

#ifdef TTK_ENABLE_TORCH
int ttkMergeTreeAutoencoder::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  runCompute(outputVector, inputTrees, inputTrees2);
  runOutput(outputVector, inputTrees, inputTrees2);
  return 1;
}

int ttkMergeTreeAutoencoder::runCompute(
  vtkInformationVector *ttkNotUsed(outputVector),
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  std::vector<ttk::ftm::MergeTree<float>> intermediateMTrees,
    intermediateMTrees2;

  bool useSecondPairsType = (mixtureCoefficient_ == 0);
  isPersistenceDiagram_ = ttk::ftm::constructTrees<float>(
    inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation,
    useSecondPairsType, DiagramPairTypes);
  // If merge trees are provided in input and normalization is not asked
  convertToDiagram_
    = (not isPersistenceDiagram_ and not normalizedWasserstein_);
  if(not isPersistenceDiagram_
     or (mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1)) {
    auto &inputTrees2ToUse
      = (not isPersistenceDiagram_ ? inputTrees2 : inputTrees);
    ttk::ftm::constructTrees<float>(inputTrees2ToUse, intermediateMTrees2,
                                    treesNodes2, treesArcs2, treesSegmentation2,
                                    !useSecondPairsType, DiagramPairTypes);
  }
  isPersistenceDiagram_ |= (not normalizedWasserstein_);

  const int numInputs = intermediateMTrees.size();
  const int numInputs2 = intermediateMTrees2.size();
  setDataVisualization(numInputs, numInputs2);

  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  execute(intermediateMTrees, intermediateMTrees2);

  ttk::ftm::mergeTreesTemplateToDouble<float>(
    intermediateMTrees, intermediateDTrees);

  return 1;
}

// TODO manage double input
int ttkMergeTreeAutoencoder::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &ttkNotUsed(inputTrees2)) {
  if(not createOutput_)
    return 1;
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  auto output_origins = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_vectors = vtkMultiBlockDataSet::GetData(outputVector, 1);
  auto output_coef = vtkMultiBlockDataSet::GetData(outputVector, 2);
  auto output_data = vtkMultiBlockDataSet::GetData(outputVector, 3);

  // ------------------------------------------
  // --- Tracking information
  // ------------------------------------------
  std::vector<std::vector<ttk::ftm::idNode>> originsMatchingVectorT,
    invOriginsMatchingVectorT;
  std::vector<std::vector<std::vector<ttk::ftm::idNode>>>
    invDataMatchingVectorT;
  std::vector<std::vector<ttk::ftm::idNode>> invReconstMatchingVectorT;
  ttk::wnn::makeMatchingVectors(
    originsMatchings_, originsCopy_, originsPrimeCopy_, originsMatchingVectorT,
    invOriginsMatchingVectorT, dataMatchings_, recs_, invDataMatchingVectorT,
    reconstMatchings_, invReconstMatchingVectorT);

  std::vector<std::vector<ttk::ftm::idNode>> originsMatchingVector;
  std::vector<std::vector<double>> originsPersPercent, originsPersDiff;
  std::vector<double> originPersPercent, originPersDiff;
  std::vector<int> originPersistenceOrder;
  ttk::wnn::computeTrackingInformation(
    originsCopy_, originsPrimeCopy_, originsMatchingVectorT,
    invOriginsMatchingVectorT, isPersistenceDiagram_, originsMatchingVector,
    originsPersPercent, originsPersDiff, originPersPercent, originPersDiff,
    originPersistenceOrder);

  // ------------------------------------------
  // --- Data
  // ------------------------------------------
  ttk::wnn::makeDataOutput(
    output_data, recs_, 1, treesSegmentation, persCorrelationMatrix_,
    invDataMatchingVectorT, invReconstMatchingVectorT, originsMatchingVectorT,
    originsMatchingVector, originsPersPercent, originsPersDiff,
    originPersistenceOrder, treesNodes, treesNodeCorr_, bestLoss_,
    mixtureCoefficient_, isPersistenceDiagram_, convertToDiagram_,
    this->debugLevel_);

  // ------------------------------------------
  // --- Origins
  // ------------------------------------------
  ttk::wnn::makeOriginsOutput(
    output_origins, originsCopy_, originsPrimeCopy_, originPersPercent,
    originPersDiff, originPersistenceOrder, originsMatchingVector,
    originsPersPercent, originsPersDiff, mixtureCoefficient_,
    isPersistenceDiagram_, convertToDiagram_, this->debugLevel_);

  // ------------------------------------------
  // --- Coefficients
  // ------------------------------------------
  ttk::wnn::makeCoefficientsOutput(output_coef, allAlphas_, allScaledAlphas_,
                                   allActAlphas_, allActScaledAlphas_,
                                   clusterAsgn_, recs_, inputTrees);

  // Field Data Input Parameters
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  for(auto paramName : paramNames) {
    vtkNew<vtkDoubleArray> array{};
    array->SetName(paramName.c_str());
    array->InsertNextTuple1(getParamValueFromName(paramName));
    output_coef->GetFieldData()->AddArray(array);
  }
  vtkNew<vtkDoubleArray> arrayActivate{};
  arrayActivate->SetName("activate");
  arrayActivate->InsertNextTuple1(activate_);
  output_coef->GetFieldData()->AddArray(arrayActivate);
  vtkNew<vtkDoubleArray> arrayActivateFunction{};
  arrayActivateFunction->SetName("activationFunction");
  arrayActivateFunction->InsertNextTuple1(activationFunction_);
  output_coef->GetFieldData()->AddArray(arrayActivateFunction);

  vtkNew<vtkIntArray> diagramPairTypesArray{};
  diagramPairTypesArray->SetName("DiagramPairTypes");
  diagramPairTypesArray->InsertNextTuple1(DiagramPairTypes);
  output_coef->GetFieldData()->AddArray(diagramPairTypesArray);

  // ------------------------------------------
  // --- Axes Vectors
  // ------------------------------------------
  std::vector<std::vector<std::vector<ttk::ftm::idNode>>> dataMatchingVectorT(
    dataMatchings_.size());
  for(unsigned int l = 0; l < dataMatchingVectorT.size(); ++l) {
    dataMatchingVectorT[l].resize(dataMatchings_[l].size());
    for(unsigned int i = 0; i < dataMatchingVectorT[l].size(); ++i) {
      auto &origin = (l == 0 ? originsCopy_[0] : originsPrimeCopy_[l - 1]);
      ttk::axa::getMatchingVector(origin.mTree, recs_[i][l].mTree,
                                  dataMatchings_[l][i],
                                  dataMatchingVectorT[l][i]);
    }
  }
  output_vectors->SetNumberOfBlocks(2);
  vtkSmartPointer<vtkMultiBlockDataSet> vectors
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  vectors->SetNumberOfBlocks(noLayers_);
  vtkSmartPointer<vtkMultiBlockDataSet> vectorsPrime
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  vectorsPrime->SetNumberOfBlocks(noLayers_);
  for(unsigned int l = 0; l < noLayers_; ++l) {
    vtkSmartPointer<vtkTable> vectorsTable = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkTable> vectorsPrimeTable
      = vtkSmartPointer<vtkTable>::New();
    for(unsigned int v = 0; v < layers_[l].getVSTensor().sizes()[1]; ++v) {
      // Vs
      vtkNew<vtkFloatArray> vectorArray{};
      std::string name = ttk::axa::getTableVectorName(
        layers_[l].getVSTensor().sizes()[1], v, 0, 0, false);
      vectorArray->SetName(name.c_str());
      vectorArray->SetNumberOfTuples(layers_[l].getVSTensor().sizes()[0]);
      for(unsigned int i = 0; i < layers_[l].getVSTensor().sizes()[0]; ++i)
        vectorArray->SetTuple1(i, layers_[l].getVSTensor()[i][v].item<float>());
      vectorsTable->AddColumn(vectorArray);
      // Vs Prime
      vtkNew<vtkFloatArray> vectorPrimeArray{};
      std::string name2 = ttk::axa::getTableVectorName(
        layers_[l].getVSTensor().sizes()[1], v, 0, 0, false);
      vectorPrimeArray->SetName(name2.c_str());
      vectorPrimeArray->SetNumberOfTuples(
        layers_[l].getVSPrimeTensor().sizes()[0]);
      for(unsigned int i = 0; i < layers_[l].getVSPrimeTensor().sizes()[0]; ++i)
        vectorPrimeArray->SetTuple1(
          i, layers_[l].getVSPrimeTensor()[i][v].item<float>());
      vectorsPrimeTable->AddColumn(vectorPrimeArray);
    }
    // Rev node corr
    vtkNew<vtkUnsignedIntArray> revNodeCorrArray{};
    revNodeCorrArray->SetName("revNodeCorr");
    revNodeCorrArray->SetNumberOfTuples(layers_[l].getVSTensor().sizes()[0]);
    std::vector<unsigned int> revNodeCorr;
    getReverseTorchNodeCorr(originsCopy_[l], revNodeCorr);
    for(unsigned int i = 0; i < layers_[l].getVSTensor().sizes()[0]; ++i)
      revNodeCorrArray->SetTuple1(i, revNodeCorr[i]);
    vectorsTable->AddColumn(revNodeCorrArray);
    // Rev node corr prime
    vtkNew<vtkUnsignedIntArray> revNodeCorrPrimeArray{};
    revNodeCorrPrimeArray->SetNumberOfTuples(
      layers_[l].getVSPrimeTensor().sizes()[0]);
    revNodeCorrPrimeArray->SetName("revNodeCorr");
    std::vector<unsigned int> revNodeCorrPrime;
    getReverseTorchNodeCorr(originsPrimeCopy_[l], revNodeCorrPrime);
    for(unsigned int i = 0; i < layers_[l].getVSPrimeTensor().sizes()[0]; ++i)
      revNodeCorrPrimeArray->SetTuple1(i, revNodeCorrPrime[i]);
    vectorsPrimeTable->AddColumn(revNodeCorrPrimeArray);
    // Origins Matchings
    auto addOriginMatchingArray
      = [&](vtkSmartPointer<vtkTable> &vectorsTableT,
            std::vector<ttk::ftm::idNode> &originMatchingVector) {
          vtkNew<vtkIntArray> matchingArray{};
          matchingArray->SetNumberOfTuples(originMatchingVector.size());
          matchingArray->SetName("nextOriginMatching");
          for(unsigned int i = 0; i < originMatchingVector.size(); ++i)
            matchingArray->SetTuple1(i, (int)originMatchingVector[i]);
          vectorsTableT->AddColumn(matchingArray);
        };
    if(l == 0)
      addOriginMatchingArray(vectorsTable, originsMatchingVectorT[l]);
    if(l < originsMatchingVectorT.size() - 1)
      addOriginMatchingArray(vectorsPrimeTable, originsMatchingVectorT[l + 1]);
    // Data Matchings
    auto addDataMatchingArray
      = [&](vtkSmartPointer<vtkTable> &vectorsTableT,
            std::vector<std::vector<ttk::ftm::idNode>> &dataMatchingVector) {
          for(unsigned int i = 0; i < dataMatchingVector.size(); ++i) {
            vtkNew<vtkIntArray> matchingArray{};
            matchingArray->SetNumberOfTuples(dataMatchingVector[i].size());
            std::stringstream ss;
            ss << "matching"
               << ttk::axa::getTableTreeName(dataMatchingVector.size(), i);
            matchingArray->SetName(ss.str().c_str());
            for(unsigned int j = 0; j < dataMatchingVector[i].size(); ++j)
              matchingArray->SetTuple1(j, (int)dataMatchingVector[i][j]);
            vectorsTableT->AddColumn(matchingArray);
          }
        };
    if(l == 0)
      addDataMatchingArray(vectorsTable, dataMatchingVectorT[l]);
    if(l < dataMatchingVectorT.size() - 1)
      addDataMatchingArray(vectorsPrimeTable, dataMatchingVectorT[l + 1]);
    // Reconst Matchings
    if(l == noLayers_ - 1) {
      for(unsigned int i = 0; i < invReconstMatchingVectorT.size(); ++i) {
        vtkNew<vtkIntArray> matchingArray{};
        matchingArray->SetNumberOfTuples(invReconstMatchingVectorT[i].size());
        std::stringstream ss;
        ss << "reconstMatching"
           << ttk::axa::getTableTreeName(invReconstMatchingVectorT.size(), i);
        matchingArray->SetName(ss.str().c_str());
        for(unsigned int j = 0; j < invReconstMatchingVectorT[i].size(); ++j)
          matchingArray->SetTuple1(j, (int)invReconstMatchingVectorT[i][j]);
        vectorsPrimeTable->AddColumn(matchingArray);
      }
    }

    // Add new block
    vectors->SetBlock(l, vectorsTable);
    std::stringstream ss;
    ss << "Vectors" << l;
    vectors->GetMetaData(l)->Set(vtkCompositeDataSet::NAME(), ss.str());
    vectorsPrime->SetBlock(l, vectorsPrimeTable);
    ss.str("");
    ss << "VectorsPrime" << l;
    vectorsPrime->GetMetaData(l)->Set(vtkCompositeDataSet::NAME(), ss.str());
  }
  output_vectors->SetBlock(0, vectors);
  output_vectors->SetBlock(1, vectorsPrime);
  unsigned int num = 0;
  output_vectors->GetMetaData(num)->Set(vtkCompositeDataSet::NAME(), "Vectors");
  num = 1;
  output_vectors->GetMetaData(num)->Set(
    vtkCompositeDataSet::NAME(), "VectorsPrime");

  // return success
  return 1;
}
#endif
