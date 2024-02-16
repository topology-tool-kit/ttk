#include <MergeTreeAxesAlgorithmUtils.h>
#include <MergeTreeTorchUtils.h>
#include <ttkMergeTreeAutoencoderDecoding.h>
#include <ttkMergeTreeAutoencoderUtils.h>
#include <ttkMergeTreeUtils.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeAutoencoderDecoding);

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
ttkMergeTreeAutoencoderDecoding::ttkMergeTreeAutoencoderDecoding() {
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(1);
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeAutoencoderDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 or port == 1 or port == 2) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
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
int ttkMergeTreeAutoencoderDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
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
int ttkMergeTreeAutoencoderDecoding::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
#ifndef TTK_ENABLE_TORCH
  TTK_FORCE_USE(inputVector);
  TTK_FORCE_USE(outputVector);
  printErr("This filter requires Torch.");
  return 0;
#else
  // --------------------------------------------------------------------------
  // --- Read Input
  // --------------------------------------------------------------------------
  auto blockBary = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto vectors = vtkMultiBlockDataSet::GetData(inputVector[1], 0);
  auto coefficients = vtkMultiBlockDataSet::GetData(inputVector[2], 0);

  // Parameters
  printMsg("Load parameters from field data.");
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  paramNames.emplace_back("activate");
  paramNames.emplace_back("activationFunction");
  auto fd = coefficients->GetFieldData();
  if(not fd->GetArray("activate"))
    fd = coefficients->GetBlock(0)->GetFieldData();
  for(auto paramName : paramNames) {
    auto array = fd->GetArray(paramName.c_str());
    if(array) {
      double const value = array->GetTuple1(0);
      if(paramName == "activate")
        activate_ = value;
      else if(paramName == "activationFunction")
        activationFunction_ = value;
      else
        setParamValueFromName(paramName, value);
    } else
      printMsg(" - " + paramName + " was not found in the field data.");
    auto stringValue = std::to_string(
      (paramName == "activate" ? activate_
                               : (paramName == "activationFunction"
                                    ? activationFunction_
                                    : getParamValueFromName(paramName))));
    printMsg(" - " + paramName + " = " + stringValue);
  }
  if(normalizedWasserstein_)
    printMsg("Computation with normalized Wasserstein.");
  else
    printMsg("Computation without normalized Wasserstein.");

  // -----------------
  // Origins
  // -----------------
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> origins, originsPrime;
  ttk::ftm::loadBlocks(
    origins, vtkMultiBlockDataSet::SafeDownCast((blockBary->GetBlock(0))));
  ttk::ftm::loadBlocks(
    originsPrime, vtkMultiBlockDataSet::SafeDownCast(blockBary->GetBlock(1)));

  std::vector<ttk::ftm::MergeTree<float>> originsTrees, originsPrimeTrees;
  std::vector<vtkUnstructuredGrid *> originsTreeNodes, originsTreeArcs,
    originsPrimeTreeNodes, originsPrimeTreeArcs;
  std::vector<vtkDataSet *> originsTreeSegmentations,
    originsPrimeTreeSegmentations;

  bool useSadMaxPairs = (mixtureCoefficient_ == 0);
  isPersistenceDiagram_ = ttk::ftm::constructTrees<float>(
    origins, originsTrees, originsTreeNodes, originsTreeArcs,
    originsTreeSegmentations, useSadMaxPairs);
  ttk::ftm::constructTrees<float>(originsPrime, originsPrimeTrees,
                                  originsTreeNodes, originsTreeArcs,
                                  originsTreeSegmentations, useSadMaxPairs);
  // If merge trees are provided in input and normalization is not asked
  convertToDiagram_
    = (not isPersistenceDiagram_ and not normalizedWasserstein_);

  // -----------------
  // Number of axes per layer
  // -----------------
  auto vS = vtkMultiBlockDataSet::SafeDownCast(vectors->GetBlock(0));
  noLayers_ = vS->GetNumberOfBlocks();
  std::vector<unsigned int> allNoAxes(noLayers_);
  for(unsigned int l = 0; l < noLayers_; ++l) {
    auto layerVectorsTable = vtkTable::SafeDownCast(vS->GetBlock(l));
    unsigned int maxBoundNoAxes = 1;
    while(not layerVectorsTable->GetColumnByName(
      ttk::axa::getTableVectorName(maxBoundNoAxes, 0, 0, 0, false).c_str()))
      maxBoundNoAxes *= 10;
    unsigned int noAxes = 1;
    while(layerVectorsTable->GetColumnByName(
      ttk::axa::getTableVectorName(maxBoundNoAxes, noAxes, 0, 0, false)
        .c_str()))
      noAxes += 1;
    allNoAxes[l] = noAxes;
  }

  // -----------------
  // Coefficients
  // -----------------
  auto numberOfInputs
    = vtkTable::SafeDownCast(coefficients->GetBlock(0))->GetNumberOfRows();
  auto noCoefs = coefficients->GetNumberOfBlocks();
  bool customRec = (noCoefs != noLayers_);
  allAlphas_.resize(numberOfInputs, std::vector<torch::Tensor>(noCoefs));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(unsigned int l = 0; l < noCoefs; ++l) {
    auto layerCoefficientsTable
      = vtkTable::SafeDownCast(coefficients->GetBlock(l));
    auto noAxes = (customRec ? allNoAxes[getLatentLayerIndex()] : allNoAxes[l]);
    std::vector<std::vector<float>> alphas(
      numberOfInputs, std::vector<float>(noAxes));
    for(unsigned int g = 0; g < noAxes; ++g) {
      auto array = layerCoefficientsTable->GetColumnByName(
        ttk::axa::getTableCoefficientName(noAxes, g).c_str());
      for(unsigned int i = 0; i < numberOfInputs; ++i)
        alphas[i][g] = array->GetVariantValue(i).ToFloat();
    }
    for(unsigned int i = 0; i < numberOfInputs; ++i)
      allAlphas_[i][l] = torch::tensor(alphas[i]).reshape({-1, 1});
  }

  // -----------------
  // Vectors
  // -----------------
  vSTensor_.resize(noLayers_);
  vSPrimeTensor_.resize(noLayers_);
  auto vSPrime = vtkMultiBlockDataSet::SafeDownCast(vectors->GetBlock(1));
  std::vector<unsigned int *> allRevNodeCorr(noLayers_),
    allRevNodeCorrPrime(noLayers_);
  std::vector<unsigned int> allRevNodeCorrSize(noLayers_),
    allRevNodeCorrPrimeSize(noLayers_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(unsigned int l = 0; l < vSTensor_.size(); ++l) {
    auto layerVectorsTable = vtkTable::SafeDownCast(vS->GetBlock(l));
    auto layerVectorsPrimeTable = vtkTable::SafeDownCast(vSPrime->GetBlock(l));
    auto noRows = layerVectorsTable->GetNumberOfRows();
    auto noRows2 = layerVectorsPrimeTable->GetNumberOfRows();
    std::vector<float> vSTensor(noRows * allNoAxes[l]),
      vSPrimeTensor(noRows2 * allNoAxes[l]);
    for(unsigned int v = 0; v < allNoAxes[l]; ++v) {
      std::string name
        = ttk::axa::getTableVectorName(allNoAxes[l], v, 0, 0, false);
      for(unsigned int i = 0; i < noRows; ++i)
        vSTensor[i * allNoAxes[l] + v]
          = layerVectorsTable->GetColumnByName(name.c_str())
              ->GetVariantValue(i)
              .ToFloat();
      for(unsigned int i = 0; i < noRows2; ++i)
        vSPrimeTensor[i * allNoAxes[l] + v]
          = layerVectorsPrimeTable->GetColumnByName(name.c_str())
              ->GetVariantValue(i)
              .ToFloat();
    }
    vSTensor_[l] = torch::tensor(vSTensor).reshape({noRows, allNoAxes[l]});
    vSPrimeTensor_[l]
      = torch::tensor(vSPrimeTensor).reshape({noRows2, allNoAxes[l]});
    allRevNodeCorr[l]
      = ttkUtils::GetPointer<unsigned int>(vtkDataArray::SafeDownCast(
        layerVectorsTable->GetColumnByName("revNodeCorr")));
    allRevNodeCorrSize[l] = noRows;
    allRevNodeCorrPrime[l]
      = ttkUtils::GetPointer<unsigned int>(vtkDataArray::SafeDownCast(
        layerVectorsPrimeTable->GetColumnByName("revNodeCorr")));
    allRevNodeCorrPrimeSize[l] = noRows2;
  }

  // -----------------
  // Call base
  // -----------------
  execute(originsTrees, originsPrimeTrees, allRevNodeCorr, allRevNodeCorrPrime,
          allRevNodeCorrSize, allRevNodeCorrPrimeSize);

  // --------------------------------------------------------------------------
  // --- Create Output
  // --------------------------------------------------------------------------
  auto output_data = vtkMultiBlockDataSet::GetData(outputVector, 0);

  // ------------------------------------------
  // --- Read Matchings
  // ------------------------------------------
  auto originsMatchingSize = getLatentLayerIndex() + 1;
  std::vector<std::vector<ttk::ftm::idNode>> originsMatchingVectorT(
    originsMatchingSize),
    invOriginsMatchingVectorT = originsMatchingVectorT;
  for(unsigned int l = 0; l < originsMatchingVectorT.size(); ++l) {
    auto array = vtkTable::SafeDownCast(
                   (l == 0 ? vS : vSPrime)->GetBlock((l == 0 ? l : l - 1)))
                   ->GetColumnByName("nextOriginMatching");
    originsMatchingVectorT[l].clear();
    originsMatchingVectorT[l].resize(array->GetNumberOfTuples());
    for(unsigned int i = 0; i < originsMatchingVectorT[l].size(); ++i)
      originsMatchingVectorT[l][i] = array->GetVariantValue(i).ToUnsignedInt();
    reverseMatchingVector<float>(originsPrime_[l].mTree,
                                 originsMatchingVectorT[l],
                                 invOriginsMatchingVectorT[l]);
  }
  auto dataMatchingSize = getLatentLayerIndex() + 2;
  std::vector<std::vector<std::vector<ttk::ftm::idNode>>> dataMatchingVectorT(
    dataMatchingSize),
    invDataMatchingVectorT = dataMatchingVectorT;
  std::vector<std::vector<ttk::ftm::idNode>> invReconstMatchingVectorT(
    numberOfInputs);
  if(not customRec) {
    for(unsigned int l = 0; l < dataMatchingVectorT.size(); ++l) {
      dataMatchingVectorT[l].resize(numberOfInputs);
      invDataMatchingVectorT[l].resize(numberOfInputs);
      for(unsigned int i = 0; i < dataMatchingVectorT[l].size(); ++i) {
        std::stringstream ss;
        ss << "matching"
           << ttk::axa::getTableTreeName(dataMatchingVectorT[l].size(), i);
        auto array = vtkTable::SafeDownCast(
                       (l == 0 ? vS : vSPrime)->GetBlock((l == 0 ? l : l - 1)))
                       ->GetColumnByName(ss.str().c_str());
        dataMatchingVectorT[l][i].clear();
        dataMatchingVectorT[l][i].resize(array->GetNumberOfTuples());
        for(unsigned int j = 0; j < dataMatchingVectorT[l][i].size(); ++j)
          dataMatchingVectorT[l][i][j]
            = array->GetVariantValue(j).ToUnsignedInt();
        auto noNodes
          = (l == 0 ? vtkTable::SafeDownCast(coefficients->GetBlock(0))
                        ->GetColumnByName("treeNoNodes")
                        ->GetVariantValue(i)
                        .ToUnsignedInt()
                    : recs_[i][l - 1].mTree.tree.getNumberOfNodes());
        reverseMatchingVector(
          noNodes, dataMatchingVectorT[l][i], invDataMatchingVectorT[l][i]);
      }
    }
    for(unsigned int i = 0; i < invReconstMatchingVectorT.size(); ++i) {
      std::stringstream ss;
      ss << "reconstMatching"
         << ttk::axa::getTableTreeName(invReconstMatchingVectorT.size(), i);
      auto array = vtkTable::SafeDownCast(vSPrime->GetBlock(noLayers_ - 1))
                     ->GetColumnByName(ss.str().c_str());
      invReconstMatchingVectorT[i].resize(array->GetNumberOfTuples());
      for(unsigned int j = 0; j < invReconstMatchingVectorT[i].size(); ++j)
        invReconstMatchingVectorT[i][j]
          = array->GetVariantValue(j).ToUnsignedInt();
    }
  }

  // ------------------------------------------
  // --- Tracking information
  // ------------------------------------------
  std::vector<std::vector<ttk::ftm::idNode>> originsMatchingVector;
  std::vector<std::vector<double>> originsPersPercent, originsPersDiff;
  std::vector<double> originPersPercent, originPersDiff;
  std::vector<int> originPersistenceOrder;
  ttk::wae::computeTrackingInformation(
    origins_, originsPrime_, originsMatchingVectorT, invOriginsMatchingVectorT,
    isPersistenceDiagram_, originsMatchingVector, originsPersPercent,
    originsPersDiff, originPersPercent, originPersDiff, originPersistenceOrder);

  // ------------------------------------------
  // --- Data
  // ------------------------------------------
  if(!recs_.empty()) {
    output_data->SetNumberOfBlocks(1);
    vtkSmartPointer<vtkMultiBlockDataSet> data
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    data->SetNumberOfBlocks(recs_[0].size());
    vtkSmartPointer<vtkMultiBlockDataSet> dataSeg
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    dataSeg->SetNumberOfBlocks(recs_.size());
    for(unsigned int l = 0; l < recs_[0].size(); ++l) {
      vtkSmartPointer<vtkMultiBlockDataSet> out_layer_i
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      out_layer_i->SetNumberOfBlocks(recs_.size());
      std::vector<ttk::ftm::MergeTree<float> *> trees(recs_.size());
      for(unsigned int i = 0; i < recs_.size(); ++i)
        trees[i] = &(recs_[i][l].mTree);

      // Custom arrays
      std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
        customIntArrays(recs_.size());
      std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
        customDoubleArrays(recs_.size());
      unsigned int lShift = 1;
      ttk::wae::computeCustomArrays(
        recs_, persCorrelationMatrix_, invDataMatchingVectorT,
        invReconstMatchingVectorT, originsMatchingVector,
        originsMatchingVectorT, originsPersPercent, originsPersDiff,
        originPersistenceOrder, l, lShift, customIntArrays, customDoubleArrays);

      // Create output
      ttk::wae::makeManyOutput(trees, out_layer_i, customIntArrays,
                               customDoubleArrays, mixtureCoefficient_,
                               isPersistenceDiagram_, convertToDiagram_,
                               this->debugLevel_);
      data->SetBlock(l, out_layer_i);
      std::stringstream ss;
      ss << "Layer" << l;
      data->GetMetaData(l)->Set(vtkCompositeDataSet::NAME(), ss.str());
    }
    output_data->SetBlock(0, data);
    unsigned int num = 0;
    std::stringstream ss;
    ss << "layers" << (isPersistenceDiagram_ ? "Diagrams" : "Trees");
    output_data->GetMetaData(num)->Set(vtkCompositeDataSet::NAME(), ss.str());
  }

  if(!customRecs_.empty()) {
    std::vector<std::vector<std::tuple<std::string, std::vector<int>>>>
      customRecsIntArrays(customRecs_.size());
    std::vector<std::vector<std::tuple<std::string, std::vector<double>>>>
      customRecsDoubleArrays(customRecs_.size());
    std::vector<ttk::ftm::MergeTree<float> *> trees(customRecs_.size());
    for(unsigned int i = 0; i < customRecs_.size(); ++i)
      trees[i] = &(customRecs_[i].mTree);
    std::vector<std::vector<int>> customOriginPersOrder(customRecs_.size());
    for(unsigned int i = 0; i < customRecs_.size(); ++i) {
      trees[i] = &(customRecs_[i].mTree);
      std::vector<ttk::ftm::idNode> matchingVector;
      getInverseMatchingVector(origins_[0].mTree, customRecs_[i].mTree,
                               customMatchings_[i], matchingVector);
      customOriginPersOrder[i].resize(
        customRecs_[i].mTree.tree.getNumberOfNodes());
      for(unsigned int j = 0; j < matchingVector.size(); ++j) {
        if(matchingVector[j] < originPersistenceOrder.size())
          customOriginPersOrder[i][j]
            = originPersistenceOrder[matchingVector[j]];
        else
          customOriginPersOrder[i][j] = -1;
      }
      std::string name4{"OriginPersOrder"};
      customRecsIntArrays[i].emplace_back(
        std::make_tuple(name4, customOriginPersOrder[i]));
    }
    vtkSmartPointer<vtkMultiBlockDataSet> dataCustom
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    ttk::wae::makeManyOutput(trees, dataCustom, customRecsIntArrays,
                             customRecsDoubleArrays, mixtureCoefficient_,
                             isPersistenceDiagram_, convertToDiagram_,
                             this->debugLevel_);
    output_data->SetBlock(0, dataCustom);
  }

  // return success
  return 1;
#endif
}
