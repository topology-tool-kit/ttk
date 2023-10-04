#include <ttkMergeAndContourTreeUtils.h>
#include <ttkMergeTreePrincipalGeodesics.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreePrincipalGeodesics);

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
ttkMergeTreePrincipalGeodesics::ttkMergeTreePrincipalGeodesics() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(4);
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreePrincipalGeodesics::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
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
int ttkMergeTreePrincipalGeodesics::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else if(port == 1 or port == 2 or port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
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
int ttkMergeTreePrincipalGeodesics::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // ------------------------------------------------------------------------------------
  // --- Load blocks
  // ------------------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTrees, inputTrees2;
  ttk::ftm::loadBlocks(inputTrees, blocks);
  ttk::ftm::loadBlocks(inputTrees2, blocks2);

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

  return run<float>(outputVector, inputTrees, inputTrees2);
}

template <class dataType>
int ttkMergeTreePrincipalGeodesics::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(outputVector, inputTrees, inputTrees2);
  runOutput<dataType>(outputVector, inputTrees, inputTrees2);
  return 1;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesics::runCompute(
  vtkInformationVector *ttkNotUsed(outputVector),
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  std::vector<ttk::ftm::MergeTree<dataType>> intermediateMTrees,
    intermediateMTrees2;

  bool const useSadMaxPairs = (mixtureCoefficient_ == 0);
  isPersistenceDiagram_ = ttk::ftm::constructTrees<dataType>(
    inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation,
    useSadMaxPairs);
  // If merge trees are provided in input and normalization is not asked
  convertToDiagram_
    = (not isPersistenceDiagram_ and not normalizedWasserstein_);
  if(convertToDiagram_) {
    mixtureCoefficient_
      = (intermediateMTrees[0].tree.template isJoinTree<dataType>() ? 1 : 0);
  }
  if(not isPersistenceDiagram_
     or (mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1)) {
    auto &inputTrees2ToUse
      = (not isPersistenceDiagram_ ? inputTrees2 : inputTrees);
    ttk::ftm::constructTrees<dataType>(inputTrees2ToUse, intermediateMTrees2,
                                       treesNodes2, treesArcs2,
                                       treesSegmentation2, !useSadMaxPairs);
  }
  isPersistenceDiagram_ |= (not normalizedWasserstein_);

  const int numInputs = intermediateMTrees.size();
  const int numInputs2 = intermediateMTrees2.size();
  setDataVisualization(numInputs, numInputs2);

  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  execute<dataType>(intermediateMTrees, intermediateMTrees2);

  ttk::ftm::mergeTreesTemplateToDouble<dataType>(
    intermediateMTrees, intermediateDTrees);

  return 1;
}

template <class dataType>
void ttkMergeTreePrincipalGeodesics::makeBarycenterOutput(
  ttk::ftm::MergeTree<dataType> &barycenter,
  int blockId,
  vtkMultiBlockDataSet *output_barycenter) {
  vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNodeBary
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputArcBary
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  ttkMergeTreeVisualization visuMakerBary;
  visuMakerBary.setShiftMode(2); // Line
  visuMakerBary.setVtkOutputNode(vtkOutputNodeBary);
  if(not isPersistenceDiagram_)
    visuMakerBary.setVtkOutputArc(vtkOutputArcBary);
  else {
    visuMakerBary.setVtkOutputArc(vtkOutputNodeBary);
    if(mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1)
      visuMakerBary.setIsPDSadMax(blockId);
    else
      visuMakerBary.setIsPDSadMax(mixtureCoefficient_ == 0);
    visuMakerBary.setPlanarLayout(true);
  }
  visuMakerBary.setDebugLevel(this->debugLevel_);
  visuMakerBary.setOutputTreeNodeId(true);
  visuMakerBary.setIsPersistenceDiagram(isPersistenceDiagram_);
  visuMakerBary.setConvertedToDiagram(convertToDiagram_);
  visuMakerBary.makeTreesOutput<dataType>(&(barycenter.tree));

  // Construct multiblock
  if(not isPersistenceDiagram_) {
    vtkMultiBlockDataSet::SafeDownCast(output_barycenter->GetBlock(0))
      ->SetBlock(blockId, vtkOutputNodeBary);
    vtkMultiBlockDataSet::SafeDownCast(output_barycenter->GetBlock(1))
      ->SetBlock(blockId, vtkOutputArcBary);
  } else
    output_barycenter->SetBlock(blockId, vtkOutputNodeBary);
}

template <class dataType>
int ttkMergeTreePrincipalGeodesics::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  auto output_barycenter = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_coef = vtkTable::GetData(outputVector, 1);
  auto output_geodVector = vtkTable::GetData(outputVector, 2);
  auto output_corrMatrix = vtkTable::GetData(outputVector, 3);

  // ------------------------------------------
  // --- Barycenter
  // ------------------------------------------
  auto noBlockBary
    = 1
      + (not isPersistenceDiagram_
         or (mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1));
  output_barycenter->SetNumberOfBlocks(noBlockBary);
  if(not isPersistenceDiagram_) {
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockNodes
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockArcs
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    int const noBlock = (inputTrees2.size() == 0 ? 1 : 2);
    vtkBlockNodes->SetNumberOfBlocks(noBlock);
    vtkBlockArcs->SetNumberOfBlocks(noBlock);
    output_barycenter->SetBlock(0, vtkBlockNodes);
    output_barycenter->SetBlock(1, vtkBlockArcs);
  }

  ttk::ftm::MergeTree<dataType> barycenterT;
  mergeTreeDoubleToTemplate(barycenter_, barycenterT);
  makeBarycenterOutput<dataType>(barycenterT, 0, output_barycenter);

  if(inputTrees2.size() != 0
     or (isPersistenceDiagram_ and mixtureCoefficient_ != 0
         and mixtureCoefficient_ != 1)) {
    ttk::ftm::MergeTree<dataType> barycenter2T;
    mergeTreeDoubleToTemplate(barycenterInput2_, barycenter2T);
    makeBarycenterOutput<dataType>(barycenter2T, 1, output_barycenter);
  }

  // ------------------------------------------
  // --- Coefficients
  // ------------------------------------------
  std::vector<std::vector<std::vector<double>> *> allTs{&allTs_, &allScaledTs_};

  vtkNew<vtkIntArray> treeIDArray{};
  treeIDArray->SetName("TreeID");
  treeIDArray->SetNumberOfTuples(inputTrees.size());
  for(unsigned int i = 0; i < inputTrees.size(); ++i) {
    treeIDArray->SetTuple1(i, i);
  }
  output_coef->AddColumn(treeIDArray);

  for(unsigned int t = 0; t < 2; ++t) {
    for(unsigned int i = 0; i < (*allTs[t]).size(); ++i) {
      vtkNew<vtkDoubleArray> tArray{};
      auto size = (*allTs[t]).size();
      std::string const name = (t == 0 ? getTableCoefficientName(size, i)
                                       : getTableCoefficientNormName(size, i));
      tArray->SetName(name.c_str());
      tArray->SetNumberOfTuples(inputTrees.size());
      for(unsigned int j = 0; j < (*allTs[t])[i].size(); ++j)
        tArray->SetTuple1(j, (*allTs[t])[i][j]);
      output_coef->AddColumn(tArray);
    }
  }

  // Copy Field Data
  // - aggregate input field data
  for(unsigned int b = 0; b < inputTrees[0]->GetNumberOfBlocks(); ++b) {
    vtkNew<vtkFieldData> fd{};
    fd->CopyStructure(inputTrees[0]->GetBlock(b)->GetFieldData());
    fd->SetNumberOfTuples(inputTrees.size());
    for(size_t i = 0; i < inputTrees.size(); ++i) {
      fd->SetTuple(i, 0, inputTrees[i]->GetBlock(b)->GetFieldData());
    }

    // - copy input field data to output row data
    for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
      auto array = fd->GetAbstractArray(i);
      array->SetName(array->GetName());
      output_coef->AddColumn(array);
    }
  }

  // Field Data Input Parameters
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  for(auto paramName : paramNames) {
    vtkNew<vtkDoubleArray> array{};
    array->SetName(paramName.c_str());
    array->InsertNextTuple1(getParamValueFromName(paramName));
    output_coef->GetFieldData()->AddArray(array);
  }

  // ------------------------------------------
  // --- Geodesics Vectors
  // ------------------------------------------
  std::vector<std::vector<std::vector<std::vector<double>>> *> vectors{
    &vS_, &v2s_};
  unsigned int maxSize = vS_[0].size();
  if(inputTrees2.size() != 0
     or (isPersistenceDiagram_ and mixtureCoefficient_ != 0
         and mixtureCoefficient_ != 1)) {
    vectors = {&vS_, &v2s_, &trees2Vs_, &trees2V2s_};
    maxSize = std::max(vS_[0].size(), trees2Vs_[0].size());
  }
  for(unsigned int v = 0; v < vectors.size(); ++v) {
    auto vector = vectors[v];
    for(unsigned int i = 0; i < (*vector).size(); ++i)
      for(unsigned int k = 0; k < 2; ++k) {
        vtkNew<vtkDoubleArray> vectorArray{};
        bool const isSecondInput = (v >= 2);
        std::string const name
          = getTableVectorName((*vector).size(), i, v % 2, k, isSecondInput);
        vectorArray->SetName(name.c_str());
        vectorArray->SetNumberOfTuples(maxSize);
        for(unsigned int j = 0; j < (*vector)[i].size(); ++j)
          vectorArray->SetTuple1(j, (*vector)[i][j][k]);
        for(unsigned int j = (*vector)[i].size(); j < maxSize; ++j)
          vectorArray->SetTuple1(j, std::nan(""));
        output_geodVector->AddColumn(vectorArray);
      }
  }

  // ------------------------------------------
  // --- Correlation Matrix
  // ------------------------------------------
  // TODO manage second input
  // Correlation coefficients
  unsigned int const noCols = branchesCorrelationMatrix_[0].size();
  unsigned int const noRows = branchesCorrelationMatrix_.size();
  for(unsigned int j = 0; j < noCols; ++j) {
    vtkNew<vtkDoubleArray> corrArray{}, persArray{};
    std::string const name = getTableCorrelationName(noCols, j);
    corrArray->SetName(name.c_str());
    corrArray->SetNumberOfTuples(noRows);
    std::string const name2 = getTableCorrelationPersName(noCols, j);
    persArray->SetName(name2.c_str());
    persArray->SetNumberOfTuples(noRows);
    for(unsigned int i = 0; i < noRows; ++i) {
      corrArray->SetTuple1(i, branchesCorrelationMatrix_[i][j]);
      persArray->SetTuple1(i, persCorrelationMatrix_[i][j]);
    }
    output_corrMatrix->AddColumn(corrArray);
    output_corrMatrix->AddColumn(persArray);
  }

  // Tree matching
  std::vector<std::vector<ttk::ftm::idNode>> matchingMatrix;
  getMatchingMatrix(
    barycenter_, intermediateDTrees, baryMatchings_, matchingMatrix);
  if(not normalizedWasserstein_)
    for(unsigned int j = 0; j < inputTrees.size(); ++j)
      matchingMatrix[barycenter_.tree.getNode(barycenter_.tree.getRoot())
                       ->getOrigin()][j]
        = intermediateDTrees[j]
            .tree.getNode(intermediateDTrees[j].tree.getRoot())
            ->getOrigin();
  for(unsigned int j = 0; j < inputTrees.size(); ++j) {
    vtkNew<vtkIntArray> corrArray{};
    std::string const name = getTableCorrelationTreeName(inputTrees.size(), j);
    corrArray->SetName(name.c_str());
    corrArray->SetNumberOfTuples(noRows);
    for(unsigned int i = 0; i < noRows; ++i)
      corrArray->SetTuple1(i, (int)matchingMatrix[i][j]);
    output_corrMatrix->AddColumn(corrArray);
  }

  // Barycenter node id
  vtkNew<vtkIntArray> baryNodeId{};
  baryNodeId->SetName("BaryNodeId");
  baryNodeId->SetNumberOfTuples(noRows);
  for(unsigned int i = 0; i < noRows; ++i)
    baryNodeId->SetTuple1(i, i);
  output_corrMatrix->AddColumn(baryNodeId);

  // Barycenter node isLeaf
  vtkNew<vtkIntArray> baryNodeIsLeaf{};
  baryNodeIsLeaf->SetName("BaryNodeIsLeaf");
  baryNodeIsLeaf->SetNumberOfTuples(noRows);
  for(unsigned int i = 0; i < noRows; ++i)
    baryNodeIsLeaf->SetTuple1(i, barycenter_.tree.isLeaf(i));
  output_corrMatrix->AddColumn(baryNodeIsLeaf);

  // Barycenter node persistence
  vtkNew<vtkDoubleArray> baryNodePers{};
  baryNodePers->SetName("BaryNodePers");
  baryNodePers->SetNumberOfTuples(noRows);
  for(unsigned int i = 0; i < noRows; ++i) {
    auto nodePers = barycenter_.tree.template getNodePersistence<double>(i);
    baryNodePers->SetTuple1(i, nodePers);
  }
  output_corrMatrix->AddColumn(baryNodePers);

  // Trees persistence order
  std::vector<std::vector<int>> treesOrder(intermediateDTrees.size());
  for(unsigned int i = 0; i < intermediateDTrees.size(); ++i) {
    treesOrder[i].resize(intermediateDTrees[i].tree.getNumberOfNodes());
    std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>> pairs;
    intermediateDTrees[i].tree.template getPersistencePairsFromTree<double>(
      pairs, false);
    for(unsigned int j = 0; j < pairs.size(); ++j) {
      int const index = pairs.size() - 1 - j;
      treesOrder[i][std::get<0>(pairs[j])] = index;
      treesOrder[i][std::get<1>(pairs[j])] = index;
    }
  }

  vtkNew<vtkDoubleArray> treesPersOrder{};
  treesPersOrder->SetName("TreesPersOrder");
  treesPersOrder->SetNumberOfTuples(noRows);
  for(unsigned int i = 0; i < noRows; ++i) {
    int minIndex = std::numeric_limits<int>::max();
    for(unsigned int j = 0; j < intermediateDTrees.size(); ++j) {
      if(int(matchingMatrix[i][j]) != -1)
        minIndex = std::min(minIndex, treesOrder[j][matchingMatrix[i][j]]);
    }
    treesPersOrder->SetTuple1(i, minIndex);
  }
  output_corrMatrix->AddColumn(treesPersOrder);

  return 1;
}
