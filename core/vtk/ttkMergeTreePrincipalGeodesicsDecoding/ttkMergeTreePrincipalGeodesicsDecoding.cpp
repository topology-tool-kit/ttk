#include <ttkMergeAndContourTreeUtils.h>
#include <ttkMergeTreePrincipalGeodesicsDecoding.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkVariantArray.h>

#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreePrincipalGeodesicsDecoding);

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
ttkMergeTreePrincipalGeodesicsDecoding::
  ttkMergeTreePrincipalGeodesicsDecoding() {
  this->SetNumberOfInputPorts(5);
  this->SetNumberOfOutputPorts(1);
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreePrincipalGeodesicsDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else if(port == 1 or port == 2 or port == 3) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    if(port == 3) // correlation table
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else if(port == 4) {
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
int ttkMergeTreePrincipalGeodesicsDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else
    return 0;
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
int ttkMergeTreePrincipalGeodesicsDecoding::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  auto blockBary = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto tableCoefficients = vtkTable::GetData(inputVector[1], 0);
  auto tableVectors = vtkTable::GetData(inputVector[2], 0);
  auto tableCorrelation = vtkTable::GetData(inputVector[3], 0);
  auto blockInputTrees = vtkMultiBlockDataSet::GetData(inputVector[4], 0);

  // ------------------------------------------------------------------------------------
  // --- Load blocks
  // ------------------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputBary, inputTrees;
  ttk::ftm::loadBlocks(inputBary, blockBary);
  ttk::ftm::loadBlocks(inputTrees, blockInputTrees);

  // If we have already computed once but the input has changed
  if((baryTreeNodes.size() != 0
      and inputBary[0]->GetBlock(0) != baryTreeNodes[0])
     or (allTs_.size() != 0
         and int(allTs_[0].size()) != tableCoefficients->GetNumberOfRows())
     or tableCoefficientsMTime != tableCoefficients->GetMTime()
     or tableVectorsMTime != tableVectors->GetMTime()
     or (tableCorrelation
         and tableCorrelationMTime != tableCorrelation->GetMTime())
     or (blockInputTrees
         and blockInputTreesMTime != blockInputTrees->GetMTime())) {
    resetDataVisualization();
  }
  tableCoefficientsMTime = tableCoefficients->GetMTime();
  tableVectorsMTime = tableVectors->GetMTime();
  if(tableCorrelation)
    tableCorrelationMTime = tableCorrelation->GetMTime();
  if(blockInputTrees)
    blockInputTreesMTime = blockInputTrees->GetMTime();

  // Parameters
  printMsg("Load parameters from field data.");
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  for(auto paramName : paramNames) {
    auto array = tableCoefficients->GetFieldData()->GetArray(paramName.c_str());
    if(array) {
      double const value = array->GetTuple1(0);
      setParamValueFromName(paramName, value);
    } else
      printMsg(" - " + paramName + " was not found in the field data.");
    printMsg(" - " + paramName + " = "
             + std::to_string(getParamValueFromName(paramName)));
  }
  if(normalizedWasserstein_)
    printMsg("Computation with normalized Wasserstein.");
  else
    printMsg("Computation without normalized Wasserstein.");

  // ------------------------------------------------------------------------------------
  // --- Load tables
  // ------------------------------------------------------------------------------------
  auto tableVNoCols = tableVectors->GetNumberOfColumns() / inputBary.size();
  unsigned int const numberOfGeodesics = tableVNoCols / 4;
  auto numberOfInputs = tableCoefficients->GetNumberOfRows();

  auto baryNoNodes = tableVectors->GetNumberOfRows();
  auto baryNoNodes2 = baryNoNodes;
  if(inputBary.size() > 1) {
    for(unsigned int i = 0; i < 2; ++i) {
      auto &baryNoNodesT = (i == 0 ? baryNoNodes : baryNoNodes2);
      while(std::isnan(
        tableVectors
          ->GetColumnByName(
            getTableVectorName(numberOfGeodesics, 0, 0, 0, (i == 1)).c_str())
          ->GetVariantValue(baryNoNodesT - 1)
          .ToDouble()))
        --baryNoNodesT;
    }
  }

  // -----------------
  // Coefficients
  // -----------------
  // Pointers can't be used because we need to access the rows (VariantArray)
  std::vector<std::string> nonFieldDataNames{"TreeID"};
  allTs_.resize(numberOfGeodesics, std::vector<double>(numberOfInputs, 0.0));
  for(unsigned int i = 0; i < numberOfGeodesics; ++i)
    for(unsigned int j = 0; j < numberOfInputs; ++j) {
      std::string const name = getTableCoefficientName(numberOfGeodesics, i);
      allTs_[i][j] = tableCoefficients->GetColumnByName(name.c_str())
                       ->GetVariantValue(j)
                       .ToDouble();
      nonFieldDataNames.push_back(name);
      std::string const normName
        = getTableCoefficientNormName(numberOfGeodesics, i);
      nonFieldDataNames.push_back(normName);
    }
  ttk::Geometry::transposeMatrix(allTs_, allTreesTs_);

  // - aggregate input field data
  vtkNew<vtkFieldData> fd{};
  for(int i = 0; i < tableCoefficients->GetNumberOfColumns(); ++i) {
    std::string const name{tableCoefficients->GetColumnName(i)};
    if(std::count(nonFieldDataNames.begin(), nonFieldDataNames.end(), name)
       == 0) {
      auto array = tableCoefficients->GetColumn(i);
      array->SetName(name.c_str());
      fd->AddArray(array);
    }
  }
  inputFieldData = vtkFieldData::New();
  inputFieldData->ShallowCopy(fd);

  // -----------------
  // Vectors
  // -----------------
  pVS_.resize(numberOfGeodesics, std::vector<double *>(2));
  pV2s_ = pVS_;
  if(inputBary.size() > 1) {
    pTrees2Vs_.resize(numberOfGeodesics, std::vector<double *>(2));
    pTrees2V2s_ = pTrees2Vs_;
  }
  for(unsigned int h = 0; h < inputBary.size(); ++h) {
    auto &pVS = (h == 0 ? pVS_ : pTrees2Vs_);
    auto &pV2s = (h == 0 ? pV2s_ : pTrees2V2s_);
    bool const secondInput = (h == 1);
    for(unsigned int i = 0; i < numberOfGeodesics; ++i) {
      for(unsigned int k = 0; k < 2; ++k) {
        std::string const name1
          = getTableVectorName(numberOfGeodesics, i, 0, k, secondInput);
        std::string const name2
          = getTableVectorName(numberOfGeodesics, i, 1, k, secondInput);
        pVS[i][k] = ttkUtils::GetPointer<double>(vtkDataArray::SafeDownCast(
          tableVectors->GetColumnByName(name1.c_str())));
        pV2s[i][k] = ttkUtils::GetPointer<double>(vtkDataArray::SafeDownCast(
          tableVectors->GetColumnByName(name2.c_str())));
      }
    }
  }
  vSize_ = baryNoNodes;
  vSize2_ = baryNoNodes2;

  // -----------------
  // Correlation
  // -----------------
  if(tableCorrelation) {
    pBranchesCorrelationMatrix_.resize(numberOfGeodesics);
    pPersCorrelationMatrix_ = pBranchesCorrelationMatrix_;

    auto tableCorrNoRows = tableCorrelation->GetNumberOfRows();
    auto baryNodeIdArray = tableCorrelation->GetColumnByName("BaryNodeId");

    for(unsigned int j = 0; j < numberOfGeodesics; ++j) {
      std::string const name = getTableCorrelationName(numberOfGeodesics, j);
      std::string const name2
        = getTableCorrelationPersName(numberOfGeodesics, j);
      pBranchesCorrelationMatrix_[j]
        = ttkUtils::GetPointer<double>(vtkDataArray::SafeDownCast(
          tableCorrelation->GetColumnByName(name.c_str())));
      pPersCorrelationMatrix_[j]
        = ttkUtils::GetPointer<double>(vtkDataArray::SafeDownCast(
          tableCorrelation->GetColumnByName(name2.c_str())));
    }

    // Pointers can't be used because processing on these data will be done
    // later
    baryMatchings_.resize(numberOfInputs);
    for(unsigned int i = 0; i < numberOfInputs; ++i) {
      std::string const name = getTableCorrelationTreeName(numberOfInputs, i);
      auto array = tableCorrelation->GetColumnByName(name.c_str());
      if(array) {
        baryMatchings_[i].resize(tableCorrNoRows);
        for(unsigned int j = 0; j < tableCorrNoRows; ++j) {
          auto baryNodeId = baryNodeIdArray->GetVariantValue(j).ToInt();
          baryMatchings_[i][j] = std::make_tuple(
            baryNodeId, array->GetVariantValue(j).ToDouble(), 0.0);
        }
      }
    }
  } else {
    pBranchesCorrelationMatrix_.clear();
    pPersCorrelationMatrix_.clear();
    baryMatchings_.clear();
  }

  // ------------------------------------------------------------------------------------

  return run<float>(outputVector, inputBary, inputTrees);
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsDecoding::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputBary,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(outputVector, inputBary, inputTrees);
  runOutput<dataType>(outputVector, inputBary, inputTrees);
  return 1;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsDecoding::runCompute(
  vtkInformationVector *ttkNotUsed(outputVector),
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputBary,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees) {
  useDoubleInput_ = (inputBary.size() > 1);
  processFirstInput = (not useDoubleInput_ or not ProcessSecondInput);
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  const int numTrees = inputTrees.size();

  setDataVisualization(numTrees);

  std::vector<ttk::ftm::MergeTree<dataType>> baryDTree, inputDTrees;

  std::vector<bool> useSadMaxPairsVec{false, true};
  if(not useDoubleInput_ and mixtureCoefficient_ == 0)
    useSadMaxPairsVec.erase(useSadMaxPairsVec.begin()); // {true}
  ttk::ftm::constructTrees<dataType>(inputBary, baryDTree, baryTreeNodes,
                                     baryTreeArcs, baryTreeSegmentation,
                                     useSadMaxPairsVec);

  if(OutputInputTrees
     or (ReconstructInputTrees
         and (computeReconstructionError_ or transferInputTreesInformation_))) {
    bool const useSadMaxPairs
      = (useDoubleInput_ and not processFirstInput) or mixtureCoefficient_ == 0;
    bool const isInputPD = ttk::ftm::constructTrees<dataType>(
      inputTrees, inputDTrees, inputTreesNodes, inputTreesArcs,
      inputTreesSegmentation, useSadMaxPairs);
    if(not isInputPD and isPersistenceDiagram_)
      mtsFlattening(inputDTrees);
  }

  //------------------------------------------------------------------------------------
  // --- Call base
  //------------------------------------------------------------------------------------
  std::string const preFix{(processFirstInput ? "" : "Second Input ")};
  auto &baryToUse = (processFirstInput ? baryDTree[0] : baryDTree[1]);

  // Geodesics distances
  if(geodesicsDistances_.size() == 0) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Compute Geodesics Distance...");
    computeGeodesicsDistance<dataType>(baryDTree);
  }

  // Input Trees
  if(OutputInputTrees) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Process Input Trees...");
    processInputTrees<dataType>(inputDTrees);
  }

  // Reconstructed Trees
  if(ReconstructInputTrees) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Reconstruct " + preFix + "Input Trees...");
    std::vector<ttk::ftm::MergeTree<dataType>> reconstructedTreesT;
    reconstruction<dataType>(baryToUse, inputDTrees, reconstructedTreesT,
                             reconstructionErrors, recInputMatchings,
                             recBaryMatchings, !processFirstInput);
    ttk::ftm::mergeTreesTemplateToDouble<dataType>(
      reconstructedTreesT, reconstructedTrees);
  } else
    reconstructedTrees.clear();

  // Geodesics Trees
  if(ConstructGeodesicsTrees) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Construct Geodesics Trees...");
    std::vector<std::vector<ttk::ftm::MergeTree<dataType>>> geodesicsTreesT;
    constructGeodesicsTrees<dataType>(
      baryToUse, geodesicsTreesT, !processFirstInput);
    ttk::ftm::mergeTreesTemplateToDouble<dataType>(
      geodesicsTreesT, geodesicsTrees);
  } else
    geodesicsTrees.clear();

  // Ellipses Trees
  if(ConstructEllipses) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Construct Ellipses Trees...");
    std::vector<ttk::ftm::MergeTree<dataType>> geodesicsEllipsesT;
    constructGeodesicsEllipses<dataType>(
      baryToUse, geodesicsEllipsesT, !processFirstInput);
    ttk::ftm::mergeTreesTemplateToDouble<dataType>(
      geodesicsEllipsesT, geodesicsEllipses);
  } else
    geodesicsEllipses.clear();

  // Rectangle Trees
  if(ConstructRectangle) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Construct Rectangle Trees...");
    std::vector<ttk::ftm::MergeTree<dataType>> geodesicsRectangleT;
    constructGeodesicsRectangle<dataType>(
      baryToUse, geodesicsRectangleT, RectangleMultiplier, !processFirstInput);
    ttk::ftm::mergeTreesTemplateToDouble<dataType>(
      geodesicsRectangleT, geodesicsRectangle);
  } else
    geodesicsRectangle.clear();

  // Surface Trees
  if(ConstructSurface) {
    printMsg(ttk::debug::Separator::L2);
    printMsg("Construct " + preFix + "Surface Trees...");
    std::vector<ttk::ftm::MergeTree<dataType>> geodesicsSurfaceT;
    constructGeodesicsSurface<dataType>(
      baryToUse, geodesicsSurfaceT, !processFirstInput);
    ttk::ftm::mergeTreesTemplateToDouble<dataType>(
      geodesicsSurfaceT, geodesicsSurface);
  } else
    geodesicsSurface.clear();

  ttk::ftm::mergeTreesTemplateToDouble<dataType>(baryDTree, baryMTree);
  ttk::ftm::mergeTreesTemplateToDouble<dataType>(inputDTrees, inputMTrees);

  return 1;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsDecoding::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &ttkNotUsed(inputBary),
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees) {
  printMsg(ttk::debug::Separator::L2);
  printMsg("Make output...");
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  auto output = vtkMultiBlockDataSet::GetData(outputVector, 0);

  unsigned int const noInput = (OutputInputTrees ? inputTrees.size() : 0);
  unsigned int const noReconst = reconstructedTrees.size();
  unsigned int const noBary = (OutputBarycenter ? 1 : 0);
  unsigned int const noGeod
    = geodesicsTrees.size()
      * (geodesicsTrees.size() == 0 ? 0 : geodesicsTrees[0].size());
  unsigned int const noEllipses = geodesicsEllipses.size();
  unsigned int const noRectangle = geodesicsRectangle.size();
  unsigned int const noSurface = geodesicsSurface.size();

  printMsg("noInput     = " + std::to_string(noInput));
  printMsg("noReconst   = " + std::to_string(noReconst));
  printMsg("noBary      = " + std::to_string(noBary));
  printMsg("noGeod      = " + std::to_string(noGeod));
  printMsg("noEllipses  = " + std::to_string(noEllipses));
  printMsg("noRectangle = " + std::to_string(noRectangle));
  printMsg("noSurface   = " + std::to_string(noSurface));

  std::vector<unsigned int> allNo{
    noInput, noReconst, noBary, noGeod, noEllipses, noRectangle, noSurface};
  std::vector<unsigned int> steps(allNo.size());
  for(unsigned int i = 0; i < allNo.size(); ++i)
    steps[i] = (i == 0 ? 0 : steps[i - 1]) + allNo[i];

  unsigned int const noTreesTotal = steps[steps.size() - 1];

  unsigned int const noBlocks
    = 2 - isPersistenceDiagram_ + OutputInputTreesSegmentation;
  output->SetNumberOfBlocks(noBlocks);
  vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockNodes
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  vtkBlockNodes->SetNumberOfBlocks(noTreesTotal);
  output->SetBlock(0, vtkBlockNodes);
  if(not isPersistenceDiagram_) {
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockArcs
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkBlockArcs->SetNumberOfBlocks(noTreesTotal);
    output->SetBlock(1, vtkBlockArcs);
  }
  if(OutputInputTreesSegmentation) {
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockSegs
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkBlockSegs->SetNumberOfBlocks(noTreesTotal);
    output->SetBlock(2 - isPersistenceDiagram_, vtkBlockSegs);
  }

  int const geodesicsTreesOffset = std::max((int)geodesicsTrees.size() - 1, 0);

  // ------------------------------------------
  // --- Trees
  // ------------------------------------------
  std::vector<std::vector<ttk::ftm::idNode>> matchingMatrix;
  if(!baryMatchings_.empty())
    getMatchingMatrix<double>(
      baryMTree[0], inputMTrees, baryMatchings_, matchingMatrix);
  // TODO compute matching to barycenter if correlation matrix is not provided
  if(transferInputTreesInformation_
     and (inputMTrees.empty() or baryMatchings_.empty()))
    printWrn("Please provide input trees and correlation matrix to transfer "
             "input trees information.");
    // TODO fix if an interpolation is empty
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(unsigned int i = 0; i < noTreesTotal; ++i) {
    int index = 0;
    int treeType = 0;
    ttk::ftm::MergeTree<dataType> mt;
    std::vector<double> ts;
    if(i < steps[0]) {
      // Input Trees
      index = i;
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(inputMTrees[index], mt);
      treeType = -1;
      ts = allTreesTs_[index];
    } else if(i < steps[1]) {
      // Reconstructed trees
      index = i - steps[0];
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
        reconstructedTrees[index], mt);
      treeType = 0;
      ts = allTreesTs_[index];
    } else if(i < steps[2]) {
      // Barycenter
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(baryMTree[0], mt);
      treeType = 1;
      std::array<double, 2> middle;
      getGeodesicsMiddle<dataType>(mt, pVS_, pV2s_, vSize_, middle);
      ts.resize(2);
      ts[0] = middle[0];
      ts[1] = middle[1];
    } else if(i < steps[3]) {
      // Geodesics Trees
      index = i - steps[2];
      int const i0 = index / geodesicsTrees[0].size();
      int const i1 = index % geodesicsTrees[0].size();
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(geodesicsTrees[i0][i1], mt);
      treeType = 2 + i0;
      ts = tGeodesics_[i0][i1];
    } else if(i < steps[4]) {
      // Ellipses Trees
      index = i - steps[3];
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
        geodesicsEllipses[index], mt);
      treeType = 3 + geodesicsTreesOffset;
      ts = tEllipses_[index];
    } else if(i < steps[5]) {
      // Rectangle Trees
      index = i - steps[4];
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
        geodesicsRectangle[index], mt);
      treeType = 4 + geodesicsTreesOffset;
      ts = tRectangle_[index];
    } else if(i < steps[6]) {
      // Surface Trees
      index = i - steps[5];
      ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
        geodesicsSurface[index], mt);
      treeType = 5 + geodesicsTreesOffset;
      ts = tSurface_[index];
    }

    bool const isInputTree = (i < steps[0]);
    bool const isInputTreeAndGotMesh
      = (isInputTree and inputTrees[i]->GetNumberOfBlocks() >= 3);
    bool const isReconstructedTree = (i >= steps[0] and i < steps[1]);
    bool const isInputOrReconstructedTree
      = (isInputTree or isReconstructedTree);

    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputArc
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputSegmentation
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

    ttkMergeTreeVisualization visuMaker;

    if(isInputTree) {
      for(unsigned int g = 0; g < pBranchesCorrelationMatrix_.size(); ++g) {
        std::vector<double> corr(mt.tree.getNumberOfNodes(), 0.0),
          corrPers = corr;
        for(unsigned int n = 0; n < vSize_; ++n) {
          if((int)matchingMatrix[n][i] != -1) {
            corr[matchingMatrix[n][i]] = pBranchesCorrelationMatrix_[g][n];
            corrPers[matchingMatrix[n][i]] = pPersCorrelationMatrix_[g][n];
          }
        }
        std::string name
          = getTableCorrelationName(pBranchesCorrelationMatrix_.size(), g);
        visuMaker.addCustomArray(name, corr);
        std::string name2
          = getTableCorrelationPersName(pBranchesCorrelationMatrix_.size(), g);
        visuMaker.addCustomArray(name2, corrPers);
      }
      if(!baryMatchings_.empty()) {
        std::vector<double> baryNodePers(mt.tree.getNumberOfNodes(), -1);
        std::vector<int> baryNodeID(mt.tree.getNumberOfNodes(), -1);
        for(unsigned int n = 0; n < vSize_; ++n) {
          if((int)matchingMatrix[n][i] != -1) {
            baryNodeID[matchingMatrix[n][i]] = n;
            baryNodePers[matchingMatrix[n][i]]
              = baryMTree[0].tree.template getNodePersistence<double>(n);
          }
        }
        std::string nameBaryNodePers{"BaryNodePers"};
        visuMaker.addCustomArray(nameBaryNodePers, baryNodePers);
        std::string nameBaryNodeId{"BaryNodeId"};
        visuMaker.addCustomIntArray(nameBaryNodeId, baryNodeID);
      }
      // visuMaker.setShiftMode(0); // Star
      visuMaker.setShiftMode(-1); // None
      visuMaker.setISampleOffset(i);
      visuMaker.setNoSampleOffset(noInput - 1);
      visuMaker.setTreesNodes(inputTreesNodes[index]);
      visuMaker.setTreesNodeCorrMesh(treesNodeCorr_[index]);
      visuMaker.setDimensionSpacing(2);
      visuMaker.setIsPersistenceDiagram(isPersistenceDiagram_);
      if(isInputTreeAndGotMesh and OutputInputTreesSegmentation) {
        visuMaker.setOutputSegmentation(true);
        visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation);
        visuMaker.setTreesSegmentation(inputTreesSegmentation[index]);
      }
    } else
      visuMaker.setShiftMode(2); // Line
    if(isPersistenceDiagram_)
      visuMaker.setPlanarLayout(true);

    // TODO transfer other information?
    if(isReconstructedTree) {
      if(transferBarycenterInformation_) {
        ttk::ftm::MergeTree<dataType> baryMT;
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(baryMTree[0], baryMT);
        std::vector<ttk::ftm::idNode> matchingVector;
        getInverseMatchingVector(
          mt, baryMT, recBaryMatchings[index], matchingVector);
        std::vector<int> baryNodeID(mt.tree.getNumberOfNodes(), -1);
        for(unsigned int n = 0; n < vSize_; ++n) {
          auto match = matchingVector[n];
          if((int)match != -1)
            baryNodeID[match] = n;
        }
        std::string nameBaryNodeId{"BaryNodeId"};
        visuMaker.addCustomIntArray(nameBaryNodeId, baryNodeID);
      }
      if(transferInputTreesInformation_ and !inputMTrees.empty()
         and !baryMatchings_.empty()) {
        ttk::ftm::MergeTree<dataType> inputMT;
        ttk::ftm::mergeTreeDoubleToTemplate<dataType>(
          inputMTrees[index], inputMT);
        std::vector<ttk::ftm::idNode> matchingVector;
        getInverseMatchingVector(
          mt, inputMT, recInputMatchings[index], matchingVector);
        std::vector<int> baryNodeID(mt.tree.getNumberOfNodes(), -1);
        for(unsigned int n = 0; n < vSize_; ++n) {
          auto match = matchingMatrix[n][index];
          if((int)match != -1) {
            match = matchingVector[match];
            if((int)match != -1)
              baryNodeID[match] = n;
          }
        }
        std::string nameBaryNodeId{"BaryNodeId"};
        visuMaker.addCustomIntArray(nameBaryNodeId, baryNodeID);
      }
    }

    visuMaker.setVtkOutputNode(vtkOutputNode);
    if(not isPersistenceDiagram_)
      visuMaker.setVtkOutputArc(vtkOutputArc);
    else {
      visuMaker.setVtkOutputArc(vtkOutputNode);
      if(mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1)
        visuMaker.setIsPDSadMax(!processFirstInput);
      else
        visuMaker.setIsPDSadMax(mixtureCoefficient_ == 0);
    }
    visuMaker.setDebugLevel(this->debugLevel_);
    visuMaker.setIsPersistenceDiagram(isPersistenceDiagram_);
    visuMaker.makeTreesOutput<dataType>(&(mt.tree));

    // Copy Input Field Data
    if(noReconst != 0 or noInput != 0) {
      vtkNew<vtkFieldData> fd{};
      if(isInputOrReconstructedTree)
        fd->DeepCopy(inputFieldData);
      else
        fd->CopyStructure(inputFieldData);
      for(int a = 0; a < inputFieldData->GetNumberOfArrays(); ++a) {
        auto array = fd->GetAbstractArray(a);
        array->SetNumberOfTuples(1);
        if(isInputOrReconstructedTree) {
          vtkNew<vtkIdList> id;
          id->InsertNextId(i % (noInput == 0 ? noReconst : noInput));
          inputFieldData->GetAbstractArray(a)->GetTuples(id, array);
        } else {
          auto downArray = vtkDataArray::SafeDownCast(array);
          if(downArray)
            downArray->SetTuple1(0, std::nan(""));
        }
      }
      vtkOutputNode->GetFieldData()->ShallowCopy(fd);
    }

    // Field Data
    vtkNew<vtkIntArray> treeTypeArray{};
    treeTypeArray->SetName("TreeType");
    treeTypeArray->InsertNextTuple1(treeType);
    vtkOutputNode->GetFieldData()->AddArray(treeTypeArray);

    vtkNew<vtkIntArray> treeIDArray{};
    treeIDArray->SetName("TreeID");
    treeIDArray->InsertNextTuple1(i);
    vtkOutputNode->GetFieldData()->AddArray(treeIDArray);

    for(unsigned int j = 0; j < ts.size(); ++j) {
      vtkNew<vtkDoubleArray> tArray{};
      std::string const name = getTableCoefficientName(ts.size(), j);
      tArray->SetName(name.c_str());
      tArray->InsertNextTuple1(ts[j]);
      vtkOutputNode->GetFieldData()->AddArray(tArray);

      vtkNew<vtkDoubleArray> scaledTArray{};
      std::string const name2 = getTableCoefficientNormName(ts.size(), j);
      scaledTArray->SetName(name2.c_str());
      scaledTArray->InsertNextTuple1(ts[j] * geodesicsDistances_[j]);
      vtkOutputNode->GetFieldData()->AddArray(scaledTArray);
    }

    bool const isSurfaceTree = (treeType == 5 + geodesicsTreesOffset);
    vtkNew<vtkIntArray> isSurfaceArray{};
    isSurfaceArray->SetName("isSurface");
    isSurfaceArray->InsertNextTuple1(isSurfaceTree);
    vtkOutputNode->GetFieldData()->AddArray(isSurfaceArray);
    bool const isBoundary = (isSurfaceTree ? surfaceIsBoundary_[index] : false);
    vtkNew<vtkIntArray> isBoundaryArray{};
    isBoundaryArray->SetName("isBoundary");
    isBoundaryArray->InsertNextTuple1(isBoundary);
    vtkOutputNode->GetFieldData()->AddArray(isBoundaryArray);
    int const boundaryID = (isSurfaceTree ? surfaceBoundaryID_[index] : -1);
    vtkNew<vtkIntArray> BoundaryIDArray{};
    BoundaryIDArray->SetName("BoundaryID");
    BoundaryIDArray->InsertNextTuple1(boundaryID);
    vtkOutputNode->GetFieldData()->AddArray(BoundaryIDArray);

    vtkNew<vtkDoubleArray> reconstructionErrorArray{};
    reconstructionErrorArray->SetName("ReconstructionError");
    auto err = (isReconstructedTree and reconstructionErrors.size() != 0
                  ? reconstructionErrors[index]
                  : std::nan(""));
    reconstructionErrorArray->InsertNextTuple1(err);
    vtkOutputNode->GetFieldData()->AddArray(reconstructionErrorArray);

    // Add new block
    vtkMultiBlockDataSet::SafeDownCast(output->GetBlock(0))
      ->SetBlock(i, vtkOutputNode);
    if(not isPersistenceDiagram_)
      vtkMultiBlockDataSet::SafeDownCast(output->GetBlock(1))
        ->SetBlock(i, vtkOutputArc);
    if(OutputInputTreesSegmentation)
      vtkMultiBlockDataSet::SafeDownCast(
        output->GetBlock(2 - isPersistenceDiagram_))
        ->SetBlock(i, vtkOutputSegmentation);
  }

  // Field Data Input Parameters
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  for(auto paramName : paramNames) {
    vtkNew<vtkDoubleArray> array{};
    array->SetName(paramName.c_str());
    array->InsertNextTuple1(getParamValueFromName(paramName));
    output->GetFieldData()->AddArray(array);
  }

  return 1;
}
