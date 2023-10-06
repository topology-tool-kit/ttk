#include <ttkMergeTreeTemporalReductionDecoding.h>

#include <ttkMergeAndContourTreeUtils.h>
#include <ttkMergeTreeVisualization.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <ttkUtils.h>

using namespace ttk;
using namespace ttk::ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeTemporalReductionDecoding);

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
ttkMergeTreeTemporalReductionDecoding::ttkMergeTreeTemporalReductionDecoding() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

ttkMergeTreeTemporalReductionDecoding::~ttkMergeTreeTemporalReductionDecoding()
  = default;

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeTemporalReductionDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
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
int ttkMergeTreeTemporalReductionDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 or port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
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
int ttkMergeTreeTemporalReductionDecoding::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  printMsg("Get input object from input vector", debug::Priority::VERBOSE);
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto table = vtkTable::GetData(inputVector[1], 0);

  // ------------------------------------------------------------------------------------
  // --- Load blocks
  // ------------------------------------------------------------------------------------
  printMsg("Load Blocks", debug::Priority::VERBOSE);
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTrees;
  loadBlocks(inputTrees, blocks);
  printMsg("Load Blocks done.", debug::Priority::VERBOSE);
  auto assignmentSolverArray
    = blocks->GetFieldData()->GetArray("AssignmentSolver");
  if(assignmentSolverArray)
    assignmentSolverID_ = assignmentSolverArray->GetTuple1(0);

  // If we have already computed once but the input has changed
  if(treesNodes.size() != 0 and inputTrees[0]->GetBlock(0) != treesNodes[0])
    resetDataVisualization();

  std::vector<std::tuple<double, int, int, int, int>> coefs;
  int const totalInputs = inputTrees.size() + table->GetNumberOfRows();
  std::vector<bool> interpolatedTrees(totalInputs, false);
  for(int i = 0; i < table->GetNumberOfRows(); ++i) {
    double const alpha
      = vtkDataArray::SafeDownCast(table->GetColumnByName("Alpha"))
          ->GetTuple1(i);
    int const index1R
      = vtkDataArray::SafeDownCast(table->GetColumnByName("Index1_R"))
          ->GetTuple1(i);
    int const index2R
      = vtkDataArray::SafeDownCast(table->GetColumnByName("Index2_R"))
          ->GetTuple1(i);
    int const index1
      = vtkDataArray::SafeDownCast(table->GetColumnByName("Index1"))
          ->GetTuple1(i);
    int const index2
      = vtkDataArray::SafeDownCast(table->GetColumnByName("Index2"))
          ->GetTuple1(i);

    coefs.emplace_back(alpha, index1R, index2R, index1, index2);
    for(int j = index1 + 1; j < index2; ++j)
      interpolatedTrees[j] = true;
  }

  return run<float>(outputVector, inputTrees, coefs, interpolatedTrees);
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<std::tuple<double, int, int, int, int>> &coefs,
  std::vector<bool> &interpolatedTrees) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(inputTrees, coefs);
  runOutput<dataType>(outputVector, inputTrees, coefs, interpolatedTrees);
  return 1;
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::runCompute(
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<std::tuple<double, int, int, int, int>> &coefs) {
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  printMsg("Construct trees", debug::Priority::VERBOSE);

  const int numInputs = inputTrees.size();

  setDataVisualization(numInputs);

  std::vector<MergeTree<dataType>> intermediateMTrees(numInputs);
  constructTrees<dataType>(
    inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation);

  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  printMsg("Call base", debug::Priority::VERBOSE);

  std::vector<MergeTree<dataType>> allMT_T;

  execute<dataType>(intermediateMTrees, coefs, allMT_T, allMatching);
  treesNodeCorrMesh = getTreesNodeCorr();

  for(size_t i = 0; i < allMT_T.size(); ++i) {
    MergeTree<double> tree;
    mergeTreeTemplateToDouble<dataType>(allMT_T[i], tree);
    intermediateSTrees.emplace_back(tree);
  }

  return 1;
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<std::tuple<double, int, int, int, int>> &coefs,
  std::vector<bool> &interpolatedTrees) {
  bool const OutputSegmentation = false; // TODO
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  auto output_sequence = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_matchings = vtkMultiBlockDataSet::GetData(outputVector, 1);

  // ------------------------------------------
  // --- Trees
  // -----------------------------------------
  std::vector<vtkUnstructuredGrid *> treesNodesT;
  std::vector<vtkUnstructuredGrid *> treesArcsT;
  std::vector<vtkDataSet *> treesSegmentationT;
  std::vector<std::vector<int>> treesNodeCorrMeshT;
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTreesT;
  int index = 0;
  size_t cpt = 0;
  while(cpt < coefs.size()) {
    while(cpt < coefs.size() and std::get<2>(coefs[cpt]) <= index) {
      treesNodesT.emplace_back(nullptr);
      treesArcsT.emplace_back(nullptr);
      treesSegmentationT.emplace_back(nullptr);
      treesNodeCorrMeshT.emplace_back();
      inputTreesT.emplace_back(nullptr);
      ++cpt;
    }
    treesNodesT.emplace_back(treesNodes[index]);
    treesArcsT.emplace_back(treesArcs[index]);
    treesSegmentationT.emplace_back(treesSegmentation[index]);
    treesNodeCorrMeshT.emplace_back(treesNodeCorrMesh[index]);
    inputTreesT.emplace_back(inputTrees[index]);
    ++index;
  }
  treesNodes.swap(treesNodesT);
  treesArcs.swap(treesArcsT);
  treesSegmentation.swap(treesSegmentationT);
  treesNodeCorrMesh.swap(treesNodeCorrMeshT);
  inputTrees.swap(inputTreesT);

  std::vector<MergeTree<dataType>> intermediateMTrees;
  mergeTreesDoubleToTemplate<dataType>(intermediateSTrees, intermediateMTrees);
  std::vector<FTMTree_MT *> trees;
  mergeTreeToFTMTree<dataType>(intermediateMTrees, trees);

  output_sequence->SetNumberOfBlocks((OutputSegmentation ? 3 : 2));
  vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockNodes
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  vtkBlockNodes->SetNumberOfBlocks(intermediateMTrees.size());
  output_sequence->SetBlock(0, vtkBlockNodes);
  vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockArcs
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  vtkBlockArcs->SetNumberOfBlocks(intermediateMTrees.size());
  output_sequence->SetBlock(1, vtkBlockArcs);
  if(OutputSegmentation) {
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockSegs
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkBlockSegs->SetNumberOfBlocks(intermediateMTrees.size());
    output_sequence->SetBlock(2, vtkBlockSegs);
  }

  double prevXMax = 0;
  std::vector<std::vector<SimplexId>> nodeCorr(intermediateMTrees.size());
  int cptInterpolatedTree = 0;
  for(size_t i = 0; i < intermediateMTrees.size(); ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputArc1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputSegmentation1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

    ttkMergeTreeVisualization visuMaker;
    visuMaker.setPlanarLayout(PlanarLayout);
    visuMaker.setBranchDecompositionPlanarLayout(
      BranchDecompositionPlanarLayout);
    visuMaker.setRescaleTreesIndividually(RescaleTreesIndividually);
    visuMaker.setOutputSegmentation(OutputSegmentation);
    visuMaker.setDimensionSpacing(DimensionSpacing);
    visuMaker.setDimensionToShift(DimensionToShift);
    visuMaker.setImportantPairs(ImportantPairs);
    visuMaker.setMaximumImportantPairs(MaximumImportantPairs);
    visuMaker.setMinimumImportantPairs(MinimumImportantPairs);
    visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
    visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
    visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
    visuMaker.setExcludeImportantPairsHigher(ExcludeImportantPairsHigher);
    visuMaker.setExcludeImportantPairsLower(ExcludeImportantPairsLower);
    // visuMaker.setShiftMode(3); // Double Line
    visuMaker.setShiftMode(2); // Line
    visuMaker.setVtkOutputNode(vtkOutputNode1);
    visuMaker.setVtkOutputArc(vtkOutputArc1);
    visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation1);
    visuMaker.setTreesNodes(treesNodes);
    visuMaker.copyPointData(treesNodes[i], treesNodeCorrMesh[i]);
    visuMaker.setTreesNodeCorrMesh(treesNodeCorrMesh);
    visuMaker.setTreesSegmentation(treesSegmentation);
    visuMaker.setInterpolatedTrees(interpolatedTrees);
    visuMaker.setPrintTreeId(i);
    visuMaker.setPrintClusterId(0);
    visuMaker.setPrevXMaxOffset(prevXMax);
    visuMaker.setDebugLevel(this->debugLevel_);
    visuMaker.makeTreesOutput<dataType>(trees);
    prevXMax = visuMaker.getPrevXMax();
    nodeCorr[i] = visuMaker.getNodeCorr()[i];

    // Field data
    if(treesNodes[i])
      vtkOutputNode1->GetFieldData()->ShallowCopy(
        treesNodes[i]->GetFieldData());
    if(treesArcs[i])
      vtkOutputArc1->GetFieldData()->ShallowCopy(treesArcs[i]->GetFieldData());
    if(treesSegmentation[i] and OutputSegmentation)
      vtkOutputSegmentation1->GetFieldData()->ShallowCopy(
        treesSegmentation[i]->GetFieldData());

    // Field data
    if(interpolatedTrees[i]) {
      // Distance previous key frame
      vtkNew<vtkDoubleArray> vtkDistancePreviousKey{};
      vtkDistancePreviousKey->SetName("DistancePreviousKeyFrame");
      vtkDistancePreviousKey->SetNumberOfTuples(1);
      vtkDistancePreviousKey->SetTuple1(
        0, distancesToKeyFrames_[cptInterpolatedTree * 2]);
      vtkOutputNode1->GetFieldData()->AddArray(vtkDistancePreviousKey);
      // Distance next key frame
      vtkNew<vtkDoubleArray> vtkDistanceNextKey{};
      vtkDistanceNextKey->SetName("DistanceNextKeyFrame");
      vtkDistanceNextKey->SetNumberOfTuples(1);
      vtkDistanceNextKey->SetTuple1(
        0, distancesToKeyFrames_[cptInterpolatedTree * 2 + 1]);
      vtkOutputNode1->GetFieldData()->AddArray(vtkDistanceNextKey);
      // Index previous key frame
      vtkNew<vtkIntArray> vtkIndexPreviousKey{};
      vtkIndexPreviousKey->SetName("IndexPreviousKeyFrame");
      vtkIndexPreviousKey->SetNumberOfTuples(1);
      vtkIndexPreviousKey->SetTuple1(
        0, std::get<3>(coefs[cptInterpolatedTree]));
      vtkOutputNode1->GetFieldData()->AddArray(vtkIndexPreviousKey);
      // Index previous key frame
      vtkNew<vtkIntArray> vtkIndexNextKey{};
      vtkIndexNextKey->SetName("IndexNextKeyFrame");
      vtkIndexNextKey->SetNumberOfTuples(1);
      vtkIndexNextKey->SetTuple1(0, std::get<4>(coefs[cptInterpolatedTree]));
      vtkOutputNode1->GetFieldData()->AddArray(vtkIndexNextKey);
      // Increment interpolated tree counter
      ++cptInterpolatedTree;
    }

    // Construct multiblock
    vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(0))
      ->SetBlock(i, vtkOutputNode1);
    vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(1))
      ->SetBlock(i, vtkOutputArc1);
    if(OutputSegmentation)
      vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(2))
        ->SetBlock(i, vtkOutputSegmentation1);
  }

  // -----------------------------------------
  // --- Matching
  // -----------------------------------------
  output_matchings->SetNumberOfBlocks(intermediateMTrees.size() - 1);
  for(size_t i = 0; i < intermediateMTrees.size() - 1; ++i) {
    // Declare vtk objects
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputMatching
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode1
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(0))
          ->GetBlock(i));
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode2
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(0))
          ->GetBlock(i + 1));
    std::vector<std::vector<SimplexId>> nodeCorrTemp;
    nodeCorrTemp.emplace_back(nodeCorr[i]);
    nodeCorrTemp.emplace_back(nodeCorr[i + 1]);

    // Fill vtk objects
    ttkMergeTreeVisualization visuMakerMatching;
    visuMakerMatching.setVtkOutputMatching(vtkOutputMatching);
    visuMakerMatching.setOutputMatching(allMatching[i]);
    visuMakerMatching.setVtkOutputNode1(vtkOutputNode2);
    visuMakerMatching.setVtkOutputNode2(vtkOutputNode1);
    visuMakerMatching.setNodeCorr1(nodeCorrTemp);
    visuMakerMatching.setDebugLevel(this->debugLevel_);

    visuMakerMatching.makeMatchingOutput<dataType>(trees[i], trees[i + 1]);

    // Field data
    vtkNew<vtkDoubleArray> vtkDistance{};
    vtkDistance->SetName("Distance");
    vtkDistance->SetNumberOfTuples(1);
    vtkDistance->SetTuple1(0, finalDistances_[i]);
    vtkOutputMatching->GetFieldData()->AddArray(vtkDistance);

    // Construct multiblock
    output_matchings->SetBlock(i, vtkOutputMatching);
  }

  // -----------------------------------------
  treesNodes.swap(treesNodesT);
  treesArcs.swap(treesArcsT);
  treesSegmentation.swap(treesSegmentationT);
  treesNodeCorrMesh.swap(treesNodeCorrMeshT);
  inputTrees.swap(inputTreesT);

  return 1;
}
