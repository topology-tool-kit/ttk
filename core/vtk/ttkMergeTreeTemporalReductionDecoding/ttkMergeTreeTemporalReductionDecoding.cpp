#include <ttkMergeTreeTemporalReductionDecoding.h>

#include <ttkFTMTreeUtils.h>
#include <ttkMergeTreeVisualization.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkMacros.h>
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

ttkMergeTreeTemporalReductionDecoding::
  ~ttkMergeTreeTemporalReductionDecoding() {
}

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
  std::vector<vtkMultiBlockDataSet *> inputTrees;
  loadBlocks(inputTrees, blocks);
  printMsg("Load Blocks done.", debug::Priority::VERBOSE);
  int dataTypeInt
    = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
        ->GetPointData()
        ->GetArray("Scalar")
        ->GetDataType();
  auto assignmentSolverArray
    = blocks->GetFieldData()->GetArray("AssignmentSolver");
  if(assignmentSolverArray)
    assignmentSolverID_ = assignmentSolverArray->GetTuple1(0);

  // If we have already computed once but the input has changed
  if(treesNodes.size() != 0 and inputTrees[0]->GetBlock(0) != treesNodes[0])
    resetDataVisualization();

  std::vector<std::tuple<double, int, int, int, int>> coefs;
  int totalInputs = inputTrees.size() + table->GetNumberOfRows();
  std::vector<bool> interpolatedTrees(totalInputs, false);
  for(int i = 0; i < table->GetNumberOfRows(); ++i) {
    double alpha = vtkDataArray::SafeDownCast(table->GetColumnByName("Alpha"))
                     ->GetTuple1(i);
    int index1R = vtkDataArray::SafeDownCast(table->GetColumnByName("Index1_R"))
                    ->GetTuple1(i);
    int index2R = vtkDataArray::SafeDownCast(table->GetColumnByName("Index2_R"))
                    ->GetTuple1(i);
    int index1 = vtkDataArray::SafeDownCast(table->GetColumnByName("Index1"))
                   ->GetTuple1(i);
    int index2 = vtkDataArray::SafeDownCast(table->GetColumnByName("Index2"))
                   ->GetTuple1(i);

    coefs.push_back(std::make_tuple(alpha, index1R, index2R, index1, index2));
    for(int j = index1 + 1; j < index2; ++j)
      interpolatedTrees[j] = true;
  }

  int res = 0;
  switch(dataTypeInt) {
    vtkTemplateMacro(
      res = run<VTK_TT>(outputVector, inputTrees, coefs, interpolatedTrees););
  }
  return res;
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::run(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees,
  std::vector<std::tuple<double, int, int, int, int>> &coefs,
  std::vector<bool> &interpolatedTrees) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(inputTrees, coefs);
  runOutput<dataType>(outputVector, inputTrees, coefs, interpolatedTrees);
  return 1;
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::runCompute(
  std::vector<vtkMultiBlockDataSet *> &inputTrees,
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
    intermediateSTrees.push_back(tree);
  }

  return 1;
}

template <class dataType>
int ttkMergeTreeTemporalReductionDecoding::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees,
  std::vector<std::tuple<double, int, int, int, int>> &coefs,
  std::vector<bool> &interpolatedTrees) {
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
  std::vector<vtkMultiBlockDataSet *> inputTreesT;
  int index = 0;
  size_t cpt = 0;
  while(cpt < coefs.size()) {
    while(cpt < coefs.size() and std::get<2>(coefs[cpt]) <= index) {
      treesNodesT.push_back(nullptr);
      treesArcsT.push_back(nullptr);
      treesSegmentationT.push_back(nullptr);
      treesNodeCorrMeshT.push_back(std::vector<int>());
      inputTreesT.push_back(nullptr);
      ++cpt;
    }
    treesNodesT.push_back(treesNodes[index]);
    treesArcsT.push_back(treesArcs[index]);
    treesSegmentationT.push_back(treesSegmentation[index]);
    treesNodeCorrMeshT.push_back(treesNodeCorrMesh[index]);
    inputTreesT.push_back(inputTrees[index]);
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
  output_sequence->SetNumberOfBlocks(intermediateMTrees.size());
  double prevXMax = 0;
  bool OutputSegmentation = false;
  std::vector<std::vector<SimplexId>> nodeCorr(intermediateMTrees.size());
  int cptInterpolatedTree = 0;
  for(size_t i = 0; i < intermediateMTrees.size(); ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkMultiBlockDataSet> vtkBlock1
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();

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
    // visuMaker.setShiftMode(3); // Double Line
    visuMaker.setShiftMode(2); // Line
    visuMaker.setVtkOutputNode(vtkOutputNode1);
    visuMaker.setVtkOutputArc(vtkOutputArc1);
    visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation1);
    visuMaker.setTreesNodes(treesNodes);
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
    if(inputTrees[i] != nullptr)
      vtkBlock1->GetFieldData()->ShallowCopy(inputTrees[i]->GetFieldData());

    // Construct multiblock
    vtkBlock1->SetNumberOfBlocks((OutputSegmentation ? 3 : 2));
    vtkBlock1->SetBlock(0, vtkOutputNode1);
    vtkBlock1->SetBlock(1, vtkOutputArc1);
    if(OutputSegmentation)
      vtkBlock1->SetBlock(2, vtkOutputSegmentation1);

    // Field data
    if(interpolatedTrees[i]) {
      // Distance previous key frame
      vtkNew<vtkDoubleArray> vtkDistancePreviousKey{};
      vtkDistancePreviousKey->SetName("DistancePreviousKeyFrame");
      vtkDistancePreviousKey->SetNumberOfTuples(1);
      vtkDistancePreviousKey->SetTuple1(
        0, distancesToKeyFrames_[cptInterpolatedTree * 2]);
      vtkBlock1->GetFieldData()->AddArray(vtkDistancePreviousKey);
      // Distance next key frame
      vtkNew<vtkDoubleArray> vtkDistanceNextKey{};
      vtkDistanceNextKey->SetName("DistanceNextKeyFrame");
      vtkDistanceNextKey->SetNumberOfTuples(1);
      vtkDistanceNextKey->SetTuple1(
        0, distancesToKeyFrames_[cptInterpolatedTree * 2 + 1]);
      vtkBlock1->GetFieldData()->AddArray(vtkDistanceNextKey);
      // Index previous key frame
      vtkNew<vtkIntArray> vtkIndexPreviousKey{};
      vtkIndexPreviousKey->SetName("IndexPreviousKeyFrame");
      vtkIndexPreviousKey->SetNumberOfTuples(1);
      vtkIndexPreviousKey->SetTuple1(
        0, std::get<3>(coefs[cptInterpolatedTree]));
      vtkBlock1->GetFieldData()->AddArray(vtkIndexPreviousKey);
      // Index previous key frame
      vtkNew<vtkIntArray> vtkIndexNextKey{};
      vtkIndexNextKey->SetName("IndexNextKeyFrame");
      vtkIndexNextKey->SetNumberOfTuples(1);
      vtkIndexNextKey->SetTuple1(0, std::get<4>(coefs[cptInterpolatedTree]));
      vtkBlock1->GetFieldData()->AddArray(vtkIndexNextKey);
      // Increment interpolated tree counter
      ++cptInterpolatedTree;
    }

    // Add to all multiblocks
    output_sequence->SetBlock(i, vtkBlock1);
  }

  // -----------------------------------------
  // --- Matching
  // -----------------------------------------
  output_matchings->SetNumberOfBlocks(intermediateMTrees.size() - 1);
  for(size_t i = 0; i < intermediateMTrees.size() - 1; ++i) {
    // Declare vtk objects
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputMatching
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(i))
          ->GetBlock(0));
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode2
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(output_sequence->GetBlock(i + 1))
          ->GetBlock(0));
    std::vector<std::vector<SimplexId>> nodeCorrTemp;
    nodeCorrTemp.push_back(nodeCorr[i]);
    nodeCorrTemp.push_back(nodeCorr[i + 1]);

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
