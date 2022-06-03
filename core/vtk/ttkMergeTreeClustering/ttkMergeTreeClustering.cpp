#include <FTMStructures.h>
#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <MergeTreeUtils.h>
#include <MergeTreeVisualization.h>
#include <ttkFTMTreeUtils.h>
#include <ttkMacros.h>
#include <ttkMergeTreeClustering.h>
#include <ttkMergeTreeVisualization.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <string>

using namespace ttk;
using namespace ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeClustering);

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
ttkMergeTreeClustering::ttkMergeTreeClustering() {
  this->setDebugMsgPrefix("MergeTreeClustering");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(3);
}

ttkMergeTreeClustering::~ttkMergeTreeClustering() = default;

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeClustering::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
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
 *      info->Set( ttkAlgorithm::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMergeTreeClustering::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
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
int ttkMergeTreeClustering::RequestData(vtkInformation *ttkNotUsed(request),
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
  loadBlocks(inputTrees, blocks);
  loadBlocks(inputTrees2, blocks2);
  int dataTypeInt
    = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
        ->GetPointData()
        ->GetArray("Scalar")
        ->GetDataType();

  // If we have already computed once but the input has changed
  if(treesNodes.size() != 0 and inputTrees[0]->GetBlock(0) != treesNodes[0])
    resetDataVisualization();

  int res = 0;
  switch(dataTypeInt) {
    vtkTemplateMacro(res = run<VTK_TT>(outputVector, inputTrees, inputTrees2););
  }
  return res;
}

template <class dataType>
int ttkMergeTreeClustering::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(outputVector, inputTrees, inputTrees2);
  runOutput<dataType>(outputVector, inputTrees, inputTrees2);
  return 1;
}

template <class dataType>
int ttkMergeTreeClustering::runCompute(
  vtkInformationVector *ttkNotUsed(outputVector),
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------

  const int numInputs = inputTrees.size();
  const int numInputs2 = inputTrees2.size();

  setDataVisualization(numInputs, numInputs2);

  std::vector<MergeTree<dataType>> intermediateMTrees(numInputs),
    intermediateMTrees2(numInputs2), barycenters(NumberOfBarycenters);
  std::vector<FTMTree_MT *> intermediateTrees(numInputs),
    intermediateTrees2(numInputs2);

  constructTrees<dataType>(
    inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation);
  constructTrees<dataType>(inputTrees2, intermediateMTrees2, treesNodes2,
                           treesArcs2, treesSegmentation2);

  mergeTreeToFTMTree<dataType>(intermediateMTrees, intermediateTrees);
  mergeTreeToFTMTree<dataType>(intermediateMTrees2, intermediateTrees2);

  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  bool AddNodes = true;
  // Classical distance
  double distance = 0;

  // Verify parameters
  if(Backend == 0) {
    BranchDecomposition = true;
    NormalizedWasserstein = true;
    KeepSubtree = false;
  } else if(Backend == 1) {
    BranchDecomposition = false;
    NormalizedWasserstein = false;
    KeepSubtree = true;
    ComputeBarycenter = false;
  }
  if(ComputeBarycenter) {
    if(not BranchDecomposition)
      printMsg("BranchDecomposition is set to true since the barycenter "
               "computation is asked.");
    BranchDecomposition = true;
    if(KeepSubtree)
      printMsg("KeepSubtree is set to false since the barycenter computation "
               "is asked.");
    KeepSubtree = false;
  }
  if(not BranchDecomposition) {
    if(NormalizedWasserstein)
      printMsg("NormalizedWasserstein is set to false since branch "
               "decomposition is not asked.");
    NormalizedWasserstein = false;
  }
  EpsilonTree2 = EpsilonTree1;
  Epsilon2Tree2 = Epsilon2Tree1;
  Epsilon3Tree2 = Epsilon3Tree1;
  printMsg("BranchDecomposition: " + std::to_string(BranchDecomposition));
  printMsg("NormalizedWasserstein: " + std::to_string(NormalizedWasserstein));
  printMsg("KeepSubtree: " + std::to_string(KeepSubtree));

  // Call base
  if(not ComputeBarycenter) {
    MergeTreeDistance mergeTreeDistance;
    mergeTreeDistance.setAssignmentSolver(AssignmentSolver);
    mergeTreeDistance.setEpsilonTree1(EpsilonTree1);
    mergeTreeDistance.setEpsilonTree2(EpsilonTree2);
    mergeTreeDistance.setEpsilon2Tree1(Epsilon2Tree1);
    mergeTreeDistance.setEpsilon2Tree2(Epsilon2Tree2);
    mergeTreeDistance.setEpsilon3Tree1(Epsilon3Tree1);
    mergeTreeDistance.setEpsilon3Tree2(Epsilon3Tree2);
    mergeTreeDistance.setProgressiveComputation(ProgressiveComputation);
    mergeTreeDistance.setBranchDecomposition(BranchDecomposition);
    mergeTreeDistance.setPersistenceThreshold(PersistenceThreshold);
    mergeTreeDistance.setNormalizedWasserstein(NormalizedWasserstein);
    mergeTreeDistance.setNormalizedWassersteinReg(NormalizedWassersteinReg);
    mergeTreeDistance.setRescaledWasserstein(RescaledWasserstein);
    mergeTreeDistance.setKeepSubtree(KeepSubtree);
    mergeTreeDistance.setUseMinMaxPair(UseMinMaxPair);
    mergeTreeDistance.setCleanTree(true);
    mergeTreeDistance.setPostprocess(OutputTrees);
    mergeTreeDistance.setDeleteMultiPersPairs(DeleteMultiPersPairs);
    mergeTreeDistance.setEpsilon1UseFarthestSaddle(Epsilon1UseFarthestSaddle);
    mergeTreeDistance.setThreadNumber(this->threadNumber_);
    mergeTreeDistance.setDebugLevel(this->debugLevel_);

    distance = mergeTreeDistance.execute<dataType>(
      intermediateMTrees[0], intermediateMTrees[1], outputMatching);
    trees1NodeCorrMesh = mergeTreeDistance.getTreesNodeCorr();
    finalDistances = std::vector<double>{distance};
  } else {
    if(NumberOfBarycenters == 1) { // and numInputs2==0){
      MergeTreeBarycenter mergeTreeBarycenter;
      mergeTreeBarycenter.setAssignmentSolver(AssignmentSolver);
      mergeTreeBarycenter.setEpsilonTree1(EpsilonTree1);
      mergeTreeBarycenter.setEpsilonTree2(EpsilonTree2);
      mergeTreeBarycenter.setEpsilon2Tree1(Epsilon2Tree1);
      mergeTreeBarycenter.setEpsilon2Tree2(Epsilon2Tree2);
      mergeTreeBarycenter.setEpsilon3Tree1(Epsilon3Tree1);
      mergeTreeBarycenter.setEpsilon3Tree2(Epsilon3Tree2);
      mergeTreeBarycenter.setProgressiveComputation(ProgressiveComputation);
      mergeTreeBarycenter.setBranchDecomposition(BranchDecomposition);
      mergeTreeBarycenter.setPersistenceThreshold(PersistenceThreshold);
      mergeTreeBarycenter.setNormalizedWasserstein(NormalizedWasserstein);
      mergeTreeBarycenter.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      mergeTreeBarycenter.setRescaledWasserstein(RescaledWasserstein);
      mergeTreeBarycenter.setKeepSubtree(KeepSubtree);
      mergeTreeBarycenter.setUseMinMaxPair(UseMinMaxPair);
      mergeTreeBarycenter.setTol(Tol);
      mergeTreeBarycenter.setAddNodes(AddNodes);
      mergeTreeBarycenter.setDeterministic(Deterministic);
      mergeTreeBarycenter.setBarycenterSizeLimitPercent(
        BarycenterSizeLimitPercent);
      mergeTreeBarycenter.setProgressiveBarycenter(ProgressiveBarycenter);
      mergeTreeBarycenter.setProgressiveSpeedDivisor(ProgressiveSpeedDivisor);
      mergeTreeBarycenter.setAlpha(Alpha);
      mergeTreeBarycenter.setPostprocess(OutputTrees);
      mergeTreeBarycenter.setDeleteMultiPersPairs(DeleteMultiPersPairs);
      mergeTreeBarycenter.setEpsilon1UseFarthestSaddle(
        Epsilon1UseFarthestSaddle);
      mergeTreeBarycenter.setThreadNumber(this->threadNumber_);
      mergeTreeBarycenter.setDebugLevel(this->debugLevel_);

      mergeTreeBarycenter.execute<dataType>(
        intermediateMTrees, outputMatchingBarycenter[0], barycenters[0]);
      trees1NodeCorrMesh = mergeTreeBarycenter.getTreesNodeCorr();
      finalDistances = mergeTreeBarycenter.getFinalDistances();
    } else {
      MergeTreeClustering<dataType> mergeTreeClustering;
      mergeTreeClustering.setAssignmentSolver(AssignmentSolver);
      mergeTreeClustering.setEpsilonTree1(EpsilonTree1);
      mergeTreeClustering.setEpsilonTree2(EpsilonTree2);
      mergeTreeClustering.setEpsilon2Tree1(Epsilon2Tree1);
      mergeTreeClustering.setEpsilon2Tree2(Epsilon2Tree2);
      mergeTreeClustering.setEpsilon3Tree1(Epsilon3Tree1);
      mergeTreeClustering.setEpsilon3Tree2(Epsilon3Tree2);
      mergeTreeClustering.setProgressiveComputation(ProgressiveComputation);
      mergeTreeClustering.setBranchDecomposition(BranchDecomposition);
      mergeTreeClustering.setPersistenceThreshold(PersistenceThreshold);
      mergeTreeClustering.setNormalizedWasserstein(NormalizedWasserstein);
      mergeTreeClustering.setNormalizedWassersteinReg(NormalizedWassersteinReg);
      mergeTreeClustering.setRescaledWasserstein(RescaledWasserstein);
      mergeTreeClustering.setKeepSubtree(KeepSubtree);
      mergeTreeClustering.setUseMinMaxPair(UseMinMaxPair);
      mergeTreeClustering.setTol(Tol);
      mergeTreeClustering.setAddNodes(AddNodes);
      mergeTreeClustering.setDeterministic(Deterministic);
      mergeTreeClustering.setNoCentroids(NumberOfBarycenters);
      mergeTreeClustering.setBarycenterSizeLimitPercent(
        BarycenterSizeLimitPercent);
      mergeTreeClustering.setProgressiveBarycenter(ProgressiveBarycenter);
      mergeTreeClustering.setProgressiveSpeedDivisor(ProgressiveSpeedDivisor);
      mergeTreeClustering.setPostprocess(OutputTrees);
      mergeTreeClustering.setDeleteMultiPersPairs(DeleteMultiPersPairs);
      mergeTreeClustering.setEpsilon1UseFarthestSaddle(
        Epsilon1UseFarthestSaddle);
      mergeTreeClustering.setMixtureCoefficient(JoinSplitMixtureCoefficient);
      mergeTreeClustering.setThreadNumber(this->threadNumber_);
      mergeTreeClustering.setDebugLevel(this->debugLevel_);

      mergeTreeClustering.template execute<dataType>(
        intermediateMTrees, outputMatchingBarycenter, clusteringAssignment,
        intermediateMTrees2, outputMatchingBarycenter2, barycenters);
      trees1NodeCorrMesh = mergeTreeClustering.getTreesNodeCorr();
      trees2NodeCorrMesh = mergeTreeClustering.getTrees2NodeCorr();
      finalDistances = mergeTreeClustering.getFinalDistances();
    }
  }

  mergeTreesTemplateToDouble<dataType>(intermediateMTrees, intermediateSTrees);
  if(numInputs2 != 0)
    mergeTreesTemplateToDouble<dataType>(
      intermediateMTrees2, intermediateSTrees2);
  if(ComputeBarycenter)
    mergeTreesTemplateToDouble<dataType>(barycenters, barycentersS);

  return 1;
}

template <class dataType>
int ttkMergeTreeClustering::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &ttkNotUsed(inputTrees2)) {
  std::vector<MergeTree<dataType>> intermediateMTrees;
  mergeTreesDoubleToTemplate<dataType>(intermediateSTrees, intermediateMTrees);
  std::vector<FTMTree_MT *> intermediateTrees;
  mergeTreeToFTMTree<dataType>(intermediateMTrees, intermediateTrees);

  std::vector<MergeTree<dataType>> barycenters;
  if(ComputeBarycenter)
    mergeTreesDoubleToTemplate<dataType>(barycentersS, barycenters);

  const int numInputs = inputTrees.size();
  // const int numInputs2 = inputTrees2.size();
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  Timer t_makeTreesOutput;

  auto output_clusters = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_centroids = vtkMultiBlockDataSet::GetData(outputVector, 1);
  auto output_matchings = vtkMultiBlockDataSet::GetData(outputVector, 2);

  // Declare internal arrays
  std::vector<std::vector<SimplexId>> nodeCorr(numInputs);

  if(not ComputeBarycenter) {
    // ---------------------------------------------------------------
    // - Classical output
    // ---------------------------------------------------------------
    if(OutputTrees) {
      FTMTree_MT *tree1 = intermediateTrees[0], *tree2 = intermediateTrees[1];
      // ------------------------------------------
      // --- Input trees
      // ------------------------------------------
      // Declare vtk objects
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
        = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1
        = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1
        = vtkSmartPointer<vtkUnstructuredGrid>::New();

      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode2
        = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc2
        = vtkSmartPointer<vtkUnstructuredGrid>::New();
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation2
        = vtkSmartPointer<vtkUnstructuredGrid>::New();

      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockNodes
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockArcs
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockSegs
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();

      // Fill vtk objects
      ttkMergeTreeVisualization visuMaker;
      visuMaker.setPlanarLayout(PlanarLayout);
      visuMaker.setBranchDecompositionPlanarLayout(
        BranchDecompositionPlanarLayout);
      visuMaker.setBranchSpacing(BranchSpacing);
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

      nodeCorr.clear();
      // First tree
      visuMaker.setNoSampleOffset(1);
      visuMaker.setTreesNodes(treesNodes[0]);
      visuMaker.copyPointData(treesNodes[0], trees1NodeCorrMesh[0]);
      visuMaker.setTreesNodeCorrMesh(trees1NodeCorrMesh[0]);
      visuMaker.setTreesSegmentation(treesSegmentation[0]);
      visuMaker.setVtkOutputNode(vtkOutputNode1);
      visuMaker.setVtkOutputArc(vtkOutputArc1);
      visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation1);
      visuMaker.setDebugLevel(this->debugLevel_);

      visuMaker.makeTreesOutput<dataType>(tree1);
      nodeCorr.push_back(visuMaker.getNodeCorr()[0]);

      // Second tree
      visuMaker.setISampleOffset(1);
      visuMaker.setTreesNodes(treesNodes[1]);
      visuMaker.clearAllCustomArrays();
      visuMaker.copyPointData(treesNodes[1], trees1NodeCorrMesh[1]);
      visuMaker.setTreesNodeCorrMesh(trees1NodeCorrMesh[1]);
      visuMaker.setTreesSegmentation(treesSegmentation[1]);
      visuMaker.setVtkOutputNode(vtkOutputNode2);
      visuMaker.setVtkOutputArc(vtkOutputArc2);
      visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation2);
      visuMaker.setDebugLevel(this->debugLevel_);

      visuMaker.makeTreesOutput<dataType>(tree2);
      nodeCorr.push_back(visuMaker.getNodeCorr()[0]);

      // Field data
      vtkOutputNode1->GetFieldData()->ShallowCopy(
        treesNodes[0]->GetFieldData());
      vtkOutputNode2->GetFieldData()->ShallowCopy(
        treesNodes[1]->GetFieldData());
      vtkOutputArc1->GetFieldData()->ShallowCopy(treesArcs[0]->GetFieldData());
      vtkOutputArc2->GetFieldData()->ShallowCopy(treesArcs[1]->GetFieldData());
      if(OutputSegmentation) {
        vtkOutputSegmentation1->GetFieldData()->ShallowCopy(
          treesSegmentation[0]->GetFieldData());
        vtkOutputSegmentation2->GetFieldData()->ShallowCopy(
          treesSegmentation[1]->GetFieldData());
      }

      // Construct multiblock
      vtkBlockNodes->SetNumberOfBlocks(2);
      vtkBlockNodes->SetBlock(0, vtkOutputNode1);
      vtkBlockNodes->SetBlock(1, vtkOutputNode2);

      vtkBlockArcs->SetNumberOfBlocks(2);
      vtkBlockArcs->SetBlock(0, vtkOutputArc1);
      vtkBlockArcs->SetBlock(1, vtkOutputArc2);

      output_clusters->SetNumberOfBlocks((OutputSegmentation ? 3 : 2));
      output_clusters->SetBlock(0, vtkBlockNodes);
      output_clusters->SetBlock(1, vtkBlockArcs);
      if(OutputSegmentation) {
        vtkBlockSegs->SetNumberOfBlocks(2);
        vtkBlockSegs->SetBlock(0, vtkOutputSegmentation1);
        vtkBlockSegs->SetBlock(1, vtkOutputSegmentation2);
        output_clusters->SetBlock(2, vtkBlockSegs);
      }

      // ------------------------------------------
      // --- Matching
      // ------------------------------------------
      // Declare vtk objects
      vtkSmartPointer<vtkUnstructuredGrid> vtkOutputMatching
        = vtkSmartPointer<vtkUnstructuredGrid>::New();

      // Fill vtk objects
      ttkMergeTreeVisualization visuMakerMatching;
      visuMakerMatching.setVtkOutputMatching(vtkOutputMatching);
      visuMakerMatching.setOutputMatching(outputMatching);
      visuMakerMatching.setVtkOutputNode2(vtkOutputNode1);
      visuMakerMatching.setVtkOutputNode1(vtkOutputNode2);
      visuMakerMatching.setNodeCorr1(nodeCorr);
      visuMakerMatching.setDebugLevel(this->debugLevel_);

      visuMakerMatching.makeMatchingOutput<dataType>(tree1, tree2);

      // Field data
      vtkNew<vtkDoubleArray> vtkDistance{};
      vtkDistance->SetName("Distance");
      vtkDistance->SetNumberOfTuples(1);
      vtkDistance->SetTuple1(0, finalDistances[0]);
      vtkOutputMatching->GetFieldData()->AddArray(vtkDistance);

      // Construct multiblock
      output_matchings->SetNumberOfBlocks(1);
      output_matchings->SetBlock(0, vtkOutputMatching);
    }
  } else {
    // ---------------------------------------------------------------
    // - Barycenter/Clustering output
    // ---------------------------------------------------------------
    if(OutputTrees) {
      // --- Declare internal arrays
      std::vector<std::vector<SimplexId>> nodeCorrBary(NumberOfBarycenters);
      std::vector<std::vector<float>> allBaryPercentMatch(NumberOfBarycenters);
      std::vector<FTMTree_MT *> barycentersTree;
      mergeTreeToFTMTree<dataType>(barycenters, barycentersTree);

      // ------------------------------------------
      // --- Input trees
      // ------------------------------------------
      output_clusters->SetNumberOfBlocks((OutputSegmentation ? 3 : 2));
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockNodes
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkBlockNodes->SetNumberOfBlocks(numInputs);
      output_clusters->SetBlock(0, vtkBlockNodes);
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockArcs
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkBlockArcs->SetNumberOfBlocks(numInputs);
      output_clusters->SetBlock(1, vtkBlockArcs);
      if(OutputSegmentation) {
        vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockSegs
          = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        vtkBlockSegs->SetNumberOfBlocks(numInputs);
        output_clusters->SetBlock(2, vtkBlockSegs);
      }
      for(unsigned int c = 0; c < NumberOfBarycenters; ++c) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
        for(int i = 0; i < numInputs; ++i) {
          if(clusteringAssignment[i] != (int)c)
            continue;

          // Declare vtk objects
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
            = vtkSmartPointer<vtkUnstructuredGrid>::New();
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1
            = vtkSmartPointer<vtkUnstructuredGrid>::New();
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1
            = vtkSmartPointer<vtkUnstructuredGrid>::New();

          // Fill vtk objects
          ttkMergeTreeVisualization visuMaker;
          visuMaker.setPlanarLayout(PlanarLayout);
          visuMaker.setBranchDecompositionPlanarLayout(
            BranchDecompositionPlanarLayout);
          visuMaker.setBranchSpacing(BranchSpacing);
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
          visuMaker.setTreesNodes(treesNodes);
          visuMaker.copyPointData(treesNodes[i], trees1NodeCorrMesh[i]);
          visuMaker.setTreesNodeCorrMesh(trees1NodeCorrMesh);
          visuMaker.setTreesSegmentation(treesSegmentation);
          visuMaker.setVtkOutputNode(vtkOutputNode1);
          visuMaker.setVtkOutputArc(vtkOutputArc1);
          visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation1);
          visuMaker.setClusteringAssignment(clusteringAssignment);
          visuMaker.setOutputMatchingBarycenter(outputMatchingBarycenter);
          visuMaker.setPrintTreeId(i);
          visuMaker.setPrintClusterId(c);
          visuMaker.setDebugLevel(this->debugLevel_);

          visuMaker.makeTreesOutput<dataType>(
            intermediateTrees, barycentersTree);
          auto nodeCorrT = visuMaker.getNodeCorr();
          nodeCorr[i] = nodeCorrT[i];

          // Field data
          vtkOutputNode1->GetFieldData()->ShallowCopy(
            treesNodes[i]->GetFieldData());
          vtkOutputArc1->GetFieldData()->ShallowCopy(
            treesArcs[i]->GetFieldData());
          if(OutputSegmentation)
            vtkOutputSegmentation1->GetFieldData()->ShallowCopy(
              treesSegmentation[i]->GetFieldData());

          vtkNew<vtkDoubleArray> vtkClusterAssignment{};
          vtkClusterAssignment->SetName("ClusterAssignment");
          vtkClusterAssignment->SetNumberOfTuples(1);
          vtkClusterAssignment->SetTuple1(0, clusteringAssignment[i]);
          vtkOutputNode1->GetFieldData()->AddArray(vtkClusterAssignment);

          // Construct multiblock
          vtkMultiBlockDataSet::SafeDownCast(output_clusters->GetBlock(0))
            ->SetBlock(i, vtkOutputNode1);
          vtkMultiBlockDataSet::SafeDownCast(output_clusters->GetBlock(1))
            ->SetBlock(i, vtkOutputArc1);
          if(OutputSegmentation)
            vtkMultiBlockDataSet::SafeDownCast(output_clusters->GetBlock(2))
              ->SetBlock(i, vtkOutputSegmentation1);
        }
      }

      // ------------------------------------------
      // --- Barycenter(s)
      // ------------------------------------------
      output_centroids->SetNumberOfBlocks(2);
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockNodes2
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkBlockNodes2->SetNumberOfBlocks(NumberOfBarycenters);
      output_centroids->SetBlock(0, vtkBlockNodes2);
      vtkSmartPointer<vtkMultiBlockDataSet> vtkBlockArcs2
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkBlockArcs2->SetNumberOfBlocks(NumberOfBarycenters);
      output_centroids->SetBlock(1, vtkBlockArcs2);
      for(unsigned int c = 0; c < NumberOfBarycenters; ++c) {

        // Declare vtk objects
        vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode2
          = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc2
          = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation2
          = vtkSmartPointer<vtkUnstructuredGrid>::New();

        // Fill vtk objects
        ttkMergeTreeVisualization visuMakerBary;
        visuMakerBary.setPlanarLayout(PlanarLayout);
        visuMakerBary.setBranchDecompositionPlanarLayout(
          BranchDecompositionPlanarLayout);
        visuMakerBary.setBranchSpacing(BranchSpacing);
        visuMakerBary.setRescaleTreesIndividually(RescaleTreesIndividually);
        visuMakerBary.setOutputSegmentation(false);
        visuMakerBary.setDimensionSpacing(DimensionSpacing);
        visuMakerBary.setDimensionToShift(DimensionToShift);
        visuMakerBary.setImportantPairs(ImportantPairs);
        visuMakerBary.setMaximumImportantPairs(MaximumImportantPairs);
        visuMakerBary.setMinimumImportantPairs(MinimumImportantPairs);
        visuMakerBary.setImportantPairsSpacing(ImportantPairsSpacing);
        visuMakerBary.setNonImportantPairsSpacing(NonImportantPairsSpacing);
        visuMakerBary.setNonImportantPairsProximity(NonImportantPairsProximity);
        visuMakerBary.setExcludeImportantPairsHigher(
          ExcludeImportantPairsHigher);
        visuMakerBary.setExcludeImportantPairsLower(ExcludeImportantPairsLower);
        visuMakerBary.setShiftMode(1); // Star Barycenter
        visuMakerBary.setTreesNodes(treesNodes);
        visuMakerBary.setTreesNodeCorrMesh(trees1NodeCorrMesh);
        visuMakerBary.setTreesSegmentation(treesSegmentation);
        visuMakerBary.setVtkOutputNode(vtkOutputNode2);
        visuMakerBary.setVtkOutputArc(vtkOutputArc2);
        visuMakerBary.setVtkOutputSegmentation(vtkOutputSegmentation2);
        visuMakerBary.setClusteringAssignment(clusteringAssignment);
        visuMakerBary.setOutputMatchingBarycenter(outputMatchingBarycenter);
        visuMakerBary.setPrintTreeId(c);
        visuMakerBary.setPrintClusterId(c);
        if(numInputs == 2 and NumberOfBarycenters == 1) {
          visuMakerBary.setBarycenterPositionAlpha(BarycenterPositionAlpha);
          visuMakerBary.setAlpha(Alpha);
        }
        visuMakerBary.setDebugLevel(this->debugLevel_);

        visuMakerBary.makeTreesOutput<dataType>(
          intermediateTrees, barycentersTree);
        auto nodeCorrBaryT = visuMakerBary.getNodeCorr();
        nodeCorrBary[c] = nodeCorrBaryT[c];
        auto allBaryPercentMatchT = visuMakerBary.getAllBaryPercentMatch();
        allBaryPercentMatch[c] = allBaryPercentMatchT[c];

        // Field data
        vtkNew<vtkDoubleArray> vtkClusterAssignment{};
        vtkClusterAssignment->SetName("ClusterAssignment");
        vtkClusterAssignment->SetNumberOfTuples(1);
        vtkClusterAssignment->SetTuple1(0, c);

        vtkOutputNode2->GetFieldData()->AddArray(vtkClusterAssignment);
        vtkOutputArc2->GetFieldData()->AddArray(vtkClusterAssignment);

        // Construct multiblock
        vtkMultiBlockDataSet::SafeDownCast(output_centroids->GetBlock(0))
          ->SetBlock(c, vtkOutputNode2);
        vtkMultiBlockDataSet::SafeDownCast(output_centroids->GetBlock(1))
          ->SetBlock(c, vtkOutputArc2);
      }

      // ------------------------------------------
      // --- Matching
      // ------------------------------------------
      output_matchings->SetNumberOfBlocks(numInputs);
      for(unsigned int c = 0; c < NumberOfBarycenters; ++c) {
        for(int i = 0; i < numInputs; ++i) {
          if(clusteringAssignment[i] != (int)c)
            continue;

          // Declare vtk objects
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputMatching
            = vtkSmartPointer<vtkUnstructuredGrid>::New();
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
            = vtkUnstructuredGrid::SafeDownCast(
              vtkMultiBlockDataSet::SafeDownCast(output_clusters->GetBlock(0))
                ->GetBlock(i));
          vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode2
            = vtkUnstructuredGrid::SafeDownCast(
              vtkMultiBlockDataSet::SafeDownCast(output_centroids->GetBlock(0))
                ->GetBlock(c));

          // Fill vtk objects
          ttkMergeTreeVisualization visuMakerMatching;
          visuMakerMatching.setVtkOutputMatching(vtkOutputMatching);
          visuMakerMatching.setOutputMatchingBarycenter(
            outputMatchingBarycenter);
          visuMakerMatching.setAllBaryPercentMatch(allBaryPercentMatch);
          visuMakerMatching.setVtkOutputNode1(vtkOutputNode1);
          visuMakerMatching.setVtkOutputNode2(vtkOutputNode2);
          visuMakerMatching.setNodeCorr1(nodeCorr);
          visuMakerMatching.setNodeCorr2(nodeCorrBary);
          visuMakerMatching.setPrintTreeId(i);
          visuMakerMatching.setPrintClusterId(c);
          visuMakerMatching.setDebugLevel(this->debugLevel_);

          visuMakerMatching.makeMatchingOutput<dataType>(
            intermediateTrees, barycentersTree);

          // Field data
          vtkNew<vtkDoubleArray> vtkDistance{};
          vtkDistance->SetName("Distance");
          vtkDistance->SetNumberOfTuples(1);
          vtkDistance->SetTuple1(0, finalDistances[i]);
          vtkOutputMatching->GetFieldData()->AddArray(vtkDistance);

          // Construct multiblock
          output_matchings->SetBlock(i, vtkOutputMatching);
        }
      }
    }
  }

  printMsg(
    "Trees output", 1, t_makeTreesOutput.getElapsedTime(), this->threadNumber_);

  return 1;
}
