#include <ttkMergeTreeTemporalReductionEncoding.h>

#include <ttkFTMTreeUtils.h>
#include <ttkMergeTreeVisualization.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkResampleToImage.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace ttk;
using namespace ttk::ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeTemporalReductionEncoding);

/**
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMergeTreeTemporalReductionEncoding::ttkMergeTreeTemporalReductionEncoding() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

ttkMergeTreeTemporalReductionEncoding::
  ~ttkMergeTreeTemporalReductionEncoding() {
}

/**
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeTemporalReductionEncoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

/**
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
int ttkMergeTreeTemporalReductionEncoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  } else
    return 0;
  return 1;
}

/**
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
int ttkMergeTreeTemporalReductionEncoding::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // ------------------------------------------------------------------------------------
  // --- Get input object from input vector
  // ------------------------------------------------------------------------------------
  printMsg("Get input object from input vector", debug::Priority::VERBOSE);
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);

  // ------------------------------------------------------------------------------------
  // --- Load blocks
  // ------------------------------------------------------------------------------------
  printMsg("Load blocks", ttk::debug::Priority::VERBOSE);
  std::vector<vtkMultiBlockDataSet *> inputTrees;
  loadBlocks(inputTrees, blocks);
  printMsg("Load blocks done.", debug::Priority::VERBOSE);
  int dataTypeInt
    = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
        ->GetPointData()
        ->GetArray("Scalar")
        ->GetDataType();

  // If we have already computed once but the input has changed
  if(treesNodes.size() != 0 and inputTrees[0]->GetBlock(0) != treesNodes[0])
    resetDataVisualization();

  // Run
  int res = 0;
  switch(dataTypeInt) {
    vtkTemplateMacro(res = run<VTK_TT>(outputVector, inputTrees););
  }
  return res;
}

template <class dataType>
int ttkMergeTreeTemporalReductionEncoding::run(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(inputTrees);
  runOutput<dataType>(outputVector, inputTrees);
  return 1;
}

std::vector<vtkSmartPointer<vtkDataSet>>
  preprocessSegmentation(std::vector<vtkDataSet *> &treesSegmentation,
                         bool doResampleToImage = false) {
  std::vector<vtkSmartPointer<vtkDataSet>> images;
  for(unsigned int i = 0; i < treesSegmentation.size(); ++i) {
    if(doResampleToImage) {
      auto seg = treesSegmentation[i];
      // auto gridSeg = vtkUnstructuredGrid::SafeDownCast(seg);
      auto resampleFilter = vtkSmartPointer<vtkResampleToImage>::New();
      resampleFilter->SetInputDataObject(seg);
      resampleFilter->SetUseInputBounds(true);
      // resampleFilter->SetSamplingDimensions(100, 100, 100);
      resampleFilter->Update();
      auto resampleFilterOut = resampleFilter->GetOutput();
      // images.push_back(resampleFilterOut);
      vtkSmartPointer<vtkDataSet> image = vtkSmartPointer<vtkImageData>::New();
      image->DeepCopy(resampleFilterOut);
      images.push_back(image);
    } else {
      vtkSmartPointer<vtkDataSet> image;
      image.TakeReference(treesSegmentation[i]);
      images.push_back(image);
    }
  }
  return images;
}

template <class dataType>
int ttkMergeTreeTemporalReductionEncoding::runCompute(
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  // ------------------------------------------------------------------------------------
  // --- Construct trees
  // ------------------------------------------------------------------------------------
  printMsg("Construct trees", debug::Priority::VERBOSE);

  const int numInputs = inputTrees.size();

  setDataVisualization(numInputs);

  std::vector<MergeTree<dataType>> intermediateMTrees(numInputs);
  constructTrees<dataType>(
    inputTrees, intermediateMTrees, treesNodes, treesArcs, treesSegmentation);

  // L2 distance preprocessing
  if(useL2Distance_) {
    auto images = preprocessSegmentation(treesSegmentation, DoResampleToImage);
    fieldL2_ = std::vector<std::vector<double>>(images.size());
    for(size_t i = 0; i < images.size(); ++i) {
      auto array = images[i]->GetPointData()->GetArray("Scalars");
      for(vtkIdType j = 0; j < array->GetNumberOfTuples(); ++j)
        fieldL2_[i].push_back(array->GetTuple1(j));
    }
  }

  // Load time variable if needed
  timeVariable_.clear();
  if(useCustomTimeVariable_) {
    timeVariable_ = std::vector<double>(inputTrees.size());
    for(size_t i = 0; i < inputTrees.size(); ++i) {
      for(int j = 0; j < 3; ++j) {
        if(j != 2 or inputTrees[i]->GetNumberOfBlocks() > 2) {
          auto array = inputTrees[i]->GetBlock(j)->GetFieldData()->GetArray(
            TimeVariableName.c_str());
          if(array) {
            timeVariable_[i] = array->GetTuple1(0);
            break;
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // --- Call base
  // ------------------------------------------------------------------------------------
  printMsg("Call base", debug::Priority::VERBOSE);

  std::vector<MergeTree<dataType>> allMT_T;

  removed = execute<dataType>(intermediateMTrees, emptyTreeDistances, allMT_T);
  treesNodeCorrMesh = getTreesNodeCorr();

  int indexRemoved = 0;
  for(size_t i = 0; i < intermediateMTrees.size(); ++i) {
    if((int)i != removed[indexRemoved]) {
      MergeTree<double> keyFrame;
      mergeTreeTemplateToDouble<dataType>(intermediateMTrees[i], keyFrame);
      keyFrames.push_back(keyFrame);
    } else
      ++indexRemoved;
  }

  return 1;
}

template <class dataType>
int ttkMergeTreeTemporalReductionEncoding::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  // ------------------------------------------------------------------------------------
  // --- Create output
  // ------------------------------------------------------------------------------------
  auto output_keyFrames = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_table = vtkTable::GetData(outputVector, 1);

  // ------------------------------------------
  // --- Trees output
  // -----------------------------------------
  std::vector<int> keyFramesIndex;
  int indexRemoved = 0;
  for(size_t i = 0; i < inputTrees.size(); ++i)
    if((int)i != removed[indexRemoved]) {
      keyFramesIndex.push_back(i);
    } else
      ++indexRemoved;

  std::vector<vtkUnstructuredGrid *> treesNodesT;
  std::vector<vtkUnstructuredGrid *> treesArcsT;
  std::vector<vtkDataSet *> treesSegmentationT;
  std::vector<std::vector<int>> treesNodeCorrMeshT;
  std::vector<vtkMultiBlockDataSet *> inputTreesT;
  for(size_t i = 0; i < keyFramesIndex.size(); ++i) {
    treesNodesT.push_back(treesNodes[keyFramesIndex[i]]);
    treesArcsT.push_back(treesArcs[keyFramesIndex[i]]);
    treesSegmentationT.push_back(treesSegmentation[keyFramesIndex[i]]);
    treesNodeCorrMeshT.push_back(treesNodeCorrMesh[keyFramesIndex[i]]);
    inputTreesT.push_back(inputTrees[keyFramesIndex[i]]);
  }
  treesNodes.swap(treesNodesT);
  treesArcs.swap(treesArcsT);
  treesSegmentation.swap(treesSegmentationT);
  treesNodeCorrMesh.swap(treesNodeCorrMeshT);
  inputTrees.swap(inputTreesT);

  std::vector<MergeTree<dataType>> keyFramesT;
  mergeTreesDoubleToTemplate<dataType>(keyFrames, keyFramesT);
  std::vector<FTMTree_MT *> keyFramesTree;
  mergeTreeToFTMTree<dataType>(keyFramesT, keyFramesTree);
  output_keyFrames->SetNumberOfBlocks(keyFramesT.size());
  double prevXMax = 0;
  for(size_t i = 0; i < keyFramesT.size(); ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation1
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkMultiBlockDataSet> vtkBlock1
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    ttkMergeTreeVisualization visuMaker;
    visuMaker.setShiftMode(2); // Line
    visuMaker.setVtkOutputNode(vtkOutputNode1);
    visuMaker.setVtkOutputArc(vtkOutputArc1);
    visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation1);
    visuMaker.setTreesNodes(treesNodes);
    visuMaker.setTreesNodeCorrMesh(treesNodeCorrMesh);
    visuMaker.setTreesSegmentation(treesSegmentation);
    visuMaker.setPrintTreeId(i);
    visuMaker.setPrintClusterId(0);
    visuMaker.setPrevXMaxOffset(prevXMax);
    visuMaker.setDebugLevel(this->debugLevel_);
    visuMaker.makeTreesOutput<dataType>(keyFramesTree);
    prevXMax = visuMaker.getPrevXMax();

    // Field data
    vtkBlock1->GetFieldData()->ShallowCopy(inputTrees[i]->GetFieldData());

    // Construct multiblock
    bool outputSegmentation = inputTrees[i]->GetNumberOfBlocks() == 3;
    vtkBlock1->SetNumberOfBlocks((outputSegmentation ? 3 : 2));
    vtkBlock1->SetBlock(0, vtkOutputNode1);
    vtkBlock1->SetBlock(1, vtkOutputArc1);
    if(outputSegmentation)
      vtkBlock1->SetBlock(2, vtkOutputSegmentation1);

    output_keyFrames->SetBlock(i, vtkBlock1);
  }
  // Add input parameters to field data
  vtkNew<vtkIntArray> vtkAssignmentSolver{};
  vtkAssignmentSolver->SetName("AssignmentSolver");
  vtkAssignmentSolver->InsertNextTuple1(assignmentSolverID_);
  output_keyFrames->GetFieldData()->AddArray(vtkAssignmentSolver);

  // ------------------------------------------
  // --- Table
  // ------------------------------------------
  vtkNew<vtkIntArray> treeIds{};
  treeIds->SetName("treeID");
  treeIds->SetNumberOfTuples(removed.size());
  vtkNew<vtkDoubleArray> coef{};
  coef->SetName("Alpha");
  coef->SetNumberOfTuples(removed.size());
  vtkNew<vtkIntArray> index1Id{};
  index1Id->SetName("Index1");
  index1Id->SetNumberOfTuples(removed.size());
  vtkNew<vtkIntArray> index2Id{};
  index2Id->SetName("Index2");
  index2Id->SetNumberOfTuples(removed.size());
  vtkNew<vtkIntArray> index1IdR{};
  index1IdR->SetName("Index1_R");
  index1IdR->SetNumberOfTuples(removed.size());
  vtkNew<vtkIntArray> index2IdR{};
  index2IdR->SetName("Index2_R");
  index2IdR->SetNumberOfTuples(removed.size());
  int keyFrameCpt = 0;
  for(size_t i = 0; i < removed.size(); ++i) {
    while(removed[i] > keyFramesIndex[keyFrameCpt])
      ++keyFrameCpt;
    int index1 = keyFramesIndex[keyFrameCpt - 1];
    int index2 = keyFramesIndex[keyFrameCpt];
    double alpha = computeAlpha(index1, removed[i], index2);
    coef->SetTuple1(i, alpha);
    index1Id->SetTuple1(i, index1);
    index2Id->SetTuple1(i, index2);
    treeIds->SetTuple1(i, removed[i]);
    index1IdR->SetTuple1(i, keyFrameCpt - 1);
    index2IdR->SetTuple1(i, keyFrameCpt);
  }
  output_table->AddColumn(treeIds);
  output_table->AddColumn(coef);
  output_table->AddColumn(index1Id);
  output_table->AddColumn(index2Id);
  output_table->AddColumn(index1IdR);
  output_table->AddColumn(index2IdR);

  // ------------------------------------------
  treesNodes.swap(treesNodesT);
  treesArcs.swap(treesArcsT);
  treesSegmentation.swap(treesSegmentationT);
  treesNodeCorrMesh.swap(treesNodeCorrMeshT);
  inputTrees.swap(inputTreesT);

  return 1;
}
