#include <ttkPlanarGraphLayout.h>

#include <ttkMacros.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <FTMTreePPUtils.h>
#include <ttkFTMTreeUtils.h>
#include <ttkMergeTreeVisualization.h>

vtkStandardNewMacro(ttkPlanarGraphLayout);

ttkPlanarGraphLayout::ttkPlanarGraphLayout() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkPlanarGraphLayout::~ttkPlanarGraphLayout() {
}

int ttkPlanarGraphLayout::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  } else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::planarGraphLayoutCall(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // Get input and output objects

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  if(!input || !output)
    return !this->printErr("Unable to retrieve input/output data objects.");

  // Copy input to output
  output->ShallowCopy(input);

  size_t nPoints = output->GetNumberOfPoints();
  size_t nEdges = output->GetNumberOfCells();

  // Get input arrays
  auto sequenceArray = this->GetInputArrayToProcess(0, inputVector);
  if(this->GetUseSequences() && !sequenceArray)
    return !this->printErr("Unable to retrieve sequence array.");

  auto sizeArray = this->GetInputArrayToProcess(1, inputVector);
  if(this->GetUseSizes() && !sizeArray)
    return !this->printErr("Unable to retrieve size array.");

  auto branchArray = this->GetInputArrayToProcess(2, inputVector);
  if(this->GetUseBranches() && !branchArray)
    return !this->printErr("Unable to retrieve branch array.");

  auto levelArray = this->GetInputArrayToProcess(3, inputVector);
  if(this->GetUseLevels() && !levelArray)
    return !this->printErr("Unable to retrieve level array.");

  // Initialize output array
  auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
  outputArray->SetName(this->GetOutputArrayName().data());
  outputArray->SetNumberOfComponents(2); // (x,y) position
  outputArray->SetNumberOfValues(nPoints * 2);

  vtkDataArray *cells{nullptr};
  if(auto outputAsUG = vtkUnstructuredGrid::SafeDownCast(output))
    cells = outputAsUG->GetCells()->GetConnectivityArray();
  else if(auto outputAsPD = vtkPolyData::SafeDownCast(output))
    cells = outputAsPD->GetLines()->GetConnectivityArray();

  if(!cells)
    return !this->printErr("Unable to retrieve connectivity array.");

  int status = 1;
  ttkTypeMacroAII(
    this->GetUseSequences() ? sequenceArray->GetDataType() : VTK_INT,
    this->GetUseBranches() ? branchArray->GetDataType() : VTK_INT,
    cells->GetDataType(),
    (status = this->computeLayout<T0, T1, T2>(
       // Output
       ttkUtils::GetPointer<float>(outputArray),
       // Input
       ttkUtils::GetPointer<T2>(cells), nPoints, nEdges,
       this->GetUseSequences() ? ttkUtils::GetPointer<T0>(sequenceArray)
                               : nullptr,
       this->GetUseSizes() ? ttkUtils::GetPointer<float>(sizeArray) : nullptr,
       this->GetUseBranches() ? ttkUtils::GetPointer<T1>(branchArray) : nullptr,
       this->GetUseLevels() ? ttkUtils::GetPointer<T1>(levelArray) : nullptr)));

  if(status != 1)
    return 0;

  // Add output field to output
  output->GetPointData()->AddArray(outputArray);

  return 1;
}

template <class dataType>
int ttkPlanarGraphLayout::mergeTreePlanarLayoutCallTemplate(
  vtkUnstructuredGrid *treeNodes,
  vtkUnstructuredGrid *treeArcs,
  vtkUnstructuredGrid *output) {
  auto mergeTree = ttk::ftm::makeTree<dataType>(treeNodes, treeArcs);
  ttk::ftm::FTMTree_MT *tree = &(mergeTree.tree);
  // tree->printTree2();

  ttk::ftm::computePersistencePairs<dataType>(tree);

  std::vector<std::vector<int>> treeNodeCorrMesh(1);
  treeNodeCorrMesh[0] = std::vector<int>(tree->getNumberOfNodes());
  for(unsigned int j = 0; j < tree->getNumberOfNodes(); ++j)
    treeNodeCorrMesh[0][j] = j;

  ttkMergeTreeVisualization visuMaker;
  visuMaker.setPlanarLayout(true);
  visuMaker.setOutputSegmentation(false);
  visuMaker.setBranchDecompositionPlanarLayout(BranchDecompositionPlanarLayout);
  visuMaker.setBranchSpacing(BranchSpacing);
  visuMaker.setImportantPairs(ImportantPairs);
  visuMaker.setMaximumImportantPairs(MaximumImportantPairs);
  visuMaker.setMinimumImportantPairs(MinimumImportantPairs);
  visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
  visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
  visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
  visuMaker.setVtkOutputNode(output);
  visuMaker.setVtkOutputArc(output);
  visuMaker.setTreesNodes(treeNodes);
  visuMaker.setTreesNodeCorrMesh(treeNodeCorrMesh);
  visuMaker.setDebugLevel(this->debugLevel_);
  visuMaker.makeTreesOutput<dataType>(tree);

  return 1;
}

int ttkPlanarGraphLayout::mergeTreePlanarLayoutCall(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkUnstructuredGrid *treeNodes
    = vtkUnstructuredGrid::GetData(inputVector[0], 0);
  vtkUnstructuredGrid *treeArcs
    = vtkUnstructuredGrid::GetData(inputVector[0], 1);
  auto output = vtkUnstructuredGrid::GetData(outputVector);
  int dataTypeInt
    = treeNodes->GetPointData()->GetArray("Scalar")->GetDataType();

  int res = 0;

  switch(dataTypeInt) {
    vtkTemplateMacro(res = mergeTreePlanarLayoutCallTemplate<VTK_TT>(
                       treeNodes, treeArcs, output););
  }

  return res;
}

int ttkPlanarGraphLayout::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  if(not InputIsAMergeTree)
    return planarGraphLayoutCall(request, inputVector, outputVector);
  else
    return mergeTreePlanarLayoutCall(request, inputVector, outputVector);
}
