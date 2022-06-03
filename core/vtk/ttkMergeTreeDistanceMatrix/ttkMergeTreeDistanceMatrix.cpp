#include <ttkFTMTreeUtils.h>
#include <ttkMacros.h>
#include <ttkMergeTreeDistanceMatrix.h>
#include <ttkUtils.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkTable.h>

using namespace ttk;
using namespace ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeDistanceMatrix);

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
ttkMergeTreeDistanceMatrix::ttkMergeTreeDistanceMatrix() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkMergeTreeDistanceMatrix::~ttkMergeTreeDistanceMatrix() = default;

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeDistanceMatrix::FillInputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    if(port == 1)
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
int ttkMergeTreeDistanceMatrix::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  else
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
template <class dataType>
int ttkMergeTreeDistanceMatrix::run(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2) {

  // Construct trees
  const int numInputs = inputTrees.size();
  std::vector<MergeTree<dataType>> intermediateTrees, intermediateTrees2;
  constructTrees(inputTrees, intermediateTrees);
  constructTrees(inputTrees2, intermediateTrees2);

  // Verify parameters
  if(not UseFieldDataParameters) {
    if(Backend == 0) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = true;
      keepSubtree_ = false;
    } else if(Backend == 1) {
      branchDecomposition_ = false;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
    }
  }
  if(not branchDecomposition_) {
    if(normalizedWasserstein_)
      printMsg("NormalizedWasserstein is set to false since branch "
               "decomposition is not asked.");
    normalizedWasserstein_ = false;
  }
  epsilonTree2_ = epsilonTree1_;
  epsilon2Tree2_ = epsilon2Tree1_;
  epsilon3Tree2_ = epsilon3Tree1_;
  printMsg("BranchDecomposition: " + std::to_string(branchDecomposition_));
  printMsg("NormalizedWasserstein: " + std::to_string(normalizedWasserstein_));
  printMsg("KeepSubtree: " + std::to_string(keepSubtree_));

  // --- Call base
  std::vector<std::vector<double>> treesDistMat(
    numInputs, std::vector<double>(numInputs));
  execute<dataType>(intermediateTrees, intermediateTrees2, treesDistMat);

  // --- Create output
  auto treesDistTable = vtkTable::GetData(outputVector);

  // zero-padd column name to keep Row Data columns ordered
  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  // copy trees distance matrix to output
  vtkNew<vtkIntArray> treeIds{};
  treeIds->SetName("treeID");
  treeIds->SetNumberOfTuples(numInputs);
  for(size_t i = 0; i < treesDistMat.size(); ++i) {
    treeIds->SetTuple1(i, i);

    std::string name{"Tree"};
    zeroPad(name, treesDistMat.size(), i);
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(numInputs);
    col->SetName(name.c_str());
    for(size_t j = 0; j < treesDistMat[i].size(); ++j) {
      col->SetTuple1(j, treesDistMat[i][j]);
    }
    treesDistTable->AddColumn(col);
  }

  treesDistTable->AddColumn(treeIds);

  // aggregate input field data
  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputTrees[0]->GetBlock(0)->GetFieldData());
  fd->SetNumberOfTuples(inputTrees.size());
  for(size_t i = 0; i < inputTrees.size(); ++i) {
    fd->SetTuple(i, 0, inputTrees[i]->GetBlock(0)->GetFieldData());
  }

  // copy input field data to output row data
  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    treesDistTable->AddColumn(fd->GetAbstractArray(i));
  }

  return 1;
}

int ttkMergeTreeDistanceMatrix::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // --- Get input object from input vector
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // --- Load blocks
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTrees, inputTrees2;
  loadBlocks(inputTrees, blocks);
  loadBlocks(inputTrees2, blocks2);

  auto arrayToGet
    = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
        ->GetPointData()
        ->GetArray("Scalar");
  if(arrayToGet == nullptr)
    arrayToGet = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
                   ->GetPointData()
                   ->GetArray("Birth");
  int dataTypeInt = arrayToGet->GetDataType();

  // --- Load field data parameters
  if(UseFieldDataParameters) {
    printMsg("Load parameters from field data.");
    std::vector<std::string> paramNames;
    getParamNames(paramNames);
    for(auto paramName : paramNames) {
      auto array = blocks->GetFieldData()->GetArray(paramName.c_str());
      if(array) {
        double value = array->GetTuple1(0);
        setParamValueFromName(paramName, value);
        printMsg(" - " + paramName + " = " + std::to_string(value));
      } else
        printMsg(" - " + paramName + " was not found in the field data.");
    }
  }

  int res = 0;
  switch(dataTypeInt) {
    vtkTemplateMacro(res = run<VTK_TT>(outputVector, inputTrees, inputTrees2));
  }

  return res;
}
