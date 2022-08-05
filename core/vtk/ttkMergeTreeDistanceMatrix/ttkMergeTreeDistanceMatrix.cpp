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

  // Verify parameters
  if(not UseFieldDataParameters) {
    if(Backend == 0) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = true;
      keepSubtree_ = false;
      baseModule = 0;
    } else if(Backend == 1) {
      branchDecomposition_ = false;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule = 0;
    } else if(Backend == 3) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule = 1;
    } else if(Backend == 4) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule = 2;
    } else {
      baseModule = 0;
    }
  }
  if(baseModule == 0) {
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
    printMsg("NormalizedWasserstein: "
             + std::to_string(normalizedWasserstein_));
    printMsg("KeepSubtree: " + std::to_string(keepSubtree_));
  }
  if(baseModule == 1) {
    printMsg("Using Branch Mapping Distance.");
    std::string metric;
    if(branchMetric == 0)
      metric = "Wasserstein Distance first degree";
    else if(branchMetric == 1)
      metric = "Wasserstein Distance second degree";
    else if(branchMetric == 2)
      metric = "Persistence difference";
    else if(branchMetric == 3)
      metric = "Shifting cost";
    else
      return 1;
    printMsg("BranchMetric: " + metric);
  }
  if(baseModule == 2) {
    printMsg("Using Path Mapping Distance.");
    std::string metric;
    if(pathMetric == 0)
      metric = "Persistence difference";
    else
      return 1;
    printMsg("PathMetric: " + metric);
  }

  // Construct trees
  const int numInputs = inputTrees.size();
  std::vector<MergeTree<dataType>> intermediateTrees, intermediateTrees2;
  std::vector<
    std::tuple<std::vector<dataType>, std::vector<std::vector<int>>, int>>
    adjacenyLists;
  if(baseModule == 0) {
    constructTrees(inputTrees, intermediateTrees);
    constructTrees(inputTrees2, intermediateTrees2);
  } else {
    auto start = std::chrono::high_resolution_clock::now();
    constructTrees(inputTrees, intermediateTrees);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Avg time ftm construction: " << duration.count()*0.000001*(1/(float)inputTrees.size()) << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    constructTrees(inputTrees2, intermediateTrees2);
    int treetype = 0;
    for(size_t i = 0; i < inputTrees.size(); i++) {
      vtkSmartPointer<vtkUnstructuredGrid> mt_nodes
        = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
      int nmax = 0;
      int nmin = 0;
      for(int j = 0; j < mt_nodes->GetNumberOfPoints(); j++) {
        if(mt_nodes->GetPointData()
             ->GetArray("CriticalType")
             ->GetComponent(j, 0)
           == 3) {
          nmax++;
          if(nmax > 1) {
            treetype = 1;
            break;
          }
        }
        if(mt_nodes->GetPointData()
             ->GetArray("CriticalType")
             ->GetComponent(j, 0)
           == 0) {
          nmin++;
          if(nmin > 1) {
            treetype = 0;
            break;
          }
        }
      }
    }
    for(size_t i = 0; i < inputTrees.size(); i++) {
      vtkSmartPointer<vtkUnstructuredGrid> mt_nodes
        = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
      vtkSmartPointer<vtkUnstructuredGrid> mt_arcs
        = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
      adjacenyLists.push_back(ftmToAdjList<dataType>(mt_nodes, mt_arcs, treetype, true));
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Avg time adj list construction: " << duration.count()*0.000001*(1/(float)inputTrees.size()) << " seconds" << std::endl;
  }

  // --- Call base
  std::vector<std::vector<double>> treesDistMat(
    numInputs, std::vector<double>(numInputs));
  if(baseModule == 0)
    execute<dataType>(intermediateTrees, intermediateTrees2, treesDistMat);
  else
    execute<dataType>(adjacenyLists, intermediateTrees, treesDistMat);

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

template<class dataType>
std::tuple<std::vector<dataType>, std::vector<std::vector<int>>, int>
  ttkMergeTreeDistanceMatrix::ftmToAdjList(
    vtkSmartPointer<vtkUnstructuredGrid> &mt_nodes,
    vtkSmartPointer<vtkUnstructuredGrid> &mt,
    int treeType,
    bool scalarLabels) {
  std::map<std::tuple<float, float, float>, dataType> posToScalarMap;
  std::vector<dataType> nodeScalars(mt->GetNumberOfPoints());
  std::vector<int> nodeTypes(mt->GetNumberOfPoints());
  std::vector<int> upIDs(mt->GetNumberOfCells());
  std::vector<int> downIDs(mt->GetNumberOfCells());
  std::vector<float> regionSizes(mt->GetNumberOfCells());
  std::vector<dataType> persistences(mt->GetNumberOfCells());
  std::vector<dataType> nodePersistences(mt->GetNumberOfPoints(), static_cast<dataType>(1000.));
  for(int i = 0; i < mt->GetNumberOfPoints(); i++) {
    int nid = mt_nodes->GetPointData()->GetArray("NodeId")->GetComponent(i, 0);
    int nt
      = mt_nodes->GetPointData()->GetArray("CriticalType")->GetComponent(i, 0);
    dataType ns
      = mt_nodes->GetPointData()->GetArray("Scalar")->GetComponent(nid, 0);
    nodeTypes[nid] = nt;
    nodeScalars[nid] = ns;
  }
  for(int i = 0; i < mt->GetNumberOfCells(); i++) {
    int downId
      = int(mt->GetCellData()->GetArray("downNodeId")->GetComponent(i, 0));
    int upId = int(mt->GetCellData()->GetArray("upNodeId")->GetComponent(i, 0));
    float rs = mt->GetCellData()->GetArray("RegionSize")->GetComponent(i, 0);
    dataType pers = abs(nodeScalars[upId] - nodeScalars[downId]);
    upIDs[i] = upId;
    downIDs[i] = downId;
    regionSizes[i] = rs;
    persistences[i] = pers;
  }

  int rootID = -1;
  std::vector<std::vector<int>> children;
  for(size_t i = 0; i < nodeScalars.size(); i++) {
    if(treeType == 0 and nodeTypes[i] == 3) {
      rootID = i;
    }
    if(treeType == 1 and nodeTypes[i] == 0) {
      rootID = i;
    }
    children.emplace_back(std::vector<int>());
  }
  for(size_t i = 0; i < upIDs.size(); i++) {
    int downId = downIDs[i];
    int upId = upIDs[i];
    children[upId].push_back(downId);
    nodePersistences[downId]
      = persistences[i]; // std::abs(nodeScalars[downId]-nodeScalars[upId]);
  }

  // random edge contractions for testing
  // std::vector<int> contractions(children.size());
  // for(int i=0; i<children.size(); i++) contractions[i] = i;
  // for (int i=0; i<mt->GetNumberOfCells(); i++){
  //   int downId =
  //   int(mt->GetCellData()->GetArray("downNodeId")->GetComponent(i,0)); int
  //   upId =
  //   contractions[int(mt->GetCellData()->GetArray("upNodeId")->GetComponent(i,0))];
  //   float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  //   if(r<0.6 && upId!=rootID && children[downId].size()>0){
  //     int rIdx;
  //     for(int cIdx=0; cIdx<children[upId].size(); cIdx++)
  //     if(children[upId][cIdx]==downId) rIdx = cIdx; children[upId][rIdx] =
  //     children[upId].back(); children[upId].pop_back(); for(auto cc :
  //     children[downId]) children[upId].push_back(cc); contractions[downId] =
  //     upId;
  //   }
  // }

  // transform to dfs order
  std::vector<int> dfs_reordering(children.size(), -1);
  // for(int i=0; i<children.size(); i++) dfs_reordering[i] = i;
  std::stack<int> stack;
  stack.push(rootID);
  int currIdx = children.size() - 1;
  while(!stack.empty()) {
    int nIdx = stack.top();
    stack.pop();
    dfs_reordering[nIdx] = currIdx;
    currIdx--;
    for(int cIdx : children[nIdx]) {
      stack.push(cIdx);
    }
  }
  for(size_t i = 0; i < children.size(); i++) {
    if(dfs_reordering[i] == -1) {
      dfs_reordering[i] = currIdx;
      currIdx--;
    }
  }
  std::vector<std::vector<int>> children_dfs(children.size());
  std::vector<dataType> labels_dfs(nodePersistences.size());
  std::vector<dataType> scalars_dfs(nodeScalars.size());
  for(size_t i = 0; i < children.size(); i++) {
    int dfsIdx = dfs_reordering[i];
    labels_dfs[dfsIdx] = nodePersistences[i];
    scalars_dfs[dfsIdx] = nodeScalars[i];
    for(int c : children[i]) {
      children_dfs[dfsIdx].push_back(dfs_reordering[c]);
    }
  }

  return std::make_tuple(scalarLabels ? scalars_dfs : labels_dfs, children_dfs,
                         labels_dfs.size() - 1);
}