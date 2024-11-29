#include <ttkMorseSmallComplexStability.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vector>
#include <utility>
#include <optional>
#include <queue>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMorseSmallComplexStability);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
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
ttkMorseSmallComplexStability::ttkMorseSmallComplexStability() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}


int ttkMorseSmallComplexStability::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmallComplexStability::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
	if(port == 1){
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
	}	
  return 0;
}

void ttkMorseSmallComplexStability::updateVertexData(const int &globalId, 
                                                      const vtkIdType &pointId, 
                                                      vtkPoints* points,
                                                      std::vector<std::array<double, 3>> &coords, 
                                                      std::vector<int> &localToGlobal,
                                                      int &localId){
  
  for (unsigned int i = 0 ; i < localToGlobal.size();i++ ){
    if (localToGlobal[i]==globalId){
      localId=i;
      return;
    }
  }
  if(localId==-1){
    localId=localToGlobal.size();
    localToGlobal.push_back(globalId);
    std::array<double, 3> newCoords;
    points->GetPoint(pointId, newCoords.data());
    coords.push_back(newCoords);
  }
  return;
}

void ttkMorseSmallComplexStability::updateVertexLinkList(std::vector<std::vector<int>> &vertexLinks, 
                                                          const int &v1, 
                                                          const int &v2){
  if ((unsigned)v1 >= vertexLinks.size()){
    std::vector<int> newEdgeList;
    newEdgeList.push_back(v2);                                                          
    vertexLinks.push_back(newEdgeList);
  }else{
    vertexLinks[v1].push_back(v2);
  }
}     

void ttkMorseSmallComplexStability::computePointIds(vtkCell* cell_1,
                      vtkCell* cell_2,
                      const int &sourceGlobalId,
                      const int &destinationGlobalId,
                      vtkDataSet *block,
                      vtkIdType &srcPointId,
                      vtkIdType &destPointId){
  
  vtkSmartPointer<vtkIdList> sourcePointIds = cell_1-> GetPointIds();
  vtkSmartPointer<vtkIdList> destinationPointIds = cell_2-> GetPointIds();
  ttkSimplexIdTypeArray* cellId = ttkSimplexIdTypeArray::SafeDownCast(block->GetPointData()->GetArray(ttk::MorseSmaleCellIdName));

  srcPointId = ((int)cellId->GetValue((sourcePointIds->GetId(0)))==sourceGlobalId) ? sourcePointIds->GetId(0) : sourcePointIds->GetId(1);
  destPointId = ((int)cellId->GetValue((destinationPointIds->GetId(0)))==destinationGlobalId) ? destinationPointIds->GetId(0) : destinationPointIds->GetId(1);
  
}

void ttkMorseSmallComplexStability::appendPoint(vtkPoints* points, 
                                                const int &index, 
                                                std::vector<std::array<double, 3>> &coords){
  std::array<double, 3> newCoords;
  points->GetPoint(index, newCoords.data());
  coords.append(newCoords);
}

int ttkMorseSmallComplexStability::prepareData(vtkDataSet* block, 
                                                std::vector<int> &localToGlobal, 
                                                std::vector<std::pair<int, int>> &edges,
                                                std::vector<std::array<double, 3>> &coords,
                                                std::vector<float>&sfValues){
  

  vtkPoints* points = block->GetPoints();
  vtkPointData pointData = points->GetPointData();
  vtkSignedCharArray* ttkScalarMask = vtkSignedCharArray::SafeDownCast(block->GetPointData()->GetArray(ttk::MaskScalarFieldName));
  vtkSignedCharArray* cellDimensions = vtkSignedCharArray::SafeDownCast(block->GetPointData()->GetArray(ttk::MorseSmaleCellDimensionName));
  
  std::vector<std::vector<std::array<double, 3>>> coords;
  std::vector<std::vector<int>> adjacencyMatrix;
  int n_points = block->GetNumberOfPoint();
  int currentVertex1 = 0;
  int currentVertex2;
  for (int i = 1 ; i < n_points; i++){
    if (ttkScalarMask->GetValue(i) == 1){
      continue;
    }
    currentVertex2 = i;


    
  }

  vtkCellData* cellData = block->GetCellData();
  ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleSeparatrixIdName));
  vtkSignedCharArray *separatrixTypes = vtkSignedCharArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleSeparatrixTypeName));
  ttkSimplexIdTypeArray *sourceIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleSourceIdName));
  ttkSimplexIdTypeArray *destinationIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleDestinationIdName));

  std::vector<std::vector<int>> criticalVertexType0;
  std::vector<std::vector<int>> criticalVertexType1;
  std::vector<std::vector<int>> criticalVertexType2;
  std::vector<std::vector<int>> criticalVertexType3;
  
  std::vector<int> localToGlobalType0;
  std::vector<int> localToGlobalType1;
  std::vector<int> localToGlobalType2;
  std::vector<int> localToGlobalType3;

  std::vector<std::array<double, 3>> coordsType0;
  std::vector<std::array<double, 3>> coordsType1;
  std::vector<std::array<double, 3>> coordsType2;
  std::vector<std::array<double, 3>> coordsType3;

  std::vector<float> sfValuesType0;
  std::vector<float> sfValuesType1;
  std::vector<float> sfValuesType2;
  std::vector<float> sfValuesType3;

  int n_cells = block->GetNumberOfCells();
  int currentCellId_1=0;

  while (currentCellId_1 < n_cells ){
    vtkCell* currentCell_1 = block->GetCell(currentCellId_1);
    int currentSeparatrixId = separatrixIds->GetValue(currentCellId_1);
    int currentSeparatrixType = static_cast<int>(separatrixTypes->GetValue(currentCellId_1));
    int currentCellId_2 = currentCellId_1 + 1;
    int potentialNextSeparatrixId=separatrixIds->GetValue(currentCellId_2);
    while (potentialNextSeparatrixId==currentSeparatrixId && currentCellId_2 < n_cells){
      potentialNextSeparatrixId=separatrixIds->GetValue(++currentCellId_2);
    }
    vtkCell* currentCell_2 = block->GetCell(--currentCellId_2);
    int destinationGlobalId = destinationIds->GetValue(currentCellId_2);
    int sourceGlobalId = sourceIds->GetValue(currentCellId_1);
    int destinationLocalId = -1;
    int sourceLocalId = -1;
    vtkIdType srcPointId;
    vtkIdType destPointId;

    computePointIds(currentCell_1, 
                      currentCell_2, 
                      sourceGlobalId, 
                      destinationGlobalId, 
                      block, 
                      srcPointId, 
                      destPointId);

    switch (currentSeparatrixType){
      case 0:
        updateVertexData(sourceGlobalId, srcPointId, points, coords, localToGlobalType0, sourceLocalId);
        updateVertexData(destinationGlobalId, destPointId, points, coords, localToGlobalType1, destinationLocalId);
        updateVertexLinkList(criticalVertexType0, sourceLocalId, destinationLocalId);
        updateVertexLinkList(criticalVertexType1, destinationLocalId, sourceLocalId);
        break;
      case 1:
        updateVertexData(sourceGlobalId, srcPointId, points, coords, localToGlobalType1, sourceLocalId);
        updateVertexData(destinationGlobalId, destPointId, points, coords, localToGlobalType2, destinationLocalId);
        updateVertexLinkList(criticalVertexType1, sourceLocalId, destinationLocalId);
        updateVertexLinkList(criticalVertexType2, destinationLocalId, sourceLocalId);
        break;
 
      case 2:
        updateVertexData(sourceGlobalId, srcPointId, points, coords, localToGlobalType1, sourceLocalId);
        updateVertexData(destinationGlobalId, destPointId, points, coords, localToGlobalType2, destinationLocalId);
        updateVertexLinkList(criticalVertexType2, sourceLocalId, destinationLocalId);
        updateVertexLinkList(criticalVertexType3, destinationLocalId, sourceLocalId);
        break;
      }
    currentCellId_1 = currentCellId_2 + 1;

  }
}

int ttkMorseSmallComplexStability::execute( vtkMultiBlockDataSet* &multiBlock1_Separatrices,
                                            vtkUnstructuredGrid* &minimalGraph,
                                            std::vector<vtkSmartPointer<vtkIntArray>> &edgesOccurences){
    

    int n_blocks = multiBlock1_Separatrices->GetNumberOfBlocks();
    std::vector<std::vector<int>> LocalToGlobal(n_blocks);
    std::vector<std::vector<std::pair<int, int>>> Edges(n_blocks);
    std::vector<std::vector<std::array<double, 3>>> Coords(n_blocks);
    std::vector<std::vector<float>> SfValues(n_blocks);

    //Build graph from separatrices
    //--> output : vector<vector<int>> vertexLocalToGlobal, 
    //              vector<vector<std::pair<int , int>>> edges, (id de chaque edge = separatrixId)

    //Build coords vector of vertex and scalar values of vertex
    // --> output : std::vector<std::vector<std::array<double, 3>>> coordinates
    //              std::vector<std::vector<float>> sfValues;
    
    for (int i = 0 ; i < n_blocks ; i++){
      vtkDataSet* block = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));
      ttkMorseSmallComplexStability::prepareData(block, LocalToGlobal[i], Edges[i], Coords[i], SfValues[i]);
    }

    ////Build equivalent classes with assignement : 
    ////--> output : vector<vector<int>> classIdToVertexId
//
    //std::vector<std::vector<int>> classIdToVertexId(n_blocks);
    //this->buildVertexEquivalenceClasses(coords, sfValues, classIdToVertexId);
   //
    ////Build N adjacency matrices from the other vectors
    //// --> output : for n < N vector<vector<int>> matrix_n for (int i = 0 ; i < edges_n.size(); i++); 
    ////                                                          matrix_n[classIdToVertexId[edges_n[i]].first][classIdToVertxId[edges_[i].second]]=1
//
    //using GraphMatrix = std::vector<std::vector<int>>; 
    //std::vector<GraphMatrix> adjacencyMatrices(n_blocks);
    //this->buildAdjacencyMatrices(adjacencyMatrices, classToVertexId);
//
    ////Build 1 occurence matrix 
    //// --> somme des adjacency matrices : std::vector<std::vector<int>> ocurrenceMatrix
//
    //GraphMatrix occurenceMatrix;
    //this->buildOccurenceMatrix(adjacencyMatrices, occurenceMatrix);
//
    ////Build for n < N std::vector<int> occurenceCount_n
    //// --> output : occurenceCount_n[i]=occurenceMatrix[classIdToVertexId[edges_n[i]].first][classIdToVertxId[edges_[i].second]]
//
    //std::vector<std::vector<int>> occurenceCountForEachEdge;
    //ttkMorseSmallComplexStability::countOccurencesOfEdges(occurenceMatrix, classIdToVertexId, occurenceCountForEachEdge);
    //for (int i = 0 ; i < n_blocks ; i++){
//
    //  vtkDataSet* block = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));
    //  int cellNumber = block->GetNumberOrCells();
//
    //  vtkCellData* blockCellData = block->GetCellData();
    //  vtkDataArray* separatrixIds = blockCellData->GetAray(ttk::MorseSmaleSeparatrixIdName);
    //  vtkNew<vtkIntArray> occurenceCount;
    //  occurencesCount->SetNumberOfTuples(cellNumber);
    //  occurencesCount->SetNumberOfComponents(1);
//
    //  for (int j = 0 ; j < cellNumber ; j++){
    //    int currentSeparatrixId = separatrixIds->GetValue(j);
    //    int newOccurenceCount = occurenceCountForEachEdge[i][currentSeparatrixId];
    //    cellNumber->SetTuple1(j, newOccurenceCount);
    //  }
    //}
//
    ////Build partial isomorphic graphs : 
    //// --> for index i_1, ... i_k, build graph from the adjacency matrix matrix_i_1 && ... && matrix_i_k 
    //// --> vertices are barycenters of each equivalent class of vertices
    return 1;
                                            }


int ttkMorseSmallComplexStability::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

    for (int i = 0; i < (*inputVector)->GetNumberOfInformationObjects(); ++i)
  {
      vtkInformation* info = (*inputVector)->GetInformationObject(i);
      info->Print(cout);
  }
  
  auto minimalGraph = vtkUnstructuredGrid::GetData(outputVector, 1);

  vtkMultiBlockDataSet *input1_Separatrices= vtkMultiBlockDataSet::GetData(inputVector[0]);

  int n_blocks = input1_Separatrices->GetNumberOfBlocks();
  std::vector<vtkSmartPointer<vtkIntArray>> edgesOccurences(n_blocks);

  if(input1_Separatrices == nullptr){
    this->printErr("No edges to perform calculation.");
    return -1;
  }

  int status = 0; 
  status = this->execute(input1_Separatrices,
                          minimalGraph,
                          edgesOccurences);

  if(status != 1)
    return 0;

  vtkMultiBlockDataSet *output1_Separatrices = vtkMultiBlockDataSet::GetData(outputVector, 0);

  output1_Separatrices->ShallowCopy(input1_Separatrices);

  for (size_t i = 0 ; i < output1_Separatrices->GetNumberOfBlocks(); i++){
    ((vtkDataSet*)(output1_Separatrices->GetBlock(i)))->GetPointData()->AddArray(edgesOccurences[i]);
  }

  return 1;
}
