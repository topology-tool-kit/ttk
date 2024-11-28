#include <ttkMorseSmallComplexStability.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vector>
#include <utility>

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

int ttkMorseSmallComplexStability::graphFromSeparatrices(vtkDataSet* block, 
                                                          std::vector<int> &localToGlobal, 
                                                          std::vector<std::pair<int, int>> &edges){

}

int ttkMorseSmallComplexStability::execute( vtkMultiBlockDataSet* &multiBlock1_Separatrices,
                                            vtkUnstructuredGrid* &minimalGraph,
                                            std::vector<vtkSmartPointer<vtkIntArray>> &edgesOccurences){
    

    int n_blocks = multiBlock1_Separatrices->GetNumberOfBlocks();
    std::vector<std::vector<int>> localToGlobal(n_blocks);
    std::vector<std::vector<std::pair<int, int>>> edges(n_blocks);
    std::vector<std::vector<std::array<double, 3>>> coords(n_blocks);
    std::vector<std::vector<float>> sfValues(n_blocks);

    //Build graph from separatrices
    //--> output : vector<vector<int>> vertexLocalToGlobal, 
    //              vector<vector<std::pair<int , int>>> edges, (id de chaque edge = separatrixId)

    //Build coords vector of vertex and scalar values of vertex
    // --> output : std::vector<std::vector<std::array<double, 3>>> coordinates
    //              std::vector<std::vector<float>> sfValues;
    
    for (int i = 0 ; i < n_blocks ; i++){
      vtkDataSet* block = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));
      ttkMorseSmallComplexStability::prepareData(block, localToGlobal[i], edges[i], coords[i], sfValues[i]);
    }

    //Build equivalent classes with assignement : 
    //--> output : vector<vector<int>> classIdToVertexId

    std::vector<std::vector<int>> classIdToVertexId(n_blocks);
    this->buildVertexEquivalenceClasses(coords, sfValues, classIdToVertexId);
   
    //Build N adjacency matrices from the other vectors
    // --> output : for n < N vector<vector<int>> matrix_n for (int i = 0 ; i < edges_n.size(); i++); 
    //                                                          matrix_n[classIdToVertexId[edges_n[i]].first][classIdToVertxId[edges_[i].second]]=1

    using GraphMatrix = std::vector<std::vector<int>>; 
    std::vector<GraphMatrix> adjacencyMatrices(n_blocks);
    this->buildAdjacencyMatrices(adjacencyMatrices, classToVertexId);

    //Build 1 occurence matrix 
    // --> somme des adjacency matrices : std::vector<std::vector<int>> ocurrenceMatrix

    GraphMatrix occurenceMatrix;
    this->buildOccurenceMatrix(adjacencyMatrices, occurenceMatrix);

    //Build for n < N std::vector<int> occurenceCount_n
    // --> output : occurenceCount_n[i]=occurenceMatrix[classIdToVertexId[edges_n[i]].first][classIdToVertxId[edges_[i].second]]

    std::vector<std::vector<int>> occurenceCountForEachEdge;
    ttkMorseSmallComplexStability::countOccurencesOfEdges(occurenceMatrix, classIdToVertexId, occurenceCountForEachEdge);
    for (int i = 0 ; i < n_blocks ; i++){

      vtkDataSet* block = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));
      int cellNumber = block->GetNumberOrCells();

      vtkCellData* blockCellData = block->GetCellData();
      vtkDataArray* separatrixIds = blockCellData->GetAray(ttk::MorseSmaleSeparatrixIdName);
      vtkNew<vtkIntArray> occurenceCount;
      occurencesCount->SetNumberOfTuples(cellNumber);
      occurencesCount->SetNumberOfComponents(1);

      for (int j = 0 ; j < cellNumber ; j++){
        int currentSeparatrixId = separatrixIds->GetValue(j);
        int newOccurenceCount = occurenceCountForEachEdge[i][currentSeparatrixId];
        cellNumber->SetTuple1(j, newOccurenceCount);
      }
    }

    //Build partial isomorphic graphs : 
    // --> for index i_1, ... i_k, build graph from the adjacency matrix matrix_i_1 && ... && matrix_i_k 
    // --> vertices are barycenters of each equivalent class of vertices
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
