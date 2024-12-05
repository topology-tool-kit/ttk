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
  this->SetNumberOfOutputPorts(1);
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
  return 0;
}

bool ttkMorseSmallComplexStability::updateVisitedVertices(const int &globalId, 
                                                          std::vector<int> &localToGlobal,
                                                          int &localId){
  
  for (int i = 0 ; i < localToGlobal.size();i++ ){
    if (localToGlobal[i]==globalId){
      localId=i;
      return true;
    }
  }
  if(localId==-1){
    localId=localToGlobal.size();
    localToGlobal.push_back(globalId);
    std::array<double, 3> newCoords;
    points->GetPoint(pointId, sourceCoords.data());
    coords.push_back(newCoords);
  }
  return false;
}

void ttkMorseSmallComplexStability::updateAdjacencyMatrix(const int &sourceLocalId,
                                                            const int &destinationLocalId,
                                                            const int &separatrixLocalId,
                                                            std::vector<std::vector<int>> &adjacencyMatrix){
  int n_row = adjacencyMatrix.size();
  int n_col = n_row == 0 ? 0 : adjacencyMatrix[0].size();
  if(sourceLocalId < n_row){
    if(destinationLocalId < n_col){
      if(adjacencyMatrix[sourceLocalId][destinationLocalId]!=-1)std::cout<<"multiple edge found"<<std::endl;
      adjacencyMatrix[sourceLocalId][destinationLocalId] = separatrixLocalId;

    }else{
      for (int i = 0 ; i < n_row ; i++){
        adjacencyMatrix[i].push_back(-1);
      }
      adjacencyMatrix[sourceLocalId][destinationLocalId]=separatrixLocalId;
    }
  }else{
    if(destinationLocalId < n_col){
      std::vector<int> newRow(n_col, -1);
      adjacencyMatrix.push_back(newRow);
      adjacencyMatrix[sourceLocalId][destinationLocalId]=separatrixLocalId;
    }else{
      for (int i = 0 ; i < n_row ; i++){
        adjacencyMatrix[i].push_back(-1);
      }
      std::vector<int> newRow(n_col+1, -1);
      adjacencyMatrix.push_back(newRow);
      adjacencyMatrix[sourceLocalId][destinationLocalId]=separatrixLocalId;
    }

  }
}     

void ttkMorseSmallComplexStability::computePointIds(const int &cellId_1,
                      const int &cellId_2,
                      const int &sourceGlobalId,
                      const int &destinationGlobalId,
                      vtkDataSet *block,
                      vtkIdType &sourcePointId,
                      vtkIdType &destinationPointId){

  vtkIdList* pointIds = block->GetCell(cellId_1)->GetPointIds();
  int id_0 = pointIds->GetId(0);
  int id_1 = pointIds->GetId(1);
  pointIds = block->GetCell(cellId_2)->GetPointIds();
  int id_2 = pointIds->GetId(0);
  int id_3 = pointIds->GetId(1);
  ttkSimplexIdTypeArray* cellId = ttkSimplexIdTypeArray::SafeDownCast(block->GetPointData()->GetArray(ttk::MorseSmaleCellIdName));

  if((int)cellId->GetValue(id_0) == sourceGlobalId)sourcePointId = id_0;
  else if((int)cellId->GetValue(id_1) == sourceGlobalId)sourcePointId = id_1;
  else if((int)cellId->GetValue(id_2) == sourceGlobalId)sourcePointId = id_2;
  else if((int)cellId->GetValue(id_3) == sourceGlobalId)sourcePointId = id_3;

  if((int)cellId->GetValue(id_0) == destinationGlobalId)destinationPointId = id_0;
  else if((int)cellId->GetValue(id_1) == destinationGlobalId)destinationPointId = id_1;
  else if((int)cellId->GetValue(id_2) == destinationGlobalId)destinationPointId = id_2;
  else if((int)cellId->GetValue(id_3) == destinationGlobalId)destinationPointId = id_3;

}


void ttkMorseSmallComplexStability::computeGraphMinor(const std::vector<std::vector<int>> &adjacencyMatrixFull, 
                                                        GraphMatrix &adjacencyMatrix){
  int n_row = adjacencyMatrixFull.size();
  int n_col = adjacencyMatrixFull[0].size();
  adjacencyMatrix.resize(n_row);
  for (int i = 0 ; i < n_row ; i++){
    adjacencyMatrix[i].resize(n_row);
  }

  for (int i = 0 ; i < n_row ; i++){
    for (int j = 0 ; j < n_col ; j++){
      if(adjacencyMatrixFull[i][j]!=-1){
        for (int k = i+1 ; k < n_row ; k++){
          if(adjacencyMatrixFull[k][j]!=-1){
            std::pair<int, int> newEdge = std::make_pair(adjacencyMatrixFull[k][j], adjacencyMatrixFull[i][j]);
            //if(adjacencyMatrix[i][k].has_value())std::cout<<"multiple edge in reduced matrix for indices "<<k<<", "<<j<<std::endl;
            adjacencyMatrix[i][k].push_back(newEdge);
          }
        }
      }
    }
  }                                                        
}

void ttkMorseSmallComplexStability::appendPoint(vtkPoints* points, 
                                                const int &index, 
                                                std::vector<std::array<double, 3>> &coords){
  std::array<double, 3> newCoords;
  points->GetPoint(index, newCoords.data());
  coords.push_back(newCoords);
}

int ttkMorseSmallComplexStability::prepareData(vtkDataSet* block, 
                                                std::vector<int> &localToGlobal, 
                                                GraphMatrix &adjacencyMatrix,
                                                std::vector<std::array<double, 3>> &coords,
                                                int &n_separatrices){
  


  vtkPoints* points = block->GetPoints();
  vtkPointData* pointData = block->GetPointData();

  //vtkSignedCharArray* cellDimensions = vtkSignedCharArray::SafeDownCast(block->GetPointData()->GetArray(ttk::MorseSmaleCellDimensionName));  

  vtkCellData* cellData = block->GetCellData();
  ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleSeparatrixIdName));
  ttkSimplexIdTypeArray *sourceIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleSourceIdName));
  ttkSimplexIdTypeArray *destinationIds = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray(ttk::MorseSmaleDestinationIdName));
  int n_cells = block->GetNumberOfCells();

  std::vector<std::vector<int>> adjacencyMatrixFull;
  std::vector<int> sourceLocalToGlobal;
  std::vector<int> destinationLocalToGlobal;

  int cellId_1=0;
  int separatrixLocalId=-1;

  while (cellId_1 < n_cells ){
    int separatrixId = separatrixIds->GetValue(cellId_1);
    int cellId_2 = cellId_1 + 1;
    int nextSeparatrixId=separatrixIds->GetValue(cellId_2);
    while (nextSeparatrixId==separatrixId && cellId_2 < n_cells){
      nextSeparatrixId=separatrixIds->GetValue(++cellId_2);
    }
    cellId_2--;
    separatrixLocalId++;
    int sourceGlobalId = sourceIds->GetValue(cellId_1);
    int destinationGlobalId = destinationIds->GetValue(cellId_2);
    int destinationLocalId = -1;
    int sourceLocalId = -1;
    vtkIdType sourcePointId;
    vtkIdType destinationPointId;

    computePointIds(cellId_1, 
                    cellId_2, 
                    sourceGlobalId, 
                    destinationGlobalId, 
                    block, 
                    sourcePointId, 
                    destinationPointId);

    bool foundSource = updateVisitedVertices(sourceGlobalId, sourceLocalToGlobal, sourceLocalId);
    bool foundDestination = updateVisitedVertices(destinationGlobalId, destinationLocalToGlobal, destinationLocalId);
//    std::cout<<"source point id = "<<sourcePointId<<std::endl;
//    std::cout<<"destination point id = "<<destinationPointId<<std::endl;
    if(!foundDestination)appendPoint(points, destinationPointId, coords);
    //std::cout<<"sourceLocalId = "<<sourceLocalId<<",   destinationLocalId = "<<destinationLocalId<<std::endl;
    //std::cout<<"separatrixLocalId = "<<separatrixLocalId<<std::endl;
    updateAdjacencyMatrix(destinationLocalId, sourceLocalId, separatrixLocalId, adjacencyMatrixFull);
    cellId_1 = cellId_2 + 1;
  }

  n_separatrices=separatrixLocalId+1;
  localToGlobal = std::move(destinationLocalToGlobal);

for (int i = 0 ; i < adjacencyMatrixFull.size(); i++){
  for (int j = 0 ;j < adjacencyMatrixFull[i].size(); j++){
    std::cout<<adjacencyMatrixFull[i][j]<<"  ";
  }
  std::cout<<std::endl;
}
//
 computeGraphMinor(adjacencyMatrixFull, adjacencyMatrix);
 
std::cout<<std::endl;
for (int i = 0 ; i < adjacencyMatrix.size(); i++){
  for (int j = 0 ;j < adjacencyMatrix[i].size(); j++){
    if(!adjacencyMatrix[i][j].empty()){
      for (int k = 0 ; k < adjacencyMatrix[i][j].size(); k++){
        std::cout<<"("<<adjacencyMatrix[i][j][k].first<<", "<<adjacencyMatrix[i][j][k].second<<") ";
      }
    }
    else std::cout<<"(X,  X) "; 
    //<<std::setw(12)
  }
  std::cout<<std::endl;
  }
  return 1;
}

int ttkMorseSmallComplexStability::execute( vtkMultiBlockDataSet* &multiBlock1_Separatrices, vtkMultiBlockDataSet* &output1_Separatrices){    


    int n_blocks = multiBlock1_Separatrices->GetNumberOfBlocks();
    std::vector<std::vector<int>> LocalToGlobal(n_blocks);
    std::vector<GraphMatrix> adjacencyMatrices(n_blocks);
    std::vector<std::vector<std::array<double, 3>>> coords(n_blocks);
    std::vector<int> separatrixCountForEachBlock(n_blocks);
    for (int i = 0 ; i < n_blocks ; i++){
      vtkDataSet* block = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));
      std::cout<<"block "<<i<<" processing..."<<std::endl;
      ttkMorseSmallComplexStability::prepareData(block, LocalToGlobal[i], adjacencyMatrices[i], coords[i], separatrixCountForEachBlock[i]);
    }

    std::vector<std::vector<int>>edgeOccurenceForEachBlock(n_blocks);

    if(!this->buildOccurenceArrays(adjacencyMatrices, 
                              separatrixCountForEachBlock, 
                              coords,
                              edgeOccurenceForEachBlock)){
                              
      printErr("Could not compute occurence array.");
      return 0;      
    };

    for (int i = 0 ; i < n_blocks ; i++){

      vtkDataSet* block = vtkDataSet::SafeDownCast(output1_Separatrices->GetBlock(i));
      int cellNumber = block->GetNumberOfCells();
      vtkCellData* blockCellData = block->GetCellData();
      ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(blockCellData->GetArray(ttk::MorseSmaleSeparatrixIdName));

      vtkNew<vtkIntArray> occurenceCount;
      occurenceCount->SetNumberOfComponents(1);
      occurenceCount->SetName(ttk::MorseSmaleStabilityOccurenceCount);

      vtkIdType currentCellId = 0;
      int currentSeparatrixId = separatrixIds->GetValue(currentCellId);
      int separatrixCount = 0;
      occurenceCount->InsertNextValue(edgeOccurenceForEachBlock[i][separatrixCount]);

      std::cout<<"size of edgeOccurenceForEachBlock["<<i<<"] = "<<edgeOccurenceForEachBlock[i].size()<<std::endl;

      for (int j = 1 ; j < cellNumber ; j++){
        if(separatrixIds->GetValue(j)!=currentSeparatrixId){
          currentSeparatrixId = separatrixIds->GetValue(j);
          separatrixCount++;
        }
        occurenceCount->InsertNextValue(edgeOccurenceForEachBlock[i][separatrixCount]);
      }
      std::cout<<"separatrixCount for block "<<i<<" = "<<separatrixCount<<std::endl;
      blockCellData->AddArray(occurenceCount);
    }
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
  
  vtkMultiBlockDataSet *input1_Separatrices = vtkMultiBlockDataSet::GetData(inputVector[0]);

  if(input1_Separatrices == nullptr){
    this->printErr("No edges to perform calculation.");
    return -1;
  }

  int status = 0; 
  vtkMultiBlockDataSet *output1_Separatrices = vtkMultiBlockDataSet::GetData(outputVector, 0);
  output1_Separatrices->ShallowCopy(input1_Separatrices);

  status = this->execute(input1_Separatrices, output1_Separatrices);

  std::cout<<"balise finale"<<std::endl;

  if(status != 1)
    return 0;

  //for (size_t i = 0 ; i < output1_Separatrices->GetNumberOfBlocks(); i++){
  //  ((vtkDataSet*)(output1_Separatrices->GetBlock(i)))->GetPointData()->AddArray(edgesOccurences[i]);
  //}
  return 1;
}
