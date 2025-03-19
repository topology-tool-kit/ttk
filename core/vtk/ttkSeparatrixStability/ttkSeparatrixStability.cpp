#include <ttkSeparatrixStability.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <Timer.h>
#include <string>
#include <iomanip>
#include <random>

vtkStandardNewMacro(ttkSeparatrixStability);

ttkSeparatrixStability::ttkSeparatrixStability() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkSeparatrixStability::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkSeparatrixStability::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

bool ttkSeparatrixStability::updateVisitedVertices(
  const int &globalId, std::vector<int> &localToGlobal, int &localId) {

  for(unsigned int i = 0; i < localToGlobal.size(); i++) {
    if(localToGlobal[i] == globalId) {
      localId = i;
      return true;
    }
  }
  if(localId == -1) {
    localId = localToGlobal.size();
    localToGlobal.push_back(globalId);
    std::array<double, 3> newCoords;
    points->GetPoint(pointId, sourceCoords.data());
    coords.push_back(newCoords);
  }
  return false;
}

void ttkSeparatrixStability::updateAdjacencyMatrix(
  const int &sourceLocalId,
  const int &destinationLocalId,
  const int &separatrixLocalId,
  GraphMatrixFull &adjacencyMatrix) {
  int n_row = adjacencyMatrix.size();
  int n_col = n_row == 0 ? 0 : adjacencyMatrix[0].size();
  if(sourceLocalId < n_row) {
    if(destinationLocalId < n_col) {
      adjacencyMatrix[sourceLocalId][destinationLocalId] = separatrixLocalId;
    } else {
      for(int i = 0; i < n_row; i++) {
        adjacencyMatrix[i].push_back(-1);
      }
      adjacencyMatrix[sourceLocalId][destinationLocalId] = separatrixLocalId;
    }
  } else {
    if(destinationLocalId < n_col) {
      std::vector<int> newRow(n_col, -1);
      adjacencyMatrix.push_back(newRow);
      adjacencyMatrix[sourceLocalId][destinationLocalId] = separatrixLocalId;
    } else {
      for(int i = 0; i < n_row; i++) {
        adjacencyMatrix[i].push_back(-1);
      }
      std::vector<int> newRow(n_col + 1, -1);
      adjacencyMatrix.push_back(newRow);
      adjacencyMatrix[sourceLocalId][destinationLocalId] = separatrixLocalId;
    }
  }
}

void ttkSeparatrixStability::computePointIds(
  const int &cellId_1,
  const int &cellId_2,
  const int &sourceGlobalId,
  const int &destinationGlobalId,
  vtkDataSet *block,
  vtkIdType &sourcePointId,
  vtkIdType &destinationPointId) {

  vtkIdList *pointIds = block->GetCell(cellId_1)->GetPointIds();
  int id_0 = pointIds->GetId(0);
  int id_1 = pointIds->GetId(1);
  pointIds = block->GetCell(cellId_2)->GetPointIds();
  int id_2 = pointIds->GetId(0);
  int id_3 = pointIds->GetId(1);
  ttkSimplexIdTypeArray *cellId = ttkSimplexIdTypeArray::SafeDownCast(
    block->GetPointData()->GetArray(ttk::MorseSmaleCellIdName));

  if((int)cellId->GetValue(id_0) == sourceGlobalId)
    sourcePointId = id_0;
  else if((int)cellId->GetValue(id_1) == sourceGlobalId)
    sourcePointId = id_1;
  else if((int)cellId->GetValue(id_2) == sourceGlobalId)
    sourcePointId = id_2;
  else if((int)cellId->GetValue(id_3) == sourceGlobalId)
    sourcePointId = id_3;

  if((int)cellId->GetValue(id_0) == destinationGlobalId)
    destinationPointId = id_0;
  else if((int)cellId->GetValue(id_1) == destinationGlobalId)
    destinationPointId = id_1;
  else if((int)cellId->GetValue(id_2) == destinationGlobalId)
    destinationPointId = id_2;
  else if((int)cellId->GetValue(id_3) == destinationGlobalId)
    destinationPointId = id_3;
}

void ttkSeparatrixStability::appendPoint(
  vtkPoints *points,
  const int &index,
  std::vector<std::array<double, 3>> &coords) {
  std::array<double, 3> newCoords;
  points->GetPoint(index, newCoords.data());
  coords.push_back(newCoords);
}

int ttkSeparatrixStability::prepareData(
  vtkDataSet *block,
  std::vector<int> &localToGlobal,
  GraphMatrixFull &adjacencyMatrixFull,
  std::vector<std::array<double, 3>> &coordsSource,
  std::vector<std::array<double, 3>> &coordsDestination,
  int &n_separatrices,
  std::vector<int> &globalSourcePointId,
  std::vector<int> &globalDestinationPointId) {

  vtkPoints *points = block->GetPoints();

  vtkCellData *cellData = block->GetCellData();
  ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(
    cellData->GetArray(ttk::MorseSmaleSeparatrixIdName));
  ttkSimplexIdTypeArray *sourceIds = ttkSimplexIdTypeArray::SafeDownCast(
    cellData->GetArray(ttk::MorseSmaleSourceIdName));
  ttkSimplexIdTypeArray *destinationIds = ttkSimplexIdTypeArray::SafeDownCast(
    cellData->GetArray(ttk::MorseSmaleDestinationIdName));
  int n_cells = block->GetNumberOfCells();

  std::vector<int> sourceLocalToGlobal;
  std::vector<int> destinationLocalToGlobal;

  int cellId_1 = 0;
  int separatrixLocalId = -1;

  while(cellId_1 < n_cells) {
    int separatrixId = separatrixIds->GetValue(cellId_1);
    int cellId_2 = cellId_1 + 1;
    int nextSeparatrixId = separatrixIds->GetValue(cellId_2);
    while(nextSeparatrixId == separatrixId && cellId_2 < n_cells) {
      nextSeparatrixId = separatrixIds->GetValue(++cellId_2);
    }
    cellId_2--;
    separatrixLocalId++;
    int sourceGlobalId = sourceIds->GetValue(cellId_1);
    int destinationGlobalId = destinationIds->GetValue(cellId_2);
    int destinationLocalId = -1;
    int sourceLocalId = -1;
    vtkIdType sourcePointId;
    vtkIdType destinationPointId;

    computePointIds(cellId_1, cellId_2, sourceGlobalId, destinationGlobalId,
                    block, sourcePointId, destinationPointId);

    bool foundSource = updateVisitedVertices(
      sourceGlobalId, sourceLocalToGlobal, sourceLocalId);
    bool foundDestination = updateVisitedVertices(
      destinationGlobalId, destinationLocalToGlobal, destinationLocalId);
    if(!foundDestination){
      appendPoint(points, destinationPointId, coordsDestination);
      globalDestinationPointId.push_back(destinationPointId);
    }
    if(!foundSource && !MergeEdgesOnSaddles){
      appendPoint(points, sourcePointId, coordsSource);
      globalSourcePointId.push_back(sourcePointId);
    }
    updateAdjacencyMatrix(sourceLocalId, destinationLocalId, separatrixLocalId,
                          adjacencyMatrixFull);
    cellId_1 = cellId_2 + 1;
  }

  n_separatrices = separatrixLocalId + 1;
  localToGlobal = std::move(destinationLocalToGlobal);

  return 1;
}

int ttkSeparatrixStability::execute(
  vtkMultiBlockDataSet *&multiBlock1_Separatrices,
  vtkMultiBlockDataSet *&output1_Separatrices) {

  int status;
  int n_blocks = multiBlock1_Separatrices->GetNumberOfBlocks();
  std::vector<std::vector<int>> LocalToGlobal(n_blocks);
  std::vector<GraphMatrixFull> adjacencyMatricesFull(n_blocks);
  std::vector<std::vector<std::array<double, 3>>> coordsDestination(n_blocks);
  std::vector<std::vector<std::array<double, 3>>> coordsSource(n_blocks);
  std::vector<int> separatrixCountForEachBlock(n_blocks);
  std::vector<std::vector<int>> edgeOccurenceForEachBlock(n_blocks);
  std::vector<std::vector<bool>> isomorphismsForEachBlock(n_blocks);
  std::vector<std::vector<std::vector<int>>> matchingArrayForEachBlockSource(n_blocks);
  std::vector<std::vector<std::vector<int>>> matchingArrayForEachBlockDestination(n_blocks);
  std::vector<std::vector<std::vector<int>>> matchingArraySeparatrixForEachBlock(n_blocks);
  std::vector<std::vector<int>> globalSourcePointIdForEachBlock(n_blocks);
  std::vector<std::vector<int>> globalDestinationPointIdForEachBlock(n_blocks);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; i++) {
    vtkDataSet *block
      = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));

    ttkSeparatrixStability::prepareData(
      block, LocalToGlobal[i], adjacencyMatricesFull[i],
      coordsSource[i], coordsDestination[i],
      separatrixCountForEachBlock[i],
      globalSourcePointIdForEachBlock[i],
      globalDestinationPointIdForEachBlock[i]);
  }

  this->setEpsilon(CostDeathBirth);
  this->setWeights(PX, PY, PZ, PF);

  for (int i = 0 ; i < n_blocks; i++){
    std::string destinationSizeString = std::to_string(globalDestinationPointIdForEachBlock[i].size());
    std::string sourceSizeString = std::to_string(globalSourcePointIdForEachBlock[i].size());
    std::string blockIdString = std::to_string(i);
    this->printMsg("Number of critical points (min-max) for block " + blockIdString + " : " + destinationSizeString);
    if(!MergeEdgesOnSaddles)
      this->printMsg("Number of critical points (1sad-2sad) for block " + blockIdString + " : " + sourceSizeString);
  }

  
  status = this->buildOccurenceArrays(adjacencyMatricesFull,  
                                      separatrixCountForEachBlock, 
                                      coordsSource,
                                      coordsDestination,  
                                      MergeEdgesOnSaddles,
                                      edgeOccurenceForEachBlock,
                                      isomorphismsForEachBlock,
                                      matchingArrayForEachBlockSource,
                                      matchingArrayForEachBlockDestination,
                                      matchingArraySeparatrixForEachBlock);



  if(status == 0)return status;

  std::vector<int> classId(n_blocks,-1);
  int idCount = 0;

  for (int i = 0 ; i < n_blocks ; i++){
    if (classId[i]!=-1)continue;
    classId[i]=idCount;
    for(int j =0; j < n_blocks; j++){
      if(isomorphismsForEachBlock[i][j]){
        classId[j]=idCount;
      }
    }
    idCount++;
  }

  for (int i = 0 ; i < n_blocks; i++){
    for (int j = 0 ; j < n_blocks ; j++){
      assert(globalDestinationPointIdForEachBlock[i].size() == matchingArrayForEachBlockDestination[j][i].size());
      assert(globalSourcePointIdForEachBlock[i].size() == matchingArrayForEachBlockSource[j][i].size());
    }
  }


  int n_sourceIds =  std::max_element(globalSourcePointIdForEachBlock.begin(), globalSourcePointIdForEachBlock.end(), 
    [](const std::vector<int>& a, const std::vector<int>& b) {
    return a.size() < b.size();})->size();

  int n_destinationIds =  std::max_element(globalDestinationPointIdForEachBlock.begin(), globalDestinationPointIdForEachBlock.end(), 
    [](const std::vector<int>& a, const std::vector<int>& b) {
    return a.size() < b.size();})->size();

  int n_separatrix = *(std::max_element(separatrixCountForEachBlock.begin(), separatrixCountForEachBlock.end()));

  std::vector<int> shuffleDestinationIds(n_destinationIds);
  std::vector<int> shuffleSourceIds(n_sourceIds);
  std::vector<int> shuffledSeparatrixIds(n_separatrix);

  if(ShuffleCriticalPointsIds){
      for (int i = 0 ; i < n_destinationIds; i++){
        shuffleDestinationIds[i]=i;
      }
      
      std::random_device rd;
      std::mt19937 gen(CriticalPointsShuffleSeed);
      std::shuffle(shuffleDestinationIds.begin(), shuffleDestinationIds.end(), gen);
      
      if(!MergeEdgesOnSaddles){
          for (int i = 0 ; i < n_sourceIds; i++){
            shuffleSourceIds[i]=i+n_destinationIds;
          }
          std::shuffle(shuffleSourceIds.begin(), shuffleSourceIds.end(), gen);
        }
      }

  if(ShuffleSeparatrixIds){
    for (int i = 0 ; i < n_separatrix; i++){
      shuffledSeparatrixIds[i]=i;
    }
    std::random_device rd;
    std::mt19937 gen(SeparatrixShuffleSeed);
    std::shuffle(shuffledSeparatrixIds.begin(), shuffledSeparatrixIds.end(), gen);
  }      
      
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; i++) {

    vtkDataSet *block
      = vtkDataSet::SafeDownCast(output1_Separatrices->GetBlock(i));
    int cellNumber = block->GetNumberOfCells();
    ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(
      block->GetCellData()->GetArray(ttk::MorseSmaleSeparatrixIdName));

    vtkNew<vtkFloatArray> occurenceCount;
    occurenceCount->SetNumberOfComponents(1);
    occurenceCount->SetName(ttk::SeparatrixStabilityOccurenceCount);

    vtkNew<vtkIntArray> isomorphismClassId;
    isomorphismClassId->SetNumberOfComponents(1);
    isomorphismClassId->SetName(ttk::SeparatrixStabilityIsomorphismClassId);

    vtkIdType currentCellId = 0;
    int currentSeparatrixId = separatrixIds->GetValue(currentCellId);
    int separatrixCount = 0;
    float newValue
      = (float)edgeOccurenceForEachBlock[i][separatrixCount] / n_blocks;
   
    occurenceCount->InsertNextValue(newValue);

    for(int j = 1; j < cellNumber; j++) {
      if(separatrixIds->GetValue(j) != currentSeparatrixId) {
        currentSeparatrixId = separatrixIds->GetValue(j);
        separatrixCount++;
        newValue
          = (float)edgeOccurenceForEachBlock[i][separatrixCount] / n_blocks;
      }
      occurenceCount->InsertNextValue(newValue);
    }

    block->GetCellData()->AddArray(occurenceCount);

    isomorphismClassId->InsertNextValue(classId[i]);

    block->GetFieldData()->AddArray(isomorphismClassId);

    for (int j = 0 ; j < n_blocks; j++){
      
      
      vtkNew<vtkIntArray> matchingIdForSeparatrix_j;
      matchingIdForSeparatrix_j->SetNumberOfComponents(1);
      std::string tmp_string_1 = std::string(ttk::SeparatrixStabilityMatchingIdSeparatrixName) 
      + std::to_string(j);
      
      
      const char* indexedArrayNameSeparatrix = tmp_string_1.c_str();
      matchingIdForSeparatrix_j->SetName(indexedArrayNameSeparatrix);
      
      int currentSeparatrixIdBis = separatrixIds->GetValue(0);
      int separatrixCountBis = 0;
      int newId = matchingArraySeparatrixForEachBlock[j][i][separatrixCountBis];
      if(ShuffleSeparatrixIds && newId!=-1)newId=shuffledSeparatrixIds[newId];
      matchingIdForSeparatrix_j->InsertNextValue(newId);
      
      for(int k = 1; k < cellNumber; k++) {
        if(separatrixIds->GetValue(k) != currentSeparatrixIdBis) {
          currentSeparatrixIdBis = separatrixIds->GetValue(k);
          separatrixCountBis++;
          newId = matchingArraySeparatrixForEachBlock[j][i][separatrixCountBis];
          if(ShuffleSeparatrixIds && newId!=-1)newId=shuffledSeparatrixIds[newId];
        }
        matchingIdForSeparatrix_j->InsertNextValue(newId);
      }
      block->GetCellData()->AddArray(matchingIdForSeparatrix_j);
      
      
      vtkNew<vtkIntArray> matchingIdForCriticalPoints_j;
      matchingIdForCriticalPoints_j->SetNumberOfComponents(1);
      std::string tmp_string_2 = std::string(ttk::SeparatrixStabilityMatchingIdName) 
      + std::to_string(j);
      
      const char* indexedArrayNameMatchingsId = tmp_string_2.c_str();
      matchingIdForCriticalPoints_j->SetName(indexedArrayNameMatchingsId);
      vtkPoints* points = block->GetPoints();
      
      for (int k = 0 ; k < points->GetNumberOfPoints(); k++){
        matchingIdForCriticalPoints_j->InsertNextValue(-2);
      }
      
      
      for (unsigned int k = 0 ; k < globalDestinationPointIdForEachBlock[i].size(); k++){
        int globalPointIdThisBlock = globalDestinationPointIdForEachBlock[i][k];
        int matchingIdOtherBlock = matchingArrayForEachBlockDestination[j][i][k];
        if(ShuffleCriticalPointsIds && matchingIdOtherBlock!=-1){
          matchingIdOtherBlock = shuffleDestinationIds[matchingIdOtherBlock];   
        }
        matchingIdForCriticalPoints_j->SetValue(globalPointIdThisBlock, matchingIdOtherBlock);
      }
      if(!MergeEdgesOnSaddles){
        for (unsigned int k = 0 ; k < globalSourcePointIdForEachBlock[i].size(); k++){
          int globalPointIdThisBlock = globalSourcePointIdForEachBlock[i][k];
          int matchingIdOtherBlock=matchingArrayForEachBlockSource[j][i][k];
          if(ShuffleCriticalPointsIds && matchingIdOtherBlock!=n_destinationIds - 1){
            matchingIdOtherBlock = shuffleSourceIds[matchingIdOtherBlock - n_destinationIds];   
          }
          matchingIdForCriticalPoints_j->SetValue(globalPointIdThisBlock, matchingIdOtherBlock);
        }
      }
      block->GetPointData()->AddArray(matchingIdForCriticalPoints_j);
      
  }
}
return status;
}

int ttkSeparatrixStability::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Timer t;
  vtkMultiBlockDataSet *input1_Separatrices
    = vtkMultiBlockDataSet::GetData(inputVector[0]);

  if(input1_Separatrices == nullptr) {
    this->printErr("No edges to perform calculation.");
    return -1;
  }

  int status = 0;
  vtkMultiBlockDataSet *output1_Separatrices
    = vtkMultiBlockDataSet::GetData(outputVector, 0);
  output1_Separatrices->ShallowCopy(input1_Separatrices);

  if(input1_Separatrices->GetNumberOfBlocks() < 2) {
    this->printErr(
      "At least two datasets are required to perform calculations.");
    return -1;
  }

  status = this->execute(input1_Separatrices, output1_Separatrices);

  this->printMsg("Occurence arrays calculated for "
                   + std::to_string(input1_Separatrices->GetNumberOfBlocks())
                   + " blocks",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  if(status != 1)
    return 0;

  return 1;
}
