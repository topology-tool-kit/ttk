#include <ttkMorseSmallComplexStability.h>

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

vtkStandardNewMacro(ttkMorseSmallComplexStability);

ttkMorseSmallComplexStability::ttkMorseSmallComplexStability() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkMorseSmallComplexStability::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmallComplexStability::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

bool ttkMorseSmallComplexStability::updateVisitedVertices(
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

void ttkMorseSmallComplexStability::updateAdjacencyMatrix(
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

void ttkMorseSmallComplexStability::computePointIds(
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

void ttkMorseSmallComplexStability::computeGraphMinor(
  const GraphMatrixFull &adjacencyMatrixFull,
  GraphMatrixMinor &adjacencyMatrix) {
  int n_row = adjacencyMatrixFull.size();
  int n_col = adjacencyMatrixFull[0].size();
  adjacencyMatrix.resize(n_col);
  for(int i = 0; i < n_col; i++) {
    adjacencyMatrix[i].resize(n_col);
  }

  for(int i = 0; i < n_row; i++) {
    for(int j = 0; j < n_col; j++) {
      if(adjacencyMatrixFull[i][j] != -1) {
        for(int k = 0; k < j; k++) {
          if(adjacencyMatrixFull[i][k] != -1) {
            std::pair<int, int> newEdge = std::make_pair(
              adjacencyMatrixFull[i][j], adjacencyMatrixFull[i][k]);
            adjacencyMatrix[k][j].push_back(newEdge);
          }
        }
      }
    }
  }
}

void ttkMorseSmallComplexStability::appendPoint(
  vtkPoints *points,
  const int &index,
  std::vector<std::array<double, 3>> &coords) {
  std::array<double, 3> newCoords;
  points->GetPoint(index, newCoords.data());
  coords.push_back(newCoords);
}

int ttkMorseSmallComplexStability::prepareData(
  vtkDataSet *block,
  std::vector<int> &localToGlobal,
  GraphMatrixFull &adjacencyMatrixFull,
  GraphMatrixMinor &adjacencyMatrixMinor,
  std::vector<std::array<double, 3>> &coordsSource,
  std::vector<std::array<double, 3>> &coordsDestination,
  int &n_separatrices) {

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
    if(!foundDestination)
      appendPoint(points, destinationPointId, coordsDestination);
    if(!foundSource && !MergeEdgesOnSaddles)
      appendPoint(points, sourcePointId, coordsSource);
    updateAdjacencyMatrix(sourceLocalId, destinationLocalId, separatrixLocalId,
                          adjacencyMatrixFull);
    cellId_1 = cellId_2 + 1;
  }

  n_separatrices = separatrixLocalId + 1;
  localToGlobal = std::move(destinationLocalToGlobal);

  if(MergeEdgesOnSaddles) {
    computeGraphMinor(adjacencyMatrixFull, adjacencyMatrixMinor);
  }

  return 1;
}

int ttkMorseSmallComplexStability::execute(
  vtkMultiBlockDataSet *&multiBlock1_Separatrices,
  vtkMultiBlockDataSet *&output1_Separatrices) {

  int status;
  int n_blocks = multiBlock1_Separatrices->GetNumberOfBlocks();
  std::vector<std::vector<int>> LocalToGlobal(n_blocks);
  std::vector<GraphMatrixMinor> adjacencyMatricesMinor(n_blocks);
  std::vector<GraphMatrixFull> adjacencyMatricesFull(n_blocks);
  std::vector<std::vector<std::array<double, 3>>> coordsDestination(n_blocks);
  std::vector<std::vector<std::array<double, 3>>> coordsSource(n_blocks);
  std::vector<int> separatrixCountForEachBlock(n_blocks);
  std::vector<std::vector<int>> edgeOccurenceForEachBlock(n_blocks);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; i++) {
    vtkDataSet *block
      = vtkDataSet::SafeDownCast(multiBlock1_Separatrices->GetBlock(i));

    ttkMorseSmallComplexStability::prepareData(
      block, LocalToGlobal[i], adjacencyMatricesFull[i],
      adjacencyMatricesMinor[i], coordsSource[i], coordsDestination[i],
      separatrixCountForEachBlock[i]);
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; i++) {
    if(!MergeEdgesOnSaddles) {
      status = this->buildOccurenceArraysFull(
        adjacencyMatricesFull, separatrixCountForEachBlock[i], coordsSource,
        coordsDestination, i, edgeOccurenceForEachBlock[i]);
    } else {
      status = this->buildOccurenceArraysMinor(
        adjacencyMatricesMinor, separatrixCountForEachBlock[i],
        coordsDestination, i, edgeOccurenceForEachBlock[i]);
    }
  }

  if(status == 0)
    return status;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < n_blocks; i++) {

    vtkDataSet *block
      = vtkDataSet::SafeDownCast(output1_Separatrices->GetBlock(i));
    int cellNumber = block->GetNumberOfCells();
    vtkCellData *blockCellData = block->GetCellData();
    ttkSimplexIdTypeArray *separatrixIds = ttkSimplexIdTypeArray::SafeDownCast(
      blockCellData->GetArray(ttk::MorseSmaleSeparatrixIdName));

    vtkNew<vtkFloatArray> occurenceCount;
    occurenceCount->SetNumberOfComponents(1);
    occurenceCount->SetName(ttk::MorseSmaleStabilityOccurenceCount);

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
    blockCellData->AddArray(occurenceCount);
  }
  return status;
}

int ttkMorseSmallComplexStability::RequestData(
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

  if(intpu1_Separatrices.GetNumberOfBlocks() < 2) {
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
