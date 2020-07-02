#include <ttkMacros.h>
#include <ttkMorseSmaleComplex.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

vtkStandardNewMacro(ttkMorseSmaleComplex);

ttkMorseSmaleComplex::ttkMorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
  this->SetDebugLevel(3); // bug? debug level is 0 without this statement

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);
}

int ttkMorseSmaleComplex::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmaleComplex::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename scalarType, typename offsetType, typename triangulationType>
int ttkMorseSmaleComplex::dispatch(
  vtkDataArray *const inputScalars,
  vtkDataArray *const inputOffsets,
  vtkUnstructuredGrid *const outputCriticalPoints,
  vtkUnstructuredGrid *const outputSeparatrices1,
  vtkUnstructuredGrid *const outputSeparatrices2,
  const triangulationType &triangulation) {

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;

  // critical points
  SimplexId criticalPoints_numberOfPoints{};
  vector<float> criticalPoints_points;
  vector<char> criticalPoints_points_cellDimensions;
  vector<scalarType> criticalPoints_points_cellScalars;
  vector<SimplexId> criticalPoints_points_cellIds;
  vector<char> criticalPoints_points_isOnBoundary;
  vector<SimplexId> criticalPoints_points_PLVertexIdentifiers;
  vector<SimplexId> criticalPoints_points_manifoldSize;

  // 1-separatrices
  SimplexId separatrices1_numberOfPoints{};
  vector<float> separatrices1_points;
  vector<char> separatrices1_points_smoothingMask;
  vector<char> separatrices1_points_cellDimensions;
  vector<SimplexId> separatrices1_points_cellIds;
  SimplexId separatrices1_numberOfCells{};
  vector<SimplexId> separatrices1_cells;
  vector<SimplexId> separatrices1_cells_sourceIds;
  vector<SimplexId> separatrices1_cells_destinationIds;
  vector<SimplexId> separatrices1_cells_separatrixIds;
  vector<char> separatrices1_cells_separatrixTypes;
  vector<char> separatrices1_cells_isOnBoundary;
  vector<scalarType> separatrices1_cells_separatrixFunctionMaxima;
  vector<scalarType> separatrices1_cells_separatrixFunctionMinima;
  vector<scalarType> separatrices1_cells_separatrixFunctionDiffs;

  // 2-separatrices
  SimplexId separatrices2_numberOfPoints{};
  vector<float> separatrices2_points;
  SimplexId separatrices2_numberOfCells{};
  vector<SimplexId> separatrices2_cells;
  vector<SimplexId> separatrices2_cells_sourceIds;
  vector<SimplexId> separatrices2_cells_separatrixIds;
  vector<char> separatrices2_cells_separatrixTypes;
  vector<char> separatrices2_cells_isOnBoundary;
  vector<scalarType> separatrices2_cells_separatrixFunctionMaxima;
  vector<scalarType> separatrices2_cells_separatrixFunctionMinima;
  vector<scalarType> separatrices2_cells_separatrixFunctionDiffs;

  if(ComputeCriticalPoints) {
    this->setOutputCriticalPoints(
      &criticalPoints_numberOfPoints, &criticalPoints_points,
      &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
      &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  } else {
    this->setOutputCriticalPoints(
      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  this->setOutputSeparatrices1(
    &separatrices1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &separatrices1_numberOfCells,
    &separatrices1_cells, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes,
    &separatrices1_cells_separatrixFunctionMaxima,
    &separatrices1_cells_separatrixFunctionMinima,
    &separatrices1_cells_separatrixFunctionDiffs,
    &separatrices1_cells_isOnBoundary);

  this->setOutputSeparatrices2(
    &separatrices2_numberOfPoints, &separatrices2_points,
    &separatrices2_numberOfCells, &separatrices2_cells,
    &separatrices2_cells_sourceIds, &separatrices2_cells_separatrixIds,
    &separatrices2_cells_separatrixTypes,
    &separatrices2_cells_separatrixFunctionMaxima,
    &separatrices2_cells_separatrixFunctionMinima,
    &separatrices2_cells_separatrixFunctionDiffs,
    &separatrices2_cells_isOnBoundary);

  const int ret
    = this->execute<scalarType, offsetType, triangulationType>(triangulation);

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    this->printErr("MorseSmaleComplex.execute() error");
    return -1;
  }
#endif

  // critical points
  {
    vtkNew<vtkPoints> points{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkSmartPointer<vtkDataArray> cellScalars{inputScalars->NewInstance()};
    vtkNew<vtkSignedCharArray> isOnBoundary{};
    vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
    vtkNew<ttkSimplexIdTypeArray> manifoldSizeScalars{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !cellDimensions || !cellIds || !cellScalars || !isOnBoundary
       || !PLVertexIdentifiers || !manifoldSizeScalars) {
      this->printErr("Critical points vtkDataArray allocation problem.");
      return -1;
    }
#endif

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");

    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(inputScalars->GetName());

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("IsOnBoundary");

    PLVertexIdentifiers->SetNumberOfComponents(1);
    PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);

    manifoldSizeScalars->SetNumberOfComponents(1);
    manifoldSizeScalars->SetName("ManifoldSize");

    for(SimplexId i = 0; i < criticalPoints_numberOfPoints; ++i) {
      points->InsertNextPoint(criticalPoints_points[3 * i],
                              criticalPoints_points[3 * i + 1],
                              criticalPoints_points[3 * i + 2]);

      cellDimensions->InsertNextTuple1(criticalPoints_points_cellDimensions[i]);
      cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);

      cellScalars->InsertNextTuple1(criticalPoints_points_cellScalars[i]);

      isOnBoundary->InsertNextTuple1(criticalPoints_points_isOnBoundary[i]);

      PLVertexIdentifiers->InsertNextTuple1(
        criticalPoints_points_PLVertexIdentifiers[i]);

      if(ComputeAscendingSegmentation and ComputeDescendingSegmentation)
        manifoldSizeScalars->InsertNextTuple1(
          criticalPoints_points_manifoldSize[i]);
      else
        manifoldSizeScalars->InsertNextTuple1(-1);
    }
    outputCriticalPoints->SetPoints(points);

    auto pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputCriticalPoints has no point data.");
      return -1;
    }
#endif

    pointData->AddArray(cellDimensions);
    pointData->AddArray(cellIds);
    pointData->AddArray(cellScalars);
    pointData->AddArray(isOnBoundary);
    pointData->AddArray(PLVertexIdentifiers);
    pointData->AddArray(manifoldSizeScalars);
  }

  // 1-separatrices
  if(ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1
     or ComputeSaddleConnectors) {

    vtkNew<vtkPoints> points{};
    vtkNew<vtkSignedCharArray> smoothingMask{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> destinationIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkSmartPointer<vtkDataArray> separatrixFunctionMaxima{
      inputScalars->NewInstance()};
    vtkSmartPointer<vtkDataArray> separatrixFunctionMinima{
      inputScalars->NewInstance()};
    vtkSmartPointer<vtkDataArray> separatrixFunctionDiffs{
      inputScalars->NewInstance()};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !smoothingMask || !cellDimensions || !cellIds || !sourceIds
       || !destinationIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("1-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif
    smoothingMask->SetNumberOfComponents(1);
    smoothingMask->SetName(ttk::MaskScalarFieldName);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");

    destinationIds->SetNumberOfComponents(1);
    destinationIds->SetName("DestinationId");

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

    for(SimplexId i = 0; i < separatrices1_numberOfPoints; ++i) {
      points->InsertNextPoint(separatrices1_points[3 * i],
                              separatrices1_points[3 * i + 1],
                              separatrices1_points[3 * i + 2]);

      smoothingMask->InsertNextTuple1(separatrices1_points_smoothingMask[i]);
      cellDimensions->InsertNextTuple1(separatrices1_points_cellDimensions[i]);
      cellIds->InsertNextTuple1(separatrices1_points_cellIds[i]);
    }
    outputSeparatrices1->SetPoints(points);

    outputSeparatrices1->Allocate(separatrices1_numberOfCells);
    SimplexId ptr{};
    for(SimplexId i = 0; i < separatrices1_numberOfCells; ++i) {
      vtkIdType line[2];
      line[0] = separatrices1_cells[ptr + 1];
      line[1] = separatrices1_cells[ptr + 2];

      outputSeparatrices1->InsertNextCell(VTK_LINE, 2, line);

      sourceIds->InsertNextTuple1(separatrices1_cells_sourceIds[i]);

      destinationIds->InsertNextTuple1(separatrices1_cells_destinationIds[i]);

      separatrixIds->InsertNextTuple1(separatrices1_cells_separatrixIds[i]);

      separatrixTypes->InsertNextTuple1(separatrices1_cells_separatrixTypes[i]);

      separatrixFunctionMaxima->InsertNextTuple1(
        separatrices1_cells_separatrixFunctionMaxima[i]);

      separatrixFunctionMinima->InsertNextTuple1(
        separatrices1_cells_separatrixFunctionMinima[i]);

      separatrixFunctionDiffs->InsertNextTuple1(
        separatrices1_cells_separatrixFunctionDiffs[i]);

      isOnBoundary->InsertNextTuple1(separatrices1_cells_isOnBoundary[i]);

      ptr += (separatrices1_cells[ptr] + 1);
    }

    auto pointData = outputSeparatrices1->GetPointData();
    auto cellData = outputSeparatrices1->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData || !cellData) {
      this->printErr("outputSeparatrices1 has no point or no cell data.");
      return -1;
    }
#endif

    pointData->AddArray(smoothingMask);
    pointData->AddArray(cellDimensions);
    pointData->AddArray(cellIds);

    cellData->AddArray(sourceIds);
    cellData->AddArray(destinationIds);
    cellData->AddArray(separatrixIds);
    cellData->AddArray(separatrixTypes);
    cellData->AddArray(separatrixFunctionMaxima);
    cellData->AddArray(separatrixFunctionMinima);
    cellData->AddArray(separatrixFunctionDiffs);
    cellData->AddArray(isOnBoundary);
  }

  // 2-separatrices
  if(dimensionality == 3
     and (ComputeAscendingSeparatrices2 or ComputeDescendingSeparatrices2)) {

    vtkNew<vtkPoints> points{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkSmartPointer<vtkDataArray> separatrixFunctionMaxima{
      inputScalars->NewInstance()};
    vtkSmartPointer<vtkDataArray> separatrixFunctionMinima{
      inputScalars->NewInstance()};
    vtkSmartPointer<vtkDataArray> separatrixFunctionDiffs{
      inputScalars->NewInstance()};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !sourceIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("2-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");

    for(SimplexId i = 0; i < separatrices2_numberOfPoints; ++i) {
      points->InsertNextPoint(separatrices2_points[3 * i],
                              separatrices2_points[3 * i + 1],
                              separatrices2_points[3 * i + 2]);
    }
    outputSeparatrices2->SetPoints(points);

    outputSeparatrices2->Allocate(separatrices2_numberOfCells);
    SimplexId ptr{};
    for(SimplexId i = 0; i < separatrices2_numberOfCells; ++i) {
      const int vertexNumber = separatrices2_cells[ptr];

      if(vertexNumber == 3) {
        vtkIdType triangle[3];
        triangle[0] = separatrices2_cells[ptr + 1];
        triangle[1] = separatrices2_cells[ptr + 2];
        triangle[2] = separatrices2_cells[ptr + 3];

        outputSeparatrices2->InsertNextCell(
          VTK_TRIANGLE, vertexNumber, triangle);
      } else {
        vtkIdType ids[16];
        for(int j = 1; j <= vertexNumber; ++j)
          ids[j - 1] = separatrices2_cells[ptr + j];

        outputSeparatrices2->InsertNextCell(VTK_POLYGON, vertexNumber, ids);
      }

      sourceIds->InsertNextTuple1(separatrices2_cells_sourceIds[i]);
      separatrixIds->InsertNextTuple1(separatrices2_cells_separatrixIds[i]);

      separatrixTypes->InsertNextTuple1(separatrices2_cells_separatrixTypes[i]);
      separatrixFunctionMaxima->InsertNextTuple1(
        separatrices2_cells_separatrixFunctionMaxima[i]);

      separatrixFunctionMinima->InsertNextTuple1(
        separatrices2_cells_separatrixFunctionMinima[i]);

      separatrixFunctionDiffs->InsertNextTuple1(
        separatrices2_cells_separatrixFunctionDiffs[i]);
      isOnBoundary->InsertNextTuple1(separatrices2_cells_isOnBoundary[i]);

      ptr += (separatrices2_cells[ptr] + 1);
    }

    auto cellData = outputSeparatrices2->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellData) {
      this->printErr("outputSeparatrices2 has no cell data.");
      return -1;
    }
#endif

    cellData->AddArray(sourceIds);
    cellData->AddArray(separatrixIds);
    cellData->AddArray(separatrixTypes);
    cellData->AddArray(separatrixFunctionMaxima);
    cellData->AddArray(separatrixFunctionMinima);
    cellData->AddArray(separatrixFunctionDiffs);
    cellData->AddArray(isOnBoundary);
  }

  return ret;
}

int ttkMorseSmaleComplex::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  int ret{};

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSeparatrices1 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSeparatrices2 = vtkUnstructuredGrid::GetData(outputVector, 2);
  auto outputMorseComplexes = vtkDataSet::GetData(outputVector, 3);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(input->GetNumberOfPoints() == 0) {
    this->printErr("Input has no point.");
    return -1;
  }
  if(!outputCriticalPoints or !outputSeparatrices1 or !outputSeparatrices2
     or !outputMorseComplexes) {
    this->printErr("Output pointers are NULL.");
    return -1;
  }
#endif

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is null");
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputScalars == nullptr) {
    this->printErr("wrong scalars.");
    return -1;
  }
#endif

  auto inputOffsets
    = ttkAlgorithm::GetOptionalArray(this->ForceInputOffsetScalarField, 1,
                                     ttk::OffsetScalarFieldName, inputVector);

  vtkNew<ttkSimplexIdTypeArray> offsets{};

  if(inputOffsets == nullptr) {
    // build a new offset field
    const SimplexId numberOfVertices = input->GetNumberOfPoints();
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(numberOfVertices);
    offsets->SetName(ttk::OffsetScalarFieldName);
    for(SimplexId i = 0; i < numberOfVertices; ++i) {
      offsets->SetTuple1(i, i);
    }
    inputOffsets = offsets;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets == nullptr) {
    this->printErr("wrong offsets.");
    return -1;
  }
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }
#endif

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");

  // morse complexes
  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfVertices) {
    this->printErr("Input has no vertices.");
    return -1;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> ascendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> descendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> morseSmaleManifold{};
#ifndef TTK_ENABLE_KAMIKAZE
  if(!ascendingManifold || !descendingManifold || !morseSmaleManifold) {
    this->printErr("Manifold vtkDataArray allocation problem.");
    return -1;
  }
#endif
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(numberOfVertices);
  ascendingManifold->SetName("AscendingManifold");

  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName("DescendingManifold");

  morseSmaleManifold->SetNumberOfComponents(1);
  morseSmaleManifold->SetNumberOfTuples(numberOfVertices);
  morseSmaleManifold->SetName("MorseSmaleManifold");

  this->setIterationThreshold(IterationThreshold);

  this->setComputeAscendingSeparatrices1(ComputeAscendingSeparatrices1);
  this->setComputeDescendingSeparatrices1(ComputeDescendingSeparatrices1);
  this->setComputeSaddleConnectors(ComputeSaddleConnectors);

  this->setComputeAscendingSeparatrices2(ComputeAscendingSeparatrices2);
  this->setComputeDescendingSeparatrices2(ComputeDescendingSeparatrices2);

  this->setReturnSaddleConnectors(ReturnSaddleConnectors);
  this->setSaddleConnectorsPersistenceThreshold(
    SaddleConnectorsPersistenceThreshold);

  this->setInputScalarField(inputScalars->GetVoidPointer(0));
  this->setInputOffsets(inputOffsets->GetVoidPointer(0));

  void *ascendingManifoldPtr = nullptr;
  void *descendingManifoldPtr = nullptr;
  void *morseSmaleManifoldPtr = nullptr;
  if(ComputeAscendingSegmentation)
    ascendingManifoldPtr = ascendingManifold->GetVoidPointer(0);
  if(ComputeDescendingSegmentation)
    descendingManifoldPtr = descendingManifold->GetVoidPointer(0);
  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
     and ComputeFinalSegmentation)
    morseSmaleManifoldPtr = morseSmaleManifold->GetVoidPointer(0);

  this->setOutputMorseComplexes(
    ascendingManifoldPtr, descendingManifoldPtr, morseSmaleManifoldPtr);

  if(inputOffsets->GetDataType() == VTK_INT) {
    ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                        (ret = dispatch<VTK_TT, SimplexId, TTK_TT>(
                           inputScalars, inputOffsets, outputCriticalPoints,
                           outputSeparatrices1, outputSeparatrices2,
                           *static_cast<TTK_TT *>(triangulation->getData()))))
  } else if(inputOffsets->GetDataType() == VTK_ID_TYPE) {
    ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                        (ret = dispatch<VTK_TT, ttk::LongSimplexId, TTK_TT>(
                           inputScalars, inputOffsets, outputCriticalPoints,
                           outputSeparatrices1, outputSeparatrices2,
                           *static_cast<TTK_TT *>(triangulation->getData()))))
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  outputMorseComplexes->ShallowCopy(input);
  // morse complexes
  if(ComputeAscendingSegmentation or ComputeDescendingSegmentation) {
    vtkPointData *pointData = outputMorseComplexes->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputMorseComplexes has no point data.");
      return -1;
    }
#endif

    if(ComputeDescendingSegmentation)
      pointData->AddArray(descendingManifold);
    if(ComputeAscendingSegmentation)
      pointData->AddArray(ascendingManifold);
    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
       and ComputeFinalSegmentation)
      pointData->AddArray(morseSmaleManifold);
  }

  return 1;
}
