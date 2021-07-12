#include <ttkMacros.h>
#include <ttkMorseSmaleComplex.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>

vtkStandardNewMacro(ttkMorseSmaleComplex);

ttkMorseSmaleComplex::ttkMorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename vtkArrayType, typename vectorType>
void setArray(vtkArrayType &vtkArray, vectorType &vector) {
  ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
}

template <typename scalarType, typename triangulationType>
int ttkMorseSmaleComplex::dispatch(vtkDataArray *const inputScalars,
                                   vtkDataArray *const inputOffsets,
                                   vtkPolyData *const outputCriticalPoints,
                                   vtkPolyData *const outputSeparatrices1,
                                   vtkPolyData *const outputSeparatrices2,
                                   const triangulationType &triangulation) {

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;
  const auto scalars
    = static_cast<scalarType *>(ttkUtils::GetVoidPointer(inputScalars));

  // critical points
  criticalPoints_points.clear();
  criticalPoints_points_cellDimensions.clear();
  criticalPoints_points_cellIds.clear();
  criticalPoints_points_isOnBoundary.clear();
  criticalPoints_points_PLVertexIdentifiers.clear();
  criticalPoints_points_manifoldSize.clear();

  // 1-separatrices
  SimplexId s1_numberOfPoints{};
  SimplexId s1_numberOfCells{};
  separatrices1_points.clear();
  separatrices1_points_smoothingMask.clear();
  separatrices1_points_cellDimensions.clear();
  separatrices1_points_cellIds.clear();
  separatrices1_cells_connectivity.clear();
  separatrices1_cells_sourceIds.clear();
  separatrices1_cells_destinationIds.clear();
  separatrices1_cells_separatrixIds.clear();
  separatrices1_cells_separatrixTypes.clear();
  separatrices1_cells_isOnBoundary.clear();
  std::vector<ttk::SimplexId> s1_separatrixFunctionMaxima{};
  std::vector<ttk::SimplexId> s1_separatrixFunctionMinima{};

  // 2-separatrices
  SimplexId s2_numberOfPoints{};
  SimplexId s2_numberOfCells{};
  separatrices2_points.clear();
  separatrices2_cells_offsets.clear();
  separatrices2_cells_connectivity.clear();
  separatrices2_cells_sourceIds.clear();
  separatrices2_cells_separatrixIds.clear();
  separatrices2_cells_separatrixTypes.clear();
  separatrices2_cells_isOnBoundary.clear();
  std::vector<ttk::SimplexId> s2_separatrixFunctionMaxima{};
  std::vector<ttk::SimplexId> s2_separatrixFunctionMinima{};

  if(ComputeCriticalPoints) {
    this->setOutputCriticalPoints(
      &criticalPoints_points, &criticalPoints_points_cellDimensions,
      &criticalPoints_points_cellIds, &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  } else {
    this->setOutputCriticalPoints(
      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  this->setOutputSeparatrices1(
    &s1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &s1_numberOfCells,
    &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes, &s1_separatrixFunctionMaxima,
    &s1_separatrixFunctionMinima, &separatrices1_cells_isOnBoundary);

  this->setOutputSeparatrices2(
    &s2_numberOfPoints, &separatrices2_points, &s2_numberOfCells,
    &separatrices2_cells_offsets, &separatrices2_cells_connectivity,
    &separatrices2_cells_sourceIds, &separatrices2_cells_separatrixIds,
    &separatrices2_cells_separatrixTypes, &s2_separatrixFunctionMaxima,
    &s2_separatrixFunctionMinima, &separatrices2_cells_isOnBoundary);

  const int ret = this->execute<scalarType, triangulationType>(triangulation);

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
    const auto nPoints = criticalPoints_points.size();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !cellDimensions || !cellIds || !cellScalars || !isOnBoundary
       || !PLVertexIdentifiers || !manifoldSizeScalars) {
      this->printErr("Critical points vtkDataArray allocation problem.");
      return -1;
    }
#endif

    points->SetNumberOfPoints(nPoints);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName(ttk::MorseSmaleCellDimensionName);
    setArray(cellDimensions, criticalPoints_points_cellDimensions);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName(ttk::MorseSmaleCellIdName);
    setArray(cellIds, criticalPoints_points_cellIds);

    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(inputScalars->GetName());
    cellScalars->SetNumberOfTuples(nPoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nPoints; ++i) {
      points->SetPoint(i, criticalPoints_points[i].data());
      cellScalars->SetTuple1(
        i, scalars[criticalPoints_points_PLVertexIdentifiers[i]]);
    }

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName(ttk::MorseSmaleBoundaryName);
    setArray(isOnBoundary, criticalPoints_points_isOnBoundary);

    PLVertexIdentifiers->SetNumberOfComponents(1);
    PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);
    setArray(PLVertexIdentifiers, criticalPoints_points_PLVertexIdentifiers);

    manifoldSizeScalars->SetNumberOfComponents(1);
    manifoldSizeScalars->SetName(ttk::MorseSmaleManifoldSizeName);
    if(!ComputeAscendingSegmentation or !ComputeDescendingSegmentation) {
      criticalPoints_points_manifoldSize.resize(nPoints);
      std::fill(criticalPoints_points_manifoldSize.begin(),
                criticalPoints_points_manifoldSize.end(), -1);
    }
    setArray(manifoldSizeScalars, criticalPoints_points_manifoldSize);

    outputCriticalPoints->SetPoints(points);

    vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(nPoints + 1);
    connectivity->SetNumberOfComponents(1);
    connectivity->SetNumberOfTuples(nPoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nPoints; ++i) {
      offsets->SetTuple1(i, i);
      connectivity->SetTuple1(i, i);
    }
    offsets->SetTuple1(nPoints, nPoints);
    vtkNew<vtkCellArray> cells{};
    cells->SetData(offsets, connectivity);
    outputCriticalPoints->SetVerts(cells);

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

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<vtkSignedCharArray> smoothingMask{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> destinationIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkDoubleArray> separatrixFunctionMaxima{};
    vtkNew<vtkDoubleArray> separatrixFunctionMinima{};
    vtkNew<vtkDoubleArray> separatrixFunctionDiffs{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !smoothingMask || !cellDimensions || !cellIds
       || !sourceIds || !destinationIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("1-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices1_points);

    smoothingMask->SetNumberOfComponents(1);
    smoothingMask->SetName(ttk::MaskScalarFieldName);
    setArray(smoothingMask, separatrices1_points_smoothingMask);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName(ttk::MorseSmaleCellDimensionName);
    setArray(cellDimensions, separatrices1_points_cellDimensions);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName(ttk::MorseSmaleCellIdName);
    setArray(cellIds, separatrices1_points_cellIds);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName(ttk::MorseSmaleSourceIdName);
    setArray(sourceIds, separatrices1_cells_sourceIds);

    destinationIds->SetNumberOfComponents(1);
    destinationIds->SetName(ttk::MorseSmaleDestinationIdName);
    setArray(destinationIds, separatrices1_cells_destinationIds);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName(ttk::MorseSmaleSeparatrixIdName);
    setArray(separatrixIds, separatrices1_cells_separatrixIds);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName(ttk::MorseSmaleSeparatrixTypeName);
    setArray(separatrixTypes, separatrices1_cells_separatrixTypes);

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName(ttk::MorseSmaleSeparatrixMaximumName);
    separatrixFunctionMaxima->SetNumberOfTuples(s1_numberOfCells);

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName(ttk::MorseSmaleSeparatrixMinimumName);
    separatrixFunctionMinima->SetNumberOfTuples(s1_numberOfCells);

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName(ttk::MorseSmaleSeparatrixDifferenceName);
    separatrixFunctionDiffs->SetNumberOfTuples(s1_numberOfCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < s1_numberOfCells; ++i) {
      const auto sepId = separatrices1_cells_separatrixIds[i];
      // inputScalars->GetTuple1 not thread safe...
      const auto min = scalars[s1_separatrixFunctionMinima[sepId]];
      const auto max = scalars[s1_separatrixFunctionMaxima[sepId]];
      separatrixFunctionMinima->SetTuple1(i, min);
      separatrixFunctionMaxima->SetTuple1(i, max);
      separatrixFunctionDiffs->SetTuple1(i, max - min);
    }

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName(ttk::MorseSmaleCriticalPointsOnBoundaryName);
    setArray(isOnBoundary, separatrices1_cells_isOnBoundary);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(s1_numberOfCells + 1);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices1_cells_connectivity);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < s1_numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 2 * i);
    }

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices1->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
    cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
    cells->SetData(offsets, connectivity);
    outputSeparatrices1->SetLines(cells);

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

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkDoubleArray> separatrixFunctionMaxima{};
    vtkNew<vtkDoubleArray> separatrixFunctionMinima{};
    vtkNew<vtkDoubleArray> separatrixFunctionDiffs{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !sourceIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("2-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices2_points);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName(ttk::MorseSmaleSourceIdName);
    setArray(sourceIds, separatrices2_cells_sourceIds);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName(ttk::MorseSmaleSeparatrixIdName);
    setArray(separatrixIds, separatrices2_cells_separatrixIds);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName(ttk::MorseSmaleSeparatrixTypeName);
    setArray(separatrixTypes, separatrices2_cells_separatrixTypes);

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName(ttk::MorseSmaleSeparatrixMaximumName);
    separatrixFunctionMaxima->SetNumberOfTuples(s2_numberOfCells);

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName(ttk::MorseSmaleSeparatrixMinimumName);
    separatrixFunctionMinima->SetNumberOfTuples(s2_numberOfCells);

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName(ttk::MorseSmaleSeparatrixDifferenceName);
    separatrixFunctionDiffs->SetNumberOfTuples(s2_numberOfCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < s2_numberOfCells; ++i) {
      const auto sepId = separatrices2_cells_separatrixIds[i];
      // inputScalars->GetTuple1 not thread safe...
      const auto min = scalars[s2_separatrixFunctionMinima[sepId]];
      const auto max = scalars[s2_separatrixFunctionMaxima[sepId]];
      separatrixFunctionMinima->SetTuple1(i, min);
      separatrixFunctionMaxima->SetTuple1(i, max);
      separatrixFunctionDiffs->SetTuple1(i, max - min);
    }

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName(ttk::MorseSmaleCriticalPointsOnBoundaryName);
    setArray(isOnBoundary, separatrices2_cells_isOnBoundary);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    setArray(offsets, separatrices2_cells_offsets);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices2_cells_connectivity);

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices2->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
    cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
    cells->SetData(offsets, connectivity);
    outputSeparatrices2->SetPolys(cells);

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

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkPolyData::GetData(outputVector, 0);
  auto outputSeparatrices1 = vtkPolyData::GetData(outputVector, 1);
  auto outputSeparatrices2 = vtkPolyData::GetData(outputVector, 2);
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

  auto inputOffsets = ttkAlgorithm::GetOrderArray(
    input, 0, 1, this->ForceInputOffsetScalarField);

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
  ascendingManifold->SetName(ttk::MorseSmaleAscendingName);

  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName(ttk::MorseSmaleDescendingName);

  morseSmaleManifold->SetNumberOfComponents(1);
  morseSmaleManifold->SetNumberOfTuples(numberOfVertices);
  morseSmaleManifold->SetName(ttk::MorseSmaleManifoldName);

  this->setIterationThreshold(IterationThreshold);

  this->setComputeAscendingSeparatrices1(ComputeAscendingSeparatrices1);
  this->setComputeDescendingSeparatrices1(ComputeDescendingSeparatrices1);
  this->setComputeSaddleConnectors(ComputeSaddleConnectors);

  this->setComputeAscendingSeparatrices2(ComputeAscendingSeparatrices2);
  this->setComputeDescendingSeparatrices2(ComputeDescendingSeparatrices2);

  this->setReturnSaddleConnectors(ReturnSaddleConnectors);
  this->setSaddleConnectorsPersistenceThreshold(
    SaddleConnectorsPersistenceThreshold);

  this->setInputScalarField(ttkUtils::GetVoidPointer(inputScalars));
  this->setInputOffsets(
    static_cast<SimplexId *>(ttkUtils::GetVoidPointer(inputOffsets)));

  void *ascendingManifoldPtr = nullptr;
  void *descendingManifoldPtr = nullptr;
  void *morseSmaleManifoldPtr = nullptr;
  if(ComputeAscendingSegmentation)
    ascendingManifoldPtr = ttkUtils::GetVoidPointer(ascendingManifold);
  if(ComputeDescendingSegmentation)
    descendingManifoldPtr = ttkUtils::GetVoidPointer(descendingManifold);
  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
     and ComputeFinalSegmentation)
    morseSmaleManifoldPtr = ttkUtils::GetVoidPointer(morseSmaleManifold);

  this->setOutputMorseComplexes(
    ascendingManifoldPtr, descendingManifoldPtr, morseSmaleManifoldPtr);

  int ret{};

  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    (ret = dispatch<VTK_TT, TTK_TT>(
       inputScalars, inputOffsets, outputCriticalPoints, outputSeparatrices1,
       outputSeparatrices2, *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0) {
    return -1;
  }

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

    pointData->AddArray(inputOffsets);
  }

  return !ret;
}
