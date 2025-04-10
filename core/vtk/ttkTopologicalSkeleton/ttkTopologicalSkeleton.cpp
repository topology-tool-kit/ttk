#include <ttkMacros.h>
#include <ttkTopologicalSkeleton.h>
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

vtkStandardNewMacro(ttkTopologicalSkeleton);

ttkTopologicalSkeleton::ttkTopologicalSkeleton() {
  this->setDebugMsgPrefix("TopologicalSkeleton");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(5);
}

int ttkTopologicalSkeleton::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTopologicalSkeleton::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2 || port == 4) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename triangulationType>
int ttkTopologicalSkeleton::fillVectorGlyphs(
  vtkPolyData *const outputVectorGlyphs,
  const triangulationType &triangulation) {

  ttk::Timer tm{};

  std::vector<std::array<float, 3>> vectorGlyphs_points;
  std::vector<char> vectorGlyphs_points_pairOrigins;
  std::vector<char> vectorGlyphs_cells_pairTypes;
  std::vector<SimplexId> vectorGlyphs_point_ids{};
  std::vector<char> vectorGlyphs_point_dimensions{};

  this->setVectorGlyphs(vectorGlyphs_points, vectorGlyphs_points_pairOrigins,
                        vectorGlyphs_cells_pairTypes, vectorGlyphs_point_ids,
                        vectorGlyphs_point_dimensions, triangulation);

  const auto nPoints = vectorGlyphs_points.size();

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(nPoints);
  vtkNew<vtkSignedCharArray> pairOrigins{};
  pairOrigins->SetNumberOfComponents(1);
  pairOrigins->SetName("PairOrigin");
  pairOrigins->SetNumberOfTuples(nPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nPoints; ++i) {
    points->SetPoint(i, vectorGlyphs_points[i].data());
    pairOrigins->SetTuple1(i, vectorGlyphs_points_pairOrigins[i]);
  }
  outputVectorGlyphs->SetPoints(points);

  const auto nCells = vectorGlyphs_cells_pairTypes.size();

  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nCells + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * nCells);
  vtkNew<vtkSignedCharArray> pairTypes{};
  pairTypes->SetNumberOfComponents(1);
  pairTypes->SetName("PairType");
  pairTypes->SetNumberOfTuples(nCells);
  vtkNew<ttkSimplexIdTypeArray> cellIds{};
  cellIds->SetNumberOfComponents(1);
  cellIds->SetName("CellId");
  cellIds->SetNumberOfTuples(2 * nCells);
  vtkNew<vtkSignedCharArray> cellDimensions{};
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");
  cellDimensions->SetNumberOfTuples(2 * nCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCells; ++i) {
    offsets->SetTuple1(i, 2 * i);
    // each glyph/line has unique points
    connectivity->SetTuple1(2 * i, 2 * i);
    connectivity->SetTuple1(2 * i + 1, 2 * i + 1);
    pairTypes->SetTuple1(i, vectorGlyphs_cells_pairTypes[i]);
    cellIds->SetTuple1(2 * i + 0, vectorGlyphs_point_ids[2 * i + 0]);
    cellIds->SetTuple1(2 * i + 1, vectorGlyphs_point_ids[2 * i + 1]);
    cellDimensions->SetTuple1(
      2 * i + 0, vectorGlyphs_point_dimensions[2 * i + 0]);
    cellDimensions->SetTuple1(
      2 * i + 1, vectorGlyphs_point_dimensions[2 * i + 1]);
  }
  offsets->SetTuple1(nCells, connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  outputVectorGlyphs->SetLines(cells);

  vtkPointData *pointData = outputVectorGlyphs->GetPointData();
  vtkCellData *cellData = outputVectorGlyphs->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(pointData == nullptr || cellData == nullptr) {
    this->printErr("In outputVectorGlyphs point or cell data");
    return -1;
  }
#endif

  pointData->AddArray(pairOrigins);
  pointData->AddArray(cellIds);
  pointData->AddArray(cellDimensions);
  cellData->SetScalars(pairTypes);

  this->printMsg(
    "Computed vector glyphs", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename vtkArrayType, typename vectorType>
void setArray(vtkArrayType &vtkArray, vectorType &vector) {
  ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
}

template <typename scalarType, typename triangulationType>
int ttkTopologicalSkeleton::dispatch(vtkDataArray *const inputVectors,
                                     vtkPolyData *const outputCriticalPoints,
                                     vtkPolyData *const outputSeparatrices1,
                                     vtkPolyData *const outputSeparatrices2,
                                     const triangulationType &triangulation) {

  const int dimensionality = triangulation.getDimensionality();
  const auto vectors = ttkUtils::GetPointer<scalarType>(inputVectors);

  const int ret = this->execute(criticalPoints_, separatrices1_, separatrices2_,
                                segmentations_, vectors,
                                inputVectors->GetMTime(), triangulation);

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    this->printErr("TopologicalSkeleton.execute() error");
    return -1;
  }
#endif

  // critical points
  {
    vtkNew<vtkPoints> points{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};
    vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
    vtkNew<ttkSimplexIdTypeArray> manifoldSizeScalars{};
    const auto nPoints = criticalPoints_.points_.size();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !cellDimensions || !cellIds || !isOnBoundary
       || !PLVertexIdentifiers || !manifoldSizeScalars) {
      this->printErr("Critical points vtkDataArray allocation problem.");
      return -1;
    }
#endif

    points->SetNumberOfPoints(nPoints);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    setArray(cellDimensions, criticalPoints_.cellDimensions_);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    setArray(cellIds, criticalPoints_.cellIds_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nPoints; ++i) {
      points->SetPoint(i, criticalPoints_.points_[i].data());
    }

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("IsOnBoundary");
    setArray(isOnBoundary, criticalPoints_.isOnBoundary_);

    PLVertexIdentifiers->SetNumberOfComponents(1);
    PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);
    setArray(PLVertexIdentifiers, criticalPoints_.PLVertexIdentifiers_);

    manifoldSizeScalars->SetNumberOfComponents(1);
    manifoldSizeScalars->SetName("ManifoldSize");
    if(!ComputeAscendingSegmentation or !ComputeDescendingSegmentation) {
      criticalPoints_.manifoldSize_.resize(nPoints);
      std::fill(criticalPoints_.manifoldSize_.begin(),
                criticalPoints_.manifoldSize_.end(), -1);
    }
    setArray(manifoldSizeScalars, criticalPoints_.manifoldSize_);

    ttkUtils::CellVertexFromPoints(outputCriticalPoints, points);

    auto pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputCriticalPoints has no point data.");
      return -1;
    }
#endif

    pointData->SetScalars(cellDimensions);
    pointData->AddArray(cellIds);
    pointData->AddArray(isOnBoundary);
    pointData->AddArray(PLVertexIdentifiers);
    pointData->AddArray(manifoldSizeScalars);
  }

  // 1-separatrices
  if(ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1
     or ComputeSaddleConnectors or ComputeAttractingCycles1
     or ComputeRepellingCycles1) {

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<vtkSignedCharArray> smoothingMask{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> destinationIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};
    vtkNew<vtkSignedCharArray> isOrbit{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !smoothingMask || !cellDimensions || !cellIds
       || !sourceIds || !destinationIds || !separatrixIds || !separatrixTypes
       || !isOnBoundary || !isOrbit) {
      this->printErr("1-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices1_.pt.points_);

    smoothingMask->SetNumberOfComponents(1);
    smoothingMask->SetName(ttk::MaskScalarFieldName);
    setArray(smoothingMask, separatrices1_.pt.smoothingMask_);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    setArray(cellDimensions, separatrices1_.pt.cellDimensions_);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    setArray(cellIds, separatrices1_.pt.cellIds_);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");
    setArray(sourceIds, separatrices1_.cl.sourceIds_);

    destinationIds->SetNumberOfComponents(1);
    destinationIds->SetName("DestinationId");
    setArray(destinationIds, separatrices1_.cl.destinationIds_);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");
    setArray(separatrixIds, separatrices1_.cl.separatrixIds_);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");
    setArray(separatrixTypes, separatrices1_.cl.separatrixTypes_);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");
    setArray(isOnBoundary, separatrices1_.cl.isOnBoundary_);

    isOrbit->SetNumberOfComponents(1);
    isOrbit->SetName("IsAnOrbit");
    setArray(isOrbit, separatrices1_.cl.isOrbit_);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(separatrices1_.cl.numberOfCells_ + 1);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices1_.cl.connectivity_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < separatrices1_.cl.numberOfCells_ + 1; ++i) {
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
    cellData->SetScalars(separatrixTypes);
    cellData->AddArray(isOnBoundary);
    cellData->AddArray(isOrbit);
  }

  // 2-separatrices
  if(dimensionality == 3
     and (ComputeAscendingSeparatrices2 or ComputeDescendingSeparatrices2)) {

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !sourceIds || !separatrixIds || !separatrixTypes
       || !isOnBoundary) {
      this->printErr("2-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices2_.pt.points_);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");
    setArray(sourceIds, separatrices2_.cl.sourceIds_);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");
    setArray(separatrixIds, separatrices2_.cl.separatrixIds_);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");
    setArray(separatrixTypes, separatrices2_.cl.separatrixTypes_);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");
    setArray(isOnBoundary, separatrices2_.cl.isOnBoundary_);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    setArray(offsets, separatrices2_.cl.offsets_);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices2_.cl.connectivity_);

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
    cellData->AddArray(isOnBoundary);
  }

  return ret;
}

int ttkTopologicalSkeleton::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkPolyData::GetData(outputVector, 0);
  auto outputSeparatrices1 = vtkPolyData::GetData(outputVector, 1);
  auto outputSeparatrices2 = vtkPolyData::GetData(outputVector, 2);
  auto outputMorseComplexes = vtkDataSet::GetData(outputVector, 3);
  auto outputVectorGlyphs = vtkPolyData::GetData(outputVector, 4);

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
     or !outputMorseComplexes or !outputVectorGlyphs) {
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

  const auto inputValues = this->GetInputArrayToProcess(0, inputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputValues == nullptr) {
    this->printErr("wrong vectors.");
    return -1;
  }
#endif

  this->printMsg("Launching computation on field `"
                 + std::string(inputValues->GetName()) + "'...");

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
  morseSmaleManifold->SetName("IntersectingManifold");

  this->segmentations_ = {ttkUtils::GetPointer<SimplexId>(ascendingManifold),
                          ttkUtils::GetPointer<SimplexId>(descendingManifold),
                          ttkUtils::GetPointer<SimplexId>(morseSmaleManifold)};

  this->setRunSimplification(RunSimplification);
  this->setSimplificationThreshold(SimplificationThreshold);
  this->setFullOrbits(ReverseFullOrbit);

  int ret{};

  ttkVtkTemplateMacro(
    inputValues->GetDataType(), triangulation->getType(),
    (ret = dispatch<VTK_TT, TTK_TT>(
       inputValues, outputCriticalPoints, outputSeparatrices1,
       outputSeparatrices2, *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0) {
    return -1;
  }

  if(ComputeVectorGlyphs) {
    ttkTemplateMacro(
      triangulation->getType(),
      (fillVectorGlyphs<TTK_TT>(
        outputVectorGlyphs, *static_cast<TTK_TT *>(triangulation->getData()))));
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
  }

  return !ret;
}
