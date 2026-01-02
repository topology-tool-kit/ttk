#include <ttkDiscreteVectorField.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(ttkDiscreteVectorField);

ttkDiscreteVectorField::ttkDiscreteVectorField() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

int ttkDiscreteVectorField::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkDiscreteVectorField::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkDiscreteVectorField::fillCriticalPoints(
  vtkPolyData *outputCriticalPoints, const triangulationType &triangulation) {

  ttk::Timer tm{};

  // critical points
  std::vector<std::array<float, 3>> critPoints_coords;
  std::vector<char> critPoints_cellDimensions;
  std::vector<SimplexId> critPoints_cellIds;
  std::vector<char> critPoints_isOnBoundary;
  std::vector<SimplexId> critPoints_PLVertexIdentifiers;

  this->setCriticalPoints<scalarType, triangulationType>(
    critPoints_coords, critPoints_cellDimensions, critPoints_cellIds,
    critPoints_isOnBoundary, critPoints_PLVertexIdentifiers, triangulation);
  const auto nPoints = critPoints_coords.size();

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(nPoints);

  vtkNew<vtkSignedCharArray> cellDimensions{};
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");
  cellDimensions->SetNumberOfTuples(nPoints);

  vtkNew<ttkSimplexIdTypeArray> cellIds{};
  cellIds->SetNumberOfComponents(1);
  cellIds->SetName("CellId");
  cellIds->SetNumberOfTuples(nPoints);

  vtkNew<vtkSignedCharArray> isOnBoundary{};
  isOnBoundary->SetNumberOfComponents(1);
  isOnBoundary->SetName("IsOnBoundary");
  isOnBoundary->SetNumberOfTuples(nPoints);

  vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
  PLVertexIdentifiers->SetNumberOfComponents(1);
  PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);
  PLVertexIdentifiers->SetNumberOfTuples(nPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nPoints; ++i) {
    points->SetPoint(i, critPoints_coords[i].data());
    cellDimensions->SetTuple1(i, critPoints_cellDimensions[i]);
    cellIds->SetTuple1(i, critPoints_cellIds[i]);
    isOnBoundary->SetTuple1(i, critPoints_isOnBoundary[i]);
#ifdef TTK_ENABLE_MPI
    if(ttk::hasInitializedMPI()) {
      PLVertexIdentifiers->SetTuple1(
        i, triangulation.getVertexGlobalId(critPoints_PLVertexIdentifiers[i]));
    } else {
      PLVertexIdentifiers->SetTuple1(i, critPoints_PLVertexIdentifiers[i]);
    }
#else
    PLVertexIdentifiers->SetTuple1(i, critPoints_PLVertexIdentifiers[i]);
#endif
  }

  ttkUtils::CellVertexFromPoints(outputCriticalPoints, points);

  vtkPointData *pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    this->printErr("outputCriticalPoints has no point data");
    return -1;
  }
#endif

  pointData->SetScalars(cellDimensions);
  pointData->AddArray(cellIds);
  // pointData->AddArray(cellScalars);
  pointData->AddArray(isOnBoundary);
  pointData->AddArray(PLVertexIdentifiers);

  this->printMsg(
    "Extracted critical points", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttkDiscreteVectorField::fillVectorGlyphs(
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
    "Computed Vector Glyphs", 1.0, tm.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttkDiscreteVectorField::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkPolyData::GetData(outputVector, 0);
  auto outputVectorGlyphs = vtkPolyData::GetData(outputVector, 1);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);

  int keepGoing = checkEmptyMPIInput<ttk::Triangulation>(triangulation);
  if(keepGoing < 2) {
    return keepGoing;
  }
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(!outputCriticalPoints or !outputVectorGlyphs) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }
  if(!input->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
#endif

  this->preconditionTriangulation(triangulation);

  const auto inputVectors = this->GetInputArrayToProcess(0, input);

  if(inputVectors == nullptr) {
    this->printErr("Input vector array is NULL");
    return 0;
  }
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputVectors->GetNumberOfComponents() != 3) {
    this->printErr("Input array needs to be a vector array(3D).");
    return 0;
  }

  int ret{};

  // baseCode processing
  this->setInputVectorField(
    ttkUtils::GetVoidPointer(inputVectors), inputVectors->GetMTime());
#ifdef TTK_ENABLE_MPI_TIME
  ttk::Timer t_mpi;
  ttk::startMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
#endif
  ttkVtkTemplateMacro(inputVectors->GetDataType(), triangulation->getType(),
                      (ret = this->buildField<VTK_TT, TTK_TT>(
                         *static_cast<TTK_TT *>(triangulation->getData()))));
#ifdef TTK_ENABLE_MPI_TIME
  double elapsedTime = ttk::endMPITimer(t_mpi, ttk::MPIrank_, ttk::MPIsize_);
  if(ttk::MPIrank_ == 0) {
    printMsg("Computation performed using " + std::to_string(ttk::MPIsize_)
             + " MPI processes lasted :" + std::to_string(elapsedTime));
  }
#endif
  if(ret != 0) {
    this->printErr("DiscreteVectorField.buildField() error code: "
                   + std::to_string(ret));
    return 0;
  }

  // critical points
  ttkVtkTemplateMacro(
    inputVectors->GetDataType(), triangulation->getType(),
    (fillCriticalPoints<VTK_TT, TTK_TT>(
      outputCriticalPoints, *static_cast<TTK_TT *>(triangulation->getData()))));

  // vector glyphs
  if(ComputeVectorGlyphs) {
    ttkTemplateMacro(
      triangulation->getType(),
      (fillVectorGlyphs<TTK_TT>(
        outputVectorGlyphs, *static_cast<TTK_TT *>(triangulation->getData()))));
  }

  return 1;
}
