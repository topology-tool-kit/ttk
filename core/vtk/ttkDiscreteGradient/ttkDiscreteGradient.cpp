#include <ttkDiscreteGradient.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>

vtkStandardNewMacro(ttkDiscreteGradient);

ttkDiscreteGradient::ttkDiscreteGradient() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

int ttkDiscreteGradient::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkDiscreteGradient::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkDiscreteGradient::fillCriticalPoints(
  vtkPolyData *outputCriticalPoints,
  vtkDataArray *const inputScalars,
  const triangulationType &triangulation) {

  // critical points
  std::vector<scalarType> criticalPoints_points_cellScalars;
  this->setOutputCriticalPoints(&criticalPoints_points_cellScalars);

  this->setCriticalPoints<scalarType>(triangulation);
  const auto nPoints = outputCriticalPoints_numberOfPoints_;

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

  vtkDataArray *cellScalars = inputScalars->NewInstance();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!cellScalars) {
    this->printErr("vtkDataArray allocation problem.");
    return -1;
  }
#endif
  cellScalars->SetNumberOfComponents(1);
  cellScalars->SetName(inputScalars->GetName());
  cellScalars->SetNumberOfTuples(nPoints);

  vtkNew<vtkSignedCharArray> isOnBoundary{};
  isOnBoundary->SetNumberOfComponents(1);
  isOnBoundary->SetName("IsOnBoundary");
  isOnBoundary->SetNumberOfTuples(nPoints);

  vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
  PLVertexIdentifiers->SetNumberOfComponents(1);
  PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);
  PLVertexIdentifiers->SetNumberOfTuples(nPoints);

  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(nPoints + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(nPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nPoints; ++i) {
    points->SetPoint(i, &outputCriticalPoints_points_[3 * i]);
    cellDimensions->SetTuple1(
      i, outputCriticalPoints_points_cellDimensions_[i]);
    cellIds->SetTuple1(i, outputCriticalPoints_points_cellIds_[i]);
    cellScalars->SetTuple1(i, criticalPoints_points_cellScalars[i]);
    isOnBoundary->SetTuple1(i, outputCriticalPoints_points_isOnBoundary_[i]);
    PLVertexIdentifiers->SetTuple1(
      i, outputCriticalPoints_points_PLVertexIdentifiers_[i]);
    offsets->SetTuple1(i, i);
    connectivity->SetTuple1(i, i);
  }
  offsets->SetTuple1(nPoints, nPoints);

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  outputCriticalPoints->SetPoints(points);
  outputCriticalPoints->SetVerts(cells);

  vtkPointData *pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    this->printErr("outputCriticalPoints has no point data");
    return -1;
  }
#endif

  pointData->AddArray(cellDimensions);
  pointData->AddArray(cellIds);
  pointData->AddArray(cellScalars);
  pointData->AddArray(isOnBoundary);
  pointData->AddArray(PLVertexIdentifiers);

  return 0;
}

int ttkDiscreteGradient::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkPolyData::GetData(outputVector, 0);
  auto outputGradientGlyphs = vtkPolyData::GetData(outputVector, 1);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(!outputCriticalPoints or !outputGradientGlyphs) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }
  if(!input->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
#endif

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is NULL");
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  const auto inputScalars = this->GetInputArrayToProcess(0, input);
  auto inputOffsets = ttkAlgorithm::GetOrderArray(
    input, 0, 1, this->ForceInputOffsetScalarField);

  if(inputScalars == nullptr || inputOffsets == nullptr) {
    this->printErr("Input scalar arrays are NULL");
    return 0;
  }

  int ret{};

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }
#endif

  // baseCode processing
  this->setInputScalarField(ttkUtils::GetVoidPointer(inputScalars));
  this->setInputOffsets(
    static_cast<SimplexId *>(ttkUtils::GetVoidPointer(inputOffsets)));

  ttkTemplateMacro(triangulation->getType(),
                   (ret = this->buildGradient<TTK_TT>(
                      *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0) {
    this->printErr("DiscreteGradient.buildGradient() error code: "
                   + std::to_string(ret));
    return 0;
  }

  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (fillCriticalPoints<VTK_TT, TTK_TT>(
                        outputCriticalPoints, inputScalars,
                        *static_cast<TTK_TT *>(triangulation->getData()))));

  // gradient glyphs
  if(ComputeGradientGlyphs) {
    ttk::Timer tm{};

    SimplexId gradientGlyphs_numberOfPoints{};
    std::vector<float> gradientGlyphs_points;
    std::vector<char> gradientGlyphs_points_pairOrigins;
    SimplexId gradientGlyphs_numberOfCells{};
    std::vector<SimplexId> gradientGlyphs_cells;
    std::vector<char> gradientGlyphs_cells_pairTypes;

    this->setGradientGlyphs(
      gradientGlyphs_numberOfPoints, gradientGlyphs_points,
      gradientGlyphs_points_pairOrigins, gradientGlyphs_numberOfCells,
      gradientGlyphs_cells, gradientGlyphs_cells_pairTypes, *triangulation);

    vtkNew<vtkPoints> points{};
    points->SetNumberOfPoints(gradientGlyphs_numberOfPoints);
    vtkNew<vtkSignedCharArray> pairOrigins{};
    pairOrigins->SetNumberOfComponents(1);
    pairOrigins->SetName("PairOrigin");
    pairOrigins->SetNumberOfTuples(gradientGlyphs_numberOfPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < gradientGlyphs_numberOfPoints; ++i) {
      points->SetPoint(i, &gradientGlyphs_points[3 * i]);
      pairOrigins->SetTuple1(i, gradientGlyphs_points_pairOrigins[i]);
    }
    outputGradientGlyphs->SetPoints(points);

    vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(gradientGlyphs_numberOfCells + 1);
    connectivity->SetNumberOfComponents(1);
    connectivity->SetNumberOfTuples(2 * gradientGlyphs_numberOfCells);
    vtkNew<vtkSignedCharArray> pairTypes{};
    pairTypes->SetNumberOfComponents(1);
    pairTypes->SetName("PairType");
    pairTypes->SetNumberOfTuples(gradientGlyphs_numberOfCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < gradientGlyphs_numberOfCells; ++i) {
      offsets->SetTuple1(i, 2 * i);
      connectivity->SetTuple1(2 * i, gradientGlyphs_cells[3 * i + 1]);
      connectivity->SetTuple1(2 * i + 1, gradientGlyphs_cells[3 * i + 2]);
      pairTypes->SetTuple1(i, gradientGlyphs_cells_pairTypes[i]);
    }
    offsets->SetTuple1(
      gradientGlyphs_numberOfCells, gradientGlyphs_numberOfCells);

    vtkNew<vtkCellArray> cells{};
    cells->SetData(offsets, connectivity);
    outputGradientGlyphs->SetLines(cells);

    vtkPointData *pointData = outputGradientGlyphs->GetPointData();
    vtkCellData *cellData = outputGradientGlyphs->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
    if(pointData == nullptr || cellData == nullptr) {
      this->printErr("In outputGradientGlyphs point or cell data");
      return -1;
    }
#endif

    pointData->AddArray(pairOrigins);
    cellData->AddArray(pairTypes);

    this->printMsg("Computed gradient glyphs", 1.0, tm.getElapsedTime(), 1);
  }

  return 1;
}
