#include <ttkDiscreteGradient.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkInformation.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;
using namespace dcg;

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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
int ttkDiscreteGradient::dispatch(vtkUnstructuredGrid *outputCriticalPoints,
                                  const char *sfName) {

  // critical points
  SimplexId criticalPoints_numberOfPoints{};
  vector<float> criticalPoints_points;
  vector<char> criticalPoints_points_cellDimensions;
  vector<SimplexId> criticalPoints_points_cellIds;
  vector<char> criticalPoints_points_isOnBoundary;
  vector<SimplexId> criticalPoints_points_PLVertexIdentifiers;
  vector<SimplexId> criticalPoints_points_manifoldSize;

  int ret = 0;
  vector<VTK_TT> criticalPoints_points_cellScalars;

  this->setOutputCriticalPoints(
    &criticalPoints_numberOfPoints, &criticalPoints_points,
    &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
    &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
    &criticalPoints_points_PLVertexIdentifiers,
    &criticalPoints_points_manifoldSize);

  if(inputOffsets_->GetDataType() == VTK_INT)
    ret = this->buildGradient<VTK_TT, int>();
  if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
    ret = this->buildGradient<VTK_TT, vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    this->printErr("DiscreteGradient.buildGradient() error code: "
                   + std::to_string(ret));
    return -1;
  }
#endif

  // critical points
  {
    this->setCriticalPoints<VTK_TT>();

    vtkNew<vtkPoints> points{};

    vtkNew<vtkCharArray> cellDimensions{};
    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");

    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");

    vtkDataArray *cellScalars = inputScalars_->NewInstance();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellScalars) {
      this->printErr("vtkDataArray allocation problem.");
      return -1;
    }
#endif
    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(sfName);

    vtkNew<vtkCharArray> isOnBoundary{};
    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("IsOnBoundary");

    vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
    PLVertexIdentifiers->SetNumberOfComponents(1);
    PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);

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
    }
    outputCriticalPoints->SetPoints(points);

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
  }

  return ret;
}

int ttkDiscreteGradient::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputGradientGlyphs = vtkUnstructuredGrid::GetData(outputVector, 1);

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

  this->inputScalars_ = this->GetInputArrayToProcess(0, input);
  this->inputOffsets_
    = ttkAlgorithm::GetOptionalArray(this->ForceInputOffsetScalarField, 1,
                                     ttk::OffsetScalarFieldName, inputVector);

  if(this->inputOffsets_ == nullptr) {
    // build a new offset field
    const SimplexId numberOfVertices = input->GetNumberOfPoints();
    offsets_->SetNumberOfComponents(1);
    offsets_->SetNumberOfTuples(numberOfVertices);
    offsets_->SetName(ttk::OffsetScalarFieldName);
    for(SimplexId i = 0; i < numberOfVertices; ++i) {
      offsets_->SetTuple1(i, i);
    }
    this->inputOffsets_ = offsets_;
  }

  if(this->inputScalars_ == nullptr || this->inputOffsets_ == nullptr) {
    this->printErr("Input scalar arrays are NULL");
    return 0;
  }

  int ret{};

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }
#endif

  // baseCode processing
  this->setIterationThreshold(IterationThreshold);
  this->setInputScalarField(inputScalars_->GetVoidPointer(0));
  this->setInputOffsets(inputOffsets_->GetVoidPointer(0));

  switch(inputScalars_->GetDataType()) {
    vtkTemplateMacro(
      ret = dispatch<VTK_TT>(outputCriticalPoints, inputScalars_->GetName()));
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // gradient glyphs
  if(ComputeGradientGlyphs) {
    Timer tm{};

    SimplexId gradientGlyphs_numberOfPoints{};
    vector<float> gradientGlyphs_points;
    vector<char> gradientGlyphs_points_pairOrigins;
    SimplexId gradientGlyphs_numberOfCells{};
    vector<SimplexId> gradientGlyphs_cells;
    vector<char> gradientGlyphs_cells_pairTypes;

    this->setGradientGlyphs(
      gradientGlyphs_numberOfPoints, gradientGlyphs_points,
      gradientGlyphs_points_pairOrigins, gradientGlyphs_numberOfCells,
      gradientGlyphs_cells, gradientGlyphs_cells_pairTypes);

    vtkNew<vtkPoints> points{};
    vtkNew<vtkCharArray> pairOrigins{};
    pairOrigins->SetNumberOfComponents(1);
    pairOrigins->SetName("PairOrigin");
    vtkNew<vtkCharArray> pairTypes{};
    pairTypes->SetNumberOfComponents(1);
    pairTypes->SetName("PairType");

    for(SimplexId i = 0; i < gradientGlyphs_numberOfPoints; ++i) {
      points->InsertNextPoint(gradientGlyphs_points[3 * i],
                              gradientGlyphs_points[3 * i + 1],
                              gradientGlyphs_points[3 * i + 2]);

      pairOrigins->InsertNextTuple1(gradientGlyphs_points_pairOrigins[i]);
    }
    outputGradientGlyphs->SetPoints(points);

    outputGradientGlyphs->Allocate(gradientGlyphs_numberOfCells);
    SimplexId ptr{};
    for(SimplexId i = 0; i < gradientGlyphs_numberOfCells; ++i) {
      std::array<vtkIdType, 2> line
        = {gradientGlyphs_cells[ptr + 1], gradientGlyphs_cells[ptr + 2]};
      outputGradientGlyphs->InsertNextCell(VTK_LINE, 2, line.data());
      pairTypes->InsertNextTuple1(gradientGlyphs_cells_pairTypes[i]);
      ptr += (gradientGlyphs_cells[ptr] + 1);
    }

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
