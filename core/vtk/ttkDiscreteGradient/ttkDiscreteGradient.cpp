#include <ttkDiscreteGradient.h>

using namespace std;
using namespace ttk;
using namespace dcg;

vtkStandardNewMacro(ttkDiscreteGradient)

  ttkDiscreteGradient::ttkDiscreteGradient()
  : UseAllCores{true}, ScalarField{}, InputOffsetScalarFieldName{},
    ForceInputOffsetScalarField{false}, ComputeGradientGlyphs{true},
    IterationThreshold{-1}, ScalarFieldId{}, OffsetFieldId{-1},

    triangulation_{}, inputScalars_{}, offsets_{}, inputOffsets_{},
    hasUpdatedMesh_{} {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

ttkDiscreteGradient::~ttkDiscreteGradient() {
  if(offsets_)
    offsets_->Delete();
}

int ttkDiscreteGradient::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

int ttkDiscreteGradient::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  switch(port) {
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
  }

  return 1;
}

int ttkDiscreteGradient::setupTriangulation(vtkDataSet *input) {
  triangulation_ = ttkTriangulation::getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation_) {
    cerr << "[ttkDiscreteGradient] Error : "
            "ttkTriangulation::getTriangulation() is null."
         << endl;
    return -1;
  }
#endif

  hasUpdatedMesh_
    = ttkTriangulation::hasChangedConnectivity(triangulation_, input, this);

  triangulation_->setWrapper(this);
  discreteGradient_.setWrapper(this);
  discreteGradient_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    cerr << "[ttkDiscreteGradient] Error : vtkTriangulation allocation problem."
         << endl;
    return -1;
  }
#endif
  return 0;
}

int ttkDiscreteGradient::getScalars(vtkDataSet *input) {
  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkDiscreteGradient] Error : input has no point data." << endl;
    return -1;
  }

  if(!ScalarField.length()) {
    cerr << "[ttkDiscreteGradient] Error : scalar field has no name." << endl;
    return -2;
  }
#endif

  if(ScalarField.length()) {
    inputScalars_ = pointData->GetArray(ScalarField.data());
  } else {
    inputScalars_ = pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    cerr << "[ttkDiscreteGradient] Error : input scalar field pointer is null."
         << endl;
    return -3;
  }
#endif

  return 0;
}

int ttkDiscreteGradient::getOffsets(vtkDataSet *input) {
  if(OffsetFieldId != -1) {
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets_) {
      InputOffsetScalarFieldName = inputOffsets_->GetName();
      ForceInputOffsetScalarField = true;
    }
  }

  if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
    inputOffsets_
      = input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_ = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(hasUpdatedMesh_ and offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
    }

    if(!offsets_) {
      const SimplexId numberOfVertices = input->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(SimplexId i = 0; i < numberOfVertices; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_) {
    cerr << "[ttkDiscreteGradient] Error : wrong input offset scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

template <typename VTK_TT>
int ttkDiscreteGradient::dispatch(vtkUnstructuredGrid *outputCriticalPoints) {

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

  discreteGradient_.setOutputCriticalPoints(
    &criticalPoints_numberOfPoints, &criticalPoints_points,
    &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
    &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
    &criticalPoints_points_PLVertexIdentifiers,
    &criticalPoints_points_manifoldSize);

  if(inputOffsets_->GetDataType() == VTK_INT)
    ret = discreteGradient_.buildGradient<VTK_TT, int>();
  if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
    ret = discreteGradient_.buildGradient<VTK_TT, vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDiscreteGradient] Error : DiscreteGradient.buildGradient() "
            "error code : "
         << ret << endl;
    return -1;
  }
#endif

  // critical points
  {
    if(inputOffsets_->GetDataType() == VTK_INT)
      discreteGradient_.setCriticalPoints<VTK_TT, int>();
    if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
      discreteGradient_.setCriticalPoints<VTK_TT, vtkIdType>();

    vtkNew<vtkPoints> points{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!points) {
      cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem."
           << endl;
      return -1;
    }
#endif

    vtkNew<vtkCharArray> cellDimensions{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellDimensions) {
      cerr << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
           << endl;
      return -1;
    }
#endif
    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");

    vtkNew<ttkSimplexIdTypeArray> cellIds{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellIds) {
      cerr << "[ttkDiscreteGradient] Error : ttkSimplexIdTypeArray allocation "
              "problem."
           << endl;
      return -1;
    }
#endif
    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");

    vtkDataArray *cellScalars = inputScalars_->NewInstance();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellScalars) {
      cerr << "[ttkDiscreteGradient] Error : vtkDataArray allocation problem."
           << endl;
      return -1;
    }
#endif
    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(ScalarField.data());

    vtkNew<vtkCharArray> isOnBoundary{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!isOnBoundary) {
      cerr << "[vtkMorseSmaleComplex] Error : vtkCharArray allocation problem."
           << endl;
      return -1;
    }
#endif
    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("IsOnBoundary");

    vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!PLVertexIdentifiers) {
      cerr << "[ttkMorseSmaleComplex] Error : ttkSimplexIdTypeArray allocation "
           << "problem." << endl;
      return -1;
    }
#endif
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
      cerr << "[ttkDiscreteGradient] Error : outputCriticalPoints has no point "
              "data."
           << endl;
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

int ttkDiscreteGradient::doIt(vector<vtkDataSet *> &inputs,
                              vector<vtkDataSet *> &outputs) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr << "[ttkDiscreteGradient] Error: not enough input information."
         << endl;
    return -1;
  }
#endif
  vtkDataSet *input = inputs[0];
  auto outputCriticalPoints = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  auto outputGradientGlyphs = vtkUnstructuredGrid::SafeDownCast(outputs[1]);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkDiscreteGradient] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!outputCriticalPoints or !outputGradientGlyphs) {
    cerr << "[ttkDiscreteGradient] Error: output pointer is NULL." << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkDiscreteGradient] Error: input has no point." << endl;
    return -1;
  }
#endif

  int ret{};

  ret = setupTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDiscreteGradient] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret = getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDiscreteGradient] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  ret = getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret) {
    cerr << "[ttkDiscreteGradient] Error : wrong offsets." << endl;
    return -1;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    cerr
      << "[ttkDiscreteGradient] Error : input offset field type not supported."
      << endl;
    return -1;
  }
#endif

  // baseCode processing
  discreteGradient_.setWrapper(this);
  discreteGradient_.setIterationThreshold(IterationThreshold);
  discreteGradient_.setInputScalarField(inputScalars_->GetVoidPointer(0));
  discreteGradient_.setInputOffsets(inputOffsets_->GetVoidPointer(0));

  switch(inputScalars_->GetDataType()) {
    vtkTemplateMacro(ret = dispatch<VTK_TT>(outputCriticalPoints));
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // gradient glyphs
  if(ComputeGradientGlyphs) {
    SimplexId gradientGlyphs_numberOfPoints{};
    vector<float> gradientGlyphs_points;
    vector<char> gradientGlyphs_points_pairOrigins;
    SimplexId gradientGlyphs_numberOfCells{};
    vector<SimplexId> gradientGlyphs_cells;
    vector<char> gradientGlyphs_cells_pairTypes;

    discreteGradient_.setGradientGlyphs(
      gradientGlyphs_numberOfPoints, gradientGlyphs_points,
      gradientGlyphs_points_pairOrigins, gradientGlyphs_numberOfCells,
      gradientGlyphs_cells, gradientGlyphs_cells_pairTypes);

    vtkNew<vtkPoints> points{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!points) {
      cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem."
           << endl;
      return -1;
    }
#endif

    vtkNew<vtkCharArray> pairOrigins{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pairOrigins) {
      cerr << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
           << endl;
      return -1;
    }
#endif
    pairOrigins->SetNumberOfComponents(1);
    pairOrigins->SetName("PairOrigin");

    vtkNew<vtkCharArray> pairTypes{};
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pairTypes) {
      cerr << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
           << endl;
      return -1;
    }
#endif
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
      vtkIdType line[2];
      line[0] = gradientGlyphs_cells[ptr + 1];
      line[1] = gradientGlyphs_cells[ptr + 2];

      outputGradientGlyphs->InsertNextCell(VTK_LINE, 2, line);

      pairTypes->InsertNextTuple1(gradientGlyphs_cells_pairTypes[i]);

      ptr += (gradientGlyphs_cells[ptr] + 1);
    }

    vtkPointData *pointData = outputGradientGlyphs->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      cerr << "[ttkDiscreteGradient] Error : outputGradientGlyphs has no point "
              "data."
           << endl;
      return -1;
    }
#endif

    pointData->AddArray(pairOrigins);

    vtkCellData *cellData = outputGradientGlyphs->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellData) {
      cerr << "[ttkDiscreteGradient] Error : outputGradientGlyphs has no cell "
              "data."
           << endl;
      return -1;
    }
#endif

    cellData->AddArray(pairTypes);
  }

  return ret;
}
