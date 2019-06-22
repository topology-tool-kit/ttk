#include <ttkDiscreteGradient.h>

using namespace std;
using namespace ttk;
using namespace dcg;

vtkStandardNewMacro(ttkDiscreteGradient)

  ttkDiscreteGradient::ttkDiscreteGradient()
  : UseAllCores{true}, ScalarField{}, InputOffsetScalarFieldName{},
    ForceInputOffsetScalarField{false}, ReverseSaddleMaximumConnection{true},
    ReverseSaddleSaddleConnection{true}, AllowSecondPass{true},
    AllowThirdPass{true}, ComputeGradientGlyphs{true}, IterationThreshold{-1},
    ScalarFieldId{}, OffsetFieldId{-1},

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
  vtkUnstructuredGrid *outputCriticalPoints
    = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *outputGradientGlyphs
    = vtkUnstructuredGrid::SafeDownCast(outputs[1]);

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

  // critical points
  SimplexId criticalPoints_numberOfPoints{};
  vector<float> criticalPoints_points;
  vector<char> criticalPoints_points_cellDimensions;
  vector<SimplexId> criticalPoints_points_cellIds;
  vector<char> criticalPoints_points_isOnBoundary;
  vector<SimplexId> criticalPoints_points_PLVertexIdentifiers;
  vector<SimplexId> criticalPoints_points_manifoldSize;

  // gradient pairs
  SimplexId gradientGlyphs_numberOfPoints{};
  vector<float> gradientGlyphs_points;
  vector<char> gradientGlyphs_points_pairOrigins;
  SimplexId gradientGlyphs_numberOfCells{};
  vector<SimplexId> gradientGlyphs_cells;
  vector<char> gradientGlyphs_cells_pairTypes;

  // baseCode processing
  discreteGradient_.setWrapper(this);
  discreteGradient_.setIterationThreshold(IterationThreshold);
  discreteGradient_.setReverseSaddleMaximumConnection(
    ReverseSaddleMaximumConnection);
  discreteGradient_.setReverseSaddleSaddleConnection(
    ReverseSaddleSaddleConnection);
  discreteGradient_.setInputScalarField(inputScalars_->GetVoidPointer(0));
  discreteGradient_.setInputOffsets(inputOffsets_->GetVoidPointer(0));

  discreteGradient_.setOutputGradientGlyphs(
    &gradientGlyphs_numberOfPoints, &gradientGlyphs_points,
    &gradientGlyphs_points_pairOrigins, &gradientGlyphs_numberOfCells,
    &gradientGlyphs_cells, &gradientGlyphs_cells_pairTypes);

  const int dimensionality = triangulation_->getDimensionality();

  switch(inputScalars_->GetDataType()) {
#ifndef _MSC_VER
    vtkTemplateMacro(({
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
        cerr << "[ttkDiscreteGradient] Error : "
                "DiscreteGradient.buildGradient() error code : "
             << ret << endl;
        return -1;
      }
#endif

      if(AllowSecondPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient2<VTK_TT, int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient2<VTK_TT, vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
        if(ret) {
          cerr << "[ttkDiscreteGradient] Error : "
                  "DiscreteGradient.buildGradient2() error code : "
               << ret << endl;
          return -1;
        }
#endif
      }

      if(dimensionality == 3 and AllowThirdPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient3<VTK_TT, int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient3<VTK_TT, vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
        if(ret) {
          cerr << "[ttkDiscreteGradient] Error : "
                  "DiscreteGradient.buildGradient2() error code : "
               << ret << endl;
          return -1;
        }
#endif
      }

      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = discreteGradient_.reverseGradient<VTK_TT, int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = discreteGradient_.reverseGradient<VTK_TT, vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
      if(ret) {
        cerr << "[ttkDiscreteGradient] Error : "
                "DiscreteGradient.reverseGradient() error code : "
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

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!points) {
          cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem."
               << endl;
          return -1;
        }
#endif

        vtkSmartPointer<vtkCharArray> cellDimensions
          = vtkSmartPointer<vtkCharArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!cellDimensions) {
          cerr
            << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
            << endl;
          return -1;
        }
#endif
        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<ttkSimplexIdTypeArray> cellIds
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!cellIds) {
          cerr << "[ttkDiscreteGradient] Error : ttkSimplexIdTypeArray "
                  "allocation problem."
               << endl;
          return -1;
        }
#endif
        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkDataArray *cellScalars = inputScalars_->NewInstance();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!cellScalars) {
          cerr
            << "[ttkDiscreteGradient] Error : vtkDataArray allocation problem."
            << endl;
          return -1;
        }
#endif
        cellScalars->SetNumberOfComponents(1);
        cellScalars->SetName(ScalarField.data());

        vtkSmartPointer<vtkCharArray> isOnBoundary
          = vtkSmartPointer<vtkCharArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!isOnBoundary) {
          cerr
            << "[vtkMorseSmaleComplex] Error : vtkCharArray allocation problem."
            << endl;
          return -1;
        }
#endif
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("IsOnBoundary");

        vtkSmartPointer<ttkSimplexIdTypeArray> PLVertexIdentifiers
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
        if(!PLVertexIdentifiers) {
          cerr << "[ttkMorseSmaleComplex] Error : ttkSimplexIdTypeArray "
                  "allocation "
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

          cellDimensions->InsertNextTuple1(
            criticalPoints_points_cellDimensions[i]);
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
          cerr << "[ttkDiscreteGradient] Error : outputCriticalPoints has no "
                  "point data."
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
    }));
#else
#ifndef TTK_ENABLE_KAMIKAZE
    vtkTemplateMacro({
      vector<VTK_TT> criticalPoints_points_cellScalars;

      discreteGradient_.setOutputCriticalPoints(
        &criticalPoints_numberOfPoints, &criticalPoints_points,
        &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
        &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
        &criticalPoints_points_PLVertexIdentifiers,
        &criticalPoints_points_manifoldSize);

      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = discreteGradient_.buildGradient<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = discreteGradient_.buildGradient<VTK_TT TTK_COMMA vtkIdType>();
      if(ret) {
        cerr << "[ttkDiscreteGradient] Error : "
                "DiscreteGradient.buildGradient() error code : "
             << ret << endl;
        return -1;
      }

      if(AllowSecondPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient2<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient2<VTK_TT TTK_COMMA vtkIdType>();
        if(ret) {
          cerr << "[ttkDiscreteGradient] Error : "
                  "DiscreteGradient.buildGradient2() error code : "
               << ret << endl;
          return -1;
        }
      }

      if(dimensionality == 3 and AllowThirdPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient3<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient3<VTK_TT TTK_COMMA vtkIdType>();
        if(ret) {
          cerr << "[ttkDiscreteGradient] Error : "
                  "DiscreteGradient.buildGradient2() error code : "
               << ret << endl;
          return -1;
        }
      }

      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = discreteGradient_.reverseGradient<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = discreteGradient_.reverseGradient<VTK_TT TTK_COMMA vtkIdType>();
      if(ret) {
        cerr << "[ttkDiscreteGradient] Error : "
                "DiscreteGradient.reverseGradient() error code : "
             << ret << endl;
        return -1;
      }

      // critical points
      {
        if(inputOffsets_->GetDataType() == VTK_INT)
          discreteGradient_.setCriticalPoints<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          discreteGradient_.setCriticalPoints<VTK_TT TTK_COMMA vtkIdType>();

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        if(!points) {
          cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem."
               << endl;
          return -1;
        }

        vtkSmartPointer<vtkCharArray> cellDimensions
          = vtkSmartPointer<vtkCharArray>::New();
        if(!cellDimensions) {
          cerr
            << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
            << endl;
          return -1;
        }
        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<ttkSimplexIdTypeArray> cellIds
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
        if(!cellIds) {
          cerr << "[ttkDiscreteGradient] Error : ttkSimplexIdTypeArray "
                  "allocation problem."
               << endl;
          return -1;
        }
        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkDataArray *cellScalars = inputScalars_->NewInstance();
        if(!cellScalars) {
          cerr
            << "[ttkDiscreteGradient] Error : vtkDataArray allocation problem."
            << endl;
          return -1;
        }
        cellScalars->SetNumberOfComponents(1);
        cellScalars->SetName(ScalarField.data());

        vtkSmartPointer<vtkCharArray> isOnBoundary
          = vtkSmartPointer<vtkCharArray>::New();
        if(!isOnBoundary) {
          cerr
            << "[vtkMorseSmaleComplex] Error : vtkCharArray allocation problem."
            << endl;
          return -1;
        }
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("IsOnBoundary");

        vtkSmartPointer<ttkSimplexIdTypeArray> PLVertexIdentifiers
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
        if(!PLVertexIdentifiers) {
          cerr << "[ttkMorseSmaleComplex] Error : ttkSimplexIdTypeArray "
                  "allocation "
               << "problem." << endl;
          return -1;
        }
        PLVertexIdentifiers->SetNumberOfComponents(1);
        PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);

        for(SimplexId i = 0; i < criticalPoints_numberOfPoints; ++i) {
          points->InsertNextPoint(criticalPoints_points[3 * i],
                                  criticalPoints_points[3 * i + 1],
                                  criticalPoints_points[3 * i + 2]);

          cellDimensions->InsertNextTuple1(
            criticalPoints_points_cellDimensions[i]);
          cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);
          cellScalars->InsertNextTuple1(criticalPoints_points_cellScalars[i]);
          isOnBoundary->InsertNextTuple1(criticalPoints_points_isOnBoundary[i]);
          PLVertexIdentifiers->InsertNextTuple1(
            criticalPoints_points_PLVertexIdentifiers[i]);
        }
        outputCriticalPoints->SetPoints(points);

        vtkPointData *pointData = outputCriticalPoints->GetPointData();
        if(!pointData) {
          cerr << "[ttkDiscreteGradient] Error : outputCriticalPoints has no "
                  "point data."
               << endl;
          return -1;
        }

        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);
        pointData->AddArray(cellScalars);
        pointData->AddArray(isOnBoundary);
        pointData->AddArray(PLVertexIdentifiers);
      }
    });
#else
    vtkTemplateMacro({
      vector<VTK_TT> criticalPoints_points_cellScalars;

      discreteGradient_.setOutputCriticalPoints(
        &criticalPoints_numberOfPoints, &criticalPoints_points,
        &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
        &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
        &criticalPoints_points_PLVertexIdentifiers,
        &criticalPoints_points_manifoldSize);

      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = discreteGradient_.buildGradient<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = discreteGradient_.buildGradient<VTK_TT TTK_COMMA vtkIdType>();

      if(AllowSecondPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient2<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient2<VTK_TT TTK_COMMA vtkIdType>();
      }

      if(dimensionality == 3 and AllowThirdPass) {
        if(inputOffsets_->GetDataType() == VTK_INT)
          ret = discreteGradient_.buildGradient3<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
          ret = discreteGradient_.buildGradient3<VTK_TT TTK_COMMA vtkIdType>();
      }

      if(inputOffsets_->GetDataType() == VTK_INT)
        discreteGradient_.reverseGradient<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        discreteGradient_.reverseGradient<VTK_TT TTK_COMMA vtkIdType>();

      // critical points
      {
        if(inputOffsets_->GetDataType() == VTK_INT)
          discreteGradient_.setCriticalPoints<VTK_TT TTK_COMMA int>();
        if(inputOffsets_->GetDataType() == VTK_INT)
          discreteGradient_.setCriticalPoints<VTK_TT TTK_COMMA vtkIdType>();

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkCharArray> cellDimensions
          = vtkSmartPointer<vtkCharArray>::New();
        cellDimensions->SetNumberOfComponents(1);
        cellDimensions->SetName("CellDimension");

        vtkSmartPointer<ttkSimplexIdTypeArray> cellIds
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
        cellIds->SetNumberOfComponents(1);
        cellIds->SetName("CellId");

        vtkDataArray *cellScalars = inputScalars_->NewInstance();
        cellScalars->SetNumberOfComponents(1);
        cellScalars->SetName(ScalarField.data());

        vtkSmartPointer<vtkCharArray> isOnBoundary
          = vtkSmartPointer<vtkCharArray>::New();
        isOnBoundary->SetNumberOfComponents(1);
        isOnBoundary->SetName("IsOnBoundary");

        vtkSmartPointer<ttkSimplexIdTypeArray> PLVertexIdentifiers
          = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
        PLVertexIdentifiers->SetNumberOfComponents(1);
        PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);

        for(SimplexId i = 0; i < criticalPoints_numberOfPoints; ++i) {
          points->InsertNextPoint(criticalPoints_points[3 * i],
                                  criticalPoints_points[3 * i + 1],
                                  criticalPoints_points[3 * i + 2]);

          cellDimensions->InsertNextTuple1(
            criticalPoints_points_cellDimensions[i]);
          cellIds->InsertNextTuple1(criticalPoints_points_cellIds[i]);
          cellScalars->InsertNextTuple1(criticalPoints_points_cellScalars[i]);
          isOnBoundary->InsertNextTuple1(criticalPoints_points_isOnBoundary[i]);
          PLVertexIdentifiers->InsertNextTuple1(
            criticalPoints_points_PLVertexIdentifiers[i]);
        }
        outputCriticalPoints->SetPoints(points);

        vtkPointData *pointData = outputCriticalPoints->GetPointData();

        pointData->AddArray(cellDimensions);
        pointData->AddArray(cellIds);
        pointData->AddArray(cellScalars);
        pointData->AddArray(isOnBoundary);
        pointData->AddArray(PLVertexIdentifiers);
      }
    });
#endif
#endif
  }

  // gradient glyphs
  if(ComputeGradientGlyphs) {
    discreteGradient_.setGradientGlyphs();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!points) {
      cerr << "[ttkDiscreteGradient] Error : vtkPoints allocation problem."
           << endl;
      return -1;
    }
#endif

    vtkSmartPointer<vtkCharArray> pairOrigins
      = vtkSmartPointer<vtkCharArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pairOrigins) {
      cerr << "[ttkDiscreteGradient] Error : vtkCharArray allocation problem."
           << endl;
      return -1;
    }
#endif
    pairOrigins->SetNumberOfComponents(1);
    pairOrigins->SetName("PairOrigin");

    vtkSmartPointer<vtkCharArray> pairTypes
      = vtkSmartPointer<vtkCharArray>::New();
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

  return 0;
}
