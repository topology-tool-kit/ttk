#include <ttkIntegralLines.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIntegralLines)

  ttkIntegralLines::ttkIntegralLines()
  : hasUpdatedMesh_{false}, inputScalars_{nullptr}, offsets_{nullptr},
    inputOffsets_{nullptr}, identifiers_{nullptr} {
  Direction = 0;
  SetNumberOfInputPorts(2);
  triangulation_ = NULL;

  OffsetScalarFieldId = -1;
  ForceInputOffsetScalarField = false;
  OffsetScalarFieldName = ttk::OffsetScalarFieldName;
  ForceInputVertexScalarField = false;
  InputVertexScalarFieldName = ttk::VertexScalarFieldName;
  UseAllCores = true;
}

ttkIntegralLines::~ttkIntegralLines() {
  if(offsets_)
    offsets_->Delete();
}

int ttkIntegralLines::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");

  return 1;
}

int ttkIntegralLines::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

int ttkIntegralLines::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

  if(!triangulation_)
    return -1;

  triangulation_->setWrapper(this);
  integralLines_.setWrapper(this);

  integralLines_.setupTriangulation(triangulation_);
  Modified();
  hasUpdatedMesh_ = true;

#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  if(triangulation_->isEmpty()) {
    cerr << "[ttkIntegralLines] Error : ttkTriangulation allocation problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkIntegralLines::getScalars(vtkDataSet *input) {
  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkIntegralLines] Error : input has no point data." << endl;
    return -1;
  }

  if(!ScalarField.length()) {
    cerr << "[ttkIntegralLines] Error : scalar field has no name." << endl;
    return -1;
  }
#endif

  inputScalars_ = pointData->GetArray(ScalarField.data());

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_) {
    cerr << "[ttkIntegralLines] Error : input scalar field pointer is null."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkIntegralLines::getOffsets(vtkDataSet *input) {
  if(ForceInputOffsetScalarField and OffsetScalarFieldName.length()) {
    inputOffsets_
      = input->GetPointData()->GetArray(OffsetScalarFieldName.data());
  } else if(OffsetScalarFieldId != -1
            and input->GetPointData()->GetArray(OffsetScalarFieldId)) {
    inputOffsets_ = input->GetPointData()->GetArray(OffsetScalarFieldId);
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    inputOffsets_ = input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  } else {
    if(hasUpdatedMesh_ and offsets_) {
      offsets_->Delete();
      offsets_ = nullptr;
      hasUpdatedMesh_ = false;
    }

    if(!offsets_) {
      const SimplexId numberOfPoints = input->GetNumberOfPoints();

      offsets_ = ttkSimplexIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfPoints);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(SimplexId i = 0; i < numberOfPoints; ++i)
        offsets_->SetTuple1(i, i);
    }

    inputOffsets_ = offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  if(!inputOffsets_) {
    cerr << "[ttkIntegralLines] Error : wrong input offset scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkIntegralLines::getIdentifiers(vtkPointSet *input) {
  if(ForceInputVertexScalarField and InputVertexScalarFieldName.length())
    identifiers_
      = input->GetPointData()->GetArray(InputVertexScalarFieldName.data());
  else if(input->GetPointData()->GetArray(ttk::VertexScalarFieldName))
    identifiers_ = input->GetPointData()->GetArray(ttk::VertexScalarFieldName);

#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  if(!identifiers_) {
    cerr << "[ttkIntegralLines] Error : wrong input vertex identifier scalar "
            "field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkIntegralLines::getTrajectories(vtkDataSet *input,
                                      vector<vector<SimplexId>> &trajectories,
                                      vtkUnstructuredGrid *output) {
  vtkSmartPointer<vtkUnstructuredGrid> ug
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> dist = vtkSmartPointer<vtkFloatArray>::New();
  dist->SetNumberOfComponents(1);
  dist->SetName("DistanceFromSeed");

  // here, copy the original scalars
  int numberOfArrays = input->GetPointData()->GetNumberOfArrays();

  vector<vtkDataArray *> scalarArrays;
  for(int k = 0; k < numberOfArrays; ++k) {
    auto a = input->GetPointData()->GetArray(k);

    if(a->GetNumberOfComponents() == 1)
      scalarArrays.push_back(a);
  }
  // not efficient, implicit conversion to double
  vector<vtkSmartPointer<vtkDoubleArray>> inputScalars(scalarArrays.size());
  for(unsigned int k = 0; k < scalarArrays.size(); ++k) {
    inputScalars[k] = vtkSmartPointer<vtkDoubleArray>::New();
    inputScalars[k]->SetNumberOfComponents(1);
    inputScalars[k]->SetName(scalarArrays[k]->GetName());
  }

  float p0[3];
  float p1[3];
  vtkIdType ids[2];
  for(SimplexId i = 0; i < (SimplexId)trajectories.size(); ++i) {
    if(trajectories[i].size()) {
      SimplexId vertex = trajectories[i][0];
      // init
      triangulation_->getVertexPoint(vertex, p0[0], p0[1], p0[2]);
      ids[0] = pts->InsertNextPoint(p0);
      // distanceScalars
      float distanceFromSeed{};
      dist->InsertNextTuple1(distanceFromSeed);
      // inputScalars
      for(unsigned int k = 0; k < scalarArrays.size(); ++k)
        inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

      for(SimplexId j = 1; j < (SimplexId)trajectories[i].size(); ++j) {
        vertex = trajectories[i][j];
        triangulation_->getVertexPoint(vertex, p1[0], p1[1], p1[2]);
        ids[1] = pts->InsertNextPoint(p1);
        // distanceScalars
        distanceFromSeed += Geometry::distance(p0, p1, 3);
        dist->InsertNextTuple1(distanceFromSeed);
        // inputScalars
        for(unsigned int k = 0; k < scalarArrays.size(); ++k)
          inputScalars[k]->InsertNextTuple1(scalarArrays[k]->GetTuple1(vertex));

        ug->InsertNextCell(VTK_LINE, 2, ids);

        // iteration
        ids[0] = ids[1];
        p0[0] = p1[0];
        p0[1] = p1[1];
        p0[2] = p1[2];
      }
    }
  }
  ug->SetPoints(pts);
  ug->GetPointData()->AddArray(dist);
  for(unsigned int k = 0; k < scalarArrays.size(); ++k)
    ug->GetPointData()->AddArray(inputScalars[k]);

  output->ShallowCopy(ug);

  return 0;
}

int ttkIntegralLines::doIt(vector<vtkDataSet *> &inputs,
                           vector<vtkDataSet *> &outputs) {
  //                            vtkPointSet* seeds, vtkUnstructuredGrid*
  //                            output){

  Memory m;

  vtkDataSet *domain = inputs[0];
  vtkPointSet *seeds = vtkPointSet::SafeDownCast(inputs[1]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  int ret{};

  ret = getTriangulation(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  // triangulation problem
  if(ret) {
    cerr << "[ttkIntegralLines] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret = getScalars(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  // field problem
  if(ret) {
    cerr << "[ttkIntegralLines] Error : wrong scalar field." << endl;
    return -1;
  }
#endif

  ret = getOffsets(domain);
#ifndef TTK_ENABLE_KAMIKAZE
  // field problem
  if(ret) {
    cerr << "[ttkIntegralLines] Error : wrong offsets." << endl;
    return -1;
  }

  if(inputOffsets_->GetDataType() != VTK_INT
     and inputOffsets_->GetDataType() != VTK_ID_TYPE) {
    cerr << "[ttkIntegralLines] Error : input offset field type not supported."
         << endl;
    return -1;
  }
#endif

  ret = getIdentifiers(seeds);
#ifndef TTK_ENABLE_KAMIKAZE
  // field problem
  if(ret) {
    cerr << "[ttkIntegralLines] Error : wrong identifiers." << endl;
    return -1;
  }
#endif

  const SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  // no points.
  if(numberOfPointsInDomain <= 0) {
    cerr << "[ttkIntegralLines] Error : domain has no points." << endl;
    return -1;
  }
#endif

  const SimplexId numberOfPointsInSeeds = seeds->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  // no points.
  if(numberOfPointsInSeeds <= 0) {
    cerr << "[ttkIntegralLines] Error : seeds have no points." << endl;
    return -1;
  }
#endif

  vector<vector<SimplexId>> trajectories;

  integralLines_.setVertexNumber(numberOfPointsInDomain);
  integralLines_.setSeedNumber(numberOfPointsInSeeds);
  integralLines_.setDirection(Direction);
  integralLines_.setInputScalarField(inputScalars_->GetVoidPointer(0));
  integralLines_.setInputOffsets(inputOffsets_->GetVoidPointer(0));

  integralLines_.setVertexIdentifierScalarField(
    identifiers_->GetVoidPointer(0));
  integralLines_.setOutputTrajectories(&trajectories);

  switch(inputScalars_->GetDataType()) {
    ttkTemplateMacro({
      if(inputOffsets_->GetDataType() == VTK_INT)
        ret = integralLines_.execute<VTK_TT TTK_COMMA int>();
      if(inputOffsets_->GetDataType() == VTK_ID_TYPE)
        ret = integralLines_.execute<VTK_TT TTK_COMMA vtkIdType>();
    });
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    cerr << "[ttkIntegralLines] IntegralLines.execute() error code : " << ret
         << endl;
    return -1;
  }
#endif

  // make the vtk trajectories
  ret = getTrajectories(domain, trajectories, output);
#ifndef TTK_ENABLE_KAMIKAZE
  // trajectories problem
  if(ret) {
    cerr << "[ttkIntegralLines] Error : wrong trajectories." << endl;
    return -1;
  }
#endif

  {
    stringstream msg;
    msg << "[ttkIntegralLines] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
