#include <ttkContinuousScatterPlot.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkContinuousScatterPlot)

  ttkContinuousScatterPlot::ttkContinuousScatterPlot()
  : inputScalars1_{}, inputScalars2_{} {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  ScatterplotResolution[0] = 1920;
  ScatterplotResolution[1] = 1080;
  ScatterplotResolution[2] = 0;

  ProjectImageSupport = true;
  WithDummyValue = false;
  DummyValue = 0;
  UcomponentId = 0;
  VcomponentId = 1;
  triangulation_ = NULL;
  UseAllCores = true;
}

ttkContinuousScatterPlot::~ttkContinuousScatterPlot() {
}

int ttkContinuousScatterPlot::FillInputPortInformation(int port,
                                                       vtkInformation *info) {

  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");

  return 1;
}

int ttkContinuousScatterPlot::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");

  return 1;
}

int ttkContinuousScatterPlot::getScalars(vtkDataSet *input) {
  vtkPointData *pointData = input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData) {
    cerr << "[ttkContinuousScatterPlot] Error : input has no point data."
         << endl;
    return -1;
  }
#endif

  inputScalars1_ = pointData->GetArray(ScalarField1.data());
  if(!inputScalars1_) {
    inputScalars1_ = pointData->GetArray(UcomponentId);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars1_) {
    cerr << "[ttkContinuousScatterPlot] Error : input scalar field 1 pointer "
            "is null."
         << endl;
    return -3;
  }
#endif

  if(inputScalars1_->GetName())
    ScalarField1 = inputScalars1_->GetName();

  inputScalars2_ = pointData->GetArray(ScalarField2.data());
  if(!inputScalars2_) {
    inputScalars2_ = pointData->GetArray(VcomponentId);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars2_) {
    cerr << "[ttkContinuousScatterPlot] Error : input scalar field 2 pointer "
            "is null."
         << endl;
    return -5;
  }
#endif

  if(inputScalars2_->GetName())
    ScalarField2 = inputScalars2_->GetName();

  stringstream msg;
  msg << "[ttkContinuousScatterPlot] U-component `" << inputScalars1_->GetName()
      << "'" << endl;
  msg << "[ttkContinuousScatterPlot] V-component `" << inputScalars2_->GetName()
      << "'" << endl;
  dMsg(cout, msg.str(), infoMsg);

  return 0;
}

int ttkContinuousScatterPlot::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

  if(!triangulation_)
    return -1;

  triangulation_->setWrapper(this);

  return 0;
}
int ttkContinuousScatterPlot::doIt(vector<vtkDataSet *> &inputs,
                                   vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

  int ret{};

  ret = getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  // wrong triangulation
  if(ret) {
    cerr << "[ttkContinuousScatterPlot] Error : wrong triangulation." << endl;
    return -1;
  }
#endif

  ret = getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  // wrong scalar fields
  if(ret) {
    cerr << "[ttkContinuousScatterPlot] Error : wrong scalar fields." << endl;
    return -2;
  }
#endif

  SimplexId numberOfPixels
    = ScatterplotResolution[0] * ScatterplotResolution[1];
#ifndef TTK_ENABLE_KAMIKAZE
  // no pixels
  if(!numberOfPixels) {
    cerr << "[ttkContinuousScatterPlot] Error : no pixels." << endl;
    return -3;
  }
#endif

  density_.clear();
  density_.resize(ScatterplotResolution[0]);
  validPointMask_.clear();
  validPointMask_.resize(ScatterplotResolution[0]);
  for(SimplexId k = 0; k < ScatterplotResolution[0]; ++k) {
    density_[k].resize(ScatterplotResolution[1], 0.0);
    validPointMask_[k].resize(ScatterplotResolution[1], 0);
  }

  SimplexId numberOfPoints = input->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  // no points
  if(!numberOfPoints) {
    cerr << "[ttkContinuousScatterPlot] Error : no points." << endl;
    return -4;
  }
#endif
  for(SimplexId k = 0; k < numberOfPoints; ++k) {
    double d1 = inputScalars1_->GetTuple1(k);
    double d2 = inputScalars2_->GetTuple1(k);

    if(!k or scalarMin_[0] > d1)
      scalarMin_[0] = d1;
    if(!k or scalarMin_[1] > d2)
      scalarMin_[1] = d2;
    if(!k or scalarMax_[0] < d1)
      scalarMax_[0] = d1;
    if(!k or scalarMax_[1] < d2)
      scalarMax_[1] = d2;
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // scalar fields stats problem
  if(scalarMin_[0] == scalarMax_[0] or scalarMin_[1] == scalarMax_[1]) {
    cerr << "[ttkContinuousScatterPlot] Error : scalar fields stats problem."
         << endl;
    return -5;
  }
#endif

  // calling the executing package
  ContinuousScatterPlot continuousScatterPlot;
  continuousScatterPlot.setWrapper(this);
  continuousScatterPlot.setVertexNumber(numberOfPoints);
  continuousScatterPlot.setDummyValue(WithDummyValue, DummyValue);
  continuousScatterPlot.setTriangulation(triangulation_);
  continuousScatterPlot.setResolutions(
    ScatterplotResolution[0], ScatterplotResolution[1]);
  continuousScatterPlot.setInputScalarField1(inputScalars1_->GetVoidPointer(0));
  continuousScatterPlot.setInputScalarField2(inputScalars2_->GetVoidPointer(0));
  continuousScatterPlot.setScalarMin(scalarMin_);
  continuousScatterPlot.setScalarMax(scalarMax_);
  continuousScatterPlot.setOutputDensity(&density_);
  continuousScatterPlot.setOutputMask(&validPointMask_);
  switch(vtkTemplate2PackMacro(
    inputScalars1_->GetDataType(), inputScalars2_->GetDataType())) {
    ttkTemplate2Macro(
      ret = continuousScatterPlot.execute<VTK_T1 TTK_COMMA VTK_T2>());
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in baseCode
  if(ret) {
    cerr << "[ttkContinuousScatterPlot] ContinuousScatterPlot.execute() error "
            "code : "
         << ret << endl;
    return -6;
  }
#endif

  vtkSmartPointer<vtkCharArray> maskScalars
    = vtkSmartPointer<vtkCharArray>::New();
  if(maskScalars) {
    maskScalars->SetNumberOfComponents(1);
    maskScalars->SetNumberOfTuples(numberOfPixels);
    maskScalars->SetName("ValidPointMask");
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  else {
    cerr << "[ttkContinuousScatterPlot] Error detected : vtkCharArray "
            "allocation problem."
         << endl;
    return -9;
  }
#endif

  vtkSmartPointer<vtkDoubleArray> densityScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  if(densityScalars) {
    densityScalars->SetNumberOfComponents(1);
    densityScalars->SetNumberOfTuples(numberOfPixels);
    densityScalars->SetName("Density");
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  else {
    cerr << "[ttkContinuousScatterPlot] Error detected : vtkDoubleArray "
            "allocation problem."
         << endl;
    return -10;
  }
#endif

  vtkSmartPointer<vtkDoubleArray> scalars1
    = vtkSmartPointer<vtkDoubleArray>::New();
  if(scalars1) {
    scalars1->SetNumberOfComponents(1);
    scalars1->SetNumberOfTuples(numberOfPixels);
    scalars1->SetName(ScalarField1.data());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  else {
    cerr << "[ttkContinuousScatterPlot] Error detected : vtkDoubleArray "
            "allocation problem."
         << endl;
    return -11;
  }
#endif

  vtkSmartPointer<vtkDoubleArray> scalars2
    = vtkSmartPointer<vtkDoubleArray>::New();
  if(scalars2) {
    scalars2->SetNumberOfComponents(1);
    scalars2->SetNumberOfTuples(numberOfPixels);
    scalars2->SetName(ScalarField2.data());
  }
#ifndef TTK_ENABLE_KAMIKAZE
  // allocation problem
  else {
    cerr << "[ttkContinuousScatterPlot] Error detected : vtkDoubleArray "
            "allocation problem."
         << endl;
    return -12;
  }
#endif

  double delta[2];
  delta[0] = (scalarMax_[0] - scalarMin_[0]) / (ScatterplotResolution[0] - 1);
  delta[1] = (scalarMax_[1] - scalarMin_[1]) / (ScatterplotResolution[1] - 1);

  double imageMin[2]{0.0, 0.0};
  double imageMax[2]{1.0, 1.0};
  if(ProjectImageSupport) {
    imageMin[0] = scalarMin_[0];
    imageMin[1] = scalarMin_[1];
    imageMax[0] = scalarMax_[0];
    imageMax[1] = scalarMax_[1];
  }
  double imageDelta[2];
  imageDelta[0] = (imageMax[0] - imageMin[0]) / (ScatterplotResolution[0] - 1);
  imageDelta[1] = (imageMax[1] - imageMin[1]) / (ScatterplotResolution[1] - 1);

  vtu_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
  pts_ = vtkSmartPointer<vtkPoints>::New();
  pts_->SetNumberOfPoints(numberOfPixels);

  SimplexId id{};
  vtkIdType ids[3];
  for(SimplexId i = 0; i < ScatterplotResolution[0]; i++) {
    for(SimplexId j = 0; j < ScatterplotResolution[1]; j++) {
      // positions:
      double x = imageMin[0] + i * imageDelta[0];
      double y = imageMin[1] + j * imageDelta[1];
      pts_->SetPoint(id, x, y, 0);

      // scalars:
      // valid mask
      maskScalars->SetTuple1(id, validPointMask_[i][j]);
      // density
      densityScalars->SetTuple1(id, density_[i][j]);
      // original scalar fields
      double d1 = scalarMin_[0] + i * delta[0];
      double d2 = scalarMin_[1] + j * delta[1];
      scalars1->SetTuple1(id, d1);
      scalars2->SetTuple1(id, d2);
      if(i < ScatterplotResolution[0] - 1
         and j < ScatterplotResolution[1] - 1) {
        ids[0] = id;
        ids[1] = id + 1;
        ids[2] = id + ScatterplotResolution[1];
        vtu_->InsertNextCell(VTK_TRIANGLE, 3, ids);

        ids[0] = id + 1;
        ids[1] = id + ScatterplotResolution[1];
        ids[2] = id + ScatterplotResolution[1] + 1;
        vtu_->InsertNextCell(VTK_TRIANGLE, 3, ids);
      }

      ++id;
    }
  }
  vtu_->SetPoints(pts_);
  vtu_->GetPointData()->AddArray(maskScalars);
  vtu_->GetPointData()->AddArray(densityScalars);
  vtu_->GetPointData()->AddArray(scalars1);
  vtu_->GetPointData()->AddArray(scalars2);
  output->ShallowCopy(vtu_);

  {
    stringstream msg;
    msg << "[ttkContinuousScatterPlot] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
