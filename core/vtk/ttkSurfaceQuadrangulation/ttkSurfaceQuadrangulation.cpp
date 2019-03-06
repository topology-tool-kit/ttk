#include <ttkSurfaceQuadrangulation.h>

vtkStandardNewMacro(ttkSurfaceQuadrangulation);

ttkSurfaceQuadrangulation::ttkSurfaceQuadrangulation()
  : InputOffsetIdentifiersFieldName{}, ForceInputIdentifiersField{false},
    SubdivisionLevel{5}, RelaxationIterations{100}, surfaceQuadrangulation_{},
    triangulation_{}, scalarField_{}, offsetIdentifiersField_{} {

  InputScalarFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));

  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));

  SetNumberOfInputPorts(1);
}

int ttkSurfaceQuadrangulation::FillInputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  return 0;
}

int ttkSurfaceQuadrangulation::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_ == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: input triangulation pointer "
            "is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  surfaceQuadrangulation_.setWrapper(this);
  surfaceQuadrangulation_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    cerr << "[ttkSurfaceQuadrangulation] Error: ttkTriangulation allocation "
            "problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkSurfaceQuadrangulation::getScalarField(vtkDataSet *input) {

  if(InputScalarFieldName.length() != 0) {
    scalarField_ = input->GetPointData()->GetArray(InputScalarFieldName.data());
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(scalarField_ == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong scalar field." << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkSurfaceQuadrangulation::getOffsetIdentifiersField(vtkDataSet *input) {

  if(ForceInputIdentifiersField
     && InputOffsetIdentifiersFieldName.length() != 0) {
    offsetIdentifiersField_
      = input->GetPointData()->GetArray(InputOffsetIdentifiersFieldName.data());
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(offsetIdentifiersField_ == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong offset identifiers field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkSurfaceQuadrangulation::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInDomain == 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: domain has no points." << endl;
    return -3;
  }
#endif

  res += getScalarField(domain);
  res += getOffsetIdentifiersField(domain);

  surfaceQuadrangulation_.setVertexNumber(numberOfPointsInDomain);
  surfaceQuadrangulation_.setSubdivisionLevel(SubdivisionLevel);
  surfaceQuadrangulation_.setRelaxationIterations(RelaxationIterations);
  surfaceQuadrangulation_.setInputScalarFieldPointer(scalarField_);
  surfaceQuadrangulation_.setInputOffsetIdentifiersFieldPointer(
    offsetIdentifiersField_);

  std::vector<ttk::SimplexId> outputVertices;
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> outputEdges;

  vtkSmartPointer<vtkDataArray> outputQuadrangulation{};

#ifndef TTK_ENABLE_KAMIKAZE
  if(outputQuadrangulation == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: vtkArray allocation problem."
         << endl;
    return -8;
  }
#endif

  res += surfaceQuadrangulation_.execute();

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/SurfaceQuadrangulation
  if(res != 0) {
    cerr << "[ttkSurfaceQuadrangulation] SurfaceQuadrangulation.execute() "
            "error code: "
         << res << endl;
    return -9;
  }
#endif

  // update result
  outputQuadrangulation->SetNumberOfComponents(1);
  outputQuadrangulation->SetVoidArray(
    static_cast<void *>(outputVertices.data()), outputVertices.size(), 1);

  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(outputQuadrangulation);

  std::stringstream msg;
  msg << "[ttkSurfaceQuadrangulation] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
