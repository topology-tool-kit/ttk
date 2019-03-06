#include <ttkSurfaceQuadrangulation.h>

vtkStandardNewMacro(ttkSurfaceQuadrangulation);

ttkSurfaceQuadrangulation::ttkSurfaceQuadrangulation()
  : UseCotanWeights{true}, SolvingMethod{0}, LogAlpha{5}, triangulation_{},
    identifiers_{}, constraints_{} {

  SetNumberOfInputPorts(2);

  InputScalarFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
  OutputScalarFieldName = "OutputSurfaceQuadrangulation";
  OutputScalarFieldType = SurfaceQuadrangulationType::Float;

  ForceConstraintIdentifiers = false;
  ThreadNumber = threadNumber_;
  UseAllCores = true;
  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
}

int ttkSurfaceQuadrangulation::FillInputPortInformation(int port,
                                                        vtkInformation *info) {

  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  }
  return 1;
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

int ttkSurfaceQuadrangulation::getIdentifiers(vtkPointSet *input) {
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);

  if(ForceConstraintIdentifiers && InputIdentifiersFieldName.length() != 0) {
    identifiers_
      = input->GetPointData()->GetArray(InputIdentifiersFieldName.data());
  } else if(input->GetPointData()->GetArray(vsfn) != nullptr) {
    identifiers_ = input->GetPointData()->GetArray(vsfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(identifiers_ == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong vertex identifier "
            "scalar field."
         << endl;
    return -1;
  }
#endif
  return 0;
}

int ttkSurfaceQuadrangulation::getConstraints(vtkPointSet *input) {
  auto osfn = static_cast<const char *>(ttk::OffsetScalarFieldName);

  if(InputScalarFieldName.length() != 0) {
    constraints_ = input->GetPointData()->GetArray(InputScalarFieldName.data());
  } else if(input->GetPointData()->GetArray(osfn) != nullptr) {
    constraints_ = input->GetPointData()->GetArray(osfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(constraints_ == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong constraint "
            "scalar field."
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
  vtkPointSet *identifiers = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  // set this early, since it should trigger some triangulation pre-processing
  surfaceQuadrangulation_.setUseCotanWeights(UseCotanWeights);

  res += getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  res += getIdentifiers(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong identifiers." << endl;
    return -2;
  }
#endif

  res += getConstraints(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: wrong constraints." << endl;
    return -2;
  }
#endif

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInDomain == 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: domain has no points." << endl;
    return -3;
  }
#endif

  const ttk::SimplexId numberOfPointsInSources
    = identifiers->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInSources == 0) {
    cerr << "[ttkSurfaceQuadrangulation] Error: sources have no points."
         << endl;
    return -4;
  }
#endif

  surfaceQuadrangulation_.setVertexNumber(numberOfPointsInDomain);
  surfaceQuadrangulation_.setConstraintNumber(numberOfPointsInSources);
  surfaceQuadrangulation_.setSources(identifiers_->GetVoidPointer(0));
  surfaceQuadrangulation_.setConstraints(constraints_->GetVoidPointer(0));
  surfaceQuadrangulation_.setSolvingMethod(SolvingMethod);
  surfaceQuadrangulation_.setLogAlpha(LogAlpha);

  vtkSmartPointer<vtkDataArray> harmonicScalarField{};

  switch(OutputScalarFieldType) {
    case SurfaceQuadrangulationType::Float:
      harmonicScalarField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case SurfaceQuadrangulationType::Double:
      harmonicScalarField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkSurfaceQuadrangulation] Error: Unknown scalar field type"
           << endl;
      return -7;
#endif // TTK_ENABLE_KAMIKAZE
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(harmonicScalarField == nullptr) {
    cerr << "[ttkSurfaceQuadrangulation] Error: vtkArray allocation problem."
         << endl;
    return -8;
  }
#endif

  harmonicScalarField->SetNumberOfComponents(1);
  harmonicScalarField->SetNumberOfTuples(numberOfPointsInDomain);
  harmonicScalarField->SetName(OutputScalarFieldName.data());

  surfaceQuadrangulation_.setOutputScalarFieldPointer(
    harmonicScalarField->GetVoidPointer(0));

  switch(OutputScalarFieldType) {
    case SurfaceQuadrangulationType::Float:
      res += surfaceQuadrangulation_.execute<float>();
      break;
    case SurfaceQuadrangulationType::Double:
      res += surfaceQuadrangulation_.execute<double>();
      break;
    default:
      break;
  }

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
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(harmonicScalarField);

  std::stringstream msg;
  msg << "[ttkSurfaceQuadrangulation] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
