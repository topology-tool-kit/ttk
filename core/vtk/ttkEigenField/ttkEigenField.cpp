#include <ttkEigenField.h>

vtkStandardNewMacro(ttkEigenField);

ttkEigenField::ttkEigenField()
  : UseCotanWeights{true}, SolvingMethod{0}, LogAlpha{5}, triangulation_{},
    identifiers_{}, constraints_{} {

  SetNumberOfInputPorts(2);

  InputScalarFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
  OutputScalarFieldName = "OutputEigenField";
  OutputScalarFieldType = EigenFieldType::Float;

  ForceConstraintIdentifiers = false;
  ThreadNumber = threadNumber_;
  UseAllCores = true;
  InputIdentifiersFieldName
    = std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
}

int ttkEigenField::FillInputPortInformation(int port, vtkInformation *info) {

  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  }
  return 1;
}

int ttkEigenField::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_ == nullptr) {
    cerr << "[ttkEigenField] Error: input triangulation pointer "
            "is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  eigenField_.setWrapper(this);
  eigenField_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangulation_->isEmpty()) {
    cerr << "[ttkEigenField] Error: ttkTriangulation allocation "
            "problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkEigenField::getIdentifiers(vtkPointSet *input) {
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);

  if(ForceConstraintIdentifiers && InputIdentifiersFieldName.length() != 0) {
    identifiers_
      = input->GetPointData()->GetArray(InputIdentifiersFieldName.data());
  } else if(input->GetPointData()->GetArray(vsfn) != nullptr) {
    identifiers_ = input->GetPointData()->GetArray(vsfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(identifiers_ == nullptr) {
    cerr << "[ttkEigenField] Error: wrong vertex identifier "
            "scalar field."
         << endl;
    return -1;
  }
#endif
  return 0;
}

int ttkEigenField::getConstraints(vtkPointSet *input) {
  auto osfn = static_cast<const char *>(ttk::OffsetScalarFieldName);
  auto pointData = input->GetPointData();

  if(InputScalarFieldName.length() != 0) {
    constraints_ = pointData->GetArray(InputScalarFieldName.data());
  } else if(pointData->GetArray(osfn) != nullptr) {
    constraints_ = pointData->GetArray(osfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(constraints_ == nullptr) {
    cerr << "[ttkEigenField] Error: wrong constraint "
            "scalar field."
         << endl;
    return -1;
  }
#endif

  if(constraints_->IsA("vtkDoubleArray")) {
    OutputScalarFieldType = EigenFieldType::Double;
  } else if(constraints_->IsA("vtkFloatArray")) {
    OutputScalarFieldType = EigenFieldType::Float;
  } else {
    return -2;
  }

  return 0;
}

int ttkEigenField::doIt(std::vector<vtkDataSet *> &inputs,
                           std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkPointSet *identifiers = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  // set this early, since it should trigger some triangulation pre-processing
  eigenField_.setUseCotanWeights(UseCotanWeights);

  res += getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkEigenField] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  res += getIdentifiers(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkEigenField] Error: wrong identifiers." << endl;
    return -2;
  }
#endif

  res += getConstraints(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkEigenField] Error: wrong constraints." << endl;
    return -2;
  }
#endif

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInDomain == 0) {
    cerr << "[ttkEigenField] Error: domain has no points." << endl;
    return -3;
  }
#endif

  const ttk::SimplexId numberOfPointsInSources
    = identifiers->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInSources == 0) {
    cerr << "[ttkEigenField] Error: sources have no points." << endl;
    return -4;
  }
#endif

  eigenField_.setVertexNumber(numberOfPointsInDomain);
  eigenField_.setConstraintNumber(numberOfPointsInSources);
  eigenField_.setSources(identifiers_->GetVoidPointer(0));
  eigenField_.setConstraints(constraints_->GetVoidPointer(0));
  eigenField_.setSolvingMethod(SolvingMethod);
  eigenField_.setLogAlpha(LogAlpha);

  vtkSmartPointer<vtkDataArray> harmonicScalarField{};

  switch(OutputScalarFieldType) {
    case EigenFieldType::Float:
      harmonicScalarField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case EigenFieldType::Double:
      harmonicScalarField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkEigenField] Error: Unknown scalar field type" << endl;
      return -7;
#endif // TTK_ENABLE_KAMIKAZE
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(harmonicScalarField == nullptr) {
    cerr << "[ttkEigenField] Error: vtkArray allocation problem." << endl;
    return -8;
  }
#endif

  harmonicScalarField->SetNumberOfComponents(1);
  harmonicScalarField->SetNumberOfTuples(numberOfPointsInDomain);
  harmonicScalarField->SetName(OutputScalarFieldName.data());

  eigenField_.setOutputScalarFieldPointer(
    harmonicScalarField->GetVoidPointer(0));

  switch(OutputScalarFieldType) {
    case EigenFieldType::Float:
      res += eigenField_.execute<float>();
      break;
    case EigenFieldType::Double:
      res += eigenField_.execute<double>();
      break;
    default:
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/EigenField
  if(res != 0) {
    cerr << "[ttkEigenField] EigenField.execute() "
            "error code: "
         << res << endl;
    return -9;
  }
#endif

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(harmonicScalarField);

  std::stringstream msg;
  msg << "[ttkEigenField] Memory usage: " << m.getElapsedUsage() << " MB."
      << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
