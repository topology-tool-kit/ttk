#include <ttkHarmonicFieldComputation.h>

vtkStandardNewMacro(ttkHarmonicFieldComputation);

ttkHarmonicFieldComputation::ttkHarmonicFieldComputation()
    : triangulation_{}, identifiers_{}, constraints_{} {

  SetNumberOfInputPorts(2);

  InputScalarFieldName =
      std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
  OutputScalarFieldName = "OutputHarmonicField";
  OutputScalarFieldType = HarmonicFieldType::Float;

  ForceInputScalarField = false;
  ThreadNumber = threadNumber_;
  UseAllCores = true;
  InputIdentifiersFieldName =
      std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
}

int ttkHarmonicFieldComputation::FillInputPortInformation(
    int port, vtkInformation *info) {

  if (port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  if (port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  }
  return 1;
}

int ttkHarmonicFieldComputation::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if (triangulation_ == nullptr) {
    cerr << "[ttkHarmonicFieldComputation] Error: input triangulation pointer "
            "is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation_->setWrapper(this);
  harmonicField_.setWrapper(this);
  harmonicField_.setupTriangulation(triangulation_);
  Modified();

#ifndef TTK_ENABLE_KAMIKAZE
  if (triangulation_->isEmpty()) {
    cerr << "[ttkHarmonicFieldComputation] Error: ttkTriangulation allocation "
            "problem."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::getIdentifiers(vtkPointSet *input) {
  auto vsfn = static_cast<const char *>(ttk::VertexScalarFieldName);

  if (ForceInputScalarField && InputIdentifiersFieldName.length() != 0) {
    identifiers_ =
        input->GetPointData()->GetArray(InputIdentifiersFieldName.data());
  } else if (input->GetPointData()->GetArray(vsfn) != nullptr) {
    identifiers_ = input->GetPointData()->GetArray(vsfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if (identifiers_ == nullptr) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong vertex identifier "
            "scalar field."
         << endl;
    return -1;
  }
#endif
  return 0;
}

int ttkHarmonicFieldComputation::getConstraints(vtkPointSet *input) {
  auto osfn = static_cast<const char *>(ttk::OffsetScalarFieldName);

  if (InputScalarFieldName.length() != 0) {
    constraints_ = input->GetPointData()->GetArray(InputScalarFieldName.data());
  } else if (input->GetPointData()->GetArray(osfn) != nullptr) {
    constraints_ = input->GetPointData()->GetArray(osfn);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if (constraints_ == nullptr) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong constraint "
            "scalar field."
         << endl;
    return -1;
  }
#endif

  return 0;
}

int ttkHarmonicFieldComputation::doIt(std::vector<vtkDataSet *> &inputs,
                                      std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkPointSet *identifiers = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if (res != 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  res += getIdentifiers(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if (res != 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong identifiers." << endl;
    return -2;
  }
#endif

  res += getConstraints(identifiers);

#ifndef TTK_ENABLE_KAMIKAZE
  if (res != 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong constraints." << endl;
    return -2;
  }
#endif

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if (numberOfPointsInDomain == 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: domain has no points."
         << endl;
    return -3;
  }
#endif

  const ttk::SimplexId numberOfPointsInSources =
      identifiers->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if (numberOfPointsInSources == 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: sources have no points."
         << endl;
    return -4;
  }
#endif

  harmonicField_.setVertexNumber(numberOfPointsInDomain);
  harmonicField_.setConstraintNumber(numberOfPointsInSources);
  harmonicField_.setSources(identifiers_->GetVoidPointer(0));
  harmonicField_.setConstraints(constraints_->GetVoidPointer(0));
  harmonicField_.setUseCotanWeights(UseCotanWeights);

  vtkDataArray *harmonicScalarField{};

  switch (OutputScalarFieldType) {
  case HarmonicFieldType::Float:
    harmonicScalarField = vtkFloatArray::New();
    break;
  case HarmonicFieldType::Double:
    harmonicScalarField = vtkDoubleArray::New();
    break;
  default:
#ifndef TTK_ENABLE_KAMIKAZE
    cerr << "[ttkHarmonicFieldComputation] Error: Unknown scalar field type"
         << endl;
    return -7;
#endif // TTK_ENABLE_KAMIKAZE
    break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if (harmonicScalarField == nullptr) {
    cerr << "[ttkHarmonicFieldComputation] Error: vtkArray allocation problem."
         << endl;
    return -8;
  }
#endif

  harmonicScalarField->SetNumberOfComponents(1);
  harmonicScalarField->SetNumberOfTuples(numberOfPointsInDomain);
  harmonicScalarField->SetName(OutputScalarFieldName.data());

  harmonicField_.setOutputScalarFieldPointer(
      harmonicScalarField->GetVoidPointer(0));

  switch (OutputScalarFieldType) {
  case HarmonicFieldType::Float:
    res += harmonicField_.execute<float>();
    break;
  case HarmonicFieldType::Double:
    res += harmonicField_.execute<double>();
    break;
  default:
    break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/HarmonicFieldComputation
  if (res != 0) {
    cerr << "[ttkHarmonicFieldComputation] HarmonicFieldComputation.execute() "
            "error code: "
         << res << endl;
    return -9;
  }
#endif

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(harmonicScalarField);
  harmonicScalarField->Delete();

  std::stringstream msg;
  msg << "[ttkHarmonicFieldComputation] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
