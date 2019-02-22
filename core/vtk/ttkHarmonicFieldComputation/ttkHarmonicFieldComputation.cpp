#include <ttkHarmonicFieldComputation.h>

vtkStandardNewMacro(ttkHarmonicFieldComputation);

ttkHarmonicFieldComputation::ttkHarmonicFieldComputation()
    : triangulation_{}, identifiers_{} {

  SetNumberOfInputPorts(3);

  InputScalarFieldName =
      std::string(static_cast<const char *>(ttk::VertexScalarFieldName));
  OutputScalarFieldName = "OutputHarmonicField";
  OutputScalarFieldType = HarmonicFieldType::Float;

  ThreadNumber = threadNumber_;
  UseAllCores = true;
}

int ttkHarmonicFieldComputation::FillInputPortInformation(
    int port, vtkInformation *info) {

  if (port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  if (port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  }
  if (port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
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

  if (InputScalarFieldName.length() != 0) {
    identifiers_ = input->GetPointData()->GetArray(InputScalarFieldName.data());
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

int ttkHarmonicFieldComputation::doIt(std::vector<vtkDataSet *> &inputs,
                                      std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkPointSet *sources = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *constraints = vtkDataSet::SafeDownCast(inputs[2]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  res |= getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if (res != 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  res |= getIdentifiers(sources);

#ifndef TTK_ENABLE_KAMIKAZE
  if (res) {
    cerr << "[ttkHarmonicFieldComputation] Error: wrong identifiers." << endl;
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

  const ttk::SimplexId numberOfPointsInSources = sources->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if (numberOfPointsInSources == 0) {
    cerr << "[ttkHarmonicFieldComputation] Error: sources have no points."
         << endl;
    return -4;
  }
#endif

  auto origin = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  if (origin != nullptr) {
    origin->SetNumberOfComponents(1);
    origin->SetNumberOfTuples(numberOfPointsInDomain);
    origin->SetName(ttk::VertexScalarFieldName);
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkHarmonicFieldComputation] Error: ttkSimplexIdTypeArray "
            "allocation problem."
         << endl;
    return -5;
  }
#endif

  auto seg = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  if (seg != nullptr) {
    seg->SetNumberOfComponents(1);
    seg->SetNumberOfTuples(numberOfPointsInDomain);
    seg->SetName("SeedIdentifier");
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else {
    cerr << "[ttkHarmonicFieldComputation] Error: ttkSimplexIdTypeArray "
            "allocation problem."
         << endl;
    return -6;
  }
#endif

  harmonicField_.setVertexNumber(numberOfPointsInDomain);
  harmonicField_.setConstraintNumber(numberOfPointsInSources);
  harmonicField_.setConstraints(identifiers_->GetVoidPointer(0));
  harmonicField_.setOutputIdentifiers(origin->GetVoidPointer(0));
  harmonicField_.setOutputSegmentation(seg->GetVoidPointer(0));

  vtkDataArray *harmonicScalars{};

  switch (OutputScalarFieldType) {
  case HarmonicFieldType::Float:
    harmonicScalars = vtkFloatArray::New();
    break;
  case HarmonicFieldType::Double:
    harmonicScalars = vtkDoubleArray::New();
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
  if (harmonicScalars == nullptr) {
    cerr << "[ttkHarmonicFieldComputation] Error: vtkArray allocation problem."
         << endl;
    return -8;
  }
#endif

  harmonicScalars->SetNumberOfComponents(1);
  harmonicScalars->SetNumberOfTuples(numberOfPointsInDomain);
  harmonicScalars->SetName(OutputScalarFieldName.data());
  harmonicField_.setOutputScalarFieldPointer(
      harmonicScalars->GetVoidPointer(0));

  switch (OutputScalarFieldType) {
  case HarmonicFieldType::Float:
    res |= harmonicField_.execute<float>();
    break;
  case HarmonicFieldType::Double:
    res |= harmonicField_.execute<double>();
    break;
  default:
    break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  // something wrong in base/HarmonicFieldComputation
  if (res) {
    cerr << "[ttkHarmonicFieldComputation] HarmonicFieldComputation.execute() "
            "error code: "
         << res << endl;
    return -9;
  }
#endif

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(harmonicScalars);
  output->GetPointData()->AddArray(origin);
  output->GetPointData()->AddArray(seg);
  harmonicScalars->Delete();

  std::stringstream msg;
  msg << "[ttkHarmonicFieldComputation] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
