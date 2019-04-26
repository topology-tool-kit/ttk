#include <ttkEigenField.h>

vtkStandardNewMacro(ttkEigenField);

ttkEigenField::ttkEigenField() {
  SetNumberOfInputPorts(2);
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
  baseWorker_.setWrapper(this);
  baseWorker_.setupTriangulation(triangulation_);
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

int ttkEigenField::doIt(std::vector<vtkDataSet *> &inputs,
                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkPointSet *identifiers = vtkPointSet::SafeDownCast(inputs[1]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(domain);

#ifndef TTK_ENABLE_KAMIKAZE
  if(res != 0) {
    cerr << "[ttkEigenField] Error: wrong triangulation." << endl;
    return -1;
  }
#endif

  const ttk::SimplexId numberOfPointsInDomain = domain->GetNumberOfPoints();

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPointsInDomain == 0) {
    cerr << "[ttkEigenField] Error: domain has no points." << endl;
    return -3;
  }
#endif

  baseWorker_.setEigenNumber(EigenNumber);

  vtkSmartPointer<vtkDataArray> eigenScalarField{};

  switch(OutputScalarFieldType) {
    case EigenFieldType::Float:
      eigenScalarField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case EigenFieldType::Double:
      eigenScalarField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkEigenField] Error: Unknown scalar field type" << endl;
      return -7;
#endif // TTK_ENABLE_KAMIKAZE
      break;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(eigenScalarField == nullptr) {
    cerr << "[ttkEigenField] Error: vtkArray allocation problem." << endl;
    return -8;
  }
#endif

  eigenScalarField->SetNumberOfComponents(1);
  eigenScalarField->SetNumberOfTuples(numberOfPointsInDomain);
  eigenScalarField->SetName(OutputScalarFieldName.data());

  baseWorker_.setOutputScalarFieldPointer(eigenScalarField->GetVoidPointer(0));

  switch(OutputScalarFieldType) {
    case EigenFieldType::Float:
      res += baseWorker_.execute<float>();
      break;
    case EigenFieldType::Double:
      res += baseWorker_.execute<double>();
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
  output->GetPointData()->AddArray(eigenScalarField);

  std::stringstream msg;
  msg << "[ttkEigenField] Memory usage: " << m.getElapsedUsage() << " MB."
      << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
