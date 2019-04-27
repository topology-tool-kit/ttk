#include <ttkEigenField.h>

#define MODULE_S "[ttkEigenField] "
#define MODULE_ERROR_S MODULE_S "Error: "
#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)    \
  if(COND) {                            \
    cerr << MODULE_ERROR_S MSG << endl; \
    return RET;                         \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro(ttkEigenField);

ttkEigenField::ttkEigenField() {
  SetNumberOfInputPorts(1);
}

int ttkEigenField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  }
  return 0;
}

int ttkEigenField::getTriangulation(vtkDataSet *input) {

  triangulation_ = ttkTriangulation::getTriangulation(input);

  TTK_ABORT_KK(
    triangulation_ == nullptr, "input triangulation pointer is NULL.", -1);

  triangulation_->setWrapper(this);
  baseWorker_.setWrapper(this);
  baseWorker_.setupTriangulation(triangulation_);
  Modified();

  TTK_ABORT_KK(
    triangulation_->isEmpty(), "ttkTriangulation allocation problem.", -2);

  return 0;
}

int ttkEigenField::doIt(std::vector<vtkDataSet *> &inputs,
                        std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *domain = vtkDataSet::SafeDownCast(inputs[0]);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outputs[0]);

  int res = 0;

  res += getTriangulation(domain);

  TTK_ABORT_KK(res != 0, "wrong triangulation", -1);

  const auto numberOfPointsInDomain = domain->GetNumberOfPoints();

  TTK_ABORT_KK(numberOfPointsInDomain == 0, "domain has no points", -2);

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
      TTK_ABORT_KK(true, "unknown scalar field type", -7);
      break;
  }

  TTK_ABORT_KK(eigenScalarField == nullptr, "vtkArray allocation problem", -3);

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

  TTK_ABORT_KK(res != 0, "EigenField execute error code " << res, -4);

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(eigenScalarField);

  std::stringstream msg;
  msg << "[ttkEigenField] Memory usage: " << m.getElapsedUsage() << " MB."
      << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
