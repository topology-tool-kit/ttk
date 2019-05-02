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

  const auto vertexNumber = domain->GetNumberOfPoints();

  TTK_ABORT_KK(vertexNumber == 0, "domain has no points", -2);

  baseWorker_.setEigenNumber(EigenNumber);

  // array of eigenfunctions
  vtkSmartPointer<vtkDataArray> eigenFunctions{};
  // statistics
  vtkSmartPointer<vtkDataArray> stats{};

  switch(OutputFieldType) {
    case EigenFieldType::Float:
      eigenFunctions = vtkSmartPointer<vtkFloatArray>::New();
      stats = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case EigenFieldType::Double:
      eigenFunctions = vtkSmartPointer<vtkDoubleArray>::New();
      stats = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    default:
      TTK_ABORT_KK(true, "unknown field type", -7);
      break;
  }

  TTK_ABORT_KK(eigenFunctions == nullptr, "vtkArray allocation problem", -3);

  eigenFunctions->SetNumberOfComponents(EigenNumber);
  eigenFunctions->SetNumberOfTuples(vertexNumber);
  eigenFunctions->SetName(OutputFieldName.data());

  stats->SetName("Statistics");
  const int statsComp = 3;
  stats->SetNumberOfComponents(statsComp);
  stats->SetNumberOfTuples(vertexNumber);
  stats->SetComponentName(0, "Min");
  stats->SetComponentName(1, "Max");
  stats->SetComponentName(2, "Sum");

  baseWorker_.setOutputFieldPointer(eigenFunctions->GetVoidPointer(0));

  switch(OutputFieldType) {
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

  // compute statistics
  switch(OutputFieldType) {
    case EigenFieldType::Float: {
      auto outp = static_cast<float *>(stats->GetVoidPointer(0));
      for(size_t i = 0; i < vertexNumber; ++i) {
        // outp index
        auto k = statsComp * i;
        // GetTuple returns a double* pointer
        auto inp = eigenFunctions->GetTuple(i);
        outp[k] = float(inp[0]);
        outp[k + 1] = outp[k];
        outp[k + 2] = outp[k];
        for(size_t j = 1; j < EigenNumber; ++j) {
          outp[k] = std::min<float>(outp[k], float(inp[j]));
          outp[k + 1] = std::max<float>(outp[k + 1], float(inp[j]));
          outp[k + 2] += float(inp[j]);
        }
      }
    } break;
    case EigenFieldType::Double: {
      auto outp = static_cast<double *>(stats->GetVoidPointer(0));
      for(size_t i = 0; i < statsComp * vertexNumber; i += statsComp) {
        // outp index
        auto k = statsComp * i;
        // GetTuple returns a double* pointer
        auto inp = eigenFunctions->GetTuple(i);
        outp[k] = inp[0];
        outp[k + 1] = outp[k];
        outp[k + 2] = outp[k];
        for(size_t j = 1; j < EigenNumber; ++j) {
          outp[k] = std::min<double>(outp[k], inp[j]);
          outp[k + 1] = std::max<double>(outp[k + 1], inp[j]);
          outp[k + 2] += inp[j];
        }
      }
    } break;
    default:
      break;
  }

  // update result
  output->ShallowCopy(domain);
  output->GetPointData()->AddArray(eigenFunctions);
  output->GetPointData()->AddArray(stats);

  std::stringstream msg;
  msg << "[ttkEigenField] Memory usage: " << m.getElapsedUsage() << " MB."
      << endl;
  dMsg(cout, msg.str(), memoryMsg);

  return 0;
}
