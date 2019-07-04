#include <ttkBarycentricSubdivision.h>

#define MODULE_S "[ttkBarycentricSubdivision] "
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

vtkStandardNewMacro(ttkBarycentricSubdivision);

int ttkBarycentricSubdivision::doIt(std::vector<vtkDataSet *> &inputs,
                                    std::vector<vtkDataSet *> &outputs) {

  ttk::Memory m;

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

  ttk::Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(triangulation == nullptr) {
    return -1;
  }

  triangulation->setWrapper(this);
  baseWorker_.setupTriangulation(triangulation);
  baseWorker_.setWrapper(this);

  vtkSmartPointer<vtkDataArray> inputScalarField
    = input->GetPointData()->GetArray(0);

  if(inputScalarField == nullptr) {
    return -2;
  }

  vtkSmartPointer<vtkDataArray> outputScalarField{};

  // allocate the memory for the output scalar field
  switch(inputScalarField->GetDataType()) {
    case VTK_CHAR:
      outputScalarField = vtkCharArray::New();
      break;
    case VTK_DOUBLE:
      outputScalarField = vtkDoubleArray::New();
      break;
    case VTK_FLOAT:
      outputScalarField = vtkFloatArray::New();
      break;
    case VTK_INT:
      outputScalarField = vtkIntArray::New();
      break;
    case VTK_ID_TYPE:
      outputScalarField = vtkIdTypeArray::New();
      break;
    default:
      std::stringstream msg;
      msg << MODULE_S "Unsupported data type :(" << std::endl;
      dMsg(std::cerr, msg.str(), fatalMsg);
      return -3;
  }

  outputScalarField->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField->SetName(inputScalarField->GetName());
  output->GetPointData()->AddArray(outputScalarField);

  // calling the executing package
  switch(inputScalarField->GetDataType()) {
    ttkTemplateMacro(baseWorker_.execute());
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(std::cout, msg.str(), memoryMsg);
  }

  return 0;
}
