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

  vtkSmartPointer<vtkDataArray> inputScalarField{};

  if(ScalarField.length() >= 0) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(0);
  }

  if(!inputScalarField) {
    return -2;
  }

  // allocate the memory for the output scalar field
  if(!outputScalarField_) {
    switch(inputScalarField->GetDataType()) {
      case VTK_CHAR:
        outputScalarField_ = vtkCharArray::New();
        break;
      case VTK_DOUBLE:
        outputScalarField_ = vtkDoubleArray::New();
        break;
      case VTK_FLOAT:
        outputScalarField_ = vtkFloatArray::New();
        break;
      case VTK_INT:
        outputScalarField_ = vtkIntArray::New();
        break;
      case VTK_ID_TYPE:
        outputScalarField_ = vtkIdTypeArray::New();
        break;
      default:
        std::stringstream msg;
        msg << MODULE_S "Unsupported data type :(" << std::endl;
        dMsg(std::cerr, msg.str(), fatalMsg);
        return -3;
    }
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());

  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()>=0) {
    output->GetPointData()->RemoveArray(ScalarField.data());
  } else {
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);

  // calling the executing package
  switch(inputScalarField->GetDataType()) {
    ttkTemplateMacro(baseWorker_.execute<VTK_TT>());
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(std::cout, msg.str(), memoryMsg);
  }

  return 0;
}
