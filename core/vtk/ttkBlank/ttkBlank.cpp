#include <ttkBlank.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;
using namespace blank;

vtkStandardNewMacro(ttkBlank);
int ttkBlank::doIt(vector<vtkDataSet *> &inputs,
                   vector<vtkDataSet *> &outputs) {

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

  if(!triangulation)
    return -1;

  triangulation->setWrapper(this);
  blank_.setupTriangulation(triangulation);
  blank_.setWrapper(this);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = nullptr;

  if(!ScalarField.empty()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(0);
  }

  if(!inputScalarField)
    return -2;

  // allocate the memory for the output scalar field
  if(!outputScalarField_) {
    switch(inputScalarField->GetDataType()) {
      case VTK_CHAR:
      case VTK_DOUBLE:
      case VTK_FLOAT:
      case VTK_ID_TYPE:
      case VTK_INT:
      case VTK_SHORT:
      case VTK_UNSIGNED_CHAR:
      case VTK_UNSIGNED_SHORT:
        outputScalarField_ = inputScalarField->NewInstance();
        break;
      default:
        stringstream msg;
        msg << "[ttkBlank] Unsupported data type :(" << endl;
        dMsg(cerr, msg.str(), fatalMsg);
        return -3;
    }
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());

  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()) {
    output->GetPointData()->RemoveArray(ScalarField.data());
  } else {
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);

  // calling the executing package
  blank_.setInputDataPointer(ttkUtils::GetVoidPointer(inputScalarField));
  blank_.setOutputDataPointer(ttkUtils::GetVoidPointer(outputScalarField_));

  // scalar type
  // ttk::Triangulation::Type triangulationType = triangulation->getType();
  // int dataType = inputScalarField->GetDataType();

  ttkVtkTemplateMacro(
    inputScalarField->GetDataType(), triangulation->getType(),
    (blank_.execute<VTK_TT, TTK_TT>(
      (TTK_TT *)triangulation->getData(), SomeIntegerArgument)));

  {
    stringstream msg;
    msg << "[ttkBlank] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
