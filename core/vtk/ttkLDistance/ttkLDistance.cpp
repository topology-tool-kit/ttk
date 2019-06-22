#include "ttkLDistance.h"

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkLDistance)

  int ttkLDistance::doIt(vector<vtkDataSet *> &inputs,
                         vector<vtkDataSet *> &outputs) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr << "[ttkLDistance] Error: not enough input information." << endl;
    return -1;
  }
#endif

  vtkDataSet *input1 = inputs[0];
  vtkDataSet *output = outputs[0];

  // Test validity of datasets (must present the same number of points).
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input1) {
    cerr << "[ttkLDistance] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!input1->GetNumberOfPoints()) {
    cerr << "[ttkLDistance] Error: input has no point." << endl;
    return -1;
  }

  if(!input1->GetPointData()) {
    cerr << "[ttkLDistance] Error: input has no point data." << endl;
    return -1;
  }

  if(!output) {
    cerr << "[ttkLDistance] Error: output pointer is NULL." << endl;
    return -1;
  }
#endif

  lDistance_.setWrapper(this);

  // Use a pointer-base copy for the input.
  output->ShallowCopy(input1);

  // In the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField1 = NULL;
  vtkDataArray *inputScalarField2 = NULL;

  if(ScalarField1.length())
    inputScalarField1 = input1->GetPointData()->GetArray(ScalarField1.data());
  else
    inputScalarField1 = input1->GetPointData()->GetArray(ScalarFieldId1);

  if(ScalarField2.length())
    inputScalarField2 = input1->GetPointData()->GetArray(ScalarField2.data());
  else
    inputScalarField2 = input1->GetPointData()->GetArray(ScalarFieldId2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField1 || !inputScalarField2
     || inputScalarField1->GetDataType() != inputScalarField2->GetDataType())
    cerr << "[ttkLDistance] Error: input scalar fields are NULL or have "
            "different types."
         << endl;
  return -1;
#endif

  // Allocate memory for the output scalar field, based on the first input.
  if(!outputScalarField_) {
    switch(inputScalarField1->GetDataType()) {
      case VTK_CHAR:
        outputScalarField_ = vtkSmartPointer<vtkCharArray>::New();
        break;
      case VTK_DOUBLE:
        outputScalarField_ = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case VTK_FLOAT:
        outputScalarField_ = vtkSmartPointer<vtkFloatArray>::New();
        break;
      case VTK_INT:
        outputScalarField_ = vtkSmartPointer<vtkIntArray>::New();
        break;
      case VTK_ID_TYPE:
        outputScalarField_ = vtkSmartPointer<vtkIdTypeArray>::New();
        break;
      default: {
        stringstream msg;
        msg << "[vtkLDistance] Unsupported data type :(" << endl;
        dMsg(cerr, msg.str(), fatalMsg);
      } break;
    }
  }

  const char *fieldName = DistanceFieldName.c_str();

  SimplexId numberOfPoints = (SimplexId)input1->GetNumberOfPoints();

  lDistance_.setNumberOfPoints(numberOfPoints);
  outputScalarField_->SetNumberOfTuples(numberOfPoints);
  outputScalarField_->SetName(fieldName);

  // On the output, replace the field array by a pointer to its processed
  // version

  output->GetPointData()->AddArray(outputScalarField_);

  // Calling the executing package.
  switch(inputScalarField1->GetDataType()) {
    vtkTemplateMacro({
      lDistance_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
      lDistance_.setInputDataPointer1(inputScalarField1->GetVoidPointer(0));
      lDistance_.setInputDataPointer2(inputScalarField2->GetVoidPointer(0));
      lDistance_.execute<VTK_TT>(DistanceType);
    });
  }

  result = lDistance_.getResult();

  return 0;
}
