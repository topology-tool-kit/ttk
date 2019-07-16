#include <ttkScalarFieldSmoother.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldSmoother)

  ttkScalarFieldSmoother::ttkScalarFieldSmoother() {

  // init
  NumberOfIterations = 1;
  ScalarFieldIdentifier = 0;
  MaskIdentifier = 0;
  ForceInputMaskScalarField = false;
  InputMask = ttk::MaskScalarFieldName;
  outputScalarField_ = NULL;

  UseAllCores = true;
}

ttkScalarFieldSmoother::~ttkScalarFieldSmoother() {

  if(outputScalarField_)
    outputScalarField_->Delete();
}

int ttkScalarFieldSmoother::doIt(vector<vtkDataSet *> &inputs,
                                 vector<vtkDataSet *> &outputs) {

  Memory m;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr << "[ttkScalarFieldSmoother] Error: not enough input information."
         << endl;
    return -1;
  }
#endif

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkScalarFieldSmoother] Error: input pointer is NULL." << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkScalarFieldSmoother] Error: input has no point." << endl;
    return -1;
  }

  if(!input->GetPointData()) {
    cerr << "[ttkScalarFieldSmoother] Error: input has no point data." << endl;
    return -1;
  }

  if(!output) {
    cerr << "[ttkScalarFieldSmoother] Error: output pointer is NULL." << endl;
    return -1;
  }
#endif

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    cerr << "[ttkScalarFieldSmoother] Error: input triangulation is NULL."
         << endl;
    return -1;
  }
#endif

  triangulation->setWrapper(this);
  smoother_.setupTriangulation(triangulation);
  smoother_.setWrapper(this);

  // This filter copies the input into a new data-set (smoothed)
  // let's use shallow copies, in order to only duplicate point positions
  // (before and after). the rest is not changed, pointers are sufficient.
  output->ShallowCopy(input);

  vtkDataArray *inputScalarField = NULL;
  vtkCharArray *inputMaskField = NULL;

  if(ScalarField.length()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(ScalarFieldIdentifier);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField) {
    cerr << "[ttkScalarFieldSmoother] Error: input scalar field poiner is NULL."
         << endl;
    return -2;
  }
#endif

  if(inputScalarField->GetName())
    ScalarField = inputScalarField->GetName();

  {
    stringstream msg;
    msg << "[ScalarFieldSmoother] Using field `" << inputScalarField->GetName()
        << "'..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  if(ForceInputMaskScalarField) {
    if(InputMask.length()) {
      inputMaskField = vtkCharArray::SafeDownCast(
        input->GetPointData()->GetArray(InputMask.data()));
    } else {
      inputMaskField = vtkCharArray::SafeDownCast(
        input->GetPointData()->GetArray(MaskIdentifier));
    }

    if(inputMaskField->GetName())
      InputMask = inputMaskField->GetName();

    {
      stringstream msg;
      msg << "[ScalarFieldSmoother] Using mask `" << inputMaskField->GetName()
          << "'..." << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  } else if(input->GetPointData()->GetArray(ttk::MaskScalarFieldName)) {
    inputMaskField = vtkCharArray::SafeDownCast(
      input->GetPointData()->GetArray(ttk::MaskScalarFieldName));
    InputMask = ttk::MaskScalarFieldName;
  }

  if(outputScalarField_) {
    outputScalarField_->Delete();
    outputScalarField_ = NULL;
  }

  // allocate the memory for the output scalar field
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

    case VTK_SHORT:
      outputScalarField_ = vtkShortArray::New();
      break;

    case VTK_UNSIGNED_CHAR:
      outputScalarField_ = vtkUnsignedCharArray::New();
      break;

    case VTK_UNSIGNED_SHORT:
      outputScalarField_ = vtkUnsignedShortArray::New();
      break;

    default: {
      stringstream msg;
      msg << "[ttkScalarFieldSmoother] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    } break;
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());

  // on the output, replace the field array by a pointer to its smoothed version
  if(ScalarField.length()) {
    output->GetPointData()->RemoveArray(ScalarField.data());
  } else {
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);

  void *inputMaskPtr
    = (inputMaskField) ? inputMaskField->GetVoidPointer(0) : nullptr;

  // calling the smoothing package
  switch(inputScalarField->GetDataType()) {

    vtkTemplateMacro({
      smoother_.setDimensionNumber(inputScalarField->GetNumberOfComponents());
      smoother_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
      smoother_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
      smoother_.setMaskDataPointer(inputMaskPtr);
      smoother_.smooth<VTK_TT>(NumberOfIterations);
    });
  }

  {
    stringstream msg;
    msg << "[ttkScalarFieldSmoother] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
