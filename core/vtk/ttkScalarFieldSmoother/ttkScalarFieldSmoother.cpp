#include                  <ttkScalarFieldSmoother.h>

vtkStandardNewMacro(ttkScalarFieldSmoother)

ttkScalarFieldSmoother::ttkScalarFieldSmoother(){

  // init
  NumberOfIterations = 1;
  ScalarFieldIdentifier = 0;
  MaskIdentifier = 0;
  UseInputMask = false;
  outputScalarField_ = NULL;
  
}

ttkScalarFieldSmoother::~ttkScalarFieldSmoother(){

  if(outputScalarField_)
    outputScalarField_->Delete();
}

int ttkScalarFieldSmoother::doIt(vector<vtkDataSet *> &inputs, 
  vector<vtkDataSet *> &outputs){

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
  
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  smoother_.setupTriangulation(triangulation);
  smoother_.setWrapper(this);

  // This filter copies the input into a new data-set (smoothed)
  // let's use shallow copies, in order to only duplicate point positions 
  // (before and after). the rest is not changed, pointers are sufficient.
  output->ShallowCopy(input);
 
  vtkDataArray *inputScalarField = NULL;
  vtkCharArray *inputMaskField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(ScalarFieldIdentifier);
  }
  
  if(!inputScalarField)
    return -2;

  if(inputScalarField->GetName())
    ScalarField = inputScalarField->GetName();
  
  {
    stringstream msg;
    msg << "[ScalarFieldSmoother] Using field `"
      << inputScalarField->GetName() << "'..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  if(UseInputMask){
      if(InputMask.length()){
          inputMaskField = vtkCharArray::SafeDownCast(input->GetPointData()->GetArray(InputMask.data()));
      } else {
          inputMaskField = vtkCharArray::SafeDownCast(input->GetPointData()->GetArray(MaskIdentifier));
      }

      if(inputMaskField->GetName())
          InputMask = inputMaskField->GetName();

      {
          stringstream msg;
          msg << "[ScalarFieldSmoother] Using mask `"
              << inputMaskField->GetName() << "'..." << endl;
          dMsg(cout, msg.str(), infoMsg);
      }
  }

  if(outputScalarField_){
    outputScalarField_->Delete();
    outputScalarField_ = NULL;
  }
  
  // allocate the memory for the output scalar field
  switch(inputScalarField->GetDataType()){
    
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
    
    case VTK_UNSIGNED_SHORT:
      outputScalarField_ = vtkUnsignedShortArray::New();
      break;
      
    default:
      {
        stringstream msg;
        msg << "[ttkScalarFieldSmoother] Unsupported data type :(" << endl;
        dMsg(cerr, msg.str(), fatalMsg);
      }
      break;
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());
  
  // on the output, replace the field array by a pointer to its smoothed version
  if(ScalarField.length()){
    output->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);

  void* inputMaskPtr = (inputMaskField)?inputMaskField->GetVoidPointer(0):nullptr;
  
  // calling the smoothing package
  switch(inputScalarField->GetDataType()){
    
#ifndef _MSC_VER
	  vtkTemplateMacro((
	  {
		  smoother_.setDimensionNumber(inputScalarField->GetNumberOfComponents());
	  smoother_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
	  smoother_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
	  smoother_.setMaskDataPointer(inputMaskPtr);
	  smoother_.smooth<VTK_TT>(NumberOfIterations);
	  }
	  ));
#else
	  vtkTemplateMacro(
	  {
		  smoother_.setDimensionNumber(inputScalarField->GetNumberOfComponents());
	  smoother_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
	  smoother_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
	  smoother_.setMaskDataPointer(inputMaskPtr);
	  smoother_.smooth<VTK_TT>(NumberOfIterations);
	  }
	  );
#endif
  }
 
  {
    stringstream msg;
    msg << "[ttkScalarFieldSmoother] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
 
  return 0;
}
