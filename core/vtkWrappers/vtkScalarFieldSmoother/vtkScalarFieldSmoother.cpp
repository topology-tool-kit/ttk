#include                  <vtkScalarFieldSmoother.h>

vtkStandardNewMacro(vtkScalarFieldSmoother)

vtkScalarFieldSmoother::vtkScalarFieldSmoother(){

  // init
  NumberOfIterations = 1;
  ScalarFieldIdentifier = 0;
  outputScalarField_ = NULL;
  
}

vtkScalarFieldSmoother::~vtkScalarFieldSmoother(){

  if(outputScalarField_)
    outputScalarField_->Delete();
}

int vtkScalarFieldSmoother::doIt(vector<vtkDataSet *> &inputs, 
  vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = vtkTriangulation::getTriangulation(input);
  
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
        msg << "[vtkScalarFieldSmoother] Unsupported data type :(" << endl;
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
  
  // calling the smoothing package
  switch(inputScalarField->GetDataType()){
    
    vtkTemplateMacro((
      {
        smoother_.setDimensionNumber(inputScalarField->GetNumberOfComponents());
        smoother_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
        smoother_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
        smoother_.smooth<VTK_TT>(NumberOfIterations);
      }      
    ));
  }
 
  {
    stringstream msg;
    msg << "[vtkScalarFieldSmoother] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
 
  return 0;
}
