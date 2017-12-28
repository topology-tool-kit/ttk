#include                  <ttkUncertainDataEstimator.h>

vtkStandardNewMacro(ttkUncertainDataEstimator)

ttkUncertainDataEstimator::ttkUncertainDataEstimator(){

  // init
  outputLowerBoundScalarField_ = NULL;
  outputUpperBoundScalarField_ = NULL;
  outputProbabilityScalarField_ = NULL;
  outputMeanField_ = NULL;

  UseAllCores = false;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

ttkUncertainDataEstimator::~ttkUncertainDataEstimator(){

  if(outputLowerBoundScalarField_)
    outputLowerBoundScalarField_->Delete();
  if(outputUpperBoundScalarField_)
    outputUpperBoundScalarField_->Delete();
  if(outputProbabilityScalarField_){
    for(int b=0 ; b<allocatedBinCount_ ; b++){
      outputProbabilityScalarField_[b]->Delete();
    }
    free(outputProbabilityScalarField_);
  }
  if(outputMeanField_) {
    outputMeanField_->Delete();
  }
}


// transmit abort signals -- to copy paste in other wrappers
bool ttkUncertainDataEstimator::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkUncertainDataEstimator::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}



int ttkUncertainDataEstimator::doIt(vtkDataSet **input, vtkDataSet *outputBoundFields, vtkDataSet *outputProbability, vtkDataSet *outputMean, int numInputs){


  // Use a pointer-base copy for the input data
  outputBoundFields->ShallowCopy(input[0]);
  outputProbability->ShallowCopy(input[0]);
  outputMean->ShallowCopy(input[0]);

  // Get arrays from input datas
  //vtkDataArray* inputScalarField[numInputs] = { NULL };
  vector<vtkDataArray*> inputScalarField(numInputs);
  for(int i=0 ; i<numInputs ; i++){
    if(ScalarField.length()){
      inputScalarField[i] = input[i]->GetPointData()->GetArray(ScalarField.data());
    }
    else{
      inputScalarField[i] = input[i]->GetPointData()->GetArray(0);
    }
    // Check if inputs have the same data type and the same number of points
    if(inputScalarField[i]->GetDataType() != inputScalarField[0]->GetDataType()){
      stringstream msg;
      msg << "[ttkUncertainDataEstimator] Inputs of different data types." << endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return -3;
      }
    if(inputScalarField[i]->GetNumberOfTuples() != inputScalarField[0]->GetNumberOfTuples()){
      stringstream msg;
      msg << "[ttkUncertainDataEstimator] Inputs with different number of points." << endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return -2;
      }
    // Check if all the inputs are here
    if(!inputScalarField[i])
      return -1;
  }

  // Allocate the memory for the output bound scalar fields
  if(!outputLowerBoundScalarField_
     && !outputUpperBoundScalarField_){
    switch(inputScalarField[0]->GetDataType()){

      case VTK_CHAR:
        outputLowerBoundScalarField_ = vtkCharArray::New();
        outputUpperBoundScalarField_ = vtkCharArray::New();
        break;

      case VTK_DOUBLE:
        outputLowerBoundScalarField_ = vtkDoubleArray::New();
        outputUpperBoundScalarField_ = vtkDoubleArray::New();
        break;

      case VTK_FLOAT:
        outputLowerBoundScalarField_ = vtkFloatArray::New();
        outputUpperBoundScalarField_ = vtkFloatArray::New();
        break;

      case VTK_INT:
        outputLowerBoundScalarField_ = vtkIntArray::New();
        outputUpperBoundScalarField_ = vtkIntArray::New();
        break;

      stringstream msg;
      msg << "[ttkUncertainDataEstimator] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
    outputLowerBoundScalarField_->SetName("lowerBoundField");
    outputUpperBoundScalarField_->SetName("upperBoundField");
  }



  // Allocate the memory for the output probability scalar fields
  // Delete existing arrays
  cout << "Bin Count = " << binCount_ << endl;
  if(outputProbabilityScalarField_){
    for(int b=0 ; b<allocatedBinCount_ ; b++){
      outputProbabilityScalarField_[b]->Delete();
    }
    free(outputProbabilityScalarField_);
  }
  // Allocate array of vtkDoubleArray*
  outputProbabilityScalarField_ = (vtkDoubleArray **) malloc(binCount_*sizeof(vtkDoubleArray *));
  allocatedBinCount_ = binCount_;
  // Delete the pointer to the input field
  if(ScalarField.length()){
    outputProbability->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    outputProbability->GetPointData()->RemoveArray(0);
  }
  // Create DoubleArray objects and link them to the data set
  for(int b=0 ; b<binCount_ ; b++)
  {
    outputProbabilityScalarField_[b] = vtkDoubleArray::New();
    outputProbabilityScalarField_[b]->SetNumberOfTuples(input[0]->GetNumberOfPoints());
    outputProbability->GetPointData()->AddArray(outputProbabilityScalarField_[b]);
    outputProbabilityScalarField_[b]->FillComponent(0,0.);
  }

  // Mean field data set
  // Remove Arrays
  int numberOfArrays = outputMean->GetPointData()->GetNumberOfArrays();
  for(int i=0 ; i<numberOfArrays ; i++) {
    outputMean->GetPointData()->RemoveArray(outputMean->GetPointData()->GetArrayName(i));
  }
  // Allocate new array
  if(!outputMeanField_) {
    outputMeanField_ = vtkDoubleArray::New();
    outputMeanField_->SetNumberOfTuples(input[0]->GetNumberOfPoints());
    outputMeanField_->SetName("meanField");
    outputMeanField_->FillComponent(0, 0.0);
  }
  outputMean->GetPointData()->AddArray(outputMeanField_);


  // On the output, replace the field array by a pointer to its processed version
  if(ScalarField.length()){
    outputBoundFields->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    outputBoundFields->GetPointData()->RemoveArray(0);
  }

  // Resize arrays and add them in the output if required
  if(computeLowerBound_){
    outputLowerBoundScalarField_->SetNumberOfTuples(input[0]->GetNumberOfPoints());
    outputBoundFields->GetPointData()->AddArray(outputLowerBoundScalarField_);
  }
  else{
    outputLowerBoundScalarField_->SetNumberOfTuples(0);
  }

  if(computeUpperBound_){
    outputUpperBoundScalarField_->SetNumberOfTuples(input[0]->GetNumberOfPoints());
    outputBoundFields->GetPointData()->AddArray(outputUpperBoundScalarField_);
  }
  else{
    outputUpperBoundScalarField_->SetNumberOfTuples(0);
  }


  // Calling the executing package
  switch(inputScalarField[0]->GetDataType()){

    vtkTemplateMacro(
    {
      UncertainDataEstimator uncertainDataEstimator;
      uncertainDataEstimator.setWrapper(this);
 
      uncertainDataEstimator.setVertexNumber(
        outputBoundFields->GetNumberOfPoints());

      uncertainDataEstimator.setNumberOfInputs(numInputs);
      for (int i = 0; i<numInputs; i++) {
        uncertainDataEstimator.setInputDataPointer(i, 
        inputScalarField[i]->GetVoidPointer(0));
      }

      uncertainDataEstimator.setComputeLowerBound(computeLowerBound_);
      uncertainDataEstimator.setComputeUpperBound(computeUpperBound_);

    
      uncertainDataEstimator.setOutputLowerBoundField(
        outputLowerBoundScalarField_->GetVoidPointer(0));
        
      uncertainDataEstimator.setOutputUpperBoundField(
        outputUpperBoundScalarField_->GetVoidPointer(0));

      uncertainDataEstimator.setOutputMeanField(
        outputMeanField_->GetVoidPointer(0));

      uncertainDataEstimator.setBinCount(binCount_);
      for (int b = 0; b<binCount_; b++) {
        uncertainDataEstimator.setOutputProbability(b, 
          outputProbabilityScalarField_[b]->GetVoidPointer(0));
      }

      uncertainDataEstimator.execute<VTK_TT>();

      for (int b = 0; b<binCount_; b++) {
        stringstream name;
        name << setprecision(8) << uncertainDataEstimator.getBinValue(b);
        outputProbabilityScalarField_[b]->SetName(name.str().c_str());
      }
    }
    );
  }

  return 0;
}



int ttkUncertainDataEstimator::FillInputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillInputPortInformation(port, info)){
    return 0;
  }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkUncertainDataEstimator::FillOutputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillOutputPortInformation(port, info)){
    return 0;
  }
  if(port==0 || port==1 || port==3)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}



// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkUncertainDataEstimator::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  // Output pointers informations
  vtkInformation* outInfo;
  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *boundFields = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Histograms
  outInfo = outputVector->GetInformationObject(1);
  vtkDataSet *probability = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Mean field
  outInfo = outputVector->GetInformationObject(2);
  vtkDataSet *mean = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] Number of inputs: " << numInputs << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Get input datas
  vtkDataSet* *input = new vtkDataSet*[numInputs];
  for(int i=0 ; i<numInputs ; i++)
  {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  doIt(input, boundFields, probability, mean, numInputs);

  delete input;

  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
