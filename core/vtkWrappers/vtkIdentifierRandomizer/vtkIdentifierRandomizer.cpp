#include                  <vtkIdentifierRandomizer.h>

vtkStandardNewMacro(vtkIdentifierRandomizer)

int vtkIdentifierRandomizer::doIt(vector<vtkDataSet *> &inputs, 
  vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  
  if(ScalarField.length()){
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  }
  else{
    inputScalarField = input->GetPointData()->GetArray(0);
  }
  
  if(!inputScalarField)
    return -2;
  
  // allocate the memory for the output scalar field
  if(!outputScalarField_){
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
        
      stringstream msg;
      msg << "[vtkIdentifierRandomizer] Unsupported data type :(" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
    }
  }
  outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
  outputScalarField_->SetName(inputScalarField->GetName());
  
  
  // on the output, replace the field array by a pointer to its processed
  // version
  if(ScalarField.length()){
    output->GetPointData()->RemoveArray(ScalarField.data());
  }
  else{
    output->GetPointData()->RemoveArray(0);
  }
  output->GetPointData()->AddArray(outputScalarField_);
  
  vector<pair<int, int> > identifierMap;
  
  for(int i = 0; i < inputScalarField->GetNumberOfTuples(); i++){
    double inputIdentifier = -1;
    inputScalarField->GetTuple(i, &inputIdentifier);
    
    bool isIn = false;
    for(int j = 0; j < (int) identifierMap.size(); j++){
      if(identifierMap[j].first == inputIdentifier){
        isIn = true;
        break;
      }
    }
    
    if(!isIn){
      identifierMap.push_back(pair<int, int>(inputIdentifier, -INT_MAX));
    }
  }
  
  // now let's shuffle things around
  int freeIdentifiers = identifierMap.size();
  int randomIdentifier = -1;
  
  for(int i = 0; i < (int) identifierMap.size(); i++){
    
    randomIdentifier = drand48()*(freeIdentifiers);
    
    int freeCounter = -1;
    for(int j = 0; j < (int) identifierMap.size(); j++){
      
      bool isFound = false;
      for(int k = 0; k < (int) identifierMap.size(); k++){
        
        if(identifierMap[k].second == identifierMap[j].first){
          isFound = true;
          break;
        }
      }
      if(!isFound){
        freeCounter++;
      }
        
      if(freeCounter >= randomIdentifier){
        randomIdentifier = identifierMap[j].first;
        break;
      }
    }
    
    identifierMap[i].second = randomIdentifier;
    freeIdentifiers--;
  }
  
  // now populate the output scalar field.
  for(int i = 0; i < inputScalarField->GetNumberOfTuples(); i++){
    double inputIdentifier = -1;
    double outputIdentifier = -1;
    
    inputScalarField->GetTuple(i, &inputIdentifier);
    
    for(int j = 0; j < (int) identifierMap.size(); j++){
      if(inputIdentifier == identifierMap[j].first){
        outputIdentifier = identifierMap[j].second;
        break;
      }
    }
    
    outputScalarField_->SetTuple(i, &outputIdentifier);
  }
  
  {
    stringstream msg;
    msg << "[vtkIdentifierRandomizer] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}