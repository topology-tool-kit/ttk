#include                  <ttkManifoldCheck.h>

vtkStandardNewMacro(ttkManifoldCheck)

int ttkManifoldCheck::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  manifoldCheck_.setupTriangulation(triangulation);
  manifoldCheck_.setWrapper(this);
 
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);
  

  manifoldCheck_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
  manifoldCheck_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
  manifoldCheck_.execute();

  {
    stringstream msg;
    msg << "[ttkManifoldCheck] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
