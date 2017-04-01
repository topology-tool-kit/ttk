#include                  <vtkGeometrySmoother.h>

vtkStandardNewMacro(vtkGeometrySmoother)

vtkGeometrySmoother::vtkGeometrySmoother(){

  // init
  NumberOfIterations = 1;
}

vtkGeometrySmoother::~vtkGeometrySmoother(){
  
}

int vtkGeometrySmoother::doIt(vector<vtkDataSet *> &inputs, 
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
  output->DeepCopy(input);

  // calling the smoothing package
  vtkPoints *inputPointSet = (vtkPointSet::SafeDownCast(input))->GetPoints();
  vtkPoints *outputPointSet = (vtkPointSet::SafeDownCast(output))->GetPoints();
  switch(outputPointSet->GetDataType()){
   
    vtkTemplateMacro((
      {
        smoother_.setDimensionNumber(3);
        smoother_.setInputDataPointer(inputPointSet->GetVoidPointer(0));
        smoother_.setOutputDataPointer(outputPointSet->GetVoidPointer(0));
        smoother_.smooth<VTK_TT>(NumberOfIterations);
      }
    ));
  }
  
  {
    stringstream msg;
    msg << "[vtkGeometrySmoother] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}