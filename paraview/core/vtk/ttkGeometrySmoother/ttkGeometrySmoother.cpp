#include                  <ttkGeometrySmoother.h>

#include <vtkPointData.h>
#include <vtkCharArray.h>

vtkStandardNewMacro(ttkGeometrySmoother)

ttkGeometrySmoother::ttkGeometrySmoother(){

  // init
  NumberOfIterations = 1;
  MaskIdentifier = 0;
  UseInputMask = false;
  InputMask = "";
}

ttkGeometrySmoother::~ttkGeometrySmoother(){
  
}

int ttkGeometrySmoother::doIt(vector<vtkDataSet *> &inputs, 
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
  
  vtkCharArray *inputMaskField = NULL;

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

  // This filter copies the input into a new data-set (smoothed)
  // let's use shallow copies, in order to only duplicate point positions 
  // (before and after). the rest is not changed, pointers are sufficient.
  output->DeepCopy(input);

  void* inputMaskPtr = (inputMaskField)?inputMaskField->GetVoidPointer(0):nullptr;
  // calling the smoothing package
  vtkPoints *inputPointSet = (vtkPointSet::SafeDownCast(input))->GetPoints();
  vtkPoints *outputPointSet = (vtkPointSet::SafeDownCast(output))->GetPoints();
  switch(outputPointSet->GetDataType()){
   
#ifndef _MSC_VER
	  vtkTemplateMacro((
	  {
		  smoother_.setDimensionNumber(3);
	  smoother_.setInputDataPointer(inputPointSet->GetVoidPointer(0));
	  smoother_.setOutputDataPointer(outputPointSet->GetVoidPointer(0));
	  smoother_.setMaskDataPointer(inputMaskPtr);
	  smoother_.smooth<VTK_TT>(NumberOfIterations);
	  }
	  ));
#else
	  vtkTemplateMacro(
	  {
		  smoother_.setDimensionNumber(3);
	  smoother_.setInputDataPointer(inputPointSet->GetVoidPointer(0));
	  smoother_.setOutputDataPointer(outputPointSet->GetVoidPointer(0));
	  smoother_.setMaskDataPointer(inputMaskPtr);
	  smoother_.smooth<VTK_TT>(NumberOfIterations);
	  }
	  );
#endif
  }
  
  {
    stringstream msg;
    msg << "[ttkGeometrySmoother] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
