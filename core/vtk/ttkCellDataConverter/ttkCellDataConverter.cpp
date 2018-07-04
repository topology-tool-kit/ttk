#include                  <ttkCellDataConverter.h>

#ifdef _WIN32
#include<ciso646>
#endif

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCellDataConverter)

ttkCellDataConverter::ttkCellDataConverter(){
}

ttkCellDataConverter::~ttkCellDataConverter(){
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkCellDataConverter::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkCellDataConverter::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkCellDataConverter] " << progress*100 
	<< "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

template<typename A, typename B, typename C>
int ttkCellDataConverter::convert(vtkDataArray* inputData, vtkDataSet* output){
  A* input_ptr=static_cast<A*>(inputData->GetVoidPointer(0));
  int n=inputData->GetNumberOfComponents();
  vtkIdType N=inputData->GetNumberOfTuples();
  B* output_ptr=new B[N*n];

  if(UseNormalization){
    double type_min=numeric_limits<B>::min();
    double type_max=numeric_limits<B>::max();

    for(int k=0; k<n; ++k){
      double* input_limits=inputData->GetRange();
      
      for(vtkIdType i=0; i<N; ++i){
	double d=(double)input_ptr[i*n+k];
	d=(d-input_limits[0])/(input_limits[1]-input_limits[0]);
	d=d*(type_max-type_min)+type_min;
	output_ptr[i*n+k]=(B)d;
      }
    }
  }
  else
    for(vtkIdType i=0; i<N*n; ++i) output_ptr[i]=(B)input_ptr[i];
    
  vtkSmartPointer<C> outputData=vtkSmartPointer<C>::New();
  outputData->SetName(ScalarField.data());
  outputData->SetNumberOfComponents(n);
  outputData->SetArray(output_ptr,N*n,0);
      
  if(ScalarField.length()) output->GetCellData()->RemoveArray(ScalarField.data());    
  else output->GetCellData()->RemoveArray(0);
  output->GetCellData()->AddArray(outputData);

  return 0;
}

int ttkCellDataConverter::doIt(vtkDataSet *input, vtkDataSet *output){
  output->ShallowCopy(input);
  
  vtkDataArray* inputScalarField=nullptr;
  if(ScalarField.length())
    inputScalarField=input->GetCellData()->GetArray(ScalarField.data());
  else
    inputScalarField=input->GetCellData()->GetArray(0);
  if(!inputScalarField)
    return -1;

  bool oldUseNormalization{UseNormalization};
  if(OutputType==SupportedType::Float or OutputType==SupportedType::Double)
    UseNormalization=false;

  switch(inputScalarField->GetDataType()){
#ifndef _MSC_VER
    vtkTemplateMacro(({
    if(OutputType==SupportedType::Float)
      convert<VTK_TT,float,vtkFloatArray>(inputScalarField,output);
    else if(OutputType==SupportedType::Int)
      convert<VTK_TT,int,vtkIntArray>(inputScalarField,output);
    else if(OutputType==SupportedType::IdType)
      convert<VTK_TT,vtkIdType,vtkIdTypeArray>(inputScalarField,output);
    else if(OutputType==SupportedType::UnsignedShort)
      convert<VTK_TT,unsigned short,vtkUnsignedShortArray>(inputScalarField,output);
    else if(OutputType==SupportedType::UnsignedChar)
      convert<VTK_TT,unsigned char,vtkUnsignedCharArray>(inputScalarField,output);
        }));
#else
    vtkTemplateMacro({
    if(OutputType==SupportedType::Float)
      convert<VTK_TT,float,vtkFloatArray>(inputScalarField,output);
    else if(OutputType==SupportedType::Int)
      convert<VTK_TT,int,vtkIntArray>(inputScalarField,output);
    else if(OutputType==SupportedType::IdType)
      convert<VTK_TT,vtkIdType,vtkIdTypeArray>(inputScalarField,output);
    else if(OutputType==SupportedType::UnsignedShort)
      convert<VTK_TT,unsigned short,vtkUnsignedShortArray>(inputScalarField,output);
    else if(OutputType==SupportedType::UnsignedChar)
      convert<VTK_TT,unsigned char,vtkUnsignedCharArray>(inputScalarField,output);
        });
#endif
  }

  UseNormalization=oldUseNormalization;
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkCellDataConverter::RequestData(vtkInformation *request, 
				   vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkCellDataConverter] Memory usage: " << m.getElapsedUsage() 
	<< " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 1;
}
