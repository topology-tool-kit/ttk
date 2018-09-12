#include                  <ttkDimensionReduction.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDimensionReduction)

  int ttkDimensionReduction::doIt(vtkTable* input, vtkTable* output){
    Memory m;

    const int numberOfRows=input->GetNumberOfRows();
    const int numberOfColumns=ScalarFields.size();

#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfRows<=0 or numberOfColumns<=0){
      cerr << "[ttkDimensionReduction] Error : input matrix has invalid dimensions" << endl;
      return -1;
    }
#endif

    vector<double> inputData;
    vector<vtkAbstractArray*> arrays;
    for(auto s : ScalarFields)
      arrays.push_back(input->GetColumnByName(s.data()));
    for(int i=0; i<numberOfRows; ++i){
      for(auto arr : arrays)
        inputData.push_back(arr->GetVariantValue(i).ToDouble());
    }

    outputData_->clear();

    dimensionReduction_.setWrapper(this);
    dimensionReduction_.setInputModulePath(ModulePath);
    dimensionReduction_.setInputModuleName(ModuleName);
    dimensionReduction_.setInputFunctionName(FunctionName);
    dimensionReduction_.setInputMatrixDimensions(numberOfRows, numberOfColumns);
    dimensionReduction_.setInputMatrix(inputData.data());
    dimensionReduction_.setInputMethod(Method);
    dimensionReduction_.setInputNumberOfComponents(NumberOfComponents);
    dimensionReduction_.setInputNumberOfNeighbors(NumberOfNeighbors);
    dimensionReduction_.setInputIsDeterministic(IsDeterministic);
    dimensionReduction_.setOutputComponents(outputData_);
    dimensionReduction_.execute();

    if(KeepAllDataArrays)
      output->ShallowCopy(input);

    for(int i=0; i<NumberOfComponents; ++i){
      string s = "Component_" + to_string(i);
      vtkSmartPointer<vtkDoubleArray> arr=vtkSmartPointer<vtkDoubleArray>::New();
      arr->SetVoidArray((*outputData_)[i].data(), numberOfRows, 1);
      arr->SetName(s.data());
      output->AddColumn(arr);
    }

    {
      stringstream msg;
      msg << "[ttkDimensionReduction] Memory usage: " << m.getElapsedUsage() 
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }

bool ttkDimensionReduction::needsToAbort(){
  return GetAbortExecute();
}

int ttkDimensionReduction::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkDimensionReduction] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkDimensionReduction::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkTable* input=vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkTable* output=vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkDimensionReduction] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
