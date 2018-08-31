#include                  <ttkManifoldLearning.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkManifoldLearning)

  int ttkManifoldLearning::doIt(vtkTable* input, vtkTable* output){
    Memory m;

    const int numberOfRows=input->GetNumberOfRows();
    const int numberOfColumns=ScalarFields.size();
    vector<double> inputData;
    vector<vtkAbstractArray*> arrays;
    for(auto s : ScalarFields)
      arrays.push_back(input->GetColumnByName(s.data()));
    for(int i=0; i<numberOfRows; ++i){
      for(auto arr : arrays)
        inputData.push_back(arr->GetVariantValue(i).ToDouble());
    }

    outputData_->clear();

    manifoldLearning_.setWrapper(this);
    manifoldLearning_.setInputModulePath(ModulePath);
    manifoldLearning_.setInputModuleName(ModuleName);
    manifoldLearning_.setInputFunctionName(FunctionName);
    manifoldLearning_.setInputMatrixDimensions(numberOfRows, numberOfColumns);
    manifoldLearning_.setInputMatrix(inputData.data());
    manifoldLearning_.setInputMethod(Method);
    manifoldLearning_.setInputNumberOfComponents(NumberOfComponents);
    manifoldLearning_.setInputNumberOfNeighbors(NumberOfNeighbors);
    manifoldLearning_.setOutputComponents(outputData_);
    manifoldLearning_.execute();

    output->ShallowCopy(input);
    for(int i=0; i<NumberOfComponents; ++i){
      vtkSmartPointer<vtkDoubleArray> arr=vtkSmartPointer<vtkDoubleArray>::New();
      arr->SetVoidArray((*outputData_)[i].data(), numberOfRows, 1);
      arr->SetName(to_string(i).data());
      output->AddColumn(arr);
    }

    {
      stringstream msg;
      msg << "[ttkManifoldLearning] Memory usage: " << m.getElapsedUsage() 
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }

bool ttkManifoldLearning::needsToAbort(){
  return GetAbortExecute();
}

int ttkManifoldLearning::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkManifoldLearning] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkManifoldLearning::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkTable* input=vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkTable* output=vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkManifoldLearning] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
