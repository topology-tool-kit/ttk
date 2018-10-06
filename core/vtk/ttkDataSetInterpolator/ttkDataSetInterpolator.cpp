#include<ttkDataSetInterpolator.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkDataSetInterpolator)

  int ttkDataSetInterpolator::doIt(vtkDataSet* source, vtkDataSet* target, vtkDataSet* output){

    Memory m;

    vtkSmartPointer<vtkProbeFilter> probe=vtkSmartPointer<vtkProbeFilter>::New();
    probe->SetInputData(target);
    probe->SetSourceData(source);
    probe->Update();
    output->ShallowCopy(probe->GetOutput());

    {
      stringstream msg;
      msg << "[ttkDataSetInterpolator] Memory usage: " << m.getElapsedUsage() 
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }

bool ttkDataSetInterpolator::needsToAbort(){
  return GetAbortExecute();
}

int ttkDataSetInterpolator::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkDataSetInterpolator] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkDataSetInterpolator::RequestData(vtkInformation *request, 
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector){

  Memory m;

  vtkDataSet* source = vtkDataSet::GetData(inputVector[1]);
  vtkDataSet* target = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet* output = vtkDataSet::GetData(outputVector);

  doIt(source, target, output);

  {
    stringstream msg;
    msg << "[ttkDataSetInterpolator] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
