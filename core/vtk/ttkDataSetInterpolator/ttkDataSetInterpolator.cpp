#include <ttkDataSetInterpolator.h>
#include <vtkProbeFilter.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkDataSetInterpolator);

ttkDataSetInterpolator::ttkDataSetInterpolator(){
    this->setDebugMsgPrefix("DataSetInterpolator");
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}
ttkDataSetInterpolator::~ttkDataSetInterpolator(){}

int ttkDataSetInterpolator::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  else
    return 0;

  return 1;
}

int ttkDataSetInterpolator::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkDataSetInterpolator::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector
){
  ttk::Timer t;

  auto source = vtkDataSet::GetData(inputVector[1]);
  auto target = vtkDataSet::GetData(inputVector[0]);

  this->printMsg(
    "Computing "+std::to_string(target->GetNumberOfPoints())+" locations",
    0,0,
    ttk::debug::LineMode::REPLACE
  );

  auto output = vtkDataSet::GetData(outputVector);

  auto probe = vtkSmartPointer<vtkProbeFilter>::New();
  probe->SetInputData(target);
  probe->SetSourceData(source);
  probe->Update();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!probe->GetOutput()) {
    this->printErr("Data probe failed.");
    return 0;
  }
#endif

  output->ShallowCopy(probe->GetOutput());

  // add original data arrays
  auto inputPointData = target->GetPointData();
  auto outputPointData = output->GetPointData();

  const size_t numberOfArrays = inputPointData->GetNumberOfArrays();
  for(size_t i = 0; i < numberOfArrays; ++i)
    outputPointData->AddArray(inputPointData->GetArray(i));

  this->printMsg(
    "Computing "+std::to_string(target->GetNumberOfPoints())+" locations",
    1,t.getElapsedTime()
  );

  return 1;
}