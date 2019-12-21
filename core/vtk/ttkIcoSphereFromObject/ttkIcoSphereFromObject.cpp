#include <ttkIcoSphereFromObject.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIcoSphereFromObject);

ttkIcoSphereFromObject::ttkIcoSphereFromObject() : ttkIcoSphere() {
  this->setDebugMsgPrefix("IcoSphereFromObject");
  this->SetNumberOfInputPorts(1);
}
ttkIcoSphereFromObject::~ttkIcoSphereFromObject() {
}

int ttkIcoSphereFromObject::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  else
    return 0;
  return 1;
}

int ttkIcoSphereFromObject::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  auto input = vtkDataObject::GetData(inputVector[0], 0);

  double bounds[6] = {0, 0, 0, 0, 0, 0};
  if(input->IsA("vtkMultiBlockDataSet")) {
    auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
    inputAsMB->GetBounds(bounds);
  } else if(input->IsA("vtkDataSet")) {
    auto inputAsDS = vtkDataSet::SafeDownCast(input);
    inputAsDS->GetBounds(bounds);
  } else {
    this->printErr("Unable to compute bounding box of "
                   + string(input->GetClassName()));
    return 0;
  }

  double dx = bounds[1] - bounds[0];
  double dy = bounds[3] - bounds[2];
  double dz = bounds[5] - bounds[4];
  float center[3]
    = {(float)(bounds[0] + dx * 0.5), (float)(bounds[2] + dy * 0.5),
       (float)(bounds[4] + dz * 0.5)};

  this->SetCenter(center);

  this->SetRadius(this->Scale * std::sqrt(dx * dx + dy * dy + dz * dz) / 2.0);

  return this->ttkIcoSphere::RequestData(request, inputVector, outputVector);
}
