#include <ttkGaussianPointCloud.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkGaussianPointCloud)

  int ttkGaussianPointCloud::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  Memory m;

  // Print status
  {
    stringstream msg;
    msg << "[ttkGaussianPointCloud] Generating " << NumberOfSamples
        << " samples in " << Dimension << "D..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Prepare input and output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto domain = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(NumberOfSamples);

  // Set Wrapper
  gaussianPointCloud_.setWrapper(this);

  if(points->GetDataType() == VTK_FLOAT) {
    gaussianPointCloud_.generate<float>(
      Dimension, NumberOfSamples, points->GetVoidPointer(0));
  }

  if(points->GetDataType() == VTK_DOUBLE) {
    gaussianPointCloud_.generate<double>(
      Dimension, NumberOfSamples, points->GetVoidPointer(0));
  }

  domain->SetPoints(points);

  // Print status
  {
    stringstream msg;
    msg << "[ttkGaussianPointCloud] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
