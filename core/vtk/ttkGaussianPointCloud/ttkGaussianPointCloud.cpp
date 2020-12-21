#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include <ttkGaussianPointCloud.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkGaussianPointCloud);

ttkGaussianPointCloud::ttkGaussianPointCloud() {
  SetNumberOfInputPorts(0);
  SetNumberOfOutputPorts(1);
}

int ttkGaussianPointCloud::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkGaussianPointCloud::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto domain = vtkUnstructuredGrid::GetData(outputVector);

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(NumberOfSamples);

  if(points->GetDataType() == VTK_FLOAT) {
    this->generate<float>(
      Dimension, NumberOfSamples,
      static_cast<float *>(ttkUtils::GetVoidPointer(points)));
  } else if(points->GetDataType() == VTK_DOUBLE) {
    this->generate<double>(
      Dimension, NumberOfSamples,
      static_cast<double *>(ttkUtils::GetVoidPointer(points)));
  }

  domain->SetPoints(points);

  return 1;
}
