#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPolyData.h>

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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkGaussianPointCloud::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  auto domain = vtkPolyData::GetData(outputVector);

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

  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(NumberOfSamples + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(NumberOfSamples);

  for(int i = 0; i < NumberOfSamples; i++) {
    offsets->SetTuple1(i, i);
    connectivity->SetTuple1(i, i);
  }
  offsets->SetTuple1(NumberOfSamples, NumberOfSamples);

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  domain->SetVerts(cells);

  return 1;
}
