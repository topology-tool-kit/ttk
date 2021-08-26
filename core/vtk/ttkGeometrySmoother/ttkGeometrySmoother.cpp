#include <ttkGeometrySmoother.h>

#include <ttkUtils.h>
#include <vtkInformation.h>

#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnsignedCharArray.h>

vtkStandardNewMacro(ttkGeometrySmoother);

ttkGeometrySmoother::ttkGeometrySmoother() {
  this->setDebugMsgPrefix("GeometrySmoother");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkGeometrySmoother::~ttkGeometrySmoother() {
}

int ttkGeometrySmoother::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkGeometrySmoother::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkGeometrySmoother::RequestData(vtkInformation *ttkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {

  auto inputPointSet = vtkPointSet::GetData(inputVector[0]);
  auto outputPointSet = vtkPointSet::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(inputPointSet);
  if(!triangulation)
    return 0;
  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputMaskField = ttkAlgorithm::GetOptionalArray(
    ForceInputMaskScalarField, 1, ttk::MaskScalarFieldName, inputPointSet);

  // This filter copies the input into a new data-set (smoothed)
  // let's use shallow copies, in order to only duplicate point positions
  // (before and after). the rest is not changed, pointers are sufficient.
  outputPointSet->DeepCopy(inputPointSet);

  // calling the smoothing package
  auto inputPoints = inputPointSet->GetPoints();
  auto outputPoints = outputPointSet->GetPoints();

  this->setDimensionNumber(3);
  this->setInputDataPointer(ttkUtils::GetVoidPointer(inputPoints));
  this->setOutputDataPointer(ttkUtils::GetVoidPointer(outputPoints));

  if(inputMaskField) {
    this->setMaskDataPointer(ttkUtils::GetVoidPointer(inputMaskField));
  }

  switch(outputPoints->GetDataType()) {
    vtkTemplateMacro(
      this->smooth<VTK_TT>(triangulation->getData(), this->NumberOfIterations));
  }

  return 1;
}
