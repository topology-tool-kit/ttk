#include <ttkScalarFieldSmoother.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldSmoother);

ttkScalarFieldSmoother::ttkScalarFieldSmoother() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldSmoother::~ttkScalarFieldSmoother() {
}

int ttkScalarFieldSmoother::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkScalarFieldSmoother::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkScalarFieldSmoother::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);

  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalarField)
    return 0;

  if(inputScalarField->GetNumberOfComponents() != 1) {
    printErr("Invalid scalar field ("
             + std::to_string(inputScalarField->GetNumberOfComponents())
             + " components)");
    return 0;
  }

  vtkDataArray *inputMaskField = ttkAlgorithm::GetOptionalArray(
    ForceInputMaskScalarField, 1, ttk::MaskScalarFieldName, input);

  // preparing the output
  vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputScalarField->NewInstance());
  outputArray->SetName(inputScalarField->GetName());
  outputArray->SetNumberOfComponents(1);
  outputArray->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());

  // This filter copies the input into a new data-set (smoothed)
  // let's use shallow copies, in order to only duplicate point positions
  // (before and after). the rest is not changed, pointers are sufficient.
  output->ShallowCopy(input);
  output->GetPointData()->AddArray(outputArray);

  printMsg("Starting computation...");
  printMsg(
    {{"  Scalar Array", inputScalarField->GetName()},
     {"  Mask Array", inputMaskField ? inputMaskField->GetName() : "None"},
     {"  #iterations", std::to_string(NumberOfIterations)}});

  void *inputMaskPtr
    = (inputMaskField) ? ttkUtils::GetVoidPointer(inputMaskField) : nullptr;

  this->setDimensionNumber(inputScalarField->GetNumberOfComponents());
  this->setInputDataPointer(ttkUtils::GetVoidPointer(inputScalarField));
  this->setOutputDataPointer(ttkUtils::GetVoidPointer(outputArray));
  this->setMaskDataPointer(inputMaskPtr);

  // calling the smoothing package
  ttkVtkTemplateMacro(
    inputScalarField->GetDataType(), triangulation->getType(),
    (this->smooth<VTK_TT, TTK_TT>(
      (TTK_TT *)triangulation->getData(), NumberOfIterations)));

  return 1;
}
