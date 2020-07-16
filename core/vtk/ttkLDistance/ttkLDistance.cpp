#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "ttkLDistance.h"
#include <ttkUtils.h>

vtkStandardNewMacro(ttkLDistance);

ttkLDistance::ttkLDistance() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkLDistance::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkLDistance::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkLDistance::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {

  const auto input1 = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  // Test validity of datasets (must present the same number of points).
#ifndef TTK_ENABLE_KAMIKAZE
  if(!input1) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(!input1->GetNumberOfPoints()) {
    this->printErr("Input has no point.");
    return -1;
  }
  if(!input1->GetPointData()) {
    this->printErr("Input has no point data.");
    return -1;
  }
  if(!output) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }
#endif

  // Use a pointer-base copy for the input.
  output->ShallowCopy(input1);

  // In the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField1 = NULL;
  vtkDataArray *inputScalarField2 = NULL;

  if(ScalarField1.length())
    inputScalarField1 = input1->GetPointData()->GetArray(ScalarField1.data());
  else
    inputScalarField1 = input1->GetPointData()->GetArray(ScalarFieldId1);

  if(ScalarField2.length())
    inputScalarField2 = input1->GetPointData()->GetArray(ScalarField2.data());
  else
    inputScalarField2 = input1->GetPointData()->GetArray(ScalarFieldId2);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField1 || !inputScalarField2
     || inputScalarField1->GetDataType() != inputScalarField2->GetDataType()) {
    this->printErr("Input scalar fields are NULL or have different types.");
    return -1;
  }
#endif

  // Allocate memory for the output scalar field, based on the first input.
  vtkSmartPointer<vtkDataArray> outputScalarField{
    inputScalarField1->NewInstance()};

  ttk::SimplexId numberOfPoints = input1->GetNumberOfPoints();

  outputScalarField->SetNumberOfTuples(numberOfPoints);
  outputScalarField->SetName(DistanceFieldName.data());
  output->GetPointData()->AddArray(outputScalarField);

  // Calling the executing package.
  switch(inputScalarField1->GetDataType()) {
    vtkTemplateMacro(this->execute<VTK_TT>(
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField1)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField2)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalarField)),
      DistanceType, numberOfPoints));
  }

  return 1;
}
