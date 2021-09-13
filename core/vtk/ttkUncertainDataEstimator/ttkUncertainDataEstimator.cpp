#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkPointData.h>

#include <ttkUncertainDataEstimator.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkUncertainDataEstimator);

ttkUncertainDataEstimator::ttkUncertainDataEstimator() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkUncertainDataEstimator::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkUncertainDataEstimator::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkUncertainDataEstimator::RequestData(vtkInformation *ttkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  // Unified bound fields
  auto boundFields = vtkDataSet::GetData(outputVector, 0);
  // Histograms
  auto probability = vtkDataSet::GetData(outputVector, 1);
  // Mean field
  auto mean = vtkDataSet::GetData(outputVector, 2);

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Inputs", std::to_string(numInputs)}});

  // Get input datas
  std::vector<vtkDataSet *> input(numInputs);
  for(int i = 0; i < numInputs; i++) {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  // Use a pointer-base copy for the input data
  boundFields->ShallowCopy(input[0]);
  probability->ShallowCopy(input[0]);
  mean->ShallowCopy(input[0]);

  int numFields = 0;
  int numArrays = 0;

  // Get arrays from input datas
  std::vector<vtkDataArray *> inputScalarField;

  for(int i = 0; i < numInputs; i++) {
    numArrays = input[i]->GetPointData()->GetNumberOfArrays();
    numFields += numArrays;
    for(int iarray = 0; iarray < numArrays; iarray++) {
      inputScalarField.push_back(input[i]->GetPointData()->GetArray(iarray));
    }
  }

  this->printMsg(std::vector<std::vector<std::string>>{
    {"#Fields", std::to_string(numFields)}});

  for(int i = 0; i < numFields; i++) {

    // Check if inputs have the same data type and the same number of points
    if(inputScalarField[i]->GetDataType()
       != inputScalarField[0]->GetDataType()) {
      this->printErr("Inputs of different data types.");
      return -3;
    }
    if(inputScalarField[i]->GetNumberOfTuples()
       != inputScalarField[0]->GetNumberOfTuples()) {
      this->printErr("Inputs with different number of points.");
      return -2;
    }
    // Check if all the inputs are here
    if(!inputScalarField[i])
      return -1;
  }
  // Allocate the memory for the output bound scalar fields
  vtkSmartPointer<vtkDataArray> outputLowerBoundScalarField{
    inputScalarField[0]->NewInstance()};
  vtkSmartPointer<vtkDataArray> outputUpperBoundScalarField{
    inputScalarField[0]->NewInstance()};
  outputLowerBoundScalarField->SetName("lowerBoundField");
  outputUpperBoundScalarField->SetName("upperBoundField");

  this->printMsg(
    std::vector<std::vector<std::string>>{{"#Bins", std::to_string(BinCount)}});

  std::vector<vtkNew<vtkDoubleArray>> probabilityScalarField(BinCount);

  // Create DoubleArray objects and link them to the data set
  for(int b = 0; b < BinCount; b++) {
    probabilityScalarField[b]->SetNumberOfTuples(input[0]->GetNumberOfPoints());
    probability->GetPointData()->AddArray(probabilityScalarField[b]);
    probabilityScalarField[b]->FillComponent(0, 0.);
  }

  // Mean field data set
  // Allocate new array
  vtkNew<vtkDoubleArray> meanField{};
  meanField->SetNumberOfTuples(input[0]->GetNumberOfPoints());
  meanField->SetName("meanField");
  meanField->FillComponent(0, 0.0);

  mean->GetPointData()->AddArray(meanField);

  // Resize arrays and add them in the output if required
  if(ComputeLowerBound) {
    outputLowerBoundScalarField->SetNumberOfTuples(
      input[0]->GetNumberOfPoints());
    boundFields->GetPointData()->AddArray(outputLowerBoundScalarField);
  } else {
    outputLowerBoundScalarField->SetNumberOfTuples(0);
  }

  if(ComputeUpperBound) {
    outputUpperBoundScalarField->SetNumberOfTuples(
      input[0]->GetNumberOfPoints());
    boundFields->GetPointData()->AddArray(outputUpperBoundScalarField);
  } else {
    outputUpperBoundScalarField->SetNumberOfTuples(0);
  }

  // Calling the executing package

  if(numFields > 0) {
    this->setVertexNumber(boundFields->GetNumberOfPoints());
    this->setBinCount(this->BinCount);

    this->setNumberOfInputs(numFields);
    for(int i = 0; i < numFields; i++) {
      this->setInputDataPointer(
        i, ttkUtils::GetVoidPointer(inputScalarField[i]));
    }

    this->setOutputLowerBoundField(
      ttkUtils::GetVoidPointer(outputLowerBoundScalarField));
    this->setOutputUpperBoundField(
      ttkUtils::GetVoidPointer(outputUpperBoundScalarField));
    this->setOutputMeanField(ttkUtils::GetVoidPointer(meanField));

    for(int b = 0; b < BinCount; b++) {
      this->setOutputProbability(
        b, ttkUtils::GetPointer<double>(probabilityScalarField[b]));
    }

    switch(inputScalarField[0]->GetDataType()) {
      vtkTemplateMacro(this->execute<VTK_TT>());
    }

    for(int b = 0; b < BinCount; b++) {
      std::stringstream name{};
      name << std::setprecision(8) << this->getBinValue(b);
      probabilityScalarField[b]->SetName(name.str().c_str());
    }
  }

  return 1;
}
