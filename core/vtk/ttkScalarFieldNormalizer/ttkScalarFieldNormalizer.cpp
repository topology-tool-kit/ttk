#include <ttkScalarFieldNormalizer.h>

#include <Geometry.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldNormalizer);

ttkScalarFieldNormalizer::ttkScalarFieldNormalizer() {
  // init
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  setDebugMsgPrefix("ScalarFieldNormalizer");
}

ttkScalarFieldNormalizer::~ttkScalarFieldNormalizer() {
}

int ttkScalarFieldNormalizer::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkScalarFieldNormalizer::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkScalarFieldNormalizer::normalize(vtkDataArray *input,
                                        vtkDataArray *output) const {

  if(!output)
    return -1;
  if(!input)
    return -2;

  double min = 0, max = 0;
  for(SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {

    double value = input->GetTuple1(i);

    if((!i) || (value < min)) {
      min = value;
    }
    if((!i) || (value > max)) {
      max = value;
    }
  }

  for(SimplexId i = 0; i < input->GetNumberOfTuples(); i++) {
    double value = input->GetTuple1(i);

    value = (value - min) / (max - min) + Geometry::powIntTen(-FLT_DIG);

    output->SetTuple1(i, value);
  }

  return 0;
}

int ttkScalarFieldNormalizer::RequestData(vtkInformation *ttkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  // get input scalar field
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(inputArray == nullptr) {
    this->printErr("No such input scalar field");
    return 0;
  }

  vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(inputArray->GetName());
  outputArray->SetNumberOfComponents(1);
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // calling the executing package
  normalize(inputArray, outputArray);

  // prepare the output
  output->ShallowCopy(input);
  output->GetPointData()->AddArray(outputArray);

  return 1;
}
