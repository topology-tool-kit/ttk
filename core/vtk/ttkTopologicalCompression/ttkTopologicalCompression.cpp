#include "ttkTopologicalCompression.h"
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(ttkTopologicalCompression);

ttkTopologicalCompression::ttkTopologicalCompression() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTopologicalCompression::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTopologicalCompression::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTopologicalCompression::RequestData(vtkInformation *request,
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(input == nullptr) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(input->GetNumberOfPoints() == 0) {
    this->printErr("Input has no point.");
    return -1;
  }
  if(input->GetPointData() == nullptr) {
    this->printErr("Input has no point data.");
    return -1;
  }
  if(output == nullptr) {
    this->printErr("Output pointer is NULL.");
    return -1;
  }
#endif

  // Triangulate
  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = nullptr;

  if(ScalarField.length()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
    this->printMsg("Starting computation on field '" + ScalarField + "'...");
  } else {
    inputScalarField = input->GetPointData()->GetArray(ScalarFieldId);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarField) {
    this->printErr("Input scalar field pointer is NULL.");
    return -1;
  }
#endif

  const auto vertexNumber = inputScalarField->GetNumberOfTuples();

  // allocate the memory for the output scalar field
  vtkSmartPointer<vtkDataArray> outputScalarField{};

  switch(inputScalarField->GetDataType()) {
    case VTK_CHAR:
      outputScalarField = vtkSmartPointer<vtkSignedCharArray>::New();
      break;
    case VTK_DOUBLE:
      outputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
      break;
    case VTK_FLOAT:
      outputScalarField = vtkSmartPointer<vtkFloatArray>::New();
      break;
    case VTK_INT:
      outputScalarField = vtkSmartPointer<vtkIntArray>::New();
      break;
    case VTK_ID_TYPE:
      outputScalarField = vtkSmartPointer<vtkIdTypeArray>::New();
      break;
    default:
      this->printErr("Unsupported data type :(");
      return -1;
  }

  outputScalarField->SetNumberOfTuples(vertexNumber);
  outputScalarField->SetName(inputScalarField->GetName());

  vtkNew<vtkIntArray> outputOffsetField{};
  outputOffsetField->SetNumberOfTuples(vertexNumber);
  outputOffsetField->SetName(ttk::OffsetScalarFieldName);

  this->setCompressionType(CompressionType);
  this->setInputDataPointer(ttkUtils::GetVoidPointer(inputScalarField));
  this->setSQ(SQMethod);
  this->setUseTopologicalSimplification(UseTopologicalSimplification);
  this->setSubdivide(!Subdivide);
  this->setOutputDataPointer(ttkUtils::GetVoidPointer(outputScalarField));
  this->setMaximumError(MaximumError);

  // Call TopologicalCompression
  switch(inputScalarField->GetDataType()) {
    vtkTemplateMacro(this->execute<VTK_TT>(Tolerance));
    default:
      break;
  }

  std::vector<int> voidOffsets = this->getCompressedOffsets();
  for(SimplexId i = 0; i < vertexNumber; ++i)
    outputOffsetField->SetTuple1(i, voidOffsets[i]);

  output->GetPointData()->AddArray(outputScalarField);
  output->GetPointData()->AddArray(outputOffsetField);

  return 1;
}
