#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>

#include <ttkPointDataConverter.h>

#include <limits>

vtkStandardNewMacro(ttkPointDataConverter);

ttkPointDataConverter::ttkPointDataConverter() {
  this->setDebugMsgPrefix("PointDataConverter");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointDataConverter::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPointDataConverter::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename A, typename B, typename C>
int ttkPointDataConverter::convert(vtkDataArray *inputData,
                                   vtkDataSet *output) {
  A *input_ptr = static_cast<A *>(inputData->GetVoidPointer(0));
  int n = inputData->GetNumberOfComponents();
  vtkIdType N = inputData->GetNumberOfTuples();
  B *output_ptr = new B[N * n];

  if(UseNormalization) {
    auto type_min = static_cast<double>(std::numeric_limits<B>::min());
    auto type_max = static_cast<double>(std::numeric_limits<B>::max());
    for(int k = 0; k < n; ++k) {
      double *input_limits = inputData->GetRange(k);

      for(vtkIdType i = 0; i < N; ++i) {
        double d = (double)input_ptr[i * n + k];
        d = (d - input_limits[0]) / (input_limits[1] - input_limits[0]);
        d = d * (type_max - type_min) + type_min;
        output_ptr[i * n + k] = (B)d;
      }
    }
  } else
    for(vtkIdType i = 0; i < N * n; ++i)
      output_ptr[i] = (B)input_ptr[i];

  vtkNew<C> outputData;
  outputData->SetName(ScalarField.data());
  outputData->SetNumberOfComponents(n);
  outputData->SetArray(output_ptr, N * n, 0);

  if(ScalarField.length())
    output->GetPointData()->RemoveArray(ScalarField.data());
  else
    output->GetPointData()->RemoveArray(0);
  output->GetPointData()->AddArray(outputData);

  return 0;
}

int ttkPointDataConverter::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  output->ShallowCopy(input);

  // TODO
  vtkDataArray *inputScalarField = nullptr;
  if(ScalarField.length())
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  else
    inputScalarField = input->GetPointData()->GetArray(0);
  if(!inputScalarField)
    return -1;

  auto InputType = inputScalarField->GetDataType();

  bool oldUseNormalization{UseNormalization};
  if(OutputType == SupportedType::Float or OutputType == SupportedType::Double)
    UseNormalization = false;

  if(InputType == VTK_CHAR) {
    if(OutputType == SupportedType::Double)
      convert<char, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<char, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<char, vtkIdType, vtkIdTypeArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<char, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<char, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<char, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<char, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_DOUBLE) {
    if(OutputType == SupportedType::Char)
      convert<double, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<double, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<double, vtkIdType, vtkIdTypeArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<double, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<double, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<double, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<double, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_FLOAT) {
    if(OutputType == SupportedType::Char)
      convert<float, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<float, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<float, vtkIdType, vtkIdTypeArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<float, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<float, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<float, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<float, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_INT) {
    if(OutputType == SupportedType::Char)
      convert<int, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<int, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<int, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<int, vtkIdType, vtkIdTypeArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<int, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<int, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<int, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_ID_TYPE) {
    if(OutputType == SupportedType::Char)
      convert<vtkIdType, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<vtkIdType, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<vtkIdType, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<vtkIdType, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<vtkIdType, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<vtkIdType, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_LONG) {
    if(OutputType == SupportedType::Char)
      convert<long long int, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<long long int, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<long long int, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<long long int, vtkIdType, vtkIdTypeArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<long long int, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<long long int, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<long long int, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<long long int, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_SHORT) {
    if(OutputType == SupportedType::Char)
      convert<short, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<short, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<short, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<short, vtkIdType, vtkIdTypeArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<short, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<short, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<short, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_UNSIGNED_SHORT) {
    if(OutputType == SupportedType::Char)
      convert<unsigned short, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<unsigned short, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<unsigned short, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<unsigned short, vtkIdType, vtkIdTypeArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<unsigned short, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<unsigned short, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedChar)
      convert<unsigned short, unsigned char, vtkUnsignedCharArray>(
        inputScalarField, output);
  } else if(InputType == VTK_UNSIGNED_CHAR) {
    if(OutputType == SupportedType::Char)
      convert<unsigned char, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<unsigned char, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<unsigned char, float, vtkFloatArray>(inputScalarField, output);
    else if(OutputType == SupportedType::IdType)
      convert<unsigned char, vtkIdType, vtkIdTypeArray>(
        inputScalarField, output);
    else if(OutputType == SupportedType::Int)
      convert<unsigned char, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<unsigned char, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<unsigned char, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
  } else {
    this->printErr("Unsupported data type");
  }

  UseNormalization = oldUseNormalization;

  return 1;
}
