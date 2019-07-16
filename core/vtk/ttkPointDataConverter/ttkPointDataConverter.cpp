#include <ttkPointDataConverter.h>

#ifdef _WIN32
#include <ciso646>
#endif

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPointDataConverter)

  ttkPointDataConverter::ttkPointDataConverter() {
  OutputType = 0;
  UseNormalization = true;
  UseAllCores = true;
}

ttkPointDataConverter::~ttkPointDataConverter() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkPointDataConverter::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPointDataConverter::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkPointDataConverter] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
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
    double type_min = numeric_limits<B>::min();
    double type_max = numeric_limits<B>::max();
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

  vtkSmartPointer<C> outputData = vtkSmartPointer<C>::New();
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

int ttkPointDataConverter::doIt(vtkDataSet *input, vtkDataSet *output) {
  output->ShallowCopy(input);

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
  } else if(InputType == VTK_SHORT) {
    if(OutputType == SupportedType::Char)
      convert<short, char, vtkCharArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Double)
      convert<short, double, vtkDoubleArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Float)
      convert<short, float, vtkFloatArray>(inputScalarField, output);
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
    else if(OutputType == SupportedType::Int)
      convert<unsigned char, int, vtkIntArray>(inputScalarField, output);
    else if(OutputType == SupportedType::Short)
      convert<unsigned char, short, vtkShortArray>(inputScalarField, output);
    else if(OutputType == SupportedType::UnsignedShort)
      convert<unsigned char, unsigned short, vtkUnsignedShortArray>(
        inputScalarField, output);
  } else {
    stringstream msg;
    msg << "[ttkCellDataConverter] Unsupported data type :(" << endl;
    dMsg(cerr, msg.str(), fatalMsg);
  }

  UseNormalization = oldUseNormalization;
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPointDataConverter::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkPointDataConverter] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
