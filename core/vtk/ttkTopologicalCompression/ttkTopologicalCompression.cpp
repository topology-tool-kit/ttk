#include                  "ttkTopologicalCompression.h"

vtkStandardNewMacro(ttkTopologicalCompression)

int ttkTopologicalCompression::doIt(
    std::vector<vtkDataSet *> &inputs,
    std::vector<vtkDataSet *> &outputs)
{

  // Prepare IO
  vtkDataSet *input1 = inputs[0];
  vtkDataSet *output1 = outputs[0];

  triangulation_.setWrapper(this);
  topologicalCompression_.setWrapper(this);

  // Triangulate
  triangulation_.setInputData(input1);
  internalTriangulation_ = ttkTriangulation::getTriangulation(input1);
  topologicalCompression_.setupTriangulation(internalTriangulation_);
  Modified();

  if (!input1) return -1;

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output1->ShallowCopy(input1);
  
  // in the following, the target scalar field of the input is replaced in the 
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you 
  // should proceed in the same way.
  vtkDataArray *inputScalarField = nullptr;
  
  if (ScalarField.length()){
    inputScalarField = input1->GetPointData()->GetArray(ScalarField.data());
    std::stringstream msg;
    msg << "[ttkTopologicalCompression] Starting computation on field '" 
      << ScalarField << "'..." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }
  else
    inputScalarField = input1->GetPointData()->GetArray(ScalarFieldId);

  if (!inputScalarField) return -1;
  
  // allocate the memory for the output scalar field
  if (!outputScalarField_) {
    switch(inputScalarField->GetDataType()) {
      case VTK_CHAR:    outputScalarField_ = vtkSmartPointer<vtkCharArray>::New(); // vtkCharArray::New();
        break;
      case VTK_DOUBLE:  outputScalarField_ = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case VTK_FLOAT:   outputScalarField_ = vtkSmartPointer<vtkFloatArray>::New();
        break;
      case VTK_INT:     outputScalarField_ = vtkSmartPointer<vtkIntArray>::New();
        break;

      default:
        {
          std::stringstream msg;
          msg << "[ttkTopologicalCompression] Unsupported data type :(" 
            << std::endl;
          dMsg(std::cerr, msg.str(), fatalMsg);
        }
        break;
    }
  }

  if (!outputOffsetField_) {
    outputOffsetField_ = vtkSmartPointer<vtkIntArray>::New();
  }

  ttk::SimplexId vertexNumber = (ttk::SimplexId) inputScalarField->GetNumberOfTuples();
  outputOffsetField_->SetNumberOfTuples(vertexNumber);
  outputOffsetField_->SetName("Offsets");

  outputScalarField_->SetNumberOfTuples(vertexNumber);
  outputScalarField_->SetName(inputScalarField->GetName());

  topologicalCompression_.setCompressionType(CompressionType);
  topologicalCompression_.setInputDataPointer(
    inputScalarField->GetVoidPointer(0));
  topologicalCompression_.setSQ(SQMethod);
  topologicalCompression_.setUseTopologicalSimplification(
    UseTopologicalSimplification);
  topologicalCompression_.setSubdivide(!Subdivide);
  topologicalCompression_.setOutputDataPointer(
    outputScalarField_->GetVoidPointer(0));
  topologicalCompression_.setMaximumError(MaximumError);
  
  // Call TopologicalCompression
  switch (inputScalarField->GetDataType()) {
    vtkTemplateMacro(
        topologicalCompression_.execute<VTK_TT>(Tolerance));
    default: break;
  }

  std::vector<int> voidOffsets = topologicalCompression_.getCompressedOffsets();
  for (ttk::SimplexId i = 0; i < vertexNumber; ++i)
    outputOffsetField_->SetTuple1(i, voidOffsets[i]);

  output1->GetPointData()->AddArray(outputScalarField_);
  output1->GetPointData()->AddArray(outputOffsetField_);

  return 0;
}
