#include <ttkMacros.h>
#include <ttkTopologicalCompressionWriter.h>
#include <ttkUtils.h>

#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>

vtkStandardNewMacro(ttkTopologicalCompressionWriter);

ttkTopologicalCompressionWriter::ttkTopologicalCompressionWriter() {
  SetNumberOfInputPorts(1);
  this->setDebugMsgPrefix("TopologicalCompressionWriter");
}

int ttkTopologicalCompressionWriter::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkTopologicalCompressionWriter::Write() {

  this->printMsg("New writing task.");

  if(ZFPOnly && ZFPTolerance < 0.0) {
    this->printErr(
      "Wrong ZFP absolute error tolerance for ZFP-only use, aborting.");
    return 0;
  }

  vtkDataObject *input = GetInput();
  vtkImageData *vti = vtkImageData::SafeDownCast(input);

  auto inputScalarField = this->GetInputArrayToProcess(0, input);
  if(inputScalarField == nullptr) {
    return 0;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(vti);
  if(triangulation == nullptr) {
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  const auto inputOffsets
    = this->GetOrderArray(vtkDataSet::SafeDownCast(input), 0, 1, false);

  vtkSmartPointer<vtkDataArray> outputScalarField;

  switch(inputScalarField->GetDataType()) {
    case VTK_SIGNED_CHAR:
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
      // Do nothing.
      return 0;
  }

  outputScalarField->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());
  outputScalarField->SetName(inputScalarField->GetName());
  Modified();

  // manage tolerance (relative % -> absolute)
  std::array<double, 2> sfRange{};
  inputScalarField->GetRange(sfRange.data());
  this->relToAbsZFPTolerance(this->ZFPTolerance, sfRange);

  ttkVtkTemplateMacro(
    inputScalarField->GetDataType(), triangulation->getType(),
    this->execute(
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField)),
      static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(inputOffsets)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalarField)),
      *static_cast<TTK_TT *>(triangulation->getData())));

  this->printMsg("Compression successful.");

  // Open file.
  FILE *fp;
  if((fp = fopen(FileName, "wb")) == nullptr) {
    this->printErr("System IO error while opening the file.");
    return 0;
  }

  int dt = inputScalarField->GetDataType();
  auto vp = static_cast<double *>(ttkUtils::GetVoidPointer(inputScalarField));

  this->setFileName(FileName);
  this->WriteToFile(fp, CompressionType, ZFPOnly, SQMethod.c_str(), dt,
                    vti->GetExtent(), vti->GetSpacing(), vti->GetOrigin(), vp,
                    Tolerance, ZFPTolerance, inputScalarField->GetName());

  this->printMsg("Wrote to " + std::string{FileName} + ".");
  return 1;
}

vtkDataObject *ttkTopologicalCompressionWriter::GetInput() {
  // copied from ParaView's vtkWriter::GetInput()
  if(this->GetNumberOfInputConnections(0) < 1) {
    return nullptr;
  }
  return this->GetExecutive()->GetInputData(0, 0);
}

void ttkTopologicalCompressionWriter::SetInputData(vtkDataObject *input) {
  // copied from ParaView's vtkWriter::SetInputData()
  this->SetInputDataInternal(0, input);
}
