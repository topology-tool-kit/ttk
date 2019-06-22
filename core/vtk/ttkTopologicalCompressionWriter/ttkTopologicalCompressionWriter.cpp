#include <ttkTopologicalCompressionWriter.h>

vtkStandardNewMacro(ttkTopologicalCompressionWriter);

ttkTopologicalCompressionWriter::ttkTopologicalCompressionWriter() {
  // GUI parameters must be initialized.
  CompressionType = (int)ttk::CompressionType::PersistenceDiagram;
  FileName = nullptr;
  ZFPBitBudget = 0;
  ZFPOnly = false;
  Tolerance = 1;
  SQMethod = "";
  Subdivide = false;
  UseTopologicalSimplification = true;
  // ScalarField = "";
  ScalarFieldId = 0;
  SetUseAllCores(true);
}

ttkTopologicalCompressionWriter::~ttkTopologicalCompressionWriter() {
}

int ttkTopologicalCompressionWriter::FillInputPortInformation(
  int, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

void ttkTopologicalCompressionWriter::ComputeTriangulation(vtkImageData *vti) {
  triangulation.setInputData(vti);
  topologicalCompression.setupTriangulation(triangulation.getTriangulation());
}

vtkDataArray *
  ttkTopologicalCompressionWriter::GetInputScalarField(vtkImageData *vti) {
  return ScalarField.length()
           ? vti->GetPointData()->GetArray(ScalarField.data())
           : vti->GetPointData()->GetArray(ScalarFieldId);
}

int ttkTopologicalCompressionWriter::AllocateOutput(
  vtkDataArray *inputScalarField) {
  switch(inputScalarField->GetDataType()) {
    case VTK_CHAR:
      outputScalarField = vtkSmartPointer<vtkCharArray>::New();
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
    default: {
      ttk::Debug d;
      std::stringstream msg;
      msg << "[vtkTopologicalCompression] Unsupported data type :(" << endl;
      d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    }

      // Do nothing.
      return -1;
  }

  outputScalarField->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());
  outputScalarField->SetName(inputScalarField->GetName());
  Modified();

  return 0;
}

void ttkTopologicalCompressionWriter::PerformCompression(
  vtkDataArray *inputScalarField) {
  topologicalCompression.setInputDataPointer(
    inputScalarField->GetVoidPointer(0));
  topologicalCompression.setSQ(SQMethod);
  topologicalCompression.setSubdivide(!Subdivide);
  topologicalCompression.setUseTopologicalSimplification(
    UseTopologicalSimplification);
  topologicalCompression.setZFPOnly(ZFPOnly);
  topologicalCompression.setCompressionType(CompressionType);
  topologicalCompression.setOutputDataPointer(
    outputScalarField->GetVoidPointer(0));
  topologicalCompression.setMaximumError(MaximumError);
  switch(inputScalarField->GetDataType()) {
    vtkTemplateMacro(topologicalCompression.execute<VTK_TT>(Tolerance));
    default: {
      ttk::Debug d;
      std::stringstream msg;
      msg << "[ttkCompressionWriter] Unsupported data type." << std::endl;
      d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    } break;
  }
}

void ttkTopologicalCompressionWriter::WriteData() {
  bool zfpOnly = ZFPOnly;
  double zfpBitBudget = ZFPBitBudget;
  topologicalCompression.setThreadNumber(threadNumber_);

  {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkCompressionWriter] New writing task." << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  if(zfpOnly && (zfpBitBudget > 64 || zfpBitBudget < 1)) {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkTopologicalCompressionReader] Wrong ZFP bit budget for "
           "ZFP-only use."
        << " Aborting" << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    return;
  }

  vtkDataObject *input = GetInput();
  vtkImageData *vti = vtkImageData::SafeDownCast(input);

  vtkDataArray *inputScalarField = GetInputScalarField(vti);

  ComputeTriangulation(vti);

  int res = AllocateOutput(inputScalarField);
  if(res < 0)
    return;

  PerformCompression(inputScalarField);

  {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkCompressionWriter] Compression successful." << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  // Open file.
  FILE *fp;
  if((fp = fopen(FileName, "wb")) == nullptr) {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkCompressionWriter] System IO error while opening the file."
        << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
    return;
  }

  int dt
    = vti->GetPointData()->GetArray(inputScalarField->GetName())->GetDataType();
  double *vp = (double *)vti->GetPointData()
                 ->GetArray(inputScalarField->GetName())
                 ->GetVoidPointer(0);

  topologicalCompression.setFileName(FileName);
  topologicalCompression.WriteToFile<double>(
    fp, CompressionType, ZFPOnly, SQMethod.c_str(), dt, vti->GetExtent(),
    vti->GetSpacing(), vti->GetOrigin(), vp, Tolerance, ZFPBitBudget);

  {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkTopologicalCompression] Wrote to " << FileName << "."
        << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }
}
