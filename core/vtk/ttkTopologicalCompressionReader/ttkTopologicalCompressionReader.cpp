#include <ttkMacros.h>
#include <ttkTopologicalCompressionReader.h>
#include <ttkUtils.h>

#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>

vtkStandardNewMacro(ttkTopologicalCompressionReader);

ttkTopologicalCompressionReader::ttkTopologicalCompressionReader() {
  SetNumberOfInputPorts(0);
  SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("TopologicalCompressionReader");
}

int ttkTopologicalCompressionReader::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkTopologicalCompressionReader::RequestInformation(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  if(FileName == nullptr) {
    return 1;
  }

  FILE *fp = fopen(FileName, "rb"); // binary mode
  if(fp == nullptr) {
    return 1;
  }

  // Fill spacing, origin, extent, scalar type
  // L8 tolerance, ZFP factor
  const auto res = this->ReadMetaData(fp);
  if(res != 0) {
    return 1;
  }
  DataScalarType = this->getDataScalarType();
  for(int i = 0; i < 3; ++i) {
    DataSpacing[i] = this->getDataSpacing()[i];
    DataOrigin[i] = this->getDataOrigin()[i];
    DataExtent[i] = this->getDataExtent()[i];
    DataExtent[3 + i] = this->getDataExtent()[3 + i];
  }

  //  ReadMetaData(fp);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkDataObject::SPACING(), DataSpacing.data(), 3);
  outInfo->Set(vtkDataObject::ORIGIN(), DataOrigin.data(), 3);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), DataExtent.data(), 6);

  int numberOfVertices = 1;
  for(int i = 0; i < 3; ++i)
    numberOfVertices *= (1 + DataExtent[2 * i + 1] - DataExtent[2 * i]);
  outInfo->Set(vtkDataObject::FIELD_NUMBER_OF_TUPLES(), numberOfVertices);

  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, DataScalarType, 1);

  rewind(fp);
  fclose(fp);
  fp = nullptr;

  return 1;
}

int ttkTopologicalCompressionReader::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  // Initialize
  if(FileName == nullptr) {
    return 1;
  }
  FILE *fp = fopen(FileName, "rb"); // binary mode

  if(fp == nullptr) {
    return 1;
  }

  this->setFileName(FileName);
  const auto res = this->ReadMetaData(fp);
  if(res != 0) {
    return 1;
  }
  DataScalarType = this->getDataScalarType();
  for(int i = 0; i < 3; ++i) {
    DataSpacing[i] = this->getDataSpacing()[i];
    DataOrigin[i] = this->getDataOrigin()[i];
    DataExtent[i] = this->getDataExtent()[i];
    DataExtent[3 + i] = this->getDataExtent()[3 + i];
  }
  int nx = 1 + DataExtent[1] - DataExtent[0];
  int ny = 1 + DataExtent[3] - DataExtent[2];
  int nz = 1 + DataExtent[5] - DataExtent[4];
  int vertexNumber = nx * ny * nz;
  ZFPOnly = this->getZFPOnly();

  vtkNew<vtkImageData> mesh{};
  BuildMesh(mesh);

  auto triangulation = ttkAlgorithm::GetTriangulation(mesh);
  this->preconditionTriangulation(triangulation);

  int status{0};
  ttkTemplateMacro(triangulation->getType(),
                   status = this->ReadFromFile(
                     fp, *static_cast<TTK_TT *>(triangulation->getData())));
  if(status != 0) {
    vtkWarningMacro("Failure when reading compressed TTK file");
  }

  mesh->GetPointData()->RemoveArray(0);
  mesh->GetPointData()->SetNumberOfTuples(vertexNumber);

  vtkNew<vtkDoubleArray> decompressed{};
  decompressed->SetNumberOfTuples(vertexNumber);
  const auto &name = this->getDataArrayName();
  if(!name.empty()) {
    decompressed->SetName(name.data());
  } else {
    decompressed->SetName("Decompressed");
  }
  for(int i = 0; i < vertexNumber; ++i) {
    decompressed->SetTuple1(i, decompressedData_[i]);
  }

  // decompressed->SetVoidArray(, vertexNumber, 0);
  mesh->GetPointData()->AddArray(decompressed);

  if(SQMethodInt != 1 && SQMethodInt != 2 && !ZFPOnly) {
    vtkNew<vtkIntArray> vertexOffset{};
    vertexOffset->SetNumberOfTuples(vertexNumber);
    vertexOffset->SetName(this->GetOrderArrayName(decompressed).data());
    const auto &voidOffsets = this->getDecompressedOffsets();
    for(size_t i = 0; i < voidOffsets.size(); ++i)
      vertexOffset->SetTuple1(i, voidOffsets[i]);
    //    vertexOffset->SetVoidArray(
    //      topologicalCompression.getDecompressedOffsets(), vertexNumber, 0);
    mesh->GetPointData()->AddArray(vertexOffset);
  }

  this->printMsg("Read " + std::to_string(mesh->GetNumberOfPoints())
                 + " vertice(s), " + std::to_string(mesh->GetNumberOfCells())
                 + " cell(s).");

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

  // Set the output
  auto output = vtkImageData::GetData(outputVector);
  output->ShallowCopy(mesh);

  return 1;
}

vtkImageData *ttkTopologicalCompressionReader::GetOutput() {
  // copied from ParaView's vtkImageAlgorithm::GetOutput(int port)
  return vtkImageData::SafeDownCast(this->GetOutputDataObject(0));
}

void ttkTopologicalCompressionReader::BuildMesh(vtkImageData *mesh) const {
  int nx = 1 + DataExtent[1] - DataExtent[0];
  int ny = 1 + DataExtent[3] - DataExtent[2];
  int nz = 1 + DataExtent[5] - DataExtent[4];
  mesh->SetDimensions(nx, ny, nz);
  mesh->SetSpacing(DataSpacing[0], DataSpacing[1], DataSpacing[2]);
  mesh->SetOrigin(DataOrigin[0], DataOrigin[1], DataOrigin[2]);
  mesh->AllocateScalars(DataScalarType, 2);
  mesh->GetPointData()->SetNumberOfTuples(nx * ny * nz);
}
