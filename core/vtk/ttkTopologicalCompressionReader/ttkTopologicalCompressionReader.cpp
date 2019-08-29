#include <ttkTopologicalCompressionReader.h>

vtkStandardNewMacro(ttkTopologicalCompressionReader);

/** Override **/

ttkTopologicalCompressionReader::ttkTopologicalCompressionReader() {

  FileName = nullptr;
  ZFPOnly = false;
  fp = nullptr;

  DataScalarType = VTK_DOUBLE;
  DataExtent[0] = 0;
  DataExtent[1] = 0;
  DataExtent[2] = 0;
  DataExtent[3] = 0;
  DataExtent[4] = 0;
  DataExtent[5] = 0;
  DataOrigin[0] = 0.0;
  DataOrigin[1] = 0.0;
  DataOrigin[2] = 0.0;
  DataSpacing[0] = 1.0;
  DataSpacing[1] = 1.0;
  DataSpacing[2] = 1.0;

  SetNumberOfInputPorts(0);
  SetNumberOfOutputPorts(1);
}

int ttkTopologicalCompressionReader::FillOutputPortInformation(
  int port, vtkInformation *info) {
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

int ttkTopologicalCompressionReader::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  if(FileName == nullptr) {
    return 1;
  }
  if(fp != nullptr) {
    return 1;
  }
  fp = fopen(FileName, "rb"); // binary mode
  if(fp == nullptr) {
    return 1;
  }

  // Fill spacing, origin, extent, scalar type
  // L8 tolerance, ZFP factor
  topologicalCompression.ReadMetaData<double>(fp);
  DataScalarType = topologicalCompression.getDataScalarType();
  for(int i = 0; i < 3; ++i) {
    DataSpacing[i] = topologicalCompression.getDataSpacing()[i];
    DataOrigin[i] = topologicalCompression.getDataOrigin()[i];
    DataExtent[i] = topologicalCompression.getDataExtent()[i];
    DataExtent[3 + i] = topologicalCompression.getDataExtent()[3 + i];
  }

  //  ReadMetaData(fp);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkDataObject::SPACING(), DataSpacing, 3);
  outInfo->Set(vtkDataObject::ORIGIN(), DataOrigin, 3);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), DataExtent, 6);

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
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // Initialize
  if(FileName == nullptr) {
    return 1;
  }
  if(fp != nullptr) {
    return 1;
  }
  fp = fopen(FileName, "rb"); // binary mode
  if(fp == nullptr) {
    return 1;
  }

  topologicalCompression.setFileName(FileName);
  topologicalCompression.ReadMetaData<double>(fp);
  DataScalarType = topologicalCompression.getDataScalarType();
  for(int i = 0; i < 3; ++i) {
    DataSpacing[i] = topologicalCompression.getDataSpacing()[i];
    DataOrigin[i] = topologicalCompression.getDataOrigin()[i];
    DataExtent[i] = topologicalCompression.getDataExtent()[i];
    DataExtent[3 + i] = topologicalCompression.getDataExtent()[3 + i];
  }
  int nx = 1 + DataExtent[1] - DataExtent[0];
  int ny = 1 + DataExtent[3] - DataExtent[2];
  int nz = 1 + DataExtent[5] - DataExtent[4];
  int vertexNumber = nx * ny * nz;
  ZFPOnly = topologicalCompression.getZFPOnly();

  BuildMesh();

  triangulation.setInputData(mesh);
  topologicalCompression.setupTriangulation(triangulation.getTriangulation());

  topologicalCompression.ReadFromFile<double>(fp);

  mesh->GetPointData()->RemoveArray(0);
  mesh->GetPointData()->SetNumberOfTuples(vertexNumber);

  decompressed = vtkSmartPointer<vtkDoubleArray>::New();
  decompressed->SetNumberOfTuples(vertexNumber);
  decompressed->SetName("Decompressed");
  std::vector<double> decompressdeData
    = topologicalCompression.getDecompressedData();
  for(int i = 0; i < vertexNumber; ++i)
    decompressed->SetTuple1(i, decompressdeData[i]);
  // decompressed->SetVoidArray(, vertexNumber, 0);
  mesh->GetPointData()->AddArray(decompressed);

  if(SQMethod != 1 && SQMethod != 2 && !ZFPOnly) {
    vertexOffset = vtkSmartPointer<vtkIntArray>::New();
    vertexOffset->SetNumberOfTuples(vertexNumber);
    vertexOffset->SetName(ttk::OffsetScalarFieldName);
    std::vector<int> voidOffsets
      = topologicalCompression.getDecompressedOffsets();
    for(int i = 0; i < vertexNumber; ++i)
      vertexOffset->SetTuple1(i, voidOffsets[i]);
    //    vertexOffset->SetVoidArray(
    //      topologicalCompression.getDecompressedOffsets(), vertexNumber, 0);
    mesh->GetPointData()->AddArray(vertexOffset);
  }

  {
    ttk::Debug d;
    std::stringstream msg;
    msg << "[ttkCompressionReader] Read " << mesh->GetNumberOfPoints()
        << " vertice(s)" << std::endl;
    msg << "[ttkCompressionReader] Read " << mesh->GetNumberOfCells()
        << " cell(s)" << std::endl;
    d.dMsg(std::cout, msg.str(), ttk::Debug::infoMsg);
  }

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);

  // Set the output
  vtkImageData *output
    = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  output->ShallowCopy(mesh);

  return 1;
}

void ttkTopologicalCompressionReader::BuildMesh() {
  int nx = 1 + DataExtent[1] - DataExtent[0];
  int ny = 1 + DataExtent[3] - DataExtent[2];
  int nz = 1 + DataExtent[5] - DataExtent[4];
  mesh = vtkSmartPointer<vtkImageData>::New();
  mesh->SetDimensions(nx, ny, nz);
  mesh->SetSpacing(DataSpacing[0], DataSpacing[1], DataSpacing[2]);
  mesh->SetOrigin(DataOrigin[0], DataOrigin[1], DataOrigin[2]);
  mesh->AllocateScalars(DataScalarType, 2);
  mesh->GetPointData()->SetNumberOfTuples(nx * ny * nz);

  int numberOfVertices = 1;
  for(int i = 0; i < 3; ++i)
    numberOfVertices *= (1 + DataExtent[2 * i + 1] - DataExtent[2 * i]);
}
