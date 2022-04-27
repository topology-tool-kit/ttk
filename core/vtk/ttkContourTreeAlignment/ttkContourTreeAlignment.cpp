#include <ttkContourTreeAlignment.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformationVector.h>
#include <vtkLongLongArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkContourTreeAlignment)

  ttkContourTreeAlignment::ttkContourTreeAlignment() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

// Specify the input data type of each input port
int ttkContourTreeAlignment::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

// Specify the data object type of each output port
int ttkContourTreeAlignment::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkContourTreeAlignment::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  //==================================================================================================================
  // Print status
  { this->printMsg("RequestData"); }

  //==================================================================================================================
  // Prepare input
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto inputMB = vtkMultiBlockDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  size_t n = inputMB->GetNumberOfBlocks();

  //==================================================================================================================
  // print input info
  this->printMsg(
    {{"#Trees", std::to_string(n)}, {"Seed", std::to_string(RandomSeed)}});

  //==================================================================================================================
  // extract topologies
  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Extracting topologies for " + std::to_string(n) + " trees",
                 debug::Priority::VERBOSE);
  if(this->debugLevel_ < static_cast<int>(debug::Priority::VERBOSE)) {
    this->printMsg("Extracting topologies for " + std::to_string(n) + " trees.",
                   0, debug::LineMode::REPLACE);
  }

  int scalarType = -1;
  vector<void *> scalars(n); // scalar type will be determined dynamically
  vector<int *> regionSizes(n);
  vector<int *> segmentationIds(n);
  vector<int *> segmentations;
  vector<size_t> segSizes;
  vector<long long *> topologies(n);
  vector<size_t> nVertices(n);
  vector<size_t> nEdges(n);

  for(size_t i = 0; i < n; i++) {
    auto contourTree = vtkUnstructuredGrid::SafeDownCast(inputMB->GetBlock(i));
    vtkDataSet *segmentation = nullptr;

    if(!contourTree) {
      auto pair = vtkMultiBlockDataSet::SafeDownCast(inputMB->GetBlock(i));

      if(!pair || pair->GetNumberOfBlocks() != 2) {
        this->printErr("block " + std::to_string(i)
                       + " is not an unstructured grid or a pair of "
                         "unstructured grid and segmentation.");
        return 0;
      }

      contourTree = vtkUnstructuredGrid::SafeDownCast(pair->GetBlock(0));
      segmentation = vtkDataSet::SafeDownCast(pair->GetBlock(1));

      if(!contourTree) {
        this->printErr("block " + std::to_string(i)
                       + " is not an unstructured grid or a pair of "
                         "unstructured grid and segmentation.");
        return 0;
      }
    }

    this->printMsg("Block " + std::to_string(i) + " read from multiblock.",
                   debug::Priority::VERBOSE);
    if(this->debugLevel_ < static_cast<int>(debug::Priority::VERBOSE)) {
      this->printMsg("Extracting topologies for " + std::to_string(n)
                       + " trees (" + std::to_string(i) + "/"
                       + std::to_string(n) + ")",
                     (float)i / (float)n, debug::LineMode::REPLACE);
    }

    nVertices[i] = contourTree->GetNumberOfPoints();
    nEdges[i] = contourTree->GetNumberOfCells();

    // auto pData = contourTree->GetPointData();
    // auto scalarArray = pData->GetArray( "Scalar" );

    if(this->GetInputArrayAssociation(0, contourTree) != 0) {
      printErr("Scalar array not point data.");
      return 0;
    }
    auto scalarArray = this->GetInputArrayToProcess(0, contourTree);
    if(!scalarArray) {
      printErr("No Point Array \"Scalar\" found in contour tree.");
      return 0;
    }
    scalars[i] = ttkUtils::GetVoidPointer(scalarArray);
    scalarType = scalarArray->GetDataType();

    this->printMsg(
      "Scalar Array read from point data.", debug::Priority::VERBOSE);

    // auto cData = contourTree->GetCellData();

    // auto regionArray = cData->GetArray( "RegionSize" );
    if(this->GetInputArrayAssociation(1, contourTree) != 1) {
      printErr("Region size array not cell data.");
      return 0;
    }
    auto regionArray = this->GetInputArrayToProcess(1, contourTree);
    if(!regionArray) {
      printErr("No Cell Array \"RegionSize\" found in contour tree.");
      return 0;
    }
    regionSizes[i] = (int *)ttkUtils::GetVoidPointer(regionArray);
    this->printMsg(
      "RegionSize Array read from cell data.", debug::Priority::VERBOSE);

    // auto segArray = cData->GetArray( "SegmentationId" );
    if(this->GetInputArrayAssociation(2, contourTree) != 1) {
      printErr("Segmentation Id array not cell data.");
      return 0;
    }
    auto segArray = this->GetInputArrayToProcess(2, contourTree);
    if(!segArray) {
      printErr("No Cell Array \"SegmentationId\" found in contour tree.");
      return 0;
    }
    segmentationIds[i] = (int *)ttkUtils::GetVoidPointer(segArray);
    this->printMsg(
      "SegmentationId Array read from cell data.", debug::Priority::VERBOSE);

    if(segmentation) {
      auto segArray_segmentation
        = this->GetInputArrayToProcess(3, segmentation);
      // auto segArray_segmentation =
      // segmentation->GetPointData()->GetArray("SegmentationId");
      if(!segArray_segmentation) {
        printErr("No Cell Array \"SegmentationId\" found in segmentation.");
        return 0;
      }
      segmentations.push_back(
        (int *)ttkUtils::GetVoidPointer(segArray_segmentation));
      this->printMsg("Segmentation read.", debug::Priority::VERBOSE);
      segSizes.push_back(segArray_segmentation->GetNumberOfValues());
    }

    auto cells = contourTree->GetCells();
    auto cellSizes
      = (vtkIdType *)ttkUtils::GetVoidPointer(cells->GetOffsetsArray());
    for(int cIdx = 0; cIdx < contourTree->GetNumberOfCells(); cIdx++) {
      if(cellSizes[cIdx + 1] - cellSizes[cIdx] != 2) {
        printErr("cell " + std::to_string(cIdx) + " of block "
                 + std::to_string(i)
                 + " not a line (input not a contour tree).");
        return 0;
      }
    }
    topologies[i]
      = (long long *)ttkUtils::GetVoidPointer(cells->GetConnectivityArray());
  }

  if(scalarType < 0 || n < 1)
    return 1;

  //==================================================================================================================
  // print status
  this->printMsg(
    "For " + std::to_string(n) + " trees topologies and data extracted.",
    debug::Priority::VERBOSE);
  if(this->debugLevel_ < static_cast<int>(debug::Priority::VERBOSE)) {
    this->printMsg(
      "Extracting topologies for " + std::to_string(n) + " trees", 1);
  }

  for(auto s : segSizes) {
    this->printMsg(std::to_string(s));
  }

  this->printMsg("Starting alignment computation.");

  //==================================================================================================================
  // Compute Tree Alignment
  vector<float> outputVertices;
  vector<int> outputEdges;
  vector<long long> outputFrequencies;
  vector<long long> outputVertexIds;
  vector<long long> outputBranchIds;
  vector<long long> outputSegmentationIds;
  vector<long long> outputArcIds;

  this->setArcMatchMode(ArcMatchMode);
  if(MatchTime)
    this->setAlignmenttreeType(ttk::cta::lastMatchedValue);
  else
    this->setAlignmenttreeType(AlignmenttreeType);
  this->setWeightArcMatch(WeightArcMatch);
  this->setWeightCombinatorialMatch(WeightCombinatorialMatch);
  this->setWeightScalarValueMatch(WeightScalarValueMatch);
  this->setDebugLevel(this->debugLevel_);
  this->setThreadNumber(this->threadNumber_);

  int success = false;
  switch(scalarType) {
    vtkTemplateMacro({
      success = this->execute<VTK_TT>(
        scalars, regionSizes, segmentationIds, topologies, nVertices, nEdges,
        segmentations, segSizes,

        outputVertices, outputFrequencies, outputVertexIds, outputBranchIds,
        outputSegmentationIds, outputArcIds, outputEdges,

        RandomSeed);
    });
  }
  if(!success) {
    printErr("base layer execution failed.");
    return 0;
  }

  //==================================================================================================================
  // print status
  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Alignment computed.");
  this->printMsg("Writing paraview/vtk output.", 0, debug::LineMode::REPLACE);

  //==================================================================================================================
  // write paraview output

  // Output
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto alignmentTree = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Prep Array
  auto prepArray = [](vtkAbstractArray *array, const string &name,
                      size_t nComponents, size_t nTuples) {
    array->SetName(name.data());
    array->SetNumberOfComponents(nComponents);
    array->SetNumberOfTuples(nTuples);
  };

  // Points
  size_t nOutputVertices = outputVertices.size();
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nOutputVertices);
  alignmentTree->SetPoints(points);

  // Point Data
  auto freq = vtkSmartPointer<vtkLongLongArray>::New();
  prepArray(freq, "Frequency", 1, nOutputVertices);
  auto freqData = (long long *)ttkUtils::GetVoidPointer(freq);

  auto scalar = vtkSmartPointer<vtkFloatArray>::New();
  prepArray(scalar, "Scalar", 1, nOutputVertices);
  auto scalarData = (float *)ttkUtils::GetVoidPointer(scalar);

  auto vertexIDs = vtkSmartPointer<vtkLongLongArray>::New();
  prepArray(vertexIDs, "VertexIDs", n, nOutputVertices);
  auto vertexIDsData = (long long *)ttkUtils::GetVoidPointer(vertexIDs);

  auto branchIDs = vtkSmartPointer<vtkLongLongArray>::New();
  prepArray(branchIDs, "BranchIDs", 1, nOutputVertices);
  auto branchIDsData = (long long *)ttkUtils::GetVoidPointer(branchIDs);

  auto segmentationIDs = vtkSmartPointer<vtkLongLongArray>::New();
  prepArray(segmentationIDs, "segmentationIDs", n, nOutputVertices);
  auto segmentationIDsData
    = (long long *)ttkUtils::GetVoidPointer(segmentationIDs);

  auto arcIDs = vtkSmartPointer<vtkLongLongArray>::New();
  prepArray(arcIDs, "arcIDs", n, nOutputVertices - 1);
  auto arcIDsData = (long long *)ttkUtils::GetVoidPointer(arcIDs);

  auto pointCoords = (float *)ttkUtils::GetVoidPointer(points);
  for(size_t i = 0; i < nOutputVertices; i++) {
    pointCoords[i * 3] = i;
    pointCoords[i * 3 + 1] = outputVertices[i];
    pointCoords[i * 3 + 2] = 0;

    freqData[i] = outputFrequencies[i];
    scalarData[i] = outputVertices[i];
    branchIDsData[i] = outputBranchIds[i];
  }
  for(size_t i = 0; i < outputVertexIds.size(); i++)
    vertexIDsData[i] = outputVertexIds[i];
  for(size_t i = 0; i < outputSegmentationIds.size(); i++)
    segmentationIDsData[i] = outputSegmentationIds[i];
  for(size_t i = 0; i < outputArcIds.size(); i++)
    arcIDsData[i] = outputArcIds[i];

  // Add Data to Points
  auto pointData = alignmentTree->GetPointData();
  pointData->AddArray(freq);
  pointData->AddArray(scalar);
  pointData->AddArray(vertexIDs);
  pointData->AddArray(branchIDs);
  pointData->AddArray(segmentationIDs);

  // Cells
  size_t nOutputEdges = outputEdges.size() / 2;
  auto cellArray = vtkSmartPointer<vtkCellArray>::New();
  // cellArray->SetCells(nOutputEdges, cells);

  // auto cellConnectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  auto cellConnectivity = cellArray->GetConnectivityArray();
  cellConnectivity->SetNumberOfTuples(nOutputEdges * 2);
  auto cellIds = (vtkIdType *)ttkUtils::GetVoidPointer(cellConnectivity);

  // auto cellOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
  auto cellOffsets = cellArray->GetOffsetsArray();
  cellOffsets->SetNumberOfTuples(nOutputEdges + 1);
  auto cellOffsetsData = (vtkIdType *)ttkUtils::GetVoidPointer(cellOffsets);

  for(size_t i = 0; i < nOutputEdges; i++) {
    cellOffsetsData[i] = i * 2;
    cellIds[i * 2] = (vtkIdType)outputEdges[i * 2];
    cellIds[i * 2 + 1] = (vtkIdType)outputEdges[i * 2 + 1];
  }
  cellOffsetsData[nOutputEdges] = 2 * nOutputEdges;

  alignmentTree->SetCells(VTK_LINE, cellArray);

  alignmentTree->GetCellData()->AddArray(arcIDs);

  this->printMsg("Writing paraview/vtk output", 1);

  if(!ExportJSON || ExportPath.length() < 1) {

    //==================================================================================================================
    // print status
    this->printMsg("No JSON output.");

  } else {

    //==================================================================================================================
    // print status
    this->printMsg("Writing JSON alignment to \"" + ExportPath + "\"", 0,
                   debug::LineMode::REPLACE);

    //==================================================================================================================
    // output alignment tree as JSON file

    std::ofstream fileJSON;
    fileJSON.open(ExportPath + "/alignment.json");

    std::vector<std::vector<int>> upEdges(nOutputVertices);
    std::vector<std::vector<int>> downEdges(nOutputVertices);

    std::vector<std::vector<int>> alignmentIDs;
    for(size_t i = 0; i < n; i++) {
      std::vector<int> vertices_i(nVertices[i], -1);
      alignmentIDs.push_back(vertices_i);
    }

    for(size_t i = 0; i < nOutputEdges; i++) {
      int id1 = outputEdges[i * 2];
      int id2 = outputEdges[i * 2 + 1];
      float v1 = scalar->GetValue(id1);
      float v2 = scalar->GetValue(id2);
      if(v1 > v2) {
        downEdges[id1].push_back(i);
        upEdges[id2].push_back(i);
      } else {
        upEdges[id1].push_back(i);
        downEdges[id2].push_back(i);
      }
    }

    fileJSON << "{\n";

    fileJSON << "  \"nodes\": [\n";

    for(size_t i = 0; i < nOutputVertices; i++) {
      fileJSON << "    {";
      fileJSON << "\"scalar\": " << scalar->GetValue(i) << ", ";
      fileJSON << "\"frequency\": " << freq->GetValue(i) << ", ";
      for(size_t j = 0; j < n; j++) {
        if(vertexIDs->GetComponent(i, j) >= 0)
          alignmentIDs[j][vertexIDs->GetComponent(i, j)] = i;
      }
      fileJSON << "\"segmentationIDs\": [";
      for(size_t j = 0; j < n; j++) {
        fileJSON << segmentationIDs->GetComponent(i, j);
        if(j < n - 1)
          fileJSON << ",";
      }
      fileJSON << "], ";
      fileJSON << "\"upEdgeIDs\": [";
      for(size_t j = 0; j < upEdges[i].size(); j++) {
        fileJSON << upEdges[i][j];
        if(j < upEdges[i].size() - 1)
          fileJSON << ",";
      }
      fileJSON << "], ";
      fileJSON << "\"downEdgeIDs\": [";
      for(size_t j = 0; j < downEdges[i].size(); j++) {
        fileJSON << downEdges[i][j];
        if(j < downEdges[i].size() - 1)
          fileJSON << ",";
      }
      fileJSON << "]";
      fileJSON << (i == nOutputVertices - 1 ? "}\n" : "},\n");
    }

    fileJSON << "  ],\n";

    fileJSON << "  \"edges\": [\n";

    for(size_t i = 0; i < nOutputEdges; i++) {
      fileJSON << "    {";
      int id1 = outputEdges[i * 2];
      int id2 = outputEdges[i * 2 + 1];
      fileJSON << "\"node1\": " << id1 << ", \"node2\": " << id2;
      fileJSON << (i == nOutputEdges - 1 ? "}\n" : "},\n");
    }

    fileJSON << "  ]\n";

    fileJSON << "}\n";

    fileJSON.close();

    //==================================================================================================================
    // print status
    this->printMsg("Writing JSON alignment to \"" + ExportPath + "\"", 1);
    this->printMsg("Writing JSON input trees to \"" + ExportPath + "\"", 0,
                   debug::LineMode::REPLACE);

    //==================================================================================================================
    // output original trees as JSON

    for(size_t t = 0; t < n; t++) {

      fileJSON.open(ExportPath + "/tree" + std::to_string(t) + ".json");

      upEdges = std::vector<std::vector<int>>(nVertices[t]);
      downEdges = std::vector<std::vector<int>>(nVertices[t]);

      for(size_t i = 0; i < nEdges[t]; i++) {
        int id1 = topologies[t][i * 2 + 0];
        int id2 = topologies[t][i * 2 + 1];
        float v1 = ((float *)scalars[t])[id1];
        float v2 = ((float *)scalars[t])[id2];
        if(v1 > v2) {
          downEdges[id1].push_back(i);
          upEdges[id2].push_back(i);
        } else {
          upEdges[id1].push_back(i);
          downEdges[id2].push_back(i);
        }
      }

      fileJSON << "{\n";

      fileJSON << "  \"nodes\": [\n";

      bool first = true;
      for(size_t k = 0; k < nOutputVertices; k++) {
        int i = vertexIDs->GetComponent(k, t);
        if(i < 0)
          continue;
        fileJSON << (first ? "    {" : ",\n    {");
        fileJSON << "\"scalar\": " << ((float *)scalars[t])[i] << ", ";
        fileJSON << "\"id\": " << alignmentIDs[t][i] << ", ";
        fileJSON << "\"upEdgeIDs\": [";
        for(size_t j = 0; j < upEdges[i].size(); j++) {
          fileJSON << upEdges[i][j];
          if(j < upEdges[i].size() - 1)
            fileJSON << ",";
        }
        fileJSON << "], ";
        fileJSON << "\"downEdgeIDs\": [";
        for(size_t j = 0; j < downEdges[i].size(); j++) {
          fileJSON << downEdges[i][j];
          if(j < downEdges[i].size() - 1)
            fileJSON << ",";
        }
        fileJSON << "]";
        fileJSON << "}";
        first = false;
      }

      fileJSON << "\n  ],\n";

      fileJSON << "  \"edges\": [\n";

      for(size_t i = 0; i < nEdges[t]; i++) {
        fileJSON << "    {";
        int id1 = alignmentIDs[t][topologies[t][i * 2 + 0]];
        int id2 = alignmentIDs[t][topologies[t][i * 2 + 1]];
        fileJSON << "\"node1\": " << id1 << ", \"node2\": " << id2;
        fileJSON << (i == nEdges[t] - 1 ? "}\n" : "},\n");
      }

      fileJSON << "  ]\n";

      fileJSON << "}\n";

      fileJSON.close();
    }
    this->printMsg("Writing JSON input trees to \"" + ExportPath + "\"", 1);
  }

  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
