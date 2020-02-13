#include <ttkExtract.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkExtract);

ttkExtract::ttkExtract() {

  this->setDebugMsgPrefix("Extract");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkExtract::~ttkExtract(){};

int ttkExtract::GetVtkDataTypeName(int outputType, std::string &dataTypeName) {
  switch(outputType) {
    case 0: // auto
    case VTK_MULTIBLOCK_DATA_SET: {
      dataTypeName = "vtkMultiBlockDataSet";
      break;
    }
    case VTK_UNSTRUCTURED_GRID: {
      dataTypeName = "vtkUnstructuredGrid";
      break;
    }
    case VTK_IMAGE_DATA: {
      dataTypeName = "vtkImageData";
      break;
    }
    default:
      return 0;
  }

  return 1;
}

int ttkExtract::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  else
    return 0;
  return 1;
}

int ttkExtract::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    if(this->Mode == 0) {
      std::string outputDataTypeName = "";
      if(!this->GetVtkDataTypeName(this->OutputType, outputDataTypeName)) {
        this->printErr("Unsupported output type");
        return 0;
      }
      info->Set(vtkDataObject::DATA_TYPE_NAME(), outputDataTypeName.data());
    } else
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  } else
    return 0;
  return 1;
}

// =============================================================================
// RequestInformation
// =============================================================================
int ttkExtract::RequestInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector *outputVector) {
  if(this->Mode == 0 && this->GetOutputType() == VTK_IMAGE_DATA) {
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // Bounds
    auto imageBounds = this->GetImageBounds();
    outInfo->Set(vtkStreamingDemandDrivenPipeline::BOUNDS(), imageBounds, 6);

    // Extent
    int wholeExtent[6]
      = {(int)imageBounds[0], (int)(imageBounds[1] - imageBounds[0]),
         (int)imageBounds[2], (int)(imageBounds[3] - imageBounds[2]),
         (int)imageBounds[4], (int)(imageBounds[5] - imageBounds[4])};
    outInfo->Set(
      vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent, 6);
  }

  return 1;
}

template <class dataType>
int markVertices(std::vector<int> &oIndex_to_mIndex_map,
                 vtkAbstractArray *inputArray,
                 const size_t &nValues,
                 const std::vector<double> &valuesAsD) {
  auto data = (dataType *)inputArray->GetVoidPointer(0);

  std::vector<dataType> valuesAsDT(nValues);
  for(size_t i = 0; i < nValues; i++)
    valuesAsDT[i] = (dataType)valuesAsD[i];

  for(size_t i = 0, nPoints = inputArray->GetNumberOfTuples(); i < nPoints;
      i++) {
    const auto &dataValue = data[i];
    bool contained = false;
    for(size_t j = 0; j < nValues; j++)
      if(valuesAsDT[j] == dataValue) {
        contained = true;
        break;
      }
    if(contained)
      oIndex_to_mIndex_map[i] = -2;
  }
  return 1;
}

// =============================================================================
// RequestData
// =============================================================================
int ttkExtract::RequestData(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  // Get Input to Output
  auto input = vtkDataObject::GetData(inputVector[0]);
  auto output = vtkDataObject::GetData(outputVector);

  int mode = this->Mode;

  // -------------------------------------------------------------------------
  // Replace Variables in ExpressionString (e.g. {time[2]})
  // -------------------------------------------------------------------------
  std::string finalExpressionString;
  {
    std::string errorMsg;
    if(!ttkUtils::replaceVariables(this->GetExpressionString(),
                                   input->GetFieldData(), finalExpressionString,
                                   errorMsg)) {
      this->printErr(errorMsg);
      return 0;
    }
  }

  std::vector<std::string> valuesAsStrings;
  ttkUtils::stringListToVector(finalExpressionString, valuesAsStrings);

  const size_t nValues = valuesAsStrings.size();
  std::vector<double> valuesAsD(nValues);

  for(size_t i = 0; i < nValues; i++) {
    try {
      double value = stod(valuesAsStrings[i]);
      valuesAsD[i] = value;
    } catch(std::invalid_argument &e) {
      this->printErr("Unable to convert element '" + valuesAsStrings[i]
                     + "' of input std::string '" + finalExpressionString
                     + "' to double");
      return 0;
    }
  }

  if(this->Mode == 0) {
    // print state
    {
      std::string outputDataTypeName = "";
      this->GetVtkDataTypeName(this->OutputType, outputDataTypeName);
      this->printMsg(ttk::debug::Separator::L1);
      this->printMsg({{"Extraction Mode", "Block"},
                      {"Output Type", outputDataTypeName},
                      {"Indicies", "[" + finalExpressionString + "]"}});
      this->printMsg(ttk::debug::Separator::L2);
    }

    this->printMsg("Extracting blocks", 0, ttk::debug::LineMode::REPLACE);
    auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
    if(!inputAsMB) {
      this->printErr("Block mode requires vtkMultiBlockDataSet input");
      return 0;
    }

    if(this->OutputType == 0) {
      auto outputAsMB = vtkMultiBlockDataSet::SafeDownCast(output);

      for(size_t i = 0; i < nValues; i++) {
        size_t blockIndex = (size_t)valuesAsD[i];
        if(0 <= blockIndex && blockIndex < inputAsMB->GetNumberOfBlocks()) {
          auto block = inputAsMB->GetBlock(blockIndex);
          auto copy
            = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());
          copy->ShallowCopy(block);
          outputAsMB->SetBlock(i, copy);
        } else {
          this->printErr("Index " + std::to_string(blockIndex)
                         + " out of range");
          return 0;
        }
      }
    } else {
      if(nValues != 1) {
        this->printErr(
          "If OutputType is specified then only one block can be extracted");
        return 0;
      }

      size_t blockIndex = (size_t)valuesAsD[0];
      if(0 <= blockIndex && blockIndex < inputAsMB->GetNumberOfBlocks()) {
        auto block = inputAsMB->GetBlock(blockIndex);
        if(block->GetDataObjectType() != this->GetOutputType()) {
          this->printErr("BlockType does not match OutputType");
          return 0;
        }
        output->ShallowCopy(block);

        // copy iteration info
        auto iterationInfo = vtkDoubleArray::SafeDownCast(
          inputAsMB->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));
        if(iterationInfo)
          output->GetFieldData()->AddArray(iterationInfo);
      } else {
        this->printErr("Index " + std::to_string(blockIndex) + " out of range");
        return 0;
      }
    }

    this->printMsg("Extracting blocks", 1);

  } else if(mode == 1) {
    int cellMode = this->GetCellMode();

    vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
    std::string inputArrayName = inputArray->GetName();

    if(this->GetInputArrayAssociation(0, inputVector) != 0) {
      this->printErr(
        "Extraction is currently only supported based on point data");
      return 0;
    }

    this->printMsg(
      {{"Ext. Mode", "Geometry"},
       {"Cell Mode",
        std::string(cellMode == 0 ? "All" : cellMode == 1 ? "Any" : "Sub")},
       {"Condition",
        "'" + inputArrayName + "' in [" + finalExpressionString + "]"}});
    this->printMsg(ttk::debug::Separator::L1);

    auto inputAsUG = vtkUnstructuredGrid::SafeDownCast(input);
    auto outputAsUG = vtkUnstructuredGrid::SafeDownCast(output);
    if(!inputAsUG || !outputAsUG) {
      this->printErr(
        "Geometry mode requires 'vtkUnstructuredGrid' input/output");
      return 0;
    }

    if(inputArray->GetNumberOfComponents() != 1) {
      this->printErr("Data array '" + inputArrayName
                     + "' must have only one component");
      return 0;
    }
    const size_t nPoints = inputArray->GetNumberOfTuples();
    auto inputPD = inputAsUG->GetPointData();
    auto inputCD = inputAsUG->GetCellData();

    // Input Topo
    size_t nInCells = inputAsUG->GetNumberOfCells();
    auto inTopologyData = static_cast<vtkIdType *>(
      ttkUtils::GetVoidPointer(inputAsUG->GetCells()->GetData()));

    // Marked Points:
    //     -1: does not satisfy condition
    //     -2: satisfies condition
    //     -3: satisfies condition extended to cell mode
    //    >=0: mIndex
    std::vector<int> oIndex_to_mIndex_map(nPoints, -1);
    std::vector<bool> markedCells(nInCells, false);
    int nMarkedPoints = 0;

    // Output Topo Stats
    size_t outTopologyDataSize = 0;
    size_t nOutCells = 0;

    // ---------------------------------------------------------------------
    // Marking Vertices
    // ---------------------------------------------------------------------
    {
      this->printMsg(
        "Marking vertices/cells", 0, ttk::debug::LineMode::REPLACE);
      ttk::Timer t;

      // Mark vertices that satisfy condition
      switch(inputArray->GetDataType()) {
        vtkTemplateMacro(markVertices<VTK_TT>(
          oIndex_to_mIndex_map, inputArray, nValues, valuesAsD));
        // vtkTemplateMacro({
        //     auto data = (VTK_TT*) inputArray->GetVoidPointer(0);

        //     std::vector<VTK_TT> valuesAsVTK_TT( nValues );
        //     for(size_t i=0; i<nValues; i++)
        //         valuesAsVTK_TT[i] = (VTK_TT)valuesAsD[i];

        //     auto contains = []( const std::vector<VTK_TT>& array, const
        //     size_t& arraySize, const VTK_TT& value ){
        //         for(size_t i=0; i<arraySize; i++)
        //             if(array[i] == value)
        //                 return true;
        //         return false;
        //     };

        //     for(size_t i=0; i<nPoints; i++)
        //         if( contains( valuesAsVTK_TT, nValues, data[i] ) )
        //             oIndex_to_mIndex_map[i] = -2;
        // });
      }

      // Mark vertices based on cell mode and determine outTopo stats
      {

        if(cellMode == 0) {
          // All
          for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
            size_t nVertices = inTopologyData[inTopoIndex];

            bool all = true;
            for(size_t j = 1; j <= nVertices; j++)
              if(oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]] == -1) {
                all = false;
                break;
              }

            if(all) {
              nOutCells++;
              markedCells[i] = true;
              for(size_t j = 1; j <= nVertices; j++)
                oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]] = -3;
              outTopologyDataSize += 1 + nVertices;
            }

            inTopoIndex += nVertices + 1;
          }

        } else if(cellMode == 1) {
          // Any

          // Mark Border Vertices
          for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
            size_t nVertices = inTopologyData[inTopoIndex];
            size_t mVertices = 0;

            for(size_t j = 1; j <= nVertices; j++) {
              auto &mIndex
                = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
              if(mIndex == -2 || mIndex == -3)
                mVertices++;
            }

            // Mark border vertices
            if(mVertices > 0) {
              nOutCells++;
              markedCells[i] = true;
              for(size_t j = 1; j <= nVertices; j++) {
                auto &mIndex
                  = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
                if(mIndex == -1)
                  mIndex = -4; // Border Vertex
                else if(mIndex == -2)
                  mIndex = -3; // Inner Vertex
              }
              outTopologyDataSize += 1 + nVertices;
            }

            inTopoIndex += nVertices + 1;
          }

          for(size_t i = 0; i < nPoints; i++) {
            auto &mIndex = oIndex_to_mIndex_map[i];
            if(mIndex == -4)
              mIndex = -3;
          }

        } else if(cellMode == 2) {
          // Sub
          for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
            size_t nVertices = inTopologyData[inTopoIndex];
            size_t mVertices = 0;

            for(size_t j = 1; j <= nVertices; j++) {
              const auto &mIndex
                = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
              if(mIndex == -2 || mIndex == -3)
                mVertices++;
            }

            if(mVertices > 0) {
              nOutCells++;
              markedCells[i] = true;
              for(size_t j = 1; j <= nVertices; j++) {
                auto &mIndex
                  = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
                if(mIndex == -2)
                  mIndex = -3;
              }
              outTopologyDataSize += 1 + mVertices;
            }

            inTopoIndex += nVertices + 1;
          }
        }
      }

      for(size_t i = 0; i < nPoints; i++)
        if(oIndex_to_mIndex_map[i] == -2 || oIndex_to_mIndex_map[i] == -3)
          oIndex_to_mIndex_map[i] = nMarkedPoints++;

      this->printMsg("Marking vertices/cells", 1, t.getElapsedTime());
    }

    // ---------------------------------------------------------------------
    // Extracting Points
    // ---------------------------------------------------------------------
    {
      ttk::Timer t;
      this->printMsg("Extracting vertices", 0, ttk::debug::LineMode::REPLACE);

      auto inPoints = inputAsUG->GetPoints();
      auto inPointCoords = (float *)inPoints->GetVoidPointer(0);

      auto outPoints = vtkSmartPointer<vtkPoints>::New();
      outPoints->SetNumberOfPoints(nMarkedPoints);
      auto outPointCoords = (float *)outPoints->GetVoidPointer(0);
      outputAsUG->SetPoints(outPoints);

      // Extract point coordinates
      {
        for(size_t i = 0, j = 0; i < nPoints; i++) {
          const auto &mIndex = oIndex_to_mIndex_map[i];
          if(mIndex >= 0) {
            int iOffset = i * 3;
            outPointCoords[j++] = inPointCoords[iOffset++];
            outPointCoords[j++] = inPointCoords[iOffset++];
            outPointCoords[j++] = inPointCoords[iOffset];
          }
        }
      }

      // Extract point data
      {
        size_t nArrays = inputPD->GetNumberOfArrays();
        auto outputPD = outputAsUG->GetPointData();
        for(size_t arrayIndex = 0; arrayIndex < nArrays; arrayIndex++) {
          auto iArray = inputPD->GetAbstractArray(arrayIndex);
          auto oArray
            = vtkSmartPointer<vtkAbstractArray>::Take(iArray->NewInstance());
          oArray->SetName(iArray->GetName());
          oArray->SetNumberOfComponents(iArray->GetNumberOfComponents());
          oArray->SetNumberOfTuples(nMarkedPoints);

          switch(iArray->GetDataType()) {
            vtkTemplateMacro({
              auto iArrayData = (VTK_TT *)iArray->GetVoidPointer(0);
              auto oArrayData = (VTK_TT *)oArray->GetVoidPointer(0);

              for(size_t i = 0, j = 0; i < nPoints; i++) {
                const auto &mIndex = oIndex_to_mIndex_map[i];
                if(mIndex >= 0)
                  oArrayData[j++] = iArrayData[i];
              }
            });
          }

          outputPD->AddArray(oArray);
        }
      }

      this->printMsg(
        "Extracting vertices (#" + std::to_string(nMarkedPoints) + ")", 1,
        t.getElapsedTime());
    }

    // ---------------------------------------------------------------------
    // Extracting Cells
    // ---------------------------------------------------------------------
    {
      ttk::Timer t;
      this->printMsg("Extracting cells", 0, ttk::debug::LineMode::REPLACE);

      const int types[20] = {VTK_EMPTY_CELL,
                             VTK_VERTEX,
                             VTK_LINE,
                             VTK_TRIANGLE,
                             VTK_TETRA,
                             VTK_CONVEX_POINT_SET, // 5
                             VTK_CONVEX_POINT_SET, // 6
                             VTK_CONVEX_POINT_SET, // 7
                             VTK_VOXEL,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET,
                             VTK_CONVEX_POINT_SET};

      std::vector<int> markedCellTypes(nOutCells);
      auto outTopology = vtkSmartPointer<vtkIdTypeArray>::New();
      outTopology->SetNumberOfValues(outTopologyDataSize);
      auto outTopologyData = (vtkIdType *)outTopology->GetVoidPointer(0);

      if(cellMode == 0 || cellMode == 1) {
        for(size_t i = 0, inTopoIndex = 0, outTopoIndex = 0, outCellIndex = 0;
            i < nInCells; i++) {
          const size_t nVertices = inTopologyData[inTopoIndex];
          if(markedCells[i]) {
            markedCellTypes[outCellIndex++] = types[nVertices];
            outTopologyData[outTopoIndex] = nVertices;
            for(size_t j = 1; j <= nVertices; j++)
              outTopologyData[outTopoIndex + j]
                = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
            outTopoIndex += nVertices + 1;
          }

          inTopoIndex += nVertices + 1;
        }
      } else {
        // Sub
        for(size_t i = 0, inTopoIndex = 0, outTopoIndex = 0, outCellIndex = 0;
            i < nInCells; i++) {
          const size_t nVertices = inTopologyData[inTopoIndex];

          if(markedCells[i]) {
            size_t mVertices = 0;
            for(size_t j = 1; j <= nVertices; j++) {
              const auto &mIndex
                = oIndex_to_mIndex_map[inTopologyData[inTopoIndex + j]];
              if(mIndex > -1) {
                outTopologyData[outTopoIndex + 1 + mVertices] = mIndex;
                mVertices++;
              }
            }
            outTopologyData[outTopoIndex] = mVertices;
            markedCellTypes[outCellIndex++] = types[mVertices];

            outTopoIndex += mVertices + 1;
          }

          inTopoIndex += nVertices + 1;
        }
      }

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetCells(nOutCells, outTopology);
      outputAsUG->SetCells(markedCellTypes.data(), cellArray);

      // Extract cell data
      {
        size_t nArrays = inputCD->GetNumberOfArrays();
        auto outputCD = outputAsUG->GetCellData();
        for(size_t arrayIndex = 0; arrayIndex < nArrays; arrayIndex++) {
          auto iArray = inputCD->GetAbstractArray(arrayIndex);
          auto oArray
            = vtkSmartPointer<vtkAbstractArray>::Take(iArray->NewInstance());
          oArray->SetName(iArray->GetName());
          oArray->SetNumberOfComponents(iArray->GetNumberOfComponents());
          oArray->SetNumberOfTuples(nOutCells);

          switch(iArray->GetDataType()) {
            vtkTemplateMacro({
              auto iArrayData = (VTK_TT *)iArray->GetVoidPointer(0);
              auto oArrayData = (VTK_TT *)oArray->GetVoidPointer(0);

              for(size_t i = 0, j = 0; i < nInCells; i++) {
                if(markedCells[i])
                  oArrayData[j++] = iArrayData[i];
              }
            });
          }

          outputCD->AddArray(oArray);
        }
      }

      this->printMsg("Extracting cells (#" + std::to_string(nOutCells) + ")", 1,
                     t.getElapsedTime());
    }

    outputAsUG->GetFieldData()->ShallowCopy(inputAsUG->GetFieldData());
  }

  // Performance
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete", 1, globalTimer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
