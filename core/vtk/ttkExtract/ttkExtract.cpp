#include <ttkExtract.h>

#include <vtkInformationVector.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>

#include <ttkUtils.h>

#include <set>

vtkStandardNewMacro(ttkExtract);

ttkExtract::ttkExtract() {
  this->setDebugMsgPrefix("Extract");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkExtract::~ttkExtract(){};

std::string ttkExtract::GetVtkDataTypeName(const int outputType) const {
  switch(outputType) {
    case -1: // in case of auto return vtkMultiBlockDataSet
    case VTK_MULTIBLOCK_DATA_SET:
      return "vtkMultiBlockDataSet";
    case VTK_UNSTRUCTURED_GRID:
      return "vtkUnstructuredGrid";
    case VTK_POLY_DATA:
      return "vtkPolyData";
    case VTK_IMAGE_DATA:
      return "vtkImageData";
    case VTK_TABLE:
      return "vtkTable";
    default:
      return "";
  }
}

int ttkExtract::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  } else
    return 0;
  return 1;
}

int ttkExtract::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    if(info->Has(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT()))
      info->Remove(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());

    if(this->OutputType != -1) {
      std::string dataTypeName = this->GetVtkDataTypeName(this->OutputType);
      if(dataTypeName.length() < 1) {
        this->printErr("Unsupported output type");
        return 0;
      }
      info->Set(vtkDataObject::DATA_TYPE_NAME(), dataTypeName.data());
    } else
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  } else
    return 0;
  return 1;
}

// =============================================================================
// RequestInformation
// =============================================================================
int ttkExtract::RequestInformation(vtkInformation *,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {

  if(this->ExtractionMode == 0 && this->GetOutputType() == VTK_IMAGE_DATA) {
    auto outInfo = outputVector->GetInformationObject(0);
    outInfo->Set(
      vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), this->ImageExtent, 6);
  }

  return 1;
}

int doubleVectorToString(std::string &str, const std::vector<double> &vec) {
  std::stringstream ss;
  for(auto &v : vec)
    ss << v << ",";
  ss.seekp(-1, ss.cur);
  str = ss.str().substr(0, ss.str().size() - 1);
  return 1;
}

int ttkExtract::ExtractBlocks(vtkDataObject *output,
                              vtkDataObject *input,
                              const std::vector<double> &indices,
                              const bool &extractTuples) const {
  // print state
  std::string indicesString = "";
  doubleVectorToString(indicesString, indices);
  {
    std::string outputDataTypeName = this->GetVtkDataTypeName(this->OutputType);
    std::string extentString = "";
    if(this->OutputType == VTK_IMAGE_DATA) {
      extentString = "[" + std::to_string(this->ImageExtent[0]);
      for(int i = 1; i < 6; i++)
        extentString += "," + std::to_string(this->ImageExtent[i]);
      extentString += "]";
    }

    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg({{"Extraction Mode",
                     "Block" + std::string(extractTuples ? " Tuples" : "")},
                    {"Output Type", outputDataTypeName + extentString},
                    {"Indicies", "[" + indicesString + "]"}});
    this->printMsg(ttk::debug::Separator::L2);
  }

  this->printMsg("Extracting blocks [" + indicesString + "]", 0,
                 ttk::debug::LineMode::REPLACE);

  auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
  if(!inputAsMB) {
    this->printErr("Block mode requires 'vtkMultiBlockDataSet' input.");
    return 0;
  }

  if(this->OutputType != -1 && extractTuples) {
    this->printErr("Block Tuple mode requires auto output.");
    return 0;
  }

  if(this->OutputType == -1) {
    // extract multiple blocks (vtkMultiBlockDataSet input/output only)
    auto outputAsMB = vtkMultiBlockDataSet::SafeDownCast(output);

    if(extractTuples) {
      int nComponents = inputAsMB->GetNumberOfBlocks();

      for(size_t i = 0; i < indices.size(); i++) {
        auto tuple = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        outputAsMB->SetBlock(i, tuple);

        size_t tupleIndex = (size_t)indices[i];
        for(int c = 0; c < nComponents; c++) {
          auto blockAsMB
            = vtkMultiBlockDataSet::SafeDownCast(inputAsMB->GetBlock(c));
          if(!blockAsMB) {
            this->printErr("Block Tuple Mode requires a vtkMultiBlockDataSet "
                           "that contains vtkMultiBlockDataSets as input.");
            return 0;
          }
          if(tupleIndex < 0 || tupleIndex >= blockAsMB->GetNumberOfBlocks()) {
            this->printErr(
              "Index out of range (" + std::to_string(tupleIndex) + "/"
              + std::to_string(blockAsMB->GetNumberOfBlocks()) + ").");
            return 0;
          }

          auto block = blockAsMB->GetBlock(tupleIndex);
          auto copy
            = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());
          copy->ShallowCopy(block);
          tuple->SetBlock(c, copy);
        }
      }
    } else {
      for(size_t i = 0; i < indices.size(); i++) {
        size_t blockIndex = (size_t)indices[i];
        if(blockIndex < 0 || blockIndex >= inputAsMB->GetNumberOfBlocks()) {
          this->printErr("Index out of range (" + std::to_string(blockIndex)
                         + "/" + std::to_string(inputAsMB->GetNumberOfBlocks())
                         + ").");
          return 0;
        }

        auto block = inputAsMB->GetBlock(blockIndex);
        auto copy = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());
        copy->ShallowCopy(block);
        outputAsMB->SetBlock(i, copy);
      }
    }
  } else {
    // extract a single block of specified output type
    if(indices.size() != 1) {
      this->printErr(
        "If OutputType is specified then only one block can be extracted.");
      return 0;
    }

    size_t blockIndex = (size_t)indices[0];
    if(blockIndex < inputAsMB->GetNumberOfBlocks()) {
      auto block = inputAsMB->GetBlock(blockIndex);
      if(output->GetDataObjectType() != block->GetDataObjectType()) {
        this->printErr("BlockType does not match OutputType");
        return 0;
      }
      output->ShallowCopy(block);

    } else {
      this->printErr("Index out of range (" + std::to_string(blockIndex) + "/"
                     + std::to_string(inputAsMB->GetNumberOfBlocks()) + ").");
      return 0;
    }
  }

  // pass field data
  auto inFD = input->GetFieldData();
  auto outFD = output->GetFieldData();
  for(int i = 0; i < inFD->GetNumberOfArrays(); i++)
    outFD->AddArray(inFD->GetArray(i));

  this->printMsg("Extracting blocks [" + indicesString + "]", 1);

  return 1;
}

int ttkExtract::ExtractRows(vtkDataObject *output,
                            vtkDataObject *input,
                            const std::vector<double> &indices) const {
  // print state
  std::string indicesString = "";
  doubleVectorToString(indicesString, indices);
  {
    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg(
      {{"Extraction Mode", "Rows"}, {"Indicies", "[" + indicesString + "]"}});
    this->printMsg(ttk::debug::Separator::L2);
  }

  ttk::Timer t;
  this->printMsg("Extracting rows [" + indicesString + "]", 0,
                 ttk::debug::LineMode::REPLACE);

  size_t nValues = indices.size();
  auto inputAsT = vtkTable::SafeDownCast(input);
  auto outputAsT = vtkTable::SafeDownCast(output);
  if(!inputAsT || !outputAsT) {
    this->printErr("Row mode requires 'vtkTable' input/output.");
    return 0;
  }

  size_t nRows = inputAsT->GetNumberOfRows();
  size_t nCols = inputAsT->GetNumberOfColumns();

  for(size_t j = 0; j < nValues; j++)
    if(((size_t)indices[j]) >= nRows || indices[j] < 0) {
      this->printErr("Index out of range (" + std::to_string((size_t)indices[j])
                     + "/" + std::to_string(nRows) + ").");
      return 0;
    }

  // Extract row at index
  for(size_t i = 0; i < nCols; i++) {
    auto iColumn = inputAsT->GetColumn(i);

    auto oColumn
      = vtkSmartPointer<vtkAbstractArray>::Take(iColumn->NewInstance());
    oColumn->SetName(iColumn->GetName());
    oColumn->SetNumberOfComponents(iColumn->GetNumberOfComponents());
    oColumn->SetNumberOfTuples(nValues);
    for(size_t j = 0; j < nValues; j++)
      oColumn->SetTuple(j, indices[j], iColumn);

    outputAsT->AddColumn(oColumn);
  }

  outputAsT->GetFieldData()->ShallowCopy(inputAsT->GetFieldData());

  this->printMsg(
    "Extracting rows [" + indicesString + "]", 1, t.getElapsedTime());

  return 1;
}

template <typename dataType>
struct ComperatorLessOrEqual {
  static bool Test(const dataType &d0, const dataType &d1) {
    return d0 <= d1;
  }
};
template <typename dataType>
struct ComperatorLess {
  static bool Test(const dataType &d0, const dataType &d1) {
    return d0 < d1;
  }
};
template <typename dataType>
struct ComperatorGreaterOrEqual {
  static bool Test(const dataType &d0, const dataType &d1) {
    return d0 >= d1;
  }
};
template <typename dataType>
struct ComperatorGreater {
  static bool Test(const dataType &d0, const dataType &d1) {
    return d0 > d1;
  }
};
template <typename dataType>
struct ComperatorEqual {
  static bool Test(const dataType &d0, const dataType &d1) {
    return d0 == d1;
  }
};

template <typename dataType, typename ComperatorType>
int testPointDataArray(std::vector<int> &oIndex_to_mIndex_map,

                       const size_t &nPoints,
                       const dataType *inputPointDataArray,
                       const std::vector<dataType> pivotValues,
                       const size_t &threadNumber) {
  const size_t nPivotValues = pivotValues.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nPoints; i++) {
    bool hasToBeMarked = false;
    const dataType &v = inputPointDataArray[i];
    for(size_t j = 0; j < nPivotValues; j++)
      if(ComperatorType::Test(v, pivotValues[j]))
        hasToBeMarked = true;
    if(hasToBeMarked)
      oIndex_to_mIndex_map[i] = -2;
  }

  return 1;
}

template <typename dataType, typename ComperatorType>
int testCellDataArray(std::vector<int> &oIndex_to_mIndex_map,

                      const size_t &nCells,
                      const vtkIdType *inConnectivityList,
                      const dataType *inputCellDataArray,
                      const std::vector<dataType> pivotValues,
                      const size_t &threadNumber) {
  const size_t nPivotValues = pivotValues.size();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nCells; i++) {
    const dataType &v = inputCellDataArray[i];

    bool hasToBeMarked = false;
    for(size_t j = 0; j < nPivotValues; j++)
      if(ComperatorType::Test(v, pivotValues[j]))
        hasToBeMarked = true;
    if(hasToBeMarked) {
      const auto &u = inConnectivityList[i * 3 + 1];
      const auto &w = inConnectivityList[i * 3 + 2];
      oIndex_to_mIndex_map[u] = -2;
      oIndex_to_mIndex_map[w] = -2;
    }
  }

  return 1;
}

template <typename dataType>
int markVerticesBasedOnPointData(std::vector<int> &oIndex_to_mIndex_map,

                                 const std::vector<double> pivotValues,
                                 const size_t &nPoints,
                                 const dataType *inputPointDataArray,
                                 const size_t &validationMode,
                                 const size_t &threadNumber) {
  const size_t nPivotValues = pivotValues.size();
  std::vector<dataType> pivotValuesDT(nPivotValues);
  for(size_t i = 0; i < nPivotValues; i++)
    pivotValuesDT[i] = (dataType)pivotValues[i];

  if(validationMode == 0)
    return testPointDataArray<dataType, ComperatorLess<dataType>>(
      oIndex_to_mIndex_map, nPoints, inputPointDataArray, pivotValuesDT,
      threadNumber);
  else if(validationMode == 1)
    return testPointDataArray<dataType, ComperatorLessOrEqual<dataType>>(
      oIndex_to_mIndex_map, nPoints, inputPointDataArray, pivotValuesDT,
      threadNumber);
  else if(validationMode == 2)
    return testPointDataArray<dataType, ComperatorEqual<dataType>>(
      oIndex_to_mIndex_map, nPoints, inputPointDataArray, pivotValuesDT,
      threadNumber);
  else if(validationMode == 3)
    return testPointDataArray<dataType, ComperatorGreaterOrEqual<dataType>>(
      oIndex_to_mIndex_map, nPoints, inputPointDataArray, pivotValuesDT,
      threadNumber);
  else if(validationMode == 4)
    return testPointDataArray<dataType, ComperatorGreater<dataType>>(
      oIndex_to_mIndex_map, nPoints, inputPointDataArray, pivotValuesDT,
      threadNumber);

  return 0;
}

template <typename dataType>
int markVerticesBasedOnCellData(std::vector<int> &oIndex_to_mIndex_map,

                                const std::vector<double> pivotValues,
                                const size_t &nCells,
                                const vtkIdType *inConnectivityList,
                                const dataType *inputCellDataArray,
                                const size_t &validationMode,
                                const size_t &threadNumber) {
  const size_t nPivotValues = pivotValues.size();
  std::vector<dataType> pivotValuesDT(nPivotValues);
  for(size_t i = 0; i < nPivotValues; i++)
    pivotValuesDT[i] = (dataType)pivotValues[i];

  if(validationMode == 0)
    return testCellDataArray<dataType, ComperatorLess<dataType>>(
      oIndex_to_mIndex_map, nCells, inConnectivityList, inputCellDataArray,
      pivotValuesDT, threadNumber);
  else if(validationMode == 1)
    return testCellDataArray<dataType, ComperatorLessOrEqual<dataType>>(
      oIndex_to_mIndex_map, nCells, inConnectivityList, inputCellDataArray,
      pivotValuesDT, threadNumber);
  else if(validationMode == 2)
    return testCellDataArray<dataType, ComperatorEqual<dataType>>(
      oIndex_to_mIndex_map, nCells, inConnectivityList, inputCellDataArray,
      pivotValuesDT, threadNumber);
  else if(validationMode == 3)
    return testCellDataArray<dataType, ComperatorGreaterOrEqual<dataType>>(
      oIndex_to_mIndex_map, nCells, inConnectivityList, inputCellDataArray,
      pivotValuesDT, threadNumber);
  else if(validationMode == 4)
    return testCellDataArray<dataType, ComperatorGreater<dataType>>(
      oIndex_to_mIndex_map, nCells, inConnectivityList, inputCellDataArray,
      pivotValuesDT, threadNumber);

  return 0;
}

int ttkExtract::ExtractGeometry(vtkDataObject *output,
                                vtkDataObject *input,
                                const std::vector<double> &expressionValues) {
  ttk::Timer globalTimer;

  auto inputArray = this->GetInputArrayToProcess(0, input);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }
  std::string inputArrayName = inputArray->GetName();

  std::string expressionValuesString = "";
  doubleVectorToString(expressionValuesString, expressionValues);

  // print status
  const std::string CellModeS[3] = {"All", "Any", "Sub"};
  const std::string CellModeSLC[3] = {"all", "any", "sub"};
  const std::string ValidationModeS[5] = {"<", "<=", "==", ">=", ">"};
  {
    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg({{"Ext. Mode", "Geometry"},
                    {"Cell Mode", CellModeS[this->CellMode]},
                    {"Condition", "'" + inputArrayName + "' "
                                    + ValidationModeS[this->ValidationMode]
                                    + " [" + expressionValuesString + "]"}});
    this->printMsg(ttk::debug::Separator::L2);
  }

  // check input/output object validity
  auto inputAsUG = vtkUnstructuredGrid::SafeDownCast(input);
  auto outputAsUG = vtkUnstructuredGrid::SafeDownCast(output);
  if(!inputAsUG || !outputAsUG) {
    this->printErr(
      "Geometry mode requires 'vtkUnstructuredGrid' input/output.");
    return 0;
  }

  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Data array '" + inputArrayName
                   + "' must have only one component.");
    return 0;
  }

  // get points and cells
  const size_t nPoints = inputAsUG->GetNumberOfPoints();
  auto inputPD = inputAsUG->GetPointData();
  auto inputCD = inputAsUG->GetCellData();

  // Input Topo
  size_t nInCells = inputAsUG->GetNumberOfCells();
  //   vtkIdType *inConnectivityList = inputAsUG->GetCells()->GetPointer();
  vtkIdType *inConnectivityList
    = inputAsUG->GetCells()->GetData()->GetPointer(0);
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
    this->printMsg("Marking " + CellModeSLC[this->CellMode] + " cells with '"
                     + inputArrayName + "' "
                     + ValidationModeS[this->ValidationMode] + " ["
                     + expressionValuesString + "]",
                   0, ttk::debug::LineMode::REPLACE);
    ttk::Timer t;

    // Mark vertices that satisfy condition
    if(this->GetInputArrayAssociation(0, input) == 0) {
      switch(inputArray->GetDataType()) {
        vtkTemplateMacro(markVerticesBasedOnPointData<VTK_TT>(
          oIndex_to_mIndex_map, expressionValues, nPoints,
          (VTK_TT *)inputArray->GetVoidPointer(0), this->ValidationMode,
          this->threadNumber_));
      }
    } else if(this->GetInputArrayAssociation(0, input) == 1) {
      switch(inputArray->GetDataType()) {
        vtkTemplateMacro(markVerticesBasedOnCellData<VTK_TT>(
          oIndex_to_mIndex_map, expressionValues, nInCells, inConnectivityList,
          (VTK_TT *)inputArray->GetVoidPointer(0), this->ValidationMode,
          this->threadNumber_));
      }
    } else {
      this->printErr(
        "Geomerty extraction is only supported based on point and cell data.");
      return 0;
    }

    // Mark vertices based on cell mode and determine outTopo stats
    {
      if(this->CellMode == 0) {
        // All
        for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
          size_t nVertices = inConnectivityList[inTopoIndex];

          bool all = true;
          for(size_t j = 1; j <= nVertices; j++)
            if(oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]]
               == -1) {
              all = false;
              break;
            }

          if(all) {
            nOutCells++;
            markedCells[i] = true;
            for(size_t j = 1; j <= nVertices; j++)
              oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]] = -3;
            outTopologyDataSize += 1 + nVertices;
          }

          inTopoIndex += nVertices + 1;
        }

      } else if(this->CellMode == 1) {
        // Any

        // Mark Border Vertices
        for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
          size_t nVertices = inConnectivityList[inTopoIndex];
          size_t mVertices = 0;

          for(size_t j = 1; j <= nVertices; j++) {
            auto &mIndex
              = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
            if(mIndex == -2 || mIndex == -3)
              mVertices++;
          }

          // Mark border vertices
          if(mVertices > 0) {
            nOutCells++;
            markedCells[i] = true;
            for(size_t j = 1; j <= nVertices; j++) {
              auto &mIndex
                = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
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

      } else if(this->CellMode == 2) {
        // Sub
        for(size_t i = 0, inTopoIndex = 0; i < nInCells; i++) {
          size_t nVertices = inConnectivityList[inTopoIndex];
          size_t mVertices = 0;

          for(size_t j = 1; j <= nVertices; j++) {
            const auto &mIndex
              = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
            if(mIndex == -2 || mIndex == -3)
              mVertices++;
          }

          if(mVertices > 0) {
            nOutCells++;
            markedCells[i] = true;
            for(size_t j = 1; j <= nVertices; j++) {
              auto &mIndex
                = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
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

    this->printMsg("Marking " + CellModeSLC[this->CellMode] + " cells with '"
                     + inputArrayName + "' "
                     + ValidationModeS[this->ValidationMode] + " ["
                     + expressionValuesString + "]",
                   1, t.getElapsedTime());
  }

  // ---------------------------------------------------------------------
  // Extracting Points
  // ---------------------------------------------------------------------
  {
    ttk::Timer t;
    this->printMsg(
      "Extracting marked vertices", 0, ttk::debug::LineMode::REPLACE);

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
      "Extracting marked vertices (#" + std::to_string(nMarkedPoints) + ")", 1,
      t.getElapsedTime());
  }

  // -------------------------------------------------------------------------
  // Extracting Cells
  // -------------------------------------------------------------------------
  {
    ttk::Timer t;
    this->printMsg("Extracting marked cells", 0, ttk::debug::LineMode::REPLACE);

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

    if(this->CellMode == 0 || this->CellMode == 1) {
      for(size_t i = 0, inTopoIndex = 0, outTopoIndex = 0, outCellIndex = 0;
          i < nInCells; i++) {
        const size_t nVertices = inConnectivityList[inTopoIndex];
        if(markedCells[i]) {
          markedCellTypes[outCellIndex++] = types[nVertices];
          outTopologyData[outTopoIndex] = nVertices;
          for(size_t j = 1; j <= nVertices; j++)
            outTopologyData[outTopoIndex + j]
              = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
          outTopoIndex += nVertices + 1;
        }

        inTopoIndex += nVertices + 1;
      }
    } else {
      // Sub
      for(size_t i = 0, inTopoIndex = 0, outTopoIndex = 0, outCellIndex = 0;
          i < nInCells; i++) {
        const size_t nVertices = inConnectivityList[inTopoIndex];

        if(markedCells[i]) {
          size_t mVertices = 0;
          for(size_t j = 1; j <= nVertices; j++) {
            const auto &mIndex
              = oIndex_to_mIndex_map[inConnectivityList[inTopoIndex + j]];
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

    this->printMsg(
      "Extracting marked cells (#" + std::to_string(nOutCells) + ")", 1,
      t.getElapsedTime());
  }

  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (" + std::to_string(nMarkedPoints) + " vertices, "
                   + std::to_string(nOutCells) + " cells)",
                 1, globalTimer.getElapsedTime());

  outputAsUG->GetFieldData()->ShallowCopy(inputAsUG->GetFieldData());

  return 1;
}

template <class dataType>
int createUniqueValueArray(vtkDataArray *uniqueValueArray,
                           vtkDataArray *valueArray) {
  std::set<dataType> uniqueValues;

  if(uniqueValueArray->GetDataType() != valueArray->GetDataType())
    return 0;

  size_t nValues
    = valueArray->GetNumberOfTuples() * valueArray->GetNumberOfComponents();
  auto valueArrayData = (dataType *)valueArray->GetVoidPointer(0);
  for(size_t i = 0; i < nValues; i++)
    uniqueValues.insert(valueArrayData[i]);

  size_t nUniqueValues = uniqueValues.size();

  uniqueValueArray->SetNumberOfComponents(1);
  uniqueValueArray->SetNumberOfTuples(nUniqueValues);

  auto uniqueValueArrayData = (dataType *)uniqueValueArray->GetVoidPointer(0);
  auto it = uniqueValues.begin();
  for(size_t i = 0; i < nUniqueValues; i++) {
    uniqueValueArrayData[i] = *it;
    it++;
  }

  return 1;
}

int ttkExtract::ExtractArrayValues(vtkDataObject *output,
                                   vtkDataObject *input,
                                   const std::vector<double> &indices) {
  size_t nValues = indices.size();

  auto inputArray = this->GetInputArrayToProcess(0, input);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }
  std::string inputArrayName = inputArray->GetName();

  output->ShallowCopy(input);

  if(this->ExtractUniqueValues) {
    ttk::Timer t;
    this->printMsg("Extracting unique values from '" + inputArrayName + "'", 0,
                   ttk::debug::LineMode::REPLACE);

    auto uniqueValueArray
      = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
    uniqueValueArray->SetName(
      ("Unique" + std::string(inputArrayName.data())).data());

    int status = 0;
    switch(inputArray->GetDataType()) {
      vtkTemplateMacro(
        status = createUniqueValueArray<VTK_TT>(uniqueValueArray, inputArray));
    }
    if(!status) {
      this->printErr("Unable to compute unique values.");
      return 0;
    }

    output->GetFieldData()->AddArray(uniqueValueArray);

    this->printMsg("Extracting unique values from '" + inputArrayName + "'", 1,
                   t.getElapsedTime());
  } else {
    ttk::Timer t;
    std::string indicesString = "";
    doubleVectorToString(indicesString, indices);

    this->printMsg("Extracting values at [" + indicesString + "] from '"
                     + inputArrayName + "'",
                   0, ttk::debug::LineMode::REPLACE);

    auto outputArray
      = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
    outputArray->SetName(
      ("Extracted" + std::string(inputArray->GetName())).data());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(indices.size());

    for(size_t i = 0; i < nValues; i++) {
      size_t index = (size_t)indices[i];
      size_t inputArraySize = inputArray->GetNumberOfTuples();
      if(index < inputArraySize) {
        outputArray->SetTuple(i, index, inputArray);
      } else {
        this->printErr("Index out of range (" + std::to_string(i) + "/"
                       + std::to_string(inputArraySize) + ").");
        return 0;
      }
    }

    output->GetFieldData()->AddArray(outputArray);

    this->printMsg("Extracting values at [" + indicesString + "] from '"
                     + inputArrayName + "'",
                   1, t.getElapsedTime());
  }

  return 1;
}

int ttkExtract::ExtractArray(vtkDataObject *output,
                             vtkDataObject *input,
                             const std::vector<double> &indices) {

  ttk::Timer t;
  std::string indicesString;
  doubleVectorToString(indicesString, indices);
  this->printMsg("Extracting array with idx [" + indicesString + "] from "
                   + std::string(this->ArrayAttributeType == 0   ? "point"
                                 : this->ArrayAttributeType == 1 ? "cell"
                                                                 : "field")
                   + " data",
                 0, 0, ttk::debug::LineMode::REPLACE);

  output->ShallowCopy(input);

  auto outputAsDS = vtkDataSet::SafeDownCast(output);
  if(!outputAsDS) {
    this->printErr("Array extraction requires vtkDataSet input.");
    return 0;
  }

  vtkFieldData *inputAttribute
    = this->ArrayAttributeType == 0   ? outputAsDS->GetPointData()
      : this->ArrayAttributeType == 1 ? outputAsDS->GetCellData()
                                      : outputAsDS->GetFieldData();

  if(indices.size() != 1) {
    this->printErr("Array extraction can only extract exactly one array.");
    return 0;
  }

  if(indices[0] < 0 || indices[0] >= inputAttribute->GetNumberOfArrays()) {
    this->printErr("Index out of bounds.");
    return 0;
  }

  auto inputArray = inputAttribute->GetArray(indices[0]);
  auto copy = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  copy->ShallowCopy(inputArray);
  copy->SetName(this->OutputArrayName.data());

  auto outputAttribute = vtkSmartPointer<vtkFieldData>::New();
  outputAttribute->AddArray(copy);

  this->ArrayAttributeType == 0
    ? outputAsDS->GetPointData()->ShallowCopy(outputAttribute)
  : this->ArrayAttributeType == 1
    ? outputAsDS->GetCellData()->ShallowCopy(outputAttribute)
    : outputAsDS->GetFieldData()->ShallowCopy(outputAttribute);

  this->printMsg("Extracting array with indices [" + indicesString + "] from "
                   + std::string(this->ArrayAttributeType == 0   ? "point"
                                 : this->ArrayAttributeType == 1 ? "cell"
                                                                 : "field")
                   + " data",
                 1, t.getElapsedTime());

  return 1;
}

// =============================================================================
// RequestData
// =============================================================================
int ttkExtract::RequestData(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  // Get Input to Output
  auto input = vtkDataObject::GetData(inputVector[0]);
  auto output = vtkDataObject::GetData(outputVector);

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

  std::vector<double> values;
  ttkUtils::stringListToDoubleVector(finalExpressionString, values);

  auto mode = this->ExtractionMode;
  if(mode < 0) {
    if(input->IsA("vtkMultiBlockDataSet"))
      mode = 0;
    else if(input->IsA("vtkTable"))
      mode = 1;
    else {
      this->printErr("Unable to automatically determine extraction mode.");
      return 0;
    }
  }

  // in case of array or geometry extraction iterate over vtkMultiBlockDataSet
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  size_t nBlocks;
  if(mode != 0 && mode != 5 && input->IsA("vtkMultiBlockDataSet")) {
    inputAsMB->ShallowCopy(input);
    nBlocks = inputAsMB->GetNumberOfBlocks();

    for(size_t b = 0; b < nBlocks; b++) {
      auto inputBlock = inputAsMB->GetBlock(b);
      auto outputBlock
        = vtkSmartPointer<vtkDataObject>::Take(inputBlock->NewInstance());
      outputAsMB->SetBlock(b, outputBlock);
    }
    output->ShallowCopy(outputAsMB);
  } else {
    inputAsMB->SetBlock(0, input);
    outputAsMB->SetBlock(0, output);
    nBlocks = 1;
  }

  if(mode == 0) {
    if(!this->ExtractBlocks(output, input, values, false))
      return 0;
  } else if(mode == 1) {
    if(!this->ExtractRows(output, input, values))
      return 0;
  } else if(mode == 2) {
    for(size_t b = 0; b < nBlocks; b++)
      if(!this->ExtractGeometry(
           outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
        return 0;
  } else if(mode == 3) {
    for(size_t b = 0; b < nBlocks; b++)
      if(!this->ExtractArrayValues(
           outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
        return 0;
  } else if(mode == 4) {
    for(size_t b = 0; b < nBlocks; b++)
      if(!this->ExtractArray(
           outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
        return 0;
  } else if(mode == 5) {
    if(!this->ExtractBlocks(output, input, values, true))
      return 0;
  } else {
    this->printErr("Unsupported Extraction Mode: " + std::to_string(mode));
    return 0;
  }

  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
