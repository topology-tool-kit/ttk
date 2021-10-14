#include "BaseClass.h"
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
#include <vtkSignedCharArray.h>
#include <vtkThreshold.h>

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
  if(port != 0)
    return 0;

  if(info->Has(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT()))
    info->Remove(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());

  if(this->OutputType != -1) {
    std::string DTName = this->GetVtkDataTypeName(this->OutputType);
    if(DTName.length() < 1) {
      this->printErr("Unsupported output type");
      return 0;
    }
    info->Set(vtkDataObject::DATA_TYPE_NAME(), DTName.data());
  } else
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);

  return 1;
}

// =============================================================================
// RequestInformation
// =============================================================================
int ttkExtract::RequestInformation(
  vtkInformation *,
  vtkInformationVector **ttkNotUsed(inputVector),
  vtkInformationVector *outputVector) {

  if(this->ExtractionMode == EXTRACTION_MODE::BLOCKS
     && this->GetOutputType() == VTK_IMAGE_DATA) {
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
          if(tupleIndex >= blockAsMB->GetNumberOfBlocks()) {
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
        if(blockIndex >= inputAsMB->GetNumberOfBlocks()) {
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

template <typename DT>
int computeMask_(signed char *mask,

                 const size_t &nValues,
                 const DT *values,
                 const std::vector<DT> &min,
                 const std::vector<DT> &max,
                 const size_t &threadNumber) {
  const size_t nPivotValues = min.size();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nValues; i++) {
    const DT &v = values[i];

    bool hasToBeMarked = false;
    for(size_t j = 0; j < nPivotValues; j++) {
      if(min[j] <= v && v <= max[j]) {
        hasToBeMarked = true;
        break;
      }
    }

    mask[i] = hasToBeMarked ? 1 : 0;
  }

  TTK_FORCE_USE(threadNumber);
  return 1;
}

template <typename DT>
int computeMask(signed char *mask,

                const std::vector<double> &pivotValues,
                const size_t &nValues,
                const DT *values,
                const ttkExtract::VALIDATION_MODE &validationMode,
                const size_t &threadNumber) {

  const size_t nPivotValues = pivotValues.size();
  std::vector<DT> pivotValuesMin(nPivotValues);
  std::vector<DT> pivotValuesMax(nPivotValues);

  const DT delta = std::numeric_limits<DT>::is_integer
                     ? 1
                     : std::numeric_limits<DT>::epsilon();

  for(size_t i = 0; i < nPivotValues; i++) {

    switch(validationMode) {
      case ttkExtract::VALIDATION_MODE::LESS_THEN: {
        pivotValuesMin[i] = std::numeric_limits<DT>::lowest();
        pivotValuesMax[i] = ((DT)pivotValues[i]) - delta;
        break;
      }

      case ttkExtract::VALIDATION_MODE::LESS_EQUAL_THEN: {
        pivotValuesMin[i] = std::numeric_limits<DT>::lowest();
        pivotValuesMax[i] = ((DT)pivotValues[i]);
        break;
      }

      case ttkExtract::VALIDATION_MODE::EQUAL:
      case ttkExtract::VALIDATION_MODE::UNEQUAL: {
        pivotValuesMin[i] = ((DT)pivotValues[i]);
        pivotValuesMax[i] = ((DT)pivotValues[i]);
        break;
      }

      case ttkExtract::VALIDATION_MODE::GREATER_EQUAL_THEN: {
        pivotValuesMin[i] = ((DT)pivotValues[i]);
        pivotValuesMax[i] = std::numeric_limits<DT>::max();
        break;
      }

      case ttkExtract::VALIDATION_MODE::GREATER_THEN: {
        pivotValuesMin[i] = ((DT)pivotValues[i]) + delta;
        pivotValuesMax[i] = std::numeric_limits<DT>::max();
        break;
      }
    }
  }

  int status = computeMask_<DT>(
    mask, nValues, values, pivotValuesMin, pivotValuesMax, threadNumber);

  if(validationMode == ttkExtract::VALIDATION_MODE::UNEQUAL)
    for(size_t i = 0; i < nValues; i++)
      mask[i] = mask[i] == 0 ? 1 : 0;

  return status;
}

int ttkExtract::AddMaskArray(vtkDataObject *output,
                             vtkDataObject *input,
                             const std::vector<double> &expressionValues) {
  ttk::Timer timer;

  this->printMsg(
    "Computing Mask", 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // check if input/output are of correct type
  auto inputAsDS = vtkDataSet::SafeDownCast(input);
  auto outputAsDS = vtkDataSet::SafeDownCast(output);
  if(!inputAsDS || !outputAsDS) {
    this->printErr("Masks can only be computed on vtkDataSet inputs.");
    return 0;
  }

  // retrieve input array
  auto inputArray = this->GetInputArrayToProcess(0, input);
  if(!inputArray || inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Unable to retrieve input scalar array.");
    return 0;
  }
  std::string inputArrayName = inputArray->GetName();
  const int inputArrayAssociation = this->GetInputArrayAssociation(0, input);
  if(inputArrayAssociation != 0 && inputArrayAssociation != 1) {
    this->printErr("Geometry extraction requires point or cell data.");
    return 0;
  }
  const bool isPointDataArray = this->GetInputArrayAssociation(0, input) == 0;

  // print updated status
  std::string expressionValuesString = "";
  doubleVectorToString(expressionValuesString, expressionValues);
  const std::string ValidationModeS[6] = {"<", "<=", "==", "!=", ">=", ">"};
  std::string msg = "Computing Mask: '" + inputArrayName + "' "
                    + ValidationModeS[static_cast<int>(this->ValidationMode)]
                    + " [" + expressionValuesString + "]";
  ;
  this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // initialize mask array
  const size_t nInPoints = inputAsDS->GetNumberOfPoints();
  const size_t nInCells = inputAsDS->GetNumberOfCells();
  const size_t nOutValus = isPointDataArray ? nInPoints : nInCells;

  auto maskArray = vtkSmartPointer<vtkSignedCharArray>::New();
  maskArray->SetName("Mask");
  maskArray->SetNumberOfTuples(nOutValus);
  auto maskArrayData = ttkUtils::GetPointer<signed char>(maskArray);

  // compute mask
  int status = 0;
  switch(inputArray->GetDataType()) {
    vtkTemplateMacro((
      status = computeMask<VTK_TT>(maskArrayData, expressionValues, nOutValus,
                                   ttkUtils::GetPointer<VTK_TT>(inputArray),
                                   this->ValidationMode, this->threadNumber_)));
  }
  if(!status) {
    this->printErr("Unable to compute mask");
    return 0;
  }

  // add to output
  outputAsDS->ShallowCopy(inputAsDS);
  if(isPointDataArray)
    outputAsDS->GetPointData()->AddArray(maskArray);
  else
    outputAsDS->GetCellData()->AddArray(maskArray);

  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}

int ttkExtract::ExtractGeometry(vtkDataObject *output,
                                vtkDataObject *input,
                                const std::vector<double> &expressionValues) {

  auto inputAsDS = vtkDataSet::SafeDownCast(input);
  auto outputAsDS = vtkDataSet::SafeDownCast(output);
  if(!inputAsDS || !outputAsDS) {
    this->printErr("Geometry mode requires vtkDataSet input.");
    return 0;
  }

  auto maskOutput = vtkSmartPointer<vtkDataSet>::Take(inputAsDS->NewInstance());

  if(!this->AddMaskArray(maskOutput, inputAsDS, expressionValues))
    return 0;

  if(this->MaskOnly) {
    outputAsDS->ShallowCopy(maskOutput);
  } else {
    ttk::Timer timer;
    this->printMsg(
      "Extracting Geometry based on Mask", 0, 0, ttk::debug::LineMode::REPLACE);

    auto outputAsUG = vtkUnstructuredGrid::SafeDownCast(output);
    if(!outputAsUG) {
      this->printErr(
        "Geometry Extraction requires vtkUnstructuredGrid input/output");
      return 0;
    }

    const bool isPointDataArray = this->GetInputArrayAssociation(0, input) == 0;

    auto threshold = vtkSmartPointer<vtkThreshold>::New();
    threshold->SetInputDataObject(maskOutput);
    threshold->SetInputArrayToProcess(
      0, 0, 0, isPointDataArray ? 0 : 1, "Mask");
    threshold->ThresholdByUpper(0.5);
    threshold->SetAllScalars(this->CellMode == CELL_MODE::ALL);
    threshold->Update();

    outputAsDS->ShallowCopy(threshold->GetOutput());

    if(isPointDataArray)
      outputAsDS->GetPointData()->RemoveArray("Mask");
    else
      outputAsDS->GetCellData()->RemoveArray("Mask");

    this->printMsg(
      "Extracting Geometry based on Mask", 1, timer.getElapsedTime());
  }

  return 1;
}

template <class DT>
int createUniqueValueArray(vtkDataArray *uniqueValueArray,
                           vtkDataArray *valueArray) {
  std::set<DT> uniqueValues;

  if(uniqueValueArray->GetDataType() != valueArray->GetDataType())
    return 0;

  size_t nValues
    = valueArray->GetNumberOfTuples() * valueArray->GetNumberOfComponents();
  auto valueArrayData = ttkUtils::GetPointer<DT>(valueArray);
  for(size_t i = 0; i < nValues; i++)
    uniqueValues.insert(valueArrayData[i]);

  size_t nUniqueValues = uniqueValues.size();

  uniqueValueArray->SetNumberOfComponents(1);
  uniqueValueArray->SetNumberOfTuples(nUniqueValues);

  auto uniqueValueArrayData = ttkUtils::GetPointer<DT>(uniqueValueArray);
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
int ttkExtract::RequestData(vtkInformation *ttkNotUsed(request),
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
  if(mode == EXTRACTION_MODE::AUTO) {
    if(input->IsA("vtkMultiBlockDataSet"))
      mode = EXTRACTION_MODE::BLOCKS;
    else if(input->IsA("vtkTable"))
      mode = EXTRACTION_MODE::ROWS;
    else {
      this->printErr("Unable to automatically determine extraction mode.");
      return 0;
    }
  }

  // in case of array or geometry extraction iterate over vtkMultiBlockDataSet
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  size_t nBlocks;
  if(mode != EXTRACTION_MODE::BLOCKS && mode != EXTRACTION_MODE::BLOCK_TUPLES
     && input->IsA("vtkMultiBlockDataSet")) {
    inputAsMB->ShallowCopy(input);
    nBlocks = inputAsMB->GetNumberOfBlocks();

    for(size_t b = 0; b < nBlocks; b++) {
      auto inputBlock = inputAsMB->GetBlock(b);

      vtkSmartPointer<vtkDataObject> outputBlock;
      if(mode == EXTRACTION_MODE::BLOCKS && !this->MaskOnly)
        outputBlock = vtkSmartPointer<vtkUnstructuredGrid>::New();
      else
        outputBlock
          = vtkSmartPointer<vtkDataObject>::Take(inputBlock->NewInstance());

      outputAsMB->SetBlock(b, outputBlock);
    }
    output->ShallowCopy(outputAsMB);
  } else {
    inputAsMB->SetBlock(0, input);
    outputAsMB->SetBlock(0, output);
    nBlocks = 1;
  }

  switch(mode) {
    case EXTRACTION_MODE::BLOCKS: {
      if(!this->ExtractBlocks(output, input, values, false))
        return 0;
      break;
    }
    case EXTRACTION_MODE::BLOCK_TUPLES: {
      if(!this->ExtractBlocks(output, input, values, true))
        return 0;
      break;
    }
    case EXTRACTION_MODE::ROWS: {
      if(!this->ExtractRows(output, input, values))
        return 0;
      break;
    }
    case EXTRACTION_MODE::GEOMETRY: {
      for(size_t b = 0; b < nBlocks; b++)
        if(!this->ExtractGeometry(
             outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
          return 0;
      break;
    }
    case EXTRACTION_MODE::ARRAY_VALUES: {
      for(size_t b = 0; b < nBlocks; b++)
        if(!this->ExtractArrayValues(
             outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
          return 0;
      break;
    }
    case EXTRACTION_MODE::ARRAYS: {
      for(size_t b = 0; b < nBlocks; b++)
        if(!this->ExtractArray(
             outputAsMB->GetBlock(b), inputAsMB->GetBlock(b), values))
          return 0;
      break;
    }
    default: {
      this->printErr("Unsupported Extraction Mode");
      return 0;
    }
  }

  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
