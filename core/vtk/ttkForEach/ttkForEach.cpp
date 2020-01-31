#include <ttkForEach.h>

#include <vtkObjectFactory.h> // for new macro

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkForEach);

ttkForEach::ttkForEach() {
  this->setDebugMsgPrefix("ForEach");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkForEach::~ttkForEach(){};

int ttkForEach::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  else
    return 0;
  return 1;
}

int ttkForEach::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkForEach::RequestInformation(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // These values need to exists to automatically enable temporal streaming
  double dummy[2] = {0, 1};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), dummy, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), dummy, 2);

  return this->Superclass::RequestInformation(
    request, inputVector, outputVector);
}

int ttkForEach::RequestData(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  // Get current row index
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  double iIteration
    = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  // Get Input and Output
  auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());
  auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
  auto inputAsT = vtkTable::SafeDownCast(input);
  auto inputAsUG = vtkUnstructuredGrid::SafeDownCast(input);

  if(!inputAsMB && !inputAsT && !inputAsUG) {
    this->printErr("Unsupported input type");
    return 0;
  }

  auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
  auto outputAsMB = vtkMultiBlockDataSet::SafeDownCast(output);
  auto outputAsT = vtkTable::SafeDownCast(output);
  auto outputAsUG = vtkUnstructuredGrid::SafeDownCast(output);

  // Determine Mode
  std::string modes[6] = {"", "B", "R", "P", "C", "F"};
  int mode = this->GetMode();
  if(mode == 0) { // Automatic
    if(inputAsMB)
      mode = 1; // For each block
    else if(inputAsT)
      mode = 2; // For each row
    else
      mode = 5; // For each field data element
  }

  // Iteration info
  auto iterationInformation = vtkSmartPointer<vtkDoubleArray>::New();
  iterationInformation->SetName("_ttk_IterationInfo");
  iterationInformation->SetNumberOfValues(2);
  iterationInformation->SetValue(0, iIteration);

  switch(mode) {
    case 1: { // Block
      if(!inputAsMB || !outputAsMB) {
        this->printErr("Block iteration requires vtkMultiBlockDataSet input");
        return 0;
      }
      iterationInformation->SetValue(1, inputAsMB->GetNumberOfBlocks());

      outputAsMB->SetBlock(0, inputAsMB->GetBlock(iIteration));
      input->GetFieldData()->ShallowCopy(output->GetFieldData());
      outputAsMB->GetFieldData()->AddArray(iterationInformation);
      break;
    }
    case 2: { // Row
      if(!inputAsT || !outputAsT) {
        this->printErr("Row iteration requires vtkTable input");
        return 0;
      }
      iterationInformation->SetValue(1, inputAsT->GetNumberOfRows());

      // Extract row at index
      size_t n = inputAsT->GetNumberOfColumns();
      for(size_t i = 0; i < n; i++) {
        auto iColumn = inputAsT->GetColumn(i);
        auto oColumn
          = vtkSmartPointer<vtkAbstractArray>::Take(iColumn->NewInstance());
        oColumn->SetName(iColumn->GetName());
        oColumn->SetNumberOfTuples(1);
        oColumn->SetTuple(0, iIteration, iColumn);
        outputAsT->AddColumn(oColumn);
      }
      input->GetFieldData()->ShallowCopy(output->GetFieldData());
      outputAsT->GetFieldData()->AddArray(iterationInformation);
      break;
    }
    case 3: { // Point
      if(!inputAsUG || !outputAsUG) {
        this->printErr("Point iteration requires 'vtkUnstructuredGrid'");
        return 0;
      }
      iterationInformation->SetValue(1, inputAsUG->GetNumberOfPoints());

      this->printErr("Point iteration is currently not supported");
      return 0;

      break;
    }
    case 4: { // Cell
      if(!inputAsUG || !outputAsUG) {
        this->printErr("Cell iteration requires 'vtkUnstructuredGrid'");
        return 0;
      }
      iterationInformation->SetValue(1, inputAsUG->GetNumberOfCells());

      this->printErr("Cell iteration is currently not supported");
      return 0;

      break;
    }
    case 5: { // Field
      output->ShallowCopy(input);
      auto fd = output->GetFieldData()->GetAbstractArray(
        this->GetFieldDataName().data());
      if(!fd) {
        this->printErr("Field data '" + this->GetFieldDataName()
                       + "' not found");
        return 0;
      }
      size_t nTuples = fd->GetNumberOfTuples();
      size_t nComponents = fd->GetNumberOfComponents();
      iterationInformation->SetValue(1, nComponents * nTuples);

      output->GetFieldData()->AddArray(iterationInformation);
      break;
    }
    default:
      this->printErr("Unsupported Mode");
      return 0;
  }

  // Print status
  this->printMsg("[" + modes[mode] + "] Iteration: ( "
                   + std::to_string((int)(iIteration + 1)) + " / "
                   + std::to_string((int)(iterationInformation->GetValue(1)))
                   + " ) ",
                 ttk::debug::Separator::SLASH);

  return 1;
}