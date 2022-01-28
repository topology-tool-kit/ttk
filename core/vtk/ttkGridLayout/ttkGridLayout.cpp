#include <ttkGridLayout.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkSmartPointer.h>

#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

vtkStandardNewMacro(ttkGridLayout);

ttkGridLayout::ttkGridLayout() {
  this->setDebugMsgPrefix("GridLayout");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkGridLayout::~ttkGridLayout() {
}

int ttkGridLayout::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else
    return 0;
  return 1;
}

int ttkGridLayout::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkGridLayout::CopyObject(vtkDataObject *output, vtkDataObject *input) {
  if(input->IsA("vtkImageData")) {
    output->ShallowCopy(input);
  } else if(input->IsA("vtkPointSet")) {
    auto outputAsPS = (vtkPointSet *)output;
    outputAsPS->ShallowCopy(input);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(outputAsPS->GetPoints());

    outputAsPS->SetPoints(points);
  } else {
    output->DeepCopy(input);
  }

  return 1;
}

int ttkGridLayout::TranslateObject(vtkDataObject *input,
                                   const size_t &colAxis,
                                   const size_t &rowAxis,
                                   const double &dw,
                                   const double &dh) {
  auto inputMB = vtkMultiBlockDataSet::SafeDownCast(input);
  auto inputAsPS = vtkPointSet::SafeDownCast(input);
  auto inputAsID = vtkImageData::SafeDownCast(input);

  if(inputMB) {
    size_t n = inputMB->GetNumberOfBlocks();
    for(size_t i = 0; i < n; i++)
      if(!this->TranslateObject(inputMB->GetBlock(i), colAxis, rowAxis, dw, dh))
        return 0;

  } else if(inputAsPS) {
    size_t nCoords = inputAsPS->GetNumberOfPoints() * 3;
    auto points = inputAsPS->GetPoints();
    auto pointCoords = (float *)points->GetVoidPointer(0);
    for(size_t i = 0; i < nCoords; i += 3) {
      pointCoords[i + colAxis] += dw;
      pointCoords[i + rowAxis] += dh;
    }

  } else if(inputAsID) {
    double origin[3] = {0, 0, 0};
    origin[colAxis] = dw;
    origin[rowAxis] = dh;
    inputAsID->SetOrigin(origin);
  } else {
    return 0;
  }

  return 1;
}

int ttkGridLayout::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  // Get Input and Output
  auto inputMB = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  // Determine Grid Axis
  int colAxis = this->GetColAxis();
  int rowAxis = this->GetRowAxis();

  // Compute width and height of grid cells
  double bounds[6];
  double maxWidth = 0;
  double maxHeight = 0;

  size_t nBlocks = inputMB->GetNumberOfBlocks();

  this->printMsg("Translating " + std::to_string(nBlocks) + " object(s)", 0,
                 ttk::debug::LineMode::REPLACE);

  for(size_t i = 0; i < nBlocks; i++) {
    auto block = inputMB->GetBlock(i);
    if(block->IsA("vtkMultiBlockDataSet")) {
      auto blockAsMB = vtkMultiBlockDataSet::SafeDownCast(block);
      blockAsMB->GetBounds(bounds);
    } else if(block->IsA("vtkDataSet")) {
      auto blockAsDS = vtkDataSet::SafeDownCast(block);
      blockAsDS->GetBounds(bounds);
    } else {
      this->printErr("Unable to determine bounding box of block #"
                     + std::to_string(i) + " with type "
                     + std::string(block->GetClassName()) + ".");
      return 0;
    }

    double blockWidth = bounds[colAxis * 2 + 1] - bounds[colAxis * 2];
    double blockHeight = bounds[rowAxis * 2 + 1] - bounds[rowAxis * 2];
    if(maxWidth < blockWidth)
      maxWidth = blockWidth;
    if(maxHeight < blockHeight)
      maxHeight = blockHeight;
  }

  // apply gap
  maxWidth += maxWidth * this->GetColGap() / 100.;
  maxHeight += maxHeight * this->GetRowGap() / 100.;

  // Determine Grid Structure
  const size_t nRows
    = this->GetNumberOfRows() < 1 ? 0 : (size_t)this->GetNumberOfRows();
  const size_t nColumns
    = nRows == 0 ? std::ceil(std::sqrt(nBlocks)) : std::ceil(nBlocks / nRows);

  for(size_t i = 0; i < nBlocks; i++) {
    // get block
    auto block = inputMB->GetBlock(i);

    auto outBlock = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());
    this->CopyObject(outBlock, block);

    const size_t row = std::floor(i / nColumns);
    const size_t col = i % nColumns;

    if(!this->TranslateObject(
         outBlock, colAxis, rowAxis, col * maxWidth, row * maxHeight)) {
      this->printErr("Unable to translate block #" + std::to_string(i)
                     + " of type '" + std::string(outBlock->GetClassName())
                     + "'.");
      return 0;
    }

    outputMB->SetBlock(i, outBlock);
  }

  this->printMsg("Translating " + std::to_string(nBlocks) + " object(s)", 1,
                 globalTimer.getElapsedTime());

  return 1;
}
