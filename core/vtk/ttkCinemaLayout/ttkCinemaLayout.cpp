#include <ttkCinemaLayout.h>

#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaLayout)

  int ttkCinemaLayout::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  Timer t;
  Memory m;

  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkCinemaLayout] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Get Input and Output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto inputMB = vtkMultiBlockDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto outputMB = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Determine Grid Size
  size_t nBlocks = inputMB->GetNumberOfBlocks();

  int columnAxis = this->GetColumnAxis();
  int rowAxis = this->GetRowAxis();

  // Iterate over blocks and arrange images on the
  size_t nRows
    = this->GetNumberOfRows() < 1 ? 0 : (size_t)this->GetNumberOfRows();
  size_t nColumns = nRows == 0 ? ceil(sqrt(nBlocks))
                               : (nBlocks % nRows) == 0 ? nBlocks / nRows
                                                        : ceil(nBlocks / nRows);

  double width = 0;
  double height = 0;

  for(size_t i = 0; i < nBlocks; i++) {
    auto block = vtkDataSet::SafeDownCast(inputMB->GetBlock(i));
    if(block == nullptr) {
      stringstream msg;
      msg << "[ttkCinemaLayout] ERROR: This filter can only arrange vtkDataSet "
             "objects"
          << endl;
      dMsg(cout, msg.str(), fatalMsg);
      return 0;
    }
    double bounds[6];
    block->GetBounds(bounds);
    width = max(width, bounds[2 * columnAxis + 1] - bounds[2 * columnAxis]);
    height = max(height, bounds[2 * rowAxis + 1] - bounds[2 * rowAxis]);
  }

  width += width * (this->GetColumnGap() / 100);
  height += height * (this->GetRowGap() / 100);

  size_t i = 0;
  double y = 0;
  while(i < nBlocks) {
    for(size_t x = 0; x < nColumns && i < nBlocks; x++) {

      double translation[3] = {0, 0, 0};
      translation[rowAxis] = y;
      translation[columnAxis] = x * width;

      auto transform = vtkSmartPointer<vtkTransform>::New();
      transform->Translate(translation);

      auto block = vtkDataSet::SafeDownCast(inputMB->GetBlock(i));

      auto transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
      transformFilter->SetTransform(transform);
      transformFilter->SetInputData(block);
      transformFilter->Update();

      outputMB->SetBlock(i, transformFilter->GetOutput());

      i++;
    }
    y += height;
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkCinemaLayout] "
           "--------------------------------------------------------------"
        << endl;
    msg << "[ttkCinemaLayout]   time: " << t.getElapsedTime() << " s." << endl;
    msg << "[ttkCinemaLayout] memory: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
