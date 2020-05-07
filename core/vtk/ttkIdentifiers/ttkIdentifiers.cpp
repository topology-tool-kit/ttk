#include <ttkIdentifiers.h>

#include <ttkMacros.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkIdentifiers);

ttkIdentifiers::ttkIdentifiers() {
  this->setDebugMsgPrefix("Identifiers");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkIdentifiers::~ttkIdentifiers() {
}

int ttkIdentifiers::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifiers::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkIdentifiers::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  ttk::Timer t;

  this->printMsg("Computing Identifiers", 0, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  auto vertexIdentifiers = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexIdentifiers->SetName(VertexFieldName.data());
  vertexIdentifiers->SetNumberOfComponents(1);
  vertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());

  auto cellIdentifiers = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  cellIdentifiers->SetName(CellFieldName.data());
  cellIdentifiers->SetNumberOfComponents(1);
  cellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());

  ttk::SimplexId vertexNumber = input->GetNumberOfPoints();
  ttk::SimplexId cellNumber = input->GetNumberOfCells();
  ttk::SimplexId count = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < vertexNumber; i++) {
    // avoid any processing if the abort signal is sent
    if((!wrapper_) || ((wrapper_) && (!wrapper_->needsToAbort()))) {

      vertexIdentifiers->SetTuple1(i, i);

      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if((wrapper_) && (!(count % ((vertexNumber) / 10)))) {
            wrapper_->updateProgress((count + 1.0) / (2 * vertexNumber));
          }

          count++;
        }
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < cellNumber; i++) {
    // avoid any processing if the abort signal is sent
    if((!wrapper_) || ((wrapper_) && (!wrapper_->needsToAbort()))) {

      cellIdentifiers->SetTuple1(i, i);

      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          if((wrapper_) && (!(count % ((cellNumber) / 10)))) {
            wrapper_->updateProgress(1 + (count + 1.0) / (2 * cellNumber));
          }

          count++;
        }
      }
    }
  }

  output->GetPointData()->AddArray(vertexIdentifiers);
  output->GetCellData()->AddArray(cellIdentifiers);

  this->printMsg(
    "Computing Identifiers", 1, t.getElapsedTime(), this->threadNumber_);

  return 1;
}