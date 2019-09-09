#include <ttkIdentifiers.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifiers)

  ttkIdentifiers::ttkIdentifiers() {

  // init
  CellFieldName = "CellIdentifiers";
  VertexFieldName = ttk::VertexScalarFieldName;
  UseAllCores = true;
}

ttkIdentifiers::~ttkIdentifiers() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkIdentifiers::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkIdentifiers::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkIdentifiers] " << progress * 100 << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkIdentifiers::doIt(vtkDataSet *input, vtkDataSet *output) {

  Timer t;

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vtkSmartPointer<ttkSimplexIdTypeArray> cellIdentifiers
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();

  vertexIdentifiers->SetName(VertexFieldName.data());
  vertexIdentifiers->SetNumberOfComponents(1);
  vertexIdentifiers->SetNumberOfTuples(input->GetNumberOfPoints());

  cellIdentifiers->SetName(CellFieldName.data());
  cellIdentifiers->SetNumberOfComponents(1);
  cellIdentifiers->SetNumberOfTuples(input->GetNumberOfCells());

  SimplexId vertexNumber = input->GetNumberOfPoints();
  SimplexId cellNumber = input->GetNumberOfCells();
  SimplexId count = 0;

  //   // see also vtkOriginalCellIds
  //   vtkDataArray *original =
  //     input->GetPointData()->GetArray("vtkOriginalPointIds");
  //   printf("original: %d\n", original);
  //   if(original){
  //     printf("\t%d entries...\n", original->GetNumberOfTuples());
  //   }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < vertexNumber; i++) {
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
  for(SimplexId i = 0; i < cellNumber; i++) {
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

  {
    stringstream msg;
    msg << "[ttkIdentifiers] Identifiers generated in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkIdentifiers::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkIdentifiers] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
