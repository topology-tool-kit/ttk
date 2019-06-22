#include <ttkComponentSize.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkComponentSize)

  ttkComponentSize::ttkComponentSize() {
  UseAllCores = true;

  connectivityFilter_ = vtkSmartPointer<vtkConnectivityFilter>::New();
  vertexNumbers_ = vtkSmartPointer<vtkDoubleArray>::New();
  cellNumbers_ = vtkSmartPointer<vtkDoubleArray>::New();
}

ttkComponentSize::~ttkComponentSize() {
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkComponentSize::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkComponentSize::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkComponentSize] " << progress * 100 << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkComponentSize::doIt(vtkPointSet *input, vtkUnstructuredGrid *output) {

  Timer t;

  connectivityFilter_->SetInputData(input);
  connectivityFilter_->SetExtractionModeToAllRegions();
  connectivityFilter_->ColorRegionsOn();
  connectivityFilter_->Update();

  output->ShallowCopy(connectivityFilter_->GetOutput());

  vector<double> vertexNumbers(
    connectivityFilter_->GetNumberOfExtractedRegions(), 0);
  vector<double> cellNumbers(
    connectivityFilter_->GetNumberOfExtractedRegions(), 0);
  vtkDataArray *cellIds = output->GetCellData()->GetArray("RegionId");
  vtkDataArray *vertexIds = output->GetPointData()->GetArray("RegionId");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < output->GetNumberOfPoints(); i++) {

    double regionId = 0;

    vertexIds->GetTuple(i, &regionId);

    vertexNumbers[(SimplexId)regionId]++;
  }

  vertexNumbers_->SetNumberOfTuples(output->GetNumberOfPoints());
  vertexNumbers_->SetNumberOfComponents(1);
  vertexNumbers_->SetName("VertexNumber");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < output->GetNumberOfPoints(); i++) {

    double regionId = 0;
    vertexIds->GetTuple(i, &regionId);

    vertexNumbers_->SetTuple1(i, vertexNumbers[(SimplexId)regionId]);
  }
  output->GetPointData()->AddArray(vertexNumbers_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < output->GetNumberOfCells(); i++) {

    double regionId = 0;

    cellIds->GetTuple(i, &regionId);

    cellNumbers[(SimplexId)regionId]++;
  }

  cellNumbers_->SetNumberOfTuples(output->GetNumberOfCells());
  cellNumbers_->SetNumberOfComponents(1);
  cellNumbers_->SetName("CellNumber");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(SimplexId i = 0; i < output->GetNumberOfCells(); i++) {

    double regionId = 0;
    cellIds->GetTuple(i, &regionId);

    cellNumbers_->SetTuple1(i, cellNumbers[(SimplexId)regionId]);
  }
  output->GetCellData()->AddArray(cellNumbers_);

  {
    stringstream msg;
    msg << "[ttkComponentSize] Connected component sizes computed in "
        << t.getElapsedTime() << " s. (" << input->GetNumberOfPoints()
        << " points)." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkComponentSize::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  Memory m;

  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkPointSet *input = vtkPointSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkComponentSize] Memory usage: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
