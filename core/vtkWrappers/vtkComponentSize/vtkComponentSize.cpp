#include                  <vtkComponentSize.h>

vtkStandardNewMacro(vtkComponentSize)

vtkComponentSize::vtkComponentSize(){

  connectivityFilter_ = vtkSmartPointer<vtkConnectivityFilter>::New();
  vertexNumbers_ = vtkSmartPointer<vtkDoubleArray>::New();
  cellNumbers_ = vtkSmartPointer<vtkDoubleArray>::New();
}

vtkComponentSize::~vtkComponentSize(){

}


// transmit abort signals -- to copy paste in other wrappers
bool vtkComponentSize::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int vtkComponentSize::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[vtkComponentSize] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

int vtkComponentSize::doIt(vtkPointSet *input, vtkUnstructuredGrid *output){

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
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < output->GetNumberOfPoints(); i++){
    
    double regionId = 0;
    
    vertexIds->GetTuple(i, &regionId);
    
    vertexNumbers[(int) regionId]++;
  }
  
  vertexNumbers_->SetNumberOfTuples(output->GetNumberOfPoints());
  vertexNumbers_->SetNumberOfComponents(1);
  vertexNumbers_->SetName("VertexNumber");

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
  for(int i = 0; i < output->GetNumberOfPoints(); i++){
    
    double regionId = 0;
    vertexIds->GetTuple(i, &regionId);
    
    vertexNumbers_->SetTuple1(i, vertexNumbers[(int) regionId]);
  }
#endif
  output->GetPointData()->AddArray(vertexNumbers_);
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < output->GetNumberOfCells(); i++){
    
    double regionId = 0;
    
    cellIds->GetTuple(i, &regionId);
    
    cellNumbers[(int) regionId]++;
  }
  
  cellNumbers_->SetNumberOfTuples(output->GetNumberOfCells());
  cellNumbers_->SetNumberOfComponents(1);
  cellNumbers_->SetName("CellNumber");

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
  for(int i = 0; i < output->GetNumberOfCells(); i++){
    
    double regionId = 0;
    cellIds->GetTuple(i, &regionId);
    
    cellNumbers_->SetTuple1(i, cellNumbers[(int) regionId]);
  }
#endif
  output->GetCellData()->AddArray(cellNumbers_);
  
  {
    stringstream msg;
    msg << "[vtkComponentSize] Connected component sizes computed in " 
      << t.getElapsedTime() << " s. (" 
      << input->GetNumberOfPoints() << " points)." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int vtkComponentSize::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkPointSet *input = vtkPointSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[vtkComponentSize] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 1;
}