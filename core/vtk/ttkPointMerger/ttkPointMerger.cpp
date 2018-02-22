#include                  <ttkPointMerger.h>

vtkStandardNewMacro(ttkPointMerger)

int ttkPointMerger::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Timer t;
  Memory m;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
  
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  if(BoundaryOnly)
    triangulation->preprocessBoundaryVertices();
  

  
  int vertexNumber = input->GetNumberOfPoints();
  vector<int> candidateVertices;
  
  if(BoundaryOnly){
    for(int i = 0; i < vertexNumber; i++){
      if(triangulation->isVertexOnBoundary(i)){
        candidateVertices.push_back(i);
      }
    }
  }
  else{
    candidateVertices.resize(vertexNumber);
    for(int i = 0; i < vertexNumber; i++)
      candidateVertices[i] = i;
  }
  
  vector<vector<int> > closePoints(vertexNumber);
  
  {
    stringstream msg;
    msg << "[ttkPointMerger] Computing pointwise distances ("
      << candidateVertices.size() << " candidates)..." << endl;
    dMsg(cout, msg.str(), Debug::timeMsg);
  }
  
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) candidateVertices.size(); i++){
    double distance = -1;
    vector<double> p0(3);
    
    input->GetPoint(candidateVertices[i], p0.data());
    for(int j = 0; j < (int) candidateVertices.size(); j++){
      if(i != j){
        vector<double> p1(3);
        input->GetPoint(candidateVertices[j], p1.data());
        distance = Geometry::distance(p0.data(), p1.data());
        if(distance < DistanceThreshold){
          closePoints[candidateVertices[i]].push_back(candidateVertices[j]);
        }
      }
    }
  }
  
  vector<double> minMergeDistance(vertexNumber, -1);
  vector<double> maxMergeDistance(vertexNumber, -1);
  vector<int> mergeCount(vertexNumber, 0);
  vector<int> mergeMap(vertexNumber, -1);
  for(int i = 0; i < vertexNumber; i++){
    if(closePoints[i].size()){
      
      int minIdentifier = i;
      do{
        for(int j = 0; j < (int) closePoints[minIdentifier].size(); j++){
          if(closePoints[minIdentifier][j] < minIdentifier)
            minIdentifier = closePoints[minIdentifier][j];
        }
        mergeMap[i] = minIdentifier;
      }while(mergeMap[minIdentifier] != minIdentifier);
      mergeCount[minIdentifier]++;
      vector<double> p0(3), p1(3);
      input->GetPoint(i, p0.data());
      input->GetPoint(minIdentifier, p1.data());
      double distance = Geometry::distance(p0.data(), p1.data());
      if((minMergeDistance[minIdentifier] == -1)
        ||(distance < minMergeDistance[minIdentifier])){
        minMergeDistance[minIdentifier] = distance;
      }
      if((maxMergeDistance[minIdentifier] == -1)
        ||(distance > maxMergeDistance[minIdentifier])){
        maxMergeDistance[minIdentifier] = distance;
        }
    }
    else{
      mergeMap[i] = i;
    }
  }

  int vertexIdGen = 0;
  vector<int> old2new(vertexNumber, -1);
  for(int i = 0; i < vertexNumber; i++){
    if(mergeMap[i] == i){
      old2new[i] = vertexIdGen;
      vertexIdGen++;
    }
  }
  
  // now create the output
  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cellArray = 
    vtkSmartPointer<vtkCellArray>::New();
    
  vector<vtkSmartPointer<vtkDoubleArray> > pointData;
  pointData.resize(input->GetPointData()->GetNumberOfArrays());
 
  {
    stringstream msg;
    msg << "[ttkPointMerger] Merge performed in "
      << t.getElapsedTime() << " s." << endl;
    msg << "[ttkPointMerger] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
