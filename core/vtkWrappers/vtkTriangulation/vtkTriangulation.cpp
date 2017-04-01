#include                  <ttkWrapper.h>
#include                  <vtkTriangulation.h>

vtkTriangulation::vtkTriangulation(){
  
  inputDataSet_ = NULL;
 
  hasAllocated_ = false;
  triangulation_ = NULL;
}

vtkTriangulation::~vtkTriangulation(){
  
  if((triangulation_)&&(hasAllocated_))
    delete triangulation_;
}

int vtkTriangulation::allocate(){

  triangulation_ = new Triangulation;
  hasAllocated_ = true;
  
  return 0;
}

int vtkTriangulation::deepCopy(vtkDataObject *other){

  if((triangulation_)&&(hasAllocated_)){
    delete triangulation_;
  }
  
  allocate();
 
  if(other->GetDataObjectType() == TTK_UNSTRUCTURED_GRID){
    // copy the triangulation object from the other
    (*triangulation_) = 
      *(((ttkUnstructuredGrid *) other)->triangulation_);
  }
  else if(other->GetDataObjectType() == TTK_IMAGE_DATA){
    // copy the triangulation object from the other
    (*triangulation_) = 
      *(((ttkImageData *) other)->triangulation_);
  }
  else if(other->GetDataObjectType() == TTK_POLY_DATA){
    // copy the triangulation object from the other
    (*triangulation_) = 
      *(((ttkPolyData *) other)->triangulation_);
  }
  
  // populate the data-structure
  if((triangulation_)&&(triangulation_->isEmpty())){
    // populate the triangulation data-structure
    setInputData((vtkDataSet *) other);
  }
  
  return 0;
}

Triangulation* vtkTriangulation::getTriangulation(vtkDataSet *other){

  if(!other)
    return NULL;

  string dataType = other->GetClassName();
  
  if((dataType == "ttkUnstructuredGrid")
    ||(dataType == "vtkUnstructuredGrid")){
    return ((ttkUnstructuredGrid *) other)->getTriangulation();
  }
  else if((dataType == "ttkImageData")
    ||(dataType == "vtkImageData")){
    return ((ttkImageData *) other)->getTriangulation();
  }
  else if((dataType == "ttkPolyData")
    ||(dataType == "vtkPolyData")){
    return ((ttkPolyData *) other)->getTriangulation();
  }
  else{
    Debug d;
    stringstream msg;
    msg << "[vtkTriangulation] Unsupported input VTK class `"
      << other->GetClassName() << "' (ref=" 
      << other->GetDataObjectType() << ")" << endl;
    d.dMsg(cerr, msg.str(), Debug::fatalMsg);
  }
  
  return NULL;
}

vtkUnstructuredGrid* vtkTriangulation::getVtkUnstructuredGrid(){
  
  if(!triangulation_)
    return NULL;
  
  if((inputDataSet_->GetDataObjectType() == VTK_IMAGE_DATA)
    ||(inputDataSet_->GetDataObjectType() == TTK_IMAGE_DATA)
    ||(inputDataSet_->GetDataObjectType() == VTK_POLY_DATA)
    ||(inputDataSet_->GetDataObjectType() == TTK_POLY_DATA)){
 
    Timer t;
  
    vtkUnstructuredGrid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkPoints_ = vtkSmartPointer<vtkPoints>::New();
  
    // here we need to create an explicit representation of the triangulation
    int vertexNumber = triangulation_->getNumberOfVertices();
    vtkPoints_->SetNumberOfPoints(vertexNumber);
    
    float p[3];
    for(int i = 0; i < vertexNumber; i++){
      triangulation_->getVertexPoint(i, p[0], p[1], p[2]);
      vtkPoints_->SetPoint(i, p[0], p[1], p[2]);
    }
    
    vtkUnstructuredGrid_->SetPoints(vtkPoints_);
    
    int cellNumber = triangulation_->getNumberOfCells();
    int dimensionality = triangulation_->getDimensionality();
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(dimensionality + 1);
    
    
    for(int i = 0; i < cellNumber; i++){
      for(int j = 0; j <= dimensionality; j++){
        int vertexId = -1;
        triangulation_->getCellVertex(i, j, vertexId);
        idList->SetId(j, vertexId);
      }
      if(dimensionality == 2){
        vtkUnstructuredGrid_->InsertNextCell(VTK_TRIANGLE, idList);
      }
      else if(dimensionality == 3){
        vtkUnstructuredGrid_->InsertNextCell(VTK_TETRA, idList);
      }
    }
    
    {
      stringstream msg;
      msg << "[vtkTriangulation] Explicit triangulation computed in "
        << t.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  else if((inputDataSet_->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    ||(inputDataSet_->GetDataObjectType() == TTK_UNSTRUCTURED_GRID)){
    return vtkUnstructuredGrid::SafeDownCast(inputDataSet_);
  }
  else{
    stringstream msg;
    msg << "[vtkTriangulation] Unsupported input VTK class `"
      << inputDataSet_->GetClassName() << "' (ref=" 
      << inputDataSet_->GetDataObjectType() << ")" << endl;
    msg << "[vtkTriangulation] Leaving an empty triangulation..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
  }
  
  return vtkUnstructuredGrid_.GetPointer();
}


bool vtkTriangulation::hasChangedConnectivity(
  Triangulation *triangulation,
  vtkDataSet *dataSet, vtkObject *callingObject)
{

#ifndef withKamikaze
  if(!triangulation)
    return false;
  if(!dataSet)
    return false;
  if(!callingObject)
    return false;
#endif
    
  if((dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    ||(dataSet->GetDataObjectType() == TTK_UNSTRUCTURED_GRID)){
    return (((vtkUnstructuredGrid *) dataSet)->GetCells()->GetMTime() 
      > callingObject->GetMTime());
  }
  else if((dataSet->GetDataObjectType() == VTK_POLY_DATA)
    ||(dataSet->GetDataObjectType() == TTK_POLY_DATA)){
    return (((vtkPolyData *) dataSet)->GetPolys()->GetMTime() 
      > callingObject->GetMTime());
  }
  else if((dataSet->GetDataObjectType() == VTK_IMAGE_DATA)
    ||(dataSet->GetDataObjectType() == TTK_IMAGE_DATA)){
    
    int vtkDimensions[3];
    vector<int> ttkDimensions;
    
    ((vtkImageData *) dataSet)->GetDimensions(vtkDimensions);
    triangulation->getGridDimensions(ttkDimensions);
    
    return !((vtkDimensions[0] == ttkDimensions[0])
      &&(vtkDimensions[1] == ttkDimensions[1])
      &&(vtkDimensions[2] == ttkDimensions[2]));
  }
  else{
    stringstream msg;
    Debug d;
    msg << "[vtkTriangulation] Unsupported input VTK class `"
      << dataSet->GetClassName() << "' (ref=" 
      << dataSet->GetDataObjectType() << ")" << endl;
    d.dMsg(cerr, msg.str(), Debug::fatalMsg);
  }
  
  return true;
}

int vtkTriangulation::setInputData(vtkDataSet* dataSet){

  if(!triangulation_){
    allocate();
  }
  
  if((dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    ||(dataSet->GetDataObjectType() == TTK_UNSTRUCTURED_GRID)){
    
    if(((vtkUnstructuredGrid *) dataSet)->GetPoints()){
      triangulation_->setInputPoints(dataSet->GetNumberOfPoints(),
        (float *)
        ((vtkUnstructuredGrid *) dataSet)->GetPoints()->GetVoidPointer(0));
    }
    if(((vtkUnstructuredGrid *) dataSet)->GetCells()){
      triangulation_->setInputCells(dataSet->GetNumberOfCells(),
        ((vtkUnstructuredGrid *) dataSet)->GetCells()->GetPointer());
    }
    inputDataSet_ = dataSet;
  }
  else if((dataSet->GetDataObjectType() == VTK_POLY_DATA)
    ||(dataSet->GetDataObjectType() == TTK_POLY_DATA)){
    
    if(((vtkPolyData *) dataSet)->GetPoints()){
      triangulation_->setInputPoints(dataSet->GetNumberOfPoints(),
        (float *)
        ((vtkPolyData *) dataSet)->GetPoints()->GetVoidPointer(0));
    }
    if(((vtkPolyData *) dataSet)->GetPolys()){
      triangulation_->setInputCells(dataSet->GetNumberOfCells(),
        ((vtkPolyData *) dataSet)->GetPolys()->GetPointer());
    }
    inputDataSet_ = dataSet;
  }
  else if((dataSet->GetDataObjectType() == VTK_IMAGE_DATA)
    ||(dataSet->GetDataObjectType() == TTK_IMAGE_DATA)){
    vtkImageData *imageData = (vtkImageData *) dataSet;

    int extents[6];
    imageData->GetExtent(extents);

    double origin[3];
    imageData->GetOrigin(origin);

    double spacing[3];
    imageData->GetSpacing(spacing);

    int gridDimensions[3];
    imageData->GetDimensions(gridDimensions);

    double firstPoint[3];
    firstPoint[0]=origin[0]+extents[0]*spacing[0];
    firstPoint[1]=origin[1]+extents[2]*spacing[1];
    firstPoint[2]=origin[2]+extents[4]*spacing[2];

    triangulation_->setInputGrid(
      firstPoint[0], firstPoint[1], firstPoint[2],
      spacing[0], spacing[1], spacing[2],
      gridDimensions[0], gridDimensions[1], gridDimensions[2]);
    
    inputDataSet_ = dataSet;
  }
  else{
    stringstream msg;
    msg << "[vtkTriangulation] Unsupported input VTK class `"
      << dataSet->GetClassName() << "' (ref=" 
      << dataSet->GetDataObjectType() << ")" << endl;
    msg << "[vtkTriangulation] Leaving an empty triangulation..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
  }
  
  return 0;
}

int vtkTriangulation::shallowCopy(vtkDataObject *other){

  if((triangulation_)&&(hasAllocated_)){
    delete triangulation_;
  }
  triangulation_ = NULL;
  
  if((other)&&(((vtkDataSet *) other)->GetNumberOfPoints())){
    if(other->GetDataObjectType() == TTK_UNSTRUCTURED_GRID){
      triangulation_ = 
        ((ttkUnstructuredGrid *) other)->triangulation_;
      hasAllocated_ = false;
    }
    else if(other->GetDataObjectType() == TTK_IMAGE_DATA){
      triangulation_ = 
        ((ttkImageData *) other)->triangulation_;
      hasAllocated_ = false;
    }
    else if(other->GetDataObjectType() == TTK_POLY_DATA){
      triangulation_ = 
        ((ttkPolyData *) other)->triangulation_;
      hasAllocated_ = false;
    }
    else{
      // let's create the object
      allocate();
    }

    // populate the data-structure
    if((triangulation_)&&(triangulation_->isEmpty())){
      setInputData((vtkDataSet *) other);
    }
  }
  
  return 0;
}

vtkStandardNewMacro(ttkUnstructuredGrid)

ttkUnstructuredGrid::ttkUnstructuredGrid(){
}

ttkUnstructuredGrid::~ttkUnstructuredGrid(){
}

void ttkUnstructuredGrid::CopyStructure(vtkDataSet *other){
  
  vtkUnstructuredGrid::CopyStructure(other);
  vtkTriangulation::shallowCopy(other);
}


void ttkUnstructuredGrid::DeepCopy(vtkDataObject *other){
  
  vtkUnstructuredGrid::DeepCopy(other);
  vtkTriangulation::deepCopy(other);
}

void ttkUnstructuredGrid::ShallowCopy(vtkDataObject *other){

  vtkUnstructuredGrid::ShallowCopy(other);
  vtkTriangulation::shallowCopy(other );
}

vtkStandardNewMacro(ttkImageData)

ttkImageData::ttkImageData(){
}

ttkImageData::~ttkImageData(){
}

void ttkImageData::CopyStructure(vtkDataSet *other){
  
  vtkImageData::CopyStructure(other);
  vtkTriangulation::shallowCopy(other);
}


void ttkImageData::DeepCopy(vtkDataObject *other){
  
  vtkImageData::DeepCopy(other);
  vtkTriangulation::deepCopy(other);
}

void ttkImageData::ShallowCopy(vtkDataObject *other){

  vtkImageData::ShallowCopy(other);
  vtkTriangulation::shallowCopy(other );
}

vtkStandardNewMacro(ttkPolyData)

ttkPolyData::ttkPolyData(){
}

ttkPolyData::~ttkPolyData(){
}

void ttkPolyData::CopyStructure(vtkDataSet *other){
  
  vtkPolyData::CopyStructure(other);
  vtkTriangulation::shallowCopy(other);
}

void ttkPolyData::DeepCopy(vtkDataObject *other){
  
  vtkPolyData::DeepCopy(other);
  vtkTriangulation::deepCopy(other);
}

void ttkPolyData::ShallowCopy(vtkDataObject *other){

  vtkPolyData::ShallowCopy(other);
  vtkTriangulation::shallowCopy(other);
}
