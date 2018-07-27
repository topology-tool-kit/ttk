#include                  <ttkReebSpace.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkReebSpace)

ttkReebSpace::ttkReebSpace(){

  // init
  SetNumberOfOutputPorts(4);
  
  ZeroSheetId = true;
  ZeroSheetType = true;
  ZeroSheetValue = true;
  ZeroSheetVertexId = true;
  
  OneSheetId = true;
  OneSheetType = true;
  OneSheetValue = true;
  OneSheetVertexId = true;
  OneSheetEdgeId = true;
  
  TwoSheets = true;
  TwoSheetCaseId = true;
  TwoSheetEdgeId = true;
  TwoSheetEdgeType = true;
  TwoSheetId = true;
  TwoSheetParameterization = true;
  TwoSheetTetId = true;
  TwoSheetValue = true;
  
  ThreeSheetTetNumber = true;
  ThreeSheetVertexNumber = true;
  ThreeSheetExpansion = true;
  ThreeSheetDomainVolume = true;
  ThreeSheetRangeArea = true;
  ThreeSheetHyperVolume = true;
  
  SimplificationThreshold = 0;
  SimplificationCriterion = 1;
  
  UseAllCores = true;
  ThreadNumber = 1;
  debugLevel_ = 3;
  ForceInputOffsetScalarField = false;
  
  uComponent_ = NULL;
  vComponent_ = NULL;
  offsetFieldU_ = NULL;
  offsetFieldV_ = NULL;

  UcomponentId = 0;
  VcomponentId = 1;

  UseOctreeAcceleration = true;
}

ttkReebSpace::~ttkReebSpace(){

}


template<class dataTypeU, class dataTypeV> int ttkReebSpace::baseCall(
  vtkDataSet *input, 
  vtkDataArray *uField, vtkDataArray *offsetFieldU,
  vtkDataArray *vField, vtkDataArray *offsetFieldV){

  Timer t;
  
  reebSpace_.setWrapper(this);
  
  bool VaryingValues = 
    (uField->GetMTime() > GetMTime())
    ||((offsetFieldU)&&(offsetFieldU->GetMTime() > GetMTime()))
    ||(vField->GetMTime() > GetMTime())
    ||((offsetFieldV)&&(offsetFieldV->GetMTime() > GetMTime()))
    ||(reebSpace_.setRangeDrivenOctree(UseOctreeAcceleration));
 
  // first time or the values changed
  reebSpace_.setInputField(uField->GetVoidPointer(0), 
    vField->GetVoidPointer(0));
  reebSpace_.setSosOffsetsU(&sosOffsetsU_);
  reebSpace_.setSosOffsetsV(&sosOffsetsV_);
    
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
  
  if(!triangulation)
    return -1;
  
  bool VaryingTriangulation = false;
  if(triangulation->isEmpty()) 
    VaryingTriangulation = true;
  if(ttkTriangulation::hasChangedConnectivity(triangulation, input, this))
    VaryingTriangulation = true;
  
  triangulation->setWrapper(this);
  reebSpace_.setupTriangulation<dataTypeU, dataTypeV>(triangulation);
 
  {
    stringstream msg;
    msg << "[ttkReebSpace] U-component: '"
      << uField->GetName() << "'" << endl;
    msg << "[ttkReebSpace] V-component: '"
      << vField->GetName() << "'" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // go!
  if((reebSpace_.empty())||(VaryingValues)||(VaryingTriangulation)){
    {
      stringstream msg;
      msg << "[ttkReebSpace] Starting computation..." << endl;
      
      dMsg(cout, msg.str(), timeMsg);
    }
    reebSpace_.execute<dataTypeU, dataTypeV>();
  }
  
  if(SimplificationThreshold > 0){
    reebSpace_.simplify<dataTypeU, dataTypeV>(
      SimplificationThreshold, 
      (ReebSpace::SimplificationCriterion) SimplificationCriterion);
  }
 
  Modified();
 
  return 0;
}
#ifdef _MSC_VER
#define COMMA ,
#endif 
int ttkReebSpace::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){

  Memory m;
  
  if(inputs.size() != 1)
    return -1;
  
  if(outputs.size() != 4)
    return -2;
  
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *sheet0 = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *sheet1 = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
  vtkUnstructuredGrid *sheet2 = vtkUnstructuredGrid::SafeDownCast(outputs[2]);
  vtkDataSet *sheet3 = outputs[3];
  
  // check data components
  if(Ucomponent.length()){
    string oldName = Ucomponent;
    if(uComponent_){
      oldName = uComponent_->GetName();
    }
    uComponent_ = input->GetPointData()->GetArray(Ucomponent.data());
    if(oldName != Ucomponent)
      uComponent_->Modified();
  }
  else{
    // default
    uComponent_ = input->GetPointData()->GetArray(UcomponentId);
  }
  if(!uComponent_)
    return -1;
  
  if(Vcomponent.length()){
    string oldName = Vcomponent;
    if(vComponent_){
      oldName = vComponent_->GetName();
    }
    vComponent_ = input->GetPointData()->GetArray(Vcomponent.data());
    if(oldName != Vcomponent)
      vComponent_->Modified();
  }
  else{
    // default
    vComponent_ = input->GetPointData()->GetArray(VcomponentId);
  }
  if(!vComponent_)
    return -2;
 
  if(ForceInputOffsetScalarField){
    if(OffsetFieldU.length()){
      
      string oldName = OffsetFieldU;
      if(offsetFieldU_){
        oldName = offsetFieldU_->GetName();
      }
      offsetFieldU_ = input->GetPointData()->GetArray(OffsetFieldU.data());
      if(oldName != OffsetFieldU)
        offsetFieldU_->Modified();
      
      if(offsetFieldU_){
        sosOffsetsU_.resize(offsetFieldU_->GetNumberOfTuples());
        for(ttkIdType i = 0; i < offsetFieldU_->GetNumberOfTuples(); i++){
          sosOffsetsU_[i] = offsetFieldU_->GetTuple1(i);
        }
      }
    }
    if(OffsetFieldV.length()){
      
      string oldName = OffsetFieldV;
      if(offsetFieldV_){
        oldName = offsetFieldV_->GetName();
      }
      offsetFieldV_ = input->GetPointData()->GetArray(OffsetFieldV.data());
      if(oldName != OffsetFieldV)
        offsetFieldV_->Modified();
      
      if(offsetFieldV_){
        sosOffsetsV_.resize(offsetFieldV_->GetNumberOfTuples());
        for(ttkIdType i = 0; i < offsetFieldV_->GetNumberOfTuples(); i++){
          sosOffsetsV_[i] = offsetFieldV_->GetTuple1(i);
        }
      }
    }
  }
  else{
    if(input->GetPointData()->GetArray(ttk::OffsetFieldUName)){
      offsetFieldU_ = input->GetPointData()->GetArray(ttk::OffsetFieldUName);
      if(offsetFieldU_){
        sosOffsetsU_.resize(offsetFieldU_->GetNumberOfTuples());
        for(ttkIdType i = 0; i < offsetFieldU_->GetNumberOfTuples(); i++){
          sosOffsetsU_[i] = offsetFieldU_->GetTuple1(i);
        }
      }
    }
    if(input->GetPointData()->GetArray(ttk::OffsetFieldVName)){
      offsetFieldV_ = input->GetPointData()->GetArray(ttk::OffsetFieldVName);
      if(offsetFieldV_){
        sosOffsetsV_.resize(offsetFieldV_->GetNumberOfTuples());
        for(ttkIdType i = 0; i < offsetFieldV_->GetNumberOfTuples(); i++){
          sosOffsetsV_[i] = offsetFieldV_->GetTuple1(i);
        }
      }
    }
  }
  
  // set the Reeb space functor
  switch(uComponent_->GetDataType()){
    
    case VTK_CHAR:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<char, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<char COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    case VTK_DOUBLE:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<double, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<double COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    case VTK_FLOAT:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<float, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<float COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    case VTK_INT:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<int, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<int COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;

    case VTK_ID_TYPE:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<vtkIdType, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<vtkIdType COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    case VTK_UNSIGNED_CHAR:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<unsigned char, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<unsigned char COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    case VTK_UNSIGNED_SHORT:
      switch(vComponent_->GetDataType()){
#ifndef _MSC_VER
		  vtkTemplateMacro((
		  {
			  baseCall<unsigned short, VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  ));
#else
		  vtkTemplateMacro(
		  {
			  baseCall<unsigned short COMMA VTK_TT>(input,
			  uComponent_, offsetFieldU_,
			  vComponent_, offsetFieldV_);
		  }
		  );
#endif
      }
      break;
      
    default:
      {
        stringstream msg;
        msg << "[ttkReebSpace] Unsupported U-component data type :( ["
          << uComponent_->GetDataType() << "]" << endl;
        dMsg(cerr, msg.str(), 1);
      }
      break;
  }
  
  // prepare the output
  {
    stringstream msg;
    msg << "[ttkReebSpace] Preparing the VTK-output..." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
  
  if(!triangulation)
    return -3;
  
  // 0-sheets - 
  // Optional additional fields:
  // PointData; u, v, vertexId, type, sheetId
  const vector<ttkIdType> *sheet0segmentation = reebSpace_.get0sheetSegmentation();
  ttkIdType vertexNumber = 0;
  for(ttkIdType i = 0; i < (ttkIdType) sheet0segmentation->size(); i++){
    ttkIdType sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1){
      const ReebSpace::Sheet0 *sheet = reebSpace_.get0sheet(sheet0Id);
      if(!sheet->pruned_){
        vertexNumber++;
      }
    }
  }
  vtkSmartPointer<vtkPoints> sheet0Points = 
    vtkSmartPointer<vtkPoints>::New();
  
  if(!sheet0->GetPoints())
    sheet0->SetPoints(sheet0Points);
    
  sheet0->GetPoints()->SetNumberOfPoints(vertexNumber);
  
  vtkSmartPointer<vtkDoubleArray> vertexScalarsU = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> vertexScalarsV = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ZeroSheetValue){
    vertexScalarsU->SetNumberOfTuples(vertexNumber);
    vertexScalarsU->SetName(uComponent_->GetName());
    vertexScalarsV->SetNumberOfTuples(vertexNumber);
    vertexScalarsV->SetName(vComponent_->GetName());
  }
  else{
    sheet0->GetPointData()->RemoveArray(uComponent_->GetName());
    sheet0->GetPointData()->RemoveArray(vComponent_->GetName());
  }
  
  vtkSmartPointer<ttkIdTypeArray> vertexIds = 
    vtkSmartPointer<ttkIdTypeArray>::New();
  if(ZeroSheetVertexId){
    vertexIds->SetNumberOfTuples(vertexNumber);
    vertexIds->SetName("VertexIds");
  }
  else{
    sheet0->GetPointData()->RemoveArray("VertexIds");
  }
  
  vtkSmartPointer<vtkCharArray> vertexTypes = 
    vtkSmartPointer<vtkCharArray>::New();
  if(ZeroSheetType){
    vertexTypes->SetNumberOfTuples(vertexNumber);
    vertexTypes->SetName("SheetType");
  }
  else{
    sheet0->GetPointData()->RemoveArray("SheetType");
  }
  
  vtkSmartPointer<ttkIdTypeArray> vertexSheetId =
    vtkSmartPointer<ttkIdTypeArray>::New();
  if(ZeroSheetId){
    vertexSheetId->SetNumberOfTuples(vertexNumber);
    vertexSheetId->SetName("0-SheetId");
  }
  else{
    sheet0->GetPointData()->RemoveArray("0-SheetId");
  }
  
  vertexNumber = 0;
  double *p = NULL;
  for(ttkIdType i = 0; i < (ttkIdType) sheet0segmentation->size(); i++){
    ttkIdType sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1){
      
      const ReebSpace::Sheet0 *sheet = reebSpace_.get0sheet(sheet0Id);
      
      if(!sheet->pruned_){
        p = input->GetPoint(i);
        
        sheet0->GetPoints()->SetPoint(vertexNumber, p);
        
        if(ZeroSheetId){
          vertexSheetId->SetTuple1(vertexNumber, 
            (*sheet0segmentation)[i]);
        }
        if(ZeroSheetVertexId){
          vertexIds->SetTuple1(vertexNumber, i);
        }
        if(ZeroSheetValue){
          double u, v;
          
          uComponent_->GetTuple(i, &u);
          vComponent_->GetTuple(i, &v);
          
          vertexScalarsU->SetTuple1(vertexNumber, u);
          vertexScalarsV->SetTuple1(vertexNumber, v);
        }
        if(ZeroSheetType){
          const ReebSpace::Sheet0* sheet = 
            reebSpace_.get0sheet((*sheet0segmentation)[i]);
          vertexTypes->SetTuple1(vertexNumber, sheet->type_);
        }
      
        vertexNumber++;
      }
    }
  }
  if(ZeroSheetId)
    sheet0->GetPointData()->AddArray(vertexSheetId);
  if(ZeroSheetVertexId)
    sheet0->GetPointData()->AddArray(vertexIds);
  if(ZeroSheetValue){
    sheet0->GetPointData()->AddArray(vertexScalarsU);
    sheet0->GetPointData()->AddArray(vertexScalarsV);
  }
  if(ZeroSheetType)
    sheet0->GetPointData()->AddArray(vertexTypes);
    
  // 1-sheets
  // Optional additional fields:
  // PointData: u, v, vertexId, 
  // CellData: edgeId, type, sheetId
  const vector<ttkIdType> *sheet1segmentation = reebSpace_.get1sheetSegmentation();
  
  vtkSmartPointer<vtkPoints> sheet1Points = 
    vtkSmartPointer<vtkPoints>::New();
  sheet1->SetPoints(sheet1Points);
  
  vtkSmartPointer<vtkCellArray> sheet1Edges =
    vtkSmartPointer<vtkCellArray>::New();
  
  vtkSmartPointer<vtkDoubleArray> edgeScalarsU = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> edgeScalarsV = 
    vtkSmartPointer<vtkDoubleArray>::New();
    
  if(OneSheetValue){
    edgeScalarsU->SetName(uComponent_->GetName());
    edgeScalarsV->SetName(vComponent_->GetName());
  }
  else{
    sheet1->GetPointData()->RemoveArray(uComponent_->GetName());
    sheet1->GetPointData()->RemoveArray(vComponent_->GetName());
  }
 
  vtkSmartPointer<ttkIdTypeArray> edgeVertexIds = 
    vtkSmartPointer<ttkIdTypeArray>::New();
  if(OneSheetVertexId){
    edgeVertexIds->SetName("VertexIds");
  }
  else{
    sheet1->GetPointData()->RemoveArray("VertexIds");
  }
  
  vtkSmartPointer<vtkIntArray> edgeType = 
    vtkSmartPointer<vtkIntArray>::New();
  if(OneSheetType){
    edgeType->SetName("EdgeType");
  }
  else{
    sheet1->GetCellData()->RemoveArray("EdgeType");
  }
  
  vtkSmartPointer<ttkIdTypeArray> edgeIds = 
    vtkSmartPointer<ttkIdTypeArray>::New();
  if(OneSheetEdgeId){
    edgeIds->SetName("EdgeIds");
  }
  else{
    sheet1->GetCellData()->RemoveArray("EdgeIds");
  }
  
  vtkSmartPointer<ttkIdTypeArray> edgeSheetIds = 
    vtkSmartPointer<ttkIdTypeArray>::New();
  if(OneSheetId){
    edgeSheetIds->SetName("1-SheetId");
  }
  else{
    sheet1->GetCellData()->RemoveArray("1-SheetId");
  }
  
  vertexNumber = 0;
  double p0[3], p1[3];
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);
  const vector<ttkIdType> *edgeTypes = reebSpace_.getEdgeTypes();
  
  for(ttkIdType i = 0; i < (ttkIdType) sheet1segmentation->size(); i++){
    
    ttkIdType sheet1Id = (*sheet1segmentation)[i];
    
    if(sheet1Id != -1){
      
      const ReebSpace::Sheet1 *sheet = reebSpace_.get1sheet(sheet1Id);
      
      if((sheet)&&(!sheet->pruned_)){
      
        ttkIdType vertexId0 = -1, vertexId1 = -1;
        triangulation->getEdgeVertex(i, 0, vertexId0);
        triangulation->getEdgeVertex(i, 1, vertexId1);
        
        input->GetPoint(vertexId0, p0);
        input->GetPoint(vertexId1, p1);
        
        sheet1->GetPoints()->InsertNextPoint(p0);
        sheet1->GetPoints()->InsertNextPoint(p1);
        
        if(OneSheetValue){
          double u, v;
          uComponent_->GetTuple(vertexId0, &u);
          vComponent_->GetTuple(vertexId0, &v);
          
          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);
          
          uComponent_->GetTuple(vertexId1, &u);
          vComponent_->GetTuple(vertexId1, &v);
          
          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);
        }
        
        if(OneSheetVertexId){
          edgeVertexIds->InsertNextTuple1(vertexId0);
          edgeVertexIds->InsertNextTuple1(vertexId1);
        }
      
        idList->SetId(0, vertexNumber);
        idList->SetId(1, vertexNumber + 1);
        sheet1Edges->InsertNextCell(idList);
        vertexNumber += 2;
        
        if(OneSheetEdgeId){
          edgeIds->InsertNextTuple1(i);
        }
        if(OneSheetType){
          edgeType->InsertNextTuple1((*edgeTypes)[i]);
        }
        if(OneSheetId){
          edgeSheetIds->InsertNextTuple1((*sheet1segmentation)[i]);
        }
      }
    }
  }
  
  sheet1->SetCells(VTK_LINE, sheet1Edges);
  
  if(OneSheetValue){
    sheet1->GetPointData()->AddArray(edgeScalarsU);
    sheet1->GetPointData()->AddArray(edgeScalarsV);
  }
  if(OneSheetVertexId){
    sheet1->GetPointData()->AddArray(edgeVertexIds);
  }
  if(OneSheetId){
    sheet1->GetCellData()->AddArray(edgeSheetIds);
  }
  if(OneSheetEdgeId){
    sheet1->GetCellData()->AddArray(edgeIds);
  }
  if(OneSheetType){
    sheet1->GetCellData()->AddArray(edgeType);
  }
  
  // 2-sheets
  // optional fields:
  // pointdata: twoSheetValues, twoSheetParameterization
  if(TwoSheets){
    const vector<FiberSurface::Vertex> *vertexList = 
      reebSpace_.getFiberSurfaceVertices();
      
    vtkSmartPointer<vtkPoints> sheet2Points = vtkSmartPointer<vtkPoints>::New();
    sheet2Points->SetNumberOfPoints(vertexList->size());
    sheet2->SetPoints(sheet2Points);
      
    vtkSmartPointer<vtkDoubleArray> triangleScalarsU = 
      vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> triangleScalarsV =
      vtkSmartPointer<vtkDoubleArray>::New();
      
    if(TwoSheetValue){
      triangleScalarsU->SetName(uComponent_->GetName());
      triangleScalarsU->SetNumberOfTuples(vertexList->size());
      triangleScalarsV->SetName(vComponent_->GetName());
      triangleScalarsV->SetNumberOfTuples(vertexList->size());
    }
    else{
      sheet2->GetPointData()->RemoveArray(uComponent_->GetName());
      sheet2->GetPointData()->RemoveArray(vComponent_->GetName());
    }
    
    vtkSmartPointer<vtkDoubleArray> triangleParameterization = 
      vtkSmartPointer<vtkDoubleArray>::New();
    if(TwoSheetParameterization){
      triangleParameterization->SetName("EdgeParameterization");
      triangleParameterization->SetNumberOfTuples(vertexList->size());
    }
    else{
      sheet2->GetPointData()->RemoveArray("EdgeParameterization");
    }
    
    ttkIdType sheet2TriangleNumber = 0;
    for(ttkIdType i = 0; i < reebSpace_.getNumberOf2sheets(); i++){
      const ReebSpace::Sheet2 *sheet = reebSpace_.get2sheet(i);
      
      if(!sheet->pruned_){
        for(ttkIdType j = 0; j  < (ttkIdType) sheet->triangleList_.size(); j++){
          sheet2TriangleNumber += sheet->triangleList_[j].size();
        }
      }
    }
    
    vtkSmartPointer<vtkCellArray> sheet2Triangles = 
    vtkSmartPointer<vtkCellArray>::New();

    // celldata: twoSheetId, twoSheetEdgeId, twoSheetTetId
    vtkSmartPointer<ttkIdTypeArray> triangleSheetIds = 
      vtkSmartPointer<ttkIdTypeArray>::New();
    if(TwoSheetId){
      triangleSheetIds->SetName("2-SheetId");
      triangleSheetIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("2-SheetId");
    }
    
    vtkSmartPointer<ttkIdTypeArray> triangleEdgeIds = 
      vtkSmartPointer<ttkIdTypeArray>::New();
    if(TwoSheetEdgeId){
      triangleEdgeIds->SetName("EdgeIds");
      triangleEdgeIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("EdgeIds");
    }
    
    vtkSmartPointer<vtkIntArray> triangleEdgeType = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetEdgeType){
      triangleEdgeType->SetName("EdgeType");
      triangleEdgeType->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("EdgeType");
    }
    
    vtkSmartPointer<ttkIdTypeArray> triangleTetIds = 
      vtkSmartPointer<ttkIdTypeArray>::New();
    if(TwoSheetTetId){
      triangleTetIds->SetName("TetIds");
      triangleTetIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("TetIds");
    }
    
    vtkSmartPointer<ttkIdTypeArray> triangleCaseIds = 
      vtkSmartPointer<ttkIdTypeArray>::New();
    if(TwoSheetCaseId){
      triangleCaseIds->SetName("CaseIds");
      triangleCaseIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("CaseIds");
    }
    
    for(ttkIdType i = 0; i < (ttkIdType) vertexList->size(); i++){
      sheet2->GetPoints()->SetPoint(i, 
        (*vertexList)[i].p_[0], (*vertexList)[i].p_[1], (*vertexList)[i].p_[2]);
      
      if(TwoSheetValue){
        triangleScalarsU->SetTuple1(i, (*vertexList)[i].uv_.first);
        triangleScalarsV->SetTuple1(i, (*vertexList)[i].uv_.second);
      }
      if(TwoSheetParameterization){
        triangleParameterization->SetTuple1(i, (*vertexList)[i].t_);
      }
    }
    if(TwoSheetValue){
      sheet2->GetPointData()->AddArray(triangleScalarsU);
      sheet2->GetPointData()->AddArray(triangleScalarsV);
    }
    if(TwoSheetParameterization){
      sheet2->GetPointData()->AddArray(triangleParameterization);
    }

    ttkIdType triangleNumber = 0;
    idList->SetNumberOfIds(3);
    for(ttkIdType i = 0; i < reebSpace_.getNumberOf2sheets(); i++){
      const ReebSpace::Sheet2 *sheet = reebSpace_.get2sheet(i);
     
      if(!sheet->pruned_){
        for(ttkIdType j = 0; j < (ttkIdType) sheet->triangleList_.size(); j++){
          
          for(ttkIdType k = 0; k < (ttkIdType) sheet->triangleList_[j].size(); k++){
            
            for(int l = 0; l < 3; l++){
              idList->SetId(l, sheet->triangleList_[j][k].vertexIds_[l]);
            }
            
            sheet2Triangles->InsertNextCell(idList);
            
            if(TwoSheetId){
              triangleSheetIds->SetTuple1(triangleNumber, i);
            }
            
            if(TwoSheetEdgeId){
              const ReebSpace::Sheet1 *sheet1 = 
                reebSpace_.get1sheet(sheet->sheet1Id_);
              triangleEdgeIds->SetTuple1(triangleNumber,
                sheet1->edgeList_[j]);
            }
            
            if(TwoSheetEdgeType){
              ttkIdType polygonEdgeId = sheet->triangleList_[j][k].polygonEdgeId_;
              ttkIdType edgeId = reebSpace_.getJacobi2Edge(polygonEdgeId);
              triangleEdgeType->SetTuple1(triangleNumber, (*edgeTypes)[edgeId]);
            }
            
            if(TwoSheetTetId){
              triangleTetIds->SetTuple1(triangleNumber, 
                sheet->triangleList_[j][k].tetId_);
            }
            if(TwoSheetCaseId){
              triangleCaseIds->SetTuple1(triangleNumber,
                sheet->triangleList_[j][k].caseId_);
            }
            
            triangleNumber++;
          }
        }
      }
    }
    sheet2->SetCells(VTK_TRIANGLE, sheet2Triangles);
    
    if(TwoSheetId){
      sheet2->GetCellData()->AddArray(triangleSheetIds);
    }
    if(TwoSheetEdgeId){
      sheet2->GetCellData()->AddArray(triangleEdgeIds);
    }
    if(TwoSheetEdgeType){
      sheet2->GetCellData()->AddArray(triangleEdgeType);
    }
    if(TwoSheetTetId){
      sheet2->GetCellData()->AddArray(triangleTetIds);
    }
    if(TwoSheetCaseId){
      sheet2->GetCellData()->AddArray(triangleCaseIds);
    }
  }
  
  // now take care of the 3 sheets
//   vector<float> *triangulationPoints 
//     = reebSpace_.getSheetTriangulationPoints();
//   vector<long long int> *triangulationCells
//     = reebSpace_.getSheetTriangulationCells();
//     
//   vtkSmartPointer<vtkPoints> sheet3Points = vtkSmartPointer<vtkPoints>::New();
//   vtkSmartPointer<vtkFloatArray> pointData = 
//     vtkSmartPointer<vtkFloatArray>::New();
//   pointData->SetNumberOfComponents(3);
//   pointData->SetVoidArray(
//     triangulationPoints->data(), triangulationPoints->size(), 1);
//   sheet3Points->SetData(pointData);
//   sheet3->SetPoints(sheet3Points);
//   
//   vtkSmartPointer<vtkCellArray> sheet3Cells 
//     = vtkSmartPointer<vtkCellArray>::New();
//   vtkSmartPointer<ttkIdTypeArray> idArray 
//     = vtkSmartPointer<ttkIdTypeArray>::New();
//   idArray->SetVoidArray(
//     triangulationCells->data(), triangulationCells->size(), 1);
//   sheet3Cells->SetCells(triangulationCells->size()/5, idArray);
//   sheet3->SetCells(VTK_TETRA, sheet3Cells);

  // now take care of the 3 sheets
  sheet3->ShallowCopy(input);
  const vector<ttkIdType> *vertex3sheets = reebSpace_.get3sheetVertexSegmentation();
  
  vtkSmartPointer<ttkIdTypeArray> vertexNumberField
    = vtkSmartPointer<ttkIdTypeArray>::New();
  vtkSmartPointer<ttkIdTypeArray> tetNumberField
    = vtkSmartPointer<ttkIdTypeArray>::New();
  
  if(ThreeSheetTetNumber){
    tetNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    tetNumberField->SetName("3-SheetTetNumber");
    for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3
        = reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_))
        tetNumberField->SetTuple1(i, sheet3->tetList_.size());
      else
        tetNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(tetNumberField);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }
  
  if(ThreeSheetVertexNumber){
    vertexNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    vertexNumberField->SetName("3-SheetVertexNumber");
    for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3
        = reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_))
        vertexNumberField->SetTuple1(i, sheet3->vertexList_.size());
      else
        vertexNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(vertexNumberField);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }
  
  vtkSmartPointer<vtkDoubleArray> domainVolume = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    domainVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    domainVolume->SetName("3-SheetDomainVolume");
    
    for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        domainVolume->SetTuple1(i, sheet3->domainVolume_);
      }
      else{
        domainVolume->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(domainVolume);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetDomainVolume");
  }
  
  vtkSmartPointer<vtkDoubleArray> rangeArea = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    rangeArea->SetNumberOfTuples(input->GetNumberOfPoints());
    rangeArea->SetName("3-SheetRangeArea");
    
    for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        rangeArea->SetTuple1(i, sheet3->rangeArea_);
      }
      else{
        rangeArea->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(rangeArea);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetRangeArea");
  }
  
  vtkSmartPointer<vtkDoubleArray> hyperVolume = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    hyperVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    hyperVolume->SetName("3-SheetHyperVolume");
    
    for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        hyperVolume->SetTuple1(i, sheet3->hyperVolume_);
      }
      else{
        hyperVolume->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(hyperVolume);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetHyperVolume");
  }
  
  vtkSmartPointer<ttkIdTypeArray> vertexSegmentation 
    = vtkSmartPointer<ttkIdTypeArray>::New();
  vertexSegmentation->SetName("3-SheetId");
  vertexSegmentation->SetNumberOfTuples(input->GetNumberOfPoints());
  for(ttkIdType i = 0; i < input->GetNumberOfPoints(); i++){
    const ReebSpace::Sheet3 *sheet = reebSpace_.get3sheet((*vertex3sheets)[i]);
    if(sheet){
      vertexSegmentation->SetTuple1(i, sheet->simplificationId_);
    }
    else{
      vertexSegmentation->SetTuple1(i, (*vertex3sheets)[i]);
    }
  }
  sheet3->GetPointData()->AddArray(vertexSegmentation);
  
  const vector<ttkIdType> *tet3sheets = reebSpace_.get3sheetTetSegmentation();
  vtkSmartPointer<ttkIdTypeArray> tetSegmentation 
    = vtkSmartPointer<ttkIdTypeArray>::New();
  tetSegmentation->SetName("3-SheetId");
  tetSegmentation->SetNumberOfTuples(input->GetNumberOfCells());
  for(ttkIdType i = 0; i < input->GetNumberOfCells(); i++){
    const ReebSpace::Sheet3 *sheet = reebSpace_.get3sheet((*tet3sheets)[i]);
    if(sheet){
      tetSegmentation->SetTuple1(i, sheet->simplificationId_);
    }
    else{
      tetSegmentation->SetTuple1(i, (*tet3sheets)[i]);
    }
  }
  sheet3->GetCellData()->AddArray(tetSegmentation);

  {
    stringstream msg;
    msg << "[ttkReebSpace] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
