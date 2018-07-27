#include <ttkPersistenceDiagram.h>

using namespace std;
using namespace ttk;
using namespace dcg;

vtkStandardNewMacro(ttkPersistenceDiagram)

  ttkPersistenceDiagram::ttkPersistenceDiagram():
    UseAllCores{},
    inputScalars_{},
    CTPersistenceDiagram_{vtkUnstructuredGrid::New()},
    offsets_{},
    inputOffsets_{},
    varyingMesh_{}
{
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  ComputeSaddleConnectors = false;
  InputOffsetScalarFieldName = ttk::OffsetScalarFieldName;
  ForceInputOffsetScalarField = false;
  ComputeSaddleConnectors = false;
  UseAllCores = true;
  ThreadNumber = 1;
  debugLevel_ = 3;
  ShowInsideDomain = false;
  computeDiagram_= true;

  triangulation_ = nullptr;
  CTDiagram_ = nullptr;
}

ttkPersistenceDiagram::~ttkPersistenceDiagram(){
  if(CTPersistenceDiagram_)
    CTPersistenceDiagram_->Delete();
  if(offsets_)
    offsets_->Delete();

  deleteDiagram();
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port, 
  vtkInformation *info){
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
  }

  return 1;
}

int ttkPersistenceDiagram::getScalars(vtkDataSet* input){
  vtkPointData* pointData=input->GetPointData();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!pointData){
    cerr << "[ttkPersistenceDiagram] Error : input has no point data." << endl;
    return -1;
  }
#endif

  if(ScalarField.length()){
    inputScalars_=pointData->GetArray(ScalarField.data());
  }
  else{
    inputScalars_=pointData->GetArray(ScalarFieldId);
    if(inputScalars_)
      ScalarField = inputScalars_->GetName();
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars_){
    cerr 
      << "[ttkPersistenceDiagram] Error : input scalar field pointer is null." 
      << endl;
    return -3;
  }
#endif

  stringstream msg;
  msg << "[ttkPersistenceDiagram] Starting computation on field `"
    << inputScalars_->GetName() << "'..." << endl;
  dMsg(cout, msg.str(), infoMsg);

  return 0;
}

int ttkPersistenceDiagram::getTriangulation(vtkDataSet* input){
  varyingMesh_=false;

  triangulation_ = ttkTriangulation::getTriangulation(input);
  if(!triangulation_)
    return 0;
  
  triangulation_->setWrapper(this);
  persistenceDiagram_.setupTriangulation(triangulation_);
  
  if(triangulation_->isEmpty() 
    or ttkTriangulation::hasChangedConnectivity(triangulation_, input, this)){
    Modified();
    varyingMesh_=true;
    computeDiagram_=true;
  }

  return 0;
}

int ttkPersistenceDiagram::getOffsets(vtkDataSet* input){
  if(OffsetFieldId != -1){
    inputOffsets_ = input->GetPointData()->GetArray(OffsetFieldId);
    if(inputOffsets_){
      InputOffsetScalarFieldName = inputOffsets_->GetName();
      ForceInputOffsetScalarField = true;
    }
  }

  if(ForceInputOffsetScalarField and InputOffsetScalarFieldName.length()){
    inputOffsets_=
      input->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
  }
  else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)){
    inputOffsets_=
      input->GetPointData()->GetArray(ttk::OffsetScalarFieldName);
  }
  else{
    if(varyingMesh_ and offsets_){
      offsets_->Delete();
      offsets_=nullptr;
    }

    if(!offsets_){
      const ttkIdType numberOfVertices=input->GetNumberOfPoints();

      offsets_=ttkIdTypeArray::New();
      offsets_->SetNumberOfComponents(1);
      offsets_->SetNumberOfTuples(numberOfVertices);
      offsets_->SetName(ttk::OffsetScalarFieldName);
      for(ttkIdType i=0; i<numberOfVertices; ++i)
        offsets_->SetTuple1(i,i);
    }

    inputOffsets_=offsets_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOffsets_){
    cerr << "[ttkPersistenceDiagram] Error : wrong input offset scalar field." 
      << endl;
    return -1;
  }
#endif

  return 0;
}
#ifdef _MSC_VER
#define COMMA ,
#endif 
int ttkPersistenceDiagram::deleteDiagram(){
  if(CTDiagram_ and inputScalars_){
    switch(inputScalars_->GetDataType()){
#ifndef _MSC_VER
		vtkTemplateMacro(({
			using tuple_t = tuple<ttkIdType,ftm::NodeType,ttkIdType,ftm::NodeType,VTK_TT,ttkIdType>;
		vector<tuple_t>* CTDiagram = (vector<tuple_t>*) CTDiagram_;
		delete CTDiagram;
		}));
#else
		vtkTemplateMacro({
			using tuple_t = tuple<ttkIdType COMMA ftm::NodeType COMMA ttkIdType COMMA ftm::NodeType COMMA VTK_TT COMMA ttkIdType>;
		vector<tuple_t>* CTDiagram = (vector<tuple_t>*) CTDiagram_;
		delete CTDiagram;
		});
#endif
    }
  }
  return 0;
}

int ttkPersistenceDiagram::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){
  
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *outputCTPersistenceDiagram = 
    vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  
  Memory m;
  
  int ret{};

  ret=getScalars(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong scalars." << endl;
    return -1;
  }
#endif

  ret=getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong triangulation." << endl;
    return -2;
  }
#endif

  ret=getOffsets(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(ret){
    cerr << "[ttkPersistenceDiagram] Error : wrong offsets." << endl;
    return -3;
  }
  if(inputOffsets_->GetDataType()!=VTK_INT and inputOffsets_->GetDataType()!=VTK_ID_TYPE){
    cerr << "[ttkPersistenceDiagram] Error : input offset field type not supported." << endl;
    return -1;
  }
#endif

  vector<tuple<Cell,Cell>> dmt_pairs;
  persistenceDiagram_.setDMTPairs(&dmt_pairs);

  persistenceDiagram_.setWrapper(this);
  persistenceDiagram_.setInputScalars(inputScalars_->GetVoidPointer(0));
  persistenceDiagram_.setInputOffsets(inputOffsets_->GetVoidPointer(0));
  persistenceDiagram_.setComputeSaddleConnectors(ComputeSaddleConnectors);
  switch(inputScalars_->GetDataType()){
#ifndef _MSC_VER
	  vtkTemplateMacro(({
		  using tuple_t = tuple<ttkIdType,
		  ftm::NodeType,
		  ttkIdType,
		  ftm::NodeType,
		  VTK_TT,
		  ttkIdType>;

	  if (CTDiagram_ and computeDiagram_) {
		  vector<tuple_t>* tmpDiagram = (vector<tuple_t>*)CTDiagram_;
		  delete tmpDiagram;
		  CTDiagram_ = new vector<tuple_t>();
	  }
	  else if (!CTDiagram_) {
		  CTDiagram_ = new vector<tuple_t>();
		  computeDiagram_ = true;
	  }

	  vector<tuple_t>* CTDiagram = (vector<tuple_t>*)CTDiagram_;

	  if (computeDiagram_) {
		  persistenceDiagram_.setOutputCTDiagram(CTDiagram);
          if(inputOffsets_->GetDataType()==VTK_INT)
            ret = persistenceDiagram_.execute<VTK_TT,int>();
          if(inputOffsets_->GetDataType()==VTK_ID_TYPE)
            ret = persistenceDiagram_.execute<VTK_TT,vtkIdType>();
#ifndef TTK_ENABLE_KAMIKAZE
		  if (ret) {
			  cerr << "[ttkPersistenceDiagram] PersistenceDiagram.execute() "
				  << "error code : " << ret << endl;
			  return -4;
		  }
#endif
	  }

	  if (ShowInsideDomain)
		  ret = getPersistenceDiagramInsideDomain<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
	  else
		  ret = getPersistenceDiagram<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
#ifndef TTK_ENABLE_KAMIKAZE
	  if (ret) {
		  cerr << "[ttkPersistenceDiagram] Error : "
			  << "build of contour tree persistence diagram has failed." << endl;
		  return -5;
	  }
#endif
	  }));
#else
#ifndef TTK_ENABLE_KAMIKAZE
	  vtkTemplateMacro({
		  using tuple_t = tuple<ttkIdType COMMA
		  ftm::NodeType COMMA
		  ttkIdType COMMA
		  ftm::NodeType COMMA
		  VTK_TT COMMA
		  ttkIdType>;

	  if (CTDiagram_ and computeDiagram_) {
		  vector<tuple_t>* tmpDiagram = (vector<tuple_t>*)CTDiagram_;
		  delete tmpDiagram;
		  CTDiagram_ = new vector<tuple_t>();
	  }
	  else if (!CTDiagram_) {
		  CTDiagram_ = new vector<tuple_t>();
		  computeDiagram_ = true;
	  }

	  vector<tuple_t>* CTDiagram = (vector<tuple_t>*)CTDiagram_;

	  if (computeDiagram_) {
		  persistenceDiagram_.setOutputCTDiagram(CTDiagram);
          if(inputOffsets_->GetDataType()==VTK_INT)
            ret = persistenceDiagram_.execute<VTK_TT COMMA int>();
          if(inputOffsets_->GetDataType()==VTK_ID_TYPE)
            ret = persistenceDiagram_.execute<VTK_TT COMMA vtkIdType>();
		  if (ret) {
			  cerr << "[ttkPersistenceDiagram] PersistenceDiagram.execute() "
				  << "error code : " << ret << endl;
			  return -4;
		  }
	  }

	  if (ShowInsideDomain)
		  ret = getPersistenceDiagramInsideDomain<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
	  else
		  ret = getPersistenceDiagram<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
	  if (ret) {
		  cerr << "[ttkPersistenceDiagram] Error : "
			  << "build of contour tree persistence diagram has failed." << endl;
		  return -5;
	  }
	  });
#else
	  vtkTemplateMacro({
		  using tuple_t = tuple<ttkIdType COMMA
		  ftm::NodeType COMMA
		  ttkIdType COMMA
		  ftm::NodeType COMMA
		  VTK_TT COMMA
		  ttkIdType>;

	  if (CTDiagram_ and computeDiagram_) {
		  vector<tuple_t>* tmpDiagram = (vector<tuple_t>*)CTDiagram_;
		  delete tmpDiagram;
		  CTDiagram_ = new vector<tuple_t>();
	  }
	  else if (!CTDiagram_) {
		  CTDiagram_ = new vector<tuple_t>();
		  computeDiagram_ = true;
	  }

	  vector<tuple_t>* CTDiagram = (vector<tuple_t>*)CTDiagram_;

	  if (computeDiagram_) {
		  persistenceDiagram_.setOutputCTDiagram(CTDiagram);
          if(inputOffsets_->GetDataType()==VTK_INT)
            ret = persistenceDiagram_.execute<VTK_TT COMMA int>();
          if(inputOffsets_->GetDataType()==VTK_ID_TYPE)
            ret = persistenceDiagram_.execute<VTK_TT COMMA vtkIdType>();
	  }

	  if (ShowInsideDomain)
		  ret = getPersistenceDiagramInsideDomain<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
	  else
		  ret = getPersistenceDiagram<VTK_TT>(ftm::TreeType::Contour, *CTDiagram);
	  });
#endif
#endif
  }

  outputCTPersistenceDiagram->ShallowCopy(CTPersistenceDiagram_);
  computeDiagram_=false;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagram] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
