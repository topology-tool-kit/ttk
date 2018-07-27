#include                  <ttkPersistenceDiagramsClustering.h>

#ifndef macroDiagramTuple
#define macroDiagramTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::NodeType, ttk::ftm::idVertex, \
  ttk::ftm::NodeType, VTK_TT, ttk::ftm::idVertex, \
  VTK_TT, float, float, float, VTK_TT, float, float, float>
#endif

#ifndef macroMatchingTuple
#define macroMatchingTuple std::tuple<ttk::ftm::idVertex, ttk::ftm::idVertex, VTK_TT>
#endif

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramsClustering)

ttkPersistenceDiagramsClustering::ttkPersistenceDiagramsClustering(){
  UseAllCores = false;
  WassersteinMetric = "2";
  UseOutputMatching = true;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  
}

ttkPersistenceDiagramsClustering::~ttkPersistenceDiagramsClustering(){}


// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceDiagramsClustering::needsToAbort(){
	return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceDiagramsClustering::updateProgress(const float &progress){
	{
		stringstream msg;
		msg << "[ttkPersistenceDiagramsClustering] " << progress*100
			<< "% processed...." << endl;
		dMsg(cout, msg.str(), advancedInfoMsg);
	}

	UpdateProgress(progress);
	return 0;
}



int ttkPersistenceDiagramsClustering::doIt(vtkDataSet** input, int numInputs){
	// Get arrays from input datas
	//vtkDataArray* inputDiagram[numInputs] = { NULL };
	vector<vtkUnstructuredGrid*> inputDiagram(numInputs);
	for(int i=0 ; i<numInputs ; ++i){
		inputDiagram[i] = vtkUnstructuredGrid::SafeDownCast(input[i]);
	}
  // Calling the executing package
	
	int dataType = inputDiagram[0]->GetCellData()->GetArray("Persistence")->GetDataType();
	
	// TODO If Windows, we need to get rid of one pair of parenthesis
	switch(dataType){

		vtkTemplateMacro((
		{
			PersistenceDiagramsClustering<VTK_TT> persistenceDiagramsClustering;
			persistenceDiagramsClustering.setWrapper(this);
			
			string wassersteinMetric = WassersteinMetric;
			persistenceDiagramsClustering.setWasserstein(wassersteinMetric);

			persistenceDiagramsClustering.setNumberOfInputs(numInputs);
			persistenceDiagramsClustering.setTimeLimit(TimeLimit);
			persistenceDiagramsClustering.setUseProgressive(UseProgressive);
			persistenceDiagramsClustering.setThreadNumber(ThreadNumber);
			persistenceDiagramsClustering.setAlpha(Alpha);
			persistenceDiagramsClustering.setNumberOfClusters(NumberOfClusters);
			persistenceDiagramsClustering.setUseAccelerated(UseAccelerated);
			persistenceDiagramsClustering.setUseKmeansppInit(UseKmeansppInit);
			
			std::vector<std::vector<macroDiagramTuple> > 
				intermediateDiagrams(numInputs);
			for(int i = 0; i < numInputs; i++){
				double Spacing = 0;
				getPersistenceDiagram<VTK_TT>(
				&(intermediateDiagrams[i]), inputDiagram[i], Spacing, 0);
			}
			persistenceDiagramsClustering.setDiagrams((void *) &intermediateDiagrams);
		
				
			std::vector<macroDiagramTuple> barycenter;
			std::vector<std::vector<macroMatchingTuple>> matchings = persistenceDiagramsClustering.execute(&barycenter);
		}
		));
	}

  return 0;
}


int ttkPersistenceDiagramsClustering::FillInputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillInputPortInformation(port, info)){
    return 0;
  }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkPersistenceDiagramsClustering::FillOutputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillOutputPortInformation(port, info)){
    return 0;
  }
  if(port==0 || port==1 || port==3)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}


// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramsClustering::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){
  Memory m;
  
  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    stringstream msg;
    dMsg(cout, msg.str(), infoMsg);
  }
  // Get input datas
  vtkDataSet* *input = new vtkDataSet*[numInputs];
  for(int i=0 ; i<numInputs ; ++i)
  {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
	if(!input[i]){
		std::cout<<"No data in input["<<i<<"]"<<std::endl;
	}
  }
  // TODO Set output
  doIt(input, numInputs);

  delete[] input;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramsClustering] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
