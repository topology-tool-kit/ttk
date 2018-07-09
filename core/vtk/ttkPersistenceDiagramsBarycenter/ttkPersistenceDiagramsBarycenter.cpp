#include                  <ttkPersistenceDiagramsBarycenter.h>

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

vtkStandardNewMacro(ttkPersistenceDiagramsBarycenter)

ttkPersistenceDiagramsBarycenter::ttkPersistenceDiagramsBarycenter(){
  UseAllCores = false;
  WassersteinMetric = "2";
  UseOutputMatching = true;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
  
}

ttkPersistenceDiagramsBarycenter::~ttkPersistenceDiagramsBarycenter(){}


// transmit abort signals -- to copy paste in other wrappers
bool ttkPersistenceDiagramsBarycenter::needsToAbort(){
	return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkPersistenceDiagramsBarycenter::updateProgress(const float &progress){
	{
		stringstream msg;
		msg << "[ttkPersistenceDiagramsBarycenter] " << progress*100
			<< "% processed...." << endl;
		dMsg(cout, msg.str(), advancedInfoMsg);
	}

	UpdateProgress(progress);
	return 0;
}



int ttkPersistenceDiagramsBarycenter::doIt(vtkDataSet** input, vtkUnstructuredGrid *outputBarycenter, vtkUnstructuredGrid *outputMatchings, vtkUnstructuredGrid *outputDiagrams, int numInputs){
	// Get arrays from input datas
	//vtkDataArray* inputDiagram[numInputs] = { NULL };
	vector<vtkUnstructuredGrid*> inputDiagram(numInputs);
	for(int i=0 ; i<numInputs ; i++){
		inputDiagram[i] = vtkUnstructuredGrid::SafeDownCast(input[i]);
	}
  // Calling the executing package
	
	int dataType = inputDiagram[0]->GetCellData()->GetArray("Persistence")->GetDataType();
	
	// TODO If Windows, we need to get rid of one pair of parenthesis
	switch(dataType){

		vtkTemplateMacro((
		{
			PersistenceDiagramsBarycenter<VTK_TT> persistenceDiagramsBarycenter;
			persistenceDiagramsBarycenter.setWrapper(this);
			
			string wassersteinMetric = WassersteinMetric;
			persistenceDiagramsBarycenter.setWasserstein(wassersteinMetric);

			persistenceDiagramsBarycenter.setNumberOfInputs(numInputs);
			persistenceDiagramsBarycenter.setTimeLimit(TimeLimit);
			persistenceDiagramsBarycenter.setUseProgressive(UseProgressive);
			persistenceDiagramsBarycenter.setThreadNumber(ThreadNumber);
			persistenceDiagramsBarycenter.setAlpha(Alpha);
			
      std::vector<std::vector<macroDiagramTuple> > 
        intermediateDiagrams(numInputs);
      for(int i = 0; i < numInputs; i++){
        double Spacing = 0;
        getPersistenceDiagram<VTK_TT>(
          &(intermediateDiagrams[i]), inputDiagram[i], Spacing, 0);
        persistenceDiagramsBarycenter.setDiagram(i, 
          (void*) &(intermediateDiagrams[i]));
      }
			
			std::vector<macroDiagramTuple> barycenter;
			std::vector<std::vector<macroMatchingTuple>> matchings = persistenceDiagramsBarycenter.execute(&barycenter);
			
			outputBarycenter->ShallowCopy(createPersistenceDiagram<VTK_TT>(&barycenter));
			outputMatchings->ShallowCopy(createMatchings(&matchings, &barycenter, 
        intermediateDiagrams));
			outputDiagrams->ShallowCopy(createOutputDiagrams(intermediateDiagrams));
		}
		));
	}

  return 0;
}


int ttkPersistenceDiagramsBarycenter::FillInputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillInputPortInformation(port, info)){
    return 0;
  }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkPersistenceDiagramsBarycenter::FillOutputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillOutputPortInformation(port, info)){
    return 0;
  }
  if(port==0 || port==1 || port==3)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}


// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramsBarycenter::RequestData(vtkInformation *request,
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
  for(int i=0 ; i<numInputs ; i++)
  {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
	if(!input[i]){
		std::cout<<"No data in input["<<i<<"]"<<std::endl;
	}
  }
  
  // TODO Set output
  
  /*vtkInformation* outInfo;
  outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *outputCT1 = vtkUnstructuredGrid::SafeDownCast(outInfo);*/
  
  vtkInformation* outInfo1;
  outInfo1 = outputVector->GetInformationObject(0);
  vtkDataSet *output1 = vtkDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *outputCT1 = vtkUnstructuredGrid::SafeDownCast(output1);
  
  vtkInformation* outInfo2;
  outInfo2 = outputVector->GetInformationObject(1);
  vtkDataSet *output2 = vtkDataSet::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *matchings = vtkUnstructuredGrid::SafeDownCast(output2);
  
  vtkInformation* outInfo3;
  outInfo3 = outputVector->GetInformationObject(2);
  vtkDataSet *output3 = vtkDataSet::SafeDownCast(outInfo3->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *outputDiagrams = vtkUnstructuredGrid::SafeDownCast(output3);
  
  
  
  doIt(input, outputCT1, matchings, outputDiagrams, numInputs);
  
  
  
  
  /*vtkUnstructuredGrid *outputCT1 = vtkUnstructuredGrid::SafeDownCast(output1);
  outputCT1->ShallowCopy(outputBarycenter);*/

  delete[] input;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramsBarycenter] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
