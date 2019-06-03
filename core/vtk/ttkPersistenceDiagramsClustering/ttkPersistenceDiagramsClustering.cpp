#include                  <ttkPersistenceDiagramsClustering.h>

#ifndef macroDiagramTuple
#define macroDiagramTuple std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId, \
  ttk::CriticalType, VTK_TT, ttk::SimplexId, \
  VTK_TT, float, float, float, VTK_TT, float, float, float>
#endif

#ifndef macroMatchingTuple
#define macroMatchingTuple std::tuple<ttk::SimplexId, ttk::SimplexId, VTK_TT>
#endif

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramsClustering)

ttkPersistenceDiagramsClustering::ttkPersistenceDiagramsClustering(){
// std::cout<<"constructor"<<std::endl;
  UseAllCores = false;
  WassersteinMetric = "2";
  UseOutputMatching = true;
  TimeLimit=9999;
  NumberOfClusters=1;
  Deterministic = 1;
  ThreadNumber = 1;
  PairTypeClustering  = -1;
  numberOfInputsFromCommandLine=1;
  UseProgressive = 1;
  UseAccelerated = 0;
  UseKmeansppInit = 0;
  Alpha = 1;
  Lambda = 1;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);

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



int ttkPersistenceDiagramsClustering::doIt(vtkDataSet** input, vtkUnstructuredGrid *outputClusters, vtkUnstructuredGrid *outputCentroids, int numInputs){
	// Get arrays from input datas
    // std::cout<<"STARTING doIt"<<std::endl;
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

            // std::cout<<"setting the PersistenceDiagramsClustering parameters"<<std::endl;

			persistenceDiagramsClustering.setWasserstein(wassersteinMetric);
            persistenceDiagramsClustering.setDeterministic(Deterministic);
            persistenceDiagramsClustering.setPairTypeClustering(PairTypeClustering);
			persistenceDiagramsClustering.setNumberOfInputs(numInputs);
			persistenceDiagramsClustering.setDebugLevel(debugLevel_);
			persistenceDiagramsClustering.setTimeLimit(TimeLimit);
			persistenceDiagramsClustering.setUseProgressive(UseProgressive);
			persistenceDiagramsClustering.setThreadNumber(ThreadNumber);
			persistenceDiagramsClustering.setAlpha(Alpha);
            persistenceDiagramsClustering.setLambda(Lambda);
			persistenceDiagramsClustering.setNumberOfClusters(NumberOfClusters);
			persistenceDiagramsClustering.setUseAccelerated(UseAccelerated);
			persistenceDiagramsClustering.setUseKmeansppInit(UseKmeansppInit);

			std::vector<std::vector<macroDiagramTuple> >
				intermediateDiagrams(numInputs);
			double max_dimension_total=0;
			for(int i = 0; i < numInputs; i++){
				double Spacing = 0;
				double max_dimension = getPersistenceDiagram<VTK_TT>(
				&(intermediateDiagrams[i]), inputDiagram[i], Spacing, 0);
				if(max_dimension_total<max_dimension){
                    max_dimension_total=max_dimension;
				}
			}
			persistenceDiagramsClustering.setDiagrams((void *) &intermediateDiagrams);


            std::vector<std::vector<macroDiagramTuple>> final_centroids;
			
            
            // std::cout<<"launching execute PerssistenceDiagramsCLustering"<<std::endl;
			std::vector<int> inv_clustering = 
                persistenceDiagramsClustering.execute(&final_centroids);
			outputClusters->ShallowCopy(createOutputClusteredDiagrams(intermediateDiagrams, inv_clustering, max_dimension_total));
			outputCentroids->ShallowCopy(createOutputCentroids<VTK_TT>(&final_centroids, inv_clustering, max_dimension_total));
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
  // std::cout<<"requestData"<<std::endl;
  int numInputs = numberOfInputsFromCommandLine;
  if(numInputs ==1){
    numInputs = inputVector[0]->GetNumberOfInformationObjects();
    // std::cout<<"1 inputs with "<< numInputs <<" objects"<<std::endl;
  }
  {
    stringstream msg;
    dMsg(cout, msg.str(), infoMsg);
  }
  // Get input datas
  vtkDataSet* *input = new vtkDataSet*[numInputs];
  for(int i=0 ; i<numInputs ; ++i)
  {
    if(numberOfInputsFromCommandLine>1){
        input[i] = vtkDataSet::GetData(inputVector[i], 0);
    }
    else{
        input[i] = vtkDataSet::GetData(inputVector[0], i);
    }
	if(!input[i]){
		std::cout<<"No data in input["<<i<<"]"<<std::endl;
	}
  }
  // TODO Set output
  // std::cout<<"settings outputs"<<std::endl;
  vtkInformation* outInfo1;
  outInfo1 = outputVector->GetInformationObject(0);
  vtkDataSet* output1 = vtkDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid* output_clusters = vtkUnstructuredGrid::SafeDownCast(output1);

  vtkInformation* outInfo2;
  outInfo2 = outputVector->GetInformationObject(1);
  vtkDataSet* output2 = vtkDataSet::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid* output_centroids = vtkUnstructuredGrid::SafeDownCast(output2);

  doIt(input, output_clusters, output_centroids, numInputs);
  delete[] input;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramsClustering] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
