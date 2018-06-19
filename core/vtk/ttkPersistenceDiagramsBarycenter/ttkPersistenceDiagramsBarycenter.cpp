#include                  <ttkPersistenceDiagramsBarycenter.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceDiagramsBarycenter)

ttkPersistenceDiagramsBarycenter::ttkPersistenceDiagramsBarycenter(){

  UseAllCores = false;

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



int ttkPersistenceDiagramsBarycenter::doIt(vtkDataSet** input, int numInputs){

	std::cout<<"Hello world ! " << std::endl;

	// Use a pointer-base copy for the input data
	/*outputBoundFields->ShallowCopy(input[0]);
	outputProbability->ShallowCopy(input[0]);
	outputMean->ShallowCopy(input[0]);*/

	// Get arrays from input datas
	//vtkDataArray* inputScalarField[numInputs] = { NULL };
	vector<vtkDataArray*> inputScalarField(numInputs);
	for(int i=0 ; i<numInputs ; i++){
		if(ScalarField.length()){
			inputScalarField[i] = input[i]->GetPointData()->GetArray(ScalarField.data());
		}
		else{
			inputScalarField[i] = input[i]->GetPointData()->GetArray(0);
		}
		// Check if inputs have the same data type and the same number of points
		if(inputScalarField[i]->GetDataType() != inputScalarField[0]->GetDataType()){
			stringstream msg;
			msg << "[ttkPersistenceDiagramsBarycenter] Inputs of different data types." << endl;
			dMsg(cerr, msg.str(), fatalMsg);
			return -3;
		}
		if(inputScalarField[i]->GetNumberOfTuples() != inputScalarField[0]->GetNumberOfTuples()){
			stringstream msg;
			msg << "[ttkPersistenceDiagramsBarycenter] Inputs with different number of points." << endl;
			dMsg(cerr, msg.str(), fatalMsg);
			return -2;
		}
		// Check if all the inputs are here
		if(!inputScalarField[i]){
			return -1;
		}
	}

	std::cout<<"Hello world2 ! " << std::endl;
  // Calling the executing package
	switch(inputScalarField[0]->GetDataType()){

		vtkTemplateMacro(
		{
			PersistenceDiagramsBarycenter persistenceDiagramsBarycenter;
			persistenceDiagramsBarycenter.setWrapper(this);

			persistenceDiagramsBarycenter.setNumberOfInputs(numInputs);
			for (int i = 0; i<numInputs; i++) {
				persistenceDiagramsBarycenter.setInputDataPointer(i, 
				inputScalarField[i]->GetVoidPointer(0));
			}

			persistenceDiagramsBarycenter.execute<VTK_TT>();
		});
	}

  return 0;
}



// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramsBarycenter::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){
	std::cout<<"Hello world ? " << std::endl;
  Memory m;

  /*
  // Output pointers informations
  vtkInformation* outInfo;
  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *boundFields = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Histograms
  outInfo = outputVector->GetInformationObject(1);
  vtkDataSet *probability = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Mean field
  outInfo = outputVector->GetInformationObject(2);
  vtkDataSet *mean = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));*/
  
  
  std::cout<<"Hello world2 ? " << std::endl;

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramsBarycenter] Number of inputs: " << numInputs << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Get input datas
  vtkDataSet* *input = new vtkDataSet*[numInputs];
  for(int i=0 ; i<numInputs ; i++)
  {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  doIt(input, numInputs);

  delete[] input;

  {
    stringstream msg;
    msg << "[ttkPersistenceDiagramsBarycenter] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
