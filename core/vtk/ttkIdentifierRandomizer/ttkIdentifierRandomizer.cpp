#include <ttkIdentifierRandomizer.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

vtkStandardNewMacro(ttkIdentifierRandomizer);

ttkIdentifierRandomizer::ttkIdentifierRandomizer(){
    this->setDebugMsgPrefix("IdentifierRandomizer");

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}
ttkIdentifierRandomizer::~ttkIdentifierRandomizer(){
}

int ttkIdentifierRandomizer::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifierRandomizer::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkIdentifierRandomizer::RequestData(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector
) {
  ttk::Timer t;

  this->printMsg("Randomizing Identifiers",0,0,
    this->threadNumber_, ttk::debug::LineMode::REPLACE);

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  auto inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray){
    this->printErr("Unable to retrieve input array 0.");
    return 0;
  }

  auto inputArrayAssociation = this->GetInputArrayAssociation(0, inputVector);

  auto outputArray = vtkSmartPointer<vtkDataArray>::Take( inputArray->NewInstance() );
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());
  outputArray->SetName(inputArray->GetName());

  // on the output, replace the field array by a pointer to its processed
  // version
  if(inputArrayAssociation==0) {
    output->GetPointData()->AddArray(outputArray);
  } else {
    output->GetCellData()->AddArray(outputArray);
  }

  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> identifierMap;

  for(ttk::SimplexId i = 0; i < inputArray->GetNumberOfTuples(); i++) {
    double inputIdentifier = -1;
    inputArray->GetTuple(i, &inputIdentifier);

    bool isIn = false;
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)identifierMap.size(); j++) {
      if(identifierMap[j].first == inputIdentifier) {
        isIn = true;
        break;
      }
    }

    if(!isIn) {
      identifierMap.push_back(
        std::pair<ttk::SimplexId, ttk::SimplexId>(inputIdentifier, -INT_MAX));
    }
  }

  // now let's shuffle things around
  ttk::SimplexId freeIdentifiers = identifierMap.size();
  ttk::SimplexId randomIdentifier = -1;

  for(ttk::SimplexId i = 0; i < (ttk::SimplexId)identifierMap.size(); i++) {

    randomIdentifier = drand48() * (freeIdentifiers);

    ttk::SimplexId freeCounter = -1;
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)identifierMap.size(); j++) {

      bool isFound = false;
      for(ttk::SimplexId k = 0; k < (ttk::SimplexId)identifierMap.size(); k++) {

        if(identifierMap[k].second == identifierMap[j].first) {
          isFound = true;
          break;
        }
      }
      if(!isFound) {
        freeCounter++;
      }

      if(freeCounter >= randomIdentifier) {
        randomIdentifier = identifierMap[j].first;
        break;
      }
    }

    identifierMap[i].second = randomIdentifier;
    freeIdentifiers--;
  }

  // now populate the output scalar field.
  for(ttk::SimplexId i = 0; i < inputArray->GetNumberOfTuples(); i++) {
    double inputIdentifier = -1;
    double outputIdentifier = -1;

    inputArray->GetTuple(i, &inputIdentifier);

    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)identifierMap.size(); j++) {
      if(inputIdentifier == identifierMap[j].first) {
        outputIdentifier = identifierMap[j].second;
        break;
      }
    }

    outputArray->SetTuple(i, &outputIdentifier);
  }

  this->printMsg("Randomizing Identifiers",1,t.getElapsedTime());

  return 1;
}
