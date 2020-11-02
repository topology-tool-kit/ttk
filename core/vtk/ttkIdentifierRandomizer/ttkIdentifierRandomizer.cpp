#include <ttkIdentifierRandomizer.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkIdentifierRandomizer);

ttkIdentifierRandomizer::ttkIdentifierRandomizer() {

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->setDebugMsgPrefix("IdentifierRandomizer");
}

int ttkIdentifierRandomizer::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifierRandomizer::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkIdentifierRandomizer::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  ttk::Timer t;

  bool isPointData = false;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, inputVector);

  if(!inputScalarField) {
    printErr("Could not retrieve mandatory input array :(");
    return 0;
  }

  if(input->GetPointData()->GetArray(inputScalarField->GetName())
     == inputScalarField) {
    isPointData = true;
  }

  this->printMsg("Shuffling " + std::string{isPointData ? "vertex" : "cell"}
                 + " field `" + std::string{inputScalarField->GetName()}
                 + "'...");

  // allocate the memory for the output scalar field
  vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputScalarField->NewInstance());
  outputArray->SetName(inputScalarField->GetName());
  outputArray->SetNumberOfComponents(1);
  outputArray->SetNumberOfTuples(inputScalarField->GetNumberOfTuples());

  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> identifierMap;

  for(int i = 0; i < inputScalarField->GetNumberOfTuples(); i++) {
    double inputIdentifier = -1;
    inputScalarField->GetTuple(i, &inputIdentifier);

    bool isIn = false;
    for(size_t j = 0; j < identifierMap.size(); j++) {
      if(identifierMap[j].first == inputIdentifier) {
        isIn = true;
        break;
      }
    }

    if(!isIn) {
      identifierMap.emplace_back(inputIdentifier, -INT_MAX);
    }
  }

  // now let's shuffle things around
  ttk::SimplexId freeIdentifiers = identifierMap.size();
  ttk::SimplexId randomIdentifier = -1;

  for(size_t i = 0; i < identifierMap.size(); i++) {

    randomIdentifier = drand48() * (freeIdentifiers);

    ttk::SimplexId freeCounter = -1;
    for(size_t j = 0; j < identifierMap.size(); j++) {

      bool isFound = false;
      for(size_t k = 0; k < identifierMap.size(); k++) {

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
  for(int i = 0; i < inputScalarField->GetNumberOfTuples(); i++) {
    double inputIdentifier = -1;
    double outputIdentifier = -1;

    inputScalarField->GetTuple(i, &inputIdentifier);

    for(size_t j = 0; j < identifierMap.size(); j++) {
      if(inputIdentifier == identifierMap[j].first) {
        outputIdentifier = identifierMap[j].second;
        break;
      }
    }

    outputArray->SetTuple(i, &outputIdentifier);
  }

  if(isPointData)
    output->GetPointData()->AddArray(outputArray);
  else
    output->GetCellData()->AddArray(outputArray);

  printMsg("Processed " + std::to_string(outputArray->GetNumberOfTuples())
             + (isPointData ? " vertices." : " cells."),
           1, t.getElapsedTime(), 1);

  printMsg(ttk::debug::Separator::L1);

  return 1;
}
