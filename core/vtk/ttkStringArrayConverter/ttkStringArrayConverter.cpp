#include <ttkStringArrayConverter.h>

#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>

#include <ttkUtils.h>

#include <map>
#include <set>

vtkStandardNewMacro(ttkStringArrayConverter);

ttkStringArrayConverter::ttkStringArrayConverter() {
  this->setDebugMsgPrefix("StringArrayConverter");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkStringArrayConverter::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkStringArrayConverter::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkStringArrayConverter::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {
  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    this->printErr("Null input data, aborting");
    return 0;
  }

  // string array
  const auto sa = vtkStringArray::SafeDownCast(
    this->GetInputAbstractArrayToProcess(0, inputVector));

  if(sa == nullptr) {
    this->printErr("Cannot find the required string array");
    return 0;
  }

  const auto nvalues = sa->GetNumberOfTuples();

  std::set<std::string> values{};

  for(vtkIdType i = 0; i < nvalues; ++i) {
    values.emplace(sa->GetValue(i));
  }

  std::map<std::string, size_t> valInd{};

  {
    size_t i = 0;
    for(const auto &el : values) {
      valInd[el] = i;
      i++;
    }
  }

  // shallow-copy input
  output->ShallowCopy(input);

  const auto pdo = output->GetPointData();

  // array of indices
  vtkNew<vtkIdTypeArray> ia{};
  // array name
  std::string colname{sa->GetName()};
  colname += "Int";
  ia->SetName(colname.data());
  ia->SetNumberOfTuples(nvalues);
  // fill array with corresponding string values indices
  for(vtkIdType i = 0; i < nvalues; ++i) {
    ia->SetValue(i, valInd[sa->GetValue(i)]);
  }

  pdo->AddArray(ia);

  return 1;
}
