#include <ttkTrackingFromFields.h>

vtkStandardNewMacro(ttkTrackingFromFields)

  constexpr unsigned int str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

int ttkTrackingFromFields::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    default:
      break;
  }

  return 1;
}

int ttkTrackingFromFields::doIt(std::vector<vtkDataSet *> &inputs,
                                std::vector<vtkDataSet *> &outputs) {
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  internalTriangulation_ = ttkTriangulation::getTriangulation(input);
  internalTriangulation_->setWrapper(this);

  // Test validity of datasets (must present the same number of points).
  if(!input)
    return -1;

  // Get number and list of inputs.
  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  int numberOfInputFields = input->GetPointData()->GetNumberOfArrays();
  if(numberOfInputFields < 3) {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] not enough input fields to perform tracking."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  vtkDataArray *firstScalarField = input->GetPointData()->GetArray(0);

  for(int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray *currentScalarField = input->GetPointData()->GetArray(i);
    if(!currentScalarField
       || firstScalarField->GetDataType()
            != currentScalarField->GetDataType()) {
      std::stringstream msg;
      msg << "[ttkTrackingFromField] inconsistent field data type or size ("
          << i << ")." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(inputScalarFieldsRaw.begin(), inputScalarFieldsRaw.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              std::string s1 = a->GetName();
              std::string s2 = b->GetName();
              return std::lexicographical_compare(
                s1.begin(), s1.end(), s2.begin(), s2.end());
            });

  int end = EndTimestep <= 0 ? numberOfInputFields
                             : std::min(numberOfInputFields, EndTimestep);
  for(int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if(pvalg >= 0) {
    switch(pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromFieldsFromField] Unrecognized tracking method."
            << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  } else {
    switch(str2int(algorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
      case str2int("1"):
      case str2int("legacy"):
      case str2int("2"):
      case str2int("geometric"):
      case str2int("3"):
      case str2int("parallel"):
        useTTKMethod = true;
        break;
      case str2int("4"):
      case str2int("greedy"):
        break;
      default:
        std::stringstream msg;
        msg << "[ttkTrackingFromField] Unrecognized tracking method."
            << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  }

  int res = 0;
  if(useTTKMethod)
    res
      = trackWithPersistenceMatching<double>(input, output, inputScalarFields);
  else {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] The specified matching method does not."
        << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return res;
}
