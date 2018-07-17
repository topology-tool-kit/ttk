#include                  <ttkTrackingFromFields.h>

vtkStandardNewMacro(ttkTrackingFromFields)

constexpr unsigned int str2int(const char* str, int h = 0)
{
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

int ttkTrackingFromFields::FillOutputPortInformation(int port, vtkInformation* info)
{
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    // default:
    //   break;
  }

  return 1;
}

int ttkTrackingFromFields::doIt(
  std::vector<vtkDataSet *> &inputs,
  std::vector<vtkDataSet *> &outputs)
{
  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  triangulation_.setWrapper(this);
  triangulation_.setInputData(input);
  internalTriangulation_ = ttkTriangulation::getTriangulation(input);
  internalTriangulation_->preprocessVertexNeighbors();

  // Test validity of datasets (must present the same number of points).
  if (!input)
    return -1;

  // Get number and list of inputs.
  std::vector<vtkDataArray*> inputScalarFieldsRaw;
  std::vector<vtkDataArray*> inputScalarFields;
  int numberOfInputFields = input->GetPointData()->GetNumberOfArrays();
  if (numberOfInputFields < 3) {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] not enough input fields to perform tracking." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  vtkDataArray* firstScalarField = input->GetPointData()->GetArray(0);
  // inputScalarFields.push_back(firstScalarField);

  for (int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray* currentScalarField = input->GetPointData()->GetArray(i);
    if (!currentScalarField ||
        firstScalarField->GetDataType() != currentScalarField->GetDataType())
    {
      std::stringstream msg;
      msg << "[ttkTrackingFromField] inconsistent field data type or size (" << i << ")." << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(
      inputScalarFieldsRaw.begin(),
      inputScalarFieldsRaw.end(),
      [] (vtkDataArray* a,
          vtkDataArray* b)
      {
          std::string s1 = a->GetName();
          std::string s2 = b->GetName();
          return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
      });

  int end = EndTimestep <= 0 ? numberOfInputFields : std::min(numberOfInputFields, EndTimestep);
  for (int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray* currentScalarField = inputScalarFieldsRaw[i];
    std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if (pvalg >= 0) {
    switch (pvalg) {
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
        msg << "[ttkTrackingFromFieldsFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  } else {
    switch (str2int(algorithm.c_str())) {
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
        msg << "[ttkTrackingFromField] Unrecognized tracking method." << std::endl;
        dMsg(std::cout, msg.str(), timeMsg);
        break;
    }
  }

  int res = 0;
  if (useTTKMethod)
    res = trackWithPersistenceMatching(
        input, output, inputScalarFields);
  else
  {
    std::stringstream msg;
    msg << "[ttkTrackingFromField] The specified matching method does not." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  return res;
}

// (*) Persistence-driven approach
int ttkTrackingFromFields::trackWithPersistenceMatching(
    vtkDataSet *input,
    vtkUnstructuredGrid *output,
    std::vector<vtkDataArray*> inputScalarFields)
{
  unsigned long fieldNumber = inputScalarFields.size();

  // 1. get persistence diagrams.
  std::vector<vtkUnstructuredGrid*> persistenceDiagrams(fieldNumber);
  // #pragma omp parallel for num_threads(ThreadNumber)
  for (int i = 0; i < (int) fieldNumber; ++i)
  {
    ttkPersistenceDiagram *ttkPD = ttkPersistenceDiagram::New();

    ttkPD->AddInputData(input);
    ttkPD->SetUseInputOffsetScalarField(false);
    ttkPD->SetComputeSaddleConnectors(false);
    ttkPD->SetShowInsideDomain(true);
    ttkPD->SetThreadNumber(1);
    // ttkPD->SetScalarFieldId(i);
    ttkPD->SetScalarField(inputScalarFields[i]->GetName());
    ttkPD->Update();

    vtkUnstructuredGrid *pd =
        vtkUnstructuredGrid::SafeDownCast(ttkPD->GetOutput(0));

    persistenceDiagrams[i] = pd;
  }

  // 2. call feature tracking with threshold.
  ttkTrackingFromPersistenceDiagrams *ttkFT = ttkTrackingFromPersistenceDiagrams::New();
  for (int i = 0; i < (int) fieldNumber; ++i) {
    vtkUnstructuredGrid *pd = persistenceDiagrams[i];
    vtkIntArray *via = vtkIntArray::New();
    via->SetNumberOfComponents(1);
    via->SetNumberOfTuples(1);
    via->SetTuple1(0, i);
    pd->GetFieldData()->AddArray(via);
    ttkFT->AddInputData(pd);
  }

  ttkFT->SetTolerance(Tolerance);
  ttkFT->SetPX(PX);
  ttkFT->SetPY(PY);
  ttkFT->SetPZ(PZ);
  ttkFT->SetPE(PE);
  ttkFT->SetPS(PS);
  ttkFT->SetDoPostProc(DoPostProc);
  ttkFT->SetPostProcThresh(PostProcThresh);
  ttkFT->SetThreadNumber(ThreadNumber);
  ttkFT->SetDistanceAlgorithm(DistanceAlgorithm);
  ttkFT->SetPVAlgorithm(PVAlgorithm);
  ttkFT->SetAlpha(Alpha);
  ttkFT->SetSpacing(Spacing);
  ttkFT->SetWassersteinMetric(WassersteinMetric);
  ttkFT->SetUseGeometricSpacing(UseGeometricSpacing);
  ttkFT->SetIs3D(Is3D);
  ttkFT->Update();

  vtkUnstructuredGrid *oo = vtkUnstructuredGrid::SafeDownCast(ttkFT->GetOutput(0));
  outputMesh_ = vtkUnstructuredGrid::New();
  outputMesh_->ShallowCopy(oo);
  output->ShallowCopy(outputMesh_);

  return 0;
}
