#include <ttkIdentifyByScalarField.h>

#include <ttkUtils.h>
#include <ttkMacros.h>

#include <numeric>
#include <vtkCellData.h>

vtkStandardNewMacro(ttkIdentifyByScalarField);

ttkIdentifyByScalarField::ttkIdentifyByScalarField() {
  this->setDebugMsgPrefix("IdentifyByScalarField");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkIdentifyByScalarField::~ttkIdentifyByScalarField() {
}

int ttkIdentifyByScalarField::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkIdentifyByScalarField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
int dispatch(
    std::vector<ttk::SimplexId> &inputIds,
    vtkDataArray* inputArray,
    bool increasingOrder
) {
  auto arrayData = (VTK_TT*) ttkUtils::GetVoidPointer(inputArray);

  auto greater_cmp = [=](int a, int b) { return arrayData[a] > arrayData[b]; };
  auto lower_cmp = [=](int a, int b) { return arrayData[a] < arrayData[b]; };

  if(increasingOrder) {
    std::sort(inputIds.begin(), inputIds.end(), lower_cmp);
  } else {
    std::sort(inputIds.begin(), inputIds.end(), greater_cmp);
  }

  return 1;
}

int ttkIdentifyByScalarField::RequestData(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector
) {
  ttk::Timer t;

  this->printMsg(
    "Computing Identifiers",
    0, 0,
    ttk::debug::LineMode::REPLACE
  );

  auto input = vtkDataSet::GetData( inputVector[0] );
  auto output = vtkDataSet::GetData( outputVector );

  auto inputArray = this->GetInputArrayToProcess(0, inputVector);

  const ttk::SimplexId numberOfCells = input->GetNumberOfCells();
  std::vector<ttk::SimplexId> inputIds(numberOfCells);
  std::iota(inputIds.begin(), inputIds.end(), 0);
  switch(inputArray->GetDataType()) {
    vtkTemplateMacro(dispatch<VTK_TT>(inputIds,inputArray,IncreasingOrder));
  }

  auto ids = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  ids->SetNumberOfComponents(1);
  ids->SetNumberOfTuples(numberOfCells);
  ids->SetName("CellScalarFieldName");
  auto outputIds = (ttk::SimplexId*) ttkUtils::GetVoidPointer(ids);

  for(int i = 0; i < numberOfCells; ++i)
    outputIds[inputIds[i]] = i;
  if(StartByOne) {
    for(int i = 0; i < numberOfCells; ++i)
      outputIds[i] += 1;
  }

  output->ShallowCopy(input);
  output->GetCellData()->AddArray(ids);

  this->printMsg(
    "Computing Identifiers",
    1, t.getElapsedTime()
  );

  return 1;
}