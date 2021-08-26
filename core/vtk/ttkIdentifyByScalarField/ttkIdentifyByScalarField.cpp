#include <numeric>

#include <ttkIdentifyByScalarField.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkIdentifyByScalarField);

ttkIdentifyByScalarField::ttkIdentifyByScalarField() {
  this->setDebugMsgPrefix("IdentifyByScalarField");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkIdentifyByScalarField::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }

  return 0;
}

int ttkIdentifyByScalarField::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }

  return 0;
}

template <typename VTK_TT>
int ttkIdentifyByScalarField::dispatch(vector<SimplexId> &inputIds) {
  VTK_TT *arr = static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars_));

  auto greater_cmp = [=](int a, int b) { return arr[a] > arr[b]; };
  auto lower_cmp = [=](int a, int b) { return arr[a] < arr[b]; };

  if(IncreasingOrder) {
    std::sort(inputIds.begin(), inputIds.end(), lower_cmp);
  } else {
    std::sort(inputIds.begin(), inputIds.end(), greater_cmp);
  }

  return 0;
}

int ttkIdentifyByScalarField::RequestData(vtkInformation *ttkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  inputScalars_ = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalars_)
    return 0;

  int inputArrayAssociation = this->GetInputArrayAssociation(0, inputVector);

  if(inputArrayAssociation > 1 || inputArrayAssociation < 0) {
    printErr("input array has to be cell data or point data.");
    return 0;
  }

  ttk::Timer t;
  this->printMsg("Computing Identifiers", 0, t.getElapsedTime(), 1,
                 ttk::debug::LineMode::REPLACE);

  const SimplexId numberOfValues = inputArrayAssociation == 0
                                     ? input->GetNumberOfPoints()
                                     : input->GetNumberOfCells();

  vector<SimplexId> inputIds(numberOfValues);
  std::iota(inputIds.begin(), inputIds.end(), 0);
  switch(inputScalars_->GetDataType()) {
    vtkTemplateMacro(dispatch<VTK_TT>(inputIds));
  }

  this->printMsg("Computing Identifiers", 1, t.getElapsedTime(), 1);

  vtkSmartPointer<ttkSimplexIdTypeArray> ids
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  ids->SetNumberOfComponents(1);
  ids->SetNumberOfTuples(numberOfValues);
  ids->SetName(inputArrayAssociation == 0 ? "PointScalarFieldName"
                                          : "CellScalarFieldName");

  SimplexId *outputIds
    = static_cast<SimplexId *>(ttkUtils::GetVoidPointer(ids));

  for(int i = 0; i < numberOfValues; ++i)
    outputIds[inputIds[i]] = i;
  if(StartByOne) {
    for(int i = 0; i < numberOfValues; ++i)
      outputIds[i] += 1;
  }

  output->ShallowCopy(input);

  if(inputArrayAssociation == 0) {
    output->GetPointData()->AddArray(ids);
  } else {
    output->GetCellData()->AddArray(ids);
  }

  return 1;
}
