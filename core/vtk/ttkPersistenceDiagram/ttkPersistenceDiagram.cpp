#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkPersistenceDiagram.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagram);

ttkPersistenceDiagram::ttkPersistenceDiagram() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistenceDiagram::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkPersistenceDiagram::dispatch(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  vtkDataArray *const inputScalarsArray,
  const scalarType *const inputScalars,
  scalarType *outputScalars,
  SimplexId *outputOffsets,
  int *outputMonotonyOffsets,
  const SimplexId *const inputOrder,
  const triangulationType *triangulation) {

  int status{};
  ttk::DiagramType CTDiagram{};

  if(BackEnd == BACKEND::APPROXIMATE_TOPOLOGY) {
    std::cout << "Chosen approx" << std::endl;
    double *range = inputScalarsArray->GetRange(0);
    this->setDeltaApproximate(range[1] - range[0]);
    this->setOutputScalars(outputScalars);
    this->setOutputOffsets(outputOffsets);
    this->setOutputMonotonyOffsets(outputMonotonyOffsets);
  }

  status = this->execute(CTDiagram, inputScalars, inputScalarsArray->GetMTime(),
                         inputOrder, triangulation);

  // something wrong in baseCode
  if(status != 0) {
    this->printErr("PersistenceDiagram::execute() error code: "
                   + std::to_string(status));
    return 0;
  }

  if(CTDiagram.empty()) {
    this->printErr("Empty diagram!");
    return 0;
  }

  vtkNew<vtkUnstructuredGrid> const vtu{};

  // convert CTDiagram to vtkUnstructuredGrid
  DiagramToVTU(vtu, CTDiagram, inputScalarsArray, *this,
               triangulation->getDimensionality(), this->ShowInsideDomain);

  outputCTPersistenceDiagram->ShallowCopy(vtu);

  if(this->ClearDGCache && this->BackEnd == BACKEND::DISCRETE_MORSE_SANDWICH) {
    this->printMsg("Clearing DiscreteGradient cache...");
    ttk::dcg::DiscreteGradient::clearCache(*triangulation);
  }

  return 1;
}

int ttkPersistenceDiagram::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputCTPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars) {
    this->printErr("Wrong input scalars");
    return 0;
  }
#endif

  vtkDataArray *offsetField
    = this->GetOrderArray(input, 0, 1, ForceInputOffsetScalarField);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!offsetField) {
    this->printErr("Wrong input offsets");
    return 0;
  }
  if(offsetField->GetDataType() != VTK_INT
     and offsetField->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported");
    return 0;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  outputOffsets->SetNumberOfComponents(1);
  outputOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputOffsets->SetName("outputOffsets");

  vtkNew<vtkIntArray> outputMonotonyOffsets{};
  outputMonotonyOffsets->SetNumberOfComponents(1);
  outputMonotonyOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputMonotonyOffsets->SetName("outputMonotonyffsets");
  outputMonotonyOffsets->FillComponent(0, 0);

  vtkSmartPointer<vtkDataArray> const outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->SetNumberOfComponents(1);
  outputScalars->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputScalars->DeepCopy(inputScalars);
  outputScalars->SetName("Cropped");

  int status{};
  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    status = this->dispatch(
      outputCTPersistenceDiagram, inputScalars,
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
      static_cast<int *>(ttkUtils::GetVoidPointer(outputMonotonyOffsets)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsetField)),
      static_cast<TTK_TT *>(triangulation->getData())));

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  return status;
}
