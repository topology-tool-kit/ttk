#include <ttkTopologicalSimplificationByPersistence.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTopologicalSimplificationByPersistence);

ttkTopologicalSimplificationByPersistence::
  ttkTopologicalSimplificationByPersistence() {
  this->setDebugMsgPrefix("PLTS");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkTopologicalSimplificationByPersistence::
  ~ttkTopologicalSimplificationByPersistence() {
}

int ttkTopologicalSimplificationByPersistence::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTopologicalSimplificationByPersistence::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTopologicalSimplificationByPersistence::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // retrieve input data object
  auto inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  // retrieve input arrays
  if(this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Input array needs to be a point data array.");

  auto inputScalars = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalars)
    return !this->printErr("Unable to retrieve input array.");

  if(inputScalars->GetNumberOfComponents() != 1)
    return !this->printErr("Input array needs to be a scalar array.");

  auto inputOrder = this->GetOrderArray(inputDataSet, 0);
  if(!inputOrder)
    return 0;

  // create output arrays
  auto outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->DeepCopy(inputScalars);
  auto outputOrder
    = vtkSmartPointer<vtkDataArray>::Take(inputOrder->NewInstance());
  outputOrder->DeepCopy(inputOrder);

  // retrieve and precondition triangulation
  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;
  this->preconditionTriangulation(triangulation);

  double persistenceThreshold = this->PersistenceThreshold;
  if(!this->ThresholdIsAbsolute) {
    double range[2];
    outputScalars->GetRange(range);
    persistenceThreshold = (range[1] - range[0]) * persistenceThreshold;
  }

  // perform simplification
  int status = 0;
  ttkVtkTemplateMacro(
    outputScalars->GetDataType(), triangulation->getType(),
    (status = this->removeNonPersistentExtrema<VTK_TT, ttk::SimplexId, TTK_TT>(
       ttkUtils::GetPointer<VTK_TT>(outputScalars),
       ttkUtils::GetPointer<ttk::SimplexId>(outputOrder),

       static_cast<TTK_TT *>(triangulation->getData()), persistenceThreshold,
       this->ComputePerturbation, this->PairType)));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  auto outputDataSet = vtkDataSet::GetData(outputVector, 0);
  outputDataSet->ShallowCopy(inputDataSet);

  auto outputDataSetPD = outputDataSet->GetPointData();
  outputDataSetPD->AddArray(outputScalars);
  outputDataSetPD->AddArray(outputOrder);

  return 1;
}
