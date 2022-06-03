#include <ttkPathCompression.h>

#include <vtkInformation.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPathCompression);

ttkPathCompression::ttkPathCompression() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPathCompression::~ttkPathCompression() = default {
}

int ttkPathCompression::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPathCompression::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }

  return 0;
}

int ttkPathCompression::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  auto order = ttkAlgorithm::GetOrderArray(inputDataSet, 0);
  if(!order) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }

  if(order->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation!");
  this->printMsg("  Scalar Array: " + std::string(order->GetName()));

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation) {
    return 0;
  }

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  vtkSmartPointer<vtkDataArray> descendingManifold
    = vtkSmartPointer<vtkDataArray>::Take(order->NewInstance());
  descendingManifold->SetName("DescendingManifold"); // set array name
  descendingManifold->SetNumberOfComponents(1); // only one component per tuple
  descendingManifold->SetNumberOfTuples(order->GetNumberOfTuples());

  vtkSmartPointer<vtkDataArray> ascendingManifold
    = vtkSmartPointer<vtkDataArray>::Take(order->NewInstance());
  ascendingManifold->SetName("AscendingManifold"); // set array name
  ascendingManifold->SetNumberOfComponents(1); // only one component per tuple
  ascendingManifold->SetNumberOfTuples(order->GetNumberOfTuples());

  this->printMsg("  Output Array 1: "
                 + std::string(descendingManifold->GetName()));
  this->printMsg("  Output Array 2: "
                 + std::string(ascendingManifold->GetName()));
  int status = 0; // this integer checks if the base code returns an error

#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    auto pointData = inputDataSet->GetPointData();
    auto rankArray = pointData->GetArray("RankArray");
    auto globalIds = pointData->GetGlobalIds();

    // Templatize over the different input array data types and call the base
    // code
    ttkTypeMacroIT(
      order->GetDataType(), triangulation->getType(),
      (status = this->computeCompression<T0, T1>(
         ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
         ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
         ttkUtils::GetPointer<T0>(order), (T1 *)triangulation->getData(),
         ttkUtils::GetPointer<int>(rankArray),
         ttkUtils::GetPointer<ttk::SimplexId>(globalIds))));
  }

#else
  ttkTypeMacroIT(
    order->GetDataType(), triangulation->getType(),
    (status = this->computeCompression<T0, T1>(
       ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
       ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
       ttkUtils::GetPointer<T0>(order), (T1 *)triangulation->getData())));
#endif // TTK_ENABLE_MPI

  // On error cancel filter execution
  if(status != 1)
    return 0;

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(descendingManifold);
  outputDataSet->GetPointData()->AddArray(ascendingManifold);

  // return success
  return 1;
}
