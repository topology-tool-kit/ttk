#include <ttkPeriodicGrid.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <Triangulation.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPeriodicGrid);

ttkPeriodicGrid::ttkPeriodicGrid() {
  this->setDebugMsgPrefix("PeriodicGrid");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPeriodicGrid::~ttkPeriodicGrid() {
}

int ttkPeriodicGrid::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  else
    return 0;

  return 1;
}

int ttkPeriodicGrid::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkPeriodicGrid::RequestData(vtkInformation *ttkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);

  // Get ttk::triangulation of the input vtkDataSet (will create one if does
  // not exist already)
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);

  bool before = triangulation->hasPeriodicBoundaries();

  triangulation->setPeriodicBoundaryConditions(Periodicity);

  printMsg("Switching regular grid periodicity from "
           + (before ? std::string("ON") : std::string("OFF")) + " to "
           + (Periodicity ? std::string("ON") : std::string("OFF")));

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  return 1;
}
