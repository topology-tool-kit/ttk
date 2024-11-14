#include <ttkMorseSmallComplexStability.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMorseSmallComplexStability);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMorseSmallComplexStability::ttkMorseSmallComplexStability() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMorseSmallComplexStability::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMorseSmallComplexStability::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkMorseSmallComplexStability::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  
  // Loop through all vtkInformation objects in the vtkInformationVector
  for (int i = 0; i < (*inputVector)->GetNumberOfInformationObjects(); ++i)
  {
      // Get the vtkInformation object at the current index
      vtkInformation* info = (*inputVector)->GetInformationObject(i);
      // Print all key-value pairs in the vtkInformation object
      info->Print(cout);
  }

  // Get input data (pointer to vtkDataObject, parent block id)
  //std::vector<std::pair<vtkDataObject *, size_t>> blocks{};
//
 
  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkMultiBlockDataSet* multiBlockInput = vtkMultiBlockDataSet::SafeDownCast(input);
  if(multiBlockInput == nullptr){
    this->printErr("Only one dataset");
    return -1;
  }

  int n_blocks = multiBlockInput->GetNumberOfBlocks();

  std::vector<vtkDataSet*> inputs;
  for (int i = 0; i < n_blocks ; i++){
    vtkDataSet* newInput = multiBlockInput.GetBlock(i);
    inputs.push_back(newInput);
  } 

  
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  vtkSmartPointer<vtkDataArray> const outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                      (status = this->computeAverages<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         (TTK_TT *)triangulation->getData())));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}
