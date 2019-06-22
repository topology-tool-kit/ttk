#include <ttkPlanarGraphLayout.h>

#include <vtkAbstractArray.h>
#include <vtkFloatArray.h>
#include <vtkLongArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPlanarGraphLayout)

  int ttkPlanarGraphLayout::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  // Print status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkPlanarGraphLayout] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Set Wrapper
  planarGraphLayout.setWrapper(this);

  // Prepare input and output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto input = vtkUnstructuredGrid::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Copy input to output
  output->ShallowCopy(input);
  size_t nPoints = output->GetNumberOfPoints();
  size_t nEdges = output->GetNumberOfCells();

  // Check input fields
  auto outputPointData = output->GetPointData();

  auto getErrorMsg = [](string arrayType, string arrayName) {
    stringstream msg;
    msg << "[ttkPlanarGraphLayout] ERROR: Input point data does not have "
        << arrayType << " '" << arrayName << "'" << endl;
    return msg.str();
  };

  auto sequences
    = outputPointData->GetAbstractArray(this->GetSequenceFieldName().data());
  if(this->GetUseSequences() && !sequences) {
    dMsg(cout, getErrorMsg("array", this->GetSequenceFieldName()), fatalMsg);
    return 0;
  }

  auto sizes = vtkFloatArray::SafeDownCast(
    outputPointData->GetAbstractArray(this->GetSizeFieldName().data()));
  if(this->GetUseSizes() && !sizes) {
    dMsg(
      cout, getErrorMsg("vtkFloatArray", this->GetSizeFieldName()), fatalMsg);
    return 0;
  }

  auto branches
    = outputPointData->GetAbstractArray(this->GetBranchFieldName().data());
  if(this->GetUseBranches() && !branches) {
    dMsg(cout, getErrorMsg("vtkIdTypeArray", this->GetBranchFieldName()),
         fatalMsg);
    return 0;
  }

  auto levels
    = outputPointData->GetAbstractArray(this->GetLevelFieldName().data());
  if(this->GetUseLevels() && !levels) {
    dMsg(
      cout, getErrorMsg("vtkIdTypeArray", this->GetLevelFieldName()), fatalMsg);
    return 0;
  }

  // Initialize output field
  vtkSmartPointer<vtkFloatArray> outputField
    = vtkSmartPointer<vtkFloatArray>::New();
  outputField->SetName(this->GetOutputFieldName().data());
  outputField->SetNumberOfComponents(2); // (x,y) position
  outputField->SetNumberOfValues(nPoints * 2);

  auto sequenceType
    = this->GetUseSequences() ? sequences->GetDataType() : VTK_CHAR;
  auto branchType = this->GetUseBranches()
                      ? branches->GetDataType()
                      : this->GetUseLevels() ? levels->GetDataType() : VTK_CHAR;
  auto levelType
    = this->GetUseLevels()
        ? levels->GetDataType()
        : this->GetUseBranches() ? branches->GetDataType() : VTK_CHAR;

  if(branchType != levelType) {
    dMsg(cout,
         "[ttkPlanarGraphLayout] ERROR: Branch and Level array must have the "
         "same type.\n",
         fatalMsg);
    return 0;
  }

  // Compute layout with base code
  switch(vtkTemplate2PackMacro(branchType, sequenceType)) {
    ttkTemplate2Macro({
      int status
        = planarGraphLayout
            .execute<vtkIdType TTK_COMMA VTK_T1 TTK_COMMA VTK_T2>(
              // Input
              !this->GetUseSequences() ? nullptr
                                       : (VTK_T2 *)sequences->GetVoidPointer(0),
              !this->GetUseSizes() ? nullptr
                                   : (float *)sizes->GetVoidPointer(0),
              !this->GetUseBranches() ? nullptr
                                      : (VTK_T1 *)branches->GetVoidPointer(0),
              !this->GetUseLevels() ? nullptr
                                    : (VTK_T1 *)levels->GetVoidPointer(0),
              output->GetCells()->GetPointer(), nPoints, nEdges,

              // Output
              (float *)outputField->GetVoidPointer(0));

      if(status != 1)
        return 0;
    });
  }

  // Add output field to output
  outputPointData->AddArray(outputField);

  return 1;
}
