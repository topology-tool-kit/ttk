#include <ttkForEachRow.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>
#include <vtkUnsignedLongLongArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkForEachRow)

  int ttkForEachRow::RequestInformation(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // These values need to exists to automatically enable temporal streaming
  double dummy[2] = {0, 1};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), dummy, 2);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), dummy, 2);

  return this->Superclass::RequestInformation(
    request, inputVector, outputVector);
}

int ttkForEachRow::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  // Get current row index
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  double index
    = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  // Print status
  {
    string divider = "========================================================="
                     "=======================";

    stringstream msg;
    msg << "[ttkForEachRow]  Iteration: " << index << "   ";
    size_t n = divider.length() - msg.str().length();
    for(size_t i = 0; i < n; i++)
      msg << "/";
    msg << endl;
    dMsg(cout, divider + "\n" + msg.str(), infoMsg);
  }

  // Get Input and Output
  auto inputTable
    = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  auto outputTable
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Propagate number of rows downstream
  auto iterationInformation = vtkSmartPointer<vtkDoubleArray>::New();
  iterationInformation->SetName("_ttk_IterationInfo");
  iterationInformation->SetNumberOfValues(2);
  iterationInformation->SetValue(0, index);
  iterationInformation->SetValue(1, inputTable->GetNumberOfRows());
  outputTable->GetFieldData()->AddArray(iterationInformation);

  outputTable->GetInformation()->Set(
    vtkDataObject::DATA_TIME_STEP(), inputTable->GetNumberOfRows());

  // Extract row at index
  size_t n = inputTable->GetNumberOfColumns();
  for(size_t i = 0; i < n; i++) {
    auto column = inputTable->GetColumn(i);
    auto newColumn
      = vtkSmartPointer<vtkAbstractArray>::Take(column->NewInstance());
    newColumn->SetName(column->GetName());
    newColumn->SetNumberOfTuples(1);
    newColumn->SetVariantValue(0, column->GetVariantValue(index));
    outputTable->AddColumn(newColumn);
  }

  return 1;
}