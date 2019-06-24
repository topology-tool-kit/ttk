#include <ttkEndFor.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkEndFor)

  int ttkEndFor::RequestInformation(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // Reset index
  this->nextIndex = 0;
  return this->Superclass::RequestInformation(
    request, inputVector, outputVector);
}

int ttkEndFor::RequestUpdateExtent(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {

  // Request next index
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(
    vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->nextIndex);

  return 1;
}

int ttkEndFor::RequestData(vtkInformation *request,
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {
  // Print status
  string divider = "==========================================================="
                   "=====================";
  dMsg(cout, divider + "\n", infoMsg);

  // Get Input
  vtkInformation *inInfo0 = inputVector[1]->GetInformationObject(0);
  auto inputFor = inInfo0->Get(vtkDataObject::DATA_OBJECT());

  // Get iteration information
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    inputFor->GetFieldData()->GetAbstractArray("_ttk_IterationInfo"));
  if(!iterationInformation) {
    dMsg(cout,
         "[ttkEndFor] ERROR: Unable to retrieve iteration information (is the "
         "input the start of the for loop?).\n",
         fatalMsg);
    return 0;
  }
  this->nextIndex = iterationInformation->GetValue(0) + 1;
  double n = iterationInformation->GetValue(1);

  if(this->nextIndex < n) {
    // Request Next Element
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);

    // Print status
    {
      stringstream msg;
      msg << "[ttkEndFor] Next Iteration: " << this->nextIndex << "   ";
      size_t padd = divider.length() - msg.str().length();
      for(size_t i = 0; i < padd; i++)
        msg << "\\";
      msg << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  } else {
    // Stop iterations
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    this->nextIndex = 0;

    // Print status
    {
      stringstream msg;
      msg << "[ttkEndFor] Iteration complete" << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    // Copy input to output
    vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);
    auto input = inInfo1->Get(vtkDataObject::DATA_OBJECT());

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
    output->ShallowCopy(input);

    output->GetFieldData()->RemoveArray("_ttk_IterationInfo");
  }

  return 1;
}
