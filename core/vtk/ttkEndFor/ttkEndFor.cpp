#include <ttkEndFor.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkEndFor)

int ttkEndFor::RequestInformation(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Reset index
    this->nextIndex = 0;
    return this->Superclass::RequestInformation(request, inputVector, outputVector);
}

int ttkEndFor::RequestUpdateExtent(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    stringstream msg;
    msg << "[ttkEndFor] Next Iteration: "<< this->nextIndex << endl;
    dMsg(cout, msg.str(), infoMsg);

    // Request next index
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    inInfo->Set( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->nextIndex);

    return 1;
}

int ttkEndFor::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Get Input
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = inInfo->Get( vtkDataObject::DATA_OBJECT() );

    // Get iteration information
    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    double n = input->GetInformation()->Get( vtkDataObject::DATA_TIME_STEP() );

    this->nextIndex = i+1;

    if(this->nextIndex<n){
        // Request Next Element
        request->Set( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );
    } else {
        // Stop iterations
        request->Remove( vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING() );

        // Print status
        {
            stringstream msg;
            msg << "[ttkEndFor] Iteration complete" << endl;
            dMsg(cout, msg.str(), infoMsg);
        }

        // Copy input to output
        vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
        auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());

        vtkInformation* outInfo = outputVector->GetInformationObject(0);
        auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
        output->ShallowCopy(input);
    }

    return 1;
}