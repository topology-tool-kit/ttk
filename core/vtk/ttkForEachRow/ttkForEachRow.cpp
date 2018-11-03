#include <ttkForEachRow.h>

#include <vtkTable.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkForEachRow)

int ttkForEachRow::RequestInformation(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    // These values need to exists to automatically enable temporal streaming
    double dummy[2] = {0, 1};
    outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_STEPS(), dummy, 2);
    outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_RANGE(), dummy, 2);

    return this->Superclass::RequestInformation(request, inputVector, outputVector);
}

int ttkForEachRow::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Get current row index
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    double index = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

    // Print status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkForEachRow]  Iteration: " << index << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Get Input and Output
    auto inputTable = vtkTable::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );
    auto outputTable = vtkTable::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Propagate number of rows downstream
    outputTable->GetInformation()->Set( vtkDataObject::DATA_TIME_STEP(), inputTable->GetNumberOfRows() );

    // Extract row at index
    size_t n = inputTable->GetNumberOfColumns();
    for(size_t i=0; i<n; i++){
        auto column = inputTable->GetColumn(i);
        auto newColumn = vtkAbstractArray::CreateArray( column->GetDataType() );
        newColumn->SetName( column->GetName() );
        newColumn->SetNumberOfValues( 1 );
        newColumn->SetVariantValue(0, column->GetVariantValue( index ) );
        outputTable->AddColumn( newColumn );
    }

    return 1;
}