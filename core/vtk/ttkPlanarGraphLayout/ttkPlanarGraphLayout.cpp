#include <ttkPlanarGraphLayout.h>

#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPlanarGraphLayout)

int ttkPlanarGraphLayout::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkPlanarGraphLayout] RequestData"<<endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Set Wrapper
    planarGraphLayout.setWrapper(this);

    // Prepare input and output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = vtkUnstructuredGrid::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Copy input to output
    output->ShallowCopy(input);

    size_t nPoints = output->GetNumberOfPoints();
    size_t nEdges = output->GetNumberOfCells();

    // Check input fields
    auto outputPointData = output->GetPointData();

    auto sequences = outputPointData->GetAbstractArray( this->GetSequenceFieldName().data() );
    if(this->GetUsePointSequence() && !sequences){
        stringstream msg;
        msg<<"[ttkPlanarGraphLayout] ERROR: Input point data does not have array '" << this->GetSequenceFieldName() << "'" <<endl;
        dMsg(cout, msg.str(), fatalMsg);
        return 0;
    }

    auto sizes = outputPointData->GetAbstractArray( this->GetSizeFieldName().data() );
    if(this->GetUsePointSize() && !sizes){
        stringstream msg;
        msg<<"[ttkPlanarGraphLayout] ERROR: Input point data does not have array '" << this->GetSizeFieldName() << "'" <<endl;
        dMsg(cout, msg.str(), fatalMsg);
        return 0;
    }

    auto branches = outputPointData->GetAbstractArray( this->GetBranchFieldName().data() );
    if(this->GetUsePointBranch() && !branches){
        stringstream msg;
        msg<<"[ttkPlanarGraphLayout] ERROR: Input point data does not have array '" << this->GetSizeFieldName() << "'" <<endl;
        dMsg(cout, msg.str(), fatalMsg);
        return 0;
    }

    auto levels = outputPointData->GetAbstractArray( this->GetLevelFieldName().data() );
    if(this->GetUsePointLevel() && !levels){
        stringstream msg;
        msg<<"[ttkPlanarGraphLayout] ERROR: Input point data does not have array '" << this->GetLevelFieldName() << "'" <<endl;
        dMsg(cout, msg.str(), fatalMsg);
        return 0;
    }

    // Initialize output field
    vtkSmartPointer<vtkFloatArray> outputField = vtkSmartPointer<vtkFloatArray>::New();
    outputField->SetName( this->GetOutputFieldName().data() );
    outputField->SetNumberOfComponents(1);
    outputField->SetNumberOfValues(nPoints);

    // Compute layout
    switch(sequences->GetDataType()){
        ttkTemplateMacro({
            int status = planarGraphLayout.execute<vtkIdType TTK_COMMA VTK_TT>(
                // Input
                (VTK_TT*) sequences->GetVoidPointer(0),
                !this->GetUsePointSize()    ? nullptr : (float*)     sizes->GetVoidPointer(0),
                !this->GetUsePointBranch()  ? nullptr : (vtkIdType*) branches->GetVoidPointer(0),
                !this->GetUsePointLevel()   ? nullptr : (vtkIdType*) levels->GetVoidPointer(0),
                output->GetCells()->GetPointer(),
                nPoints,
                nEdges,

                // Output
                (float*) outputField->GetVoidPointer(0)
            );

            if(status!=1) return 0;
        });
    }

    // Add output field to output
    outputPointData->AddArray( outputField );

    return 1;
}
