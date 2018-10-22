#include <ttkOverlapTracking.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkOverlapTracking)

int ttkOverlapTracking::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkOverlapTracking] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Get Output Graph
    // vtkInformation* outInfo = outputVector->GetInformationObject(0);
    // vtkUnstructuredGrid* outTable = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* outMB = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    outMB->ShallowCopy(inMB);

    vtkPointSet* block1 = vtkPointSet::SafeDownCast( inMB->GetBlock(0) );
    vtkPointSet* block2 = vtkPointSet::SafeDownCast( inMB->GetBlock(1) );

    overlapTracking.reset();

    {
        auto block = block1;
        overlapTracking.processTimestep(
            (float*) block->GetPoints()->GetVoidPointer(0),
            (labelType*) block->GetPointData()->GetAbstractArray("RegionId")->GetVoidPointer(0),
            block->GetNumberOfPoints()
        );
    }
    {
        auto block = block2;

        // cout<< block->GetPointData()->GetAbstractArray("RegionId")->GetTypeAsString() <<endl;
        block->GetPointData()->GetAbstractArray("RegionId")->Print(cout);

        auto xxx = vtkIdTypeArray::SafeDownCast( block->GetPointData()->GetAbstractArray("RegionId") );
        cout<<"... "<<xxx->GetSize()<<" "<<block->GetNumberOfPoints()<<endl;

        size_t count = 0;
        for(size_t i=0; i<xxx->GetSize(); i++)
            if(xxx->GetValue(i)==0) count++;
        cout<<count<<endl;

        overlapTracking.processTimestep(
            (float*) block->GetPoints()->GetVoidPointer(0),
            (labelType*) block->GetPointData()->GetAbstractArray("RegionId")->GetVoidPointer(0),
            block->GetNumberOfPoints()
        );
    }

    // vtkSmartPointer<vtkDoubleArray> indices = vtkSmartPointer<vtkDoubleArray>::New();
    // indices->SetNumberOfComponents(1);
    // indices->SetNumberOfValues( n );
    // indices->SetName( "Indices" );
    // pd->AddArray( indices );

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkOverlapTracking] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
