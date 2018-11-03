#include <ttkOverlapTracking.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkOverlapTracking)

int ttkOverlapTracking::processTimestep(vtkDataObject* dataObject){

    auto block = vtkPointSet::SafeDownCast( dataObject );
    if(block==nullptr){
        dMsg(cerr, "[ttkOverlapTracking] ERROR: Input data can not be converted to vtkPointSet\n", memoryMsg);
        return 0;
    }

    auto labels = block->GetPointData()->GetAbstractArray("RegionId");
    if(labels==nullptr){
        dMsg(cout, "[ttkOverlapTracking] ERROR: No label point data found\n", timeMsg);
        return 0;
    }

    this->overlapTracking.processTimestep(
        (float*) block->GetPoints()->GetVoidPointer(0),
        (labelType*) block->GetPointData()->GetAbstractArray("RegionId")->GetVoidPointer(0),
        block->GetNumberOfPoints()
    );

    return 1;
}

int ttkOverlapTracking::finalize(vtkUnstructuredGrid* trackingGraph){
    auto& timeNodeLabelMap = this->overlapTracking.getTimeNodeLabelMap();
    auto& timeEdgesMap = this->overlapTracking.getTimeEdgesMap();
    size_t tn = timeNodeLabelMap.size();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Point Data
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(size_t t=0; t<tn; t++){
            auto& nodes = timeNodeLabelMap[t];
            for(size_t i=0; i<nodes.size(); i++){
                points->InsertNextPoint( t, i, 0 );
            }
        }

        mesh->SetPoints(points);
    }

    // Cell Data
    {
        vtkSmartPointer<vtkUnsignedLongLongArray> overlap = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        overlap->SetNumberOfComponents(1);
        overlap->SetName("Overlap");

        vtkSmartPointer<vtkIdTypeArray> r0 = vtkSmartPointer<vtkIdTypeArray>::New();
        r0->SetNumberOfComponents(1);
        r0->SetName("RegionId_t");

        vtkSmartPointer<vtkIdTypeArray> r1 = vtkSmartPointer<vtkIdTypeArray>::New();
        r1->SetNumberOfComponents(1);
        r1->SetName("RegionId_t+1");

        vector<size_t> offset(tn);
        offset[0] = 0;
        for(size_t t=1; t<tn; t++)
            offset[t] = offset[t-1] + timeNodeLabelMap[t-1].size();

        size_t nEdges = 0;
        for(size_t t=0; t<tn-1; t++){
            nEdges += timeEdgesMap[t].size();
        }

        for(size_t t=0; t<tn-1; t++){
            auto& edges = timeEdgesMap[t];
            for(auto& e: edges){
                vtkIdType ids[2] = { (vtkIdType)(e.i+offset[t]), (vtkIdType)(e.j+offset[t+1]) };
                mesh->InsertNextCell(VTK_LINE, 2, ids);
                r0->InsertNextValue( timeNodeLabelMap[t  ][e.i] );
                r1->InsertNextValue( timeNodeLabelMap[t+1][e.j] );
                overlap->InsertNextValue( e.overlap );
            }
        }
        auto cellData = mesh->GetCellData();
        cellData->AddArray( r0 );
        cellData->AddArray( r1 );
        cellData->AddArray( overlap );
    }

    trackingGraph->ShallowCopy( mesh );

    return 1;
}

int ttkOverlapTracking::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkOverlapTracking] RequestData" << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    double n = inMB->GetInformation()->Get( vtkDataObject::DATA_TIME_STEP() );

    // First Timestep
    if(i==0){
        this->overlapTracking.reset();
    }

    size_t nBlocks = inMB->GetNumberOfBlocks();

    // Process Timestep(s)
    for(size_t i=0; i<nBlocks; i++){
        int status = this->processTimestep( inMB->GetBlock(i) );
        if(status==0) return 0;
        if(nBlocks>1)
            dMsg(cout, "[ttkOverlapTracking] -----------------------------------------------------------\n", timeMsg);
    }

    // Last Timestep
    if(i==n-1 || nBlocks>1){
        if(nBlocks<2)
            dMsg(cout, "[ttkOverlapTracking] -----------------------------------------------------------\n", timeMsg);

        // Get Output
        vtkInformation* outInfo = outputVector->GetInformationObject(0);
        auto trackingGraph = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

        // Create Tracking Graph
        this->finalize( trackingGraph );

        // Print Status
        dMsg(cout, "[ttkOverlapTracking] Tracking Graph Finalized\n", timeMsg);
    }

    return 1;
}
