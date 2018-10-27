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
        stringstream msg;
        msg << "[ttkOverlapTracking] ERROR: Input data can not be converted to vtkPointSet" << endl;
        dMsg(cerr, msg.str(), memoryMsg);
        return 0;
    }

    auto labels = block->GetPointData()->GetAbstractArray("RegionId");
    if(labels==nullptr){
        stringstream msg;
        msg<<"[ttkOverlapTracking] ERROR: No label point data found"<<endl;
        dMsg(cout, msg.str(), timeMsg);
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

    this->overlapTracking.reset();

    return 1;
}

int ttkOverlapTracking::RequestInformation(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    cout<<"[ttkOverlapTracking] RequestInformation"<<endl;
    this->currentIndex = 0;
    return this->Superclass::RequestInformation(request, inputVector, outputVector);
}

int ttkOverlapTracking::RequestUpdateExtentInformation(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    cout<<"[ttkOverlapTracking] RequestUpdateExtentInformation"<<endl;
    return 1;
}

int ttkOverlapTracking::RequestUpdateExtent(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    cout<<"[ttkOverlapTracking] RequestUpdateExtent: "<<this->currentIndex<<endl;
    inInfo->Set( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->currentIndex);

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
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkOverlapTracking] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Get Output Graph
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkUnstructuredGrid* trackingGraph = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Process input based on mode
    auto nBlocks = inMB->GetNumberOfBlocks();
    if(nBlocks>1) { // Process all blocks without streaming
        this->overlapTracking.reset();
        for(size_t i=0; i<nBlocks; i++){
            int status = this->processTimestep( inMB->GetBlock(i) );
            if(status==0) return 0;
        }
        this->finalize( trackingGraph );
    } else if(nBlocks==1) {
        // Process streamed element
        int status = this->processTimestep( inMB->GetBlock(0) );
        if(status==0) return 0;

        this->currentIndex++;
        request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
        request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
    } else {
        // Streaming complete
        this->finalize( trackingGraph );
        request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    }

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
