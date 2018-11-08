#include <ttkTrackingFromOverlap.h>

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

vtkStandardNewMacro(ttkTrackingFromOverlap)

int ttkTrackingFromOverlap::processTimestep(vtkDataObject* dataObject){

    auto pointSet = vtkPointSet::SafeDownCast( dataObject );
    if(pointSet==nullptr){
        dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted to vtkPointSet\n", memoryMsg);
        return 0;
    }

    auto n = pointSet->GetNumberOfPoints();

    auto labels = pointSet->GetPointData()->GetAbstractArray( this->GetLabelScalarField().data() );
    if(labels==nullptr && n>0){
        dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels not found\n", timeMsg);
        return 0;
    }

    this->trackingFromOverlap.processTimestep(
        n>0 ? (float*) pointSet->GetPoints()->GetVoidPointer(0) : nullptr,
        n>0 ? (labelType*) labels->GetVoidPointer(0)            : nullptr,
        pointSet->GetNumberOfPoints()
    );

    return 1;
}

int ttkTrackingFromOverlap::finalize(vtkUnstructuredGrid* trackingGraph){
    auto& timeNodeLabelMap = this->trackingFromOverlap.getTimeNodeLabelMap();
    auto& timeEdgesMap = this->trackingFromOverlap.getTimeEdgesMap();
    size_t tn = timeNodeLabelMap.size();
    auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Point Data
    {
        size_t n = 0;
        for(size_t t=0; t<tn; t++)
            n += timeNodeLabelMap[t].size();

        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints( n );
        auto pointCoords = (float*) points->GetVoidPointer(0);

        vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
        time->SetNumberOfComponents(1);
        time->SetName("TimeIndex");
        time->SetNumberOfValues(n);
        auto timeData = (double*) time->GetVoidPointer(0);

        vtkSmartPointer<vtkIdTypeArray> label = vtkSmartPointer<vtkIdTypeArray>::New();
        label->SetNumberOfComponents(1);
        label->SetName( this->GetLabelScalarField().data() );
        label->SetNumberOfValues(n);
        auto labelData = (long long*) label->GetVoidPointer(0);

        size_t q=0;
        size_t q2=0;
        for(size_t t=0; t<tn; t++){
            auto& nodes = timeNodeLabelMap[t];
            for(size_t i=0; i<nodes.size(); i++){
                pointCoords[q++] = t;
                pointCoords[q++] = i;
                pointCoords[q++] = 0;

                timeData[q2] = t;
                labelData[q2] = nodes[i];
                q2++;
            }
        }

        mesh->SetPoints(points);

        auto pointData = mesh->GetPointData();
        pointData->AddArray( time );
        pointData->AddArray( label );
    }

    // Cell Data
    {
        vtkSmartPointer<vtkUnsignedLongLongArray> overlap = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        overlap->SetNumberOfComponents(1);
        overlap->SetName("Overlap");

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
                overlap->InsertNextValue( e.overlap );
            }
        }
        auto cellData = mesh->GetCellData();
        cellData->AddArray( overlap );
    }

    trackingGraph->ShallowCopy( mesh );

    return 1;
}

int ttkTrackingFromOverlap::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkTrackingFromOverlap] RequestData" << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    double n = inMB->GetInformation()->Get( vtkDataObject::DATA_TIME_STEP() );

    // First Timestep
    if(i==0){
        this->trackingFromOverlap.reset();
    }

    size_t nBlocks = inMB->GetNumberOfBlocks();

    // Process Timestep(s)
    for(size_t i=0; i<nBlocks; i++){
        int status = this->processTimestep( inMB->GetBlock(i) );
        if(status==0) return 0;
        if(nBlocks>1)
            dMsg(cout, "[ttkTrackingFromOverlap] -----------------------------------------------------------\n", timeMsg);
        this->updateProgress( ((float)i)/((float)(nBlocks-1)) );
    }

    // Last Timestep
    if(i==n-1 || nBlocks>1){
        if(nBlocks<2)
            dMsg(cout, "[ttkTrackingFromOverlap] -----------------------------------------------------------\n", timeMsg);

        // Get Output
        vtkInformation* outInfo = outputVector->GetInformationObject(0);
        auto trackingGraph = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

        // Create Tracking Graph
        this->finalize( trackingGraph );

        // Print Status
        dMsg(cout, "[ttkTrackingFromOverlap] Tracking Graph Finalized\n", timeMsg);
    }

    return 1;
}