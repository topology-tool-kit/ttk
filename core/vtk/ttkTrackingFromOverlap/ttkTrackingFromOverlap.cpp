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
        dMsg(cerr, "[ttkTrackingFromOverlap] ERROR: Input data can not be converted to 'vtkPointSet'.\n", fatalMsg);
        return 0;
    }

    auto n = pointSet->GetNumberOfPoints();

    auto labels = pointSet->GetPointData()->GetAbstractArray( this->GetLabelFieldName().data() );
    if(labels==nullptr && n>0){
        dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' not found.\n" , fatalMsg);
        return 0;
    }
    // TODO: Fix identification of VTK_ID_TYPE_IMPL for vtk6 or lower.
    if( labels->GetDataType()!=VTK_LONG_LONG && (labels->GetDataType()!=VTK_ID_TYPE || labels->GetDataTypeSize()!=8) ){
        dMsg(cout, "[ttkTrackingFromOverlap] ERROR: Point labels '" + this->GetLabelFieldName() + "' are not of type 'Long Long'\n", fatalMsg);
        return 0;
    }

    this->trackingFromOverlap.processTimestep(
        n>0 ? (float*) pointSet->GetPoints()->GetVoidPointer(0) : nullptr,
        n>0 ? (labelType*) labels->GetVoidPointer(0)            : nullptr,
        pointSet->GetNumberOfPoints()
    );

    return 1;
}

int ttkTrackingFromOverlap::finalize(vtkDataObject* trackingGraphObject){
    auto trackingGraph = vtkUnstructuredGrid::SafeDownCast( trackingGraphObject );

    auto& timeNodesMap = this->trackingFromOverlap.getTimeNodesMap();
    auto& timeEdgesMap = this->trackingFromOverlap.getTimeEdgesMap();
    size_t tn = timeNodesMap.size();

    // Add Points
    {
        size_t n = 0;
        for(size_t t=0; t<tn; t++)
            n += timeNodesMap[t].size();

        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints( n );
        auto pointCoords = (float*) points->GetVoidPointer(0);

        vtkSmartPointer<vtkDoubleArray> time = vtkSmartPointer<vtkDoubleArray>::New();
        time->SetNumberOfComponents(1);
        time->SetName("TimeIndex");
        time->SetNumberOfValues(n);
        auto timeData = (double*) time->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedLongLongArray> size = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        size->SetNumberOfComponents(1);
        size->SetName("Size");
        size->SetNumberOfValues(n);
        auto sizeData = (unsigned long long*) size->GetVoidPointer(0);

        vtkSmartPointer<vtkIdTypeArray> label = vtkSmartPointer<vtkIdTypeArray>::New();
        label->SetNumberOfComponents(1);
        label->SetName( this->GetLabelFieldName().data() );
        label->SetNumberOfValues(n);
        auto labelData = (long long*) label->GetVoidPointer(0);

        size_t q=0;
        size_t q2=0;
        for(size_t t=0; t<tn; t++){
            auto& nodes = timeNodesMap[t];
            for(size_t i=0; i<nodes.size(); i++){
                Node& n = nodes[i];

                pointCoords[q++] = n.x;
                pointCoords[q++] = n.y;
                pointCoords[q++] = n.z;

                timeData[q2] = t;
                labelData[q2] = n.label;
                sizeData[q2] =  n.size;
                q2++;
            }
        }

        trackingGraph->SetPoints(points);

        auto pointData = trackingGraph->GetPointData();
        pointData->AddArray( time );
        pointData->AddArray( size );
        pointData->AddArray( label );
    }


    // Add Cells
    {
        size_t n = 0;
        size_t tnM1 = tn>0 ? tn-1 : 0;
        for(size_t t=0; t<tnM1; t++)
            n += timeEdgesMap[t].size();

        auto cells = vtkSmartPointer<vtkIdTypeArray>::New();
        cells->SetNumberOfValues( 3*n );
        auto cellIds = (vtkIdType*) cells->GetVoidPointer(0);

        vtkSmartPointer<vtkUnsignedLongLongArray> overlap = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        overlap->SetNumberOfValues( n );
        overlap->SetName("Overlap");

        auto overlapData = (unsigned long long*) overlap->GetVoidPointer(0);

        vector<size_t> offset(tnM1+1);
        offset[0] = 0;
        for(size_t t=1; t<tn; t++)
            offset[t] = offset[t-1] + timeNodesMap[t-1].size();

        size_t q0=0;
        size_t q1=0;
        for(size_t t=0; t<tnM1; t++){
            auto& edges = timeEdgesMap[t];
            for(auto& e: edges){
                cellIds[q0++] = 2;
                cellIds[q0++] = (vtkIdType) e.i+offset[t  ];
                cellIds[q0++] = (vtkIdType) e.j+offset[t+1];

                overlapData[q1++] = e.overlap;
            }
        }

        auto cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->SetCells(n, cells);
        trackingGraph->SetCells(VTK_LINE, cellArray);

        trackingGraph->GetCellData()->AddArray( overlap );
    }

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
        dMsg(cout, msg.str(), infoMsg);
    }

    // Set Wrapper
    this->trackingFromOverlap.setWrapper(this);

    // Get Input Object
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto inMB = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    double n = inMB->GetInformation()->Get( vtkDataObject::DATA_TIME_STEP() );

    // First Timestep
    if(i==0)
        this->trackingFromOverlap.reset();

    size_t nBlocks = inMB->GetNumberOfBlocks();

    // Process Timestep(s)
    for(size_t i=0; i<nBlocks; i++){
        int status = this->processTimestep( inMB->GetBlock(i) );
        if(status==0) return 0;
        if(nBlocks>1)
            dMsg(cout, "[ttkTrackingFromOverlap] -------------------------------------------------------\n", infoMsg);
        this->updateProgress( ((float)i)/((float)(nBlocks-1)) );
    }

    // Last Timestep
    if(i==n-1 || n==0){
        if(nBlocks<2)
            dMsg(cout, "[ttkTrackingFromOverlap] -------------------------------------------------------\n", infoMsg);

        // Get Output
        vtkInformation* outInfo = outputVector->GetInformationObject(0);
        auto trackingGraph = outInfo->Get(vtkDataObject::DATA_OBJECT());

        // Create Tracking Graph
        this->finalize( trackingGraph );

        // Print Status
        dMsg(cout, "[ttkTrackingFromOverlap] Tracking Graph Finalized\n", infoMsg);
    }

    return 1;
}