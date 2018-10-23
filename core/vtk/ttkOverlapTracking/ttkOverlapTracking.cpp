#include <ttkOverlapTracking.h>

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

int ttkOverlapTracking::processTimestep(vtkPointSet* block){
    this->overlapTracking.processTimestep(
        (float*) block->GetPoints()->GetVoidPointer(0),
        (labelType*) block->GetPointData()->GetAbstractArray("RegionId")->GetVoidPointer(0),
        block->GetNumberOfPoints()
    );
    return 1;
}

int ttkOverlapTracking::finalize(vtkUnstructuredGrid* trackingGraph){
    auto& timeNodesMap = this->overlapTracking.getTimeNodesMap();
    auto& timeEdgesMap = this->overlapTracking.getTimeEdgesMap();
    size_t tn = timeNodesMap.size();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Point Data
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(size_t t=0; t<tn; t++){
            auto& nodes = timeNodesMap[t];
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
            offset[t] = offset[t-1] + timeNodesMap[t-1].size();

        size_t nEdges = 0;
        for(size_t t=0; t<tn-1; t++){
            nEdges += timeEdgesMap[t].size();
        }

        for(size_t t=0; t<tn-1; t++){
            auto& edges = timeEdgesMap[t];
            for(auto& e: edges){
                vtkIdType ids[2] = { (vtkIdType)(e.i+offset[t]), (vtkIdType)(e.j+offset[t+1]) };
                mesh->InsertNextCell(VTK_LINE, 2, ids);
                r0->InsertNextValue( timeNodesMap[t  ][e.i].label );
                r1->InsertNextValue( timeNodesMap[t+1][e.j].label );
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

    auto nBlocks = inMB->GetNumberOfBlocks();

    if(nBlocks>1) { // Not streaming
        this->overlapTracking.reset();

        for(size_t i=0; i<nBlocks; i++)
            this->processTimestep( vtkPointSet::SafeDownCast(inMB->GetBlock(i)) );

        this->finalize( trackingGraph );

    } else if(nBlocks==1) { // Streaming
        cout<<"Streaming not supported yet"<<endl;
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
