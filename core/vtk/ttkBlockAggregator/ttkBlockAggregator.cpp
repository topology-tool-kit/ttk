#include <ttkBlockAggregator.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCompositeDataReader.h>
#include <vtkDirectedGraph.h>
#include <vtkGraph.h>
#include <vtkGraphReader.h>
#include <vtkHierarchicalBoxDataSet.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMolecule.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkNonOverlappingAMR.h>
#include <vtkObjectFactory.h>
#include <vtkOverlappingAMR.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPoints.h>
#include <vtkTable.h>
#include <vtkTree.h>
#include <vtkUndirectedGraph.h>
#include <vtkUnstructuredGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkBlockAggregator)

int ttkBlockAggregator::AddBlock(vtkDataObject* dataObject, vtkMultiBlockDataSet* mb, size_t index, bool useShallowCopy){
    auto type = dataObject->GetDataObjectType();

    vtkSmartPointer<vtkDataObject> copy;
    switch (type) {
        case VTK_DIRECTED_GRAPH:
            copy = vtkSmartPointer<vtkDirectedGraph>::Take(vtkDirectedGraph::New());
            break;
        case VTK_MOLECULE:
        case VTK_UNDIRECTED_GRAPH:
            copy = vtkSmartPointer<vtkUndirectedGraph>::Take(vtkUndirectedGraph::New());
            break;
        case VTK_IMAGE_DATA:
            copy = vtkSmartPointer<vtkImageData>::Take(vtkImageData::New());
            break;
        case VTK_POLY_DATA:
            copy = vtkSmartPointer<vtkPolyData>::Take(vtkPolyData::New());
            break;
        case VTK_RECTILINEAR_GRID:
            copy = vtkSmartPointer<vtkRectilinearGrid>::Take(vtkRectilinearGrid::New());
            break;
        case VTK_STRUCTURED_GRID:
            copy = vtkSmartPointer<vtkStructuredGrid>::Take(vtkStructuredGrid::New());
            break;
        case VTK_STRUCTURED_POINTS:
            copy = vtkSmartPointer<vtkStructuredPoints>::Take(vtkStructuredPoints::New());
            break;
        case VTK_TABLE:
            copy = vtkSmartPointer<vtkTable>::Take(vtkTable::New());
            break;
        case VTK_TREE:
            copy = vtkSmartPointer<vtkTree>::Take(vtkTree::New());
            break;
        case VTK_UNSTRUCTURED_GRID:
            copy = vtkSmartPointer<vtkUnstructuredGrid>::Take(vtkUnstructuredGrid::New());
            break;
        case VTK_MULTIBLOCK_DATA_SET:
            copy = vtkSmartPointer<vtkMultiBlockDataSet>::Take(vtkMultiBlockDataSet::New());
            break;
        case VTK_MULTIPIECE_DATA_SET:
            copy = vtkSmartPointer<vtkMultiPieceDataSet>::Take(vtkMultiPieceDataSet::New());
            break;
        case VTK_HIERARCHICAL_BOX_DATA_SET:
            copy = vtkSmartPointer<vtkHierarchicalBoxDataSet>::Take(vtkHierarchicalBoxDataSet::New());
            break;
        case VTK_OVERLAPPING_AMR:
            copy = vtkSmartPointer<vtkOverlappingAMR>::Take(vtkOverlappingAMR::New());
            break;
        case VTK_NON_OVERLAPPING_AMR:
            copy = vtkSmartPointer<vtkNonOverlappingAMR>::Take(vtkNonOverlappingAMR::New());
            break;
    }
    if(copy==nullptr){
        dMsg(cout, "[ttkBlockAggregator] ERROR: Unsupported DataType\n", fatalMsg);
        return 0;
    }

    if(useShallowCopy)
        copy->ShallowCopy( dataObject );
    else
        copy->DeepCopy( dataObject );

    mb->SetBlock(index, copy);

    return 1;
}

int ttkBlockAggregator::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg << "================================================================================" << endl
            << "[ttkBlockAggregator] RequestData" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Get Input and Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    // Get iteration information
    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

    // First timestep
    if(!this->GetUseStreaming() || i==0 || this->AggregatedMultiBlockDataSet==nullptr)
        this->AggregatedMultiBlockDataSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();
    for(size_t i=0; i<5; i++){
        vtkInformationVector* inVector = inputVector[i];

        if(inVector->GetNumberOfInformationObjects()!=1) continue;

        vtkInformation* inInfo = inVector->GetInformationObject(0);
        auto input = inInfo->Get( vtkDataObject::DATA_OBJECT() );
        if(input){
            size_t index = nBlocks + i;

            // Add block
            this->AddBlock(input, this->AggregatedMultiBlockDataSet, index, true);

            // Print status
            stringstream msg;
            msg << "[ttkBlockAggregator] Added block at index " << index << endl;
            dMsg(cout, msg.str(), infoMsg);
        }
    }

    // Get Output
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outMB = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    outMB->ShallowCopy( this->AggregatedMultiBlockDataSet );

    return 1;
}