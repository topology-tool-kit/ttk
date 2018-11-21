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

    // Get Input
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = inInfo->Get( vtkDataObject::DATA_OBJECT() );

    // Get Output
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outMB = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Get iteration information
    double i = inInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );

    // First timestep
    if(i==0 || this->AggregatedMultiBlockDataSet==nullptr)
        this->AggregatedMultiBlockDataSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    size_t n = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();

    auto inputType = input->GetDataObjectType();

    vtkDataObject* copy = nullptr;
    switch (inputType) {
        case VTK_DIRECTED_GRAPH:
            copy = vtkDirectedGraph::New();
            break;
        case VTK_MOLECULE:
        case VTK_UNDIRECTED_GRAPH:
            copy = vtkUndirectedGraph::New();
            break;
        case VTK_IMAGE_DATA:
            copy = vtkImageData::New();
            break;
        case VTK_POLY_DATA:
            copy = vtkPolyData::New();
            break;
        case VTK_RECTILINEAR_GRID:
            copy = vtkRectilinearGrid::New();
            break;
        case VTK_STRUCTURED_GRID:
            copy = vtkStructuredGrid::New();
            break;
        case VTK_STRUCTURED_POINTS:
            copy = vtkStructuredPoints::New();
            break;
        case VTK_TABLE:
            copy = vtkTable::New();
            break;
        case VTK_TREE:
            copy = vtkTree::New();
            break;
        case VTK_UNSTRUCTURED_GRID:
            copy = vtkUnstructuredGrid::New();
            break;
        case VTK_MULTIBLOCK_DATA_SET:
            copy = vtkMultiBlockDataSet::New();
            break;
        case VTK_MULTIPIECE_DATA_SET:
            copy = vtkMultiPieceDataSet::New();
            break;
        case VTK_HIERARCHICAL_BOX_DATA_SET:
            copy = vtkHierarchicalBoxDataSet::New();
            break;
        case VTK_OVERLAPPING_AMR:
            copy = vtkOverlappingAMR::New();
            break;
        case VTK_NON_OVERLAPPING_AMR:
            copy = vtkNonOverlappingAMR::New();
            break;
    }
    if(copy==nullptr){
        dMsg(cout, "[ttkBlockAggregator] ERROR: Unsupported DataType\n", fatalMsg);
        return 0;
    }
    copy->DeepCopy( input );

    this->AggregatedMultiBlockDataSet->SetBlock( n, copy );
    outMB->ShallowCopy( this->AggregatedMultiBlockDataSet );

    // Print status
    {
        stringstream msg;
        msg << "[ttkBlockAggregator] Added block at index " << n << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    return 1;
}