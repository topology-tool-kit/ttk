#include <ttkMeshGraph.h>

#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMeshGraph)

int ttkMeshGraph::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkMeshGraph] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Set Wrapper
    meshGraph.setWrapper(this);

    // Prepare input and output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = vtkUnstructuredGrid::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    float* inputPoints = (float*) input->GetPoints()->GetVoidPointer(0);

    auto inputCells = input->GetCells();
    vtkIdType* inputTopology = inputCells->GetPointer();

    size_t nInputPoints = input->GetNumberOfPoints();
    size_t nInputCells = input->GetNumberOfCells();

    this->Subdivisions = 3;

    // Output Points
    auto nOutputPoints = meshGraph.computeNumberOfOutputPoints(nInputPoints, nInputCells, this->Subdivisions);
    auto outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetNumberOfPoints( nOutputPoints );
    auto outputVertices = (float*) outputPoints->GetVoidPointer(0);

    // Output Topology
    auto outputTopologySize = meshGraph.computeOutputTopologySize(nInputCells, this->Subdivisions);
    auto outputCells = vtkSmartPointer<vtkIdTypeArray>::New();
    outputCells->SetNumberOfValues( outputTopologySize );
    vtkIdType* outputTopology = (vtkIdType*) outputCells->GetVoidPointer(0);

    auto inputPointSizeArray = input->GetPointData()->GetArray( "Size2" );
    double* inputPointSizes = inputPointSizeArray ? (double*) inputPointSizeArray->GetVoidPointer(0) : nullptr;

    meshGraph.execute<vtkIdType>(
        // Input
        inputPoints,
        inputTopology,
        inputPointSizes,
        nInputPoints,
        nInputCells,
        this->Subdivisions,

        // Output
        outputVertices,
        outputTopology
    );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    output->SetPoints( outputPoints );

    auto outputCellArray = vtkSmartPointer<vtkCellArray>::New();
    outputCellArray->SetCells(nInputCells, outputCells);
    output->SetCells(VTK_POLYGON, outputCellArray);

    return 1;
}
