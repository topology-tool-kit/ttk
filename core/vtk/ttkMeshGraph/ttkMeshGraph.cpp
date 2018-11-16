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
    meshGraph_.setWrapper(this);

    // Prepare input and output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = vtkUnstructuredGrid::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    float* inputVertices = (float*) input->GetPoints()->GetVoidPointer(0);

    auto inputCells = input->GetCells();
    vtkIdType* inputTopology = inputCells->GetPointer();

    size_t nVertices = input->GetNumberOfPoints();
    size_t nEdges = input->GetNumberOfCells();

    auto outputVertexNumber = nVertices*2 + nEdges*(this->Subdivisions*2); // one vertex becomes two + subdivisons on edges
    auto outputCellVertexNumber = 5+this->Subdivisions*2; // cellDim + 4 corners + 2 per each subdivision

    auto outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetNumberOfPoints( nVertices*outputVertexNumber );
    auto outputVertices = (float*) outputPoints->GetVoidPointer(0);

    auto outputCells = vtkSmartPointer<vtkIdTypeArray>::New();
    outputCells->SetNumberOfValues( nEdges*outputCellVertexNumber );
    vtkIdType* outputTopology = (vtkIdType*) outputCells->GetVoidPointer(0);

    meshGraph_.execute<vtkIdType>(
        // Input
        inputVertices,
        inputTopology,
        nullptr,
        nVertices,
        nEdges,
        this->Subdivisions,

        // Output
        outputVertices,
        outputTopology
    );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    output->SetPoints( outputPoints );

    auto outputCellArray = vtkSmartPointer<vtkCellArray>::New();
    outputCellArray->SetCells(nEdges, outputCells);
    output->SetCells(VTK_TRIANGLE, outputCellArray);


    return 1;
}
