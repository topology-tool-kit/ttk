#include <ttkMeshGraph.h>

#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkDataSetTriangleFilter.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkMeshGraph)

int ttkMeshGraph::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    Timer t;
    Memory m;

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

    const float* inputPoints = (float*) input->GetPoints()->GetVoidPointer(0);

    vtkCellArray* inputCells = input->GetCells();

    size_t nInputPoints = input->GetNumberOfPoints();
    size_t nInputCells = input->GetNumberOfCells();

    // Output Points
    auto nOutputPoints = meshGraph.computeNumberOfOutputPoints(nInputPoints, nInputCells, this->Subdivisions);
    auto outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetNumberOfPoints( nOutputPoints );
    auto outputVertices = (float*) outputPoints->GetVoidPointer(0);

    // Output Topology
    auto outputTopologySize = meshGraph.computeOutputTopologySize(nInputCells, this->Subdivisions);
    auto outputCells = vtkSmartPointer<vtkIdTypeArray>::New();
    outputCells->SetNumberOfValues( outputTopologySize );

    auto inputPointSizes = input->GetPointData()->GetArray( this->GetSizeFieldName().data() );
    if(!inputPointSizes){
        dMsg(cout, "[ttkMeshGraph] ERROR: Point data '" + this->GetSizeFieldName() + "' not found.\n", fatalMsg);
        return 0;
    }

    int status = 0;

    switch(inputPointSizes->GetDataType()){
        ttkTemplateMacro({
            status = meshGraph.execute<vtkIdType TTK_COMMA VTK_TT>(
                // Input
                inputPoints,
                inputCells->GetPointer(),
                (VTK_TT*) inputPointSizes->GetVoidPointer(0),
                nInputPoints,
                nInputCells,
                this->GetSubdivisions(),
                this->SizeScale,

                this->PrimaryAxis,
                this->SecondaryAxis,

                // Output
                outputVertices,
                (vtkIdType*) outputCells->GetVoidPointer(0)
            );
        });
    }
    if(status!=1) return 0;


    auto meshedGraph = vtkSmartPointer<vtkUnstructuredGrid>::New();
    meshedGraph->SetPoints( outputPoints );
    auto outputCellArray = vtkSmartPointer<vtkCellArray>::New();
    outputCellArray->SetCells(nInputCells, outputCells);
    meshedGraph->SetCells(VTK_POLYGON, outputCellArray);

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    if(this->GetSubdivisions()>1 && this->GetTetrahedralize()){
        auto dataSetTriangleFilter = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
        dataSetTriangleFilter->SetInputData( meshedGraph );
        dataSetTriangleFilter->Update();

        output->ShallowCopy( dataSetTriangleFilter->GetOutput() );
    } else {
        output->ShallowCopy( meshedGraph );
    }

    // Print status
    {
        stringstream msg;
        msg << "[ttkMeshGraph] -----------------------------------------------------------------" << endl
            << "[ttkMeshGraph]   Time: " << t.getElapsedTime() << " s" << endl
            << "[ttkMeshGraph] Memory: " << m.getElapsedUsage() << " MB" << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 1;
}
