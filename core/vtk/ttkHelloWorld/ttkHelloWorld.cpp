#include <ttkHelloWorld.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkPoints.h>
#include <vtkCellArray.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkHelloWorld);

/**
* TODO 3: Implement the filter constructor and destructor in the cpp file.
*
* The constructor has to specify the number of input and output ports
* with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
* respectively. It should also set default values for all filter
* parameters.
*
* The destructor is usually empty unless you want to manage memory
* explicitly, by for example allocating memory on the heap that needs
* to be freed when the filter is destroyed.
*/
ttkHelloWorld::ttkHelloWorld(){
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

ttkHelloWorld::~ttkHelloWorld(){}

/**
* TODO 4: Specify the required input data type of each input port
*
* This method specifies the required input object data types of the
* filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
* the port information.
*/
int ttkHelloWorld::FillInputPortInformation(int port, vtkInformation* info) {
    if(port==0)
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    else
        return 0;

    return 1;
}

/**
* TODO 5: Specify the data object type of each output port
*
* This method specifies in the port information the data type of the
* corresponding output objects. It is possible to either explicitly
* specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key, or to
* pass a type of an input port to an output port by adding the
* ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key.
*
* Note: prior to the execution of the RequestData method the pipeline will
* initialize empty output data objects based on this information.
*/
int ttkHelloWorld::FillOutputPortInformation(int port, vtkInformation* info) {
    if(port==0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    else
        return 0;

    return 1;
}

/**
* TODO 6: Implement the purpose of this filter
*
* This method is called during the pipeline execution to update the
* already initialized output data objects based on the given input
* data objects and filter parameters.
*
* Note: The output objects are already initialized based on the
* information provided by the FillOutputPortInformation method.
*/
int ttkHelloWorld::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Get input
    vtkDataSet* input = vtkDataSet::GetData( inputVector[0] );

    // Prepare output point buffer
    vtkSmartPointer<vtkPoints> boundingBoxPoints = vtkSmartPointer<vtkPoints>::New();
    boundingBoxPoints->SetNumberOfPoints( 8 ); // the 8 corner vertices of the bounding box

    // Get raw pointer to point coordinates
    float* boundingBoxPointCoordinates = (float*) boundingBoxPoints->GetVoidPointer(0);

    // Prepare output cell buffer
    // In VTK cells are defined by a connectivity list, which is a flat array
    // that contains for each cell 1 + n values, where the first value specifies
    // the number of points that span a cell, and the remaining values specify
    // the indicies of these points inside the vtkPoints array
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    // Initialize a buffer for one voxel cell:
    // 1 cell that is stored via 9 values: 1 cell dimension + 8 point indicies
    vtkIdType* boundingBoxConnectivityList = cells->WritePointer( 1, 1 + 8 );

    // Get triangulation of the input object (will create one if does not exist already)
    ttk::Triangulation* triangulation = ttkAlgorithm::GetTriangulation( input );

    // Compute bounding box in the base code
    int status = this->computeBoundingBox<vtkIdType>(
        boundingBoxPointCoordinates,
        boundingBoxConnectivityList,
        triangulation,
        this->Scale
    );

    // On error cancel filter execution
    if(status==0)
        return 0;

    // Finalize output
    vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData( outputVector );
    output->SetPoints( boundingBoxPoints ); // assign the computed points to the output
    output->SetCells( VTK_VOXEL, cells ); // assign the computed cells to the output

    return 1;
}