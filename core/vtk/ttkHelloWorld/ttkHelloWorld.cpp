#include <ttkHelloWorld.h>

#include <vtkDataObject.h> // For output port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

vtkStandardNewMacro(ttkHelloWorld);

ttkHelloWorld::ttkHelloWorld(){
    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
}

ttkHelloWorld::~ttkHelloWorld(){
}

int ttkHelloWorld::FillInputPortInformation(int port, vtkInformation* info) {
    if(port!=0) return 0;
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

int ttkHelloWorld::FillOutputPortInformation(int port, vtkInformation* info) {
    if(port!=0) return 0;
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
}

int ttkHelloWorld::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    return 1;
}