#include <ttkPlanarGraphLayout.h>

#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPlanarGraphLayout)

int ttkPlanarGraphLayout::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkPlanarGraphLayout] RequestData"<<endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Set Wrapper
    planarGraphLayout_.setWrapper(this);
    planarGraphLayout_.setDebugLevel(this->debugLevel_);

    // Prepare input and output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = vtkUnstructuredGrid::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Copy input to output
    output->ShallowCopy(input);

    size_t nPoints = output->GetNumberOfPoints();
    size_t nEdges = output->GetNumberOfCells();

    auto outputPointData = output->GetPointData();

    vtkSmartPointer<vtkFloatArray> layout = vtkSmartPointer<vtkFloatArray>::New();
    layout->SetName("Layout");
    layout->SetNumberOfComponents(1);
    layout->SetNumberOfValues(nPoints);
    auto layoutData = (float*) layout->GetVoidPointer(0);

    outputPointData->AddArray( layout );

    auto cells = output->GetCells();
    long long* topology = cells->GetPointer();
    auto levels = outputPointData->GetAbstractArray( this->GetAxisFieldName().data() );
    if(!levels){
        stringstream msg;
        msg<<"[ttkPlanarGraphLayout] ERROR: Input point data does not have array '" << this->GetAxisFieldName() << "'" <<endl;
        dMsg(cout, msg.str(), fatalMsg);
        return 0;
    }

    switch(levels->GetDataType()){
        vtkTemplateMacro({
            planarGraphLayout_.execute<VTK_TT>(
                // Input
                levels->GetVoidPointer(0),
                topology,
                nPoints,
                nEdges,

                // Output
                layoutData
            );
        });
    }

    return 1;
}
