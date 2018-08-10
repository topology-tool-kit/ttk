#include                  <ttkA3Refinement.h>
#include                  <vtkImageData.h>
#include                  <vtkHyperTreeGrid.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkA3Refinement)

int ttkA3Refinement::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Welcome MSG
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkA3Refinement] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData* inGrid = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkHyperTreeGrid* outGrid = vtkHyperTreeGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Compute Refinement
    {
        a3Refinement_.execute<int>();
    }

    // Print Performance
    {
        stringstream msg;
        msg << "[ttkA3Refinement] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}

int ttkA3Refinement::ProcessTrees(vtkHyperTreeGrid*, vtkDataObject*){
    return 1;
}
