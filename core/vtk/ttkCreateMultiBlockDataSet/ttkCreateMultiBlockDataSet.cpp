#include <ttkCreateMultiBlockDataSet.h>

#include <vtkMultiBlockDataSet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCreateMultiBlockDataSet)

int ttkCreateMultiBlockDataSet::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkCreateMultiBlockDataSet] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Timer t;

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outputMB = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // Set Blocks
    for(size_t i=0; i<5; i++){
        vtkInformationVector* inInfoVector = inputVector[i];
        if(inInfoVector!=nullptr){
            auto inInfo = inInfoVector->GetInformationObject(0);
            if(inInfo!=nullptr)
                outputMB->SetBlock( i, inInfo->Get(vtkDataObject::DATA_OBJECT()) );
        }
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCreateMultiBlockDataSet] MultiBlockDataSet created in " << t.getElapsedTime() << " s." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}