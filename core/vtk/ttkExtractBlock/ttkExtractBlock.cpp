#include                  <ttkExtractBlock.h>

#include                  <vtkMultiBlockDataSet.h>
#include                  <vtkExtractBlock.h>
#include                  <vtkSmartPointer.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkExtractBlock)

int ttkExtractBlock::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkExtractBlock] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkMultiBlockDataSet* input = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Extract Block
    {
        vtkSmartPointer<vtkExtractBlock> extractor = vtkSmartPointer<vtkExtractBlock>::New();
        extractor->SetInputData( input );
        extractor->AddIndex( this->Index );
        extractor->Update();

        output->ShallowCopy( extractor->GetOutput() );
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkExtractBlock] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}