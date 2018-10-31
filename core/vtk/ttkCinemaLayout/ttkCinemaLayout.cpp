#include <ttkCinemaLayout.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaLayout)

int ttkCinemaLayout::RequestData (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    Timer t;
    Memory m;

    // Print Status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkCinemaLayout] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Copy Input to Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto inputMB = vtkMultiBlockDataSet::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outputMB = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    size_t nBlocks = inputMB->GetNumberOfBlocks();
    size_t n = ceil( sqrt(nBlocks) );

    auto firstImage = vtkImageData::SafeDownCast( inputMB->GetBlock(0) );
    double* bounds = firstImage->GetBounds();
    double width = bounds[1]-bounds[0];
    double height = bounds[3]-bounds[2];

    size_t i = 0;
    double y = 0;
    while(i<nBlocks){
        for(size_t x=0; x<n && i<nBlocks; x++){
            auto image = vtkImageData::SafeDownCast( inputMB->GetBlock(i) );

            auto transform = vtkSmartPointer<vtkTransform>::New();
            transform->Translate(x*width,y,0);

            auto transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
            transformFilter->SetTransform( transform );
            transformFilter->SetInputData( image );
            transformFilter->Update();

            outputMB->SetBlock(i, transformFilter->GetOutput());

            i++;
        }
        y += height;
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaLayout] --------------------------------------------------------------"<<endl;
        msg << "[ttkCinemaLayout]   time: " << t.getElapsedTime() << " s." << endl;
        msg << "[ttkCinemaLayout] memory: " << m.getElapsedUsage() << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
