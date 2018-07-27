#include <ttkCinemaQueryLoader.h>
#include <vtkVariantArray.h>
#include <vtkXMLImageDataReader.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQueryLoader)

int ttkCinemaQueryLoader::doIt(vtkTable* input, vtkImageData* output){

    Memory m;

    auto* paths = vtkVariantArray::SafeDownCast( input->GetColumnByName("path") );
    auto path = paths->GetValue(0).ToString();

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(path.data());
    reader->Update();

    output->DeepCopy( reader->GetOutput() );

    {
        stringstream msg;
        msg << "[ttkCinemaQueryLoader] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

  return 0;
}

int ttkCinemaQueryLoader::RequestData(vtkInformation* request,
    vtkInformationVector** inputVector, vtkInformationVector* outputVector){

    Memory m;

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* input = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkImageData *data = vtkImageData::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    doIt(
        input,
        data
    );
    // outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), data->GetExtent(), 6);
    // outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), data->GetExtent());

    return 1;
}

int ttkCinemaQueryLoader::RequestInformation (
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector *outputVector)
{
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* input = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    auto* paths = vtkVariantArray::SafeDownCast( input->GetColumnByName("path") );

    // auto maxIndex = paths.GetSize();

    auto path = paths->GetValue( 0 ).ToString();

    cout<<"Z"<<endl;

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(path.data());
    reader->UpdateInformation();

    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    outInfo->Print(cout);
    outInfo->GetInformationObject(0)->Print();

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
        reader->GetOutputInformation(0)->Get(
            vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()
        ),
        6
    );
  return 1;
}