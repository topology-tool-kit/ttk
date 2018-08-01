#include <ttkCinemaQueryReader.h>
#include <vtkVariantArray.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiBlockDataGroupFilter.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQueryReader)

int ttkCinemaQueryReader::RequestData(vtkInformation* request,
    vtkInformationVector** inputVector, vtkInformationVector* outputVector){

    Memory m;

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    auto* paths = vtkVariantArray::SafeDownCast( inputTable->GetColumnByName("CDB_Filepath") );
    size_t max = paths->GetSize();
    for(size_t i=0; i<max; i++){
        auto path = paths->GetValue(i).ToString();
        auto ext = path.substr( path.length() - 3 );

        if(ext=="vti"){
            vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
            reader->SetFileName(path.data());
            reader->Update();
            output->SetBlock(i, reader->GetOutput());
        } else if(ext=="vtu"){
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(path.data());
            reader->Update();
            output->SetBlock(i, reader->GetOutput());
        } else if(ext=="vtp"){
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(path.data());
            reader->Update();
            output->SetBlock(i, reader->GetOutput());
        } else {
            cout<<"Unknown File type: "<<ext<<endl;
        }

        // vtkNew<vtkCubeSource> cube1;
    }

  // Append geometry of the two first meshes
//   vtkNew<vtkAppendFilter> append;
//   append->AddInputConnection(id1->GetOutputPort());
//   append->Update();
//   vtkUnstructuredGrid* aug = append->GetOutput();

//     vtkNew<vtkMultiBlockDataGroupFilter> group;
//     group->AddInputData( reader->GetOutput() );

//     output->ShallowCopy( group->GetOutput() );

    // doIt(
    //     input,
    //     data
    // );
    // outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), data->GetExtent(), 6);
    // outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), data->GetExtent());

    return 1;
}

// int ttkCinemaQueryReader::RequestInformation (
//   vtkInformation* request,
//   vtkInformationVector** inputVector,
//   vtkInformationVector *outputVector)
// {
//     vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
//     vtkTable* input = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
//     auto* paths = vtkVariantArray::SafeDownCast( input->GetColumnByName("path") );
//     auto path = paths->GetValue( 0 ).ToString();

//     vtkInformation* outInfo = outputVector->GetInformationObject(0);

//     cout<<"Z"<<endl;



//     // vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
//     // reader->SetFileName(path.data());
//     // reader->UpdateInformation();

//     // outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
//     //     reader->GetOutputInformation(0)->Get(
//     //         vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()
//     //     ),
//     //     6
//     // );
//   return 1;
// }