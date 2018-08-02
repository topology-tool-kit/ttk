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

int ttkCinemaQueryReader::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto t = outInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() );
    cout<<"-------------------------------------------------------------"<<endl;
    cout<<"[ttkCinemaQueryReader] RequestData t="<<t<<endl;

    // Memory m;
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    int n = inputTable->GetNumberOfRows();
    cout<<"[ttkCinemaQueryReader] Req n: "<<n<<endl;



    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    // for(size_t i=0; i<n; i++){
    {
        size_t i = t;

        auto* paths = vtkVariantArray::SafeDownCast( inputTable->GetColumnByName("CDB_Filepath") );
        auto path = paths->GetValue(i).ToString();
        cout<<"[ttkCinemaQueryReader] Loading: "<<path<<endl;
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
            cout<<"[ttkCinemaQueryReader] Unknown File type: "<<ext<<endl;
        }
    }
    cout<<"-------------------------------------------------------------"<<endl;

    return 1;
}

int ttkCinemaQueryReader::RequestInformation (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    vtkMultiBlockDataSetAlgorithm::RequestInformation(request, inputVector, outputVector);

    cout<<"[ttkCinemaQueryReader] RequestInformation"<<endl;
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto n = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    cout<<"[ttkCinemaQueryReader] RequestInformation n "<< n <<endl;

    // vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    // vtkTable* input = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    // int n = input->GetNumberOfRows();
    // cout<<"[ttkCinemaQueryReader] n: "<<n<<endl;


    // vtkInformation* outInfo = outputVector->GetInformationObject(0);

    // vector<double> steps(n);
    // for(int i=0; i<n; i++)
    //     steps[i] = i;
    // outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_STEPS(), steps.data(), n );

    // vector<double> range(2);
    // range[0]=0;
    // range[1]=n;
    // outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_RANGE(), range.data(), 2 );

    return 1;
}

int ttkCinemaQueryReader::RequestModified(vtkInformation* info, int when){
    cout<<"[ttkCinemaQueryReader] RequestModified "<<when<<endl;
    return 0;
}

// int ttkCinemaQueryReader::RequestUpdateExtend (
//     vtkInformation* request,
//     vtkInformationVector** inputVector,
//     vtkInformationVector* outputVector
// ){
//     // vtkMultiBlockDataSetAlgorithm::RequestInformation(request, inputVector, outputVector);

//     // cout<<"ttkCinemaQueryReader: RequestInformation"<<endl;

//     // vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
//     // vtkTable* input = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
//     // int n = input->GetNumberOfRows();

//     // vtkInformation* outInfo = outputVector->GetInformationObject(0);

//     // vector<double> steps(n);
//     // for(int i=0; i<n; i++)
//     //     steps[i] = i;
//     // outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_STEPS(), steps.data(), n );

//     // vector<double> range(2);
//     // range[0]=0;
//     // range[1]=n;
//     // outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_RANGE(), range.data(), 2 );

//     // return 1;
// }