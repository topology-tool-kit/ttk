#include                <ttkCinemaQuery.h>
#include                <vtkStringArray.h>
#include                <vtkVariantArray.h>
#include                <vtkSmartPointer.h>
#include                <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQuery)

int ttkCinemaQuery::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

    cout<<"-------------------------------------------------------------"<<endl;
    cout<<"[ttkCinemaQuery] RequestData"<<endl;
    Memory m;

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkTable* outputTable = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    for(int i=outputTable->GetNumberOfColumns()-1; i>=0; i--)
        outputTable->RemoveColumn(i);

    string serverAddress = "127.0.0.1:8888";

    vector<vector<string>> rMatrix = cinemaQuery_.execute<int>( serverAddress, QueryString );

    unsigned int rows = rMatrix.size();
    unsigned int cols = rMatrix[0].size();

    if(rows<1)
        return 0;

    for ( unsigned int i = 0; i < cols; i++ ) {
        vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();

        column->SetNumberOfTuples(rows-1);
        column->SetName( rMatrix[0][i].data() );

        for (unsigned int j = 1; j < rows; j++) {
            column->SetValue(j-1, vtkVariant(rMatrix[j][i]) );
        }
        outputTable->AddColumn(column);
    }

    cout<<"-------------------------------------------------------------"<<endl;

    return 1;
}

int ttkCinemaQuery::RequestInformation (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    vtkTableReader::RequestInformation(request, inputVector, outputVector);
    cout<<"[ttkCinemaQuery] RequestInformation"<<endl;

    string serverAddress = "127.0.0.1:8888";
    vector<vector<string>> rMatrix = cinemaQuery_.execute<int>( serverAddress, QueryString );
    unsigned int n = rMatrix.size()-1;

    vector<double> steps(n);
    for(int i=0; i<n; i++)
        steps[i] = i;
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_STEPS(), steps.data(), n );

    cout<<"[ttkCinemaQuery] RequestInformation n "<< n <<endl;

    // this->RequestData(request, inputVector, outputVector);

    // vtkTable* outputTable = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    // int n = outputTable->GetNumberOfRows();

    // cout<<"[ttkCinemaQuery] Set tn "<<n<<endl;


    // vector<double> range(2);
    // range[0]=0;
    // range[1]=n;
    // outInfo->Set( vtkStreamingDemandDrivenPipeline::TIME_RANGE(), range.data(), 2 );

    return 1;
}