#include <ttkCinemaProductReader.h>

#include <vtkTable.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkVariantArray.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaProductReader)

int ttkCinemaProductReader::RequestInformation(
    vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector
){
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    cout<<"[ttkCinemaProductReader] RequestInformation"<<endl;

    double* tSteps = new double[1];
    tSteps[0] = 0;
    double* tRange = new double[2];
    tRange[0] = 0;
    tRange[1] = 0;

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), tSteps, 1);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), tRange, 2);

    return this->Superclass::RequestInformation(request, inputVector, outputVector);
}

// int ttkCinemaProductReader::RequestUpdateExtentInformation(
//     vtkInformation* request,
//     vtkInformationVector** inputVector,
//     vtkInformationVector* outputVector
// ){
//     cout<<"[ttkCinemaProductReader] RequestUpdateExtentInformation"<<endl;
//     return this->Superclass::RequestInformation(request, inputVector, outputVector);
// };

int ttkCinemaProductReader::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaProductReader] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Prepare Input and Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // inInfo->Print(cout);
    // outInfo->Print(cout);

    // Read Data
    {
        // Determine number of files
        int n = inputTable->GetNumberOfRows();
        int m = inputTable->GetNumberOfColumns();
        cout<<"[ttkCinemaProductReader] Reading "<<n<<" files:"<<endl;

        // Compute DatabasePath
        auto databasePath = inputTable->GetFieldData()->GetAbstractArray("DatabasePath")->GetVariantValue(0).ToString();

        // Get column that contains paths
        auto paths = inputTable->GetColumnByName( this->FilepathColumnName.data() );
        if(paths==nullptr){
            stringstream msg;
            msg<<"[ttkCinemaProductReader] ERROR: Table does not have FilepathColumn '"<<this->FilepathColumnName<<"'"<<endl;
            dMsg(cerr, msg.str(), timeMsg);
            return 0;
        }

        int i, limit;
        if(this->UseStreaming){
            i = (int) inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
            limit = i+1;

            if(i>=n) // Send empty vtkMultiBlock to indicate end of sequence
                return 1;
        } else {
            i = 0;
            limit = n;
        }

        // For each row
        for(; i<limit; i++){
            // Get path
            auto path = databasePath + "/" + paths->GetVariantValue(i).ToString();
            auto ext = path.substr( path.length() - 3 );

            {
                stringstream msg;
                msg<<"[ttkCinemaProductReader]    "<<i<<": "<<path<<endl;
                dMsg(cout, msg.str(), timeMsg);
            }

            // Check if file exists
            std::ifstream infile(path.data());
            bool exists = infile.good();
            if(!exists){
                stringstream msg;
                msg<<"[ttkCinemaProductReader]    ERROR: File does not exist."<<endl;
                dMsg(cerr, msg.str(), timeMsg);
                continue;
            }

            // Read any data using vtkXMLGenericDataObjectReader
            {
                vtkSmartPointer<vtkXMLGenericDataObjectReader> reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
                reader->SetFileName( path.data() );
                reader->Update();
                output->SetBlock( this->UseStreaming ? 0 : i, reader->GetOutput());
            }

            // Augment read data with row information
            // TODO: Make Optional
            auto block = output->GetBlock( this->UseStreaming ? 0 : i );
            for(int j=0; j<m; j++){
                auto columnName = inputTable->GetColumnName(j);
                auto fieldData = block->GetFieldData();
                if(!fieldData->HasArray( columnName )){
                    bool isNumeric = inputTable->GetColumn(j)->IsNumeric();

                    if(isNumeric){
                        vtkSmartPointer<vtkDoubleArray> c = vtkSmartPointer<vtkDoubleArray>::New();
                        c->SetName( columnName );
                        c->SetNumberOfValues(1);
                        c->SetValue(0, inputTable->GetValue(i,j).ToDouble());
                        fieldData->AddArray( c );
                    } else {
                        vtkSmartPointer<vtkStringArray> c = vtkSmartPointer<vtkStringArray>::New();
                        c->SetName( columnName );
                        c->SetNumberOfValues(1);
                        c->SetValue(0, inputTable->GetValue(i,j).ToString());
                        fieldData->AddArray( c );
                    }
                }
            }

            if(!this->UseStreaming)
                this->updateProgress( ((float)i)/((float)(n-1)) );
        }

    }

    // Print status
    {
        stringstream msg;
        msg << "[ttkCinemaProductReader] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}