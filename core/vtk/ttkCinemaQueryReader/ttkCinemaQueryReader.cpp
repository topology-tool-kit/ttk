#include <ttkCinemaQueryReader.h>
#include <vtkVariantArray.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQueryReader)

int ttkCinemaQueryReader::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaQueryReader] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Prepare Input and Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* inputTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Read Data
    {
        // Determine number of files
        int n = inputTable->GetNumberOfRows();
        int m = inputTable->GetNumberOfColumns();
        cout<<"[ttkCinemaQueryReader] Reading "<<n<<" files:"<<endl;

        // Compute DatabasePath
        auto databasePath = inputTable->GetFieldData()->GetAbstractArray("DatabasePath")->GetVariantValue(0).ToString();
        databasePath = databasePath.substr( 0, databasePath.find_last_of("/"));

        // Get column that contains paths
        auto paths = inputTable->GetColumnByName( this->FilepathColumnName.data() );
        if(paths==nullptr){
            stringstream msg;
            msg<<"[ttkCinemaQueryReader] ERROR: Table does not have FilepathColumn '"<<this->FilepathColumnName<<"'"<<endl;
            dMsg(cerr, msg.str(), timeMsg);
            return 0;
        }

        // For each row
        for(int i=0; i<n; i++){
            // get path
            auto path = databasePath + "/" + paths->GetVariantValue(i).ToString();
            auto ext = path.substr( path.length() - 3 );

            std::ifstream infile(path.data());
            bool exists = infile.good();

            {
                stringstream msg;
                msg<<"[ttkCinemaQueryReader]    "<<i<<": "<<path<<endl;
                dMsg(cout, msg.str(), timeMsg);
            }

            if(!exists){
                stringstream msg;
                msg<<"[ttkCinemaQueryReader]    ERROR: File does not exist."<<endl;
                dMsg(cerr, msg.str(), timeMsg);
                continue;
            }

            // load data using correct reader
            if(ext=="vti"){
                vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
                reader->SetFileName( path.data() );
                reader->Update();
                output->SetBlock(i, reader->GetOutput());
            } else if(ext=="vtu"){
                vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
                reader->SetFileName( path.data() );
                reader->Update();
                output->SetBlock(i, reader->GetOutput());
            } else if(ext=="vtp"){
                vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
                reader->SetFileName( path.data() );
                reader->Update();
                output->SetBlock(i, reader->GetOutput());
            } else {
                cout<<"[ttkCinemaQueryReader] Unknown File type: "<<ext<<endl;
            }

            // Augment read data with row information
            // TODO: Make Optional
            auto block = output->GetBlock(i);
            for(int j=0; j<m; j++){
                vtkSmartPointer<vtkVariantArray> c = vtkSmartPointer<vtkVariantArray>::New();
                c->SetName( inputTable->GetColumnName(j) );
                c->SetNumberOfValues(1);
                c->SetValue(0, inputTable->GetValue(i,j));
                block->GetFieldData()->AddArray( c );
            }
        }
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaQueryReader] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}