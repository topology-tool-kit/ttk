#include <ttkCinemaWriter.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkFieldData.h>
#include <vtkStdString.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <stdlib.h>
#include <vtkDirectory.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaWriter)

int ttkCinemaWriter::RequestData (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    Timer t;
    Memory m;
    struct stat info;

    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaWriter] RequestData"<<endl;
        msg<<"[ttkCinemaWriter] Database Path: "<<this->DatabasePath<<endl;
        msg<<"[ttkCinemaWriter] Override     : "<<(this->OverrideDatabase?"yes":"no")<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Copy Input to Output
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outputMBD = vtkMultiBlockDataSet::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // If input is vtkMultiBlockDataSet then copy, otherwise add as single block
    if(input->IsA("vtkMultiBlockDataSet")){
        outputMBD->ShallowCopy( input );
    } else {
        outputMBD->SetBlock(0, input);
    }

    // Check if Database exists
    if(this->DatabasePath.compare("")==0 ){
        stringstream msg;
        msg<<"[ttkCinemaWriter] ERROR: Invalid database path"<<endl;
        dMsg(cout, msg.str(), timeMsg);
        return 0;
    }

    // Initialize path variables
    string dataPrefix = "data/";
    string pathPrefix = this->DatabasePath+"/" + dataPrefix;
    string pathSuffix = ".vtm";

    // Create directory if does not already exist
    {
        vtkSmartPointer<vtkDirectory> directory = vtkSmartPointer<vtkDirectory>::New();
        int opened = directory->Open( this->DatabasePath.data() );
        if(!opened)
            directory->MakeDirectory( this->DatabasePath.data() );
    }

    // If OverrideDatabase then delete old data products
    if( this->OverrideDatabase ){
        vtkSmartPointer<vtkDirectory> directory = vtkSmartPointer<vtkDirectory>::New();
        int opened = directory->Open( pathPrefix.data() );
        if(opened){
            int status = directory->DeleteDirectory( pathPrefix.data() );

            stringstream msg;
            if(status==0)
                msg<<"[ttkCinemaWriter] ERROR: Unable to delete existing data products"<<endl;
            else
                msg<<"[ttkCinemaWriter] Old data products deleted"<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    // Determine unique path to new products (for now just generate random number)
    string id;
    string path;
    bool unique = false;

    while(!unique){
        id = to_string( rand()%1000000 );
        path = pathPrefix + id + pathSuffix;
        unique = stat( path.data(), &info ) != 0;
    }

    // Write input to disk
    {
        stringstream msg;
        msg<<"[ttkCinemaWriter] Writing new data products to disk ... ";
        dMsg(cout, msg.str(), timeMsg);
    }
    vtkSmartPointer<vtkXMLMultiBlockDataWriter> mbWriter = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    mbWriter->SetFileName( path.data() );
    mbWriter->SetDataModeToAppended();
    mbWriter->SetCompressorTypeToZLib();
    vtkZLibDataCompressor::SafeDownCast( mbWriter->GetCompressor() )->SetCompressionLevel(9);
    mbWriter->SetInputData( outputMBD );
    mbWriter->Write();
    {
        stringstream msg;
        msg<<"Done"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    string dataCsvPath = this->DatabasePath+"/data.csv";

    // Create data.csv file if does not already exist
    if( stat( dataCsvPath.data(), &info ) != 0 ){
        string csv = "FILE\n";
        ofstream csvFile;
        csvFile.open( dataCsvPath.data() );
        csvFile << csv;
        csvFile.close();
    }

    // Update data.csv file
    {
        stringstream msg;
        msg<<"[ttkCinemaWriter] Updating data.csv file ... "<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }
    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetFileName( dataCsvPath.data() );
    reader->DetectNumericColumnsOff();
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(",");
    reader->Update();
    auto table = vtkTable::SafeDownCast( reader->GetOutput() );

    int n = outputMBD->GetNumberOfBlocks();
    int offset = this->OverrideDatabase ? 0 : table->GetNumberOfRows();
    if(this->OverrideDatabase)
        table->SetNumberOfRows(n);
    else
        for(int i=0; i<n; i++) // SetNumberOfRows clears existing rows -> insert blank rows
            table->InsertNextBlankRow();

    for(int i=0; i<n; i++){
        auto block = outputMBD->GetBlock(i);
        string blockExtension = "vtk";
        #if VTK_MAJOR_VERSION <= 7
            stringstream msg;
            msg << "Failed" << endl
                << "[ttkCinemaQuery] ERROR: VTK version too old."<<endl
                << "[ttkCinemaQuery]        This filter requires vtkXMLPMultiBlockDataWriter"<<endl
                << "[ttkCinemaQuery]        of version 7.0 or higher."<<endl;
            dMsg(cout, msg.str(), memoryMsg);
            return 0;
        #else
            blockExtension = this->GetDefaultFileExtensionForDataSet( block->GetDataObjectType() );
        #endif

        // auto fieldData = block->GetFieldData();
        for(int j=0; j<table->GetNumberOfColumns(); j++){
            auto columnName = table->GetColumnName(j);
            // auto columnFD = fieldData->GetAbstractArray(columnName);
            auto columnCSV = vtkStringArray::SafeDownCast( table->GetColumn(j) );
            if(string(columnName).compare("FILE")==0){
                columnCSV->SetValue(offset+i, vtkStdString( dataPrefix+id+"/"+id+"_"+to_string(i)+"."+blockExtension ));
            // TODO store field data in CSV
            // } else if(columnFD!=nullptr){
            //     columnCSV->SetValue(offset+i, columnFD->GetVariantValue(0).ToString());
            } else {
                columnCSV->SetValue(offset+i, "");
            }
        }
    }

    // Write data.csv file
    vtkSmartPointer<vtkDelimitedTextWriter> csvWriter = vtkSmartPointer<vtkDelimitedTextWriter>::New();
    csvWriter->SetUseStringDelimiter(false);
    csvWriter->SetFileName( (this->DatabasePath+"/data.csv").data() );
    csvWriter->SetInputData( table );
    csvWriter->Write();
    {
        stringstream msg;
        msg<<"Done"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaWriter]   time: " << t.getElapsedTime() << " s." << endl;
        msg << "[ttkCinemaWriter] memory: " << m.getElapsedUsage() << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
