#include  <ttkCinemaQuery.h>
#include  <vtkSmartPointer.h>
#include  <vtkStringArray.h>
#include  <vtkFieldData.h>
#include  <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQuery)

int ttkCinemaQuery::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"-------------------------------------------------------------"<<endl;
        msg<<"[ttkCinemaQuery] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    string result;
    string sqlTableDefinition, sqlTableRows;

    // Get Input Table
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkTable* inTable = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Get Output Table
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkTable* outTable = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    // Convert Input Table to SQL Table
    {
        int nc = inTable->GetNumberOfColumns();
        int nr = inTable->GetNumberOfRows();

        sqlTableDefinition = "CREATE TABLE InputTable (";
        vector<bool> isNumeric(nc);
        for(int i=0; i<nc; i++){
            auto c = inTable->GetColumn(i);
            isNumeric[i] = c->IsNumeric();
            sqlTableDefinition += (i>0 ? "," : "") + string(c->GetName()) + " " + (isNumeric[i] ? "REAL" : "TEXT");
        }
        sqlTableDefinition+=")";

        sqlTableRows = "INSERT INTO InputTable VALUES ";
        for(int j=0; j<nr; j++){
            if(j>0) sqlTableRows+=",";
            sqlTableRows+="(";
            for(int i=0; i<nc; i++){
                if(i>0) sqlTableRows+=",";
                if(isNumeric[i])
                    sqlTableRows += inTable->GetValue(j,i).ToString();
                else
                    sqlTableRows += "'" + inTable->GetValue(j,i).ToString() + "'";
            }
            sqlTableRows+=")";
        }
    }

    // Compute Query Result
    {
        result = cinemaQuery.execute(
            sqlTableDefinition,
            sqlTableRows,
            this->QueryString
        );
    }

    // Process Result
    {
        if(result!=""){
            #if VTK_MAJOR_VERSION <= 7
                stringstream msg;
                msg << "[ttkCinemaQuery] ERROR: VTK version too old."<<endl
                    << "[ttkCinemaQuery]        This filter requires vtkDelimitedTextReader"<<endl
                    << "[ttkCinemaQuery]        of version 7.0 or higher."<<endl;
                dMsg(cout, msg.str(), memoryMsg);
                return 0;
            #else
                vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
                reader->SetReadFromInputString( true );
                reader->SetInputString( result );
                reader->DetectNumericColumnsOn();
                reader->SetHaveHeaders(true);
                reader->SetFieldDelimiterCharacters(",");
                reader->Update();

                outTable->ShallowCopy( reader->GetOutput() );
            #endif
        } else {
            vtkSmartPointer<vtkTable> emptyTable = vtkSmartPointer<vtkTable>::New();

            int nc = inTable->GetNumberOfColumns();
            for(int i=0; i<nc; i++){
                auto c = inTable->GetColumn(i);
                vtkSmartPointer<vtkStringArray> emptyColumn = vtkSmartPointer<vtkStringArray>::New();
                emptyColumn->SetNumberOfValues(1);
                emptyColumn->SetValue(0,"NULL");
                emptyColumn->SetName( c->GetName() );
                emptyTable->AddColumn( emptyColumn );
            }
            outTable->ShallowCopy( emptyTable );
        }

        outTable->GetFieldData()->ShallowCopy( inTable->GetFieldData() );
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkCinemaQuery] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
