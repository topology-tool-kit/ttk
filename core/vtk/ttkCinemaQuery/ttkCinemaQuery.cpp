#include <ttkCinemaQuery.h>

#include <vtkTable.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>
#include <vtkDelimitedTextReader.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

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
        msg<<"================================================================================"<<endl;
        msg<<"[ttkCinemaQuery] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    string result;
    string sqlTableDefinition, sqlTableRows;

    // -------------------------------------------------------------------------
    // Get Input Table
    // -------------------------------------------------------------------------
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto inTable = vtkTable::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // -------------------------------------------------------------------------
    // Get Output Table
    // -------------------------------------------------------------------------
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto outTable = vtkTable::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    // -------------------------------------------------------------------------
    // Convert Input Table to SQL Table
    // -------------------------------------------------------------------------
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

    // -------------------------------------------------------------------------
    // Replace Variables in QueryString
    // -------------------------------------------------------------------------
    string finalQueryString = this->QueryString;
    {
        vector<vtkFieldData*> candidates;
        candidates.push_back( inTable->GetFieldData() );

        // Add candidates from second optional port
        auto inVector = inputVector[1];
        if(inVector!=nullptr){
            auto inInfo2 = inVector->GetInformationObject(0);

            if(inInfo2!=nullptr){
                auto inObject = inInfo2->Get(vtkDataObject::DATA_OBJECT());
                if(inObject!=nullptr){
                    auto inMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                    if(inObject->IsA("vtkMultiBlockDataSet"))
                        inMB->ShallowCopy( vtkMultiBlockDataSet::SafeDownCast( inObject ) );
                    else
                        inMB->SetBlock(0, inObject );

                    // Iterate over vtkMultiBlockDataSet
                    size_t n = inMB->GetNumberOfBlocks();
                    for(size_t i=0; i<n; i++){
                        auto block = inMB->GetBlock(i);
                        candidates.push_back( block->GetFieldData() );
                        auto blockAsDataSet = vtkDataSet::SafeDownCast( block );
                        if( blockAsDataSet!=nullptr ){
                            candidates.push_back( (vtkFieldData*) blockAsDataSet->GetPointData() );
                            candidates.push_back( (vtkFieldData*) blockAsDataSet->GetCellData() );
                        }
                    }
                }
            }
        }

        // Replace Variables in Query String (e.g. {time[2]})
        {
            vector<string> temp; // vector of substrings of SQL string separated by '{'
            boost::split(temp, finalQueryString, boost::is_any_of("{"));

            for(size_t i=1; i<temp.size(); i++){
                vector<string> temp2; // vector of substrings of temp separated by '}'
                boost::split(temp2, temp[i], boost::is_any_of("}"));

                string varToken = temp2[0]; // e.g. time[2]
                string varName  = varToken; // time
                size_t varIndex = 0;        // 2

                // If varIndex specified in string then update varName and varIndex
                size_t indexDelimiter0 = varToken.find("[");
                size_t indexDelimiter1 = varToken.find("]");
                if( indexDelimiter0!=string::npos && indexDelimiter1!=string::npos ){
                    varName = varToken.substr(0,indexDelimiter0);
                    varIndex = stoi( varToken.substr(indexDelimiter0+1,indexDelimiter1-indexDelimiter0-1) );
                }

                // Search in candidates
                bool found = false;
                for(auto fieldData: candidates){
                    auto column = fieldData->GetAbstractArray( varName.data() );
                    if( column!=nullptr && varIndex<((size_t)column->GetNumberOfTuples()) ){
                        found = true;

                        // Replace Variable
                        boost::replace_all(
                            finalQueryString,
                            "{"+varToken+"}",
                            "\""+column->GetVariantValue( varIndex ).ToString()+"\""
                        );
                    }
                }
                if(!found){
                    // Print Error
                    stringstream msg;
                    msg<<"[ttkCinemaQuery] ERROR: Variable {"<<varToken<<"} not found in field data."<<endl;
                    msg<<"[ttkCinemaQuery] ---------------------------------------------------------------"<<endl;
                    dMsg(cout, msg.str(), timeMsg);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Compute Query Result
    // -------------------------------------------------------------------------
    {
        result = cinemaQuery.execute(
            sqlTableDefinition,
            sqlTableRows,
            finalQueryString
        );
    }

    // -------------------------------------------------------------------------
    // Process Result
    // -------------------------------------------------------------------------
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
                auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
                reader->SetReadFromInputString( true );
                reader->SetInputString( result );
                reader->DetectNumericColumnsOn();
                reader->SetHaveHeaders(true);
                reader->SetFieldDelimiterCharacters(",");
                reader->Update();

                outTable->ShallowCopy( reader->GetOutput() );
            #endif
        } else {
            auto emptyTable = vtkSmartPointer<vtkTable>::New();

            int nc = inTable->GetNumberOfColumns();
            for(int i=0; i<nc; i++){
                auto c = inTable->GetColumn(i);
                auto emptyColumn = vtkSmartPointer<vtkStringArray>::New();
                emptyColumn->SetNumberOfValues(1);
                emptyColumn->SetValue(0,"NULL");
                emptyColumn->SetName( c->GetName() );
                emptyTable->AddColumn( emptyColumn );
            }
            outTable->ShallowCopy( emptyTable );
        }

        outTable->GetFieldData()->ShallowCopy( inTable->GetFieldData() );
    }

    // -------------------------------------------------------------------------
    // Output Performance
    // -------------------------------------------------------------------------
    {
        stringstream msg;
        msg << "[ttkCinemaQuery] Memory usage: "
            << m.getElapsedUsage()
            << " MB." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}
