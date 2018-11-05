#include <ttkTableSource.h>

#include <vtkArrayData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkFieldData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkTableSource)

int ttkTableSource::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkTableSource] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Timer t;

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkTable::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );
    auto outputFD = output->GetFieldData();

    // If FieldDataString is not empty then parse and add as new field data
    if( this->FieldDataString.length()>0 ){

        // Initialize reader that will be used to parse columns
        auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
        reader->SetReadFromInputString( true );
        reader->SetHaveHeaders( true );
        reader->SetDetectNumericColumns( true );
        reader->SetFieldDelimiterCharacters( "," );

        // Iterate over lines
        stringstream ss0( this->FieldDataString );
        string line;
        string element;

        // For each line
        while( getline(ss0, line, '\n') ){

            // Generate csv string that has only one column
            stringstream result;
            stringstream ss1( line );
            while( getline(ss1, element, ',') )
                result << element << endl;

            // Parse csv string with reader
            reader->SetInputString( result.str() );
            reader->Update();
            auto csvTable = reader->GetOutput();

            // Add column to field data
            size_t n = csvTable->GetNumberOfColumns(); // n should always be 1
            for(size_t i=0; i<n; i++){
                auto column = csvTable->GetColumn(i);
                if(column->GetNumberOfTuples()>0)
                    outputFD->AddArray( csvTable->GetColumn(i) );
            }
        }
    }

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkTableSource] ArrayData created in " << t.getElapsedTime() << " s." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}