#include <ttkAddFieldData.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkDelimitedTextReader.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkAddFieldData)

int ttkAddFieldData::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkAddFieldData] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Memory m;

    // Prepare Input and Output
    auto source = inputVector[0]->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT());
    auto target = inputVector[1]->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT());


    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = vtkTable::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

    output->ShallowCopy(target);

    auto outputFieldData = output->GetFieldData();

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
            while( getline(ss1, element, ',') ){
                result << element << endl;
            }

            // Parse csv string with reader
            reader->SetInputString( result.str() );
            reader->Update();
            auto csvTable = reader->GetOutput();

            // Add column to field data
            size_t n = csvTable->GetNumberOfColumns(); // n should always be 1
            for(size_t i=0; i<n; i++){
                auto column = csvTable->GetColumn(i);
                if(column->GetNumberOfValues()>0)
                    outputFieldData->AddArray( csvTable->GetColumn(i) );
            }
        }
    }

    // If source is a vtkDataSet and is not the target then add source point/cell/field data to target field data
    if( source->IsA("vtkDataSet") && source!=target ){
        auto sourceAsVtkDataSet = vtkDataSet::SafeDownCast( source );
        auto pointData = sourceAsVtkDataSet->GetPointData();

        size_t n = pointData->GetNumberOfArrays(); // n should always be 1
        for(size_t i=0; i<n; i++)
            outputFieldData->AddArray( pointData->GetArray(i) );
    }

    return 1;
}