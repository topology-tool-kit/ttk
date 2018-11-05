#include <ttkAddFieldData.h>

#include <vtkDataSet.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkFieldData.h>
#include <vtkDelimitedTextReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkAddFieldData)

void ttkAddFieldData::AddArray(int fieldType, const char* name){
    if(!name){
        vtkErrorMacro("name cannont be null.");
        return;
    }
    string n = name;
    this->ArraySelection.push_back(
        make_pair( fieldType, n )
    );
    this->Modified();
}

void ttkAddFieldData::AddPointDataArray(const char* name){
  this->AddArray(vtkDataObject::POINT, name);
}

void ttkAddFieldData::AddCellDataArray(const char* name){
  this->AddArray(vtkDataObject::CELL, name);
}

void ttkAddFieldData::AddFieldDataArray(const char* name){
  this->AddArray(vtkDataObject::FIELD, name);
}

void ttkAddFieldData::RemoveArray(int fieldType, const char* name, bool deleteType){
    bool found = true;
    while( found ){
        found = false;
        auto iter = this->ArraySelection.begin();
        while( iter!=this->ArraySelection.end() ){
            if( iter->first==fieldType && (deleteType || iter->second==name) ){
                found = true;
                iter = this->ArraySelection.erase(iter);
                this->Modified();
            } else ++iter;
        }
    }
}

void ttkAddFieldData::RemovePointDataArray(const char* name){
  this->RemoveArray(vtkDataObject::POINT, name, false);
}

void ttkAddFieldData::RemoveCellDataArray(const char* name){
  this->RemoveArray(vtkDataObject::CELL, name, false);
}

void ttkAddFieldData::RemoveFieldDataArray(const char* name){
  this->RemoveArray(vtkDataObject::FIELD, name, false);
}

void ttkAddFieldData::ClearArrays(){
  if(this->ArraySelection.size()>0)
    this->Modified();
  this->ArraySelection.clear();
}

void ttkAddFieldData::ClearPointDataArrays(){
    this->RemoveArray(vtkDataObject::POINT, "", true);
    this->Modified();
}

void ttkAddFieldData::ClearCellDataArrays(){
    this->RemoveArray(vtkDataObject::CELL, "", true);
    this->Modified();
}

void ttkAddFieldData::ClearFieldDataArrays(){
    this->RemoveArray(vtkDataObject::FIELD, "", true);
    this->Modified();
}


int ttkAddFieldData::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Print Status
    {
        stringstream msg;
        msg<<"================================================================================"<<endl;
        msg<<"[ttkAddFieldData] RequestData"<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    Timer t;

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());

    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
    output->ShallowCopy(input);

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

    // Optionally: copy data from second port
    {
        if(this->ArraySelection.size>0){
            auto inInfo = inputVector[1]->GetInformationObject(0);
            if(inInfo){
                auto input = inInfo1->Get(vtkDataObject::DATA_OBJECT());
                if(input){

                }
            }

        }
    }
    // for(auto x: this->ArraySelection)
    //     cout<<x.first<<" "<<x.second<<endl;

    // Output Performance
    {
        stringstream msg;
        msg << "[ttkAddFieldData] ArrayData created in " << t.getElapsedTime() << " s." << endl;
        dMsg(cout, msg.str(), memoryMsg);
    }

    return 1;
}