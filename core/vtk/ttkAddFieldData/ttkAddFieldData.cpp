#include <ttkAddFieldData.h>

#include <vtkVersion.h>

#include <vtkDataSet.h>
#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkAddFieldData)

  // =============================================================================
  // GUI functions to add/remove/clear array selections
  // =============================================================================

  // Add Array Selections
  void ttkAddFieldData::AddArray(int fieldType, const char *name) {
  if(!name) {
    vtkErrorMacro("name cannont be null.");
    return;
  }
  string n = name;
  this->ArraySelection.push_back(make_pair(fieldType, n));
  this->Modified();
}
void ttkAddFieldData::AddPointDataArray(const char *name) {
  this->AddArray(vtkDataObject::POINT, name);
}
void ttkAddFieldData::AddCellDataArray(const char *name) {
  this->AddArray(vtkDataObject::CELL, name);
}
void ttkAddFieldData::AddFieldDataArray(const char *name) {
  this->AddArray(vtkDataObject::FIELD, name);
}

// Remove Array Selections
void ttkAddFieldData::RemoveArray(int fieldType,
                                  const char *name,
                                  bool deleteType) {
  bool found = true;
  while(found) {
    found = false;
    auto iter = this->ArraySelection.begin();
    while(iter != this->ArraySelection.end()) {
      if(iter->first == fieldType && (deleteType || iter->second == name)) {
        found = true;
        iter = this->ArraySelection.erase(iter);
        this->Modified();
      } else
        ++iter;
    }
  }
}
void ttkAddFieldData::RemovePointDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::POINT, name, false);
}
void ttkAddFieldData::RemoveCellDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::CELL, name, false);
}
void ttkAddFieldData::RemoveFieldDataArray(const char *name) {
  this->RemoveArray(vtkDataObject::FIELD, name, false);
}

// Clear Array Selections
void ttkAddFieldData::ClearArrays() {
  if(this->ArraySelection.size() > 0)
    this->Modified();
  this->ArraySelection.clear();
}
void ttkAddFieldData::ClearPointDataArrays() {
  this->RemoveArray(vtkDataObject::POINT, "", true);
}
void ttkAddFieldData::ClearCellDataArrays() {
  this->RemoveArray(vtkDataObject::CELL, "", true);
}
void ttkAddFieldData::ClearFieldDataArrays() {
  this->RemoveArray(vtkDataObject::FIELD, "", true);
}

// =============================================================================
// RequestData
// =============================================================================
int ttkAddFieldData::RequestData(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkAddFieldData] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  Timer t;

  // Pass Input to Output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto output = outInfo->Get(vtkDataObject::DATA_OBJECT());
  output->ShallowCopy(input);

  // Get FieldData of output
  auto outputFD = output->GetFieldData();

  // If FieldDataString is not empty then parse and add as new field data
  if(this->FieldDataString.length() > 0) {
    // Print Status
    dMsg(
      cout,
      "[ttkAddFieldData] Adding field data parsed from string           ... ",
      timeMsg);

#if VTK_MAJOR_VERSION <= 7
    stringstream msg;
    msg
      << "failed\n[ttkAddFieldData] ERROR: VTK version too old." << endl
      << "[ttkAddFieldData]        This filter requires vtkDelimitedTextReader"
      << endl
      << "[ttkAddFieldData]        of version 7.0 or higher." << endl;
    dMsg(cout, msg.str(), fatalMsg);
    return 0;
#else

    // Initialize reader that will be used to parse columns
    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetReadFromInputString(true);
    reader->SetHaveHeaders(true);
    reader->SetDetectNumericColumns(true);
    reader->SetFieldDelimiterCharacters(",");

    // Iterate over lines
    stringstream ss0(this->FieldDataString);
    string line;
    string element;

    // For each line
    while(getline(ss0, line, '\n')) {

      // Generate csv string that has only one column
      stringstream result;
      stringstream ss1(line);
      while(getline(ss1, element, ','))
        result << element << endl;

      // Parse csv string with reader
      reader->SetInputString(result.str());
      reader->Update();
      auto csvTable = reader->GetOutput();

      // Add column to field data
      size_t n = csvTable->GetNumberOfColumns(); // n should always be 1
      for(size_t i = 0; i < n; i++) {
        auto column = csvTable->GetColumn(i);
        if(column->GetNumberOfTuples() > 0)
          outputFD->AddArray(csvTable->GetColumn(i));
      }
    }

#endif

    dMsg(cout, "done\n", timeMsg);
  }

  // Optionally: copy data from second port
  if(this->ArraySelection.size() > 0) {
    // Print Status
    dMsg(
      cout,
      "[ttkAddFieldData] Adding point/cell/field data from second input ... ",
      timeMsg);

    // Check if second input exists
    auto inInfoObj = inputVector[1]->GetInformationObject(0);
    if(inInfoObj) {
      // Vector that stores all vtkFieldData sources
      vector<vtkFieldData *> candidates(3);

      // Get Second Input
      auto input2 = inInfoObj->Get(vtkDataObject::DATA_OBJECT());

      // If second input exists add its field data to candidates
      if(input2) {
        candidates[2] = input2->GetFieldData();

        // If second input is a vtkDataSet then also add point/cell data to
        // candidates
        if(input2->IsA("vtkDataSet")) {
          auto inputAsDataSet = vtkDataSet::SafeDownCast(input2);
          candidates[0] = (vtkFieldData *)inputAsDataSet->GetPointData();
          candidates[1] = (vtkFieldData *)inputAsDataSet->GetCellData();
        }
      }

      // Iterate over ArraySelection and select array from candidates
      // accordingly
      for(auto &x : this->ArraySelection) {
        bool found = false;
        vtkFieldData *fd = candidates[x.first];
        if(fd != nullptr) {
          auto array = fd->GetAbstractArray(x.second.data());
          if(array != nullptr) {
            outputFD->AddArray(array);
            found = true;
          }
        }
        if(!found) {
          stringstream msg;
          msg << "failed\n[ttkAddFieldData] ERROR: Point/Cell/Field data not "
                 "found: ("
              << x.first << " -> " << x.second << ")" << endl;
          dMsg(cout, msg.str(), fatalMsg);
          return 0;
        }
      }
    }

    dMsg(cout, "done\n", timeMsg);
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkAddFieldData] "
           "--------------------------------------------------------------"
        << endl;
    msg << "[ttkAddFieldData] time: " << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
