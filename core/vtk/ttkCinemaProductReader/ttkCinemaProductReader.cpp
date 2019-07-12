#include <ttkCinemaProductReader.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>
#include <vtkXMLGenericDataObjectReader.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaProductReader)

  int ttkCinemaProductReader::RequestData(vtkInformation *request,
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
  // Print status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkCinemaProductReader] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  Timer t;
  Memory mem;

  // Prepare Input and Output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto inputTable
    = vtkTable::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Read Data
  {
    // Determine number of files
    size_t n = inputTable->GetNumberOfRows();
    size_t m = inputTable->GetNumberOfColumns();

    {
      stringstream msg;
      msg << "[ttkCinemaProductReader] Reading " << n << " files:" << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    // Compute DatabasePath
    auto databasePath = inputTable->GetFieldData()
                          ->GetAbstractArray("DatabasePath")
                          ->GetVariantValue(0)
                          .ToString();

    // Get column that contains paths
    auto paths = inputTable->GetColumnByName(this->FilepathColumnName.data());
    if(paths == nullptr) {
      stringstream msg;
      msg << "[ttkCinemaProductReader] ERROR: Table does not have "
             "FilepathColumn '"
          << this->FilepathColumnName << "'" << endl;
      dMsg(cerr, msg.str(), fatalMsg);
      return 0;
    }

    // For each row
    for(size_t i = 0; i < n; i++) {
      // Get path
      auto path = databasePath + "/" + paths->GetVariantValue(i).ToString();
      auto ext = path.substr(path.length() - 3);

      {
        stringstream msg;
        msg << "[ttkCinemaProductReader]    " << i << ": " << path << endl;
        dMsg(cout, msg.str(), infoMsg);
      }

      // Check if file exists
      std::ifstream infile(path.data());
      bool exists = infile.good();
      if(!exists) {
        stringstream msg;
        msg << "[ttkCinemaProductReader]    ERROR: File does not exist."
            << endl;
        dMsg(cerr, msg.str(), fatalMsg);
        continue;
      }

      // Read any data using vtkXMLGenericDataObjectReader
      {
        auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
        reader->SetFileName(path.data());
        reader->Update();
        output->SetBlock(i, reader->GetOutput());
      }

      // Augment read data with row information
      // TODO: Make Optional
      auto block = output->GetBlock(i);
      for(size_t j = 0; j < m; j++) {
        auto columnName = inputTable->GetColumnName(j);
        auto fieldData = block->GetFieldData();
        if(!fieldData->HasArray(columnName)) {
          bool isNumeric = inputTable->GetColumn(j)->IsNumeric();

          if(isNumeric) {
            auto c = vtkSmartPointer<vtkDoubleArray>::New();
            c->SetName(columnName);
            c->SetNumberOfValues(1);
            c->SetValue(0, inputTable->GetValue(i, j).ToDouble());
            fieldData->AddArray(c);
          } else {
            auto c = vtkSmartPointer<vtkStringArray>::New();
            c->SetName(columnName);
            c->SetNumberOfValues(1);
            c->SetValue(0, inputTable->GetValue(i, j).ToString());
            fieldData->AddArray(c);
          }
        }
      }

      this->updateProgress(((float)i) / ((float)(n - 1)));
    }
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkCinemaProductReader] "
           "-------------------------------------------------------"
        << endl;
    msg << "[ttkCinemaProductReader]   time: " << t.getElapsedTime() << " s."
        << endl;
    msg << "[ttkCinemaProductReader] memory: " << mem.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
