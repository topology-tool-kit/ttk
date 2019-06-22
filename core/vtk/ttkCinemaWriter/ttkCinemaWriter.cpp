#include <ttkCinemaWriter.h>

#include <vtkVersion.h>

#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStdString.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vtkDirectory.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaWriter)

  int ttkCinemaWriter::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  Timer t;
  double t0 = 0;
  Memory m;
  struct stat info;

  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkCinemaWriter] RequestData" << endl;
    msg << "[ttkCinemaWriter]     Path: " << this->DatabasePath << endl;
    msg << "[ttkCinemaWriter] Override: "
        << (this->OverrideDatabase ? "yes" : "no") << endl;
    msg << "[ttkCinemaWriter] C. Level: " << this->GetCompressLevel() << endl;
    msg << "[ttkCinemaWriter] "
           "--------------------------------------------------------------"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Copy Input to Output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  auto input = inInfo->Get(vtkDataObject::DATA_OBJECT());

  auto inputIsAlreadyMB = input->IsA("vtkMultiBlockDataSet");
  auto inputMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  // If input is vtkMultiBlockDataSet then copy, otherwise add as single block
  if(inputIsAlreadyMB)
    inputMB->ShallowCopy(input);
  else
    inputMB->SetBlock(0, input);

  // -------------------------------------------------------------------------
  // Prepare Database
  // -------------------------------------------------------------------------

  // Check if database path exists and if it has the correct extension
  if(this->DatabasePath.length() < 4
     || this->DatabasePath.substr(this->DatabasePath.length() - 4, 4)
            .compare(".cdb")
          != 0) {
    dMsg(cout,
         "[ttkCinemaWriter] ERROR: Database path has to end with '.cdb'.\n",
         fatalMsg);
    return 0;
  }

  // Initialize path variables
  string dataPrefix = "data/";
  string pathPrefix = this->DatabasePath + "/" + dataPrefix;
  string dataCsvPath = this->DatabasePath + "/data.csv";
  string pathSuffix = ".vtm";

  // Create directory if it does not already exist
  {
    auto directory = vtkSmartPointer<vtkDirectory>::New();
    int opened = directory->Open(this->DatabasePath.data());
    if(!opened) {
      int status = directory->MakeDirectory(this->DatabasePath.data());
      if(status == 1)
        dMsg(cout, "[ttkCinemaWriter] - Directory created\n", infoMsg);
      else {
        dMsg(cout,
             "[ttkCinemaWriter] ERROR: Unable to create database directory.\n",
             fatalMsg);
        return 0;
      }
    }
  }

  // If OverrideDatabase then delete old data products
  if(this->OverrideDatabase) {
    dMsg(cout, "[ttkCinemaWriter] - Deleting old data products        ... ",
         timeMsg);

    t0 = t.getElapsedTime();

    // Delete data.csv
    remove(dataCsvPath.data());

    // Delete data folder
    auto directory = vtkSmartPointer<vtkDirectory>::New();
    if(directory->Open(pathPrefix.data())
       && directory->DeleteDirectory(pathPrefix.data()) == 0)
      dMsg(cout,
           "failed.\n[ttkCinemaWriter] ERROR: Unable to delete existing data "
           "products.\n",
           fatalMsg);
    else {
      stringstream msg;
      msg << "done (" << (t.getElapsedTime() - t0) << " s).\n";
      dMsg(cout, msg.str(), timeMsg);
    }
  }

  // -------------------------------------------------------------------------
  // Store Data products
  // -------------------------------------------------------------------------

  // Determine unique path to new products (for now just generate random number)
  string id;
  string path;
  bool unique = false;

  while(!unique) {
    id = to_string(rand() % 1000000);
    path = pathPrefix + id + pathSuffix;
    unique = stat(path.data(), &info) != 0;
  }

  // Write input to disk
  dMsg(cout, "[ttkCinemaWriter] - Writing new data products to disk ... ",
       timeMsg);

  t0 = t.getElapsedTime();

  auto mbWriter = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
  mbWriter->SetFileName(path.data());
  mbWriter->SetDataModeToAppended();
  mbWriter->SetCompressorTypeToZLib();
  vtkZLibDataCompressor::SafeDownCast(mbWriter->GetCompressor())
    ->SetCompressionLevel(this->GetCompressLevel());
  mbWriter->SetInputData(inputMB);
  mbWriter->Write();

  {
    stringstream msg;
    msg << "done (" << (t.getElapsedTime() - t0) << " s).\n";
    dMsg(cout, msg.str(), timeMsg);
  }

  // -------------------------------------------------------------------------
  // Update 'data.csv' File
  // -------------------------------------------------------------------------

  // Update data.csv file
  dMsg(cout, "[ttkCinemaWriter] - Updating data.csv file            ... ",
       timeMsg);
  t0 = t.getElapsedTime();

  // Create data.csv file if it does not already exist
  if(stat(dataCsvPath.data(), &info) != 0) {
    vtkDataObject *firstBlock = inputMB->GetBlock(0);

    if(firstBlock != nullptr) {
      // Open file
      ofstream csvFile;
      csvFile.open(dataCsvPath.data());
      if(!csvFile.is_open()) {
        dMsg(cout,
             "failed.\n[ttkCinemaWriter] ERROR: Unable to create 'data.csv' "
             "file.\n",
             fatalMsg);
        return 0;
      }

      // Get column names
      auto fieldData = firstBlock->GetFieldData();
      size_t n = fieldData->GetNumberOfArrays();
      for(size_t i = 0; i < n; i++) {
        auto array = fieldData->GetAbstractArray(i);
        string name = array->GetName();
        if(array->GetNumberOfTuples() == 1 && name.compare("FILE") != 0)
          csvFile << name << ",";
      }
      csvFile << "FILE\n";

      // Close file
      csvFile.close();
    }
  }

  // Read csv file
  auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
  reader->SetFileName(dataCsvPath.data());
  reader->DetectNumericColumnsOff();
  reader->SetHaveHeaders(true);
  reader->SetFieldDelimiterCharacters(",");
  reader->Update();
  auto table = vtkTable::SafeDownCast(reader->GetOutput());

  int n = inputMB->GetNumberOfBlocks();
  int offset = this->OverrideDatabase ? 0 : table->GetNumberOfRows();
  if(this->OverrideDatabase)
    table->SetNumberOfRows(n);
  else
    for(int i = 0; i < n;
        i++) // SetNumberOfRows clears existing rows -> insert blank rows
      table->InsertNextBlankRow();

  for(int i = 0; i < n; i++) {
    auto block = inputMB->GetBlock(i);
    string blockExtension = "vtk";
#if VTK_MAJOR_VERSION <= 7
    stringstream msg;
    msg << "failed." << endl
        << "[ttkCinemaQuery] ERROR: VTK version too old." << endl
        << "[ttkCinemaQuery]        This filter requires "
           "vtkXMLPMultiBlockDataWriter"
        << endl
        << "[ttkCinemaQuery]        of version 7.0 or higher." << endl;
    dMsg(cout, msg.str(), fatalMsg);
    return 0;
#else
    blockExtension
      = this->GetDefaultFileExtensionForDataSet(block->GetDataObjectType());
#endif

    auto fieldData = block->GetFieldData();

    for(int j = 0; j < table->GetNumberOfColumns(); j++) {
      auto columnName = table->GetColumnName(j);
      auto columnCSV = vtkStringArray::SafeDownCast(table->GetColumn(j));
      if(string(columnName).compare("FILE") == 0)
        columnCSV->SetValue(
          offset + i, vtkStdString(dataPrefix + id + "/" + id + "_"
                                   + to_string(i) + "." + blockExtension));
      else {
        auto array = fieldData->GetAbstractArray(columnName);
        string value = array->GetVariantValue(0).ToString();
        columnCSV->SetValue(offset + i, value);
      }
    }
  }

  // Write data.csv file
  auto csvWriter = vtkSmartPointer<vtkDelimitedTextWriter>::New();
  csvWriter->SetUseStringDelimiter(false);
  csvWriter->SetFileName((this->DatabasePath + "/data.csv").data());
  csvWriter->SetInputData(table);
  csvWriter->Write();

  {
    stringstream msg;
    msg << "done (" << (t.getElapsedTime() - t0) << " s).\n";
    dMsg(cout, msg.str(), timeMsg);
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkCinemaWriter] "
           "--------------------------------------------------------------"
        << endl;
    msg << "[ttkCinemaWriter]   Time: " << t.getElapsedTime() << " s." << endl;
    msg << "[ttkCinemaWriter] Memory: " << m.getElapsedUsage() << " MB."
        << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
