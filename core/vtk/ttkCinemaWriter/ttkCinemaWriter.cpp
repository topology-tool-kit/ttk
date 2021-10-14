#include <ttkCinemaWriter.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDirectory.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStdString.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

// writers common
#include <vtkZLibDataCompressor.h>

// CSV writers
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkMultiBlockDataSet.h>

// product writers
#include <vtkPNGWriter.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>

// file lock
#include <boost/interprocess/sync/file_lock.hpp>
#include <sys/stat.h>

vtkStandardNewMacro(ttkCinemaWriter);

ttkCinemaWriter::ttkCinemaWriter() {
  this->setDebugMsgPrefix("CinemaWriter");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaWriter::~ttkCinemaWriter() {
}

int ttkCinemaWriter::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  } else {
    return 0;
  }
  return 1;
}

int ttkCinemaWriter::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  } else {
    return 0;
  }
  return 1;
}

int ensureFolder(const std::string &path) {
  auto directory = vtkSmartPointer<vtkDirectory>::New();
  if(directory->Open(path.data()) == 1
     || vtkDirectory::MakeDirectory(path.data()) == 1)
    return 1;
  else
    return 0;
}

int ttkCinemaWriter::ValidateDatabasePath() {
  if(this->DatabasePath.length() < 4
     || this->DatabasePath.substr(this->DatabasePath.length() - 4, 4)
            .compare(".cdb")
          != 0) {
    this->printErr("Database path has to end with '.cdb'.");
    return 0;
  }

  return 1;
}

int ttkCinemaWriter::DeleteDatabase() {
  ttk::Timer t;
  this->printMsg("Deleting CDB: " + this->DatabasePath, 0,
                 ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

  this->Modified();
  if(this->ValidateDatabasePath() == 0)
    return 0;
  int status = vtkDirectory::DeleteDirectory(this->DatabasePath.data());

  this->printMsg("Deleting CDB: " + this->DatabasePath, 1, t.getElapsedTime());

  return status;
}

int ttkCinemaWriter::GetLockFilePath(std::string &path) {
  if(!this->ValidateDatabasePath())
    return 0;

  path = this->DatabasePath + ".lockfile";

  return 1;
}

int ttkCinemaWriter::InitializeLockFile() {
  std::string lockFilePath;
  if(!this->GetLockFilePath(lockFilePath))
    return 0;

  std::ofstream output(lockFilePath);
  output.close();

  return 1;
}

// =============================================================================
// Process Request
// =============================================================================
int ttkCinemaWriter::ProcessDataProduct(vtkDataObject *input) {

  // ---------------------------------------------------------------------------
  // Get Correct Data Product Extension
  // ---------------------------------------------------------------------------
  vtkSmartPointer<vtkXMLWriter> xmlWriter;
  if(input->IsA("vtkDataSet")) {
    xmlWriter = vtkSmartPointer<vtkXMLWriter>::Take(
      vtkXMLDataObjectWriter::NewWriter(input->GetDataObjectType()));
  } else if(input->IsA("vtkMultiBlockDataSet")) {
    xmlWriter = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
  } else {
    return 0;
  }

  xmlWriter->SetDataModeToAppended();
  xmlWriter->SetCompressorTypeToZLib();
  const auto compressor
    = vtkZLibDataCompressor::SafeDownCast(xmlWriter->GetCompressor());
  if(compressor != nullptr)
    compressor->SetCompressionLevel(this->CompressionLevel);

  std::string productExtension = this->Format == FORMAT::VTK
                                   ? xmlWriter->GetDefaultFileExtension()
                                 : this->Format == FORMAT::PNG ? "png"
                                                               : "ttk";

  // -------------------------------------------------------------------------
  // Prepare Field Data
  // -------------------------------------------------------------------------
  auto inputFD = vtkSmartPointer<vtkFieldData>::New();
  inputFD->ShallowCopy(input->GetFieldData());
  size_t nFields = inputFD->GetNumberOfArrays();

  // -------------------------------------------------------------------------
  // Ignore Meta Fields
  // -------------------------------------------------------------------------
  {
    std::vector<std::string> toIgnore;

    // ignore file column
    toIgnore.emplace_back("FILE");

    // remove temporary columns
    for(size_t i = 0; i < nFields; i++) {
      std::string name(inputFD->GetArrayName(i));
      if(name.substr(0, 4).compare("_ttk") == 0)
        toIgnore.emplace_back(name);
    }

    // delete columns from fd
    for(const auto &name : toIgnore)
      inputFD->RemoveArray(name.data());

    nFields = inputFD->GetNumberOfArrays();
  }

  // ===========================================================================
  // Determine ProductId and collect values
  // ===========================================================================
  std::string productId;
  std::string rDataProductPath;
  std::vector<std::string> fields;
  std::vector<std::string> values;
  {

    if(nFields < 1) {
      productId = "FILE";
    } else {
      for(size_t i = 0; i < nFields; i++) {
        auto array = inputFD->GetAbstractArray(i);
        auto n = array->GetNumberOfTuples();
        auto m = array->GetNumberOfComponents();
        std::string value;
        if(n < 1) {
          value = "null";
        } else {
          value = "";
          for(int j = 0; j < n; j++) {
            for(int k = 0; k < m; k++) {
              value += array->GetVariantValue(j * m + k).ToString() + "|";
            }
            value = value.substr(0, value.size() - 1);
            value += ";";
          }
          value = value.substr(0, value.size() - 1);
        }

        fields.emplace_back(array->GetName());
        values.emplace_back(value);
      }

      productId = values[0];
      for(size_t i = 1; i < nFields; i++)
        productId += "_" + values[i];
    }

    rDataProductPath = "data/" + productId + "." + productExtension;
  }

  // print keys
  {
    std::vector<std::vector<std::string>> rows(nFields + 1);
    for(size_t i = 0; i < nFields; i++)
      rows[i] = {fields[i], values[i]};
    rows[nFields] = {"Key", productId};

    this->printMsg(rows, ttk::debug::Priority::VERBOSE);
  }

  // ===========================================================================
  // Update database
  // ===========================================================================
  {
    // Initilize file lock for remaining operations
    std::string lockFilePath;
    if(!this->GetLockFilePath(lockFilePath))
      return 0;

    boost::interprocess::file_lock flock;
    try {
      flock = boost::interprocess::file_lock(lockFilePath.data());
      flock.lock();
    } catch(boost::interprocess::interprocess_exception &) {
    }

    std::string csvPath = this->DatabasePath + "/data.csv";
    struct stat info;

    // -------------------------------------------------------------------------
    // If data.csv file does not exsist create it
    // -------------------------------------------------------------------------
    if(stat(csvPath.data(), &info) != 0) {
      ttk::Timer t;
      this->printMsg("Creating data.csv file", 0, ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);

      std::ofstream csvFile;
      csvFile.open(csvPath.data());
      if(!csvFile.is_open()) {
        this->printErr("Unable to create 'data.csv' file.");
        return 0;
      }

      std::string header;
      std::string firstRow;
      for(size_t i = 0; i < nFields; i++) {
        header += fields[i] + ",";
        firstRow += values[i] + ",";
      }
      header += "FILE";
      firstRow += rDataProductPath;

      csvFile << header << endl << firstRow << endl;

      // Close file
      csvFile.close();

      this->printMsg("Creating data.csv file", 1, t.getElapsedTime(),
                     ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
    }

    // -------------------------------------------------------------------------
    // Update data.csv file
    // -------------------------------------------------------------------------

    // read data.csv file
    auto csvTable = vtkSmartPointer<vtkTable>::New();
    {
      ttk::Timer t;
      this->printMsg("Reading data.csv file", 0, ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);

      auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
      reader->SetFileName(csvPath.data());
      reader->DetectNumericColumnsOff();
      reader->SetHaveHeaders(true);
      reader->SetFieldDelimiterCharacters(",");
      reader->Update();
      auto readerOutput = vtkTable::SafeDownCast(reader->GetOutput());
      if(!readerOutput) {
        this->printErr("Unable to read 'data.csv' file.");
        return 0;
      }

      csvTable->ShallowCopy(readerOutput);

      this->printMsg("Reading data.csv file", 1, t.getElapsedTime(),
                     ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
    }

    // check CSV file integrity
    std::vector<vtkStringArray *> fieldToCSVColumnMap(nFields);
    size_t nRows = csvTable->GetNumberOfRows();
    size_t nColumns = csvTable->GetNumberOfColumns();
    {
      // Check If CSV file is empty
      if(nColumns == 0) {
        this->printErr(
          "Empty 'data.csv' file (vtkDelimitedTextReader limitation).");
        return 0;
      }

      // Check If CSV file contains columns not present in fields
      for(size_t i = 0; i < nColumns; i++) {
        std::string columnName = csvTable->GetColumnName(i);

        // skip FILE column
        if(columnName.compare("FILE") == 0)
          continue;

        bool exsist = false;
        for(size_t j = 0; j < nFields; j++) {
          if(fields[j].compare(columnName) == 0)
            exsist = true;
        }
        if(!exsist) {
          this->printErr("'data.csv' file contains column '" + columnName
                         + "' not present in field data.");
          this->printErr("Unable to insert data product into cinema database.");
          return 0;
        }
      }

      // check if data product has fields not present in the CSV file
      for(size_t i = 0; i < nFields; i++) {
        auto column = vtkStringArray::SafeDownCast(
          csvTable->GetColumnByName(fields[i].data()));
        if(!column) {
          this->printErr("Data product has field data array '" + fields[i]
                         + "' no recorded in the data.csv file.");
          return 0;
        }
        fieldToCSVColumnMap[i] = column;
      }
    }

    // TODO: make dynamic
    auto fileColumn
      = vtkStringArray::SafeDownCast(csvTable->GetColumnByName("FILE"));
    if(!fileColumn) {
      this->printErr("'data.csv' file has no 'FILE' column");
      return 0;
    }

    // -------------------------------------------------------------------------
    // delete products of the database with same keys
    // -------------------------------------------------------------------------
    {
      std::vector<int> rowsToDelete;
      for(size_t i = 0; i < nRows; i++) {
        auto equal = true;
        for(size_t j = 0; j < nFields; j++)
          if(values[j].compare(fieldToCSVColumnMap[j]->GetValue(i)) != 0)
            equal = false;
        if(equal)
          rowsToDelete.emplace_back(i);
      }

      if(rowsToDelete.size() > 0) {
        ttk::Timer t;
        this->printMsg("Deleting products with same keys", 0,
                       ttk::debug::LineMode::REPLACE,
                       ttk::debug::Priority::DETAIL);

        for(int i = rowsToDelete.size() - 1; i >= 0; i--) {
          auto path = fileColumn->GetValue(rowsToDelete[i]);

          // Remove DataProduct
          remove((this->DatabasePath + "/" + path).data());

          // Remove Row from CSV
          csvTable->RemoveRow(rowsToDelete[i]);
        }

        this->printMsg("Deleting products with same keys", 1,
                       t.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
      }
    }

    // -----------------------------------------------------------------
    // Update data.csv file
    // -----------------------------------------------------------------
    {
      ttk::Timer t;
      this->printMsg("Updating data.csv file", 0, ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);

      size_t rowIndex = csvTable->GetNumberOfRows();
      csvTable->InsertNextBlankRow();

      for(size_t j = 0; j < nFields; j++)
        fieldToCSVColumnMap[j]->SetValue(rowIndex, values[j]);

      fileColumn->SetValue(rowIndex, rDataProductPath);

      // Write data.csv file
      auto csvWriter = vtkSmartPointer<vtkDelimitedTextWriter>::New();
      csvWriter->SetUseStringDelimiter(false);
      csvWriter->SetFileName((this->DatabasePath + "/data.csv").data());
      csvWriter->SetInputData(csvTable);
      csvWriter->Write();

      this->printMsg("Updating data.csv file", 1, t.getElapsedTime(),
                     ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
    }
  }

  // =========================================================================
  // Store Data products
  // =========================================================================
  {
    // Write input to disk
    ttk::Timer t;
    this->printMsg("Writing data product to disk", 0,
                   ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

    switch(this->Format) {

      case FORMAT::VTK: {
        xmlWriter->SetFileName(
          (this->DatabasePath + "/" + rDataProductPath).data());
        xmlWriter->SetInputData(input);
        xmlWriter->Write();
        break;
      }
      case FORMAT::PNG: {
        auto inputAsID = vtkImageData::SafeDownCast(input);
        if(!inputAsID) {
          this->printErr("PNG format requires input of type 'vtkImageData'.");
          return 0;
        }

        // search color array
        {
          bool found = false;
          auto inputPD = inputAsID->GetPointData();
          for(int i = 0; i < inputPD->GetNumberOfArrays(); i++) {
            auto array = inputPD->GetAbstractArray(i);
            if(array->IsA("vtkUnsignedCharArray")) {
              inputPD->SetActiveScalars(inputPD->GetArrayName(i));
              found = true;
              break;
            }
          }

          if(!found) {
            this->printErr("Input image does not have any color array.");
            return 0;
          }
        }

        auto imageWriter = vtkSmartPointer<vtkPNGWriter>::New();
        imageWriter->SetCompressionLevel(this->CompressionLevel);
        imageWriter->SetFileName(
          (this->DatabasePath + "/" + rDataProductPath).data());
        imageWriter->SetInputData(inputAsID);
        imageWriter->Write();
        break;
      }
      case FORMAT::TTK: {
        // Topological Compression
        if(!input->IsA("vtkImageData")) {
          vtkErrorMacro(
            "Cannot use Topological Compression without a vtkImageData");
          return 0;
        }

        const auto inputData = vtkImageData::SafeDownCast(input);
        const auto sf = this->GetInputArrayToProcess(0, inputData);

        vtkNew<ttkTopologicalCompressionWriter> topologicalCompressionWriter{};
        topologicalCompressionWriter->SetInputArrayToProcess(
          0, 0, 0, 0, sf->GetName());

        topologicalCompressionWriter->SetTolerance(this->Tolerance);
        topologicalCompressionWriter->SetMaximumError(this->MaximumError);
        topologicalCompressionWriter->SetZFPTolerance(this->ZFPTolerance);
        topologicalCompressionWriter->SetCompressionType(this->CompressionType);
        topologicalCompressionWriter->SetSQMethodPV(this->SQMethodPV);
        topologicalCompressionWriter->SetZFPOnly(this->ZFPOnly);
        topologicalCompressionWriter->SetSubdivide(this->Subdivide);
        topologicalCompressionWriter->SetUseTopologicalSimplification(
          this->UseTopologicalSimplification);

        // Check that input scalar field is indeed scalar
        if(sf->GetNumberOfComponents() != 1) {
          vtkErrorMacro("Input scalar field should have only 1 component");
          return 0;
        }
        topologicalCompressionWriter->SetDebugLevel(this->debugLevel_);
        topologicalCompressionWriter->SetFileName(
          (this->DatabasePath + "/" + rDataProductPath).data());
        topologicalCompressionWriter->SetInputData(inputData);
        topologicalCompressionWriter->Write();
        break;
      }
      default:
        this->printErr("Unsupported Format");
        return 0;
    }

    this->printMsg("Writing data product to disk", 1, t.getElapsedTime(),
                   ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
  }
  this->printMsg("Wrote " + productId + "." + productExtension);
  this->printMsg(ttk::debug::Separator::L2, ttk::debug::Priority::DETAIL);
  return 1;
}

int ttkCinemaWriter::RequestData(vtkInformation *ttkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  ttk::Timer timer;

  // Print Status
  {
    std::string format = this->Format == FORMAT::VTK   ? "VTK"
                         : this->Format == FORMAT::PNG ? "PNG"
                                                       : "TTK";
    this->printMsg({{"Database", this->DatabasePath},
                    {"C. Level", std::to_string(this->CompressionLevel)},
                    {"Format", format},
                    {"Iterate", this->IterateMultiBlock ? "Yes" : "No"}});
    this->printMsg(ttk::debug::Separator::L1);
  }

  // -------------------------------------------------------------------------
  // Copy Input to Output
  // -------------------------------------------------------------------------
  auto input = vtkDataObject::GetData(inputVector[0]);
  auto output = vtkDataObject::GetData(outputVector);
  if(this->ForwardInput)
    output->ShallowCopy(input);

  // -------------------------------------------------------------------------
  // Prepare Database
  // -------------------------------------------------------------------------
  {
    // Initilize file lock for remaining operations
    std::string lockFilePath;
    if(!this->GetLockFilePath(lockFilePath))
      return 0;

    boost::interprocess::file_lock flock;
    try {
      flock = boost::interprocess::file_lock(lockFilePath.data());
      flock.lock();
    } catch(boost::interprocess::interprocess_exception &) {
    }

    if(this->ValidateDatabasePath() == 0)
      return 0;

    if(ensureFolder(this->DatabasePath) == 0) {
      this->printErr("Unable to open/create cinema database.");
      return 0;
    }

    if(ensureFolder(this->DatabasePath + "/data") == 0) {
      this->printErr("Unable to open/create cinema database.");
      return 0;
    }
  }

  auto inputAsMB = vtkMultiBlockDataSet::SafeDownCast(input);
  if(this->IterateMultiBlock && inputAsMB) {
    size_t n = inputAsMB->GetNumberOfBlocks();
    for(size_t i = 0; i < n; i++)
      if(!this->ProcessDataProduct(inputAsMB->GetBlock(i)))
        return 0;
  } else if(!this->ProcessDataProduct(input))
    return 0;

  // Output Performance
  {
    std::string resultString = "Complete (#products: ";
    resultString += !this->IterateMultiBlock || !inputAsMB
                      ? "1"
                      : std::to_string(inputAsMB->GetNumberOfBlocks());
    resultString += ")";

    // print stats
    this->printMsg(resultString, 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
  }

  return 1;
}
