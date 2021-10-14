#include <ttkCinemaProductReader.h>

#include <vtkInformation.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkXMLGenericDataObjectReader.h>

vtkStandardNewMacro(ttkCinemaProductReader);

ttkCinemaProductReader::ttkCinemaProductReader() {
  this->setDebugMsgPrefix("CinemaProductReader");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkCinemaProductReader::~ttkCinemaProductReader() {
}

int ttkCinemaProductReader::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
  else
    return 0;
  return 1;
}

int ttkCinemaProductReader::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  else
    return 0;
  return 1;
}

template <class readerT>
vtkSmartPointer<vtkDataObject> readFileLocal_(const std::string &pathToFile,
                                              vtkNew<readerT> &reader) {
  reader->SetFileName(pathToFile.data());
  reader->Update();
  if(reader->GetErrorCode() != 0)
    return nullptr;

  auto result
    = vtkSmartPointer<vtkDataObject>::Take(reader->GetOutput()->NewInstance());
  result->ShallowCopy(reader->GetOutput());
  return result;
}

vtkSmartPointer<vtkDataObject>
  ttkCinemaProductReader::readFileLocal(const std::string &pathToFile) {

  if(pathToFile.substr(pathToFile.length() - 4, 4).compare(".ttk") == 0) {
    this->topologicalCompressionReader->SetDebugLevel(this->debugLevel_);
    return readFileLocal_(pathToFile, this->topologicalCompressionReader);
  } else if(pathToFile.substr(pathToFile.size() - 4) == ".tif"
            || pathToFile.substr(pathToFile.size() - 5) == ".tiff") {
    return readFileLocal_(pathToFile, this->tiffReader);
  } else {
    // Check if dataset is XML encoded
    std::ifstream is(pathToFile.data());
    char prefix[10] = "";
    is.get(prefix, 10);
    bool isXML = std::string(prefix).compare("<VTKFile ") == 0
                 || std::string(prefix).compare("<?xml ver") == 0;

    if(isXML)
      // If isXML use vtkXMLGenericDataObjectReader
      return readFileLocal_(pathToFile, this->xmlGenericDataObjectReader);
    else
      // Otherwise use vtkGenericDataObjectReader
      return readFileLocal_(pathToFile, this->genericDataObjectReader);
  }

  return nullptr;
}

int ttkCinemaProductReader::addFieldDataRecursively(vtkDataObject *object,
                                                    vtkFieldData *fd) {
  auto objectAsMB = vtkMultiBlockDataSet::SafeDownCast(object);
  if(objectAsMB) {
    for(size_t i = 0, j = objectAsMB->GetNumberOfBlocks(); i < j; i++)
      addFieldDataRecursively(objectAsMB->GetBlock(i), fd);
  } else {
    auto ofd = object->GetFieldData();
    for(size_t i = 0, j = fd->GetNumberOfArrays(); i < j; i++)
      ofd->AddArray(fd->GetAbstractArray(i));
  }

  return 1;
}

int ttkCinemaProductReader::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  ttk::Timer timer;

  // get input
  auto inputTable = vtkTable::GetData(inputVector[0]);

  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  size_t n = inputTable->GetNumberOfRows();
  size_t m = inputTable->GetNumberOfColumns();
  // Determine number of files
  this->printMsg(
    {{"#Files", std::to_string(n)}, {"FILE Column", this->FilepathColumnName}});
  this->printMsg(ttk::debug::Separator::L1);

  // Read Data
  {
    // Get FILE column
    auto paths = inputTable->GetColumnByName(this->FilepathColumnName.data());
    if(!paths) {
      this->printErr("Table does not have column '" + this->FilepathColumnName
                     + "'.");
      return 0;
    }

    // For each row
    for(size_t i = 0; i < n; i++) {

      // initialize timer for individual file
      ttk::Timer fileTimer;

      // get filepath
      auto path = paths->GetVariantValue(i).ToString();
      auto file = path.substr(path.find_last_of("/") + 1);

      // print progress
      this->printMsg("Reading (" + std::to_string(i + 1) + "/"
                       + std::to_string(n) + "): \"" + file + "\"",
                     0, ttk::debug::LineMode::REPLACE);

      // read local file
      {
        std::ifstream infile(path.data());
        bool exists = infile.good();
        if(!exists) {
          this->printErr("File does not exist.");
          return 0;
        }

        auto readerOutput = this->readFileLocal(path);
        if(!readerOutput) {
          this->printErr("Unable to read file.");
          return 0;
        }

        outputMB->SetBlock(i, readerOutput);
      }

      // augment data products with row data
      {
        auto block = outputMB->GetBlock(i);
        auto fieldData = block->GetFieldData();
        for(size_t j = 0; j < m; j++) {
          auto columnName = inputTable->GetColumnName(j);

          // always write FILE column
          if(!fieldData->HasArray(columnName)
             || columnName == this->FilepathColumnName) {
            if(inputTable->GetColumn(j)->IsNumeric()) {
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

        if(this->AddFieldDataRecursively)
          this->addFieldDataRecursively(block, fieldData);

        this->printMsg("Reading (" + std::to_string(i + 1) + "/"
                         + std::to_string(n) + "): \"" + file + "\"",
                       1, fileTimer.getElapsedTime());
      }
    }
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#products: " + std::to_string(n) + ")", 1,
                 timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
