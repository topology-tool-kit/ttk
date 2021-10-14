#include <ttkCinemaReader.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaReader);

ttkCinemaReader::ttkCinemaReader() {
  this->setDebugMsgPrefix("CinemaReader");

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaReader::~ttkCinemaReader() {
}

int ttkCinemaReader::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  else
    return 0;
  return 1;
}

int ttkCinemaReader::validateDatabasePath() {
  if(this->DatabasePath.length() < 4
     || this->DatabasePath.substr(this->DatabasePath.length() - 4, 4)
            .compare(".cdb")
          != 0) {
    this->printErr("Database path has to end with '.cdb'.");
    return 0;
  }

  return 1;
}

int ttkCinemaReader::RequestData(vtkInformation *ttkNotUsed(request),
                                 vtkInformationVector **ttkNotUsed(inputVector),
                                 vtkInformationVector *outputVector) {
  ttk::Timer timer;

  // print input
  this->printMsg({
    {"Database", this->GetDatabasePath()},
    {"FILE Columns", this->GetFilePathColumnNames()},
  });
  this->printMsg(ttk::debug::Separator::L1);

  if(!this->validateDatabasePath())
    return 0;

  this->printMsg("Reading CSV file", 0, ttk::debug::LineMode::REPLACE);

  // get output
  auto outTable = vtkTable::GetData(outputVector);

  // read CSV file which is in Spec D format
  {
    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetFileName((this->GetDatabasePath() + "/data.csv").data());
    reader->DetectNumericColumnsOn();
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(",");
    reader->Update();

    if(reader->GetLastError().compare("") != 0)
      return 0;

    // copy information to output
    outTable->ShallowCopy(reader->GetOutput());

    // prepend database path to FILE column
    std::vector<std::string> filePathColumnNames;
    ttkUtils::stringListToVector(
      this->GetFilePathColumnNames(), filePathColumnNames);

    for(size_t j = 0; j < filePathColumnNames.size(); j++) {
      auto filepathColumn = vtkStringArray::SafeDownCast(
        outTable->GetColumnByName(filePathColumnNames[j].data()));
      if(filepathColumn) {
        size_t n = filepathColumn->GetNumberOfValues();
        for(size_t i = 0; i < n; i++)
          filepathColumn->SetValue(
            i, this->GetDatabasePath() + "/" + filepathColumn->GetValue(i));
      } else {
        this->printErr("Input table does not have column '"
                       + filePathColumnNames[j]
                       + "' or is not of type 'vtkStringArray'.");
        return 0;
      }
    }

    this->printMsg("Reading CSV file", 1, timer.getElapsedTime());
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg(
    "Complete (#rows: " + std::to_string(outTable->GetNumberOfRows()) + ")", 1,
    timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
