#include <ttkCinemaQuery.h>

#include <vtkVersion.h>

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkCinemaQuery)

  int ttkCinemaQuery::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  // Print Status
  {
    stringstream msg;
    msg << "==================================================================="
           "============="
        << endl;
    msg << "[ttkCinemaQuery] RequestData" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  Timer t;
  Memory m;

  // -------------------------------------------------------------------------
  // Get Input Tables
  // -------------------------------------------------------------------------
  auto nTables = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkTable *> inTables(nTables);
  for(int i = 0; i < nTables; ++i) {
    inTables[i] = vtkTable::GetData(inputVector[0], i);
  }

  // -------------------------------------------------------------------------
  // Get Output Table
  // -------------------------------------------------------------------------
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto outTable
    = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // -------------------------------------------------------------------------
  // Convert Input Table to SQL Table
  // -------------------------------------------------------------------------
  std::vector<std::pair<string, string>> sqlTablesDefinitionAndRows(nTables);

  for(int k = 0; k < nTables; ++k) {
    const auto &inTable = inTables[k];
    string sqlTableDefinition, sqlTableRows;
    int nc = inTable->GetNumberOfColumns();
    int nr = inTable->GetNumberOfRows();

    sqlTableDefinition = "CREATE TABLE InputTable" + std::to_string(k) + " (";
    vector<bool> isNumeric(nc);
    for(int i = 0; i < nc; i++) {
      auto c = inTable->GetColumn(i);
      isNumeric[i] = c->IsNumeric();
      sqlTableDefinition += (i > 0 ? "," : "") + string(c->GetName()) + " "
                            + (isNumeric[i] ? "REAL" : "TEXT");
    }
    sqlTableDefinition += ")";

    sqlTableRows = "INSERT INTO InputTable" + std::to_string(k) + " VALUES ";
    for(int j = 0; j < nr; j++) {
      if(j > 0)
        sqlTableRows += ",";
      sqlTableRows += "(";
      for(int i = 0; i < nc; i++) {
        if(i > 0)
          sqlTableRows += ",";
        if(isNumeric[i])
          sqlTableRows += inTable->GetValue(j, i).ToString();
        else
          sqlTableRows += "'" + inTable->GetValue(j, i).ToString() + "'";
      }
      sqlTableRows += ")";
    }
    sqlTablesDefinitionAndRows[k] = {sqlTableDefinition, sqlTableRows};
  }

  // -------------------------------------------------------------------------
  // Backward compatibility: replace "InputTable" with "InputTable0"
  // in query string
  // -------------------------------------------------------------------------
  const std::string needle{"InputTable"};
  auto index = this->QueryString.find(needle);
  while(index != std::string::npos) {
    // move index to the end of needle
    index += needle.size();
    if(this->QueryString.size() == index) {
      // add "0" to the end of the query
      this->QueryString.append("0");
    } else if(!std::isdigit(this->QueryString[index])) {
      // replace "e" in "InputTable" with "e0"
      this->QueryString.replace(index - 1, 1, "e0");
    }
    index = this->QueryString.find(needle, index);
  }

  // -------------------------------------------------------------------------
  // Replace Variables in QueryString (e.g. {time[2]})
  // -------------------------------------------------------------------------
  string finalQueryString = this->QueryString;

  // TODO: test replace variables with several input tables
  for(const auto &inTable : inTables) {
    vtkFieldData *fieldData = inTable->GetFieldData();

    vector<string> temp; // vector of substrings of SQL string separated by '{'
    boost::split(temp, finalQueryString, boost::is_any_of("{"));

    for(size_t i = 1; i < temp.size(); i++) {
      vector<string> temp2; // vector of substrings of temp separated by '}'
      boost::split(temp2, temp[i], boost::is_any_of("}"));

      string varToken = temp2[0]; // e.g. time[2]
      string varName = varToken; // time
      size_t varIndex = 0; // 2

      // If varIndex specified in string then update varName and varIndex
      size_t indexDelimiter0 = varToken.find("[");
      size_t indexDelimiter1 = varToken.find("]");
      if(indexDelimiter0 != string::npos && indexDelimiter1 != string::npos) {
        varName = varToken.substr(0, indexDelimiter0);
        varIndex = stoi(varToken.substr(
          indexDelimiter0 + 1, indexDelimiter1 - indexDelimiter0 - 1));
      }

      // Search in candidates
      auto column = fieldData->GetAbstractArray(varName.data());
      if(column != nullptr
         && varIndex < ((size_t)column->GetNumberOfTuples())) {
        // Replace Variable
        boost::replace_all(
          finalQueryString, "{" + varToken + "}",
          "\"" + column->GetVariantValue(varIndex).ToString() + "\"");
      } else {
        // Print Error
        stringstream msg;
        msg << "[ttkCinemaQuery] ERROR: Variable {" << varToken
            << "} not found in field data." << endl;
        dMsg(cout, msg.str(), fatalMsg);
      }
    }
  }

  // -------------------------------------------------------------------------
  // Compute Query Result
  // -------------------------------------------------------------------------
  string result = "";
  {
    int status = cinemaQuery.execute(
      sqlTablesDefinitionAndRows, finalQueryString, result);
    if(status != 1)
      return 0;
    if(result.compare("") == 0) {
      stringstream msg;
      msg << "[ttkCinemaQuery] ERROR: Empty result (VTK does not support empty "
             "tables)."
          << endl;
      dMsg(cout, msg.str(), fatalMsg);
    }
  }

  // -------------------------------------------------------------------------
  // Process Result
  // -------------------------------------------------------------------------
  {
#if VTK_MAJOR_VERSION <= 7
    stringstream msg;
    msg << "[ttkCinemaQuery] ERROR: VTK version too old." << endl
        << "[ttkCinemaQuery]        This filter requires vtkDelimitedTextReader"
        << endl
        << "[ttkCinemaQuery]        of version 7.0 or higher." << endl;
    dMsg(cout, msg.str(), fatalMsg);
    return 0;
#else
    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetReadFromInputString(true);
    reader->SetInputString(result);
    reader->DetectNumericColumnsOn();
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(",");
    reader->Update();

    outTable->ShallowCopy(reader->GetOutput());

    // only copy first available Field Data
    for(const auto &inTable : inTables) {
      const auto fd = inTable->GetFieldData();
      // pass database name
      if(fd != nullptr && fd->HasArray("DatabasePath")) {
        outTable->GetFieldData()->ShallowCopy(fd);
        break;
      }
    }

#endif
  }

  // Output Performance
  {
    stringstream msg;
    msg << "[ttkCinemaQuery] "
           "---------------------------------------------------------------"
        << endl;
    msg << "[ttkCinemaQuery]   time: " << t.getElapsedTime() << " s." << endl;
    msg << "[ttkCinemaQuery] memory: " << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
