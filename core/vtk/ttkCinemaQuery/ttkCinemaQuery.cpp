#include <ttkCinemaQuery.h>

#include <vtkVersion.h>
#include <vtkInformation.h>

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVersionMacros.h>

#include <ttkUtils.h>

#include <numeric>
#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

vtkStandardNewMacro(ttkCinemaQuery);

ttkCinemaQuery::ttkCinemaQuery() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkCinemaQuery::~ttkCinemaQuery() {
}

int ttkCinemaQuery::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkCinemaQuery::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkCinemaQuery::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  ttk::Timer timer;

  // get output table
  auto outTable = vtkTable::GetData(outputVector);

  // ===========================================================================
  // Convert Input Tables to SQL Tables
  auto nTables = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkTable *> inTables(nTables);
  for(int i = 0; i < nTables; ++i) {
    inTables[i] = vtkTable::GetData(inputVector[0], i);
  }

  auto firstTable = inTables[0];

  std::vector<std::string> sqlTableDefinitions;
  std::vector<std::string> sqlInsertStatements;
  {
    ttk::Timer conversionTimer;
    this->printMsg("Converting input VTK tables to SQL tables", 0,
                   ttk::debug::LineMode::REPLACE);

    for(int i = 0; i < nTables; i++) {
      auto inTable = inTables[i];

      size_t nc = inTable->GetNumberOfColumns();
      size_t nr = inTable->GetNumberOfRows();
      std::vector<bool> isNumeric(nc);

      // select all input columns whose name is NOT matching the regexp
      std::vector<size_t> includeColumns{};
      if(this->ExcludeColumnsWithRegexp) {
        for(size_t j = 0; j < nc; ++j) {
          const auto &name = inTable->GetColumnName(j);
          if(!std::regex_match(name, std::regex(RegexpString))) {
            includeColumns.emplace_back(j);
          }
        }
      } else {
        includeColumns.resize(nc);
        std::iota(includeColumns.begin(), includeColumns.end(), 0);
      }

      // -----------------------------------------------------------------------
      // Table Definition
      std::string sqlTableDefinition
        = "CREATE TABLE InputTable" + std::to_string(i) + " (";
      bool firstCol{true};
      for(const auto j : includeColumns) {
        auto c = inTable->GetColumn(j);
        isNumeric[j] = c->IsNumeric();
        sqlTableDefinition += (firstCol ? "" : ",") + std::string(c->GetName())
                              + " " + (isNumeric[j] ? "REAL" : "TEXT");
        if(firstCol) {
          firstCol = false;
        }
      }
      sqlTableDefinition += ")";
      sqlTableDefinitions.emplace_back(sqlTableDefinition);

      // -----------------------------------------------------------------------
      // Insert Statements
      size_t q = 0;
      while(q < nr) {
        std::string sqlInsertStatement
          = "INSERT INTO InputTable" + std::to_string(i) + " VALUES ";
        for(size_t j = 0; j < 500 && q < nr; j++) {
          if(j > 0)
            sqlInsertStatement += ",";

          sqlInsertStatement += "(";
          firstCol = true;
          for(const auto k : includeColumns) {
            sqlInsertStatement += (firstCol ? "" : ",");
            if(firstCol) {
              firstCol = false;
            }
            if(isNumeric[k])
              sqlInsertStatement += inTable->GetValue(q, k).ToString();
            else
              sqlInsertStatement
                += "'" + inTable->GetValue(q, k).ToString() + "'";
          }
          sqlInsertStatement += ")";
          q++;
        }
        sqlInsertStatements.emplace_back(sqlInsertStatement);
      }
    }

    this->printMsg("Converting input VTK tables to SQL tables", 1,
                   conversionTimer.getElapsedTime());
  }

  // ===========================================================================
  // Replace Variables in Querystd::string (e.g. {time[2]})
  std::string finalQueryString;
  {
    std::string errorMsg;
    if(!ttkUtils::replaceVariables(this->GetSQLStatement(),
                                   firstTable->GetFieldData(), finalQueryString,
                                   errorMsg)) {
      this->printErr(errorMsg);
      return 0;
    }
  }

  // ===========================================================================
  // Compute Query Result
  std::stringstream csvResult;
  int csvNColumns = 0;
  int csvNRows = 0;

  int status
    = this->execute(sqlTableDefinitions, sqlInsertStatements, finalQueryString,
                    csvResult, csvNColumns, csvNRows);

  // ===========================================================================
  // Process Result
  {
#if VTK_MAJOR_VERSION < 7
    this->printErr("This filter requires at least VTK 7.0");
    return 0;
#else
    if(csvNRows < 1 || status != 1) {
      csvResult << "NULL";
      for(int i = 0; i < csvNColumns; i++)
        csvResult << ",NULL";
    }

    ttk::Timer conversionTimer;

    this->printMsg(
      "Converting SQL result to VTK table", 0, ttk::debug::LineMode::REPLACE);

    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetReadFromInputString(true);
    reader->SetInputString(csvResult.str().data());
    reader->DetectNumericColumnsOn();
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(",");
    reader->Update();

    outTable->ShallowCopy(reader->GetOutput());

    auto outFD = outTable->GetFieldData();
    for(const auto inTable : inTables) {
      if(!inTable)
        continue;
      auto inFD = inTable->GetFieldData();

      size_t n = inFD->GetNumberOfArrays();
      for(size_t i = 0; i < n; i++) {
        auto iArray = inFD->GetAbstractArray(i);
        if(!outFD->GetAbstractArray(iArray->GetName())) {
          outFD->AddArray(iArray);
        }
      }
    }

    this->printMsg("Converting SQL result to VTK table", 1,
                   conversionTimer.getElapsedTime());
#endif
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#rows: " + std::to_string(csvNRows) + ")", 1,
                 timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
