#include <ttkCinemaQuery.h>

#include <vtkInformation.h>

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

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

int ttkCinemaQuery::RequestData(vtkInformation *ttkNotUsed(request),
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

    for(int i = 0; i < nTables; i++) {
      auto inTable = inTables[i];

      size_t nc = inTable->GetNumberOfColumns();
      size_t nr = inTable->GetNumberOfRows();
      std::vector<bool> isNumeric(nc);

      std::vector<size_t> includeColumns(nc);
      std::iota(includeColumns.begin(), includeColumns.end(), 0);

      // exclude input columns whose name matches the regexp
      if(this->ExcludeColumnsWithRegexp) {
        const auto oldSize{includeColumns.size()};
        includeColumns.erase(
          std::remove_if(includeColumns.begin(), includeColumns.end(),
                         [inTable, this](const size_t a) {
                           const auto name{inTable->GetColumnName(a)};
                           return std::regex_match(
                             name, std::regex(this->RegexpString));
                         }),
          includeColumns.end());
        if(includeColumns.size() < oldSize) {
          this->printMsg("Removed "
                         + std::to_string(oldSize - includeColumns.size())
                         + " columns with regexp `" + this->RegexpString + "'");
        }
      }

      // exclude non-scalar columnns
      {
        const auto oldSize{includeColumns.size()};
        includeColumns.erase(
          std::remove_if(includeColumns.begin(), includeColumns.end(),
                         [inTable](const size_t a) {
                           const auto col{inTable->GetColumn(a)};
                           return col->GetNumberOfComponents() != 1;
                         }),
          includeColumns.end());
        if(includeColumns.size() < oldSize) {
          this->printWrn("Removed "
                         + std::to_string(oldSize - includeColumns.size())
                         + " non-scalar columns");
        }
      }

      // exclude columns with a space in their name
      {
        const auto oldSize{includeColumns.size()};
        includeColumns.erase(
          std::remove_if(includeColumns.begin(), includeColumns.end(),
                         [inTable](const size_t a) {
                           const std::string colName{inTable->GetColumnName(a)};
                           return colName.find(' ') != std::string::npos;
                         }),
          includeColumns.end());
        if(includeColumns.size() < oldSize) {
          this->printWrn("Removed "
                         + std::to_string(oldSize - includeColumns.size())
                         + " columns with a space in their name");
        }
      }

      if(includeColumns.empty()) {
        this->printWrn("No columns to process!");
        return 1;
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
            if(isNumeric[k]) {
              const auto var = inTable->GetValue(q, k);
              if(var.IsChar() || var.IsSignedChar()) {
                // convert char/signed char to int to get its value
                // instead of its char representation
                sqlInsertStatement += std::to_string(var.ToInt());
              } else {
                sqlInsertStatement += var.ToString();
              }
            } else {
              sqlInsertStatement
                += "'" + inTable->GetValue(q, k).ToString() + "'";
            }
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
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#rows: " + std::to_string(csvNRows) + ")", 1,
                 timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
