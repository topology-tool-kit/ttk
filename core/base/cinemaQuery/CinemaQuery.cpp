#include <CinemaQuery.h>
#include <iostream>

#if TTK_ENABLE_SQLITE3
#include <sqlite3.h>
#endif

ttk::CinemaQuery::CinemaQuery() {
  this->setDebugMsgPrefix("CinemaQuery");
}
ttk::CinemaQuery::~CinemaQuery() {
}

int ttk::CinemaQuery::execute(
  const std::vector<std::string> &sqlTableDefinitions,
  const std::vector<std::string> &sqlInsertStatements,
  const std::string &sqlQuery,
  std::stringstream &resultCSV,
  int &csvNColumns,
  int &csvNRows) const {

#if TTK_ENABLE_SQLITE3
  // print input
  {
    std::vector<std::string> sqlLines;
    {
      std::stringstream ss(sqlQuery);
      std::string line;
      while(std::getline(ss, line))
        sqlLines.push_back(line);
    }
    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg(sqlLines);
    this->printMsg(ttk::debug::Separator::L1);
  }

  // SQLite Variables
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;

  // Create Temporary Database
  {
    Timer timer;
    this->printMsg(
      "Creating inmemory database", 0, ttk::debug::LineMode::REPLACE);

    // Initialize DB in memory
    rc = sqlite3_open(":memory:", &db);
    if(rc != SQLITE_OK) {
      this->printErr("Creating database: " + std::string{sqlite3_errmsg(db)});
      return 0;
    }

    // Create table
    for(auto &sqlTableDefinition : sqlTableDefinitions) {
      rc = sqlite3_exec(db, sqlTableDefinition.data(), nullptr, 0, &zErrMsg);
      if(rc != SQLITE_OK) {
        this->printErr("Create table: " + std::string{zErrMsg});

        sqlite3_free(zErrMsg);
        sqlite3_close(db);

        return 0;
      }
    }

    // Fill table
    for(auto &sqlInsertStatement : sqlInsertStatements) {
      rc = sqlite3_exec(db, sqlInsertStatement.data(), nullptr, 0, &zErrMsg);
      if(rc != SQLITE_OK) {
        this->printErr("Insert values: " + std::string{zErrMsg});

        sqlite3_free(zErrMsg);
        sqlite3_close(db);
        return 0;
      }
    }

    this->printMsg("Creating inmemory database", 1, timer.getElapsedTime());
  }

  // Run SQL statement on temporary database
  {
    this->printMsg("Querying database", 0, ttk::debug::LineMode::REPLACE);
    Timer timer;

    sqlite3_stmt *sqlStatement;

    if(sqlite3_prepare_v2(db, sqlQuery.data(), -1, &sqlStatement, NULL)
       != SQLITE_OK) {
      this->printErr("Query: " + std::string{sqlite3_errmsg(db)});

      sqlite3_close(db);
      return 0;
    }
    csvNColumns = sqlite3_column_count(sqlStatement);

    // Get Header
    {
      if(csvNColumns < 1) {
        this->printErr("Query result has no columns.");

        sqlite3_close(db);
        return 0;
      }

      resultCSV << sqlite3_column_name(sqlStatement, 0);
      for(int i = 1; i < csvNColumns; i++)
        resultCSV << "," << sqlite3_column_name(sqlStatement, i);

      resultCSV << "\n";
    }

    // Get Content
    {
      while((rc = sqlite3_step(sqlStatement)) == SQLITE_ROW) {
        csvNRows++;

        resultCSV << sqlite3_column_text(sqlStatement, 0);
        for(int i = 1; i < csvNColumns; i++)
          resultCSV << "," << sqlite3_column_text(sqlStatement, i);
        resultCSV << "\n";
      }

      if(rc != SQLITE_DONE) {
        this->printErr("Fetching result: " + std::string{sqlite3_errmsg(db)});

        sqlite3_close(db);
        return 0;
      } else {
        this->printMsg("Querying database", 1, timer.getElapsedTime());
      }
    }

    sqlite3_finalize(sqlStatement);
  }

  // Close database
  {
    Timer timer;
    this->printMsg("Closing database", 0, ttk::debug::LineMode::REPLACE);

    // Delete DB
    rc = sqlite3_close(db);

    // Print status
    if(rc != SQLITE_OK) {
      this->printErr("Closing database:" + std::string{sqlite3_errmsg(db)});
      return 0;
    } else {
      this->printMsg("Closing database", 1, timer.getElapsedTime());
    }
  }

  return 1;

#else
  TTK_FORCE_USE(sqlTableDefinitions);
  TTK_FORCE_USE(sqlInsertStatements);
  TTK_FORCE_USE(sqlQuery);
  TTK_FORCE_USE(resultCSV);
  TTK_FORCE_USE(csvNColumns);
  TTK_FORCE_USE(csvNRows);

  this->printErr("This filter requires Sqlite3");
  return 0;
#endif // TTK_ENABLE_SQLITE3
}
