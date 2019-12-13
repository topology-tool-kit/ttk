/// \ingroup base
/// \class ttk::CinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %cinemaQuery processing package.
///
/// %CinemaQuery is a TTK processing package that generates a temporary SQLite3
/// Database to perform a SQL query which is returned as a CSV String

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk {
  class CinemaQuery : public Debug {
  public:
    CinemaQuery();
    ~CinemaQuery();

    /** @brief Create a temporary database to execute the SQL query
     *
     * Create a in-memory SQLite database, create tables in it and
     * fill them with input data. Execute the provided SQL query onto
     * those tables.
     *
     * @param[in] sqlTablesDefinitionAndRows Vector of string
     * pairs. The first string is a CREATE TABLE query, the second is
     * an INSERT INTO to fill data in the corresponding table.
     * @param[in] sqlQuery SQL query that will be executed on the
     * filled tables
     * @param[out] resultCSV SQL query output in a CSV format
     *
     * @return 1 in case of success
     */
    int execute(
      const std::vector<std::pair<string, string>> &sqlTablesDefinitionAndRows,
      const string &sqlQuery,
      string &resultCSV) const;
  };
} // namespace ttk