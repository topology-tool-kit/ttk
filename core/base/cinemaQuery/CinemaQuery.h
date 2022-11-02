/// \ingroup base
/// \class ttk::CinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %cinemaQuery processing package.
///
/// %CinemaQuery is a TTK processing package that generates a temporary SQLite3
/// Database to perform a SQL query which is returned as a CSV String
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">Contour
///   Tree Alignment example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/nestedTrackingFromOverlap/">Nested
///   Tracking from Overlap example</a> \n

#pragma once

// base code includes
#include <Debug.h>
#include <string>
#include <vector>

namespace ttk {
  class CinemaQuery : virtual public Debug {
  public:
    CinemaQuery();
    ~CinemaQuery() override;

    /** Creates a temporary database based on a SQL table definition and
     *  and table content to subsequentually return a query result.
     */
    int execute(const std::vector<std::string> &sqlTableDefinitions,
                const std::vector<std::string> &sqlInsertStatements,
                const std::string &sqlQuery,
                std::stringstream &resultCSV,
                int &csvNColumns,
                int &csvNRows) const;
  };
} // namespace ttk
