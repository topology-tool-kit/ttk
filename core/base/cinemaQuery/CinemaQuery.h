/// \ingroup base
/// \class ttk::CinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %cinemaQuery processing package.
///
/// %CinemaQuery is a TTK processing package that generates a temporary SQLite3 Database to perform a SQL query which is returned as a CSV String

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk{
    class CinemaQuery : public Debug{
        public:
            CinemaQuery();
            ~CinemaQuery();

            /*
            Creates a temporary database based on a SQL table definition and table content
            in order to subsequentually return the result of a query.
            */
            string execute(
                const string& sqlTableDefinition,
                const string& sqlTableRows,
                const string& sqlQuery
            ) const;
    };
}