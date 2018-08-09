/// \ingroup base
/// \class ttk::CinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %cinemaQuery processing package.
///
/// %CinemaQuery is a TTK processing package that generates a temporary SQLite3 Database to perform a SQL query which is returned as a CSV String
///
/// \sa ttkCinemaQuery.cpp %for a usage example.

#pragma once

// base code includes
#include <Wrapper.h>
#include <vector>
#include <sqlite3.h>

using namespace std;

namespace ttk{

    class CinemaQuery : public Debug{

        public:

            CinemaQuery(){};
            ~CinemaQuery(){};

            template <class dataType> string execute(
                const string& sqlTableDefinition,
                const string& sqlTableRows,
                const string& sqlQuery
            ) const;
    };
}

static int processRow(void *data, int argc, char **argv, char **azColName){
    int i;
    string& result = *((string*)data);

    if(result==""){
        for(i = 0; i<argc; i++)
            result+= (i>0?",":"") + string(azColName[i]);
        result+="\n";
    }

    for(i = 0; i<argc; i++)
        result+= (i>0?",":"") + string(argv[i]);

    result+="\n";

    return 0;
}

// template functions
template <class dataType> string ttk::CinemaQuery::execute(
    const string& sqlTableDefinition,
    const string& sqlTableRows,
    const string& sqlQuery
) const{

    string result;

    sqlite3* db;
    char* zErrMsg = 0;
    int rc;

    {
        Timer t;

        {
            stringstream msg;
            msg << "[ttkCinemaQuery] Create Temporary Database"<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }

        rc = sqlite3_open(":memory:", &db);

        if(rc){
            fprintf(stderr, "[ttkCinemaQuery] Unable to create database: %s\n", sqlite3_errmsg(db));
            return result;
        }

        rc = sqlite3_exec(db, sqlTableDefinition.data(), nullptr, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
            sqlite3_close(db);
            return result;
        }

        rc = sqlite3_exec(db, sqlTableRows.data(), nullptr, 0, &zErrMsg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
            sqlite3_close(db);
            return result;
        }

        {
            stringstream msg;
            msg << "[ttkCinemaQuery] Created in "
                << t.getElapsedTime() << " s."<< endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    {
        Timer t;

        rc = sqlite3_exec(db, sqlQuery.data(), processRow, (void*)(&result), &zErrMsg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
            sqlite3_close(db);
            return result;
        }
        sqlite3_close(db);

        {
            stringstream msg;
            msg << "[ttkCinemaQuery] Query processed in "
                << t.getElapsedTime() << " s." << endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    return result;
}


