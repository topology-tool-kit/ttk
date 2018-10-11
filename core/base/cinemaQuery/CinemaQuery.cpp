#include <CinemaQuery.h>

#if TTK_ENABLE_SQLITE3

    #include <sqlite3.h>

    // Process a row of a query result and append it to a string
    static int processRow(void *data, int argc, char **argv, char **azColName){
        int i;

        // Get output string as reference
        string& result = *((string*)data);

        // If string is empty then add a first row that records column names
        if(result==""){
            for(int i = 0; i<argc; i++)
                result+= (i>0?",":"") + string(azColName[i]);
            result+="\n";
        }

        // Append row content to string
        for(i = 0; i<argc; i++)
            result+= (i>0?",":"") + string(argv[i]);
        result+="\n";

        return 0;
    }
#endif

ttk::CinemaQuery::CinemaQuery(){}
ttk::CinemaQuery::~CinemaQuery(){}

string ttk::CinemaQuery::execute(
    const string& sqlTableDefinition,
    const string& sqlTableRows,
    const string& sqlQuery
) const{

    string result="";

    #if TTK_ENABLE_SQLITE3
        // SQLite Variables
        sqlite3* db;
        char* zErrMsg = 0;
        int rc;

        // Create Temporary Database
        {
            Timer t;

            // Print status
            {
                stringstream msg;
                msg << "[ttkCinemaQuery] Create Temporary Database"<<endl;
                dMsg(cout, msg.str(), timeMsg);
            }

            // Initialize DB in memory
            rc = sqlite3_open(":memory:", &db);
            if(rc){
                fprintf(stderr, "[ttkCinemaQuery] Unable to create database: %s\n", sqlite3_errmsg(db));
                return result;
            }

            // Create table
            rc = sqlite3_exec(db, sqlTableDefinition.data(), nullptr, 0, &zErrMsg);
            if( rc != SQLITE_OK ){
                fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
                sqlite3_free(zErrMsg);
                sqlite3_close(db);
                return result;
            }

            // Fill table
            rc = sqlite3_exec(db, sqlTableRows.data(), nullptr, 0, &zErrMsg);
            if( rc != SQLITE_OK ){
                fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
                sqlite3_free(zErrMsg);
                sqlite3_close(db);
                return result;
            }

            // Print status
            {
                stringstream msg;
                msg << "[ttkCinemaQuery] Created in "
                    << t.getElapsedTime() << " s."<< endl;
                dMsg(cout, msg.str(), timeMsg);
            }
        }

        // Run SQL on temporary database
        {
            Timer t;

            // Perform query
            rc = sqlite3_exec(db, sqlQuery.data(), processRow, (void*)(&result), &zErrMsg);
            if( rc != SQLITE_OK ){
                fprintf(stderr, "[ttkCinemaQuery] SQL error: %s\n", zErrMsg);
                sqlite3_free(zErrMsg);
                sqlite3_close(db);
                return result;
            }

            // Delete DB
            sqlite3_close(db);

            // Print status
            {
                stringstream msg;
                msg << "[ttkCinemaQuery] Query processed in "
                    << t.getElapsedTime() << " s." << endl;
                dMsg(cout, msg.str(), timeMsg);
            }
        }

    #else
        // Print Error
        {
            stringstream msg;
            msg << "[ttkCinemaQuery] ERROR: TTK is build without Sqlite3 support. Query is not executed."<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    #endif

    return result;
}
