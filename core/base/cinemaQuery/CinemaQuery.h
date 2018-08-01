/// \ingroup base
/// \class ttk::CinemaQuery
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %cinemaQuery processing package.
///
/// %CinemaQuery is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkCinemaQuery.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;
#include <iostream>
#include <vector>

using namespace std;

namespace ttk{

  class CinemaQuery : public Debug{

    public:

      CinemaQuery(){};
      ~CinemaQuery(){};

      template <class dataType>
        vector<vector<string>> execute(const string &serverAddress, const string &sqlQuery) const;
      //This generates a valid url from a string
      inline string generateURL(string input) const{
        string resultString;
        //iterator over blanks
        for(string::iterator it = input.begin(); it != input.end(); ++it) {
          if(*it == ' ' || *it == '\n') {
            resultString += "%20";
          } else if(*it == '%'){
            resultString += "%25";
          } else {
            resultString += *it;
          }
        }

        return resultString;
      }

      inline string queryDataBase(string serverAddress, string sqlQuery) const{
        //change this to io_service for boost version < 1.66
        boost::asio::io_context context;
        boost::asio::ip::tcp::socket socket(context);


        string ip = serverAddress.substr(0,serverAddress.find_last_of(":"));
        string port = serverAddress.substr(serverAddress.find_last_of(":")+1);

        {
            stringstream msg;
            msg << "[CinemaQuery] Server: "
                << serverAddress << endl
                << "[CinemaQuery]  Query:\n" << sqlQuery << endl;
            dMsg(cout, msg.str(), timeMsg);
        }

        boost::asio::ip::tcp::endpoint endpoint(
            boost::asio::ip::address::from_string( ip ),
            stoi(port)
        );

        boost::system::error_code ec;
        socket.connect(endpoint, ec);
        if (ec) {
          cout << "Connection Error: ";
          cout << ec.message() << endl;
          //TODO error handling in paraview
        }

        // create requestQuery
        string validQuery = ttk::CinemaQuery::generateURL(sqlQuery);

        //Concatenation of requestString
        stringstream requestStringStream;
        // requestStringStream << "GET " << validQuery << " HTTP/1.1\r\n\r\n";
        requestStringStream << "GET /" << validQuery << " HTTP/1.1\r\n";
        requestStringStream << "Connection: close\r\n";
        requestStringStream << "\r\n";
        string requestString = requestStringStream.str();

        // send requestQuery
        socket.send(boost::asio::buffer(requestString));

        // Handle the serverResponse
        // Response Buffer
        boost::asio::streambuf response;
        boost::asio::read_until(socket, response, "\r\n");

        // Check that response is Valid.
        istream response_stream(&response);

        // response_stream = HTTP/1.1 200 OK
        string http_version;
        response_stream >> http_version;
        unsigned int status_code;
        response_stream >> status_code;
        string status_message;
        getline(response_stream, status_message);
        if (!response_stream || http_version.substr(0, 5) != "HTTP/") {
          cout << "Invalid response\n";
          return "ERROR";
          //TODO error handling in paraview
        }
        if (status_code != 200) {
          cout << "Response returned with status code " << status_code << "\n";
          return "ERROR";
          //TODO error handling in paraview
        }

        boost::system::error_code error;
        while (boost::asio::read(socket, response, boost::asio::transfer_at_least(1), error)){
        //   cout << "Query Successful" << endl;
        }

        if (error != boost::asio::error::eof) {
          throw boost::system::system_error(error);
        }

        //Convert streambuffer to string
        //Creates String with elements from first element to last element of buffer
        boost::asio::streambuf::const_buffers_type bufs = response.data();
        string resultStr(boost::asio::buffers_begin(bufs),
                        boost::asio::buffers_begin(bufs) + response.size());

        // Todo check what happens, if json is multirow
        istream is(&response);
        string line;
        string testJ;
        // parse the result and return it
        while (getline(is, line)) {
            if (line.compare(0, 2, "[{") == 0){
              if(line.compare("[{ERROR}]") == 0){
                cout << "Error in parse" << line <<endl;
                return "ERROR";
              } else {
                testJ = line;
              }
            }
        }
        return testJ;
      };

      inline vector<string> getVecOfAttributes(istringstream &input) const{
        // Root Node
        pt::ptree root;

        // Fill ptree with the input
        pt::read_json(input, root);
        pt::ptree children = root.get_child("");

        //result vector
        vector<string> vec;
        for (const auto& kv : children) {
          pt::ptree temp = kv.second;
          string test = temp.get<string>("COLUMN_NAME");
          vec.push_back(test);
        }
        return vec;
      };

      inline vector<vector<string>> parseQueryResults(istringstream &input) const{

        vector<vector<string>> matrix;
        if(input.str().length()==0)
            return matrix;

        // Fill ptree with the input
        pt::ptree root;
        pt::read_json(*&input, root);

        pt::ptree children = root.get_child("");

        // create result matrix
        int no_of_rows = children.size();
        int no_of_cols=0;
        vector<string> colNames;
        for (const auto& kv : children) {
            pt::ptree temp = kv.second;
            no_of_cols = temp.size();
            for(const auto& c: temp)
                colNames.push_back(c.first);

            break;
        }

        string initial_value = string("NULL");

        matrix.resize(no_of_rows+1, vector<string>(no_of_cols, initial_value));
        for(int i=0; i < no_of_cols; i++) {
            matrix[0][i] = colNames[i];
        }

        // Fill matrix and return it
        int rowCount = 1;
        for (const auto& kv : children) {
          pt::ptree temp = kv.second;

          for(int i=0; i < no_of_cols; i++) {
              // Crashes here when you use select statement
              //TODO Filter selected cols and adjust the matrix
              matrix[rowCount][i] =  temp.get<string>(colNames[i]);
          }
          rowCount++;
        }

        return matrix;
      }

    protected:
  };
}


// template functions
template <class dataType> vector<vector<string>> ttk::CinemaQuery::execute(
  const string &serverAddress, const string &sqlQuery) const{

    Timer t;

    string resultJson = ttk::CinemaQuery::queryDataBase(serverAddress, sqlQuery);

    string er_str= "ERROR";
    if(resultJson.compare(er_str) == 0){
        cout << "Error" << endl;
        string initial_value = string("ERROR");
        vector<vector<string>> matrix;
        matrix.resize(1, vector<string>(1, initial_value));
        return matrix;
    }

    istringstream streamQueryResult(resultJson);

    vector<vector<string>> resultMartix = parseQueryResults(streamQueryResult);

    int n = resultMartix.size()>1 ? resultMartix.size()-1 : 0;
    int m = n>0 ? resultMartix[0].size() : 0;

    {
        stringstream msg;
        msg << "[CinemaQuery] ("<< n << ","<<m<<") tuples fetched in "
            << t.getElapsedTime() << " s."
            << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return resultMartix;
}


