/// \ingroup base
/// \class ttk::WebSocketIO
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 13.04.2021
///
/// This module provides functions for bi-directional communication between two
/// websockets, where the heavy lifting is done by the websocketpp library. The
/// primary use case for this module is to send data that was computed in a c++
/// environment to a browser client, which can then freely process/render the
/// data via JavaScript and HTML. To this end, this module needs to run a
/// websocket server, to which a client can connect to via the ttkWebSocketIO.js
/// library. The server can also receive data from the client, which is then fed
/// into the processEvent function. In order to add custom processing of events,
/// one needs to inherit from the WebSocketIO class and override the
/// processEvent method. The vtk wrapper ttkWebSocketIO provides an example on
/// how this abstract base module can be specialized to the application of
/// sending/receiving vtkDataObjects. Specifically, the ttkWebSocketIO filter
/// will send its input vtkDataObject as a serialized JSON object to all
/// connected clients every time the filter is called with a new input. When the
/// server receives a serialized JSON object form the client, then the filter
/// will instantiate a vtkDataObject and pass it as the filter output.

#pragma once

#include <Debug.h>
#include <list>
#include <set>

#include <iostream>

#include <functional>

#if TTK_ENABLE_WEBSOCKETPP

// fix undefined reference to boost::system link error
#define BOOST_ERROR_CODE_HEADER_ONLY

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

using websocketpp::connection_hdl;
using websocketpp::lib::thread;

using con_list = std::set<connection_hdl, std::owner_less<connection_hdl>>;
using WSServer = websocketpp::server<websocketpp::config::asio>;

#endif

namespace ttk {

  class WebSocketIO : virtual public Debug {
  public:
    struct Message {
      size_t binaryPayloadSize{0}; // Used for Binary Data
      const void *binaryPayload{nullptr}; // Used for Binary Data
      std::string stringPayload; // Used for String Data
      Message(const std::string &msg) : stringPayload(msg) {
      }

      Message(const size_t &sizeInBytes, const void *data)
        : binaryPayloadSize(sizeInBytes), binaryPayload(data) {
      }
      Message(const Message &msg) {
        binaryPayloadSize = msg.binaryPayloadSize;
        binaryPayload = msg.binaryPayload;
        stringPayload = msg.stringPayload;
      }
    };

    WebSocketIO();
    ~WebSocketIO();

    int isListening();

    int virtual processEvent(const std::string &eventName,
                             const std::string &eventData = "");

    int startServer(int PortNumber);
    int stopServer();

    int getPortNumber();

    int sendString(const std::string &msg) const;
    int sendBinary(const size_t &sizeInBytes, const void *data) const;
    int sendMessage(const Message &msg) const;

    int queueMessage(const std::string &msg);
    int queueMessage(const size_t &sizeInBytes, const void *data);
    int queueMessage(const Message &msg);
    int sendNextQueuedMessage();
    int processMessageQueue();
    int clearMessageQueue();

  private:
#if TTK_ENABLE_WEBSOCKETPP
    mutable WSServer server;

    ttk::Timer msgTimer;
    size_t nMessages;

    std::thread *serverThread = nullptr;
    con_list connections;
    std::mutex mutex;
    websocketpp::lib::error_code ec;

    // keep the state of the object sending process
    bool serverThreadRunning = false;
    int portNumber = 0;

    int on_open(const websocketpp::connection_hdl &hdl);
    int on_close(const websocketpp::connection_hdl &hdl);
    int on_message(const websocketpp::connection_hdl &hdl,
                   const WSServer::message_ptr &msg);

    std::list<Message> messageQueue;
#endif
  };
} // namespace ttk
