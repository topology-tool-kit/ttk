/// \ingroup base
/// \class ttk::WebSocketIO
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2019
///
/// TODO

#pragma once

#include <Debug.h>
#include <list>
#include <set>

#include <iostream>

#include <functional>

#if TTK_ENABLE_WEBSOCKETPP

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

typedef websocketpp::server<websocketpp::config::asio> server;
using websocketpp::connection_hdl;
using websocketpp::lib::bind;
using websocketpp::lib::thread;

typedef std::set<connection_hdl, std::owner_less<connection_hdl>> con_list;
typedef websocketpp::server<websocketpp::config::asio> WSServer;

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

    int on_open(websocketpp::connection_hdl hdl);
    int on_close(websocketpp::connection_hdl hdl);
    int on_message(websocketpp::connection_hdl hdl, server::message_ptr msg);

    int packageIndex = 0;
    std::list<Message> messageQueue;
#endif
  };
} // namespace ttk