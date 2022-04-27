#include <WebSocketIO.h>

#include <exception>

ttk::WebSocketIO::WebSocketIO() {
  this->setDebugMsgPrefix("WebSocketIO");

#if TTK_ENABLE_WEBSOCKETPP
  this->server.set_error_channels(websocketpp::log::elevel::none);
  this->server.set_access_channels(websocketpp::log::alevel::none);

  this->server.set_reuse_addr(true);

  // Initialize Asio
  this->server.init_asio();

  // Set the default message handler to the echo handler
  this->server.set_message_handler(bind(&WebSocketIO::on_message, this,
                                        websocketpp::lib::placeholders::_1,
                                        websocketpp::lib::placeholders::_2));
  this->server.set_open_handler(
    bind(&WebSocketIO::on_open, this, websocketpp::lib::placeholders::_1));
  this->server.set_close_handler(
    bind(&WebSocketIO::on_close, this, websocketpp::lib::placeholders::_1));

#else
  this->printErr("WebSocketIO requires websocketpp header only library!");
#endif
}

ttk::WebSocketIO::~WebSocketIO() {
  try {
    this->stopServer();
  } catch(const std::exception &e) {
    this->printErr(e.what());
  }
}

#if TTK_ENABLE_WEBSOCKETPP

int ttk::WebSocketIO::isListening() {
  return this->server.is_listening();
}

int ttk::WebSocketIO::getPortNumber() {
  return this->portNumber;
}

int ttk::WebSocketIO::startServer(int PortNumber) {
  ttk::Timer t;

  this->portNumber = PortNumber;

  this->printMsg("Starting Server at Port: " + std::to_string(this->portNumber),
                 0, 0, ttk::debug::LineMode::REPLACE);

  this->server.reset();
  this->server.listen(this->portNumber);

  // Queues a connection accept operation
  this->server.start_accept();

  // Start the Asio io_service run loop
  this->serverThread = new thread([this]() {
    try {
      {
        std::lock_guard<std::mutex> guard(this->mutex);
        this->serverThreadRunning = true;
      }
      this->server.run();
      {
        std::lock_guard<std::mutex> guard(this->mutex);
        this->serverThreadRunning = false;
      }
    } catch(websocketpp::exception const &e) {
      this->printErr("Unable to start server: " + std::string(e.what()));
    }
  });
  this->serverThread->detach();

  this->printMsg("Starting Server at Port: " + std::to_string(this->portNumber),
                 1, t.getElapsedTime());

  return 1;
}

int ttk::WebSocketIO::stopServer() {
  ttk::Timer t;

  this->printMsg("Stopping Server", 1, 0, ttk::debug::LineMode::REPLACE);

  if(this->server.is_listening()) {

    // Stopping the Websocket listener and closing outstanding connections.
    this->server.stop_listening(this->ec);
    if(this->ec) {
      this->printErr(this->ec.message());
      return 0;
    }

    // Close all existing websocket connections.
    {
      this->printMsg("Closing Connections", 0, 0, ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);
      // initiate closing
      {
        std::lock_guard<std::mutex> guard(this->mutex);
        for(con_list::iterator it = this->connections.begin();
            it != this->connections.end(); ++it) {
          this->server.close(*it, websocketpp::close::status::normal,
                             "Terminating connection ...", this->ec);
          if(this->ec) {
            this->printErr(this->ec.message());
            return 0;
          }
        }
        this->connections.clear();
      }

      // wait until all closed
      size_t nC = 1;
      while(nC > 0) {
        std::lock_guard<std::mutex> guard(this->mutex);
        nC = this->connections.size();
      }

      this->printMsg("Closing Connections", 1, 0, ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);
    }
    // Stop the endpoint.
    {
      this->printMsg("Terminating Server Thread", 0, 0,
                     ttk::debug::LineMode::REPLACE,
                     ttk::debug::Priority::DETAIL);

      this->server.stop();
      // wait until thread terminated
      bool con = true;
      while(con) {
        std::lock_guard<std::mutex> guard(this->mutex);
        con = this->serverThreadRunning;
      }
      delete this->serverThread;
      this->serverThread = nullptr;

      this->printMsg("Terminating Server Thread", 1, 0,
                     ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
    }
  }
  this->printMsg("Stopping Server", 1, t.getElapsedTime());
  return 1;
}

int ttk::WebSocketIO::sendString(const std::string &msg) const {
  if(this->connections.size() > 0) {
    this->server.send(
      *this->connections.begin(), msg, websocketpp::frame::opcode::text);
    return 1;
  }
  return 0;
}

int ttk::WebSocketIO::sendBinary(const size_t &sizeInBytes,
                                 const void *data) const {
  if(this->connections.size() > 0) {
    this->server.send(*this->connections.begin(), data, sizeInBytes,
                      websocketpp::frame::opcode::binary);
    return 1;
  }
  return 0;
}

int ttk::WebSocketIO::sendMessage(const Message &msg) const {
  if(msg.binaryPayload != nullptr)
    return this->sendBinary(msg.binaryPayloadSize, msg.binaryPayload);
  else
    return this->sendString(msg.stringPayload);
}

int ttk::WebSocketIO::queueMessage(const std::string &msg) {
  this->messageQueue.emplace_back(msg);
  return 1;
}
int ttk::WebSocketIO::queueMessage(const size_t &sizeInBytes,
                                   const void *data) {
  this->messageQueue.emplace_back(sizeInBytes, data);
  return 1;
}
int ttk::WebSocketIO::queueMessage(const Message &msg) {
  this->messageQueue.emplace_back(msg);
  return 1;
}
int ttk::WebSocketIO::clearMessageQueue() {
  this->messageQueue = std::list<Message>();
  return 1;
}

int ttk::WebSocketIO::sendNextQueuedMessage() {
  const size_t nRemainingMessages = this->messageQueue.size();
  if(nRemainingMessages < 1) {
    this->printWrn("Empty message queue.");
    return 0;
  }

  if(nRemainingMessages > 1) {
    float progress = ((float)(this->nMessages - nRemainingMessages))
                     / ((float)(this->nMessages - 1));
    this->printMsg(
      "Sending " + std::to_string(this->nMessages) + " queued messages",
      progress, this->msgTimer.getElapsedTime(), ttk::debug::LineMode::REPLACE);
  } else {
    this->printMsg(
      "Sending " + std::to_string(this->nMessages) + " queued messages", 1,
      this->msgTimer.getElapsedTime());
  }

  auto msg = this->messageQueue.front();
  this->messageQueue.pop_front();
  return this->sendMessage(msg);
}

int ttk::WebSocketIO::processMessageQueue() {
  this->nMessages = this->messageQueue.size() + 2;
  msgTimer.reStart();
  this->printMsg(
    "Sending " + std::to_string(this->nMessages) + " queued messages", 0, 0,
    ttk::debug::LineMode::REPLACE);
  if(this->nMessages < 1) {
    this->printWrn("Empty message list.");
    return 0;
  }

  this->messageQueue.emplace_front(
    ttk::WebSocketIO::Message("ttk_WSIO_BeginMessageSequence"));
  this->messageQueue.emplace_back(
    ttk::WebSocketIO::Message("ttk_WSIO_EndMessageSequence"));

  return this->sendNextQueuedMessage();
}

int ttk::WebSocketIO::processEvent(const std::string &eventName,
                                   const std::string &eventData) {
  this->printMsg("processEventBase:" + eventName + " -> " + eventData,
                 ttk::debug::Priority::VERBOSE);

  if(eventName.compare("on_message") == 0) {
    if(eventData.compare("ttk_WSIO_RequestNextMessage") == 0)
      return this->sendNextQueuedMessage();
  }

  return 1;
}

int ttk::WebSocketIO::on_open(const websocketpp::connection_hdl &hdl) {
  std::lock_guard<std::mutex> lock(this->mutex);

  ttk::Timer t;
  this->printMsg(
    "Establishing Connection", 0, 0, ttk::debug::LineMode::REPLACE);

  if(!this->connections.empty()) {
    this->printErr("One client is already connected.");
    this->server.close(hdl, websocketpp::close::status::normal,
                       "Terminating connection ...", ec);
    return 0;
  }

  this->connections.insert(hdl);
  this->printMsg("Establishing Connection", 1, t.getElapsedTime());

  this->processEvent("on_open");

  return 1;
}

int ttk::WebSocketIO::on_close(const websocketpp::connection_hdl &hdl) {
  std::lock_guard<std::mutex> lock(this->mutex);

  ttk::Timer t;
  this->printMsg("Closing Connection", 0, 0, ttk::debug::LineMode::REPLACE);
  this->connections.erase(hdl);
  this->printMsg("Closing Connection", 1, t.getElapsedTime());

  return 1;
}

int ttk::WebSocketIO::on_message(
  const websocketpp::connection_hdl &ttkNotUsed(hdl),
  const WSServer::message_ptr &msg) {

  const auto &eventData = msg->get_payload();
  if(eventData.rfind("ttk_WSIO_", 9) != 0)
    this->printMsg("Custom Message Received", 1, 0);
  this->processEvent("on_message", eventData);
  return 1;
}

#else

int ttk::WebSocketIO::isListening() {
  return 0;
}
int ttk::WebSocketIO::getPortNumber() {
  return -1;
}
int ttk::WebSocketIO::startServer(int ttkNotUsed(PortNumber)) {
  return 0;
}
int ttk::WebSocketIO::stopServer() {
  return 0;
}
int ttk::WebSocketIO::sendString(const std::string &ttkNotUsed(msg)) const {
  return 0;
}
int ttk::WebSocketIO::sendBinary(const size_t &ttkNotUsed(sizeInBytes),
                                 const void *ttkNotUsed(data)) const {
  return 0;
}
int ttk::WebSocketIO::sendMessage(const Message &ttkNotUsed(msg)) const {
  return 0;
}
int ttk::WebSocketIO::queueMessage(const std::string &) {
  return 0;
}
int ttk::WebSocketIO::queueMessage(const size_t &, const void *) {
  return 0;
}
int ttk::WebSocketIO::queueMessage(const Message &) {
  return 0;
}
int ttk::WebSocketIO::clearMessageQueue() {
  return 0;
}
int ttk::WebSocketIO::sendNextQueuedMessage() {
  return 0;
}
int ttk::WebSocketIO::processMessageQueue() {
  return 0;
}
int ttk::WebSocketIO::processEvent(const std::string &ttkNotUsed(eventName),
                                   const std::string &ttkNotUsed(eventData)) {
  return 0;
}

#endif
