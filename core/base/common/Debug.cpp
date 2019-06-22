#include <Debug.h>

bool ttk::welcomeMsg_ = true;
bool ttk::goodbyeMsg_ = true;
int ttk::globalDebugLevel_ = 0;

using namespace std;
using namespace ttk;

Debug::Debug() {

  debugLevel_ = ttk::globalDebugLevel_;

  // avoid warnings
  if(goodbyeMsg_)
    goodbyeMsg_ = true;
}

Debug::~Debug() {
  if((lastObject_) && (ttk::goodbyeMsg_)) {
    stringstream msg;
    msg << "[Common] Goodbye :)" << endl;
    dMsg(cout, msg.str(), 1);
    ttk::goodbyeMsg_ = false;
  }
}

int Debug::dMsg(ostream &stream, string msg, const int &debugLevel) const {

  if((ttk::welcomeMsg_) && (debugLevel_)) {
    ttk::welcomeMsg_ = false;
    stringstream s;

    s << "[Common] "
         " _____ _____ _  __                       __  __    ____   ___  _  ___"
      << endl
      << "[Common] "
         "|_   _|_   _| |/ /                      / /__\\ \\  |___ \\ / _ \\/ "
         "|/ _ \\"
      //"|_   _|_   _| |/ /                      / /__\\ \\  |___ \\ / _ \\/ |(
      //_ )"
      << endl
      << "[Common] "
         "  | |   | | | ' /                      | |/ __| |   __) | | | | | "
         "(_) |"
      //"  | |   | | | ' /                      | |/ __| |   __) | | | | |/ _
      //\\"
      << endl
      << "[Common] "
         "  | |   | | | . \\                      | | (__| |  / __/| |_| | "
         "|\\__, |"
      //"  | |   | | | . \\                      | | (__| |  / __/| |_| | | (_)
      //|"
      << endl
      << "[Common] "
         "  |_|   |_| |_|\\_\\                     | |\\___| | "
         "|_____|\\___/|_|  /_/"
      //"  |_|   |_| |_|\\_\\                     | |\\___| |
      //|_____|\\___/|_|\\___/"
      << endl
      << "[Common] "
         "                                        \\_\\  /_/"
      //"                                        \\_\\  /_/"
      << endl;
    s << "[Common] Welcome!" << endl;
    dMsg(cout, s.str(), 1);
  }

  if((debugLevel_ >= debugLevel) || (globalDebugLevel_ >= debugLevel))
    stream << msg.data() << flush;
  return 0;
}

int Debug::err(const string msg, const int &debugLevel) const {
  return dMsg(cerr, msg, 0);
}

int Debug::msg(const char *msg, const int &debugLevel) const {
  return dMsg(cout, string(msg), debugLevel);
}

int Debug::setDebugLevel(const int &debugLevel) {
  debugLevel_ = debugLevel;
  return 0;
}

int Debug::setWrapper(const Wrapper *wrapper) {

  BaseClass::setWrapper(wrapper);

  setDebugLevel(reinterpret_cast<Debug *>(wrapper_)->debugLevel_);

  return 0;
}
