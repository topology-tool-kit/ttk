#include <Debug.h>

ttk::debug::LineMode ttk::Debug::lastLineMode = ttk::debug::LineMode::NEW;

bool ttk::welcomeMsg_ = true;
bool ttk::goodbyeMsg_ = true;
int ttk::globalDebugLevel_ = 0;

using namespace std;
using namespace ttk;

Debug::Debug() {

  setDebugMsgPrefix("Common");

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

  welcomeMsg(stream);

  if((debugLevel_ >= debugLevel) || (globalDebugLevel_ >= debugLevel))
    stream << msg.data() << flush;

  return 0;
}

int Debug::welcomeMsg(ostream &stream) const {

  if((ttk::welcomeMsg_) && (debugLevel_)) {
    ttk::welcomeMsg_ = false;
    stringstream s;

#include <welcomeLogo.inl>
#include <welcomeMsg.inl>

#ifndef NDEBUG
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "TTK has been built in debug mode!",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "DEVELOPERS ONLY!",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(
      debug::output::YELLOW + "Expect important performance degradation.",
      debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
#endif
#ifndef TTK_ENABLE_KAMIKAZE
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(
      debug::output::YELLOW + "TTK has *NOT* been built in performance mode!",
      debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "DEVELOPERS ONLY!",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(
      debug::output::YELLOW + "Expect important performance degradation.",
      debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);

    printMsg(debug::output::YELLOW
               + "To enable the performance mode, rebuild TTK with:",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "  -DTTK_ENABLE_KAMIKAZE=ON",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
#endif
  }

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
