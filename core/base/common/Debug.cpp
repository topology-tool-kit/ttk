#include <Debug.h>

COMMON_EXPORTS ttk::debug::LineMode ttk::Debug::lastLineMode
  = ttk::debug::LineMode::NEW;

COMMON_EXPORTS bool ttk::welcomeMsg_ = true;
COMMON_EXPORTS bool ttk::goodbyeMsg_ = true;
COMMON_EXPORTS int ttk::globalDebugLevel_ = 0;

using namespace std;
using namespace ttk;

Debug::Debug() {

  setDebugMsgPrefix("Debug");

  debugLevel_ = ttk::globalDebugLevel_;

  // avoid warnings
  if(goodbyeMsg_)
    goodbyeMsg_ = true;
}

Debug::~Debug() {
  if((lastObject_) && (ttk::goodbyeMsg_)) {

    printMsg(
      "Goodbye :)", debug::Priority::PERFORMANCE, debug::LineMode::NEW, cout);

    ttk::goodbyeMsg_ = false;
  }
}

int Debug::welcomeMsg(ostream &stream) {

  int priorityAsInt = (int)debug::Priority::PERFORMANCE;

  if((ttk::welcomeMsg_) && (debugLevel_ > priorityAsInt)) {
    ttk::welcomeMsg_ = false;

    string currentPrefix = debugMsgPrefix_;
    debugMsgPrefix_ = "[Common] ";

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
#ifndef TTK_ENABLE_DOUBLE_TEMPLATING
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW
               + "TTK has *NOT* been built in double-templating mode!",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "DEVELOPERS ONLY!",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(
      debug::output::YELLOW + "Expect unsupported types for bivariate data.",
      debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);

    printMsg(debug::output::YELLOW
               + "To enable the double-templating mode, rebuild TTK with:",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg(debug::output::YELLOW + "  -DTTK_ENABLE_DOUBLE_TEMPLATING=ON",
             debug::Priority::WARNING, debug::LineMode::NEW, stream);
    printMsg("", debug::Priority::WARNING, debug::LineMode::NEW, stream);
#endif

    debugMsgPrefix_ = currentPrefix;
  }

  return 0;
}

int Debug::setDebugLevel(const int &debugLevel) {
  debugLevel_ = debugLevel;

  welcomeMsg(cout);

  return 0;
}

int Debug::setWrapper(const Wrapper *wrapper) {

  BaseClass::setWrapper(wrapper);

  setDebugLevel(reinterpret_cast<Debug *>(wrapper_)->debugLevel_);

  return 0;
}
