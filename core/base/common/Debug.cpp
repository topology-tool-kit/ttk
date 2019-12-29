#include <Debug.h>

ttk::debug::LineMode ttk::Debug::lastLineMode = ttk::debug::LineMode::NEW;

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

    s << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "                       "
      << debug::output::ENDCOLOR << debug::output::GREEN << ","
      << debug::output::ENDCOLOR << debug::output::GREEN
      << ",⌂µµ▒▒▒▒▒▒▒▒▒▒▒▒µµµ," << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN
      << "                  ,µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒µµ,"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "              " << debug::output::ENDCOLOR
      << debug::output::GREEN << "," << debug::output::ENDCOLOR
      << debug::output::GREEN << "µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒µ,"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "           ,µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░░░░▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒▒▒▒▒▒µ"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "         ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░░░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░  "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░░░░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░░░"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒▒▒▒µ"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "       ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░                    " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒µ" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "      ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░  ⌂µ▒▒▒▒▒▒▒▒▒▒▒▒∩         " << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒µ" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "    ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░ ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░▒       " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░▒░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒µ" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "   ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::DARKGREY << ",▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒▒▒▒▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "░░▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "▒▒░░▒     "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "  ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░  "
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▒"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓▓▓█████"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "█▓▄"
      << debug::output::ENDCOLOR << debug::output::GREY << "▄"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░▒   "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░░▒ "
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "  ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "╣"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓╜^"
      << debug::output::ENDCOLOR << debug::output::GREY << "^"
      << debug::output::ENDCOLOR << debug::output::GREY << ""
                                                           "^"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "^"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▀"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▀█"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "███▓▓▓"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓"
      << debug::output::ENDCOLOR << debug::output::GREY << "µ"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░    "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░ "
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "╙"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░           "
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░░"
      << debug::output::ENDCOLOR << debug::output::GREY << "T"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▀"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "██▓▓▓"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓"
      << debug::output::ENDCOLOR << debug::output::GREY << "µ"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "▒   "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "'                 " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▀" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "█▓▓▓╪" << debug::output::ENDCOLOR
      << debug::output::GREY << "▓" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░   " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░▒ " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░ " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "1" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "∩                     "
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "█"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓▓╪@"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "@    "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░▒▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "▒∩                       ░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "µ" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "█" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "@@" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "@" << debug::output::ENDCOLOR
      << debug::output::GREY << "h" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░  " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░▒ " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " ░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░                       " << debug::output::ENDCOLOR
      << debug::output::GREEN << "░ " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "▐" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "@" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "@" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒ " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░▒▒µ" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << " " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░                      ░░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "@" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▒" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░░▒▒µ" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "  ░░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY
      << "░░               " << debug::output::ENDCOLOR << debug::output::GREEN
      << "░" << debug::output::ENDCOLOR << debug::output::GREEN << "░░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "j"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "▓"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "@"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░░░░▒▒∩"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "  " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░         ▒░░░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "{"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "@"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░░▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "   " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░░░░░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "@"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░░▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "▒"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "    " << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░░░░░▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "▓"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "     " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "║▓" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓▓@▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::BRIGHTWHITE << "@"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░ "
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "      " << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "ß" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓" << debug::output::ENDCOLOR
      << debug::output::BRIGHTWHITE << "▓▓" << debug::output::ENDCOLOR
      << debug::output::GREEN << "@" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "a"
      << debug::output::ENDCOLOR << debug::output::GREY << "▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "M  "
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒▒▒▒▒▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREEN << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "      " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░▒░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "e" << debug::output::ENDCOLOR
      << debug::output::GREY << "M" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░  " << debug::output::ENDCOLOR
      << debug::output::BLUE << "," << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒▒▒▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "        " << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒" << debug::output::ENDCOLOR
      << debug::output::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "╓╜"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "º"
      << debug::output::ENDCOLOR << debug::output::BLUE << ".   "
      << debug::output::ENDCOLOR << debug::output::BLUE << "¿▒▒▒▒▒▒▒▒▒▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "         " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░░░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░░▒▒▒▒▒▒▒▒░░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR
      << debug::output::GREY << "╓" << debug::output::ENDCOLOR
      << debug::output::GREY << "m" << debug::output::ENDCOLOR
      << debug::output::LIGHTBLUE << "\"" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░      " << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒▒▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "           " << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░░░░░░░░" << debug::output::ENDCOLOR
      << debug::output::GREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░░░░░░░∞⌐*\"░" << debug::output::ENDCOLOR
      << debug::output::BLUE << ".        " << debug::output::ENDCOLOR
      << debug::output::BLUE << "╓▒▒▒▒▒▒▒▒" << debug::output::ENDCOLOR
      << debug::output::LIGHTBLUE << "▒▒" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::GREEN << "░" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "             " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░░░░░░░░░░░░░░░              ░░,▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "               " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "▒░░░░░░░░░░░░░░░░░░░░░░░░░░╓▒"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "                  " << debug::output::ENDCOLOR
      << debug::output::BLUE << "░" << debug::output::ENDCOLOR
      << debug::output::BLUE << "░▒░░░░░░░░░░░░░░░░░░░░α"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::LIGHTBLUE << "▒▒▒▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒"
      << debug::output::ENDCOLOR << debug::output::BLUE << "▒░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "                        "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::BLUE
      << "░░▒░░░░░░░░░░▒▒▒▒▒░░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░" << debug::output::ENDCOLOR
      << debug::output::DARKGREY << "░░" << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR
      << debug::output::GREEN << "                                 "
      << debug::output::ENDCOLOR << debug::output::BLUE << "░"
      << debug::output::ENDCOLOR << debug::output::GREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░░░░"
      << debug::output::ENDCOLOR << debug::output::DARKGREY << "░"
      << debug::output::ENDCOLOR << endl
      << debug::output::BOLD << "[Common] " << debug::output::ENDCOLOR << endl;

    s << debug::output::BOLD << "[Common] ";
    s << " _____ _____ _  __                    __  __    ____   ___ ____   ___"
      << endl
      << "[Common] "
         "|_   _|_   _| |/ /                   / /__"
         "\\ \\  |___ \\ / _ \\___ \\ / _ \\"
      << endl
      << "[Common] "
         "  | |   | | | ' /                   | |/ __| |   __) | | | |__) | | "
         "| |"
      << endl
      << "[Common] "
         "  | |   | | | . \\                   | | (__| |  / __/| |_| / __/| "
         "|_| "
         "|"
      << endl
      << "[Common] "
         "  |_|   |_| |_|\\_\\                  "
         "| |\\___| | |_____|\\___/_____|\\___/"
      << endl
      << "[Common] "
         "                                     \\_\\  /_/"
      << debug::output::ENDCOLOR << endl;
    s << debug::output::BOLD << "[Common] Welcome!" << debug::output::ENDCOLOR
      << endl;
#ifndef NDEBUG
    s << debug::output::YELLOW << "[Common]" << endl;
    s << "[Common] WARNING:" << endl;
    s << "[Common] TTK has been built in debug mode! (developers only)" << endl;
    s << "[Common] Expect important performance degradation." << endl;
    s << "[Common]" << debug::output::ENDCOLOR << endl;
#endif
#ifndef TTK_ENABLE_KAMIKAZE
    s << debug::output::YELLOW << "[Common]" << endl;
    s << "[Common] WARNING:" << endl;
    s << "[Common] TTK has *NOT* been built in performance mode!"
      << " (developers only)" << endl;
    s << "[Common] Expect important performance degradation." << endl;
    s << "[Common] " << endl;
    s << "[Common] To enable the performance mode, rebuild TTK with:" << endl;
    s << "[Common]   -DTTK_ENABLE_KAMIKAZE=ON" << endl;
    s << "[Common]" << debug::output::ENDCOLOR << endl;
#endif
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
