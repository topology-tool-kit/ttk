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

    s << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "                       " << debug::ENDCOLOR << debug::GREEN << ","
      << debug::ENDCOLOR << debug::GREEN << ",⌂µµ▒▒▒▒▒▒▒▒▒▒▒▒µµµ,"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "                  ,µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒µµ," << debug::ENDCOLOR
      << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "              " << debug::ENDCOLOR << debug::GREEN << ","
      << debug::ENDCOLOR << debug::GREEN
      << "µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒µ," << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "           ,µ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::GREY << "░░░░▒▒▒▒▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒▒▒▒▒▒µ" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "         ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░░░"
      << debug::ENDCOLOR << debug::BLUE << "░  " << debug::ENDCOLOR
      << debug::BLUE << "░░░░░░░░░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░░░" << debug::ENDCOLOR
      << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒µ" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "       ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░                    " << debug::ENDCOLOR
      << debug::BLUE << "░" << debug::ENDCOLOR << debug::BLUE << "░░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒▒µ" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "      ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY
      << "░" << debug::ENDCOLOR << debug::DARKGREY
      << "░  ⌂µ▒▒▒▒▒▒▒▒▒▒▒▒∩         " << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "░" << debug::ENDCOLOR
      << debug::GREEN << "▒▒µ" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "    ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░ ,▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░▒       "
      << debug::ENDCOLOR << debug::BLUE << "░" << debug::ENDCOLOR << debug::BLUE
      << "░▒░" << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒µ" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "   ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░ "
      << debug::ENDCOLOR << debug::DARKGREY << ",▒" << debug::ENDCOLOR
      << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY << "▒▒▒▒▒▒▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY
      << "░░▒░" << debug::ENDCOLOR << debug::DARKGREY << "▒" << debug::ENDCOLOR
      << debug::DARKGREY << "▒▒░░▒     " << debug::ENDCOLOR << debug::BLUE
      << "░" << debug::ENDCOLOR << debug::BLUE << "░░▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "  ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░  "
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY
      << "▒▒▒▒▒" << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::GREY << "▒▒▒" << debug::ENDCOLOR << debug::BRIGHTWHITE << "▒"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "▓▓▓█████" << debug::ENDCOLOR << debug::BRIGHTWHITE << "█▓▄"
      << debug::ENDCOLOR << debug::GREY << "▄" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░▒   "
      << debug::ENDCOLOR << debug::BLUE << "░" << debug::ENDCOLOR << debug::BLUE
      << "░░░░▒ " << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "░" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "  ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░ " << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::GREY << "▒▒▒" << debug::ENDCOLOR
      << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY << "▒▒▒▒"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "╣" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "▓╜^" << debug::ENDCOLOR << debug::GREY << "^"
      << debug::ENDCOLOR << debug::GREY
      << ""
         "^"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "^" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "▀" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "▀█" << debug::ENDCOLOR << debug::BRIGHTWHITE << "███▓▓▓"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR
      << debug::GREY << "µ" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░    " << debug::ENDCOLOR
      << debug::BLUE << "░" << debug::ENDCOLOR << debug::BLUE << "░░░ "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░ "
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY
      << "▒▒" << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::GREY << "▒▒▒" << debug::ENDCOLOR << debug::GREY << "╙"
      << debug::ENDCOLOR << debug::DARKGREY << "░           " << debug::ENDCOLOR
      << debug::DARKGREY << "░░" << debug::ENDCOLOR << debug::GREY << "T"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▀" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "██▓▓▓" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "▓" << debug::ENDCOLOR << debug::GREY << "µ" << debug::ENDCOLOR
      << debug::DARKGREY << "▒   " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░░░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "░" << debug::ENDCOLOR
      << debug::GREEN << "░" << debug::ENDCOLOR << debug::DARKGREY << "░ "
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY
      << "▒▒▒▒" << debug::ENDCOLOR << debug::DARKGREY << "'                 "
      << debug::ENDCOLOR << debug::DARKGREY << "░░" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "▀" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "█▓▓▓╪" << debug::ENDCOLOR << debug::GREY << "▓" << debug::ENDCOLOR
      << debug::DARKGREY << "░   " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░░▒ " << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "░░" << debug::ENDCOLOR << debug::DARKGREY << "░ "
      << debug::ENDCOLOR << debug::DARKGREY << "1" << debug::ENDCOLOR
      << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::DARKGREY << "∩                     "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::BRIGHTWHITE << "█"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓▓╪@" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@    " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░▒▒" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << " ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░" << debug::ENDCOLOR << debug::DARKGREY << "░ "
      << debug::ENDCOLOR << debug::DARKGREY << "▒∩                       ░"
      << debug::ENDCOLOR << debug::DARKGREY << "µ" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "█" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "▓" << debug::ENDCOLOR << debug::BRIGHTWHITE << "@@" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::GREY << "h"
      << debug::ENDCOLOR << debug::BLUE << "░  " << debug::ENDCOLOR
      << debug::BLUE << "░░░░▒ " << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << " ░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░                       "
      << debug::ENDCOLOR << debug::GREEN << "░ " << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREY << "▐"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "@" << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::GREY << "▒ " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░░░▒▒µ" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN << " "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << debug::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY
      << "░                      ░░" << debug::ENDCOLOR << debug::DARKGREY
      << "░░" << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::BRIGHTWHITE
      << "▒" << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::DARKGREY << "▒" << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░░░░▒▒µ" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "  ░░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY
      << "░░               " << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << debug::GREEN << "░░░░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "j"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR << debug::BLUE
      << "░" << debug::ENDCOLOR << debug::BLUE << "░░░░░░▒▒∩" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN << "  "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::GREEN << "░" << debug::ENDCOLOR << debug::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░         ▒░░░░░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "{" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░░░▒▒▒▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN << "   "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::GREEN << "░" << debug::ENDCOLOR << debug::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░░░░░░░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR
      << debug::GREY << "▒" << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::BLUE << "▒" << debug::ENDCOLOR << debug::BLUE << "░░░▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN << "    "
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░░░░░░░░▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "▒" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "▓" << debug::ENDCOLOR << debug::GREY
      << "▒" << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "░▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "     " << debug::ENDCOLOR << debug::BLUE << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << debug::GREEN << "║▓"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓▓@▒" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << debug::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::BRIGHTWHITE << "@" << debug::ENDCOLOR << debug::GREY << "▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░ " << debug::ENDCOLOR
      << debug::BLUE << "▒" << debug::ENDCOLOR << debug::BLUE << "▒▒▒▒▒▒▒▒▒░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "      " << debug::ENDCOLOR << debug::BLUE << "▒" << debug::ENDCOLOR
      << debug::BLUE << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "▒" << debug::ENDCOLOR
      << debug::GREEN << "░" << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "ß" << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓"
      << debug::ENDCOLOR << debug::BRIGHTWHITE << "▓▓" << debug::ENDCOLOR
      << debug::GREEN << "@" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "a" << debug::ENDCOLOR << debug::GREY
      << "▒" << debug::ENDCOLOR << debug::DARKGREY << "M  " << debug::ENDCOLOR
      << debug::BLUE << "▒" << debug::ENDCOLOR << debug::BLUE << "▒▒▒▒▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "      " << debug::ENDCOLOR << debug::BLUE << "░" << debug::ENDCOLOR
      << debug::BLUE << "░▒░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "░" << debug::ENDCOLOR
      << debug::GREEN << "░" << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::GREEN << "▒" << debug::ENDCOLOR
      << debug::GREEN << "▒▒▒▒" << debug::ENDCOLOR << debug::GREEN << "▒"
      << debug::ENDCOLOR << debug::GREEN << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░"
      << debug::ENDCOLOR << debug::GREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREY << "e" << debug::ENDCOLOR << debug::GREY
      << "M" << debug::ENDCOLOR << debug::DARKGREY << "░  " << debug::ENDCOLOR
      << debug::BLUE << "," << debug::ENDCOLOR << debug::BLUE << "▒▒▒▒▒▒▒▒▒▒▒"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "        " << debug::ENDCOLOR << debug::BLUE << "▒" << debug::ENDCOLOR
      << debug::BLUE << "░░░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░░" << debug::ENDCOLOR
      << debug::GREEN << "▒" << debug::ENDCOLOR << debug::GREEN
      << "▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREY << "╓╜"
      << debug::ENDCOLOR << debug::LIGHTBLUE << "º" << debug::ENDCOLOR
      << debug::BLUE << ".   " << debug::ENDCOLOR << debug::BLUE
      << "¿▒▒▒▒▒▒▒▒▒▒░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "         " << debug::ENDCOLOR << debug::BLUE << "░" << debug::ENDCOLOR
      << debug::BLUE << "░░░░░░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░░" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "░" << debug::ENDCOLOR
      << debug::GREEN << "░░▒▒▒▒▒▒▒▒░░" << debug::ENDCOLOR << debug::GREY << "░"
      << debug::ENDCOLOR << debug::GREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░░"
      << debug::ENDCOLOR << debug::GREY << "╓" << debug::ENDCOLOR << debug::GREY
      << "m" << debug::ENDCOLOR << debug::LIGHTBLUE << "\"" << debug::ENDCOLOR
      << debug::DARKGREY << "░      " << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "▒▒▒▒▒▒▒▒▒▒" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::GREEN << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "           " << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "░░░░░░░░░░░" << debug::ENDCOLOR
      << debug::GREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << debug::DARKGREY << "░░░░░░░░∞⌐*\"░"
      << debug::ENDCOLOR << debug::BLUE << ".        " << debug::ENDCOLOR
      << debug::BLUE << "╓▒▒▒▒▒▒▒▒" << debug::ENDCOLOR << debug::LIGHTBLUE
      << "▒▒" << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::GREEN << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "             " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE
      << "░░░░░░░░░░░░░░░              ░░,▒▒▒▒▒" << debug::ENDCOLOR
      << debug::LIGHTBLUE << "▒" << debug::ENDCOLOR << debug::LIGHTBLUE << "▒▒"
      << debug::ENDCOLOR << debug::BLUE << "▒" << debug::ENDCOLOR << debug::BLUE
      << "░" << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "               " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "▒░░░░░░░░░░░░░░░░░░░░░░░░░░╓▒"
      << debug::ENDCOLOR << debug::LIGHTBLUE << "▒" << debug::ENDCOLOR
      << debug::LIGHTBLUE << "▒▒▒▒▒" << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "░░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "                  " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░▒░░░░░░░░░░░░░░░░░░░░α"
      << debug::ENDCOLOR << debug::LIGHTBLUE << "▒" << debug::ENDCOLOR
      << debug::LIGHTBLUE << "▒▒▒▒" << debug::ENDCOLOR << debug::BLUE << "▒"
      << debug::ENDCOLOR << debug::BLUE << "▒░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░"
      << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "                        " << debug::ENDCOLOR << debug::BLUE << "░"
      << debug::ENDCOLOR << debug::BLUE << "░░▒░░░░░░░░░░▒▒▒▒▒░░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << debug::GREEN
      << "                                 " << debug::ENDCOLOR << debug::BLUE
      << "░" << debug::ENDCOLOR << debug::GREY << "░" << debug::ENDCOLOR
      << debug::DARKGREY << "░" << debug::ENDCOLOR << debug::DARKGREY << "░░░░"
      << debug::ENDCOLOR << debug::DARKGREY << "░" << debug::ENDCOLOR << endl
      << debug::BOLD << "[Common] " << debug::ENDCOLOR << endl;

    s << debug::BOLD << "[Common] ";
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
      << debug::ENDCOLOR << endl;
    s << debug::BOLD << "[Common] Welcome!" << debug::ENDCOLOR << endl;
#ifndef NDEBUG
    s << debug::YELLOW << "[Common]" << endl;
    s << "[Common] WARNING:" << endl;
    s << "[Common] TTK has been built in debug mode! (developers only)" << endl;
    s << "[Common] Expect important performance degradation." << endl;
    s << "[Common]" << debug::ENDCOLOR << endl;
#endif
#ifndef TTK_ENABLE_KAMIKAZE
    s << debug::YELLOW << "[Common]" << endl;
    s << "[Common] WARNING:" << endl;
    s << "[Common] TTK has *NOT* been built in performance mode!"
      << " (developers only)" << endl;
    s << "[Common] Expect important performance degradation." << endl;
    s << "[Common] " << endl;
    s << "[Common] To enable the performance mode, rebuild TTK with:" << endl;
    s << "[Common]   -DTTK_ENABLE_KAMIKAZE=ON" << endl;
    s << "[Common]" << debug::ENDCOLOR << endl;
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
