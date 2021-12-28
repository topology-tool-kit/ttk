/// \namespace ttk The Topology ToolKit

/// \mainpage TTK 1.1 Documentation
/// \image html "splash.png"
/// Useful links:
///   - TTK Home:
/// <a href="https://topology-tool-kit.github.io/"
/// target="new">http://topology-tool-kit .github.io/</a>

/// \defgroup base base
/// \brief The Topology ToolKit - Base code processing packages.
/// @{
/// \class ttk::Debug
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2011.
///
/// \brief Minimalist debugging class.
///
/// %Debug provides a few mechanisms to handle debugging messages at a global
/// and local scope, time and memory measurements, etc.
/// Each ttk class should inheritate from it.

#pragma once

#include <BaseClass.h>

#include <algorithm>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ttk {

  COMMON_EXPORTS extern bool welcomeMsg_;
  COMMON_EXPORTS extern bool goodbyeMsg_;
  COMMON_EXPORTS extern int globalDebugLevel_;

  namespace debug {
    enum class Priority : int {
      ERROR, // 0
      WARNING, // 1
      PERFORMANCE, // 2
      INFO, // 3
      DETAIL, // 4
      VERBOSE // 5
    };

    enum class Separator : char {
      L0 = '%',
      L1 = '=',
      L2 = '-',
      SLASH = '/',
      BACKSLASH = '\\'
    };

    enum class LineMode : int {
      NEW,
      APPEND, // append
      REPLACE // replace line and append
    };

    namespace output {
      const std::string BOLD = "\33[0;1m";
      const std::string GREY = "\33[2;1m";
      const std::string ITALIC = "\33[3;1m";
      const std::string UNDERLINED = "\33[4;1m";
      const std::string FLASHING = "\33[5;1m";
      const std::string INVERTED = "\33[7;1m";
      const std::string STRIKETHROUGH = "\33[9;1m";
      const std::string DARKGREY = "\33[30;1m";
      const std::string RED = "\33[31;1m";
      const std::string GREEN = "\33[32;1m";
      const std::string YELLOW = "\33[33;1m";
      const std::string BLUE = "\33[34;1m";
      const std::string PINK = "\33[35;1m";
      const std::string LIGHTBLUE = "\33[36;1m";
      const std::string BRIGHTWHITE = "\33[37;1m";
      const std::string ENDCOLOR = "\33[0m";
    } // namespace output

    const int LINEWIDTH = 80;
  } // namespace debug

  class Debug : public BaseClass {

  public:
    // 1) constructors, destructors, operators, etc.
    Debug();

    virtual ~Debug();

    // 2) functions
    /// Set the debug level of a particular object. The global variable
    /// globalDebugLevel_ will over-ride this setting if it has a lower value.
    /// \return Returns 0 upon success, negative values otherwise.
    virtual int setDebugLevel(const int &debugLevel);

    /// Specify a pointer to a calling object that wraps the current class
    /// deriving from ttk::BaseClass.
    ///
    /// This function is useful to pass the execution context (debug level,
    /// number of threads, etc.) from a wrapper to a base object.
    /// \param wrapper Pointer to the wrapping object.
    /// \return Returns 0 upon success, negative values otherwise.
    virtual int setWrapper(const Wrapper *wrapper);

    // =========================================================================
    // New Debug Methods
    // =========================================================================

    /**
     * Prints a string debug message.
     */
    inline int printMsg(const std::string &msg,
                        const debug::Priority &priority = debug::Priority::INFO,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      return this->printMsgInternal(msg, priority, lineMode, stream);
    }

    /**
     * Prints multiple string debug messages at once.
     */
    inline int printMsg(const std::vector<std::string> &msgs,
                        const debug::Priority &priority = debug::Priority::INFO,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      size_t prints = 0;
      for(auto &msg : msgs)
        prints += this->printMsgInternal(msg, priority, lineMode, stream);
      return prints == msgs.size() ? 1 : 0;
    }

    /**
     * Prints an error debug message.
     */
    inline int printErr(const std::string &msg,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cerr) const {
      return this->printMsgInternal(
        msg, debug::Priority::ERROR, lineMode, stream);
    }

    /**
     * Prints a warning debug message.
     */
    inline int printWrn(const std::string &msg,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cerr) const {
      return this->printMsgInternal(
        msg, debug::Priority::WARNING, lineMode, stream);
    }

    /**
     * Prints a performance debug message with specified progress, time,
     * nThreads, and memory (values can be omitted form the message by passing
     * the value -1)
     */
    inline int printMsg(const std::string &msg,
                        const double &progress,
                        const double &time,
                        const int &threads,
                        const double &memory,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority
                        = debug::Priority::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      std::vector<std::string> chunks(4);
      size_t q = 0;

      if(memory >= 0)
        chunks[q++] = std::to_string(static_cast<int>(memory)) + "MB";
      if(time >= 0) {
        std::stringstream sStream;
        sStream.precision(3);
        sStream << std::fixed;
        sStream << time;
        chunks[q++] = sStream.str() + "s";
      }
      if(threads >= 0)
        chunks[q++] = std::to_string(threads) + "T";
      if(progress >= 0)
        chunks[q++] = std::to_string((int)(progress * 100)) + "%";

      std::string stats = "";
      if(q > 0) {
        stats += " [";
        stats += chunks[0];
        for(size_t i = 1; i < q; i++)
          stats += "|" + chunks[i];

        stats += "]";
      }

      return this->printMsgInternal(
        msg, stats, msg.length() < 1 ? ">" : ".", priority, lineMode, stream);
    }

    /**
     * Prints a performance debug message with specified progress and time.
     */
    inline int printMsg(const std::string &msg,
                        const double &progress,
                        const double &time,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority
                        = debug::Priority::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      return this->printMsg(
        msg, progress, time, -1, -1, lineMode, priority, stream);
    }

    /**
     * Prints a performance debug message with specified progress, time, and
     * threads.
     */
    inline int printMsg(const std::string &msg,
                        const double &progress,
                        const double &time,
                        const int &threads,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority
                        = debug::Priority::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      return this->printMsg(
        msg, progress, time, threads, -1, lineMode, priority, stream);
    }

    /**
     * Prints a performance debug message with specified progress.
     */
    inline int printMsg(const std::string &msg,
                        const double &progress,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority
                        = debug::Priority::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      return this->printMsg(
        msg, progress, -1, -1, -1, lineMode, priority, stream);
    }

    /**
     * Prints a performance debug message with specified progress and custom
     * priority.
     */
    inline int printMsg(const std::string &msg,
                        const double &progress,
                        const debug::Priority &priority,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cout) const {
      return this->printMsg(msg, progress, -1, -1, lineMode, priority, stream);
    }

    /**
     * Prints a table.
     */
    inline int printMsg(const std::vector<std::vector<std::string>> &rows,
                        const debug::Priority &priority = debug::Priority::INFO,
                        const bool hasHeader = true,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      int nRows = rows.size();
      int nCols = nRows > 0 ? rows[0].size() : 0;
      if(nCols < 1)
        return 0;

      std::vector<std::string> formatedRows(nRows);
      std::vector<size_t> colSizes(nCols, 0);
      for(int i = 0; i < nRows; i++)
        for(int j = 0; j < nCols; j++)
          colSizes[j] = std::max(colSizes[j], rows[i][j].size());

      auto formatCell = [](const std::string &value, const size_t &width,
                           const std::string &fillSymbol) {
        std::string cell = value;
        int diff = width - cell.size();
        for(int i = 0; i < diff; i++)
          cell += fillSymbol;
        return cell;
      };

      // Values
      int resultIndex = 0;
      for(int i = 0; i < nRows; i++) {
        auto &row = formatedRows[resultIndex++];
        row
          = formatCell(rows[i][0], colSizes[0], " ") + (hasHeader ? ": " : "");
        if(nCols > 1)
          row += formatCell(rows[i][1], colSizes[1], " ");
        for(int j = 2; j < nCols; j++)
          row += "," + formatCell(rows[i][j], colSizes[j], " ");
      }

      return this->printMsg(formatedRows, priority, lineMode, stream);
    }

    /**
     * Prints a separator.
     */
    inline int printMsg(const debug::Separator &separator,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority = debug::Priority::INFO,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      return this->printMsgInternal(
        "", "", std::string(1, (char &)separator), priority, lineMode, stream);
    }

    /**
     * Prints a separator with custom priority.
     */
    inline int printMsg(const debug::Separator &separator,
                        const debug::Priority &priority,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        std::ostream &stream = std::cout) const {
      return this->printMsg(separator, lineMode, priority, stream);
    }

    /**
     * Prints a message and fills the remaining space with a separator.
     */
    inline int printMsg(const std::string &msg,
                        const debug::Separator &separator,
                        const debug::LineMode &lineMode = debug::LineMode::NEW,
                        const debug::Priority &priority = debug::Priority::INFO,
                        std::ostream &stream = std::cout) const {
      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      return this->printMsgInternal(
        msg, "", std::string(1, (char &)separator), priority, lineMode, stream);
    }

    /**
     * Sets the prefix that will be print at the beginning of every
     * debug message.
     */
    inline void setDebugMsgPrefix(const std::string &prefix) {
      this->debugMsgPrefix_ = prefix.length() > 0 ? "[" + prefix + "] " : "";
    }

  protected:
    mutable int debugLevel_;

    COMMON_EXPORTS static debug::LineMode lastLineMode;

    std::string debugMsgPrefix_;

    /**
     * Internal debug method that formats debug messages.
     */
    inline int
      printMsgInternal(const std::string &msg,
                       const std::string &right,
                       const std::string &filler,
                       const debug::Priority &priority = debug::Priority::INFO,
                       const debug::LineMode &lineMode = debug::LineMode::NEW,
                       std::ostream &stream = std::cout) const {

      std::string combinedMsg = msg;

      if(filler.length() > 0) {
        if(msg.length() > 0)
          combinedMsg += " ";

        int gapWidth = debug::LINEWIDTH - this->debugMsgPrefix_.length()
                       - combinedMsg.length() - right.length();
        gapWidth = std::max(gapWidth / filler.length(), (size_t)1);

        for(int i = 0; i < gapWidth; i++)
          combinedMsg += filler;

        combinedMsg += debug::output::BLUE + right + debug::output::ENDCOLOR;
      }

      return this->printMsgInternal(combinedMsg, priority, lineMode, stream);
    }

    /**
     * Internal debug method that actually prints messages.
     */
    inline int printMsgInternal(const std::string &msg,
                                const debug::Priority &priority,
                                const debug::LineMode &lineMode,
                                std::ostream &stream = std::cout) const {

      if((this->debugLevel_ < (int)priority)
         && (globalDebugLevel_ < (int)priority))
        return 0;

      // on error or warning print end of line
      if((int)priority < 2 && this->lastLineMode == debug::LineMode::REPLACE)
        stream << "\n";

      // print prefix
      if(lineMode != debug::LineMode::APPEND)
        stream << debug::output::GREEN << this->debugMsgPrefix_
               << debug::output::ENDCOLOR;

      // print error or warning prefix
      if((int)priority == 0)
        stream << debug::output::RED << "[ERROR]" << debug::output::ENDCOLOR
               << " ";
      else if((int)priority == 1)
        stream << debug::output::YELLOW << "[WARNING]"
               << debug::output::ENDCOLOR << " ";

      // print msg
      stream << msg.data();

      // go either into new line or replace current line
      if(lineMode == debug::LineMode::NEW)
        stream << "\n";
      else if(lineMode == debug::LineMode::REPLACE)
        stream << "\r";

      // flush stream
      stream.flush();

      this->lastLineMode = lineMode;

      return 1;
    }

    int welcomeMsg(std::ostream &stream);
  };

} // namespace ttk

#include <Os.h>
#include <Timer.h>

namespace ttk {
  /// \brief Legacy backward compatibility
  class DebugTimer : public Timer {};
  /// \brief Legacy backward compatibility.
  class DebugMemory : public Memory {};
} // namespace ttk

#include <OrderDisambiguation.h>

/// @}
