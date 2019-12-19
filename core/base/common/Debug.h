/// \namespace ttk The Topology ToolKit

/// \mainpage TTK 0.9.9 Documentation
/// \image html "../img/splash.png"
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

#ifndef _DEBUG_H
#define _DEBUG_H

#include <BaseClass.h>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ttk {

  extern bool welcomeMsg_;
  extern bool goodbyeMsg_;
  extern int globalDebugLevel_;

  namespace debug {
    enum class PRIORITY : int {
      ERROR, // 0
      WARNING, // 1
      PERFORMANCE, // 2
      INFO, // 3
      DETAIL, // 4
      VERBOSE // 5
    };

    enum class SEPARATOR : char {
      L0 = '%',
      L1 = '=',
      L2 = '-',
      SLASH = '/',
      BACKSLASH = '\\'
    };

    enum class LINEMODE : int {
      NEW,
      APPEND, // append
      REPLACE // replace line and append
    };

    const int LINEWIDTH = 80;
  }; // namespace DEBUG

  class Debug : virtual public BaseClass {

  public:
    // 1) constructors, destructors, operators, etc.
    Debug();

    virtual ~Debug();

    enum debugPriority {
      fatalMsg, // 0
      timeMsg, // 1
      memoryMsg, // 2
      infoMsg, // 3
      detailedInfoMsg, // 4
      advancedInfoMsg // 5
    };

    // 2) functions
    /// Send a debug message to a stream with a priority debugLevel (lower
    /// means higher priority).
    /// If the global debug level for the program is set to 0, the program
    /// should be completely quiet. So the '0' priority should only be
    /// reserved for fatal errors.
    /// \param stream Output stream.
    /// \param msg %Debug message (can contain std::endl characters).
    /// \param debugLevel Priority of the message.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa msg(), err()
    virtual int dMsg(std::ostream &stream,
                     std::string msg,
                     const int &debugLevel = infoMsg) const;

    /// Wrapper for dMsg() that sends a debug message to the standard error
    /// output stream.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa dMsg(), msg()
    int err(const std::string msg, const int &debugLevel = fatalMsg) const;

    /// Wrapper for dMsg() that sends a debug message to the standard
    /// output stream.
    /// \return Returns 0 upon success, negative values otherwise.
    /// \sa dMsg(), msg()
    int msg(const char *msg, const int &debugLevel = infoMsg) const;

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
    /// \sa ttkBlank
    virtual int setWrapper(const Wrapper *wrapper);

    // =========================================================================
    // New Debug Methods
    // =========================================================================
    inline int PrintMsg(const std::string &msg,
                        const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
        return 0;

      return this->PrintMsgInternal(msg, priority, lineMode, stream);
    }

    inline int PrintMsg(const std::vector<std::string> &msgs,
                        const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
        return 0;

      size_t prints = 0;
      for(auto &msg : msgs)
        prints += this->PrintMsgInternal(msg, priority, lineMode, stream);
      return prints == msgs.size() ? 1 : 0;
    }

    inline int PrintErr(const std::string &msg,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cerr) const {
      return this->PrintMsgInternal(
        msg, debug::PRIORITY::ERROR, lineMode, stream);
    }

    inline int PrintMsg(const std::string &msg,
                        const double &progress,
                        const double &time,
                        const double &memory,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        const debug::PRIORITY &priority
                        = debug::PRIORITY::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
        return 0;

      std::vector<std::string> chunks(3);
      size_t q = 0;

      if(memory >= 0)
        chunks[q++] = std::to_string(memory) + "mb";
      if(time >= 0)
        chunks[q++] = std::to_string(time) + "s";
      if(progress >= 0)
        chunks[q++] = std::to_string((int)(progress * 100)) + "%";

      std::string stats = "";
      if(q >= 0) {
        stats += " [";
        stats += chunks[0];
        for(size_t i = 1; i < q; i++)
          stats += "|" + chunks[i];

        stats += "]";
      }

      return this->PrintMsgInternal(
        msg, stats, msg.length() < 1 ? ">" : ".", priority, lineMode, stream);
    }

    inline int PrintMsg(const std::string &msg,
                        const double &progress,
                        const double &time,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        const debug::PRIORITY &priority
                        = debug::PRIORITY::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      return this->PrintMsg(
        msg, progress, time, -1, lineMode, priority, stream);
    }

    inline int PrintMsg(const std::string &msg,
                        const double &progress,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        const debug::PRIORITY &priority
                        = debug::PRIORITY::PERFORMANCE,
                        std::ostream &stream = std::cout) const {
      return this->PrintMsg(msg, progress, -1, -1, lineMode, priority, stream);
    }

    inline int PrintMsg(const std::string &msg,
                        const double &progress,
                        const debug::PRIORITY &priority,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cout) const {
      return this->PrintMsg(msg, progress, -1, -1, lineMode, priority, stream);
    }

    inline int PrintMsg(const std::vector<std::vector<std::string>> &rows,
                        const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                        const bool hasHeader = true,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
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
                           const std::string &fillSymbol = " ") {
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
        row = formatCell(rows[i][0], colSizes[0]) + (hasHeader ? ": " : "");
        if(nCols > 1)
          row += formatCell(rows[i][1], colSizes[1]);
        for(int j = 2; j < nCols; j++)
          row += "," + formatCell(rows[i][j], colSizes[j]);
      }

      return this->PrintMsg(formatedRows, priority, lineMode, stream);
    }

    inline int PrintMsg(const debug::SEPARATOR &separator,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
        return 0;

      return this->PrintMsgInternal(
        "", "", std::string(1, (char &)separator), priority, lineMode, stream);
    }

    inline int PrintMsg(const debug::SEPARATOR &separator,
                        const debug::PRIORITY &priority,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        std::ostream &stream = std::cout) const {
      return this->PrintMsg(separator, lineMode, priority, stream);
    }

    inline int PrintMsg(const std::string &msg,
                        const debug::SEPARATOR &separator,
                        const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                        const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                        std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);
      if((this->debugLevel_ < priorityAsInt)
         && (globalDebugLevel_ < priorityAsInt))
        return 0;

      return this->PrintMsgInternal(
        msg, "", std::string(1, (char &)separator), priority, lineMode, stream);
    }

    inline void SetDebugMsgPrefix(const std::string &prefix) {
      this->DebugMsgPrefix = prefix.length() > 0 ? "[" + prefix + "] " : "";
    }

  protected:
    mutable int debugLevel_;

    std::string DebugMsgPrefix;

    inline int
      PrintMsgInternal(const std::string &msg,
                       const std::string &right,
                       const std::string &filler,
                       const debug::PRIORITY &priority = debug::PRIORITY::INFO,
                       const debug::LINEMODE &lineMode = debug::LINEMODE::NEW,
                       std::ostream &stream = std::cout) const {
      std::string combinedMsg = msg;

      if(filler.length() > 0) {
        int gapWidth = debug::LINEWIDTH - this->DebugMsgPrefix.length()
                       - combinedMsg.length() - right.length();
        gapWidth = std::max(gapWidth / filler.length(), (size_t)1);

        for(int i = 0; i < gapWidth; i++)
          combinedMsg += filler;

        combinedMsg += "\33[34;1m" + right + "\33[0m";
      }

      return this->PrintMsgInternal(combinedMsg, priority, lineMode, stream);
    }

    inline int PrintMsgInternal(const std::string &msg,
                                const debug::PRIORITY &priority,
                                const debug::LINEMODE &lineMode,
                                std::ostream &stream = std::cout) const {
      int priorityAsInt = static_cast<int>(priority);

      // go either into new line or replace current line
      if(lineMode == debug::LINEMODE::NEW)
        stream << "\n";
      if(lineMode == debug::LINEMODE::REPLACE)
        stream << "\r";

      // print prefix
      if(lineMode != debug::LINEMODE::APPEND)
        stream << "\33[32;1m" << this->DebugMsgPrefix << "\33[0m";

      // print error or warning prefix
      if(priorityAsInt == 0)
        stream << "\33[31;1m[ERROR]\33[0m ";
      if(priorityAsInt == 1)
        stream << "\33[33m[WARNING]\33[0m ";

      // print msg
      stream << msg.data();

      // on error or warning print end of line
      if(priorityAsInt < 2)
        stream << "\n";

      // flush stream
      stream.flush();

      return 1;
    }
  };
} // namespace ttk

#include <Os.h>

namespace ttk {
  /// \brief Legacy backward compatibility
  class DebugTimer : public Timer {};
  /// \brief Legacy backward compatibility.
  class DebugMemory : public Memory {};
} // namespace ttk

#endif
/// @}
