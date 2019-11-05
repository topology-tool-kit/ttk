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

  class Debug : public BaseClass {

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

  protected:
    mutable int debugLevel_;
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
