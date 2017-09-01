/// \namespace ttk The Topology ToolKit

/// \mainpage TTK 0.9.2 Documentation
/// \image html "../img/splash.png"
/// Useful links:
///   - TTK Home: 
/// <a href="https://topology-tool-kit.github.io/" 
/// target="new">http://topology-tool-kit .github.io/</a>

/// \defgroup baseCode baseCode
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

#ifndef                 _DEBUG_H
#define                 _DEBUG_H

#ifdef withOpenMP
  #include              <omp.h>
#endif

#include                <cerrno>
#include                <fstream>
#include                <iostream>
#include                <sstream>
#include                <string>
#include                <vector>


using namespace std;

namespace ttk{

  extern bool welcomeMsg_;
  extern bool goodbyeMsg_;
  extern int globalDebugLevel_;
 
  class Wrapper;
  
  class Debug{
    
    public:
      
      // 1) constructors, destructors, operators, etc.
      Debug();
      
      virtual ~Debug();
      
      enum debugPriority{
        fatalMsg,         // 0
        timeMsg,          // 1
        memoryMsg,        // 2
        infoMsg,          // 3
        detailedInfoMsg,  // 4
        advancedInfoMsg   // 5
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
      virtual int dMsg(ostream &stream, string msg, 
        const int &debugLevel = infoMsg) const;
      
      /// Wrapper for dMsg() that sends a debug message to the standard error
      /// output stream.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa dMsg(), msg()
      int err(const string msg, const int &debugLevel = fatalMsg) const;
     
      int getThreadNumber() const { return threadNumber_;};

      /// Wrapper for dMsg() that sends a debug message to the standard 
      /// output stream.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa dMsg(), msg()
      int msg(const char *msg, const int &debugLevel = infoMsg) const;
     
      /// Set the debug level of a particular object. The global variable 
      /// globalDebugLevel_ will over-ride this setting if it has a lower value.
      /// \return Returns 0 upon success, negative values otherwise.
      virtual int setDebugLevel(const int &debugLevel);
      
      int setThreadNumber(const int threadNumber){
        
        threadNumber_ = threadNumber;
        return 0;
      }
      
      /// Specify a pointer to a calling object that wraps the current class 
      /// deriving from ttk::Debug.
      /// 
      /// This function is useful to pass the execution context (debug level,
      /// number of threads, etc.) from a wrapper to a baseCode object.
      /// \param wrapper Pointer to the wrapping object.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa vtkBlank
      int setWrapper(const Wrapper *wrapper);
      
      
    protected:
      
      bool                    lastObject_;
      mutable int             debugLevel_, threadNumber_;
      Wrapper                 *wrapper_;
  };
}

using namespace ttk;

#include                <Os.h>

namespace ttk{
  /// \brief Legacy backward compatibility
  class DebugTimer : public Timer{};
  /// \brief Legacy backward compatibility.
  class DebugMemory : public Memory{};
}

#endif
/// @}
