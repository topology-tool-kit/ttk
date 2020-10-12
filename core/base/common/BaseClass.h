/// \ingroup base
/// \class ttk::BaseClass
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date June 2018.
///
///\brief TTK base package.

#ifndef _BASECLASS_H
#define _BASECLASS_H

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
// enable the `and` and `or` keywords on MSVC
#include <ciso646>
#endif // _WIN32

#include <DataTypes.h>

#if defined(_MSC_VER) && defined(common_EXPORTS)
#define COMMON_EXPORTS __declspec(dllexport)
#elif defined(_MSC_VER)
#define COMMON_EXPORTS __declspec(dllimport)
#else
#define COMMON_EXPORTS
#endif // _MSC_VER

namespace ttk {

  COMMON_EXPORTS extern int globalThreadNumber_;

  class Wrapper;

  class BaseClass {
  public:
    BaseClass();

    virtual ~BaseClass(){};

    int getThreadNumber() const {
      return threadNumber_;
    };

    virtual int setThreadNumber(const int threadNumber) {
      threadNumber_ = threadNumber;
      return 0;
    }

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
    bool lastObject_;
    mutable int threadNumber_;
    Wrapper *wrapper_;
  };
} // namespace ttk

#endif // _BASECLASS_H
