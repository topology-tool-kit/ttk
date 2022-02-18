/// \ingroup base
/// \class ttk::BaseClass
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date June 2018.
///
///\brief TTK base package.

#ifndef _BASECLASS_H
#define _BASECLASS_H

#ifdef _WIN32
// enable the `and` and `or` keywords on MSVC
#include <ciso646>
#endif // _WIN32

#include <DataTypes.h>
#include <OpenMP.h>

#if defined(_MSC_VER) && defined(TTK_ENABLE_SHARED_BASE_LIBRARIES)
#if defined(common_EXPORTS)
// building common.dll
#define COMMON_EXPORTS __declspec(dllexport)
// building every other .dll that include this header
#else
#define COMMON_EXPORTS __declspec(dllimport)
#endif // common_EXPORTS
#else
#define COMMON_EXPORTS
#endif // _MSC_VER && TTK_ENABLE_SHARED_BASE_LIBRARIES

/**
 * @brief Mark function/method parameters that are not used in the
 * function body at all.
 *
 * This is the TTK version of vtkNotUsed. It silences compiler
 * warnings about unused parameters when those parameters are not used
 * in the current function/method body but still necessary in the
 * function/method signature (for instance when implementing a
 * particular API with inheritance).
 *
 * Basically, it removes the name of the parameter at its definition
 * site (similar to commenting the parameter name).
 */
#define ttkNotUsed(x)

/**
 * @brief Force the compiler to use the function/method parameter.
 *
 * Some function/method parameters are used under some `#ifdef`
 * preprocessor conditions. This macro introduce a no-op statement
 * that uses those parameters, to silence compiler warnings in the
 * `#else` blocks. It can be inserted anywhere in the function body.
 */
#define TTK_FORCE_USE(x) (void)(x)

namespace ttk {

  COMMON_EXPORTS extern int globalThreadNumber_;

  class Wrapper;

  class BaseClass {
  public:
    BaseClass();

    virtual ~BaseClass() = default;

    int getThreadNumber() const {
      return threadNumber_;
    }

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
    virtual int setWrapper(const Wrapper *wrapper);

  protected:
    bool lastObject_;
    mutable int threadNumber_;
    Wrapper *wrapper_;
  };
} // namespace ttk

#endif // _BASECLASS_H
