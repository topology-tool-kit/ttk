/// \ingroup base
/// \class ttk::Wrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2014.
///
/// \brief Wrapper class to wrap ttk code.

#ifndef _WRAPPER_H
#define _WRAPPER_H

#include <Debug.h>

namespace ttk {

  class Wrapper : public Debug {

  public:
    Wrapper() {
      processingProgress_ = 0;
    };

    virtual ~Wrapper(){};

    virtual bool needsToAbort() = 0;

    virtual int updateProgress(const float &progress) = 0;

  protected:
    float processingProgress_;
  };
} // namespace ttk

#endif
