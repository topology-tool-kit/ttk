/// \ingroup base
/// \class ttk::Wrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2014.
///
/// \brief Wrapper class to wrap ttk code.

#pragma once

#include <Debug.h>

namespace ttk {

  class Wrapper : public Debug {

  public:
    Wrapper() {
      processingProgress_ = 0;
    }

    ~Wrapper() override = default;

    virtual bool needsToAbort() = 0;

    virtual int updateProgress(const float &progress) = 0;

  protected:
    float processingProgress_;
  };
} // namespace ttk
