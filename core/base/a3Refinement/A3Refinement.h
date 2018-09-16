/// \ingroup base
/// \class ttk::A3Refinement
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %a3Refinement processing package.
///
/// %A3Refinement is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkA3Refinement.cpp %for a usage example.

#pragma once

#include                  <Wrapper.h>

using namespace std;

namespace ttk{

  class A3Refinement : public Debug{

    public:

      A3Refinement();
      ~A3Refinement();

      template <class dataType> int execute() const;
  };
}

// if the package is a pure template class, uncomment the following line
// #include                  <A3Refinement.cpp>

// template functions
template <class dataType> int ttk::A3Refinement::execute() const{

  Timer t;

  {
    stringstream msg;
    msg << "[A3Refinement] Refined in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
