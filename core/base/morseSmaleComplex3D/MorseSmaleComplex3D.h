/// \ingroup base
/// \class ttk::MorseSmaleComplex3D::
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex3D processing package.
///
/// %MorseSmaleComplex3D is a TTK processing package that takes a scalar field
/// on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex3D.cpp %for a usage example.

#pragma once

// base code includes
#include <AbstractMorseSmaleComplex.h>

#include <set>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 3D triangulation.
   */
  class MorseSmaleComplex3D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex3D();
  };
} // namespace ttk
