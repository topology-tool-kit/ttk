/// \ingroup base
/// \class ttk::MorseSmaleComplex2D
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK %morseSmaleComplex2D processing package.
///
/// %MorseSmaleComplex2D is a TTK processing package that takes a scalar field
/// on the input and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMorseSmaleComplex2D.cpp %for a usage example.

#pragma once

// base code includes
#include <AbstractMorseSmaleComplex.h>

namespace ttk {

  /**
   * Class specialized in building the Morse-Smale complex
   * of 2D triangulation.
   */
  class MorseSmaleComplex2D : public AbstractMorseSmaleComplex {

  public:
    MorseSmaleComplex2D();
    ~MorseSmaleComplex2D();

    /**
     * Main function for computing the whole Morse-Smale complex.
     */
    template <typename dataType, typename idType, typename triangulationType>
    int execute(const triangulationType &triangulation);

    /**
     * Compute the descending 1-separatrices by reading into the discrete
     * gradient.
     */
    template <typename triangulationType>
    int getAscendingSeparatrices1(
      const std::vector<dcg::Cell> &criticalPoints,
      std::vector<Separatrix> &separatrices,
      std::vector<std::vector<dcg::Cell>> &separatricesGeometry,
      const triangulationType &triangulation) const;
  };
} // namespace ttk
