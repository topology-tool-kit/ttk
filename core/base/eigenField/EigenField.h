/// \ingroup base
/// \class ttk::EigenField
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date April 2019
///
/// \brief TTK processing package for computing eigenfunctions of a
/// triangular mesh.
///
/// \sa ttkEigenField.cpp % for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_casting/">Persistent
///   Generators Casting example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_fertility/">Persistent
///   Generators Fertility example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_skull/">Persistent
///   Generators Skull example</a> \n

#pragma once

// base code includes
#include <Debug.h>
#include <Geometry.h>
#include <Laplacian.h>
#include <Triangulation.h>

namespace ttk {

  class EigenField : virtual public Debug {
  public:
    EigenField() {
      this->setDebugMsgPrefix("EigenField");
    }

    inline void
      preconditionTriangulation(AbstractTriangulation &triangulation) const {
      Laplacian::preconditionTriangulation(triangulation);
    }

    template <typename T, class TriangulationType = AbstractTriangulation>
    int execute(const TriangulationType &triangulation,
                T *const outputFieldPointer,
                const unsigned int eigenNumber = 500,
                bool computeStatistics = false,
                T *const outputStatistics = nullptr) const;
  };

} // namespace ttk
