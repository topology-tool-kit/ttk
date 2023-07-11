/// \ingroup base
/// \class ttk::TopologicalSimplification
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK processing package for the topological simplification of scalar
/// data.
///
/// Given an input scalar field and a list of critical points to remove, this
/// class minimally edits the scalar field such that the listed critical points
/// disappear. This procedure is useful to speedup subsequent topological data
/// analysis when outlier critical points can be easily identified. It is
/// also useful for data simplification.
///
/// \b Related \b publications \n
/// "Generalized Topological Simplification of Scalar Fields on Surfaces" \n
/// Julien Tierny, Valerio Pascucci \n
/// Proc. of IEEE VIS 2012.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2012.
///
/// "Localized Topological Simplification of Scalar Data"
/// Jonas Lukasczyk, Christoph Garth, Ross Maciejewski, Julien Tierny
/// Proc. of IEEE VIS 2020.
/// IEEE Transactions on Visualization and Computer Graphics
///
/// \sa ttkTopologicalSimplification.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/BuiltInExample1/">BuiltInExample1
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">Contour
///   Tree Alignment example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/ctBones/">CT Bones
///   example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions/">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 0 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 1 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 2 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 3 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 4 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tribute/">Tribute
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/uncertainStartingVortex/">
///   Uncertain Starting Vortex example</a> \n

#pragma once

// base code includes
#include <Debug.h>
#include <LegacyTopologicalSimplification.h>
#include <LocalizedTopologicalSimplification.h>
#include <Triangulation.h>

#include <cmath>
#include <set>
#include <tuple>
#include <type_traits>

namespace ttk {

  class TopologicalSimplification : virtual public Debug {
  public:
    TopologicalSimplification();

    enum class BACKEND { LEGACY, LTS };
    /*
     * Either execute this file "legacy" algorithm, or the
     * lts algorithm. The choice depends on the value of the variable backend_.
     * Default is lts (localized).
     */
    template <typename dataType, typename triangulationType>
    int execute(const dataType *const inputScalars,
                dataType *const outputScalars,
                const SimplexId *const identifiers,
                const SimplexId *const inputOffsets,
                SimplexId *const offsets,
                const SimplexId constraintNumber,
                const bool addPerturbation,
                const triangulationType &triangulation);

    inline void setBackend(const BACKEND arg) {
      backend_ = arg;
    }

    inline int preconditionTriangulation(AbstractTriangulation *triangulation) {
      switch(backend_) {
        case BACKEND::LEGACY:
          legacyObject_.setDebugLevel(debugLevel_);
          legacyObject_.setThreadNumber(threadNumber_);
          legacyObject_.preconditionTriangulation(triangulation);
          break;

        case BACKEND::LTS:
          ltsObject_.setDebugLevel(debugLevel_);
          ltsObject_.setThreadNumber(threadNumber_);
          ltsObject_.preconditionTriangulation(triangulation);
          break;

        default:
          this->printErr(
            "Error, the backend for topological simplification is invalid");
          return -1;
      }
      return 0;
    }

  protected:
    BACKEND backend_{BACKEND::LTS};
    LegacyTopologicalSimplification legacyObject_;
    lts::LocalizedTopologicalSimplification ltsObject_;
  };
} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::TopologicalSimplification::execute(
  const dataType *const inputScalars,
  dataType *const outputScalars,
  const SimplexId *const identifiers,
  const SimplexId *const inputOffsets,
  SimplexId *const offsets,
  const SimplexId constraintNumber,
  const bool addPerturbation,
  const triangulationType &triangulation) {
  switch(backend_) {
    case BACKEND::LTS:
      return ltsObject_
        .removeUnauthorizedExtrema<dataType, SimplexId, triangulationType>(
          outputScalars, offsets, &triangulation, identifiers, constraintNumber,
          addPerturbation);
    case BACKEND::LEGACY:
      return legacyObject_.execute(inputScalars, outputScalars, identifiers,
                                   inputOffsets, offsets, constraintNumber,
                                   triangulation);

    default:
      this->printErr(
        "Error, the backend for topological simplification is invalid");
      return -1;
  }
}
