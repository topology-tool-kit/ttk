/// \ingroup base
/// \class ttk::DistanceMatrixDistortion
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %DistanceMatrixDistortion class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa DistanceMatrixDistortion

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>

namespace ttk {

  /**
   * The DistanceMatrixDistortion class provides a method to compute the
   * distortion score between two distance matrices representing the same
   * points.
   */
  class DistanceMatrixDistortion : virtual public Debug {

  public:
    DistanceMatrixDistortion();

    int execute(const std::vector<double *> &highDistMatrix,
                const std::vector<double *> &lowDistMatrix,
                double &distortionValue,
                double *distortionVerticesValues) const;

  protected:
    // Warning: If this is set to true, then the similarity values we return
    // will lie between -inf and 1 (best similarity).
    bool DoNotNormalize{false};
  };

} // namespace ttk
