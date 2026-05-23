/// \ingroup base
/// \class ttk::VPath
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date May 2026
/// \date VPath extractor wrapping the DiscreteGradient class.
///
/// \brief TTK convenience class wrapping the DiscreteGradient class for
/// the easy extraction of vpaths.
///
/// Given a simplexId and dimension, this class returns a descending (or
/// ascending) vpath started in the given input simplex.
///
/// \sa VPath.cpp %for an alternative integral line backend.
/// \sa DiscreteGradient.cpp %for the core mechanisms.
/// \sa ttkVPath.cpp %for a usage example.
///

#pragma once

// base code includes
#include <Triangulation.h>
// std includes

namespace ttk {
  namespace vp {

  class VPath : virtual public Debug {

  public:
    VPath();
    ~VPath() override;

    // template <class triangulationType = ttk::AbstractTriangulation>
    // int execute(triangulationType *triangulation);

    /**
     * @brief Computes the integral line starting at the vertex of global id
     * seedIdentifier.
     *
     * @tparam triangulationType
     * @param triangulation
     * @param integralLine integral line to compute
     * @param offsets Order array of the scalar array
     */
    // template <class triangulationType = ttk::AbstractTriangulation>
    // void computeIntegralLine(const triangulationType *triangulation,
    //                          ttk::intgl::IntegralLine *integralLine,
    //                          const ttk::SimplexId *offsets) const;

  protected:

#endif
  };
} // namespace ttk
