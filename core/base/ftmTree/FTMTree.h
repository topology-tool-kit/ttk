/// \ingroup base
/// \class ttk::FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date December 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkFTMTree.cpp %for a usage example.

#ifndef FTMTREE_H
#define FTMTREE_H

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "FTMDataTypes.h"
#include "FTMTree_CT.h"

namespace ttk {
  namespace ftm {

    /**
     * Compute the join tree, split tree or contour tree of a function on a
     * triangulation. TTK assumes that the input dataset is made of only one
     * connected component.
     */
    class FTMTree : public FTMTree_CT {
    public:
      // -----------------
      // CONSTRUCTORS
      // -----------------

      FTMTree();
      virtual ~FTMTree();

      // -------
      // PROCESS
      // -------

      // Initialize structures then build tree
      // Need triangulation, scalars and all params set before call
      template <typename scalarType, typename idType>
      void build(void);
    };

#include "FTMTree_Template.h"

  } // namespace ftm
} // namespace ttk

#endif // TASKEDTREE_H
