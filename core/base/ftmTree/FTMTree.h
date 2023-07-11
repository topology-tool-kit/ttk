/// \ingroup base
/// \class ttk::ftm::FTMTree
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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">Contour
///   Tree Alignment example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/ctBones/">CT Bones
///   example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a>\n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeTemporalReduction/">Merge
///   Tree Temporal Reduction</a> \n

#pragma once

// base code includes
#include <Geometry.h>
#include <Triangulation.h>

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
      ~FTMTree() override = default;

      // -------
      // PROCESS
      // -------

      // Initialize structures then build tree
      // Need triangulation, scalars and all params set before call
      template <typename scalarType, class triangulationType>
      void build(const triangulationType *mesh);
    };

#include "FTMTree_Template.h"

  } // namespace ftm
} // namespace ttk
