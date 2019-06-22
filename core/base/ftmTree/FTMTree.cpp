/// \class ttk:FTMTree
/// \ingroup base
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#include "FTMTree.h"

using namespace std;
using namespace ttk;
using namespace ftm;

FTMTree::FTMTree() : FTMTree_CT(new Params, nullptr, new Scalars) {
}

FTMTree::~FTMTree() {
  delete params_;
  delete scalars_;
}
