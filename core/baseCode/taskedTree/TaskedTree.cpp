/// \ingroup baseCode
/// \class ttk:TaskedTree
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

#include "TaskedTree.h"

// ------------
// CONSTRUCTOR
// ------------
// {

TaskedTree::TaskedTree() : ContourTree(new Params, nullptr, new Scalars)
{
}

TaskedTree::~TaskedTree()
{
}
