/// \ingroup baseCode
/// \class ttk::TaskedTree
/// \author Charles Gueuent <charles.gueunet@lip6.fr>
/// \date December 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkTaskedTree.cpp %for a usage example.

#ifndef _TASKEDTREE_H_
#define _TASKEDTREE_H_

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "ContourTree.h"
#include "DataTypes.h"

namespace ttk
{
   class TaskedTree : public ContourTree
   {
     public:
      // -----------------
      // CONSTRUCTORS
      // -----------------
      // {

      TaskedTree();
      virtual ~TaskedTree();

      // }

      // -------
      // PROCESS
      // -------
      // {

      // Initialize structures then build tree
      // Need triangulation, scalars and all params set before call
      template <typename scalarType>
      void build(void);

      // }
   };

#include "TaskedTreeTemplate.h"
}

#endif  // TASKEDTREE_H
