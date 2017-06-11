/// \ingroup baseCode
/// \class ttk::ContourTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef _CONTOURTREE_H
#define _CONTOURTREE_H

#include <queue>
#include <set>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "DataTypes.h"
#include "ExtendedUF.h"
#include "MergeTree.h"

namespace ttk
{
   class ContourTree : public MergeTree
   {
     protected:
      MergeTree *jt_, *st_;

     public:
      // -----------------
      // Constructors
      // -----------------
      // {

      ContourTree(Params* const params, Triangulation* mesh, Scalars* const scalars);
      virtual ~ContourTree();

      // }
      // -----------------
      // ACCESSOR
      // -----------------
      // {

      inline MergeTree* getJoinTree(void) const
      {
         return jt_;
      }

      inline MergeTree* getSplitTree(void) const
      {
         return st_;
      }

      inline MergeTree* getTree(const TreeType& tt)
      {
         switch (tt) {
            case TreeType::Split:
               return getSplitTree();
               break;
            case TreeType::Join:
               return getJoinTree();
               break;
            case TreeType::Contour:
               return this;
               break;
         }
         return this;
      }

      inline void setupTriangulation(Triangulation* m, const bool preproc = true)
      {
         MergeTree::setupTriangulation(m, preproc);
         jt_->setupTriangulation(m, false);
         st_->setupTriangulation(m, false);
      }

      inline int setDebugLevel(const int d)
      {
         Debug::setDebugLevel(d);
         jt_->setDebugLevel(d);
         st_->setDebugLevel(d);
         return 0;
      }

      inline int setThreadNumber(const int n)
      {
         Debug::setThreadNumber(n);
         jt_->setThreadNumber(n);
         st_->setThreadNumber(n);
         return 0;
      }

      // }
      // -----------------
      // PROCESS
      // -----------------
      // {

      int vertexPrecomputation();

      void build(TreeType tt);

      void insertNodes();

      int combine();

      void updateRegion(const ArcRegion& arcRegion, idSuperArc ctArc);

      void createCTArcSegmentation(idSuperArc ctArc, const bool isJT, idSuperArc xtArc);

      void finalizeSegmentation(void);

      // }
   };
}

#include "ContourTreeTemplate.h"

#endif  // CONTOURTREE_H
