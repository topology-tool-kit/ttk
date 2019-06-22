/// \ingroup base
/// \class ttk::FTMTree_CT
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
/// \sa ttkContourForests.cpp %for a usage example.

#ifndef FTMTREE_CT_H
#define FTMTREE_CT_H

#include <queue>
#include <set>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "FTMDataTypes.h"
#include "FTMTree_MT.h"

namespace ttk {
  namespace ftm {

    class FTMTree_CT : public FTMTree_MT {
    protected:
      FTMTree_MT *jt_, *st_;

    public:
      // -----------------
      // Constructors
      // -----------------

      FTMTree_CT(Params *const params,
                 Triangulation *mesh,
                 Scalars *const scalars);
      virtual ~FTMTree_CT();

      // -----------------
      // ACCESSOR
      // -----------------

      inline FTMTree_MT *getJoinTree(void) const {
        return jt_;
      }

      inline FTMTree_MT *getSplitTree(void) const {
        return st_;
      }

      inline FTMTree_MT *getTree(const TreeType tt) {
        switch(tt) {
          case TreeType::Split:
            return getSplitTree();
            break;
          case TreeType::Join:
            return getJoinTree();
            break;
          case TreeType::Contour:
            return this;
            break;
          default:
            return this;
            break;
        }
        return this;
      }

      inline void setupTriangulation(Triangulation *m,
                                     const bool preproc = true) {
        FTMTree_MT::setupTriangulation(m, preproc);
        jt_->setupTriangulation(m, false);
        st_->setupTriangulation(m, false);
      }

      inline int setDebugLevel(const int &d) {
        Debug::setDebugLevel(d);
        jt_->setDebugLevel(d);
        st_->setDebugLevel(d);
        return 0;
      }

      inline int setThreadNumber(const int n) {
        Debug::setThreadNumber(n);
        jt_->setThreadNumber(n);
        st_->setThreadNumber(n);
        return 0;
      }

      // -----------------
      // PROCESS
      // -----------------

      int leafSearch();

      void build(TreeType tt);

      void insertNodes();

      int combine();

      void updateRegion(const ArcRegion &arcRegion, idSuperArc ctArc);

      void createCTArcSegmentation(idSuperArc ctArc,
                                   const bool isJT,
                                   idSuperArc xtArc);

      void finalizeSegmentation(void);
    };

  } // namespace ftm
} // namespace ttk

#endif // CONTOURTREE_H
