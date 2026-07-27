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

#pragma once

#include <queue>
#include <set>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>

#include "FTMDataTypes.h"
#include "FTMTree_MT.h"

namespace ttk {
  namespace ftm {

    class FTMTree_CT : public FTMTree_MT {
    protected:
      FTMTree_MT jt_, st_;

    public:
      // -----------------
      // Constructors
      // -----------------

      FTMTree_CT(const std::shared_ptr<Params> &params,
                 const std::shared_ptr<Scalars> &scalars);
      ~FTMTree_CT() override = default;

      // -----------------
      // ACCESSOR
      // -----------------

      inline FTMTree_MT *getJoinTree() {
        return &jt_;
      }

      inline FTMTree_MT *getSplitTree() {
        return &st_;
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
            break;
        }
        return this;
      }

      inline void preconditionTriangulation(AbstractTriangulation *tri,
                                            const bool preproc = true) {
        FTMTree_MT::preconditionTriangulation(tri, preproc);
        jt_.preconditionTriangulation(tri, false);
        st_.preconditionTriangulation(tri, false);
      }

      inline int setDebugLevel(const int &d) override {
        Debug::setDebugLevel(d);
        jt_.setDebugLevel(d);
        st_.setDebugLevel(d);
        return 0;
      }

      inline int setThreadNumber(const int n) override {
        Debug::setThreadNumber(n);
        jt_.setThreadNumber(n);
        st_.setThreadNumber(n);
        return 0;
      }

      // -----------------
      // PROCESS
      // -----------------

      template <class triangulationType>
      int leafSearch(const triangulationType *mesh);

      template <class triangulationType>
      void build(const triangulationType *mesh, TreeType tt);

      void insertNodes();

      int combine();

      void updateRegion(const ArcRegion &arcRegion, idSuperArc ctArc);

      void createCTArcSegmentation(idSuperArc ctArc,
                                   const bool isJT,
                                   idSuperArc xtArc);

      void finalizeSegmentation();
    };

  } // namespace ftm
} // namespace ttk

#include <FTMTree_CT_Template.h>
