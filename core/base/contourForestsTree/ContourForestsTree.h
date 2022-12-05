/// \ingroup base
/// \class ttk::ContourForestsTree
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
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa ttkContourForests.cpp %for a usage example.

#pragma once

#include <memory>
#include <queue>
#include <set>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "DeprecatedDataTypes.h"
#include "ExtendedUF.h"
#include "MergeTree.h"

namespace ttk {
  namespace cf {
    class ContourForestsTree : public MergeTree {
      friend class ContourForests;

    protected:
      MergeTree jt_, st_;

    public:
      // -----------------
      // Constructors
      // -----------------
      // {

      ContourForestsTree(const std::shared_ptr<Params> &params,
                         const std::shared_ptr<Scalars> &scalars,
                         idPartition part = nullPartition);
      ~ContourForestsTree() override;

      // }
      // -----------------
      // INITIALIZE
      // -----------------
      // {

      void flush() {
        MergeTree::flush();
        jt_.flush();
        st_.flush();
      }

      // }
      // -----------------
      // ACCESSOR
      // -----------------
      // {

      inline MergeTree *getJoinTree() {
        return &jt_;
      }

      inline MergeTree *getSplitTree() {
        return &st_;
      }

      inline MergeTree *getTree(const TreeType &tt) {
        switch(tt) {
          case TreeType::JoinAndSplit:
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

      // }
      // -----------------
      // PROCESS
      // -----------------
      // {

      /// \brief Combine tree with Natarajan's algorithm
      int combine(const SimplexId &seed0,
                  const SimplexId &seed1,
                  std::list<std::vector<std::pair<SimplexId, bool>>> &storage);

    private:
      // -----------------
      // PROCESS
      // -----------------
      // {

      /// \brief initialize data of the Merge Trees jt & st
      template <typename scalarType>
      void initDataMT();

      // }
    };
  } // namespace cf
} // namespace ttk
