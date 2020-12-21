/// \ingroup base
/// \class ttk::ContourForests
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

#include "ContourForestsTree.h"

namespace ttk {
  namespace cf {
    // Classes Interface
    //         ContourForests

    class Interface {
    private:
      // same size : number of conn. components.
      SimplexId seed_;

      // Overlap
      std::vector<SimplexId> lowerOverlap_;
      std::vector<SimplexId> upperOverlap_;

    public:
      Interface(const SimplexId &seed);

      // Getter & Setter
      // {

      inline const SimplexId &getSeed(void) const {
        return seed_;
      }

      inline std::vector<SimplexId> &getUpper(void) {
        return upperOverlap_;
      }

      inline std::vector<SimplexId> &getLower(void) {
        return lowerOverlap_;
      }
      inline SimplexId getNbUpper(void) const {
        return upperOverlap_.size();
      }

      inline SimplexId getNbLower(void) const {
        return lowerOverlap_.size();
      }

      inline void setSeed(const SimplexId &local_seed) {
        seed_ = local_seed;
      }

      inline void addUpper(const SimplexId &lb) {
        upperOverlap_.emplace_back(lb);
      }

      inline void addLower(const SimplexId &lb) {
        lowerOverlap_.emplace_back(lb);
      }

      inline void upReserve(const SimplexId &u) {
        upperOverlap_.reserve(u);
      }

      inline void loReserve(const SimplexId &l) {
        lowerOverlap_.reserve(l);
      }

      inline void appendUpper(const std::vector<SimplexId> &vertices) {
        upperOverlap_.insert(
          upperOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      inline void appendLower(const std::vector<SimplexId> &vertices) {
        lowerOverlap_.insert(
          lowerOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      inline void swapUpper(std::vector<SimplexId> &verts) {
        upperOverlap_.swap(verts);
      }

      inline void swapLower(std::vector<SimplexId> &verts) {
        lowerOverlap_.swap(verts);
      }

      // }
    };

    struct ParallelParams {
      numThread nbThreads;
      idInterface nbInterfaces;
      idPartition nbPartitions;
      int partitionNum;
      bool lessPartition;
    };

    struct ParallelData {
      std::vector<Interface> interfaces;
      std::vector<ContourForestsTree> trees;
    };

    class ContourForests : public ContourForestsTree {
    private:
      // global (one instance -> no pointer)
      ParallelParams parallelParams_;

      // local
      ParallelData parallelData_;

    public:
      ContourForests();

      virtual ~ContourForests();

      // Getters & Setters
      // {

      inline int setThreadNumber(const int nbThread) {
        if(nbThread) {
          parallelParams_.nbThreads = nbThread;
        } else {
          parallelParams_.nbThreads = OsCall::getNumberOfCores();
        }
        return 0;
      }

      inline void setPartitionNum(int p) {
        parallelParams_.partitionNum = p;
      }

      inline void setLessPartition(bool l) {
        parallelParams_.lessPartition = l;
      }

      // range of partitions, position of seeds , ...

      inline std::tuple<SimplexId, SimplexId>
        getJTRange(const idPartition &i) const {
        const SimplexId &start
          = (i == 0)
              ? 0
              : scalars_->sosOffsets[parallelData_.interfaces[i - 1].getSeed()];

        const SimplexId &end
          = (i == parallelParams_.nbInterfaces)
              ? scalars_->size
              : scalars_->sosOffsets[parallelData_.interfaces[i].getSeed()];

        return std::make_tuple(start, end);
      }

      inline std::tuple<SimplexId, SimplexId>
        getSTRange(const idPartition &i) const {
        const SimplexId &start
          = (i == parallelParams_.nbInterfaces)
              ? scalars_->size - 1
              : scalars_->sosOffsets[parallelData_.interfaces[i].getSeed()] - 1;

        const SimplexId &end
          = (i == 0)
              ? -1
              : scalars_->sosOffsets[parallelData_.interfaces[i - 1].getSeed()]
                  - 1;

        return std::make_tuple(start, end);
      }

      inline std::tuple<SimplexId, SimplexId>
        getSeedsPos(const idPartition &i) const {
        const SimplexId &seed0
          = (i == 0)
              ? -1
              : scalars_->sosOffsets[parallelData_.interfaces[i - 1].getSeed()];

        const SimplexId &seed1
          = (i == parallelParams_.nbInterfaces)
              ? nullVertex
              : scalars_->sosOffsets[parallelData_.interfaces[i].getSeed()];

        return std::make_tuple(seed0, seed1);
      }

      inline std::tuple<std::vector<SimplexId>, std::vector<SimplexId>>
        getOverlaps(const idPartition &i) {
        const std::vector<SimplexId> &lower
          = (i == 0) ? std::vector<SimplexId>()
                     : parallelData_.interfaces[i - 1].getLower();

        const std::vector<SimplexId> &upper
          = (i == parallelParams_.nbInterfaces)
              ? std::vector<SimplexId>()
              : parallelData_.interfaces[i].getUpper();

        return std::make_tuple(lower, upper);
      }

      idPartition vertex2partition(const SimplexId &v);

      // }

      // Init
      // {
      void initInterfaces(void);

      template <typename triangulationType>
      void initOverlap(const triangulationType &mesh);

      void initNbPartitions(void);

      //}
      // Process
      // {

      template <typename scalarType, typename triangulationType>
      int build(const triangulationType &mesh);

      template <typename scalarType, typename triangulationType>
      int
        parallelBuild(std::vector<std::vector<ExtendedUnionFind *>> &baseUF_JT,
                      std::vector<std::vector<ExtendedUnionFind *>> &baseUF_ST,
                      const triangulationType &mesh);

      void stitch(void);
      void stitchTree(const char tree);

      // replace distributed tree by a global one, will be removed
      void unify();
      void unifyTree(const char treetype);
      // }

      // Print
      // {
      void printDebug(DebugTimer &timer, const std::string &str);

      void printVectCT();
      // }
    };

  } // namespace cf
} // namespace ttk

#include <ContourForestsTemplate.h>
