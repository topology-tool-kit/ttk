/// \ingroup base
/// \class ttk::ProgressiveTopology
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date 2020.
///
/// \brief TTK processing package for progressive Topological Data Analysis
///
/// This package introduces a multiresolution hierarchical representation of the
/// data which allows the definition of efficient progressive algorithms for
/// TDA. It is applied to the progressive computation of Critical Points
/// and Persistence Diagrams.
///
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa PersistenceDiagram
/// \sa ScalarFieldCriticalPoints

#pragma once

// base code includes
#include <DynamicTree.h>
#include <MultiresTopology.h>

#include <limits>
#include <tuple>

namespace ttk {

  /**
   * @brief Persistence pair type (with persistence in double)
   */

  using triplet = std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>;
  using polarity = unsigned char;

  /**
   * Compute the persistence diagram of a function on a triangulation.
   * TTK assumes that the input dataset is made of only one connected component.
   */
  class ProgressiveTopology : public MultiresTopology {

  public:
    // struct PersistencePair {
    //   /** first (lower) vertex id */
    //   ttk::SimplexId birth{};
    //   /** second (higher) vertex id */
    //   ttk::SimplexId death{};
    //   /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
    //   ttk::SimplexId pairType{};

    //   PersistencePair() = default;
    //   PersistencePair(const SimplexId b,
    //                   const SimplexId d,
    //                   const SimplexId pType)
    //     : birth{b}, death{d}, pairType{pType} {
    //   }
    // };

    ProgressiveTopology() {
      this->setDebugMsgPrefix("ProgressiveTopology");
    }

    /* PROGRESSIVE MODE DECLARATIONS */
    int executeCPProgressive(int computePersistenceDiagram,
                             const SimplexId *inputOffsets);

    int resumeProgressive(int computePersistenceDiagram,
                          const SimplexId *offsets);

    inline void setAlgorithm(int data) {
      computePersistenceDiagram_ = data;
    }
    void setStartingDecimationLevel(int data) {
      if(data != startingDecimationLevel_) {
        resumeProgressive_ = false;
      }
      startingDecimationLevel_ = std::max(data, 0);
    }
    void setStoppingDecimationLevel(int data) {
      std::cout << "Stopping DL from progressive topo" << std::endl;
      if(data >= stoppingDecimationLevel_) {
        resumeProgressive_ = false;
      }
      stoppingDecimationLevel_ = std::max(data, 0);
    }
    inline void setIsResumable(const bool b) {
      this->isResumable_ = b;
      if(!b) {
        resumeProgressive_ = false;
      }
    }
    void setPreallocateMemory(const bool b) {
      if(b != preallocateMemory_) {
        resumeProgressive_ = false;
      }
      this->preallocateMemory_ = b;
    }
    inline void setTimeLimit(const double d) {
      if(d <= 0.0) {
        this->timeLimit_ = std::numeric_limits<double>::infinity();
      } else {
        this->timeLimit_ = d;
      }
    }

    int computeProgressivePD(std::vector<PersistencePair> &CTDiagram,
                             const SimplexId *offsets);

    int computeProgressiveCP(
      std::vector<std::pair<SimplexId, char>> *criticalPoints,
      const SimplexId *offsets);

  protected:
    void sortPersistenceDiagram2(std::vector<PersistencePair> &diagram,
                                 const SimplexId *const offsets) const;

    using VLBoundaryType
      = std::array<std::vector<std::pair<SimplexId, SimplexId>>, nLink_>;

    void
      initCriticalPoints(std::vector<polarity> &isNew,
                         std::vector<std::vector<std::pair<polarity, polarity>>>
                           &vertexLinkPolarity,
                         std::vector<polarity> &toProcess,
                         std::vector<DynamicTree> &link,
                         std::vector<uint8_t> &vertexLink,
                         VLBoundaryType &vertexLinkByBoundaryType,
                         std::vector<char> &vertexTypes,
                         const SimplexId *const offsets) const;

    void initSaddleSeeds(std::vector<polarity> &isNew,
                         std::vector<std::vector<std::pair<polarity, polarity>>>
                           &vertexLinkPolarity,
                         std::vector<polarity> &toPropageMin,
                         std::vector<polarity> &toPropageMax,
                         std::vector<polarity> &toProcess,
                         std::vector<DynamicTree> &link,
                         std::vector<uint8_t> &vertexLink,
                         VLBoundaryType &vertexLinkByBoundaryType,
                         std::vector<std::vector<SimplexId>> &saddleCCMin,
                         std::vector<std::vector<SimplexId>> &saddleCCMax,
                         const SimplexId *const offsets) const;

    void initPropagation(
      std::vector<polarity> &toPropageMin,
      std::vector<polarity> &toPropageMax,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
      std::vector<std::vector<SimplexId>> &saddleCCMin,
      std::vector<std::vector<SimplexId>> &saddleCCMax,
      std::vector<Lock> &vertLockMin,
      std::vector<Lock> &vertLockMax,
      std::vector<polarity> &isUpdatedMin,
      std::vector<polarity> &isUpdatedMax,
      const SimplexId *const offsets) const;

    void updatePropagation(
      std::vector<polarity> &toPropageMin,
      std::vector<polarity> &toPropageMax,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
      std::vector<std::vector<SimplexId>> &saddleCCMin,
      std::vector<std::vector<SimplexId>> &saddleCCMax,
      std::vector<Lock> &vertLockMin,
      std::vector<Lock> &vertLockMax,
      std::vector<polarity> &isUpdatedMin,
      std::vector<polarity> &isUpdatedMax,
      const SimplexId *const offsets) const;

    void
      buildVertexLinkPolarity(const SimplexId vertexId,
                              std::vector<std::pair<polarity, polarity>> &vlp,
                              const SimplexId *const offsets) const;

    void sortTriplets(std::vector<triplet> &triplets,
                      const SimplexId *const offsets,
                      const bool splitTree) const;

    void tripletsToPersistencePairs(
      std::vector<PersistencePair> &pairs,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<triplet> &triplets,
      const SimplexId *const offsets,
      const bool splitTree) const;

    void initDynamicLink(const SimplexId &vertexId,
                         std::vector<std::pair<polarity, polarity>> &vlp,
                         uint8_t &vertexLink,
                         DynamicTree &link,
                         VLBoundaryType &vlbt,
                         const SimplexId *const offsets) const;

    char getCriticalTypeFromLink(
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link) const;

    void
      updateDynamicLink(DynamicTree &link,
                        std::vector<std::pair<polarity, polarity>> &vlp,
                        std::vector<std::pair<SimplexId, SimplexId>> &vl) const;

    void updateCriticalPoints(
      std::vector<polarity> &isNew,
      std::vector<std::vector<std::pair<polarity, polarity>>>
        &vertexLinkPolarity,
      std::vector<polarity> &toProcess,
      std::vector<polarity> &toReprocess,
      std::vector<DynamicTree> &link,
      std::vector<uint8_t> &vertexLink,
      VLBoundaryType &vertexLinkByBoundaryType,
      std::vector<char> &vertexTypes,
      const SimplexId *const offsets) const;

    void
      updateSaddleSeeds(std::vector<polarity> &isNew,
                        std::vector<std::vector<std::pair<polarity, polarity>>>
                          &vertexLinkPolarity,
                        std::vector<polarity> &toPropageMin,
                        std::vector<polarity> &toPropageMax,
                        std::vector<polarity> &toProcess,
                        std::vector<polarity> &toReprocess,
                        std::vector<DynamicTree> &link,
                        std::vector<uint8_t> &vertexLink,
                        VLBoundaryType &vertexLinkByBoundaryType,
                        std::vector<std::vector<SimplexId>> &saddleCCMin,
                        std::vector<std::vector<SimplexId>> &saddleCCMax,
                        std::vector<polarity> &isUpdatedMin,
                        std::vector<polarity> &isUpdatedMax,
                        const SimplexId *const offsets) const;

    bool getMonotonyChangeByOldPointCP(
      const SimplexId vertexId,
      const std::vector<polarity> &isNew,
      std::vector<polarity> &toProcess,
      std::vector<polarity> &toReprocess,
      std::vector<std::pair<polarity, polarity>> &vlp,
      const SimplexId *const offsets) const;

    void updateLinkPolarity(const SimplexId vertexId,
                            std::vector<std::pair<polarity, polarity>> &vlp,
                            const SimplexId *const offsets) const;

    ttk::SimplexId propageFromSaddles(
      const SimplexId vertexId,
      std::vector<Lock> &vertLock,
      std::vector<polarity> &toPropage,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<std::vector<SimplexId>> &saddleCC,
      std::vector<polarity> &isUpdated,
      std::vector<SimplexId> &globalExtremum,
      const SimplexId *const offsets,
      const bool splitTree) const;

    void computePersistencePairsFromSaddles(
      std::vector<PersistencePair> &CTDiagram,
      const SimplexId *const offsets,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
      const std::vector<polarity> &toPropageMin,
      const std::vector<polarity> &toPropageMax) const;

    double predictNextIterationDuration(const double currItDuration,
                                        const size_t nCurrPairs) const;

    void stopComputationIf(const bool b);
    void clearResumableState();

    // progressive approach
    int computePersistenceDiagram_{1};

    // do some extra computations to allow to resume computation
    bool isResumable_{false};
    bool resumeProgressive_{false};
    // time limit
    double timeLimit_{0.0};
    bool preallocateMemory_{true};

    // keep state in case of resuming computation
    std::vector<std::vector<SimplexId>> vertexRepresentativesMax_{},
      vertexRepresentativesMin_{};
    std::vector<std::vector<std::pair<polarity, polarity>>>
      vertexLinkPolarity_{};
    std::vector<polarity> isNew_{};

    std::vector<uint8_t> vertexLink_{};
    VLBoundaryType vertexLinkByBoundaryType_{};
    std::vector<DynamicTree> link_{};
    std::vector<polarity> toProcess_{};
    std::vector<polarity> toReprocess_{};
    std::vector<std::vector<SimplexId>> saddleCCMin_{};
    std::vector<std::vector<SimplexId>> saddleCCMax_{};
    std::vector<char> vertexTypes_{};
  };
} // namespace ttk
