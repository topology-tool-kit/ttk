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
#include <ImplicitTriangulation.h>
#include <MultiresTriangulation.h>
#include <OpenMPLock.h>

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
  class ProgressiveTopology : public Debug {

  public:
    struct PersistencePair {
      /** first (lower) vertex id */
      ttk::SimplexId birth{};
      /** second (higher) vertex id */
      ttk::SimplexId death{};
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      ttk::SimplexId pairType{};

      PersistencePair() = default;
      PersistencePair(const SimplexId b,
                      const SimplexId d,
                      const SimplexId pType)
        : birth{b}, death{d}, pairType{pType} {
      }
    };

    ProgressiveTopology() {
      this->setDebugMsgPrefix("ProgressiveTopology");
    }

    inline void setupTriangulation(ImplicitTriangulation *const data) {
      triangulation_ = data;
      multiresTriangulation_.setTriangulation(triangulation_);
    }

    /* PROGRESSIVE MODE DECLARATIONS */
    int executeCPProgressive(int computePersistenceDiagram,
                             const SimplexId *inputOffsets);

    int resumeProgressive(int computePersistenceDiagram,
                          const SimplexId *offsets);

    inline void setAlgorithm(int data) {
      computePersistenceDiagram_ = data;
    }
    inline void setStartingDecimationLevel(int data) {
      if(data != startingDecimationLevel_) {
        resumeProgressive_ = false;
      }
      startingDecimationLevel_ = std::max(data, 0);
    }
    inline void setStoppingDecimationLevel(int data) {
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
    inline void setPreallocateMemory(const bool b) {
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
    inline int getStoppingDecimationLevel() {
      return this->stoppingDecimationLevel_;
    }

    int computeProgressivePD(std::vector<PersistencePair> &CTDiagram,
                             const SimplexId *offsets);

    int computeProgressiveCP(
      std::vector<std::pair<SimplexId, char>> *criticalPoints,
      const SimplexId *offsets);

    void setStartingResolutionLevel(int rl) {
      this->setStartingDecimationLevel(multiresTriangulation_.RL_to_DL(rl));
    }

    void setStoppingResolutionLevel(int rl) {
      this->setStoppingDecimationLevel(multiresTriangulation_.RL_to_DL(rl));
    }

  protected:
    void sortPersistenceDiagram2(std::vector<PersistencePair> &diagram,
                                 const SimplexId *const offsets) const;

    // maximum link size in 3D
    static const size_t nLink_ = 27;
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

    void buildVertexLinkByBoundary(const SimplexId vertexId,
                                   VLBoundaryType &vlbt) const;

    void initDynamicLink(const SimplexId &vertexId,
                         std::vector<std::pair<polarity, polarity>> &vlp,
                         uint8_t &vertexLink,
                         DynamicTree &link,
                         VLBoundaryType &vlbt,
                         const SimplexId *const offsets) const;

    char getCriticalTypeFromLink(
      const SimplexId vertexId,
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link) const;

    void getValencesFromLink(
      const SimplexId vertexId,
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link,
      std::vector<polarity> &toPropageMin,
      std::vector<polarity> &toPropageMax,
      std::vector<std::vector<SimplexId>> &saddleCCMin,
      std::vector<std::vector<SimplexId>> &saddleCCMax) const;

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

    void getTripletsFromSaddles(
      const SimplexId vertexId,
      std::vector<triplet> &triplets,
      const std::vector<std::vector<SimplexId>> &vertexReps) const;

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
    std::string resolutionInfoString();

    ImplicitTriangulation *triangulation_{};

    // new non-progressive approach
    MultiresTriangulation multiresTriangulation_{};

    // store the two global extrema extracted from the whole dataset vertices
    // sorting operation
    mutable SimplexId globalMax_{}, globalMin_{};

    // progressive approach
    int computePersistenceDiagram_{1};
    int decimationLevel_{};
    int startingDecimationLevel_{};
    int stoppingDecimationLevel_{};
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
    std::vector<PersistencePair> CTDiagram_;
  };
} // namespace ttk
