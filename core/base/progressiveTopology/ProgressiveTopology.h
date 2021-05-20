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
#include <MultiresTriangulation.h>
#include <Triangulation.h>

#include <tuple>

#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#endif

namespace ttk {

  /**
   * @brief Persistence pair type (with persistence in double)
   */

#if defined(_GLIBCXX_PARALLEL_FEATURES_H) && defined(TTK_ENABLE_OPENMP)
#define PARALLEL_SORT                       \
  omp_set_num_threads(this->threadNumber_); \
  __gnu_parallel::sort
#else
#define PARALLEL_SORT std::sort
#endif // _GLIBCXX_PARALLEL_FEATURES_H && TTK_ENABLE_OPENMP

  using triplet = std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>;
  using polarity = unsigned char;

  /**
   * @brief RAII wrapper around OpenMP lock
   */
  class Lock {
#ifdef TTK_ENABLE_OPENMP
  public:
    Lock() {
      omp_init_lock(&this->lock_);
    }
    ~Lock() {
      omp_destroy_lock(&this->lock_);
    }
    inline void lock() {
      omp_set_lock(&this->lock_);
    }
    inline void unlock() {
      omp_unset_lock(&this->lock_);
    }
    Lock(const Lock &) = delete;
    Lock(Lock &&) = delete;
    Lock &operator=(const Lock &) = delete;
    Lock &operator=(Lock &&) = delete;

  private:
    omp_lock_t lock_{};
#else
  public:
    inline void lock() {
    }
    inline void unlock() {
    }
#endif // TTK_ENABLE_OPENMP
  };

  /**
   * Compute the persistence diagram of a function on a triangulation.
   * TTK assumes that the input dataset is made of only one connected component.
   */
  class ProgressiveTopology : public Debug {

  public:

    struct PersistencePair {
      /** first (lower) vertex id */
      ttk::SimplexId birth{};
      /** first vertex type */
      ttk::CriticalType birthType{};
      /** second (higher) vertex id */
      ttk::SimplexId death{};
      /** second vertex type */
      ttk::CriticalType deathType{};
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

    inline int setupTriangulation(ImplicitTriangulation *data) {
      triangulation_ = data;
      multiresTriangulation_.setTriangulation(triangulation_);
      return 0;
    }
    inline int setInputScalars(void *data) {
      inputScalars_ = data;
      return 0;
    }
    inline int setInputOffsets(void *data) {
      inputOffsets_ = data;
      return 0;
    }
    // inline int setOutputCTDiagram(void *data) {
    //   CTDiagram_ = data;
    //   return 0;
    // }

  protected:
    void sortPersistenceDiagram2(std::vector<PersistencePair> &diagram,
                                 const SimplexId *const offsets) const;
    // template <typename scalarType>
    // int sortPersistenceDiagram(std::vector<PersistencePair> &diagram,
    //                            const scalarType *const scalars,
    //                            const SimplexId *const offsets) const;

    /* PROGRESSIVE MODE DECLARATIONS */
  public:
    template <class scalarType>
    int executeCPProgressive(int computePersistenceDiagram,
                             const scalarType *inputScalars,
                             const SimplexId *inputOffsets);

    template <typename scalarType>
    int resumeProgressive(int computePersistenceDiagram,
                          const scalarType *scalars,
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

    template <typename scalarType>
    int computeProgressivePD(std::vector<PersistencePair> &CTDiagram,
                             const scalarType *scalars,
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
    // maximum link size in 3D
    static const size_t nLink_ = 27;
    using VLBoundaryType
      = std::array<std::vector<std::pair<SimplexId, SimplexId>>, nLink_>;

    template <typename ScalarType, typename OffsetType>
    void
      initCriticalPoints(std::vector<polarity> &isNew,
                         std::vector<std::vector<std::pair<polarity, polarity>>>
                           &vertexLinkPolarity,
                         std::vector<polarity> &toProcess,
                         std::vector<DynamicTree> &link,
                         std::vector<uint8_t> &vertexLink,
                         VLBoundaryType &vertexLinkByBoundaryType,
                         std::vector<char> &vertexTypes,
                         const ScalarType *const scalars,
                         const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
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
                         const ScalarType *const scalars,
                         const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
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
      const ScalarType *const scalars,
      const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
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
      const ScalarType *const scalars,
      const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
    void
      buildVertexLinkPolarity(const SimplexId vertexId,
                              std::vector<std::pair<polarity, polarity>> &vlp,
                              const ScalarType *const scalars,
                              const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
    void sortTriplets(std::vector<triplet> &triplets,
                      const ScalarType *const scalars,
                      const OffsetType *const offsets,
                      const bool splitTree) const;

    template <typename scalarType, typename offsetType>
    void tripletsToPersistencePairs(
      std::vector<PersistencePair> &pairs,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<triplet> &triplets,
      const scalarType *const scalars,
      const offsetType *const offsets,
      const bool splitTree) const;

    void buildVertexLinkByBoundary(const SimplexId vertexId,
                                   VLBoundaryType &vlbt) const;

    template <typename ScalarType, typename OffsetType>
    void initDynamicLink(const SimplexId &vertexId,
                         std::vector<std::pair<polarity, polarity>> &vlp,
                         uint8_t &vertexLink,
                         DynamicTree &link,
                         VLBoundaryType &vlbt,
                         const ScalarType *const scalars,
                         const OffsetType *const offsets) const;

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

    template <typename ScalarType, typename OffsetType>
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
      const ScalarType *const scalars,
      const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
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
                        const ScalarType *const scalars,
                        const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
    bool getMonotonyChangeByOldPointCP(
      const SimplexId vertexId,
      const std::vector<polarity> &isNew,
      std::vector<polarity> &toProcess,
      std::vector<polarity> &toReprocess,
      std::vector<std::pair<polarity, polarity>> &vlp,
      const ScalarType *const scalars,
      const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
    void updateLinkPolarity(const SimplexId vertexId,
                            std::vector<std::pair<polarity, polarity>> &vlp,
                            const ScalarType *const scalars,
                            const OffsetType *const offsets) const;

    template <typename ScalarType, typename OffsetType>
    ttk::SimplexId propageFromSaddles(
      const SimplexId vertexId,
      std::vector<Lock> &vertLock,
      std::vector<polarity> &toPropage,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<std::vector<SimplexId>> &saddleCC,
      std::vector<polarity> &isUpdated,
      std::vector<SimplexId> &globalExtremum,
      const ScalarType *const scalars,
      const OffsetType *const offsets,
      const bool splitTree) const;

    void getTripletsFromSaddles(
      const SimplexId vertexId,
      std::vector<triplet> &triplets,
      const std::vector<std::vector<SimplexId>> &vertexReps) const;

    template <typename ScalarType, typename OffsetType>
    void computePersistencePairsFromSaddles(
      std::vector<PersistencePair> &CTDiagram,
      const ScalarType *const scalars,
      const OffsetType *const offsets,
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
    void *inputScalars_{};
    void *inputOffsets_{};
    // void *CTDiagram_{};

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

// template <typename scalarType>
// int ttk::ProgressiveTopology::sortPersistenceDiagram(
//   std::vector<PersistencePair> &diagram,
//   const scalarType *const scalars,
//   const SimplexId *const offsets) const {
//   auto cmp
//     = [scalars, offsets](
//         const std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,
//                          ttk::CriticalType, scalarType, ttk::SimplexId> &a,
//         const std::tuple<ttk::SimplexId, ttk::CriticalType, ttk::SimplexId,
//                          ttk::CriticalType, scalarType, ttk::SimplexId> &b) {
//         const ttk::SimplexId idA = std::get<0>(a);
//         const ttk::SimplexId idB = std::get<0>(b);
//         const ttk::SimplexId va = offsets[idA];
//         const ttk::SimplexId vb = offsets[idB];
//         const scalarType sa = scalars[idA];
//         const scalarType sb = scalars[idB];

//         if(sa != sb)
//           return sa < sb;
//         else
//           return va < vb;
//       };

//   std::sort(diagram.begin(), diagram.end(), cmp);

//   return 0;
// }
template <typename scalarType>
int ttk::ProgressiveTopology::computeProgressivePD(
  std::vector<PersistencePair> &CTDiagram,
  const scalarType *scalars,
  const SimplexId *offsets) {
  int ret = -1;
  printMsg("Progressive Persistence Diagram computation");
  ret = executeCPProgressive(1, scalars, offsets);
  CTDiagram = std::move(CTDiagram_);
  CTDiagram_.clear();

  return ret;
}

template <typename scalarType>
int ttk::ProgressiveTopology::executeCPProgressive(
  int computePersistenceDiagram,
  const scalarType *scalars,
  const SimplexId *offsets) {

  printMsg(ttk::debug::Separator::L1);

  if(resumeProgressive_) {
    resumeProgressive(computePersistenceDiagram, scalars, offsets);
    return 0;
  }

  Timer timer;

  decimationLevel_ = startingDecimationLevel_;
  // multiresTriangulation_.setTriangulation(triangulation_);
  multiresTriangulation_.setDecimationLevel(0);
  const SimplexId vertexNumber = multiresTriangulation_.getVertexNumber();

#ifdef TTK_ENABLE_KAMIKAZE
  if(vertexNumber == 0) {
    this->printErr("No points in triangulation");
    return 1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  double tm_allocation = timer.getElapsedTime();

  // clean state (in case previous operation was a resume)
  clearResumableState();

  const auto dim = multiresTriangulation_.getDimensionality();
  const size_t maxNeigh = dim == 3 ? 14 : (dim == 2 ? 6 : 0);

  std::vector<std::vector<SimplexId>> saddleCCMin{}, saddleCCMax{};
  std::vector<std::vector<SimplexId>> vertexRepresentativesMin{},
    vertexRepresentativesMax{};
  std::vector<polarity> toPropageMin{}, toPropageMax{};
  std::vector<polarity> isUpToDateMin{}, isUpToDateMax{};
  std::vector<char> vertexTypes{};

  // std::vector<char> vertexTypes{};

  if(computePersistenceDiagram) {
    saddleCCMin.resize(vertexNumber);
    saddleCCMax.resize(vertexNumber);
    vertexRepresentativesMin.resize(vertexNumber);
    vertexRepresentativesMax.resize(vertexNumber);
    toPropageMin.resize(vertexNumber, 0);
    toPropageMax.resize(vertexNumber, 0);
    isUpToDateMin.resize(vertexNumber, 0);
    isUpToDateMax.resize(vertexNumber, 0);
  } else {
    vertexTypes.resize(vertexNumber, static_cast<char>(CriticalType::Regular));
  }

  // index in vertexLinkByBoundaryType
  std::vector<uint8_t> vertexLink(vertexNumber);
  VLBoundaryType vertexLinkByBoundaryType{};
  std::vector<DynamicTree> link(vertexNumber);
  std::vector<polarity> isNew(vertexNumber, 255);
  std::vector<std::vector<std::pair<polarity, polarity>>> vertexLinkPolarity(
    vertexNumber);
  std::vector<polarity> toProcess(vertexNumber, 0), toReprocess{};

  if(this->startingDecimationLevel_ > this->stoppingDecimationLevel_
     || this->isResumable_) {
    // only needed for progressive computation
    toReprocess.resize(vertexNumber, 0);
  }

  // lock vertex thread access for firstPropage
  std::vector<Lock> vertLockMin(vertexNumber), vertLockMax(vertexNumber);

  // pre-allocate memory
  if(preallocateMemory_) {
    double tm_prealloc = timer.getElapsedTime();
    printMsg("Pre-allocating data structures", 0, 0, threadNumber_,
             ttk::debug::LineMode::REPLACE);
    for(SimplexId i = 0; i < vertexNumber; ++i) {
      vertexLinkPolarity[i].reserve(maxNeigh);
      link[i].alloc(maxNeigh);
    }
    printMsg("Pre-allocating data structures", 1,
             timer.getElapsedTime() - tm_prealloc, threadNumber_);
  }

  tm_allocation = timer.getElapsedTime() - tm_allocation;
  printMsg("Total memory allocation", 1, tm_allocation, threadNumber_);

  // computation of implicit link
  std::vector<SimplexId> boundReps{};
  multiresTriangulation_.findBoundaryRepresentatives(boundReps);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(size_t i = 0; i < boundReps.size(); i++) {
    if(boundReps[i] != -1) {
      buildVertexLinkByBoundary(boundReps[i], vertexLinkByBoundaryType);
    }
  }

  if(debugLevel_ > 4) {
    std::cout << "boundary representatives : ";
    for(auto bb : boundReps) {
      std::cout << ", " << bb;
    }
    std::cout << std::endl;
  }

  printMsg(this->resolutionInfoString(), 0,
           timer.getElapsedTime() - tm_allocation, this->threadNumber_,
           ttk::debug::LineMode::REPLACE);
  multiresTriangulation_.setDecimationLevel(decimationLevel_);

  if(computePersistenceDiagram) {
    initSaddleSeeds(isNew, vertexLinkPolarity, toPropageMin, toPropageMax,
                    toProcess, link, vertexLink, vertexLinkByBoundaryType,
                    saddleCCMin, saddleCCMax, scalars, offsets);
    initPropagation(toPropageMin, toPropageMax, vertexRepresentativesMin,
                    vertexRepresentativesMax, saddleCCMin, saddleCCMax,
                    vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                    scalars, offsets);

    // compute pairs in non-progressive mode
    computePersistencePairsFromSaddles(
      CTDiagram_, scalars, offsets, vertexRepresentativesMin,
      vertexRepresentativesMax, toPropageMin, toPropageMax);
  } else {
    initCriticalPoints(isNew, vertexLinkPolarity, toProcess, link, vertexLink,
                       vertexLinkByBoundaryType, vertexTypes, scalars, offsets);
  }

  printMsg(this->resolutionInfoString(), 1,
           timer.getElapsedTime() - tm_allocation, this->threadNumber_);

  // skip subsequent propagations if time limit is exceeded
  stopComputationIf(timer.getElapsedTime() - tm_allocation
                    > 0.9 * this->timeLimit_);

  while(decimationLevel_ > stoppingDecimationLevel_) {
    Timer tmIter{};
    decimationLevel_--;
    multiresTriangulation_.setDecimationLevel(decimationLevel_);

    printMsg(this->resolutionInfoString(), 0,
             timer.getElapsedTime() - tm_allocation, this->threadNumber_,
             ttk::debug::LineMode::REPLACE);

    if(computePersistenceDiagram) {
      updateSaddleSeeds(isNew, vertexLinkPolarity, toPropageMin, toPropageMax,
                        toProcess, toReprocess, link, vertexLink,
                        vertexLinkByBoundaryType, saddleCCMin, saddleCCMax,
                        isUpToDateMin, isUpToDateMax, scalars, offsets);
      updatePropagation(toPropageMin, toPropageMax, vertexRepresentativesMin,
                        vertexRepresentativesMax, saddleCCMin, saddleCCMax,
                        vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                        scalars, offsets);
      computePersistencePairsFromSaddles(
        CTDiagram_, scalars, offsets, vertexRepresentativesMin,
        vertexRepresentativesMax, toPropageMin, toPropageMax);
    } else {
      updateCriticalPoints(isNew, vertexLinkPolarity, toProcess, toReprocess,
                           link, vertexLink, vertexLinkByBoundaryType,
                           vertexTypes, scalars, offsets);
    }

    const auto itDuration = tmIter.getElapsedTime();
    const auto nextItDuration
      = predictNextIterationDuration(itDuration, CTDiagram_.size() + 1);

    printMsg(this->resolutionInfoString(), 1,
             timer.getElapsedTime() - tm_allocation, this->threadNumber_);

    // skip subsequent propagations if time limit is exceeded
    stopComputationIf(timer.getElapsedTime() + nextItDuration - tm_allocation
                      > this->timeLimit_);

    if(debugLevel_ > 3) {
      std::cout << "current iteration lasted " << itDuration << "s"
                << std::endl;
      std::cout << "next iteration predicted to last at most " << nextItDuration
                << "s" << std::endl;
    }
  }

  // ADD GLOBAL MIN-MAX PAIR
  if(computePersistenceDiagram) {
    CTDiagram_.emplace_back(this->globalMin_, this->globalMax_, -1);
  }

  // store state for resuming computation
  if(this->isResumable_ and decimationLevel_ > 0) {
    if(computePersistenceDiagram) {
      this->vertexRepresentativesMax_ = std::move(vertexRepresentativesMax);
      this->vertexRepresentativesMin_ = std::move(vertexRepresentativesMin);
      this->saddleCCMin_ = std::move(saddleCCMin);
      this->saddleCCMax_ = std::move(saddleCCMax);
    }
    this->isNew_ = std::move(isNew);
    this->toProcess_ = std::move(toProcess);
    this->toReprocess_ = std::move(toReprocess);
    this->vertexLinkPolarity_ = std::move(vertexLinkPolarity);
    this->link_ = std::move(link);
    this->vertexLink_ = std::move(vertexLink);
    this->vertexLinkByBoundaryType_ = std::move(vertexLinkByBoundaryType);
    this->resumeProgressive_ = true;
  }

  // prepare outputs
  vertexTypes_ = std::move(vertexTypes);

  // finally sort the diagram
  sortPersistenceDiagram2(CTDiagram_, offsets);
  // sortPersistenceDiagram(CTDiagram, scalars, offsets);
  this->printMsg(
    "Total", 1.0, timer.getElapsedTime() - tm_allocation, this->threadNumber_);
  return 0;
}

template <typename scalarType>
int ttk::ProgressiveTopology::resumeProgressive(int computePersistenceDiagram,
                                                const scalarType *scalars,
                                                const SimplexId *offsets) {

  const auto vertexNumber = multiresTriangulation_.getVertexNumber();

  this->printMsg(
    "Resuming computation from resolution level "
    + std::to_string(multiresTriangulation_.DL_to_RL(decimationLevel_))
    + " to level "
    + std::to_string(
      multiresTriangulation_.DL_to_RL(stoppingDecimationLevel_)));

  // lock vertex thread access for firstPropage
  std::vector<Lock> vertLockMin(vertexNumber), vertLockMax(vertexNumber);
  // propagation markers
  std::vector<polarity> toPropageMin{}, toPropageMax{};
  std::vector<polarity> isUpToDateMin{}, isUpToDateMax{};

  if(computePersistenceDiagram) {
    toPropageMin.resize(vertexNumber, 0);
    toPropageMax.resize(vertexNumber, 0);
    isUpToDateMin.resize(vertexNumber, 0);
    isUpToDateMax.resize(vertexNumber, 0);
  }

  Timer timer;

  while(this->decimationLevel_ > this->stoppingDecimationLevel_) {
    Timer tmIter{};
    this->decimationLevel_--;
    multiresTriangulation_.setDecimationLevel(this->decimationLevel_);
    printMsg(this->resolutionInfoString(), 0, timer.getElapsedTime(),
             this->threadNumber_, ttk::debug::LineMode::REPLACE);
    if(computePersistenceDiagram) {
      updateSaddleSeeds(isNew_, vertexLinkPolarity_, toPropageMin, toPropageMax,
                        toProcess_, toReprocess_, link_, vertexLink_,
                        vertexLinkByBoundaryType_, saddleCCMin_, saddleCCMax_,
                        isUpToDateMin, isUpToDateMax, scalars, offsets);
      updatePropagation(toPropageMin, toPropageMax, vertexRepresentativesMin_,
                        vertexRepresentativesMax_, saddleCCMin_, saddleCCMax_,
                        vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                        scalars, offsets);
      computePersistencePairsFromSaddles(
        CTDiagram_, scalars, offsets, vertexRepresentativesMin_,
        vertexRepresentativesMax_, toPropageMin, toPropageMax);
    } else {
      updateCriticalPoints(
        isNew_, vertexLinkPolarity_, toProcess_, toReprocess_, link_,
        vertexLink_, vertexLinkByBoundaryType_, vertexTypes_, scalars, offsets);
    }

    const auto itDuration = tmIter.getElapsedTime();
    const auto nextItDuration
      = predictNextIterationDuration(itDuration, CTDiagram_.size() + 1);

    printMsg(this->resolutionInfoString(), 1, timer.getElapsedTime(),
             this->threadNumber_);
    //
    // skip subsequent propagations if time limit is exceeded
    stopComputationIf(timer.getElapsedTime() + nextItDuration
                      > this->timeLimit_);
  }

  // ADD GLOBAL MIN-MAX PAIR
  if(computePersistenceDiagram) {
    CTDiagram_.emplace_back(this->globalMin_, this->globalMax_, -1);
  }
  // finally sort the diagram
  sortPersistenceDiagram2(CTDiagram_, offsets);
  // sortPersistenceDiagram(CTDiagram, scalars, offsets);
  this->printMsg("Total", 1.0, timer.getElapsedTime(), this->threadNumber_);

  // clean state (we don't need it anymore)
  if(this->decimationLevel_ == 0) {
    clearResumableState();
  }

  return 0;
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::computePersistencePairsFromSaddles(
  std::vector<PersistencePair> &CTDiagram,
  const ScalarType *const scalars,
  const OffsetType *const offsets,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
  const std::vector<polarity> &toPropageMin,
  const std::vector<polarity> &toPropageMax) const {

  Timer timer{};
  std::vector<triplet> tripletsMax{}, tripletsMin{};
  const SimplexId nbDecVert = multiresTriangulation_.getDecimatedVertexNumber();

  for(SimplexId localId = 0; localId < nbDecVert; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(toPropageMin[globalId]) {
      getTripletsFromSaddles(globalId, tripletsMin, vertexRepresentativesMin);
    }
    if(toPropageMax[globalId]) {
      getTripletsFromSaddles(globalId, tripletsMax, vertexRepresentativesMax);
    }
  }

  if(debugLevel_ > 3) {
    std::cout << "TRIPLETS " << timer.getElapsedTime() << std::endl;
  }
  double tm_pairs = timer.getElapsedTime();

  sortTriplets(tripletsMax, scalars, offsets, true);
  sortTriplets(tripletsMin, scalars, offsets, false);

  const auto tm_sort = timer.getElapsedTime();
  if(debugLevel_ > 3) {
    std::cout << "TRIPLETS SORT " << tm_sort - tm_pairs << std::endl;
  }

  typename std::remove_reference<decltype(CTDiagram)>::type CTDiagramMin{},
    CTDiagramMax{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    tripletsToPersistencePairs(CTDiagramMin, vertexRepresentativesMax,
                               tripletsMax, scalars, offsets, true);
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    tripletsToPersistencePairs(CTDiagramMax, vertexRepresentativesMin,
                               tripletsMin, scalars, offsets, false);
  }
  CTDiagram = std::move(CTDiagramMin);
  CTDiagram.insert(CTDiagram.end(), CTDiagramMax.begin(), CTDiagramMax.end());

  if(debugLevel_ > 3) {
    std::cout << "PAIRS " << timer.getElapsedTime() - tm_sort << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::sortTriplets(std::vector<triplet> &triplets,
                                            const ScalarType *const scalars,
                                            const OffsetType *const offsets,
                                            const bool splitTree) const {
  if(triplets.empty())
    return;

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
    // return (scalars[a] < scalars[b])
    //        || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
  };

  // Sorting step
  const auto cmp = [=](const triplet &t1, const triplet &t2) {
    const SimplexId s1 = std::get<0>(t1);
    const SimplexId s2 = std::get<0>(t2);
    const SimplexId m1 = std::get<2>(t1);
    const SimplexId m2 = std::get<2>(t2);
    if(s1 != s2)
      return lt(s1, s2) != splitTree;
    else // s1 == s2
      return lt(m1, m2) == splitTree;
  };

  PARALLEL_SORT(triplets.begin(), triplets.end(), cmp);
}

template <typename scalarType, typename offsetType>
void ttk::ProgressiveTopology::tripletsToPersistencePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<triplet> &triplets,
  const scalarType *const scalars,
  const offsetType *const offsets,
  const bool splitTree) const {

  Timer tm;
  if(triplets.empty())
    return;
  size_t numberOfPairs = 0;

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
    // return (scalars[a] < scalars[b])
    //        || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
  };

  const auto getRep = [&](SimplexId v) -> SimplexId {
    auto r = vertexRepresentatives[v][0];
    while(r != v) {
      v = r;
      r = vertexRepresentatives[v][0];
    }
    return r;
  };

  for(const auto &t : triplets) {
    SimplexId r1 = getRep(std::get<1>(t));
    SimplexId r2 = getRep(std::get<2>(t));
    if(r1 != r2) {
      SimplexId s = std::get<0>(t);
      numberOfPairs++;

      // Add pair
      if(splitTree) {
        // r1 = min(r1, r2), r2 = max(r1, r2)
        if(lt(r2, r1)) {
          std::swap(r1, r2);
        }
        // pair saddle-max: s -> min(r1, r2);
        pairs.emplace_back(s, r1, 2);

      } else {
        // r1 = max(r1, r2), r2 = min(r1, r2)
        if(lt(r1, r2)) {
          std::swap(r1, r2);
        }
        // pair min-saddle: max(r1, r2) -> s;
        pairs.emplace_back(r1, s, 0);
      }

      vertexRepresentatives[std::get<1>(t)][0] = r2;
      vertexRepresentatives[r1][0] = r2;
    }
  }

  if(debugLevel_ > 3) {
    std::string prefix = splitTree ? "[sad-max]" : "[min-sad]";
    std::cout << prefix << "  found all pairs in " << tm.getElapsedTime()
              << " s." << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::buildVertexLinkPolarity(
  const SimplexId vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vlp.resize(neighborNumber);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);
    const bool lower0 = offsets[neighborId0] < offsets[vertexId];
    // const bool lower0 = (scalars[neighborId0] < scalars[vertexId])
    //                     || (scalars[neighborId0] == scalars[vertexId]
    //                         && offsets[neighborId0] < offsets[vertexId]);
    const polarity isUpper0 = static_cast<polarity>(!lower0) * 255;
    vlp[i] = std::make_pair(isUpper0, 0);
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::initDynamicLink(
  const SimplexId &vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  uint8_t &vertexLink,
  DynamicTree &link,
  VLBoundaryType &vlbt,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  if(vlp.empty()) {
    buildVertexLinkPolarity(vertexId, vlp, scalars, offsets);
  }

  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  link.alloc(neighborNumber);

  // associate vertex link boundary
  vertexLink = multiresTriangulation_.getVertexBoundaryIndex(vertexId);

  // update the link polarity for old points that are processed for
  // the first time
  const auto &vl = vlbt[vertexLink];
  for(size_t edgeId = 0; edgeId < vl.size(); edgeId++) {
    const SimplexId n0 = vl[edgeId].first;
    const SimplexId n1 = vl[edgeId].second;
    if(vlp[n0].first == vlp[n1].first) {
      // the smallest id (n0) becomes the parent of n1
      link.insertEdge(n1, n0);
    }
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::updateCriticalPoints(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<char> &vertexTypes,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer tm;
  const auto nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // find breaking edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId localId = 0; localId < nDecVerts; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(isNew[globalId]) {
      if(decimationLevel_ > stoppingDecimationLevel_ || isResumable_) {
        buildVertexLinkPolarity(
          globalId, vertexLinkPolarity[globalId], scalars, offsets);
      }
    } else {
      getMonotonyChangeByOldPointCP(globalId, isNew, toProcess, toReprocess,
                                    vertexLinkPolarity[globalId], scalars,
                                    offsets);
    }
  }

  if(debugLevel_ > 3) {
    std::cout << "MONOTONY " << tm.getElapsedTime() << " s." << std::endl;
  }

  double t_critical = tm.getElapsedTime();
  // second Loop  process or reprocess
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);

    if(isNew[globalId]) { // new point
      if(toProcess[globalId]) {
        initDynamicLink(globalId, vertexLinkPolarity[globalId],
                        vertexLink[globalId], link[globalId],
                        vertexLinkByBoundaryType, scalars, offsets);
        vertexTypes[globalId] = getCriticalTypeFromLink(
          globalId, vertexLinkPolarity[globalId], link[globalId]);
      }
      isNew[globalId] = false;

    } else { // old point
      if(toReprocess[globalId]) {
        if(toProcess[globalId]) { // was already processed : need to reprocess
          updateDynamicLink(link[globalId], vertexLinkPolarity[globalId],
                            vertexLinkByBoundaryType[vertexLink[globalId]]);
        } else { // first processing
          updateLinkPolarity(
            globalId, vertexLinkPolarity[globalId], scalars, offsets);
          initDynamicLink(globalId, vertexLinkPolarity[globalId],
                          vertexLink[globalId], link[globalId],
                          vertexLinkByBoundaryType, scalars, offsets);
          toProcess[globalId] = 255; // mark as processed
        }
        vertexTypes[globalId] = getCriticalTypeFromLink(
          globalId, vertexLinkPolarity[globalId], link[globalId]);
        toReprocess[globalId] = 0;
      }
    }

  } // end for openmp

  if(debugLevel_ > 3) {
    std::cout << "CRITICAL POINTS UPDATE " << tm.getElapsedTime() - t_critical
              << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::updateSaddleSeeds(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
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
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer tm;
  const auto nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // find breaking edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId localId = 0; localId < nDecVerts; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(isNew[globalId]) {
      if(decimationLevel_ > stoppingDecimationLevel_ || isResumable_) {
        buildVertexLinkPolarity(
          globalId, vertexLinkPolarity[globalId], scalars, offsets);
      }
    } else {
      getMonotonyChangeByOldPointCP(globalId, isNew, toProcess, toReprocess,
                                    vertexLinkPolarity[globalId], scalars,
                                    offsets);
    }
  }

  if(debugLevel_ > 3) {
    std::cout << "MONOTONY " << tm.getElapsedTime() << " s." << std::endl;
  }

  double t_critical = tm.getElapsedTime();
  // second Loop  process or reprocess
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);

    if(isNew[globalId]) { // new point
      if(toProcess[globalId]) {
        initDynamicLink(globalId, vertexLinkPolarity[globalId],
                        vertexLink[globalId], link[globalId],
                        vertexLinkByBoundaryType, scalars, offsets);
        getValencesFromLink(globalId, vertexLinkPolarity[globalId],
                            link[globalId], toPropageMin, toPropageMax,
                            saddleCCMin, saddleCCMax);
      }
      isNew[globalId] = false;

    } else { // old point
      if(toReprocess[globalId]) {
        if(toProcess[globalId]) { // was already processed : need to reprocess
          updateDynamicLink(link[globalId], vertexLinkPolarity[globalId],
                            vertexLinkByBoundaryType[vertexLink[globalId]]);
        } else { // first processing
          updateLinkPolarity(
            globalId, vertexLinkPolarity[globalId], scalars, offsets);
          initDynamicLink(globalId, vertexLinkPolarity[globalId],
                          vertexLink[globalId], link[globalId],
                          vertexLinkByBoundaryType, scalars, offsets);
          toProcess[globalId] = 255; // mark as processed
        }
        getValencesFromLink(globalId, vertexLinkPolarity[globalId],
                            link[globalId], toPropageMin, toPropageMax,
                            saddleCCMin, saddleCCMax);
        toReprocess[globalId] = 0;
      }
    }

    // reset updated flag
    isUpdatedMin[globalId] = 0;
    isUpdatedMax[globalId] = 0;
  } // end for openmp

  if(debugLevel_ > 3) {
    std::cout << "CRITICAL POINTS UPDATE " << tm.getElapsedTime() - t_critical
              << std::endl;
  }
}
template <typename ScalarType, typename OffsetType>

bool ttk::ProgressiveTopology::getMonotonyChangeByOldPointCP(
  const SimplexId vertexId,
  const std::vector<polarity> &isNew,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  bool hasMonotonyChanged = false;
  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);

    // check for monotony changes
    const bool lower = offsets[neighborId] < offsets[vertexId];
    // const bool lower = (scalars[neighborId] < scalars[vertexId])
    //                    || (scalars[neighborId] == scalars[vertexId]
    //                        && offsets[neighborId] < offsets[vertexId]);
    const polarity isUpper = lower ? 0 : 255;
    const polarity isUpperOld = vlp[i].first;

    if(isUpper != isUpperOld) { // change of monotony
      hasMonotonyChanged = true;

      toReprocess[vertexId] = 255;
      toProcess[neighborId] = 255;
      const SimplexId neighborNumberNew
        = multiresTriangulation_.getVertexNeighborNumber(neighborId);
      for(SimplexId j = 0; j < neighborNumberNew; j++) {
        SimplexId neighborIdNew = -1;
        multiresTriangulation_.getVertexNeighbor(neighborId, j, neighborIdNew);
        if(isNew[neighborIdNew])
          toProcess[neighborIdNew] = 255;
      }

      vlp[i].second = 255;
    }
  }
  return hasMonotonyChanged;
}

template <typename ScalarType, typename OffsetType>
ttk::SimplexId ttk::ProgressiveTopology::propageFromSaddles(
  const SimplexId vertexId,
  std::vector<Lock> &vertLock,
  std::vector<polarity> &toPropage,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<std::vector<SimplexId>> &saddleCC,
  std::vector<polarity> &isUpdated,
  std::vector<SimplexId> &globalExtremum,
  const ScalarType *const scalars,
  const OffsetType *const offsets,
  const bool splitTree) const {

  auto &toProp = toPropage[vertexId];
  auto &reps = vertexRepresentatives[vertexId];
  auto &updated = isUpdated[vertexId];

  const auto gt = [=](const SimplexId a, const SimplexId b) {
    return (offsets[a] > offsets[b]) == splitTree;
    // return ((scalars[a] > scalars[b])
    //         || (scalars[a] == scalars[b] && offsets[a] > offsets[b]))
    //        == splitTree;
  };

  if(updated) {
    return reps[0];
  }

  if(this->threadNumber_ > 1) {
    vertLock[vertexId].lock();
  }

  if(toProp) { // SADDLE POINT
    const auto &CC = saddleCC[vertexId];
    reps.clear();
    reps.reserve(CC.size());
    for(size_t r = 0; r < CC.size(); r++) {
      SimplexId neighborId = -1;
      SimplexId localId = CC[r];
      multiresTriangulation_.getVertexNeighbor(vertexId, localId, neighborId);
      SimplexId ret = propageFromSaddles(
        neighborId, vertLock, toPropage, vertexRepresentatives, saddleCC,
        isUpdated, globalExtremum, scalars, offsets, splitTree);
      reps.emplace_back(ret);
    }

    if(reps.size() > 1) {
      // sort & remove duplicate elements
      std::sort(reps.begin(), reps.end(), gt);
      const auto last = std::unique(reps.begin(), reps.end());
      reps.erase(last, reps.end());
    }

    updated = 255;
    if(this->threadNumber_ > 1) {
      vertLock[vertexId].unlock();
    }

    return reps[0];

  } else {

    SimplexId ret = vertexId;
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    SimplexId maxNeighbor = vertexId;
    for(SimplexId i = 0; i < neighborNumber; i++) {
      SimplexId neighborId = -1;
      multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
      if(gt(neighborId, maxNeighbor)) {
        maxNeighbor = neighborId;
      }
    }
    if(maxNeighbor != vertexId) { // not an extremum
      ret = propageFromSaddles(maxNeighbor, vertLock, toPropage,
                               vertexRepresentatives, saddleCC, isUpdated,
                               globalExtremum, scalars, offsets, splitTree);
    } else {
#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
      if(gt(vertexId, globalExtremum[tid])) {
        globalExtremum[tid] = vertexId;
      }
    }
    reps.resize(1);
    reps[0] = ret;
    updated = 255;
    if(this->threadNumber_ > 1) {
      vertLock[vertexId].unlock();
    }
    return ret;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::initCriticalPoints(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<char> &vertexTypes,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // computes the critical types of all points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    buildVertexLinkPolarity(
      globalId, vertexLinkPolarity[globalId], scalars, offsets);
    initDynamicLink(globalId, vertexLinkPolarity[globalId],
                    vertexLink[globalId], link[globalId],
                    vertexLinkByBoundaryType, scalars, offsets);
    vertexTypes[globalId] = getCriticalTypeFromLink(
      globalId, vertexLinkPolarity[globalId], link[globalId]);
    toProcess[globalId] = 255;
    isNew[globalId] = 0;
  }

  if(debugLevel_ > 3) {
    std::cout << "initial critical types in " << timer.getElapsedTime() << " s."
              << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::initSaddleSeeds(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<polarity> &toProcess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // computes the critical types of all points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    buildVertexLinkPolarity(
      globalId, vertexLinkPolarity[globalId], scalars, offsets);
    initDynamicLink(globalId, vertexLinkPolarity[globalId],
                    vertexLink[globalId], link[globalId],
                    vertexLinkByBoundaryType, scalars, offsets);
    getValencesFromLink(globalId, vertexLinkPolarity[globalId], link[globalId],
                        toPropageMin, toPropageMax, saddleCCMin, saddleCCMax);
    toProcess[globalId] = 255;
    isNew[globalId] = 0;
  }

  if(debugLevel_ > 3) {
    std::cout << "initial critical types in " << timer.getElapsedTime() << " s."
              << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::initPropagation(
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
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  std::vector<SimplexId> globalMaxThr(threadNumber_, 0);
  std::vector<SimplexId> globalMinThr(threadNumber_, 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    if(toPropageMin[v]) {
      propageFromSaddles(v, vertLockMin, toPropageMin, vertexRepresentativesMin,
                         saddleCCMin, isUpdatedMin, globalMinThr, scalars,
                         offsets, false);
    }
    if(toPropageMax[v]) {
      propageFromSaddles(v, vertLockMax, toPropageMax, vertexRepresentativesMax,
                         saddleCCMax, isUpdatedMax, globalMaxThr, scalars,
                         offsets, true);
    }
  }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
    // return (scalars[a] < scalars[b])
    //        || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
  };

  globalMin_ = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
  globalMax_ = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);

  if(debugLevel_ > 3) {
    std::cout << "FIRSTPROPAGATION " << timer.getElapsedTime() << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::updatePropagation(
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
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  Timer tm{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  if(debugLevel_ > 5) {
    const auto pred = [](const polarity a) { return a > 0; };
    const auto numberOfCandidatesToPropageMax
      = std::count_if(toPropageMax.begin(), toPropageMax.end(), pred);
    std::cout << " sad-max we have " << numberOfCandidatesToPropageMax
              << " vertices to propage from outta " << nDecVerts << std::endl;
    const auto numberOfCandidatesToPropageMin
      = std::count_if(toPropageMin.begin(), toPropageMin.end(), pred);
    std::cout << " min-sad we have " << numberOfCandidatesToPropageMin
              << " vertices to propage from outta " << nDecVerts << std::endl;
  }

  std::vector<SimplexId> globalMaxThr(threadNumber_, 0);
  std::vector<SimplexId> globalMinThr(threadNumber_, 0);

  // propage along split tree
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    if(toPropageMin[v]) {
      propageFromSaddles(v, vertLockMin, toPropageMin, vertexRepresentativesMin,
                         saddleCCMin, isUpdatedMin, globalMinThr, scalars,
                         offsets, false);
    }
    if(toPropageMax[v]) {
      propageFromSaddles(v, vertLockMax, toPropageMax, vertexRepresentativesMax,
                         saddleCCMax, isUpdatedMax, globalMaxThr, scalars,
                         offsets, true);
    }
  }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return offsets[a] < offsets[b];
    // return (scalars[a] < scalars[b])
    //        || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
  };

  globalMin_ = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
  globalMax_ = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);

  if(debugLevel_ > 3) {
    std::cout << "PROPAGATION UPDATE " << tm.getElapsedTime() << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ProgressiveTopology::updateLinkPolarity(
  const SimplexId vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const ScalarType *const scalars,
  const OffsetType *const offsets) const {

  for(size_t i = 0; i < vlp.size(); i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
    // const bool lower = (scalars[neighborId] < scalars[vertexId])
    //                    || (scalars[neighborId] == scalars[vertexId]
    //                        && offsets[neighborId] < offsets[vertexId]);
    const bool lower = offsets[neighborId] < offsets[vertexId];
    const polarity isUpper = lower ? 0 : 255;
    vlp[i] = std::make_pair(isUpper, 0);
  }
}
