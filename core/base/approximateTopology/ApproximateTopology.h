/// \ingroup base
/// \class ttk::ApproximateTopology
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \date 2021.
///
/// \brief TTK processing package for progressive Topological Data Analysis
///
/// This package introduces an approximation algorithm for the
/// computation of the extremum-saddle persistence diagram of a scalar field.
/// The approximation comes with a user-controlled error on the Bottleneck
/// distance to the exact diagram.
///
/// \b Related \b publication \n
/// "Fast Approximation of Persistence Diagrams with Guarantees" \n
/// Jules Vidal, Julien Tierny\n
/// IEEE Symposium on Large Data Visualization and Analysis (LDAV), 2021
///
/// \sa PersistenceDiagram

#pragma once

// base code includes
#include <DynamicTree.h>
#include <ImplicitTriangulation.h>
#include <MultiresTopology.h>
#include <MultiresTriangulation.h>

#include <limits>
#include <numeric>
#include <tuple>

namespace ttk {

  using triplet = std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>;
  using polarity = unsigned char;

  class ApproximateTopology : public MultiresTopology {

  public:
    ApproximateTopology() {
      this->setDebugMsgPrefix("ApproximateTopology");
    }
    void setEpsilon(double data) {
      epsilon_ = data;
    }
    void setDelta(double data) {
      delta_ = data;
    }
    template <typename scalarType>
    int computeApproximatePD(std::vector<PersistencePair> &CTDiagram,
                             const scalarType *scalars,
                             scalarType *const fakeScalars,
                             SimplexId *const outputOffsets,
                             int *const outputMonotonyOffsets);

    template <typename scalarType>
    int executeApproximateTopology(const scalarType *scalars,
                                   scalarType *fakeScalars,
                                   SimplexId *outputOffsets,
                                   int *outputMonotonyOffsets);

    template <typename scalarType, typename offsetType>
    void
      initGlobalPolarity(std::vector<polarity> &isNew,
                         std::vector<std::vector<std::pair<polarity, polarity>>>
                           &vertexLinkPolarity,
                         std::vector<polarity> &toProcess,
                         const scalarType *fakeScalars,
                         const offsetType *const offsets,
                         const int *const monotonyOffsets) const;

    template <typename scalarType, typename offsetType>
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
      const scalarType *fakeScalars,
      const offsetType *const offsetField,
      const int *const monotonyOffsets);

    template <typename scalarType, typename offsetType>
    void buildVertexLinkPolarityApproximate(
      const SimplexId vertexId,
      std::vector<std::pair<polarity, polarity>> &vlp,
      const scalarType *fakeScalars,
      const offsetType *const offsetField,
      const int *const monotonyOffsets) const;

    template <typename ScalarType, typename OffsetType>
    void sortTripletsApproximate(std::vector<triplet> &triplets,
                                 const ScalarType *const scalars,
                                 const ScalarType *const fakeScalars,
                                 const OffsetType *const offsets,
                                 const int *const monotonyOffsets,
                                 const bool splitTree) const;

    template <typename scalarType, typename offsetType>
    void getCriticalTypeApproximate(
      const SimplexId &vertexId,
      std::vector<std::pair<polarity, polarity>> &vlp,
      uint8_t &vertexLink,
      DynamicTree &link,
      VLBoundaryType &vlbt,
      const scalarType *fakeScalars,
      const offsetType *const offsets,
      const int *const monotonyOffsets) const;

    template <typename ScalarType, typename offsetType>
    void computeCriticalPoints(
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
      ScalarType *fakeScalars,
      const offsetType *const offsets,
      int *monotonyOffsets);

    template <typename scalarType, typename offsetType>
    ttk::SimplexId propageFromSaddles(
      const SimplexId vertexId,
      std::vector<Lock> &vertLock,
      std::vector<polarity> &toPropage,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<std::vector<SimplexId>> &saddleCC,
      std::vector<polarity> &isUpdated,
      std::vector<SimplexId> &globalExtremum,
      const bool splitTree,
      const scalarType *fakeScalars,
      const offsetType *const offsetField,
      const int *const monotonyOffsets) const;

    template <typename ScalarType, typename offsetType>
    void computePersistencePairsFromSaddles(
      std::vector<PersistencePair> &CTDiagram,
      const ScalarType *const fakeScalars,
      const offsetType *const offsets,
      const int *const monotonyOffsets,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
      std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
      const std::vector<polarity> &toPropageMin,
      const std::vector<polarity> &toPropageMax) const;

    template <typename scalarType, typename offsetType>
    void sortVertices(const SimplexId vertexNumber,
                      std::vector<SimplexId> &sortedVertices,
                      SimplexId *vertsOrder,
                      const scalarType *const fakeScalars,
                      const offsetType *const offsetField,
                      const int *const monotonyOffsets);

    template <typename scalarType>
    int sortPersistenceDiagramApproximate(std::vector<PersistencePair> &diagram,
                                          scalarType *fakeScalars,
                                          const SimplexId *const offsets,
                                          int *monotonyOffsets) const;

  protected:
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

    template <typename scalarType, typename offsetType>
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
      const scalarType *fakeScalars,
      const offsetType *const offsets,
      const int *const monotonyOffsets) const;

    template <typename scalarType, typename offsetType>
    void
      buildVertexLinkPolarity(const SimplexId vertexId,
                              std::vector<std::pair<polarity, polarity>> &vlp,
                              const scalarType *fakeScalars,
                              const offsetType *const offsets,
                              const int *const monotonyOffsets) const;

    template <typename ScalarType, typename OffsetType>
    void sortTriplets(std::vector<triplet> &triplets,
                      const ScalarType *const fakeScalars,
                      const OffsetType *const offsets,
                      const int *const monotonyOffsets,
                      const bool splitTree) const;

    template <typename scalarType, typename offsetType>
    void tripletsToPersistencePairs(
      std::vector<PersistencePair> &pairs,
      std::vector<std::vector<SimplexId>> &vertexRepresentatives,
      std::vector<triplet> &triplets,
      const scalarType *const fakeScalars,
      const offsetType *const offsets,
      const int *const monotonyOffsets,
      const bool splitTree) const;

    template <typename scalarType, typename offsetType>
    int getMonotonyChangeByOldPointCPApproximate(
      const SimplexId vertexId,
      double eps,
      const std::vector<polarity> &isNew,
      std::vector<polarity> &toProcess,
      std::vector<polarity> &toReprocess,
      std::vector<std::pair<polarity, polarity>> &vlp,
      scalarType *fakeScalars,
      const offsetType *const offsets,
      int *monotonyOffsets) const;

    template <typename ScalarType, typename offsetType>
    int updateGlobalPolarity(
      double eps,
      std::vector<polarity> &isNew,
      std::vector<std::vector<std::pair<polarity, polarity>>>
        &vertexLinkPolarity,
      std::vector<polarity> &toProcess,
      std::vector<polarity> &toReprocess,
      ScalarType *fakeScalars,
      const offsetType *const offsets,
      int *monotonyOffsets) const;

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

    template <typename scalarType, typename offsetType>
    bool printPolarity(std::vector<polarity> &isNew,
                       const SimplexId vertexId,
                       std::vector<std::vector<std::pair<polarity, polarity>>>
                         &vertexLinkPolarity,
                       const scalarType *scalars,
                       const scalarType *fakeScalars,
                       const offsetType *const offsets,
                       const int *const monotonyOffsets,
                       bool verbose = false);

    double epsilon_{};
    double delta_{};
  };
} // namespace ttk

template <typename scalarType>
int ttk::ApproximateTopology::executeApproximateTopology(
  const scalarType *ttkNotUsed(scalars),
  scalarType *fakeScalars,
  SimplexId *outputOffsets,
  int *outputMonotonyOffsets) {

  Timer timer;

  SimplexId *const vertsOrder = static_cast<SimplexId *>(outputOffsets);

  decimationLevel_ = startingDecimationLevel_;
  multiresTriangulation_.setTriangulation(triangulation_);
  const SimplexId vertexNumber = multiresTriangulation_.getVertexNumber();

  int *monotonyOffsets = outputMonotonyOffsets;

#ifdef TTK_ENABLE_KAMIKAZE
  if(vertexNumber == 0) {
    this->printErr("No points in triangulation");
    return 1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  double tm_allocation = timer.getElapsedTime();

  const auto dim = multiresTriangulation_.getDimensionality();
  const size_t maxNeigh = dim == 3 ? 14 : (dim == 2 ? 6 : 0);

  std::vector<std::vector<SimplexId>> saddleCCMin(vertexNumber),
    saddleCCMax(vertexNumber);
  std::vector<std::vector<SimplexId>> vertexRepresentativesMin(vertexNumber),
    vertexRepresentativesMax(vertexNumber);

  std::vector<std::vector<std::pair<polarity, polarity>>> vertexLinkPolarity(
    vertexNumber);

  std::vector<polarity> isNew(vertexNumber, 255);
  std::vector<polarity> toPropageMin(vertexNumber, 0),
    toPropageMax(vertexNumber, 0);
  std::vector<polarity> isUpToDateMin(vertexNumber, 0),
    isUpToDateMax(vertexNumber, 0);

  // index in vertexLinkByBoundaryType
  std::vector<uint8_t> vertexLink(vertexNumber);
  VLBoundaryType vertexLinkByBoundaryType{};
  std::vector<DynamicTree> link(vertexNumber);
  std::vector<polarity> toProcess(vertexNumber, 0), toReprocess{};

  std::vector<SimplexId> offsetsVec(vertexNumber);
  std::iota(offsetsVec.begin(), offsetsVec.end(), 0);
  const SimplexId *const offsets = offsetsVec.data();

  if(this->startingDecimationLevel_ > this->stoppingDecimationLevel_) {
    // only needed for progressive computation
    toReprocess.resize(vertexNumber, 0);
  }

  std::vector<Lock> vertLockMin(vertexNumber), vertLockMax(vertexNumber);

  if(preallocateMemory_) {
    double tm_prealloc = timer.getElapsedTime();
    printMsg("Pre-allocating data structures", 0, 0, threadNumber_,
             ttk::debug::LineMode::REPLACE);
    for(SimplexId i = 0; i < vertexNumber; ++i) {
      // for(int d = 0; d < vertexLinkPolarity.size(); d++) {
      //   vertexLinkPolarity[d][i].reserve(maxNeigh);
      // }
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

  // if(debugLevel_ > 4) {
  //   std::cout << "boundary representatives : ";
  //   for(auto bb : boundReps) {
  //     std::cout << ", " << bb;
  //   }
  //   std::cout << std::endl;
  // }

  multiresTriangulation_.setDecimationLevel(decimationLevel_);

  initGlobalPolarity(isNew, vertexLinkPolarity, toProcess, fakeScalars, offsets,
                     monotonyOffsets);

  double delta = epsilon_ * delta_;

  while(decimationLevel_ > stoppingDecimationLevel_) {
    Timer tmIter{};
    decimationLevel_--;
    multiresTriangulation_.setDecimationLevel(decimationLevel_);
    // double progress = (double)(startingDecimationLevel_ - decimationLevel_)
    //                   / (startingDecimationLevel_);
    // this->printMsg("decimation level: " + std::to_string(decimationLevel_),
    //                progress, timer.getElapsedTime() - tm_allocation,
    //                threadNumber_);

    int ret = updateGlobalPolarity(delta, isNew, vertexLinkPolarity, toProcess,
                                   toReprocess, fakeScalars, offsets,
                                   monotonyOffsets);
    if(ret == -1) {
      std::cout << "Found ERROR - aborting" << std::endl;
      return -1;
    }
  } // end while

  computeCriticalPoints(vertexLinkPolarity, toPropageMin, toPropageMax,
                        toProcess, link, vertexLink, vertexLinkByBoundaryType,
                        saddleCCMin, saddleCCMax, fakeScalars, offsets,
                        monotonyOffsets);

  updatePropagation(toPropageMin, toPropageMax, vertexRepresentativesMin,
                    vertexRepresentativesMax, saddleCCMin, saddleCCMax,
                    vertLockMin, vertLockMax, isUpToDateMin, isUpToDateMax,
                    fakeScalars, offsets, monotonyOffsets);

  computePersistencePairsFromSaddles(
    CTDiagram_, fakeScalars, offsets, monotonyOffsets, vertexRepresentativesMin,
    vertexRepresentativesMax, toPropageMin, toPropageMax);
  // ADD GLOBAL MIN-MAX PAIR
  CTDiagram_.emplace_back(this->globalMin_, this->globalMax_, -1);
  // fakeScalars[this->globalMax_] - fakeScalars[this->globalMin_], -1);

  // std::cout << "# epsilon " << epsilon_ << " in "
  //           << timer.getElapsedTime() - tm_allocation << " s." << std::endl;

  this->printMsg("Complete", 1.0, timer.getElapsedTime() - tm_allocation,
                 this->threadNumber_);

  // finally sort the diagram
  // std::cout << "SORTING" << std::endl;
  sortPersistenceDiagramApproximate(
    CTDiagram_, fakeScalars, offsets, monotonyOffsets);
  // std::cout << "Final epsilon " << epsilon_ << std::endl;

  // sort vertices to generate correct output offset order
  std::vector<SimplexId> sortedVertices{};
  sortVertices(vertexNumber, sortedVertices, vertsOrder, fakeScalars, offsets,
               monotonyOffsets);
  return 0;
}

template <typename ScalarType, typename OffsetType>
void ttk::ApproximateTopology::computePersistencePairsFromSaddles(
  std::vector<PersistencePair> &CTDiagram,
  const ScalarType *const fakeScalars,
  const OffsetType *const offsets,
  const int *const monotonyOffsets,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMin,
  std::vector<std::vector<SimplexId>> &vertexRepresentativesMax,
  const std::vector<polarity> &toPropageMin,
  const std::vector<polarity> &toPropageMax) const {

  Timer timer{};
  // CTDiagram.clear();
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

  sortTriplets(tripletsMax, fakeScalars, offsets, monotonyOffsets, true);
  sortTriplets(tripletsMin, fakeScalars, offsets, monotonyOffsets, false);

  const auto tm_sort = timer.getElapsedTime();

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
                               tripletsMax, fakeScalars, offsets,
                               monotonyOffsets, true);
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif // TTK_ENABLE_OPENMP
    tripletsToPersistencePairs(CTDiagramMax, vertexRepresentativesMin,
                               tripletsMin, fakeScalars, offsets,
                               monotonyOffsets, false);
  }
  CTDiagram = std::move(CTDiagramMin);
  CTDiagram.insert(CTDiagram.end(), CTDiagramMax.begin(), CTDiagramMax.end());

  if(debugLevel_ > 3) {
    std::cout << "PAIRS " << timer.getElapsedTime() - tm_sort << std::endl;
  }
}

template <typename ScalarType, typename OffsetType>
void ttk::ApproximateTopology::sortTriplets(std::vector<triplet> &triplets,
                                            const ScalarType *const fakeScalars,
                                            const OffsetType *const offsets,
                                            const int *const monotonyOffsets,
                                            const bool splitTree) const {
  if(triplets.empty())
    return;

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return ((fakeScalars[a] < fakeScalars[b])
            || (fakeScalars[a] == fakeScalars[b]
                && ((monotonyOffsets[a] < monotonyOffsets[b])
                    || (monotonyOffsets[a] == monotonyOffsets[b]
                        && offsets[a] < offsets[b]))));
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

  TTK_PSORT(this->threadNumber_, triplets.begin(), triplets.end(), cmp);
}

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::tripletsToPersistencePairs(
  std::vector<PersistencePair> &pairs,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<triplet> &triplets,
  const scalarType *const fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets,

  const bool splitTree) const {
  Timer tm;
  if(triplets.empty())
    return;
  size_t numberOfPairs = 0;

  // // accelerate getRep lookup?
  // std::vector<SimplexId> firstRep(vertexRepresentatives.size());

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif // TTK_ENABLE_OPENMP
  //   for(size_t i = 0; i < firstRep.size(); ++i) {
  //     if(!vertexRepresentatives[i].empty()) {
  //       firstRep[i] = vertexRepresentatives[i][0];
  //     }
  //   }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return ((fakeScalars[a] < fakeScalars[b])
            || (fakeScalars[a] == fakeScalars[b]
                && ((monotonyOffsets[a] < monotonyOffsets[b])
                    || (monotonyOffsets[a] == monotonyOffsets[b]
                        && offsets[a] < offsets[b]))));
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
        // fakeScalars[r1] - fakeScalars[s], 2);

      } else {
        // r1 = max(r1, r2), r2 = min(r1, r2)
        if(lt(r1, r2)) {
          std::swap(r1, r2);
        }
        // pair min-saddle: max(r1, r2) -> s;
        pairs.emplace_back(r1, s, 0);
        // fakeScalars[s] - fakeScalars[r1], 0);
      }

      vertexRepresentatives[std::get<1>(t)][0] = r2;
      vertexRepresentatives[r1][0] = r2;
      // vertexRepresentatives[std::get<1>(t)][0] = r2;
      // vertexRepresentatives[r1][0] = r2;
    }
  }

  if(debugLevel_ > 3) {
    std::string prefix = splitTree ? "[sad-max]" : "[min-sad]";
    std::cout << prefix << "  found all pairs in " << tm.getElapsedTime()
              << " s." << std::endl;
  }
}

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::sortVertices(
  const SimplexId vertexNumber,
  std::vector<SimplexId> &sortedVertices,
  SimplexId *vertsOrder,
  const scalarType *const fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) {

  Timer tm;

  sortedVertices.resize(vertexNumber);

  // fill with numbers from 0 to vertexNumber - 1
  std::iota(sortedVertices.begin(), sortedVertices.end(), 0);

  // sort vertices in ascending order following scalarfield / offsets
  TTK_PSORT(this->threadNumber_, sortedVertices.begin(), sortedVertices.end(),
            [&](const SimplexId a, const SimplexId b) {
              return ((fakeScalars[a] < fakeScalars[b])
                      || (fakeScalars[a] == fakeScalars[b]
                          && ((monotonyOffsets[a] < monotonyOffsets[b])
                              || (monotonyOffsets[a] == monotonyOffsets[b]
                                  && offsets[a] < offsets[b]))));
            });

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sortedVertices.size(); ++i) {
    vertsOrder[sortedVertices[i]] = i;
  }

  // if(debugLevel_ > 2) {
  //   std::cout << "SORT " << tm.getElapsedTime() << std::endl;
  // }
}

template <typename scalarType, typename offsetType>
int ttk::ApproximateTopology::getMonotonyChangeByOldPointCPApproximate(
  const SimplexId vertexId,
  double eps,
  const std::vector<polarity> &isNew,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  std::vector<std::pair<polarity, polarity>> &vlp,
  scalarType *fakeScalars,
  const offsetType *const offsets,
  int *monotonyOffsets) const {

  int hasMonotonyChanged = 0;
  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);

    // check for monotony changes
    // const bool lowerStatic
    //   = (scalarField[neighborId] < fakeScalars[vertexId])
    //     || ((scalarField[neighborId] == fakeScalars[vertexId])
    //         && (offsets[neighborId] < offsets[vertexId]));
    const bool lowerDynamic
      = ((fakeScalars[neighborId] < fakeScalars[vertexId])
         || (fakeScalars[neighborId] == fakeScalars[vertexId]
             && ((monotonyOffsets[neighborId] < monotonyOffsets[vertexId])
                 || (monotonyOffsets[neighborId] == monotonyOffsets[vertexId]
                     && offsets[neighborId] < offsets[vertexId]))));

    const polarity isUpperDynamic = lowerDynamic ? 0 : 255;
    // const polarity isUpperStatic = lowerStatic ? 0 : 255;
    const polarity isUpperOld = vlp[i].first;

    if(isUpperDynamic != isUpperOld) { // change of monotony
      SimplexId oldNeighbor = -1;
      int oldDecimation = pow(2, decimationLevel_ + 1);
      multiresTriangulation_.getVertexNeighborAtDecimation(
        vertexId, i, oldNeighbor, oldDecimation);

      double replacementValueDynamic
        = (0.5 * (double)fakeScalars[oldNeighbor]
           + .5 * (double)fakeScalars[vertexId]); // depends on epsilon
      double deltaDynamic
        = fabs((double)fakeScalars[neighborId] - replacementValueDynamic);

      //=====================
      SimplexId oldNeighNumber = 0;
      SimplexId nnumber
        = multiresTriangulation_.getVertexNeighborNumber(neighborId);
      for(SimplexId iii = 0; iii < nnumber; iii++) {
        SimplexId neighborId2 = -1;
        multiresTriangulation_.getVertexNeighbor(neighborId, iii, neighborId2);
        if(!isNew[neighborId2]) {
          oldNeighNumber++;
        }
      }

      if(deltaDynamic > eps or !isNew[neighborId] or oldNeighNumber > 2) {
        hasMonotonyChanged = 1;

        toReprocess[vertexId] = 255;
        if(isNew[neighborId]) {
          toProcess[neighborId] = 255;
        } else {
          toReprocess[neighborId] = 255;
        }
        const SimplexId neighborNumberNew
          = multiresTriangulation_.getVertexNeighborNumber(neighborId);
        for(SimplexId j = 0; j < neighborNumberNew; j++) {
          SimplexId neighborIdNew = -1;
          multiresTriangulation_.getVertexNeighbor(
            neighborId, j, neighborIdNew);
          if(isNew[neighborIdNew])
            toProcess[neighborIdNew] = 255;
        }
        vlp[i].second = 255;
      } else {
        fakeScalars[neighborId] = replacementValueDynamic;

        // corrects rounding error when we process an integer scalar field
        if(fakeScalars[neighborId] == fakeScalars[oldNeighbor]) {
          fakeScalars[neighborId] = fakeScalars[vertexId];
        }

        // change monotony offset
        // The monotony must be preserved, meaning that
        // if we had f(vertexId)>f(oldNeighbor)
        // we should have f(vertexId)>f(neighborId)
        // which is enforced when
        // monotonyOffsets[vertexId]>monotonyOffsets[neighborId]
        if(isUpperOld) { // we should enforce f(vertexId)<f(neighborId)
          if(offsets[vertexId] > offsets[neighborId]) {
            monotonyOffsets[neighborId]
              = monotonyOffsets[vertexId] + pow(2, decimationLevel_);
            if(monotonyOffsets[vertexId] == monotonyOffsets[oldNeighbor]
               and fakeScalars[vertexId] == fakeScalars[oldNeighbor]) {
              std::cout << "THIS IS AN ISSUE" << std::endl;
            }
          } else {
            monotonyOffsets[neighborId] = monotonyOffsets[vertexId];
          }
        } else { // we should enforce f(vertexId)>f(neighborId)
          if(offsets[vertexId] < offsets[neighborId]) {
            monotonyOffsets[neighborId]
              = monotonyOffsets[vertexId] - pow(2, decimationLevel_);
            if(monotonyOffsets[vertexId] == monotonyOffsets[oldNeighbor]
               and fakeScalars[vertexId] == fakeScalars[oldNeighbor]) {
              std::cout << "THIS IS AN ISSUE" << std::endl;
            }
          } else {
            monotonyOffsets[neighborId] = monotonyOffsets[vertexId];
          }
        }
      }
    } // end if change of monotony
  } // end for neighbors
  return hasMonotonyChanged;
}

template <typename scalarType, typename offsetType>
ttk::SimplexId ttk::ApproximateTopology::propageFromSaddles(
  const SimplexId vertexId,
  std::vector<Lock> &vertLock,
  std::vector<polarity> &toPropage,
  std::vector<std::vector<SimplexId>> &vertexRepresentatives,
  std::vector<std::vector<SimplexId>> &saddleCC,
  std::vector<polarity> &isUpdated,
  std::vector<SimplexId> &globalExtremum,
  const bool splitTree,
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) const {

  auto &toProp = toPropage[vertexId];
  auto &reps = vertexRepresentatives[vertexId];
  auto &updated = isUpdated[vertexId];

  if(updated) {
    return reps[0];
  }

  // const auto gt = [=](const SimplexId v1, const SimplexId v2) {
  //   return ((fakeScalars[v1] > fakeScalars[v2])
  //           || (fakeScalars[v1] == fakeScalars[v2]
  //               && offsets[v1] > offsets[v2]))
  //          == splitTree;
  // };
  const auto gt = [=](const SimplexId v1, const SimplexId v2) {
    return ((fakeScalars[v1] > fakeScalars[v2])
            || (fakeScalars[v1] == fakeScalars[v2]
                && ((monotonyOffsets[v1] > monotonyOffsets[v2])
                    || (monotonyOffsets[v1] == monotonyOffsets[v2]
                        && offsets[v1] > offsets[v2]))))
           == splitTree;
  };

  if(this->threadNumber_ > 1) {
    vertLock[vertexId].lock();
  }
  if(saddleCC[vertexId].size()
     and !toProp) { // tis a saddle point, should have to propage on it
    printErr("ERRRROR");
  }
  if(toProp) { // SADDLE POINT
    if(debugLevel_ > 5)
      printMsg("to saddle " + std::to_string(vertexId) + " "
               + std::to_string(saddleCC[vertexId].size()));
    // for(auto ccid : saddleCC[vertexId]) {
    //   SimplexId ccV = -1;
    //   multiresTriangulation_.getVertexNeighbor(vertexId, ccid, ccV);
    //   std::cout << " " << ccid << " (" << ccV << ") , ";
    // }
    // std::cout << std::endl;
    const auto &CC = saddleCC[vertexId];
    reps.clear();
    reps.reserve(CC.size());
    for(size_t r = 0; r < CC.size(); r++) {
      SimplexId neighborId = -1;
      SimplexId localId = CC[r];
      multiresTriangulation_.getVertexNeighbor(vertexId, localId, neighborId);
      // printMsg("       CC " + std::to_string(CC[r]) + " "
      //          + std::to_string(neighborId) + " ("
      //          + std::to_string(fakeScalars[neighborId]) + " , "
      //          + std::to_string(offsets[neighborId]) + ")");
      SimplexId ret = propageFromSaddles(neighborId, vertLock, toPropage,
                                         vertexRepresentatives, saddleCC,
                                         isUpdated, globalExtremum, splitTree,
                                         fakeScalars, offsets, monotonyOffsets);
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
    if(debugLevel_ > 5)
      printMsg("to non saddle " + std::to_string(vertexId) + " "
               + std::to_string(saddleCC[vertexId].size()));

    SimplexId ret = vertexId;
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    SimplexId maxNeighbor = vertexId;
    // std::cout << "neigh number" << neighborNumber << std::endl;
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
                               globalExtremum, splitTree, fakeScalars, offsets,
                               monotonyOffsets);

    } else { // needed to find the globalExtremum per thread
#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP
      if(gt(vertexId, globalExtremum[tid])) {
        // if(splitTree)
        //   std::cout << "new global max " << std::endl;
        // else
        //   std::cout << "new global min " << std::endl;
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

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::updatePropagation(
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
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) {

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

  // reset updated flag
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    isUpdatedMin[v] = 0;
    isUpdatedMax[v] = 0;
  }

  // propage along split tree
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);
    if(toPropageMin[v]) {
      propageFromSaddles(v, vertLockMin, toPropageMin, vertexRepresentativesMin,
                         saddleCCMin, isUpdatedMin, globalMinThr, false,
                         fakeScalars, offsets, monotonyOffsets);
    }
    if(toPropageMax[v]) {
      propageFromSaddles(v, vertLockMax, toPropageMax, vertexRepresentativesMax,
                         saddleCCMax, isUpdatedMax, globalMaxThr, true,
                         fakeScalars, offsets, monotonyOffsets);
    }
  }

  const auto lt = [=](const SimplexId a, const SimplexId b) -> bool {
    return ((fakeScalars[a] < fakeScalars[b])
            || (fakeScalars[a] == fakeScalars[b]
                && ((monotonyOffsets[a] < monotonyOffsets[b])
                    || (monotonyOffsets[a] == monotonyOffsets[b]
                        && offsets[a] < offsets[b]))));
  };

  globalMin_ = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
  globalMax_ = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);

  if(globalMin_ == 0 or globalMax_ == 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDecVerts; i++) {
      SimplexId v = multiresTriangulation_.localToGlobalVertexId(i);

#ifdef TTK_ENABLE_OPENMP
      const auto tid = omp_get_thread_num();
#else
      const auto tid = 0;
#endif // TTK_ENABLE_OPENMP

      if(lt(globalMaxThr[tid], v)) {
        globalMaxThr[tid] = v;
      }
      if(lt(v, globalMinThr[tid])) {
        globalMinThr[tid] = v;
      }
    }
    globalMin_
      = *std::min_element(globalMinThr.begin(), globalMinThr.end(), lt);
    globalMax_
      = *std::max_element(globalMaxThr.begin(), globalMaxThr.end(), lt);
    // printMsg("Explicitely found global extremas");
  }
  if(debugLevel_ > 3) {
    printMsg("Propagation Update", 1, tm.getElapsedTime(), threadNumber_);
  }
}

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::buildVertexLinkPolarityApproximate(
  const SimplexId vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) const {

  const SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vlp.resize(neighborNumber);

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);

    // const bool lower0 = vertsOrder_[neighborId0] < vertsOrder_[vertexId];
    const bool lower0
      = ((fakeScalars[neighborId0] < fakeScalars[vertexId])
         || (fakeScalars[neighborId0] == fakeScalars[vertexId]
             && ((monotonyOffsets[neighborId0] < monotonyOffsets[vertexId])
                 || (monotonyOffsets[neighborId0] == monotonyOffsets[vertexId]
                     && offsets[neighborId0] < offsets[vertexId]))));

    const polarity isUpper0 = static_cast<polarity>(!lower0) * 255;
    vlp[i] = std::make_pair(isUpper0, 0);
  }
}

template <typename scalarType, typename offsetType>
bool ttk::ApproximateTopology::printPolarity(
  std::vector<polarity> &isNew,
  const SimplexId vertexId,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  const scalarType *scalars,
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets,
  bool verbose) {

  bool error = false;
  std::stringstream mymsg;
  std::vector<std::pair<polarity, polarity>> &vlp
    = vertexLinkPolarity[vertexId];
  mymsg << "POLARITY PRINT"
        << "\n";
  mymsg << "vertex " << vertexId << " has "
        << multiresTriangulation_.getVertexNeighborNumber(vertexId)
        << " neighbors"
        << "\n";
  mymsg << "\tself  f:" << fakeScalars[vertexId] << " s:" << scalars[vertexId]
        << " o:" << offsets[vertexId] << " m:" << monotonyOffsets[vertexId]
        << "  isnew: " << (int)isNew[vertexId] << "\n";
  if(vlp.empty()) {
    if(verbose) {
      std::cout << "\tpolarity not initialized for " << vertexId << std::endl;
      std::cout << mymsg.str() << std::endl;
    }
    return false;
  }
  for(SimplexId i = 0;
      i < multiresTriangulation_.getVertexNeighborNumber(vertexId); i++) {
    SimplexId nId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, nId);
    // find reverse pol
    std::vector<std::pair<polarity, polarity>> &vlp2 = vertexLinkPolarity[nId];
    bool init = false;
    polarity rpol;
    if(!vlp2.empty()) {
      for(SimplexId j = 0;
          j < multiresTriangulation_.getVertexNeighborNumber(nId); j++) {
        SimplexId nId2 = -1;
        multiresTriangulation_.getVertexNeighbor(nId, j, nId2);
        if(nId2 == vertexId) {
          rpol = vlp2[j].first;
          init = true;
        }
      }
    }

    const bool lower
      = ((fakeScalars[nId] < fakeScalars[vertexId])
         || (fakeScalars[nId] == fakeScalars[vertexId]
             && ((monotonyOffsets[nId] < monotonyOffsets[vertexId])
                 || (monotonyOffsets[nId] == monotonyOffsets[vertexId]
                     && offsets[nId] < offsets[vertexId]))));

    const polarity isUpper = lower ? 0 : 255;

    mymsg << " " << i << "th: " << nId << " f:" << fakeScalars[nId]
          << " s:" << scalars[nId] << " o:" << offsets[nId]
          << " m:" << monotonyOffsets[nId] << "  , pol:" << (bool)vlp[i].first
          << "(" << (bool)vlp[i].second << ")"
          << " rpol:" << (bool)rpol << "  true pol:" << (bool)isUpper
          << " init " << init << "  isnew: " << (int)isNew[nId] << "\n";
    if((rpol == isUpper and !vlp2.empty())
       or (isUpper != vlp[i].first and !vlp[i].second)) {
      mymsg << "POLARITY ERROR "
            << "\n";
      error = true;
    }
  }
  if(error or verbose) {
    std::cout << mymsg.str() << std::endl;
  }
  return error;
}

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::getCriticalTypeApproximate(
  const SimplexId &vertexId,
  std::vector<std::pair<polarity, polarity>> &vlp,
  uint8_t &vertexLink,
  DynamicTree &link,
  VLBoundaryType &vlbt,
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) const {

  if(vlp.empty()) {
    buildVertexLinkPolarityApproximate(
      vertexId, vlp, fakeScalars, offsets, monotonyOffsets);
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

template <typename scalarType, typename offsetType>
void ttk::ApproximateTopology::initGlobalPolarity(
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  const scalarType *fakeScalars,
  const offsetType *const offsets,
  const int *const monotonyOffsets) const {

  Timer timer{};
  const size_t nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

  // computes the critical types of all points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDecVerts; i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    buildVertexLinkPolarityApproximate(globalId, vertexLinkPolarity[globalId],
                                       fakeScalars, offsets, monotonyOffsets);
    toProcess[globalId] = 255;
    isNew[globalId] = 0;
  }
  printMsg("Polarity Init", 1, timer.getElapsedTime(), threadNumber_,
           debug::LineMode::NEW, debug::Priority::DETAIL);
}

template <typename ScalarType, typename offsetType>
int ttk::ApproximateTopology::updateGlobalPolarity(
  double eps,
  std::vector<polarity> &isNew,
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toProcess,
  std::vector<polarity> &toReprocess,
  ScalarType *fakeScalars,
  const offsetType *const offsets,
  int *monotonyOffsets) const {

  Timer tm;
  const auto nDecVerts = multiresTriangulation_.getDecimatedVertexNumber();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId localId = 0; localId < nDecVerts; localId++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(localId);
    if(!isNew[globalId]) {
      getMonotonyChangeByOldPointCPApproximate(
        globalId, eps, isNew, toProcess, toReprocess,
        vertexLinkPolarity[globalId], fakeScalars, offsets, monotonyOffsets);
    }
  }

  // second Loop  process or reprocess
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < multiresTriangulation_.getDecimatedVertexNumber(); i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
    if(isNew[globalId]) { // new point
      if(decimationLevel_ > stoppingDecimationLevel_) {
        buildVertexLinkPolarityApproximate(
          globalId, vertexLinkPolarity[globalId], fakeScalars, offsets,
          monotonyOffsets);
      }
      isNew[globalId] = 0;

    } else { // old point
      if(toReprocess[globalId]) {
        updateLinkPolarityPonctual(vertexLinkPolarity[globalId]);

        toProcess[globalId] = 255; // mark for processing
        toReprocess[globalId] = 0;
      }
    }
  } // end for openmp
  return 0;
}

template <typename ScalarType, typename offsetType>
void ttk::ApproximateTopology::computeCriticalPoints(
  std::vector<std::vector<std::pair<polarity, polarity>>> &vertexLinkPolarity,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<polarity> &toProcess,
  std::vector<DynamicTree> &link,
  std::vector<uint8_t> &vertexLink,
  VLBoundaryType &vertexLinkByBoundaryType,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax,
  ScalarType *fakeScalars,
  const offsetType *const offsets,
  int *monotonyOffsets) {

  Timer tm;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < multiresTriangulation_.getDecimatedVertexNumber(); i++) {
    SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);

    if(toProcess[globalId]) {
      // nbComputations++;
      getCriticalTypeApproximate(globalId, vertexLinkPolarity[globalId],
                                 vertexLink[globalId], link[globalId],
                                 vertexLinkByBoundaryType, fakeScalars, offsets,
                                 monotonyOffsets);
      getValencesFromLink(globalId, vertexLinkPolarity[globalId],
                          link[globalId], toPropageMin, toPropageMax,
                          saddleCCMin, saddleCCMax);
    }
  } // end for openmp

  if(debugLevel_ > 3) {
    printMsg(
      "Critical Points Computation", 1, tm.getElapsedTime(), threadNumber_);
  }
}

template <typename scalarType>
int ttk::ApproximateTopology::sortPersistenceDiagramApproximate(
  std::vector<PersistencePair> &diagram,
  scalarType *fakeScalars,
  const SimplexId *const offsets,
  int *monotonyOffsets) const {
  auto cmp = [fakeScalars, offsets, monotonyOffsets](
               const PersistencePair &pA, const PersistencePair &pB) {
    const ttk::SimplexId a = pA.birth;
    const ttk::SimplexId b = pB.birth;

    return ((fakeScalars[a] < fakeScalars[b])
            || (fakeScalars[a] == fakeScalars[b]
                && ((monotonyOffsets[a] < monotonyOffsets[b])
                    || (monotonyOffsets[a] == monotonyOffsets[b]
                        && offsets[a] < offsets[b]))));
  };

  std::sort(diagram.begin(), diagram.end(), cmp);

  return 0;
}

template <typename scalarType>
int ttk::ApproximateTopology::computeApproximatePD(
  std::vector<PersistencePair> &CTDiagram,
  const scalarType *scalars,
  scalarType *const fakeScalars,
  SimplexId *const outputOffsets,
  int *const outputMonotonyOffsets) {

  int ret = -1;
  std::stringstream ss;
  ss << "Approximate Persistence Diagram computation with "
     << debug::output::UNDERLINED << debug::output::YELLOW << epsilon_ * 100
     << "%" << debug::output::ENDCOLOR << debug::output::ENDCOLOR << " error";
  printMsg(ss.str());

  ret = executeApproximateTopology(
    scalars, fakeScalars, outputOffsets, outputMonotonyOffsets);
  CTDiagram = std::move(CTDiagram_);
  CTDiagram_.clear();

  return ret;
}
