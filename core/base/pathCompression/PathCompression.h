/// \ingroup base
/// \class ttk::PathCompression
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK processing package for the computation of Morse-Smale
/// segmentations using Path Compression.
///
/// Given an input order field, this class computes its ascending and descending
/// segmentation by assigning every vertex to its minimum or maximum in gradient
/// or inverse gradient direction. For convienience a hash (no hash collision
/// detection) of both segmentations can be created to represent the Morse-Smale
/// segmentation.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttkPathCompression.cpp %for a usage example
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleSegmentation_at/">Morse-Smale
///   segmentation example</a> \n

#pragma once

// base code includes
#include <Triangulation.h>

using ttk::SimplexId;

namespace ttk {
  namespace pcp {
#ifdef TTK_ENABLE_64BIT_IDS
    /**
     * @brief Get a hash value from two keys
     *
     * @param a First Hash key
     * @param b Second hash key
     *
     * @return Hash value
     */
    constexpr unsigned long long int getHash(const unsigned long long int a,
                                             const unsigned long long int b) {
      return (a * b + (a * a) + (b * b) + (a * a * a) * (b * b * b))
             % ULONG_LONG_MAX;
    }
#else
    /**
     * @brief Get a hash value from two keys
     *
     * @param a First Hash key
     * @param b Second hash key
     *
     * @return Hash value
     */
    constexpr unsigned int getHash(const unsigned int a, const unsigned int b) {
      return (a * b + (a * a) + (b * b) + (a * a * a) * (b * b * b)) % UINT_MAX;
    }
#endif
  } // namespace pcp
  class PathCompression : public virtual Debug {
  public:
    PathCompression();

    /** @brief Pointers to pre-allocated segmentation point data arrays */
    struct OutputSegmentation {
      SimplexId *ascending_;
      SimplexId *descending_;
      SimplexId *morseSmale_;
    };

    /**
     * Compute necessary triangulation information
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      data->preconditionVertexNeighbors();
    }

    /**
     * @brief Main function for computing the Morse-Smale complex.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] outSegmentation segmentations as a struct
     * @param[in] orderArray order field
     * @param[in] triangulation triangulation
     *
     * @pre PathCompression::preconditionTriangulation must be
     * called prior to this.
     *
     * @return 0 on success
     */
    template <typename triangulationType>
    inline int execute(OutputSegmentation &outSegmentation,
                       const SimplexId *const orderArray,
                       const triangulationType &triangulation);

    /**
     * @brief Compute the ascending and descending segmentation in one run
     *
     * This function computes the ascending and descending segmentation on the
     * order field. First, the ascending and descending segmentation is set to
     * the largest and smallest neighbor of each vertex. Then, using path
     * compression, the vertices are assigned to their minimum/maximum in
     * positive/negative gradient direction.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] ascSegmentation ascending segmentation
     * @param[out] dscSegmentation descending segmentation
     * @param[in] orderArray order array
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computePathCompression(SimplexId *const ascSegmentation,
                               SimplexId *const dscSegmentation,
                               const SimplexId *const orderArray,
                               const triangulationType &triangulation) const;

    /**
     * @brief Compute the ascending or descending segmentation
     *
     * This function computes the ascending or descending segmentation on the
     * order field. First, the ascending or descending segmentation is set to
     * the largest or smallest neighbor of each vertex. Then, using path
     * compression, the vertices are assigned to their minimum/maximum in
     * positive/negative gradient direction.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] segmentation segmentation
     * @param[in] computeAscending compute the ascending or descending
     * segmentation
     * @param[in] orderArray order array
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computePathCompressionSingle(
      SimplexId *const segmentation,
      const bool computeAscending,
      const SimplexId *const orderArray,
      const triangulationType &triangulation) const;

    /**
     * @brief Computes a MS segmentation hash
     *
     * Computes a hash from the ascending and descending segmentation as keys.
     * The function does not check for hash conflicts.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] morseSmaleSegmentation
     * @param[out] ascSegmentation ascending segmentation
     * @param[out] dscSegmentation descending segmentation
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computeMSHash(SimplexId *const morseSmaleSegmentation,
                      const SimplexId *const ascSegmentation,
                      const SimplexId *const dscSegmentation,
                      const triangulationType &triangulation) const;

  protected:
    // Compute ascending segmentation?
    bool ComputeAscendingSegmentation{true};

    // Compute descending segmentation?
    bool ComputeDescendingSegmentation{true};

    // Compute Morse-Smale segmentation hash?
    bool ComputeMSSegmentationHash{true};
  };
} // namespace ttk

template <typename triangulationType>
int ttk::PathCompression::execute(OutputSegmentation &outSegmentation,
                                  const SimplexId *const orderArray,
                                  const triangulationType &triangulation) {
  if(orderArray == nullptr)
    return this->printErr("Input offset field pointer is null.");

  Timer t;

  this->printMsg("Start computing segmentations", 0.0, t.getElapsedTime(),
                 this->threadNumber_);

  if((ComputeAscendingSegmentation && ComputeDescendingSegmentation)
     || ComputeMSSegmentationHash) {
    computePathCompression(outSegmentation.ascending_,
                           outSegmentation.descending_, orderArray,
                           triangulation);
  } else if(ComputeAscendingSegmentation) {
    computePathCompressionSingle(
      outSegmentation.ascending_, true, orderArray, triangulation);
  } else if(ComputeDescendingSegmentation) {
    computePathCompressionSingle(
      outSegmentation.descending_, false, orderArray, triangulation);
  }

  if(ComputeMSSegmentationHash) {
    computeMSHash(outSegmentation.morseSmale_, outSegmentation.ascending_,
                  outSegmentation.descending_, triangulation);
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompression(
  SimplexId *const ascSegmentation,
  SimplexId *const dscSegmentation,
  const SimplexId *const orderArray,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const SimplexId nVertices = triangulation.getNumberOfVertices();
  std::vector<SimplexId> lActiveVertices;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_) private(lActiveVertices)
  {
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  lActiveVertices.reserve(nVertices);
#endif // TTK_ENABLE_OPENMP
    // find the largest and smallest neighbor for each vertex
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId{0};
      SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &dmi = dscSegmentation[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascSegmentation[i];
      ami = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(orderArray[neighborId] < orderArray[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        } else if(orderArray[neighborId] > orderArray[dmi]) {
          dmi = neighborId;
          hasLargerNeighbor = true;
        }
      }

      if(hasLargerNeighbor || hasSmallerNeighbor) {
        lActiveVertices.push_back(i);
      }
    }

    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId const &v = lActiveVertices[i];
        SimplexId &vDsc = dscSegmentation[v];
        SimplexId &vAsc = ascSegmentation[v];

// compress paths
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vDsc = dscSegmentation[vDsc];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vAsc = ascSegmentation[vAsc];

        // check if fully compressed
        if(vDsc != dscSegmentation[vDsc] || vAsc != ascSegmentation[vAsc]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Asc. and Desc. segmentation computed", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompressionSingle(
  SimplexId *const segmentation,
  const bool computeAscending,
  const SimplexId *const orderArray,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const SimplexId nVertices = triangulation.getNumberOfVertices();
  std::vector<SimplexId> lActiveVertices;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
  {
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  lActiveVertices.reserve(nVertices);
#endif // TTK_ENABLE_OPENMP
    // find the largest neighbor for each vertex
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId{0};
      SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &mi = segmentation[i];
      mi = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(computeAscending) {
          if(orderArray[neighborId] < orderArray[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        } else {
          if(orderArray[neighborId] > orderArray[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        }
      }

      if(hasLargerNeighbor) {
        lActiveVertices.push_back(i);
      }
    }

    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId const &v = lActiveVertices[i];
        SimplexId &vMan = segmentation[v];

// compress paths
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vMan = segmentation[vMan];

        // check if fully compressed
        if(vMan != segmentation[vMan]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

  if(computeAscending) {
    this->printMsg("Ascending segmentation computed", 1.0,
                   localTimer.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  } else {
    this->printMsg("Descending segmentation computed", 1.0,
                   localTimer.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  return 0; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computeMSHash(
  SimplexId *const morseSmaleSegmentation,
  const SimplexId *const ascSegmentation,
  const SimplexId *const dscSegmentation,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const size_t nVerts = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleSegmentation[i]
      = pcp::getHash(ascSegmentation[i], dscSegmentation[i]);
  }

  this->printMsg("Morse-Smale segmentation hash computed", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0;
}
