/// \ingroup base
/// \class ttk::PathCompression
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date January 2023.
///
/// \brief TTK processing package for the computation of Morse-Smale
/// segmentations using Path Compression.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Triangulation.h>


using ttk::SimplexId;

#ifdef TTK_ENABLE_64BIT_IDS
constexpr unsigned long long int hash_max = ULONG_LONG_MAX;

constexpr unsigned long long int getHash(
  const unsigned long long int a, const unsigned long long int b)  {
  return (a*b + (a*a) + (b*b) + (a*a*a)*(b*b*b)) % hash_max;
	//return std::rotl(a,1) ^ b;
}
#else
constexpr unsigned int hash_max = UINT_MAX;

/**
 * @brief Get a hash
 * 
 * @param a 
 * @param b 
 * @return constexpr unsigned int 
 */
constexpr unsigned int getHash(
  const unsigned int a, const unsigned int b) {
  return (a*b + (a*a) + (b*b) + (a*a*a)*(b*b*b)) % hash_max;
  //return std::rotl(a,1) ^ b;
}
#endif

namespace ttk {
  class PathCompression : public virtual Debug {
  public:
    PathCompression();

    /** @brief Pointers to pre-allocated segmentation point data arrays */
    struct OutputManifold {
      SimplexId *ascending_;
      SimplexId *descending_;
      SimplexId *morseSmale_;
    };

    /**
     * Main function for computing the Morse-Smale complex.
     *
     * @pre PathCompression::preconditionTriangulation must be
     * called prior to this.
     */
    template <typename dataType, typename triangulationType>
    inline int execute(OutputManifold &outManifold,
                       const dataType *const scalars,
                       const SimplexId *const offsets,
                       const triangulationType &triangulation);

    /**
     * Enable/Disable computation of the geometrical embedding of the
     * manifolds of the critical points.
     */
    inline void setComputeSegmentation(const bool doAscending,
                                       const bool doDescending,
                                       const bool doMorseSmale) {
      this->ComputeAscendingSegmentation = doAscending;
      this->ComputeDescendingSegmentation = doDescending;
      this->ComputeFinalSegmentation = doMorseSmale;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      data->preconditionVertexNeighbors();
    }

    template <typename triangulationType>
    int computePathCompression(
      SimplexId *const ascManifold,
      SimplexId *const dscManifold,
      const SimplexId *const orderArr,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computePathCompressionSingle(
      SimplexId *const manifold,
      const bool computeAscending,
      const SimplexId *const orderArr,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      SimplexId *const morseSmaleManifold,
      const SimplexId *const ascManifold,
      const SimplexId *const desManifold,
      const triangulationType &triangulation) const;

    bool ComputeDescendingSeparatrices2{false};
    bool ComputeAscendingSegmentation{true};
    bool ComputeDescendingSegmentation{true};
    bool ComputeFinalSegmentation{true};

  };
} // namespace ttk

// ---------------- //
//  Execute method  //
// ---------------- //

template <typename dataType, typename triangulationType>
int ttk::PathCompression::execute(OutputManifold &outManifold,
                                  const dataType *const scalars,
                                  const SimplexId *const orderArray,
                                  const triangulationType &triangulation) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(scalars == nullptr) {
    this->printErr("Input scalar field pointer is null.");
    return -1;
  }

  if(orderArray == nullptr) {
    this->printErr("Input offset field pointer is null.");
    return -1;
  }
#endif

  Timer t;

  if((ComputeAscendingSegmentation && ComputeDescendingSegmentation) ||
    ComputeFinalSegmentation) {
    computePathCompression(outManifold.ascending_,
      outManifold.descending_, orderArray, triangulation);
  } else if(ComputeAscendingSegmentation) {
    computePathCompressionSingle(
      outManifold.ascending_, true, orderArray, triangulation);
  } else if(ComputeDescendingSegmentation) {
    computePathCompressionSingle(
      outManifold.descending_, false, orderArray, triangulation);
  }

  if(ComputeFinalSegmentation) {
    computeFinalSegmentation(outManifold.morseSmale_, outManifold.ascending_,
      outManifold.descending_, triangulation);
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompression(
  SimplexId *const ascManifold,
  SimplexId *const dscManifold,
  const SimplexId *const orderArr,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  //this->printWrn(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printWrn("Computing Manifolds");

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  SimplexId nActiveVertices;
  std::vector<SimplexId> activeVertices;
  activeVertices.reserve(nVertices);
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation.getVertexNeighborNumber(i);      

    bool hasLargerNeighbor = false;
    SimplexId &dmi = dscManifold[i];
    dmi = i;

    bool hasSmallerNeighbor = false;
    SimplexId &ami = ascManifold[i];
    ami = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation.getVertexNeighbor(i, n, neighborId);

      if(orderArr[neighborId] < orderArr[ami]) {
        ami = neighborId;
        hasSmallerNeighbor = true;
      } else if(orderArr[neighborId] > orderArr[dmi]) {
        dmi = neighborId;
        hasLargerNeighbor = true;
      }
    }

    if(hasLargerNeighbor || hasSmallerNeighbor) {
      activeVertices.push_back(i);
    }
  }

  nActiveVertices = activeVertices.size();
  size_t currentIndex = 0;

  // compress paths until no changes occur
  while(nActiveVertices > 0) {      
    for(SimplexId i = 0; i < nActiveVertices; i++) {
      SimplexId &v = activeVertices[i];
      SimplexId &vDsc = dscManifold[v];
      SimplexId &vAsc = ascManifold[v];

      // compress paths
      vDsc = dscManifold[vDsc];
      vAsc = ascManifold[vAsc];

      // check if not fully compressed
      if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
        activeVertices[currentIndex++] = v;
      }
    }

    nActiveVertices = currentIndex;
    currentIndex = 0;
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  #pragma omp parallel num_threads(threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices;
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &dmi = dscManifold[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascManifold[i];
      ami = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(orderArr[neighborId] < orderArr[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        } else if(orderArr[neighborId] > orderArr[dmi]) {
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
        SimplexId &v = lActiveVertices[i];
        SimplexId &vDsc = dscManifold[v];
        SimplexId &vAsc = ascManifold[v];

        // compress paths
        #pragma omp atomic read
        vDsc = dscManifold[vDsc];

        #pragma omp atomic read
        vAsc = ascManifold[vAsc];

        // check if fully compressed
        if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printWrn("Computed Manifolds");

  return 1; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompressionSingle(
  SimplexId *const manifold,
  const bool computeAscending,
  const SimplexId *const orderArr,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  //this->printWrn(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printWrn("Computing Manifolds");

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  SimplexId nActiveVertices;
  std::vector<SimplexId> activeVertices;
  activeVertices.reserve(nVertices);
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation.getVertexNeighborNumber(i);      

    bool hasLargerNeighbor = false;
    SimplexId &mi = manifold[i];
    mi = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation.getVertexNeighbor(i, n, neighborId);

      if(computeAscending) {
        if(orderArr[neighborId] < orderArr[mi]) {
          mi = neighborId;
          hasLargerNeighbor = true;
        }
      } else {
        if(orderArr[neighborId] > orderArr[mi]) {
          mi = neighborId;
          hasLargerNeighbor = true;
        }
      }
    }

    if(hasLargerNeighbor) {
      activeVertices.push_back(i);
    }
  }

  nActiveVertices = activeVertices.size();
  size_t currentIndex = 0;

  // compress paths until no changes occur
  while(nActiveVertices > 0) {      
    for(SimplexId i = 0; i < nActiveVertices; i++) {
      SimplexId &v = activeVertices[i];
      SimplexId &vMan = manifold[v];

      // compress paths
      vMan = manifold[vMan];

      // check if not fully compressed
      if(vMan != manifold[vMan]) {
        activeVertices[currentIndex++] = v;
      }
    }

    nActiveVertices = currentIndex;
    currentIndex = 0;
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  #pragma omp parallel num_threads(threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices;
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &mi = manifold[i];
      mi = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(computeAscending) {
          if(orderArr[neighborId] < orderArr[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        } else {
          if(orderArr[neighborId] > orderArr[mi]) {
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
        SimplexId &v = lActiveVertices[i];
        SimplexId &vMan = manifold[v];

        // compress paths
        #pragma omp atomic read
        vMan = manifold[vMan];

        // check if fully compressed
        if(vMan != manifold[vMan]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printWrn("Computed Manifolds");

  return 1; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computeFinalSegmentation(
  SimplexId *const morseSmaleManifold,
  const SimplexId *const ascManifold,
  const SimplexId *const dscManifold,
  const triangulationType &triangulation) const {
    ttk::Timer localTimer;

  this->printMsg("Computing MSC Manifold",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
  const size_t nVerts = triangulation.getNumberOfVertices();
  

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = getHash(ascManifold[i], dscManifold[i]);
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  #pragma omp parallel for schedule(static) num_threads(threadNumber_)
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = getHash(ascManifold[i], dscManifold[i]);
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed MSC Manifold",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}
