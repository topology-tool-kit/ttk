#pragma once

#include <BaseClass.h>
#include <iostream>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using IT = long long int;
namespace ttk {

  struct value {
    float scalar;
    IT globalId;
    IT localId;
    IT ordering = 0;

    value(float _scalar, IT _globalId, IT _localId)
      : scalar(_scalar), globalId(_globalId), localId(_localId) {
    }
  };

  // creates a vector of value structs based on the given pointers
  // only takes values of vertices which mainly belong to the current rank
  // ( so only vertices which are no ghost cells)
  inline std::
    tuple<std::vector<value>, std::vector<IT>, std::unordered_map<IT, IT>>
    populateVector(const size_t nVerts,
                   const float *const scalars,
                   const IT *const globalIds,
                   const char *const ghostCells) {
    std::vector<value> valuesToSortVector;
    std::vector<IT> gidsToGetVector;
    std::unordered_map<IT, IT> gidToLidMap;
    for(size_t i = 0; i < nVerts; i++) {
      IT globalId = globalIds[i];
      IT localId = i;
      gidToLidMap[globalId] = localId;
      if((int)ghostCells[i] == 0) {
        float scalarValue = scalars[i];
        valuesToSortVector.emplace_back(scalarValue, globalId, localId);
      } else {
        gidsToGetVector.push_back(globalId);
      }
    }
    return std::make_tuple(valuesToSortVector, gidsToGetVector, gidToLidMap);
  }

  // orders an value vector first by their scalar value and then by global id
  inline void sortVerticesDistributed(std::vector<value> &values) {
    std::sort(values.begin(), values.end(), [](value v1, value v2) {
      return (v1.scalar < v2.scalar)
             || (v1.scalar == v2.scalar && v1.globalId < v2.globalId);
    });
  }

  // send the highest burstSize values and decrease the vector by that amount
  // check if there are actually that many elements in the vector
  inline std::vector<value> returnVectorForBurstsize(std::vector<value> &values,
                                                     size_t burstSize) {
    std::vector<value> outVector;
    if(burstSize > values.size()) {
      outVector.assign(values.begin(), values.end());
      values.clear();
    } else {
      outVector.assign(values.end() - burstSize, values.end());
      values.erase(values.end() - burstSize, values.end());
    }

    return outVector;
  }

  // takes in an ordered (as defined above) vector of values and creates an
  // ordermap for each scalar value
  inline std::unordered_map<float, IT> buildOrderMap(std::vector<value> &values,
                                                     size_t totalSize) {
    std::unordered_map<float, IT> orderMap;

    // omp only creates larger overhead because orderMap needs to accessed
    // critically
    for(size_t i = 0; i < totalSize; i++) {
      float scalarVal = values[i].scalar;
      IT orderVal = totalSize - i - 1;
      orderMap[scalarVal] = orderVal;
    }
    return orderMap;
  }

  /**
   * @brief Sort vertices according to scalars disambiguated by offsets
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars array of size nVerts, the scalar values which we want to
   * order
   * @param[in] orderMap map which maps scalar values to a defined order
   * @param[out] order array of size nVerts, computed order of vertices
   * @param[in] nThreads number of parallel threads
   */
  inline void buildOrderArray(const size_t nVerts,
                              const float *const scalars,
                              std::unordered_map<float, IT> &orderMap,
                              SimplexId *const order,
                              const int nThreads) {

    TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nVerts; ++i) {
      order[i] = orderMap[scalars[i]];
    }
  }

  /**
   * @brief Sort vertices according to scalars disambiguated by offsets
   *
   * @param[in] nInts number of long ints, as globalid / order pairs
   * @param[in] orderedValuesForRank array of size nInts, the actual pairs, [i]
   * = gId, [i+1] = order
   * @param[in] gidToLidMap map which maps scalar values to a defined order
   * @param[out] order array of size nVerts, computed order of vertices, this
   * procedure doesn't fill it completely
   * @param[in] nThreads number of parallel threads
   */
  inline void buildArrayForReceivedData(const size_t nInts,
                                        const IT *const orderedValuesForRank,
                                        std::unordered_map<IT, IT> &gidToLidMap,
                                        SimplexId *const order,
                                        const int nThreads) {

    TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nInts; i += 2) {
      order[gidToLidMap[orderedValuesForRank[i]]] = orderedValuesForRank[i + 1];
    }
  }

  // gets all the neighbors of a rank based on the rankarray of the vertices
  // belonging to this rank
  inline std::unordered_set<int> getNeighbors(const size_t nVertices,
                                              const int ownRank,
                                              const int *const rankArray) {
    std::unordered_set<int> neighbors;
    for(size_t i = 0; i < nVertices; i++) {
      int value = rankArray[i];
      if(value != ownRank)
        neighbors.emplace(value);
    }
    return neighbors;
  }

  // returns a vector of gid vectors, one vector for each neighbor
  // this shows us where we need to get the order for these gids from
  inline std::vector<std::vector<IT>>
    getGIdForNeighbors(const size_t nGIds,
                       const size_t nNeighbors,
                       std::vector<IT> &gidsToGetVector,
                       std::unordered_set<int> &neighbors,
                       std::unordered_map<IT, IT> &gidToLidMap,
                       const int *const rankArray) {
    std::vector<std::vector<IT>> outVector;
    outVector.resize(nNeighbors);
    for(size_t i = 0; i < nGIds; i++) {
      IT gId = gidsToGetVector[i];
      IT lId = gidToLidMap[gId];
      int rankForGId = rankArray[lId];
      auto distance
        = std::distance(neighbors.begin(), neighbors.find(rankForGId));
      outVector[distance].push_back(gId);
    }
    return outVector;
  }

  inline std::vector<IT>
    getOrderForGIds(const size_t nGIds,
                    const IT *const gIds,
                    std::unordered_map<IT, IT> &gidToLidMap,
                    const SimplexId *const order,
                    const int nThreads) {

    TTK_FORCE_USE(nThreads);
    std::vector<IT> outVector;
    outVector.resize(nGIds * 2);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nGIds; i++) {
      IT gId = gIds[i];
      IT lId = gidToLidMap[gId];
      SimplexId orderForThisGId = order[lId];
      outVector[2 * i] = gId;
      outVector[2 * i + 1] = orderForThisGId;
    }
    return outVector;
  }
  /**
   * @brief Sort vertices according to scalars disambiguated by offsets
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars array of size nVerts, main vertex comparator
   * @param[in] offsets array of size nVerts, disambiguate scalars on plateaux
   * @param[out] order array of size nVerts, computed order of vertices
   * @param[in] nThreads number of parallel threads
   */
  template <typename scalarType, typename idType>
  void sortVertices(const size_t nVerts,
                    const scalarType *const scalars,
                    const idType *const offsets,
                    SimplexId *const order,
                    const int nThreads) {

    // array of pre-sorted vertices
    std::vector<SimplexId> sortedVertices(nVerts);

    TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      sortedVertices[i] = i;
    }

    if(offsets != nullptr) {
      TTK_PSORT(
        nThreads, sortedVertices.begin(), sortedVertices.end(),
        [&](const SimplexId a, const SimplexId b) {
          return (scalars[a] < scalars[b])
                 || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
        });
    } else {
      TTK_PSORT(nThreads, sortedVertices.begin(), sortedVertices.end(),
                [&](const SimplexId a, const SimplexId b) {
                  return (scalars[a] < scalars[b])
                         || (scalars[a] == scalars[b] && a < b);
                });
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      order[sortedVertices[i]] = i;
    }
  }

  /**
   * @brief Precondition an order array to be consumed by the base layer API
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars pointer to scalar field buffer of size @p nVerts
   * @param[out] order pointer to pre-allocated order buffer of size @p nVerts
   * @param[in] nThreads number of threads to be used
   */
  template <typename scalarType>
  inline void preconditionOrderArray(const size_t nVerts,
                                     const scalarType *const scalars,
                                     SimplexId *const order,
                                     const int nThreads
                                     = ttk::globalThreadNumber_) {
    ttk::sortVertices(
      nVerts, scalars, static_cast<int *>(nullptr), order, nThreads);
  }
} // namespace ttk
