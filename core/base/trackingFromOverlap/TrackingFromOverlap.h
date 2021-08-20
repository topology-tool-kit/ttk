/// \ingroup base
/// \class ttk::TrackingFromOverlap
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %trackingFromOverlap processing package that tracks labled point
/// sets.
///
/// %TrackingFromOverlap is a TTK processing package that provides algorithms to
/// track labled point sets across time (and optionally levels) based on spatial
/// overlap, where two points overlap iff their corresponding coordinates are
/// equal.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///

#pragma once

#include <Debug.h>
#include <algorithm>
#include <boost/variant.hpp>
#include <map>
#include <unordered_map>

using topologyType = unsigned char;
using idType = long long int;

using labelTypeVariant = boost::variant<double,
                                        float,
                                        long long,
                                        unsigned long long,
                                        long,
                                        unsigned long,
                                        int,
                                        unsigned int,
                                        short,
                                        unsigned short,
                                        char,
                                        signed char,
                                        unsigned char>;
using sizeType = float;

namespace ttk {
  class TrackingFromOverlap : virtual public Debug {
  public:
    TrackingFromOverlap() {
      this->setDebugMsgPrefix("TrackingFromOverlap");
    }
    ~TrackingFromOverlap() = default;

    struct Node {
      labelTypeVariant label{};
      sizeType size{};
      float x{};
      float y{};
      float z{};

      idType branchID{-1};
      idType maxPredID{-1};
      idType maxSuccID{-1};

      Node() = default;
    };

    using Edges = std::vector<idType>; // [index0, index1, overlap, branch,...]
    using Nodes = std::vector<Node>;

    struct CoordinateComparator {
      const float *coordinates;

      CoordinateComparator(const float *coords) : coordinates(coords) {
      }

      inline bool operator()(const size_t &i, const size_t &j) {
        size_t ic = i * 3;
        size_t jc = j * 3;
        return coordinates[ic] == coordinates[jc]
                 ? coordinates[ic + 1] == coordinates[jc + 1]
                     ? coordinates[ic + 2] < coordinates[jc + 2]
                     : coordinates[ic + 1] < coordinates[jc + 1]
                 : coordinates[ic] < coordinates[jc];
      }
    };

    // This function sorts points based on their x, y, and then z coordinate
    int sortCoordinates(const float *pointCoordinates,
                        const size_t nPoints,
                        std::vector<size_t> &sortedIndicies) const {
      printMsg("Sorting coordinates ... ", debug::Priority::PERFORMANCE);
      Timer t;

      sortedIndicies.resize(nPoints);
      for(size_t i = 0; i < nPoints; i++)
        sortedIndicies[i] = i;
      CoordinateComparator c = CoordinateComparator(pointCoordinates);
      sort(sortedIndicies.begin(), sortedIndicies.end(), c);

      std::stringstream msg;
      msg << "done (" << t.getElapsedTime() << " s).";
      printMsg(msg.str(), debug::Priority::PERFORMANCE);

      return 1;
    }

    int computeBranches(std::vector<Edges> &timeEdgesMap,
                        std::vector<Nodes> &timeNodesMap) const {
      printMsg("Computing branches  ... ", debug::Priority::PERFORMANCE);
      Timer tm;

      size_t nT = timeNodesMap.size();

      // Compute max pred and succ
      for(size_t t = 1; t < nT; t++) {
        auto &nodes0 = timeNodesMap[t - 1];
        auto &nodes1 = timeNodesMap[t];
        auto &edges = timeEdgesMap[t - 1];

        size_t nE = edges.size();

        for(size_t i = 0; i < nE; i += 4) {
          auto n0Index = edges[i];
          auto n1Index = edges[i + 1];
          auto &n0 = nodes0[n0Index];
          auto &n1 = nodes1[n1Index];

          sizeType n0MaxSuccSize
            = n0.maxSuccID != -1 ? nodes1[n0.maxSuccID].size : 0;
          sizeType n1MaxPredSize
            = n1.maxPredID != -1 ? nodes0[n1.maxPredID].size : 0;
          if(n0MaxSuccSize < n1.size)
            n0.maxSuccID = n1Index;
          if(n1MaxPredSize < n0.size)
            n1.maxPredID = n0Index;
        }
      }

      // Label first nodes of branches
      idType branchCounter = 0;

      for(size_t t = 0; t < nT; t++)
        for(auto &n : timeNodesMap[t])
          n.branchID = n.maxPredID == -1 ? branchCounter++ : -1;

      for(size_t t = 1; t < nT; t++) {
        auto &nodes0 = timeNodesMap[t - 1];
        auto &nodes1 = timeNodesMap[t];

        for(size_t i = 0; i < nodes1.size(); i++) {
          auto &n1 = nodes1[i];
          if(n1.maxPredID != -1
             && ((idType)i) != nodes0[n1.maxPredID].maxSuccID)
            n1.branchID = branchCounter++;
        }
      }

      // Propagate branch labels
      for(size_t t = 1; t < nT; t++) {
        auto &nodes0 = timeNodesMap[t - 1];
        auto &nodes1 = timeNodesMap[t];
        auto &edges = timeEdgesMap[t - 1];

        size_t nE = edges.size();

        for(size_t i = 0; i < nE; i += 4) {
          auto n0Index = edges[i];
          auto n1Index = edges[i + 1];
          auto &n0 = nodes0[n0Index];
          auto &n1 = nodes1[n1Index];

          if(n1.branchID == -1 && n0Index == n1.maxPredID)
            n1.branchID = n0.branchID;
        }
      }

      // Label edges
      for(size_t t = 1; t < nT; t++) {
        auto &nodes0 = timeNodesMap[t - 1];
        auto &nodes1 = timeNodesMap[t];
        auto &edges = timeEdgesMap[t - 1];

        size_t nE = edges.size();

        for(size_t i = 0; i < nE; i += 4) {
          auto n0Index = edges[i];
          auto n1Index = edges[i + 1];
          auto &n0 = nodes0[n0Index];
          auto &n1 = nodes1[n1Index];

          edges[i + 3] = n0.branchID == n1.branchID ? n0.branchID
                         : n0.maxSuccID == n1Index  ? n0.branchID
                                                    : n1.branchID;
        }
      }

      std::stringstream msg;
      msg << "done (" << tm.getElapsedTime() << " s).";
      printMsg(msg.str(), debug::Priority::PERFORMANCE);

      return 1;
    }

    // This function sorts all unique lables of a point set and then maps these
    // lables to their respective index in the sorted list
    template <typename labelType>
    int computeLabelIndexMap(const labelType *pointLabels,
                             const size_t nPoints,
                             std::map<labelType, size_t> &labelIndexMap) const;

    // This function computes all nodes and their properties based on a labeled
    // point set
    template <typename labelType>
    int computeNodes(const float *pointCoordinates,
                     const labelType *pointLabels,
                     const size_t nPoints,
                     Nodes &nodes) const;

    // This function computes the overlap between two labeled point sets
    template <typename labelType>
    int computeOverlap(const float *pointCoordinates0,
                       const float *pointCoordinates1,
                       const labelType *pointLabels0,
                       const labelType *pointLabels1,
                       const size_t nPoints0,
                       const size_t nPoints1,

                       Edges &edges) const;

  private:
  };
} // namespace ttk

// =============================================================================
// Compute LabelIndexMap
// =============================================================================
template <typename labelType>
int ttk::TrackingFromOverlap::computeLabelIndexMap(
  const labelType *pointLabels,
  const size_t nPoints,
  std::map<labelType, size_t> &labelIndexMap) const {
  for(size_t i = 0; i < nPoints; i++)
    labelIndexMap[pointLabels[i]] = 0;
  size_t i = 0;
  for(auto &it : labelIndexMap)
    it.second = i++;
  return 1;
}

// =============================================================================
// Identify Nodes
// =============================================================================
template <typename labelType>
int ttk::TrackingFromOverlap::computeNodes(const float *pointCoordinates,
                                           const labelType *pointLabels,
                                           const size_t nPoints,
                                           Nodes &nodes) const {
  printMsg("Identifying nodes ..... ", debug::Priority::PERFORMANCE);

  Timer t;

  std::map<labelType, size_t> labelIndexMap;
  this->computeLabelIndexMap(pointLabels, nPoints, labelIndexMap);

  size_t nNodes = labelIndexMap.size();

  nodes.resize(nNodes);
  for(size_t i = 0, q = 0; i < nPoints; i++) {
    labelType label = pointLabels[i];
    Node &n = nodes[labelIndexMap[label]];
    n.label = label;
    n.size++;
    n.x += pointCoordinates[q++];
    n.y += pointCoordinates[q++];
    n.z += pointCoordinates[q++];
  }

  for(size_t i = 0; i < nNodes; i++) {
    Node &n = nodes[i];
    float size = (float)n.size;
    n.x /= size;
    n.y /= size;
    n.z /= size;
  }

  // Print Status
  {
    std::stringstream msg;
    msg << "done (#" << nNodes << " in " << t.getElapsedTime() << " s).";
    printMsg(msg.str(), debug::Priority::PERFORMANCE);
  }

  return 1;
}

// =============================================================================
// Track Nodes
// =============================================================================
template <typename labelType>
int ttk::TrackingFromOverlap::computeOverlap(const float *pointCoordinates0,
                                             const float *pointCoordinates1,
                                             const labelType *pointLabels0,
                                             const labelType *pointLabels1,
                                             const size_t nPoints0,
                                             const size_t nPoints1,

                                             Edges &edges) const {
  // -------------------------------------------------------------------------
  // Compute labelIndexMaps
  // -------------------------------------------------------------------------
  std::map<labelType, size_t> labelIndexMap0;
  std::map<labelType, size_t> labelIndexMap1;
  this->computeLabelIndexMap<labelType>(pointLabels0, nPoints0, labelIndexMap0);
  this->computeLabelIndexMap<labelType>(pointLabels1, nPoints1, labelIndexMap1);

  // -------------------------------------------------------------------------
  // Sort coordinates
  // -------------------------------------------------------------------------
  std::vector<size_t> sortedIndicies0;
  std::vector<size_t> sortedIndicies1;
  this->sortCoordinates(pointCoordinates0, nPoints0, sortedIndicies0);
  this->sortCoordinates(pointCoordinates1, nPoints1, sortedIndicies1);

  // -------------------------------------------------------------------------
  // Track Nodes
  // -------------------------------------------------------------------------
  printMsg("Tracking .............. ", debug::Priority::PERFORMANCE);
  Timer t;

  /* Function that determines configuration of point p0 and p1:
      0: p0Coords = p1Coords
     >0: p0Coords < p1Coords
     <0: p0Coords > p1Coords
  */
  auto compare = [&](size_t p0, size_t p1) {
    size_t p0CoordIndex = p0 * 3;
    size_t p1CoordIndex = p1 * 3;

    float p0_X = pointCoordinates0[p0CoordIndex++];
    float p0_Y = pointCoordinates0[p0CoordIndex++];
    float p0_Z = pointCoordinates0[p0CoordIndex];

    float p1_X = pointCoordinates1[p1CoordIndex++];
    float p1_Y = pointCoordinates1[p1CoordIndex++];
    float p1_Z = pointCoordinates1[p1CoordIndex];

    return p0_X == p1_X  ? p0_Y == p1_Y  ? p0_Z == p1_Z  ? 0
                                           : p0_Z < p1_Z ? -1
                                                         : 1
                           : p0_Y < p1_Y ? -1
                                         : 1
           : p0_X < p1_X ? -1
                         : 1;
  };

  size_t i = 0; // iterator for 0
  size_t j = 0; // iterator for 1

  size_t nEdges = 0;
  std::unordered_map<size_t, std::unordered_map<size_t, size_t>> edgesMap;
  // Iterate over both point sets synchronously using comparison function
  while(i < nPoints0 && j < nPoints1) {
    size_t pointIndex0 = sortedIndicies0[i];
    size_t pointIndex1 = sortedIndicies1[j];

    // Determine point configuration
    int c = compare(pointIndex0, pointIndex1);

    if(c == 0) { // Points have same coordinates -> track
      labelType label0 = pointLabels0[pointIndex0];
      labelType label1 = pointLabels1[pointIndex1];

      size_t &nodeIndex0 = labelIndexMap0[label0];
      size_t &nodeIndex1 = labelIndexMap1[label1];

      // Find edge and increase overlap counter
      auto edges0 = edgesMap.find(nodeIndex0); // Edges from label0 to nodes1

      // If map does not exist then create it
      if(edges0 == edgesMap.end()) {
        edgesMap[nodeIndex0] = std::unordered_map<size_t, size_t>();
        edges0 = edgesMap.find(nodeIndex0);
      }

      // Find edge label0 -> label1
      auto edge = edges0->second.find(nodeIndex1);

      // If edge does not exist then create it
      if(edge == edges0->second.end()) {
        edges0->second[nodeIndex1] = 0;
        edge = edges0->second.find(nodeIndex1);
        nEdges++;
      }

      // Increase overlap
      edge->second++;

      i++;
      j++;
    } else if(c > 0) { // p0 in front of p1 -> let p1 catch up
      j++;
    } else { // p1 in front of p0 -> let p0 catch up
      i++;
    }
  }

  // -------------------------------------------------------------------------
  // Pack Output
  // -------------------------------------------------------------------------
  {
    edges.resize(nEdges * 4);
    size_t q = 0;
    for(auto &it0 : edgesMap) {
      for(auto &it1 : it0.second) {
        edges[q++] = it0.first;
        edges[q++] = it1.first;
        edges[q++] = it1.second;
        edges[q++] = -1;
      }
    }
  }

  // Print Status
  {
    std::stringstream msg;
    msg << "done (#" << nEdges << " in " << t.getElapsedTime() << " s).";
    printMsg(msg.str(), debug::Priority::PERFORMANCE);
  }

  return 0;
}
