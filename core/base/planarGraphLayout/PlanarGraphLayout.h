/// \ingroup base
/// \class ttk::PlanarGraphLayout
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK %planarGraphLayout processing package.
///
/// %PlanarGraphLayout is a TTK processing package that computes a planar graph
/// layout of a \b vtkUnstructuredGrid. To improve the quality of the layout it
/// is possible to pass additional field data to the algorithm:\n \b 1) \b
/// Sequences: Points are positioned along the x-axis based on a sequence (e.g.,
/// time indicies or scalar values). \b 1) \b Sizes: Points cover space on the
/// y-axis based on their size. \b 1) \b Branches: Points with the same branch
/// label are positioned on straight lines. \b 1) \b Levels: The layout of
/// points with the same level label are computed individually and afterwards
/// nested based on the level hierarchy. This makes it possible to draw nested
/// graphs where each level is a layer of the resulting graph.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///

#pragma once

#include <map>

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk {

  class PlanarGraphLayout : public Debug {

  public:
    PlanarGraphLayout(){};
    ~PlanarGraphLayout(){};

    template <typename topoType, typename idType, typename sequenceType>
    int execute(
      // Input
      const sequenceType *pointSequences,
      const float *sizes,
      const idType *branches,
      const idType *levels,
      const topoType *topology,
      const size_t &nPoints,
      const size_t &nEdges,

      // Output
      float *layout) const;

    template <typename topoType, typename idType>
    int extractLevel(
      // Input
      const idType &level,
      const idType *levels,
      const topoType *topology,
      const size_t &nPoints,
      const size_t &nEdges,

      // Output
      vector<size_t> &nodeIndicies,
      vector<size_t> &edgeIndicies) const;

    template <typename topoType, typename idType, typename sequenceType>
    int computeDotString(
      // Input
      const sequenceType *pointSequences,
      const float *sizes,
      const idType *branches,
      const topoType *topology,
      const vector<size_t> &nodeIndicies,
      const vector<size_t> &edgeIndicies,
      const map<sequenceType, size_t> &sequenceValueToIndexMap,

      // Output
      string &dotString) const;

    template <typename topoType, typename idType>
    int computeSlots(
      // Input
      const float *sizes,
      const idType *levels,
      const topoType *topology,
      const size_t &nPoints,
      const size_t &nEdges,
      const idType &nLevels,

      // Output
      float *layout) const;

    // Compute Dot Layout
    int computeDotLayout(
      // Input
      const vector<size_t> &nodeIndicies,
      const string &dotString,

      // Output
      float *layout) const;
  };
} // namespace ttk

// =============================================================================
// Extract Level
// =============================================================================
template <typename topoType, typename idType>
int ttk::PlanarGraphLayout::extractLevel(
  // Input
  const idType &level,
  const idType *levels,
  const topoType *topology,
  const size_t &nPoints,
  const size_t &nEdges,

  // Output
  vector<size_t> &nodeIndicies,
  vector<size_t> &edgeIndicies) const {

  // If levels==nullptr then return all points and edges
  if(levels == nullptr) {
    nodeIndicies.resize(nPoints);
    for(size_t i = 0; i < nPoints; i++)
      nodeIndicies[i] = i;

    edgeIndicies.resize(nEdges);
    for(size_t i = 0; i < nEdges; i++)
      edgeIndicies[i] = i;

    return 1;
  }

  // Get nodes at level
  for(size_t i = 0; i < nPoints; i++)
    if(levels[i] == level)
      nodeIndicies.push_back(i);

  // Get edges at level
  size_t nEdges3 = nEdges * 3;
  for(size_t i = 0; i < nEdges3; i += 3) {
    auto n0l = levels[topology[i + 1]];
    auto n1l = levels[topology[i + 2]];
    if(n0l == level && n0l == n1l)
      edgeIndicies.push_back(i / 3);
  }

  return 1;
}

// =============================================================================
// Compute Dot String
// =============================================================================
template <typename topoType, typename idType, typename sequenceType>
int ttk::PlanarGraphLayout::computeDotString(
  const sequenceType *pointSequences,
  const float *sizes,
  const idType *branches,
  const topoType *topology,
  const vector<size_t> &nodeIndicies,
  const vector<size_t> &edgeIndicies,
  const map<sequenceType, size_t> &sequenceValueToIndexMap,

  string &dotString) const {

  Timer t;

  bool useSequences = pointSequences != nullptr;
  bool useSizes = sizes != nullptr;
  bool useBranches = branches != nullptr;

  string headString = "digraph g {rankdir=LR;";
  string nodeString = "";
  string edgeString = "";
  string rankString = "";

  // lambda functions that generate string representations of nodes
  auto sl = [](size_t s) { return "\"s" + to_string(s) + "\""; };
  auto nl = [](size_t id) { return to_string(id); };

  // -------------------------------------------------------------------------
  // Nodes
  // -------------------------------------------------------------------------
  {
    // Set default node style
    nodeString += "node[label=\"\",shape=box,width=1,height=1];";

    // If useSizes then map size to node height
    if(useSizes)
      for(auto &i : nodeIndicies)
        nodeString += nl(i) + "[height=" + to_string(sizes[i]) + "];";
  }

  // -------------------------------------------------------------------------
  // Ranks
  // -------------------------------------------------------------------------
  if(useSequences) {
    size_t nSequenceValues = sequenceValueToIndexMap.size();

    // Sequence Chain
    {
      edgeString += sl(0);
      for(size_t s = 1; s < nSequenceValues; s++)
        edgeString += "->" + sl(s);
      edgeString += "[weight=1];";
    }

    // Collect nodes with the same sequence index
    vector<vector<size_t>> sequenceIndexToPointIndexMap(nSequenceValues);
    for(auto &i : nodeIndicies)
      sequenceIndexToPointIndexMap
        [sequenceValueToIndexMap.find(pointSequences[i])->second]
          .push_back(i);

    // Compute individual ranks
    for(size_t s = 0; s < nSequenceValues; s++) {
      rankString += "{rank=same " + sl(s);

      auto &nodes = sequenceIndexToPointIndexMap[s];
      for(auto &i : nodes)
        rankString += " " + nl(i);

      rankString += "}";
    }
  }

  // -------------------------------------------------------------------------
  // Edges
  // -------------------------------------------------------------------------
  {
    for(auto &edgeIndex : edgeIndicies) {
      size_t temp = edgeIndex * 3;
      auto &i0 = topology[temp + 1];
      auto &i1 = topology[temp + 2];
      edgeString += nl(i0) + "->" + nl(i1);

      if(useBranches) {
        auto b0 = branches[i0];
        auto b1 = branches[i1];
        edgeString += b0 == b1 ? "[weight=1]" : "[weight=0]";
      }

      edgeString += ";";
    }
  }

  // -------------------------------------------------------------------------
  // Build Dot String
  // -------------------------------------------------------------------------
  { dotString = headString + nodeString + edgeString + rankString + "}"; }

  // -------------------------------------------------------------------------
  // Print Status
  // -------------------------------------------------------------------------
  {
    stringstream msg;
    msg << "[ttkPlanarGraphLayout] Dot String generated in "
        << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);

    dMsg(cout, "\n" + dotString + "\n\n", advancedInfoMsg);
  }

  return 1;
}

// =============================================================================
// Compute Slots
// =============================================================================
template <typename topoType, typename idType>
int ttk::PlanarGraphLayout::computeSlots(
  // Input
  const float *sizes,
  const idType *levels,
  const topoType *topology,
  const size_t &nPoints,
  const size_t &nEdges,
  const idType &nLevels,

  // Output
  float *layout) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(sizes == nullptr || levels == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  // Comparator that sorts children based on layout.y
  struct ChildrenComparator {
    const float *layout_;

    ChildrenComparator(const float *layout) : layout_(layout){};

    inline bool operator()(const size_t &i, const size_t &j) {
      return layout_[i * 2 + 1] < layout_[j * 2 + 1];
    }
  };

  auto comparator = ChildrenComparator(layout);

  // -------------------------------------------------------------------------
  // Compute Children
  // -------------------------------------------------------------------------
  vector<vector<size_t>> nodeIndexChildrenIndexMap(nPoints);

  size_t nEdges3 = nEdges * 3;
  for(size_t i = 0; i < nEdges3; i += 3) {
    auto n0 = topology[i + 1];
    auto n1 = topology[i + 2];
    if((levels[n0] + 1) == levels[n1])
      nodeIndexChildrenIndexMap[n0].push_back(n1);
  }

  // -------------------------------------------------------------------------
  // Adjust positions from bottom to top (skip last level)
  // -------------------------------------------------------------------------
  for(idType l = 0; l < nLevels - 1; l++) {
    vector<size_t> nodeIndicies;
    vector<size_t> edgeIndicies;

    // get nodes at current level (parents)
    this->extractLevel<topoType, idType>(
      // Input
      l, levels, topology, nPoints, nEdges,

      // Output
      nodeIndicies, edgeIndicies);

    // for each parent adjust position of children
    for(auto &parent : nodeIndicies) {
      auto &children = nodeIndexChildrenIndexMap[parent];
      size_t nChildren = children.size();
      if(nChildren < 1)
        continue;

      // sort children
      sort(children.begin(), children.end(), comparator);

      // size of parent
      float sizeParent = sizes[parent];

      // size of child
      float sizeChildren = 0;
      for(auto &child : children)
        sizeChildren += sizes[child];

      // gap space
      float gap = sizeParent - sizeChildren;
      float gapDelta = (gap / (nChildren + 1)) / 2;

      float y = layout[parent * 2 + 1] + sizeParent * 0.5 - gapDelta;
      for(auto &child : children) {
        float temp = gapDelta + sizes[child] / 2;
        layout[child * 2 + 1] = y - temp;
        y -= 2 * temp;
      }
    }
  }

  return 1;
}

// =============================================================================
// Execute
// =============================================================================
template <typename topoType, typename idType, typename sequenceType>
int ttk::PlanarGraphLayout::execute(
  // Input
  const sequenceType *pointSequences,
  const float *sizes,
  const idType *branches,
  const idType *levels,
  const topoType *topology,
  const size_t &nPoints,
  const size_t &nEdges,

  // Output
  float *layout) const {

  Timer t;

  // Init Input
  bool useSequences = pointSequences != nullptr;
  bool useSizes = sizes != nullptr;
  bool useBranches = branches != nullptr;
  bool useLevels = levels != nullptr;

  // Print Input
  {
    stringstream msg;
    msg << "[ttkPlanarGraphLayout] Computing layout for graph with" << endl
        << "[ttkPlanarGraphLayout]  - " << nPoints << " vertices" << endl
        << "[ttkPlanarGraphLayout]  - " << nEdges << " edges" << endl;
    if(useSequences)
      msg << "[ttkPlanarGraphLayout]  - using sequences" << endl;
    if(useSizes)
      msg << "[ttkPlanarGraphLayout]  - using sizes" << endl;
    if(useBranches)
      msg << "[ttkPlanarGraphLayout]  - using branches" << endl;
    if(useLevels)
      msg << "[ttkPlanarGraphLayout]  - using levels" << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  if(useLevels && !useSizes) {
    dMsg(cout,
         "[ttkPlanarGraphLayout] ERROR: When 'UseLevels' is enabled then "
         "'UseSizes' must also be enabled.\n",
         fatalMsg);
    return 0;
  }

  // Global SequenceValue to SequenceIndex map
  map<sequenceType, size_t> sequenceValueToIndexMap;
  if(useSequences) {
    for(size_t i = 0; i < nPoints; i++)
      sequenceValueToIndexMap[pointSequences[i]] = 0;
    size_t i = 0;
    for(auto &el : sequenceValueToIndexMap)
      el.second = i++;
  }

  // Get number of levels
  idType nLevels = 1;
  if(useLevels) {
    for(size_t i = 0; i < nPoints; i++)
      if(nLevels < levels[i])
        nLevels = levels[i];
    nLevels += 1;
  }

  // -------------------------------------------------------------------------
  // Compute initial layout for each level
  // -------------------------------------------------------------------------
  for(idType l = 0; l < nLevels; l++) {
    vector<size_t> nodeIndicies;
    vector<size_t> edgeIndicies;

    // Extract nodes and edges at certain level
    {
      int status = this->extractLevel<topoType, idType>(
        // Input
        l, levels, topology, nPoints, nEdges,

        // Output
        nodeIndicies, edgeIndicies);
      if(status != 1)
        return 0;
    }

    // Compute Dot String
    string dotString;
    {
      int status = this->computeDotString<topoType, idType, sequenceType>(
        // Input
        pointSequences, sizes, branches, topology, nodeIndicies, edgeIndicies,
        sequenceValueToIndexMap,

        // Output
        dotString);
      if(status != 1)
        return 0;
    }

    // Compute Dot Layout
    {
      int status = this->computeDotLayout(nodeIndicies, dotString, layout);
      if(status != 1)
        return 0;
    }
  }

  // -------------------------------------------------------------------------
  // If nLevels>1 then compute slots
  // -------------------------------------------------------------------------
  if(nLevels > 1) {
    this->computeSlots<topoType, idType>(
      // Input
      sizes, levels, topology, nPoints, nEdges, nLevels,

      // Output
      layout);
  }

  // -------------------------------------------------------------------------
  // Print performance
  // -------------------------------------------------------------------------
  {
    stringstream msg;
    msg << "[ttkPlanarGraphLayout] Layout computed in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 1;
}
