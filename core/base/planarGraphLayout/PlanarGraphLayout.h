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
/// time indices or scalar values). \b 1) \b Sizes: Points cover space on the
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
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/contourTreeAlignment/">Contour
///   Tree Alignment example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/nestedTrackingFromOverlap/">Nested
///   Tracking from Overlap example</a> \n

#pragma once

#include <map>

// base code includes
#include <Debug.h>

namespace ttk {

  class PlanarGraphLayout : virtual public Debug {

  public:
    PlanarGraphLayout();
    ~PlanarGraphLayout() override;

    template <typename ST, typename IT, typename CT>
    int computeLayout(
      // Output
      float *layout,

      // Input
      const CT *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const ST *pointSequences,
      const float *sizes,
      const IT *branches,
      const IT *levels) const;

    template <typename IT, typename CT>
    int extractLevel(
      // Output
      std::vector<size_t> &nodeIndices,
      std::vector<size_t> &edgeIndices,

      // Input
      const CT *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const IT &level,
      const IT *levels) const;

    template <typename ST, typename IT, typename CT>
    int computeDotString(
      // Output
      std::string &dotString,

      // Input
      const CT *connectivityList,
      const ST *pointSequences,
      const float *sizes,
      const IT *branches,
      const std::vector<size_t> &nodeIndices,
      const std::vector<size_t> &edgeIndices,
      const std::map<ST, size_t> &sequenceValueToIndexMap) const;

    template <typename IT, typename CT>
    int computeSlots(
      // Output
      float *layout,

      // Input
      const CT *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const float *sizes,
      const IT *levels,
      const IT &nLevels) const;

    // Compute Dot Layout
    int computeDotLayout(
      // Output
      float *layout,

      // Input
      const std::vector<size_t> &nodeIndices,
      const std::string &dotString) const;
  };
} // namespace ttk

// =============================================================================
// Extract Level
// =============================================================================
template <typename IT, typename CT>
int ttk::PlanarGraphLayout::extractLevel(
  // Output
  std::vector<size_t> &nodeIndices,
  std::vector<size_t> &edgeIndices,

  // Input
  const CT *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const IT &level,
  const IT *levels) const {

  // If levels==nullptr then return all points and edges
  if(levels == nullptr) {
    nodeIndices.resize(nPoints);
    for(size_t i = 0; i < nPoints; i++)
      nodeIndices[i] = i;

    edgeIndices.resize(nEdges);
    for(size_t i = 0; i < nEdges; i++)
      edgeIndices[i] = i;

    return 1;
  }

  // Get nodes at level
  for(size_t i = 0; i < nPoints; i++)
    if(levels[i] == level)
      nodeIndices.push_back(i);

  // Get edges at level
  size_t const nEdges2 = nEdges * 2;
  for(size_t i = 0; i < nEdges2; i += 2) {
    auto n0l = levels[connectivityList[i + 0]];
    auto n1l = levels[connectivityList[i + 1]];
    if(n0l == level && n0l == n1l)
      edgeIndices.push_back(i / 2);
  }

  return 1;
}

// =============================================================================
// Compute Dot String
// =============================================================================
template <typename ST, typename IT, typename CT>
int ttk::PlanarGraphLayout::computeDotString(
  // Output
  std::string &dotString,

  // Input
  const CT *connectivityList,
  const ST *pointSequences,
  const float *sizes,
  const IT *branches,
  const std::vector<size_t> &nodeIndices,
  const std::vector<size_t> &edgeIndices,
  const std::map<ST, size_t> &sequenceValueToIndexMap) const {

  Timer t;

  this->printMsg("Generating DOT String", 0, debug::LineMode::REPLACE);

  bool const useSequences = pointSequences != nullptr;
  bool const useSizes = sizes != nullptr;
  bool const useBranches = branches != nullptr;

  std::string const headString = "digraph g {rankdir=LR;";
  std::string nodeString = "";
  std::string edgeString = "";
  std::string rankString = "";

  // lambda functions that generate string representations of nodes
  auto sl = [](size_t s) { return "\"s" + std::to_string(s) + "\""; };
  auto nl = [](size_t id) { return std::to_string(id); };

  // ---------------------------------------------------------------------------
  // Nodes
  // ---------------------------------------------------------------------------
  {
    // Set default node style
    nodeString += "node[label=\"\",shape=box,width=1,height=1];";

    // If useSizes then map size to node height
    if(useSizes)
      for(auto &i : nodeIndices)
        nodeString += nl(i) + "[height=" + std::to_string(sizes[i]) + "];";
  }

  // ---------------------------------------------------------------------------
  // Ranks
  // ---------------------------------------------------------------------------
  if(useSequences) {
    size_t const nSequenceValues = sequenceValueToIndexMap.size();

    // Sequence Chain
    {
      edgeString += sl(0);
      for(size_t s = 1; s < nSequenceValues; s++)
        edgeString += "->" + sl(s);
      edgeString += "[weight=1];";
    }

    // Collect nodes with the same sequence index
    std::vector<std::vector<size_t>> sequenceIndexToPointIndexMap(
      nSequenceValues);
    for(auto &i : nodeIndices)
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

  // ---------------------------------------------------------------------------
  // Edges
  // ---------------------------------------------------------------------------
  {
    for(auto &edgeIndex : edgeIndices) {
      size_t const temp = edgeIndex * 2;
      auto &i0 = connectivityList[temp + 0];
      auto &i1 = connectivityList[temp + 1];
      edgeString += nl(i0) + "->" + nl(i1);

      if(useBranches) {
        auto b0 = branches[i0];
        auto b1 = branches[i1];
        edgeString += b0 == b1 ? "[weight=1]" : "[weight=0]";
      }

      edgeString += ";";
    }
  }

  // ---------------------------------------------------------------------------
  // Finalize
  // ---------------------------------------------------------------------------

  // Build Dot String
  { dotString = headString + nodeString + edgeString + rankString + "}"; }

  // Print Status
  this->printMsg("Generating DOT string", 1, t.getElapsedTime());
  this->printMsg("\n" + dotString + "\n", debug::Priority::VERBOSE);

  return 1;
}

// =============================================================================
// Compute Slots
// =============================================================================
template <typename IT, typename CT>
int ttk::PlanarGraphLayout::computeSlots(
  // Output
  float *layout,

  // Input
  const CT *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const float *sizes,
  const IT *levels,
  const IT &nLevels) const {

  if(sizes == nullptr || levels == nullptr) {
    return -1;
  }

  Timer t;
  this->printMsg("Computing slots", 0, debug::LineMode::REPLACE);

  // Comparator that sorts children based on layout.y
  struct ChildrenComparator {
    const float *layout_;

    ChildrenComparator(const float *layout) : layout_(layout) {
    }

    inline bool operator()(const size_t &i, const size_t &j) {
      return layout_[i * 2 + 1] < layout_[j * 2 + 1];
    }
  };

  auto comparator = ChildrenComparator(layout);

  // ---------------------------------------------------------------------------
  // Compute Children
  // ---------------------------------------------------------------------------
  std::vector<std::vector<size_t>> nodeIndexChildrenIndexMap(nPoints);

  size_t const nEdges2 = nEdges * 2;
  for(size_t i = 0; i < nEdges2; i += 2) {
    auto n0 = connectivityList[i + 0];
    auto n1 = connectivityList[i + 1];
    if((levels[n0] + 1) == levels[n1])
      nodeIndexChildrenIndexMap[n0].push_back(n1);
  }

  // ---------------------------------------------------------------------------
  // Adjust positions from bottom to top (skip last level)
  // ---------------------------------------------------------------------------
  for(IT l = 0; l < nLevels - 1; l++) {
    std::vector<size_t> nodeIndices;
    std::vector<size_t> edgeIndices;

    // get nodes at current level (parents)
    this->extractLevel<IT, CT>(
      // Output
      nodeIndices, edgeIndices,

      // Input
      connectivityList, nPoints, nEdges, l, levels);

    // for each parent adjust position of children
    for(auto &parent : nodeIndices) {
      auto &children = nodeIndexChildrenIndexMap[parent];
      size_t const nChildren = children.size();
      if(nChildren < 1)
        continue;

      // sort children
      sort(children.begin(), children.end(), comparator);

      // size of parent
      float const sizeParent = sizes[parent];

      // size of child
      float sizeChildren = 0;
      for(auto &child : children)
        sizeChildren += sizes[child];

      // gap space
      float const gap = sizeParent - sizeChildren;
      float const gapDelta = (gap / (nChildren + 1)) / 2;

      float y = layout[parent * 2 + 1] + sizeParent * 0.5 - gapDelta;
      for(auto &child : children) {
        float const temp = gapDelta + sizes[child] / 2;
        layout[child * 2 + 1] = y - temp;
        y -= 2 * temp;
      }
    }
  }

  this->printMsg("Computing slots", 1, t.getElapsedTime());

  return 1;
}

// =============================================================================
// Execute
// =============================================================================
template <typename ST, typename IT, typename CT>
int ttk::PlanarGraphLayout::computeLayout(
  // Output
  float *layout,

  // Input
  const CT *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const ST *pointSequences,
  const float *sizes,
  const IT *branches,
  const IT *levels) const {

  Timer t;

  // Init Input
  bool const useSequences = pointSequences != nullptr;
  bool const useSizes = sizes != nullptr;
  bool const useBranches = branches != nullptr;
  bool const useLevels = levels != nullptr;

  // Print Input
  {
    std::string modeS = "";
    if(useSequences)
      modeS += "Sequence + ";
    if(useSizes)
      modeS += "Size + ";
    if(useBranches)
      modeS += "Branches + ";
    if(useLevels)
      modeS += "Levels + ";

    this->printMsg(debug::Separator::L1);
    this->printMsg({{"#Nodes", std::to_string(nPoints)},
                    {"#Edges", std::to_string(nEdges)},
                    {"Mode", modeS.substr(0, modeS.length() - 3)}});
    this->printMsg(debug::Separator::L2);
  }

  if(useLevels && !useSizes) {
    this->printErr("'UseLevels' requires 'UseSizes'.");
    return 0;
  }

  // Global SequenceValue to SequenceIndex map
  std::map<ST, size_t> sequenceValueToIndexMap;
  if(useSequences) {
    for(size_t i = 0; i < nPoints; i++)
      sequenceValueToIndexMap[pointSequences[i]] = 0;
    size_t i = 0;
    for(auto &el : sequenceValueToIndexMap)
      el.second = i++;
  }

  // Get number of levels
  IT nLevels = 1;
  if(useLevels) {
    for(size_t i = 0; i < nPoints; i++)
      if(nLevels < levels[i])
        nLevels = levels[i];
    nLevels += 1;
  }

  // ---------------------------------------------------------------------------
  // Compute initial layout for each level
  // ---------------------------------------------------------------------------
  for(IT l = 0; l < nLevels; l++) {
    std::vector<size_t> nodeIndices;
    std::vector<size_t> edgeIndices;

    // Extract nodes and edges at certain level
    {
      int const status = this->extractLevel<IT, CT>(
        // Output
        nodeIndices, edgeIndices,

        // Input
        connectivityList, nPoints, nEdges, l, levels);
      if(status != 1)
        return 0;
    }

    // Compute Dot String
    std::string dotString;
    {
      int const status = this->computeDotString<ST, IT, CT>(
        // Output
        dotString,

        // Input
        connectivityList, pointSequences, sizes, branches, nodeIndices,
        edgeIndices, sequenceValueToIndexMap);
      if(status != 1)
        return 0;
    }

    // Compute Dot Layout
    {
      int const status = this->computeDotLayout(layout, nodeIndices, dotString);
      if(status != 1)
        return 0;
    }
  }

  // ---------------------------------------------------------------------------
  // If nLevels>1 then compute slots
  // ---------------------------------------------------------------------------
  if(nLevels > 1) {
    this->computeSlots<IT, CT>(
      // Output
      layout,

      // Input
      connectivityList, nPoints, nEdges, sizes, levels, nLevels);
  }

  // ---------------------------------------------------------------------------
  // Print performance
  // ---------------------------------------------------------------------------
  this->printMsg(debug::Separator::L2);
  this->printMsg("Complete", 1, t.getElapsedTime());
  this->printMsg(debug::Separator::L1);

  return 1;
}
