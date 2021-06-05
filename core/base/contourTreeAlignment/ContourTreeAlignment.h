/// \ingroup base
/// \class ttk::ContourTreeAlignment
/// \author Florian Wetzels (f_wetzels13@cs.uni-kl.de), Jonas Lukasczyk
/// <jl@jluk.de> \date 31.01.2020
///
/// \brief TTK %contourTreeAlignment processing package.
///
/// %ContourTreeAlignment is a TTK processing package that computes an alignment
/// for n contour trees. To compute the alignment, use the execute function.
/// Each contour tree is represented by an integer array for the topology,  an
/// integer array for each of the edge scalars "regionSize" and
/// segementationId", a <scalarType> array for the vertex scalars and two
/// integers for the number of edges and vertices. These properties are passed
/// as vectors of arrays, where the i-th array in a vector represents the
/// corresponding array for the i-th tree. The alignment tree is written to the
/// output vectors which are passed by reference.
///
/// \b Related \b publication: \n
/// 'Fuzzy contour trees: Alignment and joint layout of multiple contour trees'
/// Anna Pia Lohfink, Florian Wetzels, Jonas Lukasczyk, Gunther H. Weber, and
/// Christoph Garth. Comput. Graph. Forum, 39(3):343-355, 2020.
///

#pragma once

// base code includes
#include "contourtree.h"
#include <Debug.h>
#include <algorithm>
#include <memory>
#include <random>

enum Type_Alignmenttree { averageValues, medianValues, lastMatchedValue };
// enum Type_Match { matchNodes, matchArcs };
enum Mode_ArcMatch { persistence, area, volume, overlap };

struct AlignmentTree {
  std::shared_ptr<AlignmentTree> child1;
  std::shared_ptr<AlignmentTree> child2;
  std::shared_ptr<BinaryTree> node1;
  std::shared_ptr<BinaryTree> node2;
  int size;
  int height;
  // int freq;
};

struct AlignmentEdge;

struct AlignmentNode {

  Type_Node type;
  int freq;
  float scalarValue;
  int branchID;

  std::vector<std::shared_ptr<AlignmentEdge>> edgeList;

  std::vector<std::pair<int, int>> nodeRefs;
};

struct AlignmentEdge {

  std::weak_ptr<AlignmentNode> node1;
  std::weak_ptr<AlignmentNode> node2;
  float scalardistance;
  float area;
  float volume;
  std::vector<int> region;
  int freq;

  std::vector<std::pair<int, int>> arcRefs;
};

using namespace std;

namespace ttk {

  class ContourTreeAlignment : virtual public Debug {

  public:
    /// constructor and destructor
    ContourTreeAlignment() {
      this->setDebugMsgPrefix("ContourTreeAlignment");
    };
    ~ContourTreeAlignment() {
      contourtrees.clear();
      nodes.clear();
      arcs.clear();
    };

    /// setter for parameters
    void setArcMatchMode(int mode) {
      arcMatchMode = static_cast<Mode_ArcMatch>(mode);
    }
    void setWeightCombinatorialMatch(float weight) {
      weightCombinatorialMatch = weight;
    }
    void setWeightArcMatch(float weight) {
      weightArcMatch = weight;
    }
    void setWeightScalarValueMatch(float weight) {
      weightScalarValueMatch = weight;
    }
    void setAlignmenttreeType(int type) {
      alignmenttreeType = static_cast<Type_Alignmenttree>(type);
    }

    /// The actual iterated n-tree-alignment algorithm
    template <class scalarType>
    int execute(const vector<void *> &scalarsVP,
                const vector<int *> &regionSizes,
                const vector<int *> &segmentationIds,
                const vector<long long *> &topologies,
                const vector<size_t> &nVertices,
                const vector<size_t> &nEdges,
                const vector<int *> &segmentations,
                const vector<size_t> &segsizes,

                vector<float> &outputVertices,
                vector<long long> &outputFrequencies,
                vector<long long> &outputVertexIds,
                vector<long long> &outputBranchIds,
                vector<long long> &outputSegmentationIds,
                vector<long long> &outputArcIds,
                vector<int> &outputEdges,
                int seed);

    /// functions for aligning single trees in iteration
    bool alignTree(const std::shared_ptr<ContourTree> &t);
    bool initialize(const std::shared_ptr<ContourTree> &t);
    bool alignTree_consistentRoot(const std::shared_ptr<ContourTree> &t);
    bool initialize_consistentRoot(const std::shared_ptr<ContourTree> &t,
                                   int rootIdx);

    /// getters for graph data structures
    std::vector<std::pair<std::vector<std::shared_ptr<CTNode>>,
                          std::vector<std::shared_ptr<CTEdge>>>>
      getGraphs();
    std::vector<std::shared_ptr<ContourTree>> getContourTrees();

    std::pair<std::vector<std::shared_ptr<AlignmentNode>>,
              std::vector<std::shared_ptr<AlignmentEdge>>>
      getAlignmentGraph();
    std::shared_ptr<BinaryTree> getAlignmentGraphRooted();
    int getAlignmentRootIdx();

    /// function for aligning two sarbitrary binary trees
    std::pair<float, std::shared_ptr<AlignmentTree>>
      getAlignmentBinary(const std::shared_ptr<BinaryTree> &t1,
                         const std::shared_ptr<BinaryTree> &t2);

    /// function that adds branch decomposition information
    void computeBranches();

  protected:
    /// filter parameters
    Type_Alignmenttree alignmenttreeType = averageValues;
    Mode_ArcMatch arcMatchMode = persistence;
    float weightArcMatch = 1;
    float weightCombinatorialMatch = 0;
    float weightScalarValueMatch = 0;

    /// alignment graph data
    std::vector<std::shared_ptr<AlignmentNode>> nodes;
    std::vector<std::shared_ptr<AlignmentEdge>> arcs;

    /// iteration variables
    std::vector<std::shared_ptr<ContourTree>> contourtrees;
    std::vector<size_t> permutation;
    std::shared_ptr<AlignmentNode> alignmentRoot;
    int alignmentRootIdx;
    float alignmentVal;

    /// functions for aligning two trees (computing the alignment value and
    /// memoization matrix)
    float alignTreeBinary(const std::shared_ptr<BinaryTree> &t1,
                          const std::shared_ptr<BinaryTree> &t2,
                          std::vector<std::vector<float>> &memT,
                          std::vector<std::vector<float>> &memF);
    float alignForestBinary(const std::shared_ptr<BinaryTree> &t1,
                            const std::shared_ptr<BinaryTree> &t2,
                            std::vector<std::vector<float>> &memT,
                            std::vector<std::vector<float>> &memF);

    /// functions for the traceback of the alignment computation (computing the
    /// actual alignment tree)
    std::shared_ptr<AlignmentTree>
      traceAlignmentTree(const std::shared_ptr<BinaryTree> &t1,
                         const std::shared_ptr<BinaryTree> &t2,
                         std::vector<std::vector<float>> &memT,
                         std::vector<std::vector<float>> &memF);
    std::vector<std::shared_ptr<AlignmentTree>>
      traceAlignmentForest(const std::shared_ptr<BinaryTree> &t1,
                           const std::shared_ptr<BinaryTree> &t2,
                           std::vector<std::vector<float>> &memT,
                           std::vector<std::vector<float>> &memF);
    std::shared_ptr<AlignmentTree>
      traceNullAlignment(const std::shared_ptr<BinaryTree> &t, bool first);

    /// function that defines the local editing costs of two nodes
    float editCost(const std::shared_ptr<BinaryTree> &t1,
                   const std::shared_ptr<BinaryTree> &t2);

    /// helper functions for tree data structures
    bool isBinary(const std::shared_ptr<Tree> &t);
    std::shared_ptr<BinaryTree>
      rootAtNode(const std::shared_ptr<AlignmentNode> &root);
    std::shared_ptr<BinaryTree>
      computeRootedTree(const std::shared_ptr<AlignmentNode> &node,
                        const std::shared_ptr<AlignmentEdge> &parent,
                        int &id);
    std::shared_ptr<BinaryTree> computeRootedDualTree(
      const std::shared_ptr<AlignmentEdge> &arc, bool parent1, int &id);
    void computeNewAlignmenttree(const std::shared_ptr<AlignmentTree> &res);

    /// helper functions for branch decomposition
    std::pair<float, std::vector<std::shared_ptr<AlignmentNode>>>
      pathToMax(const std::shared_ptr<AlignmentNode> &root,
                const std::shared_ptr<AlignmentNode> &parent);
    std::pair<float, std::vector<std::shared_ptr<AlignmentNode>>>
      pathToMin(const std::shared_ptr<AlignmentNode> &root,
                const std::shared_ptr<AlignmentNode> &parent);
  };
} // namespace ttk

template <class scalarType>
int ttk::ContourTreeAlignment::execute(const vector<void *> &scalarsVP,
                                       const vector<int *> &regionSizes,
                                       const vector<int *> &segmentationIds,
                                       const vector<long long *> &topologies,
                                       const vector<size_t> &nVertices,
                                       const vector<size_t> &nEdges,
                                       const vector<int *> &segmentations,
                                       const vector<size_t> &segsizes,
                                       vector<float> &outputVertices,
                                       vector<long long> &outputFrequencies,
                                       vector<long long> &outputVertexIds,
                                       vector<long long> &outputBranchIds,
                                       vector<long long> &outputSegmentationIds,
                                       vector<long long> &outputArcIds,
                                       vector<int> &outputEdges,
                                       int seed) {

  Timer timer;

  size_t nTrees = nVertices.size();

  vector<float *> scalars(nTrees);
  for(size_t t = 0; t < nTrees; t++) {
    scalars[t] = (float *)((scalarType *)scalarsVP[t]);
  }

  // Print Input
  {
    this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
    this->printMsg("Execute base layer");
    this->printMsg("Computing Alignment for " + std::to_string(nTrees)
                   + " trees.");
  }

  for(size_t t = 0; t < nTrees; t++) {
    this->printMsg(ttk::debug::Separator::L2, debug::Priority::VERBOSE);
    this->printMsg("Input Tree " + std::to_string(t) + " topology:",
                   debug::Priority::VERBOSE);

    std::vector<std::vector<std::string>> tableLines;
    tableLines.push_back(
      {"cellId", "vId0", "vId1", "scalar0", "scalar1", "region", "segId"});
    for(size_t i = 0; i < nEdges[t]; i++) {
      long long vertexId0 = topologies[t][i * 2 + 0];
      long long vertexId1 = topologies[t][i * 2 + 1];
      int regionSize = regionSizes[t][i];
      int segmentationId = segmentationIds[t][i];
      scalarType scalarOfVertexId0 = scalars[t][vertexId0];
      scalarType scalarOfVertexId1 = scalars[t][vertexId1];

      std::vector<std::string> tableLine;
      tableLine.push_back(std::to_string(i));
      tableLine.push_back(std::to_string(vertexId0));
      tableLine.push_back(std::to_string(vertexId1));
      tableLine.push_back(std::to_string(scalarOfVertexId0));
      tableLine.push_back(std::to_string(scalarOfVertexId1));
      tableLine.push_back(std::to_string(regionSize));
      tableLine.push_back(std::to_string(segmentationId));

      tableLines.push_back(tableLine);
    }
    this->printMsg(tableLines, debug::Priority::VERBOSE);
  }

  std::vector<std::vector<std::vector<int>>> segRegions;
  if(!segsizes.empty()) {
    for(size_t i = 0; i < nTrees; i++) {
      int maxSegId = -1;
      std::vector<int> seg(segsizes[i]);
      for(size_t j = 0; j < segsizes[i]; j++) {
        maxSegId = std::max(maxSegId, segmentations[i][j]);
        seg[j] = segmentations[i][j];
      }
      std::vector<std::vector<int>> reg(maxSegId + 1);
      for(size_t j = 0; j < segsizes[i]; j++) {
        reg[segmentations[i][j]].push_back(j);
      }
      segRegions.push_back(reg);
    }
    for(size_t i = 0; i < nTrees; i++) {
      int sum = 0;
      for(const auto &seg : segRegions[i]) {
        sum += seg.size();
      }
      this->printMsg("Tree " + std::to_string(i)
                     + ", sum of segment sizes: " + std::to_string(sum));
      sum = 0;
      for(size_t j = 0; j < nEdges[i]; j++) {
        sum += regionSizes[i][j];
      }
      this->printMsg("Tree " + std::to_string(i)
                     + ", sum of region sizes: " + std::to_string(sum));
    }
  }

  // prepare data structures
  contourtrees = std::vector<std::shared_ptr<ContourTree>>();
  nodes = std::vector<std::shared_ptr<AlignmentNode>>();
  arcs = std::vector<std::shared_ptr<AlignmentEdge>>();
  int bestRootIdx{};

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Shuffling input trees. Used seed: " + std::to_string(seed), 0,
                 debug::LineMode::REPLACE);

  // permute list input trees randomly with given seed
  permutation = std::vector<size_t>();
  for(size_t i = 0; i < nTrees; i++) {
    permutation.push_back(i);
  }
  std::srand(seed);
  if(alignmenttreeType != lastMatchedValue)
    std::random_shuffle(permutation.begin(), permutation.end());

  this->printMsg(
    "Shuffling input trees. Used seed: " + std::to_string(seed), 1);

  // print permutation
  if(this->debugLevel_ >= static_cast<int>(debug::Priority::DETAIL)) {
    std::string permutationString = "";
    for(size_t i = 0; i < permutation.size(); i++) {
      permutationString += std::to_string(permutation[i])
                           + (i == permutation.size() - 1 ? "" : ",");
    }
    this->printMsg(
      "Permutation: " + permutationString, debug::Priority::DETAIL);
  }

  this->printMsg("Starting alignment heuristic.");

  std::tuple<std::vector<std::shared_ptr<AlignmentNode>>,
             std::vector<std::shared_ptr<AlignmentEdge>>,
             std::vector<std::shared_ptr<ContourTree>>>
    bestAlignment;
  float bestAlignmentValue = FLT_MAX;

  printMsg(ttk::debug::Separator::L2);
  printMsg("Filtering input contour trees", 0, debug::LineMode::REPLACE);

  std::vector<std::shared_ptr<ContourTree>> contourtreesToAlign;
  for(size_t i = 0; i < nTrees; i++) {
    std::shared_ptr<ContourTree> ct(new ContourTree(
      scalars[permutation[i]], regionSizes[permutation[i]],
      segmentationIds[permutation[i]], topologies[permutation[i]],
      nVertices[permutation[i]], nEdges[permutation[i]],
      segRegions.empty() ? std::vector<std::vector<int>>() : segRegions[i]));
    if(ct->isBinary()) {
      contourtreesToAlign.push_back(ct);
    } else {
      this->printWrn("Input " + std::to_string(permutation[i])
                     + " not binary. Will not be aligned.");
    }
    printMsg("Filtering input contour trees (" + std::to_string(i) + "/"
               + std::to_string(nTrees) + ")",
             (float)i / (float)nTrees, debug::LineMode::REPLACE);
  }
  if(contourtreesToAlign.empty()) {

    this->printErr("No input binary.");

    return 0;
  }

  printMsg("Filtering input contour trees", 1);

  for(size_t rootIdx = 0;
      rootIdx < contourtreesToAlign[0]->getGraph().first.size(); rootIdx++) {

    contourtrees.clear();
    nodes.clear();
    arcs.clear();
    alignmentVal = 0;

    this->printMsg(ttk::debug::Separator::L2);
    this->printMsg("Starting alignment computation with root "
                   + std::to_string(rootIdx));

    // initialize alignment with first tree
    size_t i = 0;

    this->printMsg(
      "Initializing alignment with tree " + std::to_string(permutation[i]), 0,
      debug::LineMode::REPLACE, debug::Priority::DETAIL);
    initialize_consistentRoot(contourtreesToAlign[i], rootIdx);
    this->printMsg(
      "Initializing alignment with tree " + std::to_string(permutation[i]), 1,
      debug::Priority::DETAIL);

    if(alignmentRoot->type == saddleNode) {

      this->printMsg("Initialized root is saddle, alignment aborted.");

      continue;
    }

    // construct other contour tree objects and align them

    i++;
    while(i < contourtreesToAlign.size()) {

      this->printMsg("Aligning tree " + std::to_string(permutation[i]), 0,
                     debug::LineMode::REPLACE, debug::Priority::DETAIL);
      alignTree_consistentRoot(contourtreesToAlign[i]);
      this->printMsg("Aligning tree " + std::to_string(permutation[i]), 1,
                     debug::Priority::DETAIL);

      i++;
    }

    this->printMsg("All trees aligned. Total alignment value: "
                   + std::to_string(alignmentVal));

    if(alignmentVal < bestAlignmentValue) {

      bestAlignmentValue = alignmentVal;
      bestAlignment = std::make_tuple(nodes, arcs, contourtrees);
      bestRootIdx = rootIdx;
    }
  }

  nodes = std::get<0>(bestAlignment);
  arcs = std::get<1>(bestAlignment);
  contourtrees = std::get<2>(bestAlignment);

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Alignment iteration complete.");
  this->printMsg("Root of optimal alignment: " + std::to_string(bestRootIdx)
                 + ".");
  this->printMsg(
    "Value of optimal alignment: " + std::to_string(bestAlignmentValue) + ".");
  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg(
    "Computing branches", 0, timer.getElapsedTime(), debug::LineMode::REPLACE);

  computeBranches();

  this->printMsg("Computing branches", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Writing output pointers", 0, timer.getElapsedTime(),
                 debug::LineMode::REPLACE);

  for(const auto &node : nodes) {

    outputVertices.push_back(node->scalarValue);
    outputFrequencies.push_back(node->freq);
    outputBranchIds.push_back(node->branchID);
    std::vector<long long> refs(nTrees, -1);
    std::vector<long long> segRefs(nTrees, -1);
    for(std::pair<int, int> ref : node->nodeRefs) {
      refs[permutation[ref.first]] = ref.second;
      int eId
        = contourtrees[ref.first]->getGraph().first[ref.second]->edgeList[0];
      segRefs[permutation[ref.first]]
        = contourtrees[ref.first]->getGraph().second[eId]->segId;
    }
    for(int ref : refs) {
      outputVertexIds.push_back(ref);
    }
    for(int segRef : segRefs) {
      outputSegmentationIds.push_back(segRef);
    }
  }

  for(const auto &edge : arcs) {

    int i = 0;
    for(const auto &node : nodes) {

      if(node == edge->node1.lock()) {
        outputEdges.push_back(i);
      }

      if(node == edge->node2.lock()) {
        outputEdges.push_back(i);
      }

      i++;
    }

    std::vector<long long> arcRefs(nTrees, -1);
    for(std::pair<int, int> ref : edge->arcRefs) {
      arcRefs[permutation[ref.first]] = ref.second;
    }
    for(int ref : arcRefs) {
      outputArcIds.push_back(ref);
    }
  }
  this->printMsg("Writing output pointers", 1, timer.getElapsedTime());

  // Print performance
  {
    this->printMsg(ttk::debug::Separator::L1);
    this->printMsg("Alignment computed in "
                   + std::to_string(timer.getElapsedTime()) + " s. ("
                   + std::to_string(threadNumber_) + " thread(s)).");
    this->printMsg("Number of nodes in alignment: "
                   + std::to_string(nodes.size()));
  }

  return 1;
}
