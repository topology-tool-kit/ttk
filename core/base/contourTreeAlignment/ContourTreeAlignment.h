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
/// segementationId", a `<scalarType>` array for the vertex scalars and two
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
/// \sa ttk::cta::ContourTree
/// \sa ttk::cta::CTNode
/// \sa ttk::cta::CTEdge
/// \sa ttk::cta::BinaryTree
/// \sa ttk::cta::Tree
/// \sa ttk::cta::AlignmentTree
/// \sa ttk::cta::AlignmentNode
/// \sa ttk::cta::AlignmentEdge
/// \sa ttkContourTreeAlignment
///

#pragma once

// base code includes
#include <CTA_contourtree.h>
#include <Debug.h>

#include <algorithm>
#include <memory>
#include <random>

namespace ttk {

  namespace cta {

    enum Type_Alignmenttree { averageValues, medianValues, lastMatchedValue };
    // enum Type_Match { matchNodes, matchArcs };
    enum Mode_ArcMatch { persistence, area, volume, overlap };

    /**
     * \ingroup base
     * @brief Basic tree data structure for an alignment of two rooted binary
     * trees.
     *
     * This structure represents nodes of a rooted alignment tree.
     *
     * It stores the following meta information for the tree structure:
     * - height of the subtree rooted in this node
     * - size of the subtree rooted in this node
     * - two pointers to its child nodes
     *
     * It stores the following alignment node properties:
     * - pointers to the two matched nodes
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::BinaryTree
     */
    struct AlignmentTree {
      std::shared_ptr<ttk::cta::AlignmentTree> child1;
      std::shared_ptr<ttk::cta::AlignmentTree> child2;
      std::shared_ptr<ttk::cta::BinaryTree> node1;
      std::shared_ptr<ttk::cta::BinaryTree> node2;
      int size;
      int height;
      // int freq;
    };

    struct AlignmentEdge;

    /**
     * \ingroup base
     * @brief Basic data structure for a node of an unrooted alignment tree.
     *
     * This structure represents nodes of an unrooted alignment tree.
     *
     * It stores the following edge properties:
     * - critical type
     * - scalar value
     * - branchID
     *
     * It stores the following information for the tree structure:
     * - a vector with reference ids to its incident edges.
     *
     * It stores the following alignment information:
     * - frequency
     * - references to the ids of the represented nodes of the origial contour
     * trees
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::AlignmentEdge
     */
    struct AlignmentNode {

      ttk::cta::Type_Node type;
      int freq;
      float scalarValue;
      int branchID;

      std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>> edgeList;

      std::vector<std::pair<int, int>> nodeRefs;
    };

    /**
     * \ingroup base
     * @brief Basic data structure for an edge of an unrooted alignment tree.
     *
     * This structure represents edges of an unrooted alignment tree.
     *
     * It stores the following edge properties:
     * - persistence
     * - area
     * - volume
     * - the corresponding segment/region.
     *
     * It stores the following information for the tree structure:
     * - two reference ids to its incident nodes.
     *
     * It stores the following alignment information:
     * - frequency
     * - references to the ids of the represented arcs of the origial contour
     * trees
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::AlignmentNode
     */
    struct AlignmentEdge {

      std::weak_ptr<ttk::cta::AlignmentNode> node1;
      std::weak_ptr<ttk::cta::AlignmentNode> node2;
      float scalardistance;
      float area;
      float volume;
      std::vector<int> region;
      int freq;

      std::vector<std::pair<int, int>> arcRefs;
    };

  } // namespace cta

} // namespace ttk

namespace ttk {

  class ContourTreeAlignment : virtual public Debug {

  public:
    /// Constructor of the Alignment Object
    ContourTreeAlignment() {
      this->setDebugMsgPrefix("ContourTreeAlignment");
    }

    /// Destructor of the Alignment Object
    ~ContourTreeAlignment() {
      contourtrees.clear();
      nodes.clear();
      arcs.clear();
    }

    /// Setter for the matching mode based on arc properties.
    ///
    /// \param mode Determines what arc properties should be compared.
    void setArcMatchMode(int mode) {
      arcMatchMode = static_cast<ttk::cta::Mode_ArcMatch>(mode);
    }
    /// Setter for the weight of combinatorial matching.
    ///
    /// \param weight Determines the factor with which the combinatorial
    /// distance is weighted.
    void setWeightCombinatorialMatch(float weight) {
      weightCombinatorialMatch = weight;
    }
    /// Setter for the weight of arc property matching.
    ///
    /// \param weight Determines the factor with which the arc property based
    /// distance is weighted.
    void setWeightArcMatch(float weight) {
      weightArcMatch = weight;
    }
    /// Setter for the weight of node matching.
    ///
    /// \param weight Determines the factor with which the scalar value distance
    /// is weighted.
    void setWeightScalarValueMatch(float weight) {
      weightScalarValueMatch = weight;
    }
    /// Setter for the type of alignment tree.
    ///
    /// \param type Determines how the labels of the alignment tree are computed
    /// from the matched labels.
    void setAlignmenttreeType(int type) {
      alignmenttreeType = static_cast<ttk::cta::Type_Alignmenttree>(type);
    }

    /// The actual iterated n-tree-alignment algorithm. Computes the alignment
    /// of n input contour trees.
    ///
    /// \param scalarsVP Vector holding n arrays that represent the node scalars
    /// of the n input trees. scalarsVP[i][j] should be the scalar value of the
    /// jth node in the ith tree. \param regionSizes Vector holding n arrays
    /// that represent the region sizes of edges of the n input trees.
    /// regionSizes[i][j] should be the size of the segment associated with the
    /// jth edge in the ith tree. \param segmentationIds Vector holding n arrays
    /// that represent the segmentation ids of the edges the n input trees.
    /// segmentationIds[i][j] should be the segment identifier of the jth node
    /// in the ith tree. \param topologies Vector holding n arrays that
    /// represent the connectivity of the the n input trees. \param nVertices
    /// Vector holding n integers representing the number of nodes of the n
    /// input trees. \param nEdges Vector holding n integers representing the
    /// number of edges of the n input trees. \param segmentations Vector
    /// holding n arrays representing the segmentation arrays of the n input
    /// trees. \param segsizes Vector holding the sizes of the corresponding
    /// segmentations array. segSizes[i] should be the size of the number of
    /// points in the ith scalar field. This should also be the size of
    /// segmentations[i]. \param outputVertices Vector for the alignment node
    /// scalars that will be filled by this algorithm. outputVertices[i] should
    /// be the scalar value of the ith node in the alignment tree. \param
    /// outputFrequencies Vector for the alignment node frequencies that will be
    /// filled by this algorithm. outputFrequencies[i] should be the number of
    /// original nodes matched the ith node in the alignment tree. \param
    /// outputVertexIds Vector for the alignment node matching that will be
    /// filled by this algorithm. outputVertexIds[i*n+j] should be the id of the
    /// vertex from the jth input tree matched in the ith node of the alignment
    /// tree. \param outputBranchIds Vector for the alignment node branch ids
    /// that will be filled by this algorithm. outputBranchIds[i] should be the
    /// branch id of the ith node in the alignment tree. \param
    /// outputSegmentationIds Vector for the alignment edge segmentations that
    /// will be filled by this algorithm. outputSegmentationIds[i*n+j] should be
    /// the id of the segment from the jth input field associated the ith node
    /// of the alignment tree. \param outputArcIds Vector for the alignment edge
    /// maching that will be filled by this algorithm. outputArcIds[i*n+j]
    /// should be the id of the arc from the jth input tree associated the ith
    /// node of the alignment tree. \param outputEdges Vector for the alignment
    /// graph connectivity. The ith edge of the alignment connects the nodes of
    /// index outputEdges[2*i] and outputEdges[2*i+1]. \param seed seed for
    /// randomization.
    template <class scalarType>
    int execute(const std::vector<void *> &scalarsVP,
                const std::vector<int *> &regionSizes,
                const std::vector<int *> &segmentationIds,
                const std::vector<long long *> &topologies,
                const std::vector<size_t> &nVertices,
                const std::vector<size_t> &nEdges,
                const std::vector<int *> &segmentations,
                const std::vector<size_t> &segsizes,

                std::vector<float> &outputVertices,
                std::vector<long long> &outputFrequencies,
                std::vector<long long> &outputVertexIds,
                std::vector<long long> &outputBranchIds,
                std::vector<long long> &outputSegmentationIds,
                std::vector<long long> &outputArcIds,
                std::vector<int> &outputEdges,
                int seed);

    using ContourTree = cta::ContourTree;

    /// This function aligns a new tree to the current alignment.
    /// \pre The given contour tree needs to be binary.
    /// \param t The input contour tree.
    /// \return Return true if alignment was successful, false otherwise.
    bool alignTree(const std::shared_ptr<ContourTree> &t);
    /// This function initializes a new alignment graph from a given contour
    /// tree. \pre The given contour tree needs to be binary. \param t The input
    /// contour tree. \return Return true if initialization was successful,
    /// false otherwise.
    bool initialize(const std::shared_ptr<ContourTree> &t);
    /// This function aligns a new tree to the current alignment but keeps the
    /// old alignment's root. \pre The given contour tree needs to be binary.
    /// \param t The input contour tree.
    /// \return Return true if alignment was successful, false otherwise.
    bool alignTree_consistentRoot(const std::shared_ptr<ContourTree> &t);
    /// This function initializes a new alignment graph from a given contour
    /// tree and sets the fixed root of the alignment to the given input. \pre
    /// The given contour tree needs to be binary. \param t The input contour
    /// tree. \param rootIdx \return Return true if initialization was
    /// successful, false otherwise.
    bool initialize_consistentRoot(const std::shared_ptr<ContourTree> &t,
                                   int rootIdx);

    /// Getter for aligned contour trees
    /// \return Vector of aligned contour trees in a graph representation.
    std::vector<std::pair<std::vector<std::shared_ptr<ttk::cta::CTNode>>,
                          std::vector<std::shared_ptr<ttk::cta::CTEdge>>>>
      getGraphs();
    /// Getter for aligned contour trees
    /// \return Vector of aligned contour trees as ttk::cta::ContourTree data
    /// structures.
    std::vector<std::shared_ptr<ContourTree>> getContourTrees();

    /// Getter for the alignment tree in a graph representation
    /// \return The alignment graph.
    std::pair<std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>,
              std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>>
      getAlignmentGraph();
    /// Getter for the alignment tree in a rooted tree representation
    /// \return The alignment tree rooted in the fixed node determined on
    /// initialization.
    std::shared_ptr<ttk::cta::BinaryTree> getAlignmentGraphRooted();
    /// Getter for the fixed root index.
    /// \return The root index.
    int getAlignmentRootIdx();

    /// Function for aligning two arbitrary binary trees.
    /// \param t1 The first binary rooted tree to align.
    /// \param t2 The second binary rooted tree to align.
    /// \return A pair consisting of the alignment distance and the alignment
    /// tree representing the matching between the two trees.
    std::pair<float, std::shared_ptr<ttk::cta::AlignmentTree>>
      getAlignmentBinary(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                         const std::shared_ptr<ttk::cta::BinaryTree> &t2);

    /// Function that adds branch decomposition information to the alignment
    /// nodes.
    void computeBranches();

  protected:
    // filter parameters
    ttk::cta::Type_Alignmenttree alignmenttreeType = ttk::cta::averageValues;
    ttk::cta::Mode_ArcMatch arcMatchMode = ttk::cta::persistence;
    float weightArcMatch = 1;
    float weightCombinatorialMatch = 0;
    float weightScalarValueMatch = 0;

    // alignment graph data
    std::vector<std::shared_ptr<ttk::cta::AlignmentNode>> nodes;
    std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>> arcs;

    // iteration variables
    std::vector<std::shared_ptr<ContourTree>> contourtrees;
    std::vector<size_t> permutation;
    std::shared_ptr<ttk::cta::AlignmentNode> alignmentRoot;
    int alignmentRootIdx;
    float alignmentVal;

    // functions for aligning two trees (computing the alignment value and
    // memoization matrix)
    float alignTreeBinary(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                          const std::shared_ptr<ttk::cta::BinaryTree> &t2,
                          std::vector<std::vector<float>> &memT,
                          std::vector<std::vector<float>> &memF);
    float alignForestBinary(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                            const std::shared_ptr<ttk::cta::BinaryTree> &t2,
                            std::vector<std::vector<float>> &memT,
                            std::vector<std::vector<float>> &memF);

    // functions for the traceback of the alignment computation (computing the
    // actual alignment tree)
    std::shared_ptr<ttk::cta::AlignmentTree>
      traceAlignmentTree(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                         const std::shared_ptr<ttk::cta::BinaryTree> &t2,
                         std::vector<std::vector<float>> &memT,
                         std::vector<std::vector<float>> &memF);
    std::vector<std::shared_ptr<ttk::cta::AlignmentTree>>
      traceAlignmentForest(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                           const std::shared_ptr<ttk::cta::BinaryTree> &t2,
                           std::vector<std::vector<float>> &memT,
                           std::vector<std::vector<float>> &memF);
    std::shared_ptr<ttk::cta::AlignmentTree>
      traceNullAlignment(const std::shared_ptr<ttk::cta::BinaryTree> &t,
                         bool first);

    // function that defines the local editing costs of two nodes
    float editCost(const std::shared_ptr<ttk::cta::BinaryTree> &t1,
                   const std::shared_ptr<ttk::cta::BinaryTree> &t2);

    // helper functions for tree data structures
    bool isBinary(const std::shared_ptr<ttk::cta::Tree> &t);
    std::shared_ptr<ttk::cta::BinaryTree>
      rootAtNode(const std::shared_ptr<ttk::cta::AlignmentNode> &root);
    std::shared_ptr<ttk::cta::BinaryTree>
      computeRootedTree(const std::shared_ptr<ttk::cta::AlignmentNode> &node,
                        const std::shared_ptr<ttk::cta::AlignmentEdge> &parent,
                        int &id);
    std::shared_ptr<ttk::cta::BinaryTree>
      computeRootedDualTree(const std::shared_ptr<ttk::cta::AlignmentEdge> &arc,
                            bool parent1,
                            int &id);
    void computeNewAlignmenttree(
      const std::shared_ptr<ttk::cta::AlignmentTree> &res);

    // helper functions for branch decomposition
    std::pair<float, std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>>
      pathToMax(const std::shared_ptr<ttk::cta::AlignmentNode> &root,
                const std::shared_ptr<ttk::cta::AlignmentNode> &parent);
    std::pair<float, std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>>
      pathToMin(const std::shared_ptr<ttk::cta::AlignmentNode> &root,
                const std::shared_ptr<ttk::cta::AlignmentNode> &parent);
  };
} // namespace ttk

template <class scalarType>
int ttk::ContourTreeAlignment::execute(
  const std::vector<void *> &scalarsVP,
  const std::vector<int *> &regionSizes,
  const std::vector<int *> &segmentationIds,
  const std::vector<long long *> &topologies,
  const std::vector<size_t> &nVertices,
  const std::vector<size_t> &nEdges,
  const std::vector<int *> &segmentations,
  const std::vector<size_t> &segsizes,
  std::vector<float> &outputVertices,
  std::vector<long long> &outputFrequencies,
  std::vector<long long> &outputVertexIds,
  std::vector<long long> &outputBranchIds,
  std::vector<long long> &outputSegmentationIds,
  std::vector<long long> &outputArcIds,
  std::vector<int> &outputEdges,
  int seed) {

  Timer timer;

  size_t nTrees = nVertices.size();

  std::vector<float *> scalars(nTrees);
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
  nodes = std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>();
  arcs = std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>();
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
  if(alignmenttreeType != ttk::cta::lastMatchedValue)
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

  std::tuple<std::vector<std::shared_ptr<ttk::cta::AlignmentNode>>,
             std::vector<std::shared_ptr<ttk::cta::AlignmentEdge>>,
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

    if(alignmentRoot->type == ttk::cta::saddleNode) {

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
