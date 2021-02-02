#ifndef CONTOURTREE_H
#define CONTOURTREE_H

#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <stack>

///=====================================================================================================================
/// enums and structs for tree data structure
///=====================================================================================================================

///=====================================================================================================================
/// node types of a contour tree
enum Type_Node { minNode, maxNode, saddleNode };

///=====================================================================================================================
/// basic tree data structure for a rooted contour tree of unbounded degree
struct Tree {
  std::vector<std::shared_ptr<Tree>> children;
  Type_Node type;
  int vertexId;
  int id;
  int size;
  int height;
  float scalardistanceParent;
  float volume;
};

///=====================================================================================================================
/// basic tree data structure for a rooted contour tree of degree 2
struct BinaryTree {
  std::shared_ptr<BinaryTree> child1;
  std::shared_ptr<BinaryTree> child2;
  Type_Node type;
  int vertexId;
  int id;
  int size;
  int height;
  float scalardistanceParent;
  float area;
  float volume;
  float scalarValue;
  std::vector<int> region;

  // for fuzzy trees
  int freq;
  std::vector<std::pair<int, int>> nodeRefs;
  std::vector<std::pair<int, int>> arcRefs;
};

///=====================================================================================================================
/// edge data structure for an unrooted contour tree (forward declaration)
struct CTEdge;

///=====================================================================================================================
/// node data structure for an unrooted contour tree
struct CTNode {

  Type_Node type;
  float scalarValue;
  int branchID;

  std::vector<int> edgeList;
};

///=====================================================================================================================
/// edge data structure for an unrooted contour tree (definition)
struct CTEdge {

  int node1Idx;
  int node2Idx;
  float scalardistance;
  float area;
  std::vector<int> region;
  float volume;
  int segId;
};

///=======================================================================================S==============================
/// class for an unrooted contour tree
///=====================================================================================================================

class ContourTree {

public:
  ContourTree(float *scalars,
              int *regionSizes,
              int *segmentationIds,
              long long *topology,
              size_t nVertices,
              size_t nEdges,
              std::vector<std::vector<int>> regions
              = std::vector<std::vector<int>>());
  ~ContourTree();

  std::shared_ptr<BinaryTree> rootAtMax();
  std::shared_ptr<BinaryTree> rootAtNode(const std::shared_ptr<CTNode> &root);
  bool isBinary();
  void computeBranches();
  std::pair<std::vector<std::shared_ptr<CTNode>>,
            std::vector<std::shared_ptr<CTEdge>>>
    getGraph();

private:
  std::vector<std::shared_ptr<CTNode>> nodes;
  std::vector<std::shared_ptr<CTEdge>> arcs;

  bool binary;

  std::shared_ptr<Tree> computeRootedTree(const std::shared_ptr<CTNode> &node,
                                          const std::shared_ptr<CTEdge> &parent,
                                          int &id);
  std::shared_ptr<BinaryTree>
    computeRootedTree_binary(const std::shared_ptr<CTNode> &node,
                             const std::shared_ptr<CTEdge> &parent,
                             int &id);
  std::pair<float, std::vector<int>> pathToMax(int root, int parent);
  std::pair<float, std::vector<int>> pathToMin(int root, int parent);
};

#endif // CONTOURTREE_H
