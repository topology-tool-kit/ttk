#ifndef CONTOURTREE_H
#define CONTOURTREE_H

#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <stack>
#include <memory>

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

  // std::vector<CTEdge*> edgeList;
  std::vector<int> edgeList;
};

///=====================================================================================================================
/// edge data structure for an unrooted contour tree (definition)
struct CTEdge {

  int node1Idx;
  int node2Idx;
  float scalardistance;
  float area;
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
              size_t nEdges);
  ~ContourTree();

  std::shared_ptr<BinaryTree> rootAtMax();
  std::shared_ptr<BinaryTree> rootAtNode(CTNode *root);
  bool isBinary();
  void computeBranches();
  std::pair<std::vector<CTNode *>, std::vector<CTEdge *>> getGraph();

  /*static void deleteBinaryTree(std::shared_ptr<BinaryTree> t) {
    if(t->child1)
      deleteBinaryTree(t->child1);
    if(t->child2)
      deleteBinaryTree(t->child2);
    delete t;
  }*/

private:
  std::vector<CTNode *> nodes;
  std::vector<CTEdge *> arcs;

  bool binary;
  float threshold;

  std::shared_ptr<Tree> computeRootedTree(CTNode *node, CTEdge *parent, int &id);
  std::shared_ptr<BinaryTree> computeRootedTree_binary(CTNode *node, CTEdge *parent, int &id);
  std::pair<float, std::vector<int>> pathToMax(int root, int parent);
  std::pair<float, std::vector<int>> pathToMin(int root, int parent);
};

#endif // CONTOURTREE_H
