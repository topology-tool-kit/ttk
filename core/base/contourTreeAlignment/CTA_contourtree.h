#pragma once

#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <stack>

namespace ttk {

  namespace cta {

    //#####################################################################################################################
    // enums and structs for tree data structure

    //=====================================================================================================================
    // node types of a contour tree
    enum Type_Node { minNode, maxNode, saddleNode };

    //=====================================================================================================================
    // Tree Data Structure with arbitrary degree

    /**
     * \ingroup base
     * @brief Basic tree data structure for a rooted contour tree of unbounded
     * degree.
     *
     * This structure represents nodes of a rooted contour tree.
     * It stores the following node properties:
     * - critical type
     * - vertex id
     *
     * It stores the following edge properties (these always refer to its parent
     * edge and are filled with dummy values for the root vertex):
     * - persistence
     * - volume
     *
     * It stores the following meta information for the tree structure:
     * - height of the subtree rooted in this node
     * - size of the subtree rooted in this node
     * - node id
     * - a vector of pointers to its child nodes.
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::ContourTree
     */
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

    //=====================================================================================================================
    // Binary Tree Data Structure

    /**
     * \ingroup base
     * @brief Basic tree data structure for a rooted contour tree of degree 2
     * and rooted contour tree alignments of degree 2.
     *
     * This structure represents nodes of a rooted contour tree.
     *
     * It stores the following node properties:
     * - critical type
     * - vertex id
     *
     * It stores the following edge properties (these always refer to its parent
     * edge and are filled with dummy values for the root vertex):
     * - persistence
     * - area
     * - volume
     * - the corresponding segment/region
     *
     * It stores the following meta information for the tree structure:
     * - height of the subtree rooted in this node
     * - size of the subtree rooted in this node
     * - node id and two pointers to its child nodes
     *
     * It stores the following alignment node properties:
     * - frequency of the node
     * - a vector of references to matched contour tree nodes
     * - a vector of references to matched contour tree arcs.
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::ContourTree
     * \sa ttk::cta::AlignmentTree
     */
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

    //=====================================================================================================================
    // Structs for unrooted Contour Tree

    struct CTEdge;

    /**
     * \ingroup base
     * @brief Basic data structure for a node of an unrooted contour tree.
     *
     * This structure represents nodes of an unrooted contour tree.
     *
     * It stores the following node properties:
     * - critical type
     * - scalar value
     * - branch id
     *
     * It stores the following information for the tree structure:
     * - a vector of reference ids to its incident edges.
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::ContourTree
     * \sa ttk::cta::CTEdge
     */
    struct CTNode {

      Type_Node type;
      float scalarValue;
      int branchID;

      std::vector<int> edgeList;
    };

    /**
     * \ingroup base
     * @brief Basic data structure for an edge of an unrooted contour tree.
     *
     * This structure represents edges of an unrooted contour tree.
     *
     * It stores the following edge properties:
     * - persistence
     * - area
     * - volume
     * - segmentation id
     * - the corresponding segment/region.
     *
     * It stores the following information for the tree structure:
     * - two reference ids to its incident nodes.
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::ContourTree
     * \sa ttk::cta::CTNode
     */
    struct CTEdge {

      int node1Idx;
      int node2Idx;
      float scalardistance;
      float area;
      std::vector<int> region;
      float volume;
      int segId;
    };

    ///#####################################################################################################################
    // class for an unrooted contour tree

    /**
     * \ingroup base
     * @brief Contour %Tree Data Structure for an unrooted contour tree of
     * unbounded degree for internal use from the ttk:ContourTreeAlignment
     * module.
     *
     * \author Florian Wetzels <f_wetzels13@cs.uni-kl.de> \date 26.11.2021
     *
     * ToDo
     *
     * \sa ttk::ContourTreeAlignment
     * \sa ttk::cta::CTNode
     * \sa ttk::cta::CTEdge
     * \sa ttk::cta::BinaryTree
     * \sa ttk::cta::Tree
     */
    class ContourTree {

    public:
      /// Constructor for internal contour tree class. Constructs graph from vtk
      /// style array input.
      ///
      /// \param scalars Scalar values of the contour tree nodes. scalars[i] is
      /// scalar values of the node of index i. \param regionSizes Region size
      /// values of the contour tree edges. regionSizes[i] is the number of
      /// vertices in the segment associated with the edge of index i. \param
      /// segmentationIds Segmentation ids of the contour tree edges.
      /// segmentationIds[i] is an identifier of the segment associated with the
      /// edge of index i. \param topology Connectivity information describing
      /// the graph stucture. topology[i][0] and topology[i][1] are the indices
      /// of the two nodes incident to the edge of index i. \param nVertices The
      /// number of vertices in the graph and the length of scalar array. \param
      /// nEdges The number of edges in the graph and the length of regionSizes,
      /// segmentationIds and topology arrays. \param regions Segments/Regions
      /// associated with the edges of the contour tree. regions[i] is a vector
      /// of vertex identifiers that belong to the edge of index i.
      ContourTree(float *scalars,
                  int *regionSizes,
                  int *segmentationIds,
                  long long *topology,
                  size_t nVertices,
                  size_t nEdges,
                  std::vector<std::vector<int>> regions = {});

      /// Destructor of internal contour tree class.
      ~ContourTree();

      /// Get a rooted binary representation of the contour tree. The root is
      /// the vertex of highest scalar value.
      ///
      /// \pre The contour tree should be binary. If it is not, the rooted tree
      /// will not be complete. \returns Smartpointer to the binary tree object.
      std::shared_ptr<BinaryTree> rootAtMax();

      /// Get a rooted binary representation of the contour tree. The root is
      /// the vertex passed as first argument.
      ///
      /// \pre The contour tree should be binary. If it is not, the rooted tree
      /// will not be complete. \param root The contour tree node that will be
      /// the root of the rooted tree. \returns Smartpointer to the binary tree
      /// object.
      std::shared_ptr<BinaryTree>
        rootAtNode(const std::shared_ptr<CTNode> &root);

      /// Checks if the maximum node degree is at most 3.
      ///
      /// \returns True if the contour tree is binary, False otherwise.
      bool isBinary();

      /// Computes a branch decomposition of the contour tree and attaches the
      /// branch information to the node objects.
      void computeBranches();

      /// Returns a graph representaiton of the contour tree as node list and
      /// edge list.
      ///
      /// \returns Pair consisting of the node list and the edge list.
      std::pair<std::vector<std::shared_ptr<CTNode>>,
                std::vector<std::shared_ptr<CTEdge>>>
        getGraph();

    private:
      std::vector<std::shared_ptr<CTNode>> nodes;
      std::vector<std::shared_ptr<CTEdge>> arcs;

      bool binary;

      std::shared_ptr<Tree>
        computeRootedTree(const std::shared_ptr<CTNode> &node,
                          const std::shared_ptr<CTEdge> &parent,
                          int &id);
      std::shared_ptr<BinaryTree>
        computeRootedTree_binary(const std::shared_ptr<CTNode> &node,
                                 const std::shared_ptr<CTEdge> &parent,
                                 int &id);
      std::pair<float, std::vector<int>> pathToMax(int root, int parent);
      std::pair<float, std::vector<int>> pathToMin(int root, int parent);
    };

  }; // namespace cta
}; // namespace ttk
