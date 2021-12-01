/// \ingroup base
/// \class Octree
/// \author Guoxi Liu (guoxil@g.clemson.edu)
/// \date May 2021.
///
/// \brief Implementation of the point region (PR) octree.
///
/// \b Related \b publications \n
/// "The PR-star octree: A spatio-topological data structure for tetrahedral
/// meshes." Kenneth Weiss, Leila Floriani, Riccardo Fellegara, and Marcelo
/// Velloso In Proceedings of the 19th ACM SIGSPATIAL International Conference
/// on Advances in Geographic Information Systems, 2011.
///
/// \sa ttk::CompactTriangulationPreconditioning

#pragma once

#include <Triangulation.h>
#include <assert.h>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <vector>

class OctreeNode {
public:
  OctreeNode() {
    locCode_ = 0;
    childExists_ = 0;
  }
  OctreeNode(uint32_t location) {
    locCode_ = location;
    childExists_ = 0;
  }

  ~OctreeNode() {
  }

protected:
  uint32_t locCode_;
  uint8_t childExists_;
  std::vector<ttk::SimplexId> vertexIds_;
  std::vector<ttk::SimplexId> cellIds_;

  friend class Octree;
};

class Octree : public virtual ttk::Debug {
public:
  Octree(const ttk::AbstractTriangulation *t);
  Octree(const ttk::AbstractTriangulation *t, const int k);
  ~Octree();
  void initialize(const ttk::AbstractTriangulation *t, const int k);

  bool empty();

  size_t getNodeTreeDepth(const OctreeNode *node);

  OctreeNode *getParentNode(OctreeNode *node);

  void visitAll(const OctreeNode *node);

  int verifyTree(ttk::SimplexId &vertexNum);

  int insertVertex(ttk::SimplexId &vertexId);

  int insertCell(ttk::SimplexId &cellId);

  void reindex(std::vector<ttk::SimplexId> &vertices,
               std::vector<ttk::SimplexId> &nodes,
               std::vector<ttk::SimplexId> &cells);

private:
  const ttk::AbstractTriangulation *triangulation_;
  std::unordered_map<uint32_t, OctreeNode> allNodes_;
  int capacity_;
  std::array<float, 3> center_{}, size_{};

  /**
   * Loop up the octree node in the unordered map given the location code.
   * Returns nullptr if not found.
   */
  inline OctreeNode *lookupNode(uint32_t locCode) {
    const auto iter = allNodes_.find(locCode);
    return (iter == allNodes_.end() ? nullptr : &(iter->second));
  }

  /**
   * Compute the center and size with given location code.
   * Note: initialize two arrays before calling this function!
   */
  void computeCenterSize(uint32_t location,
                         std::array<float, 3> &centerArr,
                         std::array<float, 3> &sizeArr);

  /**
   * Get the location code of the child.
   */
  uint32_t getChildLocation(uint32_t parLoc,
                            ttk::SimplexId vertexId,
                            const std::array<float, 3> &centerArr);

  /**
   * Recursive subdividing function to make sure each node does not exceed the
   * capacity limit.
   */
  void subdivide(OctreeNode *node);
};