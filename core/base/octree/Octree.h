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

using namespace ttk;
using namespace std;

class OctreeNode {
public:
  OctreeNode() {
    OctreeNode(0);
  }
  OctreeNode(uint32_t location) {
    locCode = location;
    childExists = 0;
    vertexIds = new vector<SimplexId>();
    cellIds = nullptr;
  }

  ~OctreeNode() {
  }

protected:
  uint32_t locCode;
  uint8_t childExists;
  vector<SimplexId> *vertexIds;
  vector<SimplexId> *cellIds;

  friend class Octree;
};

class Octree : public virtual Debug {
public:
  Octree(const AbstractTriangulation *t);
  Octree(const AbstractTriangulation *t, const int k);
  ~Octree();

  bool empty();

  size_t getNodeTreeDepth(const OctreeNode *node);

  OctreeNode *getParentNode(OctreeNode *node);

  void visitAll(const OctreeNode *node);

  int verifyTree(SimplexId &vertexNum);

  int insertVertex(SimplexId &vertexId);

  int insertCell(SimplexId &cellId);

  void reindex(vector<SimplexId> &vertices,
               vector<SimplexId> &nodes,
               vector<SimplexId> &cells);

private:
  const AbstractTriangulation *triangulation_;
  unordered_map<uint32_t, OctreeNode> allNodes;
  int capacity;
  float center[3];
  float size[3];

  /**
   * Loop up the octree node in the unordered map given the location code.
   * Returns nullptr if not found.
   */
  OctreeNode *lookupNode(uint32_t locCode) {
    const auto iter = allNodes.find(locCode);
    return (iter == allNodes.end() ? nullptr : &(iter->second));
  }

  /**
   * Compute the center and size with given location code.
   * Note: initialize two arrays before calling this function!
   */
  void
    computeCenterSize(uint32_t location, float centerArr[3], float sizeArr[3]) {
    int leftmost = 0;
    uint32_t tmp = location;
    while(tmp >>= 1) {
      leftmost++;
    }
    if(leftmost % 3 != 0) {
      this->printMsg("Location: " + to_string(location)
                     + ", leftmost: " + to_string(leftmost));
      this->printErr("computeCenterSize(): the location seems not correct!");
      this->printErr("Please try a larger bucket capacity!");
      return;
    }

    // initialize the center and size arrays
    for(int i = 0; i < 3; i++) {
      centerArr[i] = center[i];
      sizeArr[i] = size[i];
    }

    uint32_t mask;
    for(int i = leftmost - 1; i >= 0; i -= 3) {
      sizeArr[0] /= 2.0;
      sizeArr[1] /= 2.0;
      sizeArr[2] /= 2.0;
      for(int j = 0; j < 3; j++) {
        mask = 1 << (i - j);
        if(location & mask) {
          centerArr[j] += sizeArr[j];
        } else {
          centerArr[j] -= sizeArr[j];
        }
      }
    }
  }

  /**
   * Get the location code of the child.
   */
  uint32_t
    getChildLocation(uint32_t parLoc, SimplexId vertexId, float centerArr[3]) {
    float xval, yval, zval;
    if(triangulation_->getVertexPoint(vertexId, xval, yval, zval)) {
      this->printErr("getChildLocation(): FAILED to get the coordinate values "
                     "of the vertex id "
                     + to_string(vertexId));
      return -1;
    }

    if(xval < centerArr[0]) {
      if(yval < centerArr[1]) {
        if(zval < centerArr[2]) {
          return (parLoc << 3);
        } else {
          return (parLoc << 3) + 1;
        }
      } else {
        if(zval < centerArr[2]) {
          return (parLoc << 3) + 2;
        } else {
          return (parLoc << 3) + 3;
        }
      }
    } else {
      if(yval < centerArr[1]) {
        if(zval < centerArr[2]) {
          return (parLoc << 3) + 4;
        } else {
          return (parLoc << 3) + 5;
        }
      } else {
        if(zval < centerArr[2]) {
          return (parLoc << 3) + 6;
        } else {
          return (parLoc << 3) + 7;
        }
      }
    }

    this->printErr("getChildLocation(): Shouldn't get here!");
    return -1;
  }

  /**
   * Recursive subdividing function to make sure each node does not exceed the
   * capacity limit.
   */
  void subdivide(OctreeNode *node) {
    if(node == nullptr)
      return;

    if((int)node->vertexIds->size() > capacity) {
      uint32_t childCode = 0;
      float ncenter[3] = {0.0}, nsize[3] = {0.0};

      computeCenterSize(node->locCode, ncenter, nsize);

      for(int v : *(node->vertexIds)) {
        childCode = getChildLocation(node->locCode, v, ncenter);
        OctreeNode *childNode = lookupNode(childCode);

        if(childNode == nullptr) {
          OctreeNode newnode(childCode);
          newnode.vertexIds->push_back(v);
          allNodes[childCode] = newnode;
          node->childExists |= (1 << (childCode & 7));
        } else {
          childNode->vertexIds->push_back(v);
        }
      }

      delete node->vertexIds;
      node->vertexIds = nullptr;

      for(uint32_t i = 0; i < 8; i++) {
        if(node->childExists & (1 << i)) {
          subdivide(lookupNode((node->locCode << 3) | i));
        }
      }
    }
  }
};