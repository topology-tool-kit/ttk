#include <Octree.h>

using namespace std;
using namespace ttk;

Octree::Octree(const AbstractTriangulation *t) {
  initialize(t, 1000);
}

Octree::Octree(const AbstractTriangulation *t, const int k) {
  initialize(t, k);
}

Octree::~Octree() {
}

// Initialize the octree with the given triangulation and bucket threshold.
void Octree::initialize(const AbstractTriangulation *t, const int k) {
  this->setDebugMsgPrefix("PR Octree");
  OctreeNode root(1);
  allNodes_[1] = root;
  capacity_ = k;
  triangulation_ = t;

  // find the minimum and maximum coordinate values
  float mins[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  float maxs[3] = {FLT_MIN, FLT_MIN, FLT_MIN};
  SimplexId vertexNum = t->getNumberOfVertices();
  for(int i = 0; i < vertexNum; i++) {
    float coord[3];
    t->getVertexPoint(i, coord[0], coord[1], coord[2]);
    for(int j = 0; j < 3; j++) {
      if(coord[j] < mins[j])
        mins[j] = coord[j];
      if(coord[j] > maxs[j])
        maxs[j] = coord[j];
    }
  }

  for(int j = 0; j < 3; j++) {
    center_[j] = (mins[j] + maxs[j]) / 2.0;
    size_[j] = (maxs[j] - mins[j]) / 2.0;
  }
}

/**
 * Return true if the octree is empty, otherwise false.
 */
bool Octree::empty() {
  OctreeNode *root = lookupNode(1);
  if(root->vertexIds_.empty() && !root->childExists_) {
    return true;
  }

  return false;
}

/**
 * Get the depth of the node in the octree.
 */
size_t Octree::getNodeTreeDepth(const OctreeNode *node) {
  // assert(node->locCode);
  size_t depth = 0;
  for(uint32_t lc = node->locCode_; lc != 1; lc >>= 3)
    depth++;
  return depth;
}

/**
 * Get the parent node of the current node.
 */
OctreeNode *Octree::getParentNode(OctreeNode *node) {
  assert(node->locCode_);
  const uint32_t locCodeParent = node->locCode_ >> 3;
  return lookupNode(locCodeParent);
}

/**
 * Traverse the octree from a given node.
 */
void Octree::visitAll(const OctreeNode *node) {
  if(node == nullptr)
    return;
  this->printMsg(to_string(node->locCode_));
  for(int i = 0; i < 8; i++) {
    if(node->childExists_ & (1 << i)) {
      const uint32_t locCodeChild = (node->locCode_ << 3) | i;
      const OctreeNode *child = lookupNode(locCodeChild);
      visitAll(child);
    }
  }
}

/**
 * Get the total number of vertices in the tree.
 */
int Octree::verifyTree(SimplexId &vertexNum) {
  int vertexCount = 0;
  unordered_map<uint32_t, OctreeNode>::iterator it;
  for(it = allNodes_.begin(); it != allNodes_.end(); it++) {
    if(it->second.childExists_ && it->second.vertexIds_.size() > 0) {
      this->printErr("[Octree] WRONG! The internal node "
                     + to_string(it->second.locCode_)
                     + " should not contain any vertices!");
      return -1;
    }
    if(it->second.vertexIds_.size() > 0) {
      vertexCount += it->second.vertexIds_.size();
    }
  }
  if(vertexCount != vertexNum) {
    this->printErr("The number of vertices in the tree is "
                   + to_string(vertexCount) + ", which is not equal to "
                   + to_string(vertexNum));
    return -1;
  }

  return 0;
}

/**
 * Insert the vertex into the octree by its id.
 */
int Octree::insertVertex(SimplexId &vertexId) {
  if(vertexId < 0 || vertexId >= triangulation_->getNumberOfVertices()) {
    return -1;
  }

  OctreeNode *current = lookupNode(1);
  if(current == nullptr) {
    return -2;
  }
  uint32_t location = current->locCode_;
  std::array<float, 3> ncenter{}, nsize{};

  while(current != nullptr && current->childExists_ != 0) {
    computeCenterSize(location, ncenter, nsize);
    location = getChildLocation(location, vertexId, ncenter);
    current = lookupNode(location);
  }

  if(current == nullptr) {
    OctreeNode newnode(location);
    newnode.vertexIds_ = vector<SimplexId>{vertexId};
    OctreeNode *parent = lookupNode(location >> 3);
    parent->childExists_ |= (1 << (location & 7));
    allNodes_[location] = newnode;
  } else {
    current->vertexIds_.push_back(vertexId);
    subdivide(current);
  }
  return 0;
}

/**
 * Insert the cell into the octree by its id.
 * Note: This function can only be called after inserting all vertices!
 */
int Octree::insertCell(SimplexId &cellId) {
  if(cellId < 0 || cellId >= triangulation_->getNumberOfCells()) {
    return -1;
  }

  int dim = triangulation_->getCellVertexNumber(cellId);
  std::array<float, 3> ncenter{}, nsize{};
  for(int i = 0; i < dim; i++) {
    SimplexId vertexId{};
    OctreeNode *current = lookupNode(1);

    triangulation_->getCellVertex(cellId, i, vertexId);

    while(current != nullptr && current->childExists_ != 0) {
      auto location = current->locCode_;
      computeCenterSize(location, ncenter, nsize);
      location = getChildLocation(location, vertexId, ncenter);
      current = lookupNode(location);
    }

    if(current == nullptr) {
      this->printErr("[Octree] insertCell(): Cannot find the vertex id ("
                     + to_string(vertexId) + ") in the tree!");
      return -1;
    } else {
      if(current->cellIds_.empty()) {
        current->cellIds_ = vector<SimplexId>{cellId};
      } else if(current->cellIds_.back() != cellId) {
        current->cellIds_.push_back(cellId);
      }
    }
  }

  return 0;
}

/**
 * Get reindexed vertices and cells.
 */
void Octree::reindex(vector<SimplexId> &vertices,
                     vector<SimplexId> &nodes,
                     vector<SimplexId> &cells) {
  int totalCells = triangulation_->getNumberOfCells();
  vector<int> cellMap(totalCells, -1);

  OctreeNode *root = lookupNode(1);
  int leafCount = 0;
  int cellCount = 0;

  // use the depth-first search to complete reindexing
  stack<const OctreeNode *> nodeStack;
  nodeStack.push(root);

  while(!nodeStack.empty()) {
    const OctreeNode *topNode = nodeStack.top();
    nodeStack.pop();

    if(topNode == nullptr) {
      this->printErr("[Octree] reindex(): shouldn't get here!");
      break;
    }
    if(topNode->childExists_) {
      for(int i = 0; i < 8; i++) {
        if(topNode->childExists_ & (1 << i)) {
          const uint32_t locCodeChild = (topNode->locCode_ << 3) | i;
          const OctreeNode *child = lookupNode(locCodeChild);
          nodeStack.push(child);
        }
      }
    } else {
      vertices.insert(
        vertices.end(), topNode->vertexIds_.begin(), topNode->vertexIds_.end());
      vector<SimplexId> tmp(topNode->vertexIds_.size(), leafCount);
      nodes.insert(nodes.end(), tmp.begin(), tmp.end());
      leafCount++;

      if(topNode->cellIds_.size() > 0) {
        for(auto it = topNode->cellIds_.begin(); it != topNode->cellIds_.end();
            it++) {
          if(cellMap[*it] == -1) {
            cellMap[*it] = cellCount++;
          }
        }
      }
    }
  }

  int count = 0;
  cells.resize(totalCells);
  for(int i = 0; i < totalCells; i++) {
    if(cellMap[i] == -1) {
      count++;
    }
    cells.at(cellMap[i]) = i;
  }
  this->printMsg("reindex(): There are " + to_string(count)
                 + " wrong entries!");
}

void Octree::computeCenterSize(uint32_t location,
                               std::array<float, 3> &centerArr,
                               std::array<float, 3> &sizeArr) {
  int leftmost = 0;
  uint32_t tmp = location;
  while(tmp >>= 1) {
    leftmost++;
  }
  if(leftmost % 3 != 0) {
    this->printMsg("Location: " + std::to_string(location)
                   + ", leftmost: " + std::to_string(leftmost));
    this->printErr("computeCenterSize(): the location seems not correct!");
    this->printErr("Please try a larger bucket capacity!");
    return;
  }

  // initialize the center and size arrays
  centerArr = center_;
  sizeArr = size_;

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

uint32_t Octree::getChildLocation(uint32_t parLoc,
                                  ttk::SimplexId vertexId,
                                  const std::array<float, 3> &centerArr) {
  float xval{}, yval{}, zval{};
  if(triangulation_->getVertexPoint(vertexId, xval, yval, zval)) {
    this->printErr("getChildLocation(): FAILED to get the coordinate values "
                   "of the vertex id "
                   + std::to_string(vertexId));
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

void Octree::subdivide(OctreeNode *node) {
  if(node == nullptr)
    return;

  if((int)node->vertexIds_.size() > capacity_) {
    uint32_t childCode = 0;
    std::array<float, 3> ncenter{}, nsize{};

    computeCenterSize(node->locCode_, ncenter, nsize);

    for(int v : node->vertexIds_) {
      childCode = getChildLocation(node->locCode_, v, ncenter);
      OctreeNode *childNode = lookupNode(childCode);

      if(childNode == nullptr) {
        OctreeNode newnode(childCode);
        newnode.vertexIds_.push_back(v);
        allNodes_[childCode] = newnode;
        node->childExists_ |= (1 << (childCode & 7));
      } else {
        childNode->vertexIds_.push_back(v);
      }
    }

    node->vertexIds_.clear();

    for(uint32_t i = 0; i < 8; i++) {
      if(node->childExists_ & (1 << i)) {
        subdivide(lookupNode((node->locCode_ << 3) | i));
      }
    }
  }
}
