#include <Octree.h>

using namespace ttk;

Octree::Octree(const AbstractTriangulation *t) {
  Octree(t, 1000);
}

Octree::Octree(const AbstractTriangulation *t, const int k) {
  this->setDebugMsgPrefix("PR Octree");
  OctreeNode root(1);
  allNodes[1] = root;
  capacity = k;
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
    center[j] = (mins[j] + maxs[j]) / 2.0;
    size[j] = (maxs[j] - mins[j]) / 2.0;
  }
}

Octree::~Octree() {
}

/**
 * Return true if the octree is empty, otherwise false.
 */
bool Octree::empty() {
  OctreeNode *root = lookupNode(1);
  if(root->vertexIds->size() == 0 && !root->childExists) {
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
  for(uint32_t lc = node->locCode; lc != 1; lc >>= 3)
    depth++;
  return depth;
}

/**
 * Get the parent node of the current node.
 */
OctreeNode *Octree::getParentNode(OctreeNode *node) {
  assert(node->locCode);
  const uint32_t locCodeParent = node->locCode >> 3;
  return lookupNode(locCodeParent);
}

/**
 * Traverse the octree from a given node.
 */
void Octree::visitAll(const OctreeNode *node) {
  if(node == nullptr)
    return;
  this->printMsg(to_string(node->locCode));
  for(int i = 0; i < 8; i++) {
    if(node->childExists & (1 << i)) {
      const uint32_t locCodeChild = (node->locCode << 3) | i;
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
  for(it = allNodes.begin(); it != allNodes.end(); it++) {
    if(it->second.childExists && it->second.vertexIds) {
      this->printErr("[Octree] WRONG! The internal node "
                     + to_string(it->second.locCode)
                     + " should not contain any vertices!");
      return -1;
    }
    if(it->second.vertexIds) {
      vertexCount += it->second.vertexIds->size();
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
  uint32_t location = current->locCode;
  float ncenter[3] = {0.0}, nsize[3] = {0.0};

  while(current != nullptr && current->childExists != 0) {
    computeCenterSize(location, ncenter, nsize);
    location = getChildLocation(location, vertexId, ncenter);
    current = lookupNode(location);
  }

  if(current == nullptr) {
    OctreeNode newnode(location);
    newnode.vertexIds = new vector<SimplexId>{vertexId};
    OctreeNode *parent = lookupNode(location >> 3);
    parent->childExists |= (1 << (location & 7));
    allNodes[location] = newnode;
  } else {
    current->vertexIds->push_back(vertexId);
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
  uint32_t location;
  float ncenter[3] = {0.0}, nsize[3] = {0.0};
  for(int i = 0; i < dim; i++) {
    SimplexId vertexId;
    OctreeNode *current = lookupNode(1);

    location = current->locCode;
    triangulation_->getCellVertex(cellId, i, vertexId);

    while(current != nullptr && current->childExists != 0) {
      computeCenterSize(location, ncenter, nsize);
      location = getChildLocation(location, vertexId, ncenter);
      current = lookupNode(location);
    }

    if(current == nullptr) {
      this->printErr("[Octree] insertCell(): Cannot find the vertex id ("
                     + to_string(vertexId) + ") in the tree!");
      return -1;
    } else {
      if(current->cellIds == nullptr) {
        current->cellIds = new vector<SimplexId>{cellId};
      } else if(current->cellIds->back() != cellId) {
        current->cellIds->push_back(cellId);
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
    if(topNode->childExists) {
      for(int i = 0; i < 8; i++) {
        if(topNode->childExists & (1 << i)) {
          const uint32_t locCodeChild = (topNode->locCode << 3) | i;
          const OctreeNode *child = lookupNode(locCodeChild);
          nodeStack.push(child);
        }
      }
    } else {
      vertices.insert(
        vertices.end(), topNode->vertexIds->begin(), topNode->vertexIds->end());
      vector<SimplexId> tmp(topNode->vertexIds->size(), leafCount);
      nodes.insert(nodes.end(), tmp.begin(), tmp.end());
      leafCount++;

      if(topNode->cellIds) {
        for(auto it = topNode->cellIds->begin(); it != topNode->cellIds->end();
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
