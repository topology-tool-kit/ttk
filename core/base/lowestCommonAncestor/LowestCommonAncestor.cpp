#include <LowestCommonAncestor.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <sstream>
#include <stack>
#include <string>

ttk::LowestCommonAncestor::LowestCommonAncestor() {
  this->setDebugMsgPrefix("LowestCommonAncestor");
}

int ttk::LowestCommonAncestor::preprocess() {

  Timer t;

  // Return value
  int retval;
  // Eulerian transverse of the tree
  retval = eulerianTransverse();
  if(retval != 0) {
    return retval;
  }
  // Divide the depth array in blocs and do preprocessing
  retval = computeBlocs();
  if(retval != 0) {
    return retval;
  }
  // Preprocess the range minimum queries on the blocs
  blocMinimumValueRMQ_.setVector(blocMinimumValue_);
  blocMinimumValueRMQ_.setDebugLevel(debugLevel_);
  retval = blocMinimumValueRMQ_.preprocess(true);
  if(retval != 0) {
    return retval;
  }

  this->printMsg("Preprocessed queries.", 1.0, t.getElapsedTime(), 1,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0;
}

int ttk::LowestCommonAncestor::RMQuery(const int &i, const int &j) const {
  // Bloc of i
  int blocI = i / blocSize_;
  // Bloc of j
  int blocJ = j / blocSize_;
  // If i and j are in different blocs
  if(blocI != blocJ) {
    // Positions and values of the 3 range minima to compare
    std::array<int, 3> min_pos, min_value;
    // Position of the min in the bloc containing the ith case
    min_pos[0] = blocI * blocSize_
                 + normalizedBlocTable_[blocToNormalizedBloc_[blocI]]
                                       [i % blocSize_][blocSize_ - 1];
    // Position of the min in the blocs between the bloc of i and j
    min_pos[1] = ((blocJ - blocI) > 1)
                   ? blocMinimumPosition_[blocMinimumValueRMQ_.query(
                     blocI + 1, blocJ - 1)]
                   : INT_MAX;
    // Position of the min in the bloc containing the jth case
    min_pos[2]
      = blocJ * blocSize_
        + normalizedBlocTable_[blocToNormalizedBloc_[blocJ]][0][j % blocSize_];
    // Values of the depth to compare
    min_value[0] = nodeDepth_[min_pos[0]];
    min_value[1] = (min_pos[1] != INT_MAX) ? nodeDepth_[min_pos[1]] : INT_MAX;
    min_value[2] = nodeDepth_[min_pos[2]];
    // Return the minimum of the 3
    return min_pos[min_pos_3(min_value)];
  } else {
    // i and j are in the same bloc
    return blocI * blocSize_
           + normalizedBlocTable_[blocToNormalizedBloc_[blocI]][i % blocSize_]
                                 [j % blocSize_];
  }
}

int ttk::LowestCommonAncestor::computeBlocs() {
  // Size and number of blocs
  int sizeOfArray = static_cast<int>(nodeDepth_.size());
  blocSize_ = static_cast<int>(log2(sizeOfArray) / 2.0);
  int numberOfBlocs = sizeOfArray / blocSize_ + (sizeOfArray % blocSize_ != 0);
  // Define the bloc ranges
  for(int i = 0; i < numberOfBlocs; i++) {
    std::pair<int, int> range;
    range.first = i * blocSize_;
    range.second = (i + 1) * blocSize_;
    blocPartition_.push_back(range);
  }
  if(sizeOfArray % blocSize_ != 0) {
    std::pair<int, int> range;
    range.first = (numberOfBlocs - 1) * blocSize_;
    range.second = sizeOfArray;
    blocPartition_.push_back(range);
  }
  // Find the minimum in each bloc
  blocMinimumValue_.resize(numberOfBlocs);
  blocMinimumPosition_.resize(numberOfBlocs);
  for(int i = 0; i < numberOfBlocs; i++) {
    blocMinimumValue_[i] = nodeDepth_[blocPartition_[i].first];
    blocMinimumPosition_[i] = blocPartition_[i].first;
    for(int j = blocPartition_[i].first + 1; j < blocPartition_[i].second;
        j++) {
      if(nodeDepth_[j] < blocMinimumValue_[i]) {
        blocMinimumValue_[i] = nodeDepth_[j];
        blocMinimumPosition_[i] = j;
      }
    }
  }
  // Allocate the query tables
  int numberOfTables = (1 << (blocSize_ - 1));
  normalizedBlocTable_.resize(numberOfTables);
  for(int i = 0; i < numberOfTables; i++) {
    normalizedBlocTable_[i].resize(blocSize_);
    for(int j = 0; j < blocSize_; j++) {
      normalizedBlocTable_[i][j].resize(blocSize_);
    }
  }
  // Build the query table for each possible normalized bloc
  for(int i = 0; i < numberOfTables; i++) {
    // Building of the ith possible bloc
    std::vector<bool> plusOne(blocSize_ - 1);
    int quotient = i;
    int remain;
    for(int j = 0; j < (blocSize_ - 1); j++) {
      remain = quotient % 2;
      quotient /= 2;
      plusOne[blocSize_ - 2 - j] = static_cast<bool>(remain);
    }
    if(blocSize_ < 1) {
      return -1;
    }
    std::vector<int> normalizedBloc(blocSize_);
    normalizedBloc[0] = 0;
    for(int j = 0; j < (blocSize_ - 1); j++) {
      if(plusOne[j]) {
        normalizedBloc[j + 1] = normalizedBloc[j] + 1;
      } else {
        normalizedBloc[j + 1] = normalizedBloc[j] - 1;
      }
    }
    // Give the normalizedBloc to a RangeMinimumQuery object
    RangeMinimumQuery<int> rmq(normalizedBloc);
    rmq.setDebugLevel(debugLevel_);
    rmq.preprocess(true);
    // Double loop to compute all queries
    for(int j = 0; j < blocSize_; j++) {
      normalizedBlocTable_[i][j][j] = j;
      for(int k = j + 1; k < blocSize_; k++) {
        normalizedBlocTable_[i][j][k] = rmq.query(j, k);
        normalizedBlocTable_[i][k][j] = normalizedBlocTable_[i][j][k];
      }
    }
  }
  // Determine the corresponding normalized bloc for each bloc
  blocToNormalizedBloc_.resize(numberOfBlocs, -1);
  for(int i = 0; i < numberOfBlocs; i++) {
    int level = 0;
    int tableId = 0;
    for(int j = blocPartition_[i].first + 1; j < blocPartition_[i].second;
        j++) {
      if(nodeDepth_[j] > nodeDepth_[j - 1]) {
        tableId += (1 << (blocSize_ - 2 - level));
      }
      level++;
    }
    blocToNormalizedBloc_[i] = tableId;
  }
  return 0;
}

int ttk::LowestCommonAncestor::eulerianTransverse() {
  // Find the root
  int rootId = -1;
  for(unsigned int i = 0; i < getNumberOfNodes(); i++) {
    if(node_[i].getAncestorId() == static_cast<int>(i)) {
      rootId = i;
      break;
    }
  }
  if(rootId == -1) {
    this->printErr("Tree root not found.");
    return -1;
  } else {
    this->printMsg("Rout found: node id = " + std::to_string(rootId),
                   debug::Priority::DETAIL);
  }
  if(!(node_[rootId].getNumberOfSuccessors() > 0)) {
    this->printErr("Tree root found with no successor.");
    return -2;
  }
  // Initialize the vectors
  nodeOrder_.clear();
  nodeOrder_.reserve(2 * getNumberOfNodes() + 1);
  nodeDepth_.clear();
  nodeDepth_.reserve(2 * getNumberOfNodes() + 1);
  nodeFirstAppearence_.clear();
  nodeFirstAppearence_.resize(getNumberOfNodes(), -1);
  // Transverse starting from the root
  std::stack<int> nodeStack;
  int depth = 0;
  int iteration = 0;
  std::vector<bool> isVisited(getNumberOfNodes(), false);
  nodeStack.push(rootId);
  while(!nodeStack.empty()) {
    int nodeId = nodeStack.top();
    nodeStack.pop();
    nodeOrder_.push_back(nodeId);
    nodeDepth_.push_back(depth);
    if(!isVisited[nodeId]) {
      // Add ancestor
      if(nodeId != rootId) {
        nodeStack.push(node_[nodeId].getAncestorId());
      }
      // Add successors
      int numberOfSuccessors = node_[nodeId].getNumberOfSuccessors();
      for(int i = 0; i < numberOfSuccessors; i++) {
        nodeStack.push(node_[nodeId].getSuccessorId(i));
      }
      nodeFirstAppearence_[nodeId] = iteration;
      isVisited[nodeId] = true;
    }
    // Next depth
    if(!nodeStack.empty()) {
      if(nodeStack.top() == node_[nodeId].getAncestorId()) {
        depth--;
      } else {
        depth++;
      }
    }
    iteration++;
  }
  return 0;
}
