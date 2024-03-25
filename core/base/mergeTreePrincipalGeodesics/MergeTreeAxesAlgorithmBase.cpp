#include <MergeTreeAxesAlgorithmBase.h>

void ttk::MergeTreeAxesAlgorithmBase::reverseMatchingVector(
  unsigned int noNodes,
  std::vector<ftm::idNode> &matchingVector,
  std::vector<ftm::idNode> &invMatchingVector) {
  invMatchingVector.clear();
  invMatchingVector.resize(noNodes, std::numeric_limits<ftm::idNode>::max());
  for(unsigned int i = 0; i < matchingVector.size(); ++i)
    if(matchingVector[i] < noNodes)
      invMatchingVector[matchingVector[i]] = i;
}
