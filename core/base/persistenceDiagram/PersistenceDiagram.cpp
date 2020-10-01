#include <PersistenceDiagram.h>

using namespace std;
using namespace ttk;

using namespace ftm;

PersistenceDiagram::PersistenceDiagram() {
  setDebugMsgPrefix("PersistenceDiagram");
}

CriticalType PersistenceDiagram::getNodeType(FTMTree_MT *tree,
                                             TreeType treeType,
                                             const SimplexId vertexId) const {
  const Node *node = tree->vertex2Node(vertexId);
  int upDegree{};
  int downDegree{};
  if(treeType == TreeType::Join or treeType == TreeType::Contour) {
    upDegree = node->getNumberOfUpSuperArcs();
    downDegree = node->getNumberOfDownSuperArcs();
  } else {
    upDegree = node->getNumberOfDownSuperArcs();
    downDegree = node->getNumberOfUpSuperArcs();
  }
  int degree = upDegree + downDegree;

  // saddle point
  if(degree > 1) {
    if(upDegree > 1)
      return CriticalType::Saddle2;
    else
      return CriticalType::Saddle1;
  }
  // local extremum
  else {
    if(upDegree)
      return CriticalType::Local_minimum;
    else
      return CriticalType::Local_maximum;
  }
}

void ttk::PersistenceDiagram::sortPersistenceDiagram(
  std::vector<PersistencePair> &diagram, const SimplexId *const offsets) const {

  auto cmp = [offsets](const PersistencePair &a, const PersistencePair &b) {
    return offsets[a.birth] < offsets[b.birth];
  };

  std::sort(diagram.begin(), diagram.end(), cmp);
}
