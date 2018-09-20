#include<PersistenceDiagram.h>

using namespace std;
using namespace ttk;

using namespace ftm;

PersistenceDiagram::PersistenceDiagram():
  ComputeSaddleConnectors{},

  triangulation_{},
  inputScalars_{},
  CTDiagram_{}
{}

PersistenceDiagram::~PersistenceDiagram(){
}

CriticalIndex PersistenceDiagram::getNodeType(FTMTree_MT* tree,
                                         TreeType    treeType,
                                         const SimplexId   vertexId) const
{
   const Node* node = tree->vertex2Node(vertexId);
   int         upDegree{};
   int         downDegree{};
   if (treeType == TreeType::Join or treeType == TreeType::Contour) {
      upDegree   = node->getNumberOfUpSuperArcs();
      downDegree = node->getNumberOfDownSuperArcs();
   } else {
      upDegree   = node->getNumberOfDownSuperArcs();
      downDegree = node->getNumberOfUpSuperArcs();
   }
   int degree = upDegree + downDegree;

   // saddle point
   if (degree > 1) {
      if (upDegree > 1)
         return CriticalIndex::Saddle2;
      else return CriticalIndex::Saddle1;
   }
   // local extremum
   else {
      if (upDegree)
         return CriticalIndex::Local_minimum;
      else return CriticalIndex::Local_maximum;
  }
}
