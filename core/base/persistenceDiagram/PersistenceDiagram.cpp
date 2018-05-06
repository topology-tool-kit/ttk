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

NodeType PersistenceDiagram::getNodeType(FTMTree_MT* tree,
                                         TreeType    treeType,
                                         const int   vertexId) const
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
         return NodeType::Saddle2;
      else return NodeType::Saddle1;
   }
   // local extremum
   else {
      if (upDegree)
         return NodeType::Local_minimum;
      else return NodeType::Local_maximum;
  }
}
