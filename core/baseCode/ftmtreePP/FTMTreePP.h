/// \ingroup baseCode
/// \class ttk::FTMTreePP
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2017-08-03
///
///\brief TTK processing package that add persistance pairs features
// to the contour tree package.
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkPersistenceDiagram.cpp %for a usage example.

#ifndef FTMTREE_PP_H
#define FTMTREE_PP_H

#include "FTMTree.h"

namespace ttk
{
namespace ftm
{

   class FTMTreePP : public FTMTree
   {
     public:

      FTMTreePP();
      virtual ~FTMTreePP();

      template <typename scalarType>
      void computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>>& pairs,
                                   const bool jt);

     private:
   };

}
}

template <typename scalarType>
void ftm::FTMTreePP::computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>>& pairs,
                                             const bool jt)
{
   ftm::FTMTree_MT* tree = jt ? getJoinTree() : getSplitTree();
   tree->sortLeaves();

   auto& vectL = tree->getLeaves();

   // Debug, check sort for ST (up to bottom)
   for (auto nid : vectL) {
       cout << tree->getNode(nid)->getVertexId() << endl;
   }
}

#endif  // FTMTREE_PP_H
