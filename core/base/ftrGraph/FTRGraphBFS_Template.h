#ifndef FTRGRAPHBFS_TEMPLATE_H
#define FTRGRAPHBFS_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::bfsPropagation(const idVertex saddle, const idCell seed,
                                                Propagation* const newLocalProp, const int bfsId)
      {
         if (bfsCells_[seed] != bfsId) {
            bfsCells_[seed] = bfsId;
            idEdge nbEdges = mesh_->getTriangleEdgeNumber(seed);
            for(idEdge en = 0; en < nbEdges; ++en) {
               idEdge edge;
               mesh_->getTriangleEdge(seed, en, edge);
               idVertex v0, v1;
               mesh_->getEdgeVertex(edge, 0, v0);
               mesh_->getEdgeVertex(edge, 1, v1);
               bool comp0 = newLocalProp->compare(saddle, v0);
               bool comp1 = newLocalProp->compare(saddle, v1);
               // crossing edge
               if (comp0 != comp1) {
                  // Add the lower one in propagation
                  if(comp0) {
                     if (bfsVerts_[v0] != bfsId) {
                        bfsVerts_[v0] = bfsId;
                        newLocalProp->addNewVertex(v0);
                     }
                  } else {
                     if (bfsVerts_[v1] != bfsId) {
                        bfsVerts_[v1] = bfsId;
                        newLocalProp->addNewVertex(v1);
                     }
                  }
                  // Recursively continue BFS
                  idCell nbTri = mesh_->getEdgeTriangleNumber(edge);
                  for (idCell tn = 0; tn < nbTri; ++tn) {
                     idCell neighTriangle;
                     mesh_->getEdgeTriangle(edge, tn, neighTriangle);
                     if (neighTriangle == seed) {
                        continue;
                     }
                     bfsPropagation(saddle, neighTriangle, newLocalProp, bfsId);
                  }

               }
            }
         }
      }

   }
}


#endif /* end of include guard: FTRGRAPHBFS_TEMPLATE_H */
