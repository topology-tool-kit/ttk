#ifndef FTRGRAPHBFS_TEMPLATE_H
#define FTRGRAPHBFS_TEMPLATE_H

#include "FTRGraph.h"

// BFs
// #define DEBUG_3(msg) std::cout msg
#define DEBUG_3(msg)

namespace ttk
{
   namespace ftr
   {
      template <typename ScalarType>
      void FTRGraph<ScalarType>::bfsPropagation(const idVertex saddle, const idCell seed,
                                                Propagation* const newLocalProp, const idSuperArc arc)
      {
         if (bfsCells_[seed] != saddle) {
            bfsCells_[seed] = saddle;
            idEdge nbEdges = mesh_->getTriangleEdgeNumber(seed);
            for(idEdge en = 0; en < nbEdges; ++en) {
               idEdge edge;
               mesh_->getTriangleEdge(seed, en, edge);
               // continue if already seen
               if(bfsEdges_[edge] == saddle) continue;
               bfsEdges_[edge] = saddle;
               // get bondary vertices
               idVertex v0, v1;
               mesh_->getEdgeVertex(edge, 0, v0);
               mesh_->getEdgeVertex(edge, 1, v1);
               bool comp0 = newLocalProp->compare(saddle, v0);
               bool comp1 = newLocalProp->compare(saddle, v1);
               // crossing edge
               if (comp0 != comp1) {
                  // Add the lower one in propagation
                  if(comp0) {
                     if (bfsVerts_[v0] != saddle) {
                        bfsVerts_[v0] = saddle;
                        newLocalProp->addNewVertex(v0);
                     }
                     graph_.visit(v1, arc);
                     DEBUG_3(<< "process " << printEdge(edge, newLocalProp) << " with " << arc << " : " << newLocalProp->getRpz() << std::endl);
                  } else {
                     if (bfsVerts_[v1] != saddle) {
                        bfsVerts_[v1] = saddle;
                        newLocalProp->addNewVertex(v1);
                     }
                     graph_.visit(v0, arc);
                     DEBUG_3(<< "process " << printEdge(edge, newLocalProp) << " with " << arc << " : " << newLocalProp->getRpz() << std::endl);
                  }
                  // Recursively continue BFS
                  idCell nbTri = mesh_->getEdgeTriangleNumber(edge);
                  for (idCell tn = 0; tn < nbTri; ++tn) {
                     idCell neighTriangle;
                     mesh_->getEdgeTriangle(edge, tn, neighTriangle);
                     if (neighTriangle == seed) {
                        continue;
                     }
                     bfsPropagation(saddle, neighTriangle, newLocalProp, arc);
                  }

               }
            }
         }
      }

   }
}


#endif /* end of include guard: FTRGRAPHBFS_TEMPLATE_H */
