#ifndef FTRGRAPHBFS_TEMPLATE_H
#define FTRGRAPHBFS_TEMPLATE_H

#include "FTRGraph.h"

namespace ttk
{
   namespace ftr
   {

      template <typename ScalarType>
      void FTRGraph<ScalarType>::bfsSeed(const std::size_t idt, const valence idcc,
                                         std::vector<idCell>& triangles, std::vector<valence>& cc,
                                         const Propagation* const localProp)
      {
         const idCell   curTri  = triangles[idt];
         const idVertex curVert = localProp->getCurVertex();
         if (cc[idt] == -1){
            cc[idt] = idcc;

            // for each edge
            idEdge nbEdges = mesh_->getTriangleEdgeNumber(curTri);
            // sould be 3
            for (idEdge en = 0; en < nbEdges; ++en) {
               idEdge edge;
               mesh_->getTriangleEdge(curTri, en, edge);
               idVertex v0, v1;
               mesh_->getEdgeVertex(edge, 0, v0);
               mesh_->getEdgeVertex(edge, 1, v1);
               // if edge crossed by the value
               if(localProp->compare(curVert, v0) != localProp->compare(curVert, v1)){
                  // for traingles attached to this edge
                  idCell nbTri = mesh_->getEdgeTriangleNumber(edge);
                  for (idCell tn = 0; tn < nbTri; ++tn) {
                     idCell neighTriangle;
                     mesh_->getEdgeTriangle(edge, tn, neighTriangle);
                     if (neighTriangle == curTri) {
                        continue;
                     }
                     const auto it = std::find(triangles.cbegin(), triangles.cend(), neighTriangle);
                     if (it != triangles.cend()) {
                        bfsSeed(std::distance(triangles.cbegin(), it), idcc, triangles, cc, localProp);
                     }
                  }
               }
            }
         }
      }

      template <typename ScalarType>
      void FTRGraph<ScalarType>::bfsPropagation(const idVertex saddle, const idCell seed,
                                                Propagation* const  newLocalProp,
                                                std::set<idCell>&   visitedCells,
                                                std::set<idVertex>& addedVertices)
      {
         if (visitedCells.find(seed) == end(visitedCells)) {
            visitedCells.emplace(seed);
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
                  // Add in propagation
                  if(comp0) {
                     if (!graph_.isVisited(v0) && addedVertices.find(v0) == end(addedVertices)) {
                        addedVertices.emplace(v0);
                        newLocalProp->addNewVertex(v0);
                        toVisit_[v0] = newLocalProp->getRpz();
                        // TODO remove the set, only use toVisit_
                     }
                  } else {
                     if (!graph_.isVisited(v1) && addedVertices.find(v1) == end(addedVertices)) {
                        addedVertices.emplace(v1);
                        newLocalProp->addNewVertex(v1);
                        toVisit_[v1] = newLocalProp->getRpz();
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

                     bfsPropagation(saddle, neighTriangle, newLocalProp, visitedCells, addedVertices);
                  }

               }
            }
         }
      }

   }
}


#endif /* end of include guard: FTRGRAPHBFS_TEMPLATE_H */
