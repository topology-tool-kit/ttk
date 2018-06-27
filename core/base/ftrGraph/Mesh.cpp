#include "Mesh.h"

using namespace std;
using namespace ttk;
using namespace ftr;

void Mesh::preSortEdges(VertCompFN lowerThan)
{
   // Can't be parallel on a vector of bool !!
   // #ifdef TTK_ENABLE_OPENMP
   // #pragma omp parallel for schedule(static)
   // #endif
   for(idEdge ei = 0; ei < nbEdges_; ++ei) {
      idVertex v0, v1;
      tri_->getEdgeVertex(ei, 0, v0);
      tri_->getEdgeVertex(ei, 1, v1);

      edgesSortId_[ei] = lowerThan(v0, v1);
   }
}

bool Mesh::compareEdges(const idEdge e0, const idEdge e1, VertCompFN lowerThan) const
{
   idVertex e0v0, e0v1;
   if(edgesSortId_[e0] == 0) {
      getEdgeVertex(e0, 0, e0v0);
      getEdgeVertex(e0, 1, e0v1);
   } else {
      getEdgeVertex(e0, 1, e0v0);
      getEdgeVertex(e0, 0, e0v1);
   }

   idVertex e1v0, e1v1;
   if(edgesSortId_[e1] == 0) {
      getEdgeVertex(e1, 0, e1v0);
      getEdgeVertex(e1, 1, e1v1);
   } else {
      getEdgeVertex(e1, 1, e1v0);
      getEdgeVertex(e1, 0, e1v1);
   }

   if (e0v0 == e1v0) return lowerThan(e0v1, e1v1);
   return lowerThan(e0v0, e1v0);
}

void Mesh::preSortTriangles(VertCompFN lowerThan)
{
   // require edges to be already sorted
   //
   // Not sure we can parallelie on a vector of bitfields
   // #ifdef TTK_ENABLE_OPENMP
   // #pragma omp parallel for schedule(static)
   // #endif
   for(idCell ti = 0; ti < nbTriangles_; ++ti)
   {
      idEdge   e0, e1, e2;
      getTriangleEdge(ti, 0, e0);
      getTriangleEdge(ti, 1, e1);
      getTriangleEdge(ti, 2, e2);

      if (compareEdges(e0, e1, lowerThan)) {
         // 1 2 3
         // 1 3 2
         // 2 3 1
         if (compareEdges(e1, e2, lowerThan)) {
            // 0 : 1 2 3
            trianglesSortId_[ti] = 0;
         } else if (compareEdges(e0, e2, lowerThan)) {
            // 1 : 1 3 2
            trianglesSortId_[ti] = 1;
         } else {
            // 3 : 2 3 1
            trianglesSortId_[ti] = 3;
         }
      } else {
         // 2 1 3
         // 3 1 2
         // 3 2 1
         if (compareEdges(e0, e2, lowerThan)) {
            // 2 : 2 1 3
            trianglesSortId_[ti] = 2;
         } else if (compareEdges(e1, e2, lowerThan)) {
            // 4 : 3 1 2
            trianglesSortId_[ti] = 4;
         } else {
            // 5 : 3 2 1
            trianglesSortId_[ti] = 5;
         }
      }
   }
}

orderedEdge Mesh::getOrderedEdge(const idEdge e, const bool increasingOrder) const
{
   idVertex v0, v1;
   tri_->getEdgeVertex(e, 0, v0);
   tri_->getEdgeVertex(e, 1, v1);
   if (edgesSortId_[e] == increasingOrder) {
      return std::make_tuple(v0, v1);
   } else {
      return std::make_tuple(v1, v0);
   }
}

orderedTriangle Mesh::getOrderedTriangle(const idCell ti, const bool increasingOrder) const
{
   idEdge e0, e1, e2;
   getTriangleEdge(ti, 0, e0);
   getTriangleEdge(ti, 1, e1);
   getTriangleEdge(ti, 2, e2);

   switch (trianglesSortId_[ti]) {
      case 0:
         if(increasingOrder) return std::make_tuple(e0, e1, e2);
         else                return std::make_tuple(e2, e1, e0);

      case 1:
         if(increasingOrder) return std::make_tuple(e0, e2, e1);
         else                return std::make_tuple(e1, e2, e0);

      case 2:
         if(increasingOrder) return std::make_tuple(e1, e0, e2);
         else                return std::make_tuple(e2, e0, e1);

      case 3:
         if(increasingOrder) return std::make_tuple(e2, e0, e1);
         else                return std::make_tuple(e1, e0, e2);

      case 4:
         if(increasingOrder) return std::make_tuple(e1, e2, e0);
         else                return std::make_tuple(e0, e2, e1);

      case 5:
         if(increasingOrder) return std::make_tuple(e2, e1, e0);
         else                return std::make_tuple(e0, e1, e2);
   }

#ifndef TTK_ENABLE_KAMIKAZE
   std::cerr << "[Mesh]: error, unknownk case in triangleSortId_ : "
             << static_cast<unsigned>(trianglesSortId_[ti]) << " at " << ti << std::endl;
#endif

   return std::make_tuple(nullEdge, nullEdge, nullEdge);
}

