/// \ingroup baseCode
//
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK structures for the contour tree
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include<iterator>

#include "DataTypes.h"
#include "Node.h"
#include "Structures.h"
#include "SuperArc.h"

namespace ttk
{
   // Compute parameters (global)
   struct Params {
      int           debugLevel;
      TreeType      treeType;
      SimplifMethod simplifyMethod;
      double        simplifyThreshold;
   };

   // Scalar related containers (global)
   struct Scalars {
      idVertex         size;
      void*            values;
      vector<idVertex> sosOffsets;
      vector<idVertex> sortedVertices, mirrorVertices;

      // Need vertices to be sorted : use mirrorVertices.

      bool isLower(const idVertex& a, const idVertex& b) const
      {
         return mirrorVertices[a] < mirrorVertices[b];
      }
      bool isEqLower(const idVertex& a, const idVertex& b) const
      {
         return mirrorVertices[a] <= mirrorVertices[b];
      }

      bool isHigher(const idVertex& a, const idVertex& b) const
      {
         return mirrorVertices[a] > mirrorVertices[b];
      }
      bool isEqHigher(const idVertex& a, const idVertex& b) const
      {
         return mirrorVertices[a] >= mirrorVertices[b];
      }
   };

   // Tree datas ( 1 per tree )
   struct TreeData {
      TreeType    treeType;
      idPartition partition;

      // components : tree / nodes / extrema
      vector<SuperArc> superArcs;
      vector<Node>     nodes;
      vector<idNode>   leaves, roots;

      // arc crossing an interface (one can be in both)
      vector<idSuperArc> arcsCrossingBelow, arcsCrossingAbove;

      // vertex 2 node / superarc
      vector<idCorresp> vert2tree;
   };

   // info on one vertex and CT arc in wich it is
   struct vertex {
      idVertex   id;
      idSuperArc ctArc;
   };

   using segmentIterator = vector<vertex>::iterator;
   using segmentRevIterator = vector<vertex>::reverse_iterator;

   // If we want to cross a Segment in the sorted order,
   // wich is form the end to the beginning in the case of Split Tree,
   // we can do so by using sbegin and send which use this class
   class sorted_iterator : public segmentIterator
   {
     public:
      sorted_iterator(segmentIterator base) : segmentIterator(base), forward_(true)
      {
      }
      sorted_iterator(segmentRevIterator base) : segmentIterator(base.base()), forward_(false)
      {
      }

      const sorted_iterator& operator++()
      {
         if (forward_) {
            vector<vertex>::iterator::operator++();
         } else {
            vector<vertex>::iterator::operator--();
         }

         return *this;
      }

     private:
      bool forward_;
   };

   // Segmentation data
   struct Region {
       // inverted in case of split tree
      segmentIterator segmentBegin;
      segmentIterator segmentEnd;
      bool forward;

      sorted_iterator sbegin() const{
          if(forward)
              return sorted_iterator(segmentBegin);
          return sorted_iterator(segmentRevIterator(segmentEnd));
      }

      sorted_iterator send() const{
          if(forward)
              return sorted_iterator(segmentEnd);
          return sorted_iterator(segmentRevIterator(segmentBegin));
      }
   };
}

#endif /* end of include guard: STRUCTURES_H */
