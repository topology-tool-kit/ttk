/// \ingroup baseCode
/// \class ttk::ContourForests
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef _CONTOURFOREST_H
#define _CONTOURFOREST_H

#include <typeinfo>

#include "ContourForestsTree.h"

namespace ttk
{
   // Classes Interface
   //         ContourForests

   class Interface
   {
     private:
      // same size : number of conn. components.
      idVertex seed_;

      // Overlap
      vector<idVertex> lowerOverlap_;
      vector<idVertex> upperOverlap_;

     public:
      Interface(const idVertex &seed);

      // Getter & Setter
      // {

      inline const idVertex &getSeed(void) const
      {
         return seed_;
      }

      inline vector<idVertex> &getUpper(void)
      {
         return upperOverlap_;
      }

      inline vector<idVertex> &getLower(void)
      {
         return lowerOverlap_;
      }
      inline idVertex getNbUpper(void) const
      {
         return upperOverlap_.size();
      }

      inline idVertex getNbLower(void) const
      {
         return lowerOverlap_.size();
      }

      inline void setSeed(const idVertex &local_seed)
      {
         seed_ = local_seed;
      }

      inline void addUpper(const idVertex &lb)
      {
         upperOverlap_.emplace_back(lb);
      }

      inline void addLower(const idVertex &lb)
      {
         lowerOverlap_.emplace_back(lb);
      }

      inline void upReserve(const idVertex &u)
      {
          upperOverlap_.reserve(u);
      }

      inline void loReserve(const idVertex &l)
      {
        lowerOverlap_.reserve(l);
      }

      inline void appendUpper(const vector<idVertex> &vertices)
      {
         upperOverlap_.insert(upperOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      inline void appendLower(const vector<idVertex> &vertices)
      {
         lowerOverlap_.insert(lowerOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      inline void swapUpper(vector<idVertex> &verts)
      {
         upperOverlap_.swap(verts);
      }

      inline void swapLower(vector<idVertex> &verts)
      {
         lowerOverlap_.swap(verts);
      }

      // }
   };


   struct ParallelParams {
      numThread   nbThreads;
      idInterface nbInterfaces;
      idPartition nbPartitions;
      int         partitionNum;
      bool        lessPartition;
   };

   struct ParallelData {
      vector<Interface>   interfaces;
      vector<ContourForestsTree> trees;
   };

   class ContourForests : public ContourForestsTree
   {
     private:
      // global (one instance -> no pointer)
      ParallelParams parallelParams_;

      // local
      ParallelData parallelData_;


     public:
      ContourForests();

      virtual ~ContourForests();

      // Getters & Setters
      // {

      inline void setThreadNumber(const unsigned short nbThread)
      {
         if (nbThread) {
            parallelParams_.nbThreads = nbThread;
         } else {
            parallelParams_.nbThreads = OsCall::getNumberOfCores();
         }
      }

      inline void setPartitionNum(int p)
      {
          parallelParams_.partitionNum = p;
      }

      inline void setLessPartition(bool l)
      {
          parallelParams_.lessPartition = l;
      }

      // range of partitions, position of seeds , ...

      inline tuple<idVertex, idVertex> getJTRange(const idPartition& i) const
      {
         const idVertex &start = (i == 0)
                                   ? 0
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i - 1].getSeed()];

         const idVertex &end   = (i == parallelParams_.nbInterfaces)
                                   ? scalars_->size
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i].getSeed()];

         return make_tuple(start, end);
      }

      inline tuple<idVertex, idVertex> getSTRange(const idPartition& i) const
      {
         const idVertex &start = (i == parallelParams_.nbInterfaces)
                                   ? scalars_->size -1
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i].getSeed()] -1;

         const idVertex &end   = (i == 0)
                                   ? -1
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i-1].getSeed()] -1;

         return make_tuple(start, end);
      }

      inline tuple<idVertex, idVertex> getSeedsPos(const idPartition &i) const
      {
         const idVertex &seed0 = (i == 0)
                                   ? -1
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i-1].getSeed()];

         const idVertex &seed1 = (i == parallelParams_.nbInterfaces)
                                   ? nullVertex
                                   : scalars_->mirrorVertices[parallelData_.interfaces[i].getSeed()];

         return make_tuple(seed0, seed1);
      }

      inline tuple<vector<idVertex>, vector<idVertex>> getOverlaps(const idPartition &i)
      {
         const vector<idVertex> &lower = (i == 0)
                                   ? vector<idVertex>()
                                   : parallelData_.interfaces[i-1].getLower();

         const vector<idVertex> &upper = (i == parallelParams_.nbInterfaces)
                                   ? vector<idVertex>()
                                   : parallelData_.interfaces[i].getUpper();

         return make_tuple(lower, upper);
      }

      idPartition vertex2partition(const idVertex &v);

      // }

      // Init
      // {
      void initInterfaces(void);

      void initOverlap(void);

      void initNbPartitions(void);

      //}
      // Process
      // {

      template <typename scalarType>
      int build();

      template <typename scalarType>
      int parallelBuild(vector<vector<ExtendedUnionFind *>> &baseUF_JT,
                        vector<vector<ExtendedUnionFind *>> &baseUF_ST);

      void stitch(void);
      void stitchTree(const char tree);

      // replace distributed tree by a global one, will be removed
      void unify();
      void unifyTree(const char treetype);
      // }

      // Print
      // {
      void printDebug(DebugTimer &timer, const string &str);

      void printVectCT();
      // }
   };

#include <ContourForestsTemplate.h>
}

#endif
