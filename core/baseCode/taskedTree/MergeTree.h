/// \ingroup baseCode
/// \class ttk::MergeTree
/// \author Charles Gueuent <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the
/// sublevel set tree of scalar data and more
/// (data segmentation
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///

#ifndef _MERGETREE_H
#define _MERGETREE_H

#include <functional>
#include <map>
#include <queue>
#include <set>
#include <vector>

#ifdef __APPLE__
# include <algorithm>
# include <numeric>
#else
# ifdef _WIN32
#  include <algorithm>
#  include <numeric>
# else
#  ifdef __clang__
#   include <algorithm>
#   include <numeric>
#  else
#   include <parallel/algorithm>
#  endif
# endif
#endif

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "AtomicUF.h"
#include "AtomicVector.h"
#include "DataTypes.h"
#include "ExtendedUF.h"
#include "Node.h"
#include "Structures.h"
#include "SuperArc.h"

namespace ttk
{
   using UF = AtomicUF *;

   /*
    * OpenMP use class field as thread-private, but we want to share them in the build
    * we use ptr to allow the copy'ed version to share the same data
    * (pointer private per thread on the same location shared)
    */
   // Tree datas ( 1 per tree )
   struct TreeData {
      TreeType treeType;

      // components : tree / nodes / extrema
      AtomicVector<SuperArc> *superArcs;
      AtomicVector<Node> *    nodes;
      AtomicVector<idNode> *  roots;
      vector<idNode> *        leaves;

      // vertex 2 node / superarc
      vector<idCorresp> *vert2tree;

      // uf
      vector<UF> *ufs, *propagation;
      // valences
      vector<valence> *valences;
      // opened nodes
      vector<char> *openedNodes;

#ifdef withStatsHeight
      vector<idSuperArc> *arcDepth;
      vector<idSuperArc> *arcPotential;
#endif

#ifdef withStatsTime
      vector<float> *   arcStart;
      vector<float> *   arcEnd;
      vector<idVertex> *arcOrig;
      vector<idNode>   *arcTasks;
#endif

      // current nb of tasks
      idNode activeTasks;

      // Segmentation, stay empty for Contour tree as
      // they are created by Merge Tree
      Segments segments_;
   };

   class MergeTree : virtual public Debug
   {
     protected:
      // global
      Params *const  params_;
      Triangulation *mesh_;
      Scalars *const scalars_;

      // local
      TreeData   treeData_;
      Comparison comp_;

      using sortedVertIt = vector<idVertex>::iterator;

     public:
      // -----------
      // CONSTRUCT
      // -----------
      // {

      // Tree with global data and partition number
      MergeTree(Params *const params, Triangulation *mesh, Scalars *const scalars, TreeType type);

      virtual ~MergeTree();

      //}
      // --------------------
      // Init
      // --------------------
      // {

      void initNbScalars(void)
      {
         scalars_->size = mesh_->getNumberOfVertices();
      }

      /// \brief init Simulation of Simplicity datastructure if not set
      void initSoS(void)
      {
         auto *sosVect = scalars_->sosOffsets.get();
         if (sosVect == nullptr) {
            sosVect = new vector<idVertex>(0);
            scalars_->sosOffsets.reset(sosVect);
         }

         if (!sosVect->size()) {
            sosVect->resize(scalars_->size);
            iota(sosVect->begin(), sosVect->end(), 0);
         }
      }

      /// \brief if sortedVertices_ is null, define and fill it
      /// Also fill the mirror vector
      template <typename scalarType>
      void sortInput(void);

      /// \brief clear local data for new computation
      void makeAlloc(void)
      {
         createAtomicVector<SuperArc>(treeData_.superArcs);

         // Stats alloc

         createAtomicVector<Node>(treeData_.nodes);
         treeData_.nodes->reserve(scalars_->size/2);

         createAtomicVector<idNode>(treeData_.roots);
         treeData_.roots->reserve(10);

         createVector<idNode>(treeData_.leaves);
         treeData_.leaves->reserve(scalars_->size/3);

         // Known size

         createVector<idCorresp>(treeData_.vert2tree);
         treeData_.vert2tree->resize(scalars_->size);

         createVector<UF>(treeData_.ufs);
         treeData_.ufs->resize(scalars_->size);

         createVector<UF>(treeData_.propagation);
         treeData_.propagation->resize(scalars_->size);

         createVector<valence>(treeData_.valences);
         treeData_.valences->resize(scalars_->size);

         createVector<char>(treeData_.openedNodes);
         treeData_.openedNodes->resize(scalars_->size);

         treeData_.segments_.clear();
      }

      void makeInit(void) {
          initVector<idCorresp>(treeData_.vert2tree, nullCorresp);
          initVector<UF>(treeData_.ufs, nullptr);
          initVector<UF>(treeData_.propagation, nullptr);
          initVector<valence>(treeData_.valences, 0);
          initVector<char>(treeData_.openedNodes,0);
      }

      //}
      // -------------
      // ACCESSOR
      // ------------
      // {
      //{

      // Tree info for wrapper

#ifdef withStatsHeight
      inline idSuperArc getArcDepth(const idSuperArc arcId)
      {
          return (*treeData_.arcDepth)[arcId];
      }

      inline idVertex getArcPotential(const idSuperArc arcId)
      {
          return (*treeData_.arcPotential)[arcId];
      }
#endif

#ifdef withStatsTime
      inline float getArcStart(const idSuperArc arcId)
      {
          return (*treeData_.arcStart)[arcId];
      }

      inline float getArcEnd(const idSuperArc arcId)
      {
          return (*treeData_.arcEnd)[arcId];
      }

      inline idVertex getArcOrig(const idSuperArc arcId)
      {
          return (*treeData_.arcOrig)[arcId];
      }

      inline idVertex getArcActiveTasks(const idSuperArc arcId)
      {
          return(*treeData_.arcTasks)[arcId];
      }
#endif

      inline idVertex getArcSize(const idSuperArc arcId)
      {
          return getSuperArc(arcId)->size();
      }

      inline bool isJT(void) const
      {
         return treeData_.treeType == TreeType::Join;
      }

      // global
      // called for the tree used by the wrapper (only).
      // On this implementation, the warpper communicate with ContourForest
      // A child class of this one.

      inline void setTreeType(const int &local_treeType)
      {
         params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setScalars(void *local_scalars)
      {
         scalars_->values = local_scalars;
      }

      inline void setupTriangulation(Triangulation *m, const bool preproc = true)
      {
         mesh_ = m;
         if (mesh_ && preproc) {
            // propage through vertices (build)
            mesh_->preprocessVertexNeighbors();
         }
      }

      // }
      // scalar
      // .....................{

      template <typename scalarType>
      inline const scalarType &getValue(idVertex idNode) const
      {
         return (((scalarType *)scalars_->values))[idNode];
      }

      template <typename scalarType>
      inline void setVertexScalars(scalarType *vals)
      {
         scalars_->values = (void *)vals;
      }

      // }
      // offset
      // .....................{

      inline void setVertexSoSoffsets(const vector<idVertex> &offsets)
      {
         // copy input vector
         scalars_->sosOffsets.reset(new vector<idVertex>(offsets.cbegin(), offsets.cend()));
      }

      // }
      // arcs
      // .....................{

      inline idSuperArc getNumberOfSuperArcs(void) const
      {
         return treeData_.superArcs->size();
      }

      inline SuperArc *getSuperArc(idSuperArc i)
      {
#ifndef withKamikaze
         if ((size_t)i >= treeData_.superArcs->size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << treeData_.superArcs->size() << endl;
            return nullptr;
         }
#endif
         return &((*treeData_.superArcs)[i]);
      }

      inline const SuperArc *getSuperArc(idSuperArc i) const
      {
#ifndef withKamikaze
         if ((size_t)i >= treeData_.superArcs->size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << treeData_.superArcs->size() << endl;
            return nullptr;
         }
#endif
         return &((*treeData_.superArcs)[i]);
      }

      // }
      // nodes
      // .....................{

      inline idNode getNumberOfNodes(void) const
      {
         return treeData_.nodes->size();
      }

      inline Node *getNode(idNode nodeId)
      {
         return &((*treeData_.nodes)[nodeId]);
      }

      inline void setValence(const idVertex v, const idVertex val)
      {
          (*treeData_.valences)[v] = val;
      }

      // }
      // leaves / root
      // .....................{

      inline idNode getNumberOfLeaves(void) const
      {
         return treeData_.leaves->size();
      }

      inline const vector<idNode> &getLeaves(void) const
      {
         // break encapsulation...
         return (*treeData_.leaves);
      }

      inline idNode getLeave(const idNode id) const
      {
#ifndef withKamikaze
         if ((id < 0) || (size_t)id > (treeData_.leaves->size())) {
            stringstream msg;
            msg << "[MergTree] getLeaves out of bounds : " << id << endl;
            err(msg.str(), fatalMsg);
            return (*treeData_.leaves)[0];
         }
#endif
         return (*treeData_.leaves)[id];
      }

      inline const vector<idNode> &getRoots(void) const
      {
         // break encapsulation...
         return (*treeData_.roots);
      }

      // }
      // vert2tree
      // .....................{

      inline void setVert2Tree(decltype(treeData_.vert2tree) const vect2tree)
      {
         treeData_.vert2tree = vect2tree;
      }

      // }
      // }
      // --------------------
      // VERT 2 TREE Special functions
      // --------------------
      // {

      // test vertex correpondance
      // ...........................{

      inline bool isCorrespondingArc(const idVertex val) const
      {
         return !isCorrespondingNull(val) && (*treeData_.vert2tree)[val] >= 0;
      }

      inline bool isCorrespondingNode(const idVertex val) const
      {
         return (*treeData_.vert2tree)[val] < 0;
      }

      inline bool isCorrespondingNull(const idVertex val) const
      {
         return (*treeData_.vert2tree)[val] == nullCorresp;
      }

      //   }
      // Get vertex info
      // ...........................{

      inline idNode getCorrespondingNodeId(const idVertex val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingNode(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingNode, ";
            debug << "Vertex :" << val << " is not a node :";
            debug << (*treeData_.vert2tree)[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return corr2idNode(val);
      }

      inline idSuperArc getCorrespondingSuperArcId(const idVertex val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingArc(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingSuperArcId, ";
            debug << "Vertex :" << val << " is not on an arc :";
            debug << (*treeData_.vert2tree)[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return (*treeData_.vert2tree)[val];
      }

      //   }
      // Get corresponding elemnt
      // ................................{
      inline SuperArc *vertex2SuperArc(const idVertex vert)
      {
         return &((*treeData_.superArcs)[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const idVertex vert)
      {
         return &((*treeData_.nodes)[getCorrespondingNodeId(vert)]);
      }
      //   }
      // Update vertex info
      // ................................{

      inline void updateCorrespondingArc(const idVertex vert, const idSuperArc arc)
      {
         (*treeData_.vert2tree)[vert] = arc;
      }

      inline void updateCorrespondingNode(const idVertex vert, const idNode node)
      {
         (*treeData_.vert2tree)[vert] = idNode2corr(node);
      }

      inline idCorresp idNode2corr(const idNode id) const
      {
         // transform idNode to special value for the array : -idNode -1
         return -(idCorresp)(id + 1);
      }

      inline idNode corr2idNode(const idCorresp &corr) const
      {
         return -(idNode)((*treeData_.vert2tree)[corr] + 1);
      }

      //   }
      // }
      // -------------------
      // Process
      // -------------------
      // {

      /// \brief Compute the merge
      void build(const bool ct);

      // extrema

      virtual int precompute();

      // skeleton

      void leaves();

      void processTask(const idVertex startVert, const idNode h, const idVertex orig);

      tuple<bool, bool> propage(CurrentState &currentState, UF curUF);

      void closeAndMergeOnSaddle(idVertex saddleVert);

      void closeOnBackBone(idVertex saddleVert);

      void closeArcsUF(idNode closeNode, UF uf);

      idVertex trunk();

      void assignChunkTrunk(const vector<idVertex> &pendingVerts, idNode &lastVertInRange,
                            idVertex &acc, const idVertex v);

      // stats

      void stats();

      idNode height(const idNode &node, const idNode h=0);

      void createArcPotential(void);

      void arcPotential(const idNode parentId, const idVertex pot = 0);

      // segmentation

      /// \brief use vert2tree to compute the segmentation of the fresh builded merge tree.
      void buildSegmentation();

      // }
      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // {
      // SuperArcs
      // .......................{

      idSuperArc openSuperArc(idNode downNodeId);

      idSuperArc makeSuperArc(idNode downNodeId, idNode upNodeId);

      void closeSuperArc(idSuperArc superArcId, idNode upNodeId);

      void hideArc(idSuperArc sa);

      void mergeArc(idSuperArc sa, idSuperArc recept, const bool changeConnectivity = true);

      // }
      // Nodes
      // ...........................{

      vector<idNode> sortedNodes(const bool parallel = false);

      idNode makeNode(idVertex vertexId, idVertex linked = nullVertex);

      idNode makeNode(const Node *const n, idVertex linked = nullVertex);

      idSuperArc insertNode(Node *node, const bool segm = true);

      Node *getDownNode(const SuperArc *a);

      Node *getUpNode(const SuperArc *a);

      idNode getParent(const idNode n)
      {
         return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
      }

      void delNode(idNode node);

      void hideNode(idNode node);

      // }
      // Segmentation
      // ...........................{

      tuple<segm_it, segm_it> addSimpleSegment(idVertex v)
      {
         return treeData_.segments_.addLateSimpleSegment(v);
      }

      // }
      // Update informations
      // ...........................{

      // Create the segmentation of all arcs by operating the pending operations
      void finalizeSegmentation(void);

      // }

      // }
      // ---------------------------
      // Operators : clone & print
      // ---------------------------
      // {

      // Clone
      MergeTree *clone() const;

      void clone(const MergeTree *mt);

      void doSwap(MergeTree *mt);

      // Print
      string printArc(idSuperArc a);

      string printNode(idNode n);

      void printTree2(void);

      void printParams(void) const;

      int printTime(DebugTimer &t, const string &s, idVertex nbScalars = -1,
                    const int debugLevel = 2) const;

      //}

     protected:

      // -----
      // Tools
      // -----
      // {

      idNode getVertInRange(const vector<idVertex> &range, const idVertex v,
                              const idNode last = 0) const;

      tuple<idVertex, idVertex> getBoundsFromVerts(const vector<idVertex> &nodes) const;

      idSuperArc upArcFromVert(const idVertex v)
      {
         return getNode(getCorrespondingNodeId(v))->getUpSuperArcId(0);
      }

      inline idVertex getChunkSize(const idVertex nbVerts = -1, const idVertex nbtasks = 100) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->size : nbVerts;
         return max(10000, 1 + (s / (nbtasks * threadNumber_)));
      }

      inline idVertex getChunkCount(const idVertex nbVerts = -1, const idVertex nbTasks = 100) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->size : nbVerts;
         return 1 + (s / getChunkSize(s, nbTasks));
      }

      // }
      // ------------------
      // Comparisons
      // -----------------
      // {
      // Compare using the scalar array : only for sort step

      template <typename scalarType>
      inline bool isLower(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 (*scalars_->sosOffsets)[a] < (*scalars_->sosOffsets)[b]);
      }

      template <typename scalarType>
      inline bool isHigher(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 (*scalars_->sosOffsets)[a] > (*scalars_->sosOffsets)[b]);
      }

      template <typename scalarType>
      inline bool isEqLower(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 (*scalars_->sosOffsets)[a] <= (*scalars_->sosOffsets)[b]);
      }

      template <typename scalarType>
      inline bool isEqHigher(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 (*scalars_->sosOffsets)[a] >= (*scalars_->sosOffsets)[b]);
      }

      template <typename type>
      void createVector(vector<type> *&ptr)
      {
         if(!ptr)
            ptr = new vector<type>;
         ptr->clear();
      }

      template <typename type>
      void createAtomicVector(AtomicVector<type> *&ptr)
      {
         if(!ptr)
            ptr = new AtomicVector<type>;
         ptr->clear();
      }

      template <typename type>
      void initVector(vector<type> *&vect, const type val)
      {
         int s = vect->size();
#pragma omp parallel for num_threads(threadNumber_) schedule(static)
         for (int i = 0; i < s; i++) {
            (*vect)[i] = val;
         }
      }

      // }
   };

   ostream &operator<<(ostream &o, Node const &n);
   ostream &operator<<(ostream &o, SuperArc const &a);

#include <MergeTreeTemplate.h>
}

#endif /* end of include guard: MERGETREE_H */
