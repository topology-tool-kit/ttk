/// \ingroup baseCode
/// \class ttk::MergeTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the
/// sublevel set tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.

#ifndef _MERGETREE_H
#define _MERGETREE_H

#include <queue>
#include <vector>

#ifdef __APPLE__
# include <algorithm>
# include <numeric>
#else
# ifdef _WIN32
#  include <algorithm>
#  include <numeric>
# else
#  include <parallel/algorithm>
# endif
#endif

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "DataTypes.h"
#include "ExtendedUF.h"
#include "Node.h"
#include "Structures.h"
#include "SuperArc.h"

namespace ttk
{
   class MergeTree : virtual public Debug
   {
      friend class ContourForests;
      friend class ContourForestsTree;

     protected:
      // global
      Params *const  params_;
      Triangulation *mesh_;
      Scalars *const scalars_;

      // local
      TreeData treeData_;

     public:

      // CONSTRUCT
      // -----------
      // {

      // Tree with global data and partition number
      MergeTree(Params *const params, Triangulation *mesh, Scalars *const scalars, TreeType type,
                idPartition part = nullPartition);

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
         vector<idVertex> &sosVect = scalars_->sosOffsets;
         if (!sosVect.size()) {
            sosVect.resize(scalars_->size);
            iota(sosVect.begin(), sosVect.end(), 0);
         }
      }

      /// \brief init the type of the current tree froms params
      void initTreeType(void)
      {
         treeData_.treeType = params_->treeType;
      }

      /// \brief if sortedVertices_ is null, define and fill it
      /// Also fill the mirror vector
      template <typename scalarType>
      void sortInput(void);

      /// \brief clear local data for new computation
      void flush(void)
      {
         treeData_.superArcs.clear();
         treeData_.nodes.clear();
         treeData_.leaves.clear();
         treeData_.arcsCrossingAbove.clear();
         treeData_.arcsCrossingBelow.clear();
         treeData_.vert2tree.clear();
         treeData_.vert2tree.resize(scalars_->size, nullCorresp);
      }

      //}
      // -------------
      // ACCESSOR
      // ------------
      // {
      //{

      // global
      // called for the tree used by the wrapper (only).
      // On this implementation, the warpper communicate with ContourForest
      // A child class of this one.

      inline int setDebugLevel(const int  &local_debugLevel)
      {
         Debug::setDebugLevel(local_debugLevel);
         params_->debugLevel = local_debugLevel;
         return 0;
      }

      inline void setTreeType(const int &local_treeType)
      {
         params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setSimplificationMethod(const int &local_simplifyMethod)
      {
         params_->simplifyMethod = static_cast<SimplifMethod>(local_simplifyMethod);
      }

      inline void setSimplificationThreshold(const double &local_simplificationThreshold)
      {
         params_->simplifyThreshold = local_simplificationThreshold;
      }

      inline void setScalars(void * local_scalars)
      {
        scalars_->values = local_scalars;
      }

      inline void setupTriangulation(Triangulation* m, const bool preproc = true)
      {
        mesh_ = m;
        if (mesh_ && preproc) {
           mesh_->preprocessEdges();
           mesh_->preprocessVertexNeighbors();
        }
      }

      // }
      // partition
      // .....................{

      inline const idPartition getPartition(void) const
      {
         return treeData_.partition;
      }

      // }
      // scalar
      // .....................{

      template <typename scalarType>
      inline const scalarType &getValue(const idVertex &idNode) const
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

      inline void setVertexSoSoffsets(const vector<idVertex>& offsets)
      {
        scalars_->sosOffsets = offsets;
      }

      // }
      // arcs
      // .....................{

      inline const idSuperArc getNumberOfSuperArcs(void) const
      {
         return treeData_.superArcs.size();
      }

      inline const idSuperArc getNumberOfVisibleArcs(void) const
      {
         // Costly ! for dedbug only
         idSuperArc visibleArc = 0;
         for (const SuperArc &arc : treeData_.superArcs) {
            if (arc.isVisible())
               ++visibleArc;
         }
         return visibleArc;
      }

      inline const vector<SuperArc> &getSuperArc(void) const
      {
         // break encapsulation...
         return treeData_.superArcs;
      }

      inline SuperArc *getSuperArc(const idSuperArc &i)
      {
#ifndef withKamikaze
         if ((size_t)i >= treeData_.superArcs.size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << treeData_.superArcs.size() << endl;
            return nullptr;
         }
#endif
         return &(treeData_.superArcs[i]);
      }

      inline idVertex getNumberOfVisibleRegularNode(const idSuperArc &sa)
      {
         // Costly ! for dedbug only
         idVertex   res   = 0;
         SuperArc * a     = getSuperArc(sa);
         const auto nbReg = a->getNumberOfRegularNodes();
         for (idVertex v = 0; v < nbReg; v++) {
            if (!a->isMasqued(v))
               ++res;
         }

         return res;
      }

      inline void addCrossingAbove(const idSuperArc &sa)
      {
         treeData_.arcsCrossingAbove.emplace_back(sa);
      }

      // arcsCrossingBelow is not used.

      // }
      // nodes
      // .....................{

      inline const idNode getNumberOfNodes(void) const
      {
         return treeData_.nodes.size();
      }

      inline const vector<Node>& getNodes(void) const
      {
         // break encapsulation...
         return treeData_.nodes;
      }

      inline Node *getNode(const idNode &nodeId)
      {
         return &(treeData_.nodes[nodeId]);
      }

      // }
      // leaves / root
      // .....................{

      inline const idVertex getNumberOfLeaves(void) const
      {
         return treeData_.leaves.size();
      }

      inline const vector<idNode>& getLeaves(void) const
      {
          // break encapsulation...
         return treeData_.leaves;
      }

      inline const idNode &getLeave(const idNode &id) const
      {
#ifndef withKamikaze
         if ((id < 0) || (size_t)id > (treeData_.leaves.size())) {
            stringstream msg;
            msg << "[MergTree] getLeaves out of bounds : " << id << endl;
            err(msg.str(), fatalMsg);
            return treeData_.leaves[0];
         }
#endif
         return treeData_.leaves[id];
      }

      inline const vector<idNode>& getRoots(void) const
      {
          // break encapsulation...
         return treeData_.roots;
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
      //{

      // test vertex correpondance
      // ...........................{

      inline const bool isCorrespondingArc(const idVertex &val) const
      {
         return !isCorrespondingNull(val) && treeData_.vert2tree[val] >= 0;
      }

      inline const bool isCorrespondingNode(const idVertex &val) const
      {
         return treeData_.vert2tree[val] < 0;
      }

      inline const bool isCorrespondingNull(const idVertex &val) const
      {
         return treeData_.vert2tree[val] == nullCorresp;
      }

      //}
      // Get vertex info
      // ...........................{

      inline const idNode getCorrespondingNodeId(const idVertex &val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingNode(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingNode, ";
            debug << "Vertex :" << val << " is not a node :";
            debug <<  treeData_.vert2tree[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return corr2idNode(val);
      }

      inline const idSuperArc getCorrespondingSuperArcId(const idVertex &val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingArc(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingSuperArcId, ";
            debug << "Vertex :" << val << " is not on an arc :";
            debug <<  treeData_.vert2tree[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return treeData_.vert2tree[val];
      }

      // }
      // Get vertex correponding object
      // ................................{

      inline SuperArc *vertex2SuperArc(const idVertex &vert)
      {
         return &(treeData_.superArcs[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const idVertex &vert)
      {
         return &(treeData_.nodes[getCorrespondingNodeId(vert)]);
      }

      // }
      // Update vertex info
      // ................................{

      inline void updateCorrespondingArc(const idVertex &arc, const idSuperArc &val)
      {
         treeData_.vert2tree[arc] = val;
      }

      inline void updateCorrespondingNode(const idVertex &vert, const idNode &val)
      {
         treeData_.vert2tree[vert] = idNode2corr(val);
      }

      inline const idCorresp idNode2corr(const idNode &id) const
      {
         // transform idNode to special value for the array : -idNode -1
         return -(idCorresp)(id + 1);
      }

      inline const idNode corr2idNode(const idCorresp &corr) const
      {
          return -(idNode)(treeData_.vert2tree[corr]+1);
      }


      // }

      // }
      // -------------------
      // Process
      // -------------------
      //{

      // build
      // ..........................{

      // Merge tree processing of a vertex during build
      void processVertex(const idVertex &vertex, vector<ExtendedUnionFind *> &vect_baseUF,
                         const bool overlapB, const bool overlapA, DebugTimer &begin);

      /// \brief Compute the merge tree using Carr's algorithm
      int build(vector<ExtendedUnionFind *> &vect_baseUF, const vector<idVertex> &overlapBefore,
                const vector<idVertex> &overlapAfter, idVertex start, idVertex end,
                const idVertex &posSeed0, const idVertex &posSeed1);

      // }
      // Simplify
      // ...........................{

      // BFS simplification for local CT
      template <typename scalarType>
      idEdge localSimplify(const idVertex &podSeed0, const idVertex &podSeed1);

      // BFS simpliciation for global CT
      template <typename scalarType>
      idEdge globalSimplify(const idVertex posSeed0, const idVertex posSeed1);

      // Having sorted pairs, simplify the current tree
      // in accordance with threashol, between the two seeds.
      template <typename scalarType>
      idEdge simplifyTree(const idVertex &posSeed0, const idVertex &posSeed1,
                          const vector<tuple<idVertex, idVertex, scalarType, bool>> &sortedPairs);

      // add this arc in the subtree which is in the parentNode
      void markThisArc(vector<ExtendedUnionFind *> &ufArray, const idNode &curNodeId,
                       const idSuperArc &mergingArcId, const idNode &parentNodeId);
      // }
      // PersistencePairs
      // ...........................{

      template <typename scalarType>
      int computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType>> &pairs);

      template <typename scalarType>
      int computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs);

      // Construct abstract JT / ST on a CT and fill pairs in accordance.
      // used for global simplification
      template <typename scalarType>
      void recoverMTPairs(const vector<idNode> &sortedNodes,
                          vector<tuple<idVertex, idVertex, scalarType, bool>> &pairsJT,
                          vector<tuple<idVertex, idVertex, scalarType, bool>> &pairsST);

      // }

      // }
      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // {

      // SuperArcs
      // .......................{

      idSuperArc openSuperArc(const idNode &downNodeId, const bool overlapB, const bool overlapA);

      idSuperArc makeSuperArc(const idNode &downNodeId, const idNode &upNodeId, const bool overlapB,
                              const bool overlapA, pair<idVertex, bool> *vertexList = nullptr,
                              idVertex vertexSize = -1);

      void closeSuperArc(const idSuperArc &superArcId, const idNode &upNodeId, const bool overlapB,
                         const bool overlapA);

      void hideArc(const idSuperArc &sa);

      void mergeArc(const idSuperArc &sa, const idSuperArc &recept,
                    const bool changeConnectivity = true);

      const idVertex insertNodeAboveSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed);

      const idVertex getVertBelowSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed,
                                      const vector<idCorresp> &vert2treeOther);

      // is there an external arc linkind node with treeNode in tree
      const bool alreadyExtLinked(const idNode &node, const idPartition &tree,
                                  const idNode &treeNode);

      idSuperArc getNumberOfExternalDownArcs(const idNode &node);

      // TODO Remove that

      void removeHiddenDownArcs(const idNode &n);

      void removeInternalDownArcs(const idNode& node);

      idSuperArc getNumberOfVisibleArcs(const idNode &n);

      idSuperArc getNumberOfUnmergedDownArcs(const idNode &n);


      // }
      // Nodes
      // ...........................{

      idNode makeNode(const idVertex &vertexId, const idVertex &linked = nullVertex);

      idNode makeNode(const Node *const n, const idVertex &linked = nullVertex);

      const idSuperArc insertNode(Node *node, const bool segment);

      const idSuperArc reverseInsertNode(Node *node, const bool segment);

      inline Node *getDownNode(const SuperArc *a);

      inline Node *getUpNode(const SuperArc *a);

      idNode getParent(const idNode &n);

      void delNode(const idNode &node, const pair<idVertex, bool> *mv = nullptr,
                   const idVertex &nbm = 0);

      void hideNode(const idNode &node);

      // For persistance pair on CT
      // these function allow to make a JT / ST od the CT
      const vector<idNode>&& getNodeNeighbors(const idNode &node);

      const vector<idNode>&& getNodeUpNeighbors(const idNode &n);

      const vector<idNode>&& getNodeDownNeighbors(const idNode &n);

      // Remove part not in partition

      void hideAndClearArcsAbove(const idNode &baseNode);

      void hideAndClearArcsBelow(const idNode &baseNode, const idVertex &seed);

      idSuperArc hideAndClearLeadingTo(const idNode &baseNode, const idVertex &v);

      // }
      // Update informations
      // ...........................{

      void updateSegmentation();

      void parallelUpdateSegmentation(const bool ct = false);

      // will disapear
      void parallelInitNodeValence(const int nbThreadValence);

      // }

      // }
      // ---------------------------
      // Operators : print & clone
      // ---------------------------
      // {

      // Print
      void printTree2(void);

      string printArc(const idSuperArc &a)
      {
         const SuperArc *sa = getSuperArc(a);
         stringstream    res;
         res << a << ": ";
         if (sa->getDownCT() == treeData_.partition)
            res << getNode(sa->getDownNodeId())->getVertexId() << " -- ";
         else
            res << "(extern) -- ";

         if (sa->getUpCT() == treeData_.partition)
            res << getNode(sa->getUpNodeId())->getVertexId();
         else
            res << "(extern)";

         res << " \t\t(vis:" << sa->isVisible() << ")";
         return res.str();
      }

      string printNode(const idNode &n)
      {
         const Node *node = getNode(n);
         stringstream res;
         res << n << " : (";
         res << node->getVertexId() << ") / ";

         for (idSuperArc i = 0; i < node->getNumberOfUpSuperArcs(); ++i) {
            if (getSuperArc(node->getUpSuperArcId(i))->isVisible()) {
               res << "+";
            } else {
               res << "-";
            }
            res << node->getUpSuperArcId(i) << " ";
         }

         res << " \\ ";

         for (idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); ++i) {
            if (getSuperArc(node->getDownSuperArcId(i))->isVisible()) {
               res << "+";
            } else {
               res << "-";
            }
            res << node->getDownSuperArcId(i) << " ";
         }

         res << "\t\t(vis:" << node->isVisible() << " )";
         return res.str();
      }

      // Clone
      MergeTree *clone() const;

      void clone(const MergeTree *mt);

      void doSwap(MergeTree *mt);

      //}

     protected:
      // ------------------
      // Comparisons
      // -----------------
      // {

      // Strict

      inline const bool isLower(const idVertex &a, const idVertex &b) const
      {
         return scalars_->mirrorVertices[a] < scalars_->mirrorVertices[b];
      }

      inline const bool isHigher(const idVertex &a, const idVertex &b) const
      {
         return scalars_->mirrorVertices[a] > scalars_->mirrorVertices[b];
      }

      // Large

      inline const bool isEqLower(const idVertex &a, const idVertex &b) const
      {
         return scalars_->mirrorVertices[a] <= scalars_->mirrorVertices[b];
      }

      inline const bool isEqHigher(const idVertex &a, const idVertex &b) const
      {
         return scalars_->mirrorVertices[a] >= scalars_->mirrorVertices[b];
      }

      //}

     private:

      // ------------------
      // Comparisons
      // -----------------
      // {
      // Compare using the scalar array : only for sort step

      template <typename scalarType>
      inline const bool isLower(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 scalars_->sosOffsets[a] < scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline const bool isHigher(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 scalars_->sosOffsets[a] > scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline const bool isEqLower(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 scalars_->sosOffsets[a] <= scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline const bool isEqHigher(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 scalars_->sosOffsets[a] >= scalars_->sosOffsets[b]);
      }

      // }
      // ----------------
      // Simplification
      // ----------------
      // {

      // preserve = do no hide it.
      void hideAndMerge(const idSuperArc &mergingArcId, const idSuperArc &receptacleArcId,
                        const bool preserveDownNode = false);

      // Use BFS from root to find down and up of the receptarc (maintaining segmentation information)
      const tuple<idNode, idNode, idVertex> createReceptArc(
          const idNode &root, const idSuperArc &receptArcId, vector<ExtendedUnionFind *> &arrayUF,
          const vector<pair<idSuperArc, idSuperArc>> &valenceOffsets);

      // during this BFS nodes should have only one arc up/down : find it :
      const idSuperArc newUpArc(const idNode &curNodeId, vector<ExtendedUnionFind *> &ufArray);

      const idSuperArc newDownArc(const idNode &curNodeId, vector<ExtendedUnionFind *> &ufArray);

      // }
      // --------------
      // Tool
      // --------------
      // {
      // create a pair with relative order : child vertex first

      inline pair<idVertex, idVertex> reorderEdgeRel(const pair<idVertex, idVertex> &vert)
      {
         if (treeData_.treeType == TreeType::Split) {
            if (isHigher(vert.first, vert.second)) {
               return vert;
            }

            return make_pair(vert.second, vert.first);
         }  // else

         if (isLower(vert.first, vert.second))
            return vert;

         return make_pair(vert.second, vert.first);
      }

      bool verifyTree(void);

      // Create a pair with the value corresponding to the simplification method

      template <typename scalarType>
      void addPair(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs, const idVertex &orig,
                   const idVertex &term, const bool goUp);

      template <typename scalarType>
      void addPair(vector<tuple<idVertex, idVertex, scalarType>> &pairs, const idVertex &orig,
                   const idVertex &term);

      // }
   };

   ostream &operator<<(ostream &o, Node const &n);
   ostream &operator<<(ostream &o, SuperArc const &a);

#include <MergeTreeTemplate.h>
}

#endif /* end of include guard: MERGETREE_H */
