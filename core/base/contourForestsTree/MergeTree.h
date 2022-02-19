/// \ingroup base
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

#pragma once

#include <map>
#include <numeric>
#include <queue>
#include <vector>

#include <Geometry.h>
#include <Triangulation.h>

#include "DeprecatedDataTypes.h"
#include "DeprecatedNode.h"
#include "DeprecatedStructures.h"
#include "DeprecatedSuperArc.h"
#include "ExtendedUF.h"

namespace ttk {
  namespace cf {
    class MergeTree : virtual public Debug {
      friend class ContourForests;
      friend class ContourForestsTree;

    protected:
      // global
      Params *const params_;
      Scalars *const scalars_;

      // local
      TreeData treeData_;

    public:
      // CONSTRUCT
      // -----------
      // {

      // Tree with global data and partition number
      MergeTree(Params *const params,
                Scalars *const scalars,
                TreeType type,
                idPartition part = nullPartition);

      virtual ~MergeTree();

      //}
      // --------------------
      // Init
      // --------------------
      // {

      template <typename triangulationType>
      void initNbScalars(const triangulationType &tri) {
        scalars_->size = tri->getNumberOfVertices();
      }

      /// \brief init the type of the current tree froms params
      void initTreeType(void) {
        treeData_.treeType = params_->treeType;
      }

      /// \brief if sortedVertices_ is null, define and fill it
      /// Also fill the mirror std::vector
      template <typename scalarType>
      void sortInput(void);

      /// \brief clear local data for new computation
      void flush(void) {
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

      inline int setDebugLevel(const int &local_debugLevel) {
        Debug::setDebugLevel(local_debugLevel);
        params_->debugLevel = local_debugLevel;
        return 0;
      }

      inline void setTreeType(const int &local_treeType) {
        params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setSimplificationMethod(const int &local_simplifyMethod) {
        params_->simplifyMethod
          = static_cast<SimplifMethod>(local_simplifyMethod);
      }

      inline void setSimplificationThreshold(
        const double &local_simplificationThreshold) {
        params_->simplifyThreshold = local_simplificationThreshold;
      }

      inline void setScalars(void *local_scalars) {
        scalars_->values = local_scalars;
      }

      inline void preconditionTriangulation(AbstractTriangulation *const m,
                                            const bool preproc = true) {
        if(m && preproc) {
          m->preconditionEdges();
          m->preconditionVertexNeighbors();
        }
      }

      // }
      // partition
      // .....................{

      inline idPartition getPartition(void) const {
        return treeData_.partition;
      }

      // }
      // scalar
      // .....................{

      template <typename scalarType>
      inline const scalarType &getValue(const SimplexId &nodeId) const {
        return (((scalarType *)scalars_->values))[nodeId];
      }

      template <typename scalarType>
      inline void setVertexScalars(scalarType *vals) {
        scalars_->values = (void *)vals;
      }

      // }
      // offset
      // .....................{

      /**
       * @pre For this function to behave correctly in the absence of
       * the VTK wrapper, ttk::preconditionOrderArray() needs to be
       * called to fill the @p offsets buffer prior to any
       * computation (the VTK wrapper already includes a mecanism to
       * automatically generate such a preconditioned buffer).
       * @see examples/c++/main.cpp for an example use.
       */
      inline void setVertexSoSoffsets(const SimplexId *const offsets) {
        scalars_->sosOffsets = offsets;
      }

      // }
      // arcs
      // .....................{

      inline idSuperArc getNumberOfSuperArcs(void) const {
        return treeData_.superArcs.size();
      }

      inline idSuperArc getNumberOfVisibleArcs(void) const {
        // Costly ! for dedbug only
        idSuperArc visibleArc = 0;
        for(const SuperArc &arc : treeData_.superArcs) {
          if(arc.isVisible())
            ++visibleArc;
        }
        return visibleArc;
      }

      inline const std::vector<SuperArc> &getSuperArc(void) const {
        // break encapsulation...
        return treeData_.superArcs;
      }

      inline SuperArc *getSuperArc(const idSuperArc &i) {
#ifndef TTK_ENABLE_KAMIKAZE
        if((size_t)i >= treeData_.superArcs.size()) {
          std::cout << "[Merge Tree] get superArc on bad id :" << i;
          std::cout << " / " << treeData_.superArcs.size() << std::endl;
        }
#endif
        return &(treeData_.superArcs[i]);
      }

      inline SimplexId getNumberOfVisibleRegularNode(const idSuperArc &sa) {
        // Costly ! for dedbug only
        SimplexId res = 0;
        SuperArc *a = getSuperArc(sa);
        const auto nbReg = a->getNumberOfRegularNodes();
        for(SimplexId v = 0; v < nbReg; v++) {
          if(!a->isMasqued(v))
            ++res;
        }

        return res;
      }

      inline void addCrossingAbove(const idSuperArc &sa) {
        treeData_.arcsCrossingAbove.emplace_back(sa);
      }

      // arcsCrossingBelow is not used.

      // }
      // nodes
      // .....................{

      inline idNode getNumberOfNodes(void) const {
        return treeData_.nodes.size();
      }

      inline const std::vector<Node> &getNodes(void) const {
        // break encapsulation...
        return treeData_.nodes;
      }

      inline Node *getNode(const idNode &nodeId) {
        return &(treeData_.nodes[nodeId]);
      }

      // }
      // leaves / root
      // .....................{

      inline SimplexId getNumberOfLeaves(void) const {
        return treeData_.leaves.size();
      }

      inline const std::vector<idNode> &getLeaves(void) const {
        // break encapsulation...
        return treeData_.leaves;
      }

      inline const idNode &getLeave(const idNode &id) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(id > treeData_.leaves.size()) {
          this->printErr("getLeaves out of bounds: " + std::to_string(id));
          return treeData_.leaves[0];
        }
#endif
        return treeData_.leaves[id];
      }

      inline const std::vector<idNode> &getRoots(void) const {
        // break encapsulation...
        return treeData_.roots;
      }

      // }
      // vert2tree
      // .....................{

      inline void setVert2Tree(decltype(treeData_.vert2tree) const &vect2tree) {
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

      inline bool isCorrespondingArc(const SimplexId &val) const {
        return !isCorrespondingNull(val) && treeData_.vert2tree[val] >= 0;
      }

      inline bool isCorrespondingNode(const SimplexId &val) const {
        return treeData_.vert2tree[val] < 0;
      }

      inline bool isCorrespondingNull(const SimplexId &val) const {
        return treeData_.vert2tree[val] == nullCorresp;
      }

      //}
      // Get vertex info
      // ...........................{

      inline idNode getCorrespondingNodeId(const SimplexId &val) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!isCorrespondingNode(val)) {
          this->printErr("getCorrespondingNode, Vertex: " + std::to_string(val)
                         + " is not a node: "
                         + std::to_string(treeData_.vert2tree[val]));
        }
#endif
        return corr2idNode(val);
      }

      inline idSuperArc getCorrespondingSuperArcId(const SimplexId &val) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!isCorrespondingArc(val)) {
          this->printErr(
            "getCorrespondingSuperArcId, Vertex: " + std::to_string(val)
            + " is not on an arc: " + std::to_string(treeData_.vert2tree[val]));
        }
#endif
        return treeData_.vert2tree[val];
      }

      // }
      // Get vertex correponding object
      // ................................{

      inline SuperArc *vertex2SuperArc(const SimplexId &vert) {
        return &(treeData_.superArcs[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const SimplexId &vert) {
        return &(treeData_.nodes[getCorrespondingNodeId(vert)]);
      }

      // }
      // Update vertex info
      // ................................{

      inline void updateCorrespondingArc(const SimplexId &arc,
                                         const idSuperArc &val) {
        treeData_.vert2tree[arc] = val;
      }

      inline void updateCorrespondingNode(const SimplexId &vert,
                                          const idNode &val) {
        treeData_.vert2tree[vert] = idNode2corr(val);
      }

      inline idCorresp idNode2corr(const idNode &id) const {
        // transform idNode to special value for the array : -idNode -1
        return -static_cast<idCorresp>(id + 1);
      }

      inline idNode corr2idNode(const idCorresp &corr) const {
        return static_cast<idNode>(-(treeData_.vert2tree[corr] + 1));
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
      template <typename triangulationType>
      void processVertex(const SimplexId &vertex,
                         std::vector<ExtendedUnionFind *> &vect_baseUF,
                         const bool overlapB,
                         const bool overlapA,
                         const triangulationType &mesh,
                         DebugTimer &begin);

      /// \brief Compute the merge tree using Carr's algorithm
      template <typename triangulationType>
      int build(std::vector<ExtendedUnionFind *> &vect_baseUF,
                const std::vector<SimplexId> &overlapBefore,
                const std::vector<SimplexId> &overlapAfter,
                SimplexId start,
                SimplexId end,
                const SimplexId &posSeed0,
                const SimplexId &posSeed1,
                const triangulationType &mesh);

      // }
      // Simplify
      // ...........................{

      // BFS simplification for local CT
      template <typename scalarType>
      SimplexId localSimplify(const SimplexId &podSeed0,
                              const SimplexId &podSeed1);

      // BFS simpliciation for global CT
      template <typename scalarType, typename triangulationType>
      SimplexId globalSimplify(const SimplexId posSeed0,
                               const SimplexId posSeed1,
                               const triangulationType &mesh);

      // Having sorted std::pairs, simplify the current tree
      // in accordance with threashol, between the two seeds.
      template <typename scalarType>
      SimplexId simplifyTree(
        const SimplexId &posSeed0,
        const SimplexId &posSeed1,
        const std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>>
          &sortedPairs);

      // add this arc in the subtree which is in the parentNode
      void markThisArc(std::vector<ExtendedUnionFind *> &ufArray,
                       const idNode &curNodeId,
                       const idSuperArc &mergingArcId,
                       const idNode &parentNodeId);
      // }
      // PersistencePairs
      // ...........................{

      template <typename scalarType, typename triangulationType>
      int computePersistencePairs(
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
        const triangulationType &mesh);

      template <typename scalarType, typename triangulationType>
      int computePersistencePairs(
        std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairs,
        const triangulationType &mesh);

      // Construct abstract JT / ST on a CT and fill std::pairs in accordance.
      // used for global simplification
      template <typename scalarType, typename triangulationType>
      void recoverMTPairs(
        const std::vector<idNode> &sortedNodes,
        std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>>
          &pairsJT,
        std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>>
          &pairsST,
        const triangulationType &mesh);

      // }

      // }
      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // {

      // SuperArcs
      // .......................{

      idSuperArc openSuperArc(const idNode &downNodeId,
                              const bool overlapB,
                              const bool overlapA);

      idSuperArc makeSuperArc(const idNode &downNodeId,
                              const idNode &upNodeId,
                              const bool overlapB,
                              const bool overlapA,
                              std::pair<SimplexId, bool> *vertexList = nullptr,
                              SimplexId vertexSize = -1);

      void closeSuperArc(const idSuperArc &superArcId,
                         const idNode &upNodeId,
                         const bool overlapB,
                         const bool overlapA);

      void hideArc(const idSuperArc &sa);

      void mergeArc(const idSuperArc &sa,
                    const idSuperArc &recept,
                    const bool changeConnectivity = true);

      SimplexId insertNodeAboveSeed(const idSuperArc &arc,
                                    const std::pair<SimplexId, bool> &seed);

      SimplexId getVertBelowSeed(const idSuperArc &arc,
                                 const std::pair<SimplexId, bool> &seed,
                                 const std::vector<idCorresp> &vert2treeOther);

      // is there an external arc linkind node with treeNode in tree
      bool alreadyExtLinked(const idNode &node,
                            const idPartition &tree,
                            const idNode &treeNode);

      idSuperArc getNumberOfExternalDownArcs(const idNode &node);

      // TODO Remove that

      void removeHiddenDownArcs(const idNode &n);

      void removeInternalDownArcs(const idNode &node);

      idSuperArc getNumberOfVisibleArcs(const idNode &n);

      idSuperArc getNumberOfUnmergedDownArcs(const idNode &n);

      // }
      // Nodes
      // ...........................{

      idNode makeNode(const SimplexId &vertexId,
                      const SimplexId &linked = nullVertex);

      idNode makeNode(const Node *const n,
                      const SimplexId &linked = nullVertex);

      idSuperArc insertNode(Node *node, const bool segment);

      idSuperArc reverseInsertNode(Node *node, const bool segment);

      inline Node *getDownNode(const SuperArc *a);

      inline Node *getUpNode(const SuperArc *a);

      idNode getParent(const idNode &n);

      void delNode(const idNode &node,
                   const std::pair<SimplexId, bool> *mv = nullptr,
                   const SimplexId &nbm = 0);

      void hideNode(const idNode &node);

      // For persistance std::pair on CT
      // these function allow to make a JT / ST od the CT
      std::vector<idNode> getNodeNeighbors(const idNode &node);

      std::vector<idNode> getNodeUpNeighbors(const idNode &n);

      std::vector<idNode> getNodeDownNeighbors(const idNode &n);

      // Remove part not in partition

      void hideAndClearArcsAbove(const idNode &baseNode);

      void hideAndClearArcsBelow(const idNode &baseNode, const SimplexId &seed);

      idSuperArc hideAndClearLeadingTo(const idNode &baseNode,
                                       const SimplexId &v);

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

      std::string printArc(const idSuperArc &a) {
        const SuperArc *sa = getSuperArc(a);
        std::stringstream res;
        res << a << ": ";
        if(sa->getDownCT() == treeData_.partition)
          res << getNode(sa->getDownNodeId())->getVertexId() << " -- ";
        else
          res << "(extern) -- ";

        if(sa->getUpCT() == treeData_.partition)
          res << getNode(sa->getUpNodeId())->getVertexId();
        else
          res << "(extern)";

        res << " \t\t(vis:" << sa->isVisible() << ")";
        return res.str();
      }

      std::string printNode(const idNode &n) {
        const Node *node = getNode(n);
        std::stringstream res;
        res << n << " : (";
        res << node->getVertexId() << ") / ";

        for(idSuperArc i = 0; i < node->getNumberOfUpSuperArcs(); ++i) {
          if(getSuperArc(node->getUpSuperArcId(i))->isVisible()) {
            res << "+";
          } else {
            res << "-";
          }
          res << node->getUpSuperArcId(i) << " ";
        }

        res << " \\ ";

        for(idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); ++i) {
          if(getSuperArc(node->getDownSuperArcId(i))->isVisible()) {
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

      inline bool isLower(const SimplexId &a, const SimplexId &b) const {
        return scalars_->isLower(a, b);
      }

      inline bool isHigher(const SimplexId &a, const SimplexId &b) const {
        return scalars_->isHigher(a, b);
      }

      //}

    private:
      // ------------------
      // Comparisons
      // -----------------
      // {
      // Compare using the scalar array : only for sort step

      template <typename scalarType>
      inline bool isLower(const SimplexId &a, const SimplexId &b) const {
        return ((scalarType *)scalars_->values)[a]
                 < ((scalarType *)scalars_->values)[b]
               || (((scalarType *)scalars_->values)[a]
                     == ((scalarType *)scalars_->values)[b]
                   && scalars_->sosOffsets[a] < scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline bool isHigher(const SimplexId &a, const SimplexId &b) const {
        return ((scalarType *)scalars_->values)[a]
                 > ((scalarType *)scalars_->values)[b]
               || (((scalarType *)scalars_->values)[a]
                     == ((scalarType *)scalars_->values)[b]
                   && scalars_->sosOffsets[a] > scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline bool isEqLower(const SimplexId &a, const SimplexId &b) const {
        return ((scalarType *)scalars_->values)[a]
                 < ((scalarType *)scalars_->values)[b]
               || (((scalarType *)scalars_->values)[a]
                     == ((scalarType *)scalars_->values)[b]
                   && scalars_->sosOffsets[a] <= scalars_->sosOffsets[b]);
      }

      template <typename scalarType>
      inline bool isEqHigher(const SimplexId &a, const SimplexId &b) const {
        return ((scalarType *)scalars_->values)[a]
                 > ((scalarType *)scalars_->values)[b]
               || (((scalarType *)scalars_->values)[a]
                     == ((scalarType *)scalars_->values)[b]
                   && scalars_->sosOffsets[a] >= scalars_->sosOffsets[b]);
      }

      // }
      // ----------------
      // Simplification
      // ----------------
      // {

      // preserve = do no hide it.
      void hideAndMerge(const idSuperArc &mergingArcId,
                        const idSuperArc &receptacleArcId,
                        const bool preserveDownNode = false);

      // Use BFS from root to find down and up of the receptarc (maintaining
      // segmentation information)
      std::tuple<idNode, idNode, SimplexId> createReceptArc(
        const idNode &root,
        const idSuperArc &receptArcId,
        std::vector<ExtendedUnionFind *> &arrayUF,
        const std::vector<std::pair<idSuperArc, idSuperArc>> &valenceOffsets);

      // during this BFS nodes should have only one arc up/down : find it :
      idSuperArc newUpArc(const idNode &curNodeId,
                          std::vector<ExtendedUnionFind *> &ufArray);

      idSuperArc newDownArc(const idNode &curNodeId,
                            std::vector<ExtendedUnionFind *> &ufArray);

      // }
      // --------------
      // Tool
      // --------------
      // {
      // create a std::pair with relative order : child vertex first

      inline std::pair<SimplexId, SimplexId>
        reorderEdgeRel(const std::pair<SimplexId, SimplexId> &vert) {
        if(treeData_.treeType == TreeType::Split) {
          if(isHigher(vert.first, vert.second)) {
            return vert;
          }

          return std::make_pair(vert.second, vert.first);
        } // else

        if(isLower(vert.first, vert.second))
          return vert;

        return std::make_pair(vert.second, vert.first);
      }

      template <typename triangulationType>
      bool verifyTree(const triangulationType &mesh);

      // Create a std::pair with the value corresponding to the simplification
      // method

      template <typename scalarType, typename triangulationType>
      void addPair(
        std::vector<std::tuple<SimplexId, SimplexId, scalarType, bool>> &pairs,
        const SimplexId &orig,
        const SimplexId &term,
        const triangulationType &mesh,
        const bool goUp);

      template <typename scalarType, typename triangulationType>
      void addPair(
        std::vector<std::tuple<SimplexId, SimplexId, scalarType>> &pairs,
        const SimplexId &orig,
        const SimplexId &term,
        const triangulationType &mesh);

      // }
    };

    std::ostream &operator<<(std::ostream &o, Node const &n);
    std::ostream &operator<<(std::ostream &o, SuperArc const &a);

  } // namespace cf
} // namespace ttk

#include <MergeTreeTemplate.h>
