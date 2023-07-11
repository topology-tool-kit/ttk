/// \ingroup base
/// \class ttk::FTMTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
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

#pragma once

#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <vector>

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include "FTMAtomicUF.h"
#include "FTMAtomicVector.h"
#include "FTMDataTypes.h"
#include "FTMNode.h"
#include "FTMStructures.h"
#include "FTMSuperArc.h"

static ttk::Timer _launchGlobalTime;

namespace ttk {
  namespace ftm {
    using UF = AtomicUF *;

    // Tree data ( 1 per tree )
    struct TreeData {
      TreeType treeType;

      // components : tree / nodes / extrema
      std::shared_ptr<FTMAtomicVector<SuperArc>> superArcs;
      std::shared_ptr<FTMAtomicVector<Node>> nodes;
      std::shared_ptr<FTMAtomicVector<idNode>> roots;
      std::vector<idNode> leaves;

      // vertex 2 node / superarc
      std::vector<idCorresp> vert2tree;
      std::vector<SimplexId> visitOrder;
      std::vector<std::list<std::vector<SimplexId>>> trunkSegments;

      // Track information
      std::vector<AtomicUF> storage;
      std::vector<UF> ufs;
      std::vector<UF> propagation;
      std::shared_ptr<FTMAtomicVector<CurrentState>> states;
      // valences
      std::vector<valence> valences;
      // opened nodes
      std::vector<char> openedNodes;

      // current nb of tasks
      idNode activeTasks;

      // Segmentation, stay empty for Contour tree as
      // they are created by Merge Tree
      Segments segments_;

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      std::vector<ActiveTask> activeTasksStats;
#endif

#ifdef TTK_ENABLE_OMP_PRIORITY
      // Is this MT to be computed with greater task priority than others
      bool prior = false;
#endif
    };

    class FTMTree_MT : virtual public Debug {
    protected:
      // global
      std::shared_ptr<Params> params_;
      std::shared_ptr<Scalars> scalars_;

      // local
      TreeData mt_data_;
      Comparison comp_;

    public:
      // -----------
      // CONSTRUCT
      // -----------

      // Tree with global data and partition number
      FTMTree_MT(const std::shared_ptr<Params> &params,
                 const std::shared_ptr<Scalars> &scalars,
                 TreeType type);

      ~FTMTree_MT() override;

      void clear();

      // --------------------
      // Init
      // --------------------

      inline void setParamsScalars(const std::shared_ptr<Params> &params,
                                   const std::shared_ptr<Scalars> &scalars) {
        this->scalars_ = scalars;
        this->params_ = params;
        this->mt_data_.treeType = params->treeType;
      }

      template <class triangulationType>
      void initNbScalars(const triangulationType *triangulation) {
        scalars_->size = triangulation->getNumberOfVertices();
      }

      void initComp() {
        if(isST()) {
          comp_.vertLower
            = [this](const SimplexId a, const SimplexId b) -> bool {
            return this->scalars_->isHigher(a, b);
          };
          comp_.vertHigher
            = [this](const SimplexId a, const SimplexId b) -> bool {
            return this->scalars_->isLower(a, b);
          };
        } else {
          comp_.vertLower
            = [this](const SimplexId a, const SimplexId b) -> bool {
            return this->scalars_->isLower(a, b);
          };
          comp_.vertHigher
            = [this](const SimplexId a, const SimplexId b) -> bool {
            return this->scalars_->isHigher(a, b);
          };
        }
      }

      bool compLower(const SimplexId a, const SimplexId b) {
        return comp_.vertLower(a, b);
      }

      /// \brief if sortedVertices_ is null, define and fill it
      template <typename scalarType>
      void sortInput();

      /// \brief clear local data for new computation
      void makeAlloc() {
        createAtomicVector<SuperArc>(mt_data_.superArcs);

        // Stats alloc

        createAtomicVector<Node>(mt_data_.nodes);
        mt_data_.nodes->reserve(scalars_->size / 2);

        createAtomicVector<idNode>(mt_data_.roots);
        mt_data_.roots->reserve(10);

        createVector<idNode>(mt_data_.leaves);
        mt_data_.leaves.reserve(scalars_->size / 3);

        // Known size

        createVector<idCorresp>(mt_data_.vert2tree);
        mt_data_.vert2tree.resize(scalars_->size);

        createVector<std::list<std::vector<SimplexId>>>(mt_data_.trunkSegments);

        createVector<SimplexId>(mt_data_.visitOrder);
        mt_data_.visitOrder.resize(scalars_->size);

        createVector<UF>(mt_data_.ufs);
        mt_data_.ufs.resize(scalars_->size);

        createVector<UF>(mt_data_.propagation);
        mt_data_.propagation.resize(scalars_->size);

        createVector<valence>(mt_data_.valences);
        mt_data_.valences.resize(scalars_->size);

        createVector<char>(mt_data_.openedNodes);
        mt_data_.openedNodes.resize(scalars_->size);

        mt_data_.segments_.clear();
      }

      void makeInit() {
        initVector<idCorresp>(mt_data_.vert2tree, nullCorresp);
        initVector<SimplexId>(mt_data_.visitOrder, nullVertex);
        initVector<UF>(mt_data_.ufs, nullptr);
        initVector<UF>(mt_data_.propagation, nullptr);
        initVector<valence>(mt_data_.valences, 0);
        initVector<char>(mt_data_.openedNodes, 0);
      }

      void initVectStates(const SimplexId nbLeaves) {
        if(!mt_data_.states) {
          mt_data_.states = std::make_shared<FTMAtomicVector<CurrentState>>(
            nbLeaves, comp_.vertHigher);
        }
        mt_data_.states->clear();
        mt_data_.states->reserve(nbLeaves);
      }

      // -------------------
      // Process
      // -------------------

      /// \brief Compute the merge
      template <class triangulationType>
      void build(const triangulationType *mesh, const bool ct);

      // extrema

      template <class triangulationType>
      int leafSearch(const triangulationType *mesh);

      // skeleton

      template <class triangulationType>
      void leafGrowth(const triangulationType *mesh);

      template <class triangulationType>
      void arcGrowth(const triangulationType *mesh,
                     const SimplexId startVert,
                     const SimplexId orig);

      template <class triangulationType>
      std::tuple<bool, bool> propagate(const triangulationType *mesh,
                                       CurrentState &currentState,
                                       UF curUF);

      template <class triangulationType>
      void closeAndMergeOnSaddle(const triangulationType *mesh,
                                 SimplexId saddleVert);

      template <class triangulationType>
      void closeOnBackBone(const triangulationType *mesh, SimplexId saddleVert);

      void closeArcsUF(idNode closeNode, UF uf);

      template <class triangulationType>
      SimplexId trunk(const triangulationType *mesh, const bool ct);

      virtual SimplexId
        trunkSegmentation(const std::vector<SimplexId> &pendingNodesVerts,
                          const SimplexId begin,
                          const SimplexId stop);

      // fill treedata_.trunkSegments
      SimplexId
        trunkCTSegmentation(const std::vector<SimplexId> &pendingNodesVerts,
                            const SimplexId begin,
                            const SimplexId stop);

      // segmentation

      /// \brief use vert2tree to compute the segmentation of the fresh built
      /// merge tree.
      void buildSegmentation();

      // Create the segmentation of all arcs by operating the pending operations
      void finalizeSegmentation();

      void normalizeIds();

      // -------------
      // ACCESSOR
      // ------------

      // Tree info for wrapper

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      const ActiveTask &getActiveTasks(const idSuperArc taskId) const {
        return (*mt_data_.activeTasksStats)[taskId];
      }
#endif

      inline SimplexId getArcSize(const idSuperArc arcId) {
        return getSuperArc(arcId)->size();
      }

      inline bool isJT() const {
        return mt_data_.treeType == TreeType::Join;
      }

      inline bool isST() const {
        return mt_data_.treeType == TreeType::Split;
      }

      // global
      // called for the tree used by the wrapper (only).
      // On this implementation, the warpper communicate with ContourForest
      // A child class of this one.

      inline void preconditionTriangulation(AbstractTriangulation *tri,
                                            const bool preproc = true) {
        if(tri && preproc) {
          // propagate through vertices (build)
          tri->preconditionVertexNeighbors();
        }
      }

      inline void setScalars(void *local_scalars) {
        scalars_->values = local_scalars;
      }

      inline void setTreeType(const int local_treeType) {
        params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setSegmentation(const bool segm) {
        params_->segm = segm;
      }

      inline void setNormalizeIds(const bool normalize) {
        params_->normalize = normalize;
      }

#ifdef TTK_ENABLE_OMP_PRIORITY
      inline void setPrior(void) {
        mt_data_.prior = true;
      }

      inline bool isPrior(void) const {
        return mt_data_.prior;
      }
#endif

      // scalar

      template <typename scalarType>
      inline const scalarType &getValue(SimplexId nodeId) const {
        return (((scalarType *)scalars_->values))[nodeId];
      }

      template <typename scalarType>
      inline void setVertexScalars(const scalarType *vals) {
        scalars_->values = static_cast<void *>(const_cast<scalarType *>(vals));
      }

      // offset
      /**
       * @pre For this function to behave correctly in the absence of
       * the VTK wrapper, ttk::preconditionOrderArray() needs to be
       * called to fill the @p sos buffer prior to any
       * computation (the VTK wrapper already includes a mechanism to
       * automatically generate such a preconditioned buffer).
       * @see examples/c++/main.cpp for an example use.
       */
      inline void setVertexSoSoffsets(const SimplexId *const sos) {
        scalars_->offsets = sos;
      }

      // arcs

      inline idSuperArc getNumberOfSuperArcs() const {
        return mt_data_.superArcs->size();
      }

      inline SuperArc *getSuperArc(idSuperArc i) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(i >= mt_data_.superArcs->size()) {
          std::cout << "[Merge Tree] get superArc on bad id :" << i;
          std::cout << " / " << mt_data_.superArcs->size() << std::endl;
        }
#endif
        return &((*mt_data_.superArcs)[i]);
      }

      inline const SuperArc *getSuperArc(idSuperArc i) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(i >= mt_data_.superArcs->size()) {
          std::cout << "[Merge Tree] get superArc on bad id :" << i;
          std::cout << " / " << mt_data_.superArcs->size() << std::endl;
        }
#endif
        return &((*mt_data_.superArcs)[i]);
      }

      // nodes

      inline idNode getNumberOfNodes() const {
        return mt_data_.nodes->size();
      }

      inline Node *getNode(idNode nodeId) {
        return &((*mt_data_.nodes)[nodeId]);
      }

      inline void setValence(const SimplexId v, const SimplexId val) {
        mt_data_.valences[v] = val;
      }

      // leaves / root

      inline idNode getNumberOfLeaves() const {
        return mt_data_.leaves.size();
      }

      inline const std::vector<idNode> &getLeaves() const {
        // break encapsulation...
        return mt_data_.leaves;
      }

      inline idNode getLeave(const idNode id) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(id > mt_data_.leaves.size()) {
          this->printErr("getLeaves out of bounds: " + std::to_string(id));
          return mt_data_.leaves[0];
        }
#endif
        return mt_data_.leaves[id];
      }

      inline const std::vector<idNode> &getRoots() const {
        // break encapsulation...
        return (*mt_data_.roots);
      }

      // vertices

      inline SimplexId getNumberOfVertices() const {
        return scalars_->size;
      }

      // vert2tree

      inline void setVert2Tree(decltype(mt_data_.vert2tree) const &vect2tree) {
        mt_data_.vert2tree = vect2tree;
      }

      // --------------------
      // VERT 2 TREE Special functions
      // --------------------

      // test vertex correpondance

      inline bool isCorrespondingArc(const SimplexId val) const {
        return !isCorrespondingNull(val) && mt_data_.vert2tree[val] >= 0;
      }

      inline bool isCorrespondingNode(const SimplexId val) const {
        return mt_data_.vert2tree[val] < 0;
      }

      inline bool isCorrespondingNull(const SimplexId val) const {
        return mt_data_.vert2tree[val] == nullCorresp;
      }

      // Get vertex info

      inline idNode getCorrespondingNodeId(const SimplexId val) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!isCorrespondingNode(val)) {
          this->printErr("getCorrespondingNode, Vertex: " + std::to_string(val)
                         + " is not a node: "
                         + std::to_string(mt_data_.vert2tree[val]));
        }
#endif
        return corr2idNode(val);
      }

      inline idSuperArc getCorrespondingSuperArcId(const SimplexId val) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(!isCorrespondingArc(val)) {
          this->printErr(
            "getCorrespondingSuperArcId, Vertex: " + std::to_string(val)
            + " is not on an arc: " + std::to_string(mt_data_.vert2tree[val]));
        }
#endif
        return mt_data_.vert2tree[val];
      }

      // Get corresponding element

      inline SuperArc *vertex2SuperArc(const SimplexId vert) {
        return &((*mt_data_.superArcs)[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const SimplexId vert) {
        return &((*mt_data_.nodes)[getCorrespondingNodeId(vert)]);
      }

      // Update vertex info

      inline void updateCorrespondingArc(const SimplexId vert,
                                         const idSuperArc arc) {
        mt_data_.vert2tree[vert] = arc;
      }

      inline void updateCorrespondingNode(const SimplexId vert,
                                          const idNode node) {
        mt_data_.vert2tree[vert] = idNode2corr(node);
      }

      inline idCorresp idNode2corr(const idNode id) const {
        // transform idNode to special value for the array : -idNode -1
        return -static_cast<idCorresp>(id + 1);
      }

      inline idNode corr2idNode(const idCorresp &corr) const {
        return static_cast<idNode>(-(mt_data_.vert2tree[corr] + 1));
      }

      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // SuperArcs

      idSuperArc openSuperArc(idNode downNodeId);

      idSuperArc makeSuperArc(idNode downNodeId, idNode upNodeId);

      void closeSuperArc(idSuperArc superArcId, idNode upNodeId);

      // Nodes

      std::vector<idNode> sortedNodes(const bool parallel = false);

      void sortLeaves(const bool parallel = false);

      /**
       * @brief Sort tree nodes according to vertex order
       *
       * The vertex order is the same for Join Trees and Split Trees:
       * minima first, maxima last.
       */
      void sortNodes();

      /**
       * @brief Sort tree arcs
       *
       * Arcs are sorted according to the lexicographic order (down
       * node order, up node order). The node order is the one used in
       * @ref sortNodes.
       */
      void sortArcs();

      idNode makeNode(SimplexId vertexId, SimplexId linked = nullVertex);

      idNode makeNode(const Node *const n, SimplexId linked = nullVertex);

      idSuperArc insertNode(Node *node, const bool segm = true);

      // get node starting / ending this arc
      // orientation depends on Join/Split tree
      Node *getDownNode(const SuperArc *a);
      Node *getUpNode(const SuperArc *a);
      idNode getDownNodeId(const SuperArc *a);
      idNode getUpNodeId(const SuperArc *a);

      // get node above / below this arc
      // in term of scalar value
      Node *getLowerNode(const SuperArc *a);
      Node *getUpperNode(const SuperArc *a);
      idNode getLowerNodeId(const SuperArc *a);
      idNode getUpperNodeId(const SuperArc *a);

      idNode getParent(const idNode n) {
        return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
      }

      void delNode(idNode node);

      // ---------------------------
      // Operators : clone/ move & print
      // ---------------------------

      std::shared_ptr<FTMTree_MT> clone() const;

      void move(FTMTree_MT &mt);

      // Print
      std::string printArc(idSuperArc a);

      std::string printNode(idNode n);

      void printTree2();

      void printParams() const;

      int printTime(Timer &t,
                    const std::string &s,
                    const int debugLevel = 2) const;

      // ----------------------------------------
      // Utils functions
      // Mathieu Pont (mathieu.pont@lip6.fr)
      // 2021
      // ----------------------------------------

      // --------------------
      // Is
      // --------------------
      bool isNodeOriginDefined(idNode nodeId);

      bool isRoot(idNode nodeId);

      bool isLeaf(idNode nodeId);

      bool isNodeAlone(idNode nodeId);

      bool isFullMerge();

      bool isBranchOrigin(idNode nodeId);

      template <class dataType>
      bool isJoinTree();

      template <class dataType>
      bool isImportantPair(idNode nodeId,
                           double threshold,
                           std::vector<double> &excludeLower,
                           std::vector<double> &excludeHigher);

      template <class dataType>
      bool isImportantPair(idNode nodeId, double threshold);

      bool isNodeMerged(idNode nodeId);

      bool isNodeIdInconsistent(idNode nodeId);

      bool isThereOnlyOnePersistencePair();

      // Do not normalize node is if root or son of a merged root
      bool notNeedToNormalize(idNode nodeId);

      bool isMultiPersPair(idNode nodeId);

      template <class dataType>
      bool isParentInconsistent(ftm::idNode nodeId);

      template <class dataType>
      bool verifyBranchDecompositionInconsistency();

      // --------------------
      // Get
      // --------------------
      idNode getRoot();

      idNode getParentSafe(idNode nodeId);

      void getChildren(idNode nodeId, std::vector<idNode> &res);

      void getLeavesFromTree(std::vector<idNode> &res);

      int getNumberOfLeavesFromTree();

      int getNumberOfNodeAlone();

      int getRealNumberOfNodes();

      template <class dataType>
      idNode getMergedRootOrigin();

      void getBranchOriginsFromThisBranch(
        idNode node, std::tuple<std::vector<idNode>, std::vector<idNode>> &res);

      void getTreeBranching(std::vector<idNode> &branching,
                            std::vector<int> &branchingID,
                            std::vector<std::vector<idNode>> &nodeBranching);

      void getTreeBranching(std::vector<idNode> &branching,
                            std::vector<int> &branchingID);

      void getAllRoots(std::vector<idNode> &res);

      int getNumberOfRoot();

      int getNumberOfChildren(idNode nodeId);

      int getTreeDepth();

      int getNodeLevel(idNode nodeId);

      void getAllNodeLevel(std::vector<int> &res);

      void getLevelToNode(std::vector<std::vector<idNode>> &res);

      void getBranchSubtree(std::vector<idNode> &branching,
                            idNode branchRoot,
                            std::vector<idNode> &res);

      template <class dataType>
      idNode getLowestNode(idNode nodeStart);

      // --------------------
      // Persistence
      // --------------------
      template <class dataType>
      std::tuple<dataType, dataType> getBirthDeathFromIds(idNode nodeId1,
                                                          idNode nodeId2);

      template <class dataType>
      std::tuple<dataType, dataType> getBirthDeathNodeFromIds(idNode nodeId1,
                                                              idNode nodeId2);

      template <class dataType>
      std::tuple<dataType, dataType> getBirthDeath(idNode nodeId);

      template <class dataType>
      std::tuple<ftm::idNode, ftm::idNode> getBirthDeathNode(idNode nodeId);

      template <class dataType>
      std::tuple<dataType, dataType> getMergedRootBirthDeath();

      template <class dataType>
      std::tuple<ftm::idNode, ftm::idNode> getMergedRootBirthDeathNode();

      template <class dataType>
      dataType getBirth(idNode nodeId);

      template <class dataType>
      dataType getNodePersistence(idNode nodeId);

      template <class dataType>
      dataType getMaximumPersistence();

      template <class dataType>
      ftm::idNode getSecondMaximumPersistenceNode();

      template <class dataType>
      dataType getSecondMaximumPersistence();

      template <class dataType>
      void getPersistencePairsFromTree(
        std::vector<std::tuple<ftm::idNode, ftm::idNode, dataType>> &pairs,
        bool useBD);

      template <class dataType>
      std::vector<ftm::idNode> getMultiPersOrigins(bool useBD);

      void getMultiPersOriginsVectorFromTree(
        std::vector<std::vector<idNode>> &res);

      // --------------------
      // Set
      // --------------------
      void setParent(idNode nodeId, idNode newParentNodeId);

      // --------------------
      // Delete
      // --------------------
      // Delete node by keeping subtree
      void deleteNode(idNode nodeId);

      void deleteIthUpArc(idNode nodeId, int arcIth);

      // Delete arc of the node to its parent
      void deleteParent(idNode nodeId);

      // Delete node without keeping subtree
      void deleteSubtree(idNode nodeId);

      // --------------------
      // Create/Delete/Modify Tree
      // --------------------
      void copyMergeTreeStructure(FTMTree_MT *tree);

      // --------------------
      // Utils
      // --------------------
      void printNodeSS(idNode node, std::stringstream &ss);

      template <class dataType>
      std::stringstream printNode2(idNode nodeId, bool doPrint = true);

      template <class dataType>
      std::stringstream printMergedRoot(bool doPrint = true);

      std::stringstream printTree(bool doPrint = true);

      std::stringstream printTreeStats(bool doPrint = true);

      template <class dataType>
      std::stringstream printTreeScalars(bool printNodeAlone = true,
                                         bool doPrint = true);

      template <class dataType>
      std::stringstream printPairsFromTree(bool useBD = false,
                                           bool printPairs = true,
                                           bool doPrint = true);

      std::stringstream printMultiPersOriginsVectorFromTree(bool doPrint
                                                            = true);

      template <class dataType>
      std::stringstream printMultiPersPairsFromTree(bool useBD = false,
                                                    bool printPairs = true,
                                                    bool doPrint = true);

      // ----------------------------------------
      // End of utils functions
      // ----------------------------------------

    protected:
      // -----
      // Tools
      // -----

      idNode getVertInRange(const std::vector<SimplexId> &range,
                            const SimplexId v,
                            const idNode last = 0) const;

      std::tuple<SimplexId, SimplexId>
        getBoundsFromVerts(const std::vector<SimplexId> &nodes) const;

      idSuperArc upArcFromVert(const SimplexId v) {
        return getNode(getCorrespondingNodeId(v))->getUpSuperArcId(0);
      }

      inline SimplexId getChunkSize(const SimplexId nbVerts = -1,
                                    const SimplexId nbtasks = 100) const {
        const SimplexId s = (nbVerts == -1) ? scalars_->size : nbVerts;
#ifndef NDEBUG
        // Debug mode
        static const SimplexId minWorks = 1;
#else
        // Release mode
        static const SimplexId minWorks = 10000;
#endif
        return std::max(minWorks, 1 + (s / (nbtasks * threadNumber_)));
      }

      inline SimplexId getChunkCount(const SimplexId nbVerts = -1,
                                     const SimplexId nbTasks = 100) const {
        const SimplexId s = (nbVerts == -1) ? scalars_->size : nbVerts;
        return 1 + (s / getChunkSize(s, nbTasks));
      }

      void sortUpArcs(const idNode nid) {
        auto comp = [&](const idSuperArc a, const idSuperArc b) -> bool {
          return comp_.vertLower(getUpperNode(getSuperArc(a))->getVertexId(),
                                 getUpperNode(getSuperArc(b))->getVertexId());
        };

        getNode(nid)->sortUpArcs(comp);
      }

      void sortDownArcs(const idNode nid) {
        auto comp = [&](const idSuperArc a, const idSuperArc b) -> bool {
          return comp_.vertHigher(getUpperNode(getSuperArc(a))->getVertexId(),
                                  getUpperNode(getSuperArc(b))->getVertexId());
        };

        getNode(nid)->sortDownArcs(comp);
      }

      // ------------------
      // Comparisons
      // -----------------
      // Compare using the scalar array : only for sort step

      inline bool isLower(SimplexId a, SimplexId b) const {
        return scalars_->offsets[a] < scalars_->offsets[b];
      }

      inline bool isHigher(SimplexId a, SimplexId b) const {
        return scalars_->offsets[a] > scalars_->offsets[b];
      }

      template <typename type>
      void createVector(std::vector<type> &vec) {
        vec.clear();
      }

      template <typename type>
      void createAtomicVector(std::shared_ptr<FTMAtomicVector<type>> &ptr) {
        if(!ptr)
          ptr = std::make_shared<FTMAtomicVector<type>>();
        ptr->clear();
      }

      template <typename type>
      void initVector(std::vector<type> &vect, const type val) {
        auto s = vect.size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static)
#endif
        for(typename std::vector<type>::size_type i = 0; i < s; i++) {
          vect[i] = val;
        }
      }
    }; // end of FTMTree_MT class

    std::ostream &operator<<(std::ostream &o, Node const &n);
    std::ostream &operator<<(std::ostream &o, SuperArc const &a);

    template <typename dataType>
    struct MergeTree {
      std::shared_ptr<ftm::Scalars> scalars;
      std::shared_ptr<std::vector<dataType>> scalarsValues;
      std::shared_ptr<ftm::Params> params;
      ftm::FTMTree_MT tree;

      std::shared_ptr<ftm::Scalars> emptyScalars() {
        auto scalarsT = std::make_shared<ftm::Scalars>();
        scalarsT->size = 0;
        scalarsT->values = nullptr;
        return scalarsT;
      }

      std::shared_ptr<ftm::Params> emptyParams() {
        auto paramsT = std::make_shared<ftm::Params>();
        paramsT->treeType = ftm::Join_Split;
        return paramsT;
      }

      MergeTree() : MergeTree(emptyScalars(), emptyParams()) {
      }

      template <typename T, typename U>
      MergeTree(const T scalarsT, U paramsT)
        : scalars(scalarsT), params(paramsT),
          tree(paramsT, scalarsT, params->treeType) {
        tree.makeAlloc();
        scalarsValues = std::make_shared<std::vector<dataType>>();
        for(unsigned int i = 0; i < tree.getNumberOfNodes(); ++i)
          scalarsValues->push_back(tree.getValue<dataType>(i));
        scalars->values = (void *)(scalarsValues->data());
      }

      MergeTree(const std::shared_ptr<ftm::Scalars> &scalarsT,
                const std::shared_ptr<std::vector<dataType>> &scalarValuesT,
                std::shared_ptr<ftm::Params> &paramsT)
        : scalars(scalarsT), scalarsValues(scalarValuesT), params(paramsT),
          tree(paramsT, scalarsT, params->treeType) {
        tree.makeAlloc();
        scalars->values = (void *)(scalarsValues->data());
      }

      void copy(const MergeTree<dataType> &mt) {
        // Copy scalars
        scalars = std::make_shared<ftm::Scalars>();
        scalars->size = mt.scalars->size;
        scalarsValues = mt.scalarsValues;
        scalars->values = (void *)(scalarsValues->data());

        // Copy params
        params = std::make_shared<ftm::Params>();
        params->treeType = mt.params->treeType;

        // Copy tree
        tree.clear();
        tree.setParamsScalars(params, scalars);
        tree.makeAlloc();
        tree.copyMergeTreeStructure(const_cast<FTMTree_MT *>(&(mt.tree)));
      }

      MergeTree(const MergeTree<dataType> &mt)
        : scalars(mt.scalars), scalarsValues(mt.scalarsValues),
          params(mt.params), tree(params, scalars, params->treeType) {
        copy(mt);
      }

      MergeTree<dataType> &operator=(const MergeTree<dataType> &mt) {
        if(&mt != this) {
          copy(mt);
        }
        return *this;
      }
    };

  } // namespace ftm
} // namespace ttk

#include <FTMTreeUtils_Template.h>
#include <FTMTree_MT_Template.h>
