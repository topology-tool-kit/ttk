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

#ifndef FTMTREE_MT_H
#define FTMTREE_MT_H

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
#include "FTMTree_DataTypes.h"
#include "Node.h"
#include "Structures.h"
#include "SuperArc.h"

namespace ttk
{
namespace ftm
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
      std::vector<idNode> *   leaves;

      // vertex 2 node / superarc
      std::vector<idCorresp> *vert2tree;
      std::vector<idVertex>  *visitOrder;
      std::vector<std::list<std::vector<idVertex>>> *trunkSegments;

      // Track informations
      std::vector<UF> *ufs, *propagation;
      AtomicVector<CurrentState> *states;
      // valences
      std::vector<valence> *valences;
      // opened nodes
      std::vector<char> *openedNodes;

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      std::vector<ActiveTask> *activeTasksStats;
#endif

      // current nb of tasks
      idNode activeTasks;

      // Segmentation, stay empty for Contour tree as
      // they are created by Merge Tree
      Segments segments_;
   };

   class FTMTree_MT : virtual public Debug
   {
     protected:
      // global
      Params *const  params_;
      Triangulation *mesh_;
      Scalars *const scalars_;

      // local
      TreeData   mt_data_;
      Comparison comp_;

     public:

      // -----------
      // CONSTRUCT
      // -----------

      // Tree with global data and partition number
      FTMTree_MT(Params *const params, Triangulation *mesh, Scalars *const scalars, TreeType type);

      virtual ~FTMTree_MT();

      // --------------------
      // Init
      // --------------------

      void initNbScalars(void)
      {
         scalars_->size = mesh_->getNumberOfVertices();
      }

      /// \brief init Simulation of Simplicity datastructure if not set
      template<typename idType>
      void initSoS(void)
      {
         if (scalars_->offsets == nullptr) {
            scalars_->offsets = new idType[scalars_->size];
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
            for (idVertex i = 0; i < scalars_->size; i++) {
               ((idType *)scalars_->offsets)[i] = i;
            }
         }

      }

      void initComp(void)
      {
         if (isST()) {
            comp_.vertLower = [this](const idVertex a, const idVertex b) -> bool {
               return this->scalars_->isHigher(a, b);
            };
            comp_.vertHigher = [this](const idVertex a, const idVertex b) -> bool {
               return this->scalars_->isLower(a, b);
            };
         } else {
            comp_.vertLower = [this](const idVertex a, const idVertex b) -> bool {
               return this->scalars_->isLower(a, b);
            };
            comp_.vertHigher = [this](const idVertex a, const idVertex b) -> bool {
               return this->scalars_->isHigher(a, b);
            };
         }
      }

      bool compLower(const idVertex a, const idVertex b)
      {
         return comp_.vertLower(a, b);
      }

      /// \brief if sortedVertices_ is null, define and fill it
      /// Also fill the mirror vector
      template <typename scalarType, typename idType>
      void sortInput(void);

      /// \brief clear local data for new computation
      void makeAlloc(void)
      {
         createAtomicVector<SuperArc>(mt_data_.superArcs);

         // Stats alloc

         createAtomicVector<Node>(mt_data_.nodes);
         mt_data_.nodes->reserve(scalars_->size/2);

         createAtomicVector<idNode>(mt_data_.roots);
         mt_data_.roots->reserve(10);

         createVector<idNode>(mt_data_.leaves);
         mt_data_.leaves->reserve(scalars_->size/3);

         // Known size

         createVector<idCorresp>(mt_data_.vert2tree);
         mt_data_.vert2tree->resize(scalars_->size);

         createVector<std::list<std::vector<idVertex>>>(mt_data_.trunkSegments);

         createVector<idVertex>(mt_data_.visitOrder);
         mt_data_.visitOrder->resize(scalars_->size);

         createVector<UF>(mt_data_.ufs);
         mt_data_.ufs->resize(scalars_->size);

         createVector<UF>(mt_data_.propagation);
         mt_data_.propagation->resize(scalars_->size);

         createVector<valence>(mt_data_.valences);
         mt_data_.valences->resize(scalars_->size);

         createVector<char>(mt_data_.openedNodes);
         mt_data_.openedNodes->resize(scalars_->size);

         mt_data_.segments_.clear();
      }

      void makeInit(void)
      {
         initVector<idCorresp>(mt_data_.vert2tree, nullCorresp);
         initVector<idVertex>(mt_data_.visitOrder, nullVertex);
         initVector<UF>(mt_data_.ufs, nullptr);
         initVector<UF>(mt_data_.propagation, nullptr);
         initVector<valence>(mt_data_.valences, 0);
         initVector<char>(mt_data_.openedNodes, 0);
      }

      void initVectStates(const idVertex nbLeaves)
      {
         if(!mt_data_.states) {
            mt_data_.states = new AtomicVector<CurrentState>(nbLeaves, comp_.vertHigher);
         }
         mt_data_.states->clear();
         mt_data_.states->reserve(nbLeaves);
      }

      // -------------------
      // Process
      // -------------------

      /// \brief Compute the merge
      void build(const bool ct);

      // extrema

      virtual int leafSearch();

      // skeleton

      void leafGrowth();

      void arcGrowth(const idVertex startVert, const idVertex orig);

      std::tuple<bool, bool> propage(CurrentState &currentState, UF curUF);

      void closeAndMergeOnSaddle(idVertex saddleVert);

      void closeOnBackBone(idVertex saddleVert);

      void closeArcsUF(idNode closeNode, UF uf);

      idVertex trunk(const bool ct);

      virtual idVertex trunkSegmentation(const std::vector<idVertex> &pendingNodesVerts,
                                         const idVertex begin,
                                         const idVertex stop);

      // fill treedata_.trunkSegments
      idVertex trunkCTSegmentation(const std::vector<idVertex> &pendingNodesVerts,
                                   const idVertex begin,
                                   const idVertex stop);

      // segmentation

      /// \brief use vert2tree to compute the segmentation of the fresh builded merge tree.
      void buildSegmentation();

      // Create the segmentation of all arcs by operating the pending operations
      void finalizeSegmentation(void);

      void normalizeIds();

      // -------------
      // ACCESSOR
      // ------------

      // Tree info for wrapper

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
      const ActiveTask& getActiveTasks(const idSuperArc taskId) const
      {
         return (*mt_data_.activeTasksStats)[taskId];
      }
#endif

      inline idVertex getArcSize(const idSuperArc arcId)
      {
          return getSuperArc(arcId)->size();
      }

      inline bool isJT(void) const
      {
         return mt_data_.treeType == TreeType::Join;
      }

      inline bool isST(void) const
      {
         return mt_data_.treeType == TreeType::Split;
      }

      // global
      // called for the tree used by the wrapper (only).
      // On this implementation, the warpper communicate with ContourForest
      // A child class of this one.

      inline void setupTriangulation(Triangulation *m, const bool preproc = true)
      {
         mesh_ = m;
         if (mesh_ && preproc) {
            // propage through vertices (build)
            mesh_->preprocessVertexNeighbors();
         }
      }

      inline void setScalars(void *local_scalars)
      {
         scalars_->values = local_scalars;
      }

      inline void setTreeType(const int local_treeType)
      {
         params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setSegmentation(const bool segm)
      {
          params_->segm = segm;
      }

      inline void setNormalizeIds(const bool normalize)
      {
          params_->normalize = normalize;
      }

      // scalar

      template <typename scalarType>
      inline const scalarType &getValue(idVertex idNode) const
      {
         return (((scalarType *)scalars_->values))[idNode];
      }

      inline void setVertexScalars(void *vals)
      {
         scalars_->values = vals;
      }

      // offset
      inline void setVertexSoSoffsets(void *sos)
      {
         scalars_->offsets = sos;
      }

      // arcs

      inline idSuperArc getNumberOfSuperArcs(void) const
      {
         return mt_data_.superArcs->size();
      }

      inline SuperArc *getSuperArc(idSuperArc i)
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)i >= mt_data_.superArcs->size()) {
            std::cout << "[Merge Tree] get superArc on bad id :" << i;
            std::cout << " / " << mt_data_.superArcs->size() << std::endl;
            return nullptr;
         }
#endif
         return &((*mt_data_.superArcs)[i]);
      }

      inline const SuperArc *getSuperArc(idSuperArc i) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)i >= mt_data_.superArcs->size()) {
            std::cout << "[Merge Tree] get superArc on bad id :" << i;
            std::cout << " / " << mt_data_.superArcs->size() << std::endl;
            return nullptr;
         }
#endif
         return &((*mt_data_.superArcs)[i]);
      }

      // nodes

      inline idNode getNumberOfNodes(void) const
      {
         return mt_data_.nodes->size();
      }

      inline Node *getNode(idNode nodeId)
      {
         return &((*mt_data_.nodes)[nodeId]);
      }

      inline void setValence(const idVertex v, const idVertex val)
      {
          (*mt_data_.valences)[v] = val;
      }

      // leaves / root

      inline idNode getNumberOfLeaves(void) const
      {
         return mt_data_.leaves->size();
      }

      inline const std::vector<idNode> &getLeaves(void) const
      {
         // break encapsulation...
         return (*mt_data_.leaves);
      }

      inline idNode getLeave(const idNode id) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)id > (mt_data_.leaves->size())) {
            std::stringstream msg;
            msg << "[MergTree] getLeaves out of bounds : " << id << std::endl;
            err(msg.str(), fatalMsg);
            return (*mt_data_.leaves)[0];
         }
#endif
         return (*mt_data_.leaves)[id];
      }

      inline const std::vector<idNode> &getRoots(void) const
      {
         // break encapsulation...
         return (*mt_data_.roots);
      }

      // vertices

      inline idVertex getNumberOfVertices(void) const
      {
         return scalars_->size;
      }

      // vert2tree

      inline void setVert2Tree(decltype(mt_data_.vert2tree) const vect2tree)
      {
         mt_data_.vert2tree = vect2tree;
      }

      // --------------------
      // VERT 2 TREE Special functions
      // --------------------

      // test vertex correpondance

      inline bool isCorrespondingArc(const idVertex val) const
      {
         return !isCorrespondingNull(val) && (*mt_data_.vert2tree)[val] >= 0;
      }

      inline bool isCorrespondingNode(const idVertex val) const
      {
         return (*mt_data_.vert2tree)[val] < 0;
      }

      inline bool isCorrespondingNull(const idVertex val) const
      {
         return (*mt_data_.vert2tree)[val] == nullCorresp;
      }

      // Get vertex info

      inline idNode getCorrespondingNodeId(const idVertex val) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!isCorrespondingNode(val)) {
            std::stringstream debug;
            debug << "[FTMTree_MT] : getCorrespondingNode, ";
            debug << "Vertex :" << val << " is not a node :";
            debug << (*mt_data_.vert2tree)[val] << std::endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return corr2idNode(val);
      }

      inline idSuperArc getCorrespondingSuperArcId(const idVertex val) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!isCorrespondingArc(val)) {
            std::stringstream debug;
            debug << "[FTMTree_MT] : getCorrespondingSuperArcId, ";
            debug << "Vertex :" << val << " is not on an arc :";
            debug << (*mt_data_.vert2tree)[val] << std::endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return (*mt_data_.vert2tree)[val];
      }

      // Get corresponding elemnt

      inline SuperArc *vertex2SuperArc(const idVertex vert)
      {
         return &((*mt_data_.superArcs)[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const idVertex vert)
      {
         return &((*mt_data_.nodes)[getCorrespondingNodeId(vert)]);
      }

      // Update vertex info


      inline void updateCorrespondingArc(const idVertex vert, const idSuperArc arc)
      {
         (*mt_data_.vert2tree)[vert] = arc;
      }

      inline void updateCorrespondingNode(const idVertex vert, const idNode node)
      {
         (*mt_data_.vert2tree)[vert] = idNode2corr(node);
      }

      inline idCorresp idNode2corr(const idNode id) const
      {
         // transform idNode to special value for the array : -idNode -1
         return -(idCorresp)(id + 1);
      }

      inline idNode corr2idNode(const idCorresp &corr) const
      {
         return -(idNode)((*mt_data_.vert2tree)[corr] + 1);
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

      idNode makeNode(idVertex vertexId, idVertex linked = nullVertex);

      idNode makeNode(const Node *const n, idVertex linked = nullVertex);

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
      Node * getUpperNode(const SuperArc *a);
      idNode getLowerNodeId(const SuperArc *a);
      idNode getUpperNodeId(const SuperArc *a);

      idNode getParent(const idNode n)
      {
         return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
      }

      void delNode(idNode node);

      // ---------------------------
      // Operators : clone/ move & print
      // ---------------------------

      FTMTree_MT *clone() const;

      void move(FTMTree_MT *mt);

      // Print
      std::string printArc(idSuperArc a);

      std::string printNode(idNode n);

      void printTree2(void);

      void printParams(void) const;

      int printTime(DebugTimer &t, const std::string &s, idVertex nbScalars = -1,
                    const int debugLevel = 2) const;

     protected:

      // -----
      // Tools
      // -----

      idNode getVertInRange(const std::vector<idVertex> &range,
                            const idVertex v,
                            const idNode last = 0) const;

      std::tuple<idVertex, idVertex> getBoundsFromVerts(const std::vector<idVertex> &nodes) const;

      idSuperArc upArcFromVert(const idVertex v)
      {
         return getNode(getCorrespondingNodeId(v))->getUpSuperArcId(0);
      }

      inline idVertex getChunkSize(const idVertex nbVerts = -1, const idVertex nbtasks = 100) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->size : nbVerts;
#ifndef NDEBUG
         // Debug mode
         static const idVertex minWorks = 1;
#else
         // Release mode
         static const idVertex minWorks = 10000;
#endif
         return std::max(minWorks, 1 + (s / (nbtasks * threadNumber_)));
      }

      inline idVertex getChunkCount(const idVertex nbVerts = -1, const idVertex nbTasks = 100) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->size : nbVerts;
         return 1 + (s / getChunkSize(s, nbTasks));
      }

      void sortUpArcs(const idNode nid)
      {
         auto comp = [&](const idSuperArc a, const idSuperArc b) -> bool {
            return comp_.vertLower(getUpperNode(getSuperArc(a))->getVertexId(),
                                   getUpperNode(getSuperArc(b))->getVertexId());
         };

         getNode(nid)->sortUpArcs(comp);
      }

      void sortDownArcs(const idNode nid)
      {
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

      template <typename scalarType,typename idType>
      inline bool isLower(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 ((idType *)scalars_->offsets)[a] < ((idType *)scalars_->offsets)[b]);
      }

      template <typename scalarType,typename idType>
      inline bool isHigher(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 ((idType *)scalars_->offsets)[a] > ((idType *)scalars_->offsets)[b]);
      }

      template <typename scalarType,typename idType>
      inline bool isEqLower(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] < ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 ((idType *)scalars_->offsets)[a] <= ((idType *)scalars_->offsets)[b]);
      }

      template <typename scalarType,typename idType>
      inline bool isEqHigher(idVertex a, idVertex b) const
      {
         return ((scalarType *)scalars_->values)[a] > ((scalarType *)scalars_->values)[b] ||
                (((scalarType *)scalars_->values)[a] == ((scalarType *)scalars_->values)[b] &&
                 ((idType *)scalars_->offsets)[a] >= ((idType*)scalars_->offsets)[b]);
      }

      template <typename type>
      void createVector(std::vector<type> *&ptr)
      {
         if(!ptr)
            ptr = new std::vector<type>;
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
      void initVector(std::vector<type> *&vect, const type val)
      {
         int s = vect->size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static)
#endif
         for (int i = 0; i < s; i++) {
            (*vect)[i] = val;
         }
      }

   };

   std::ostream &operator<<(std::ostream &o, Node const &n);
   std::ostream &operator<<(std::ostream &o, SuperArc const &a);


}
}

#include <FTMTree_MT_Template.h>

#endif /* end of include guard: MERGETREE_H */
