/// \ingroup base
/// \class ttk::ftr::Graph
/// \author Gueunet Charles <charles.gueunet+ttk@gmail.com>
/// \date 2018-01-25
///
/// \brief TTK %FTRGraph graph skeleton
///
/// This class manage nodes and arcs of the graph structure
/// along with their Segmentation
///
/// \sa ttk::FTRGraph

#ifndef GRAPH_H
#define GRAPH_H

#include "FTRAtomicVector.h"
#include "FTRCommon.h"
#include "FTRDataTypes.h"
#include "FTRNode.h"
#include "FTRPropagation.h"
#include "FTRScalars.h"
#include "FTRSuperArc.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif

#include <map>
#include <random>
#include <vector>

namespace ttk {
  namespace ftr {
    struct SegmInfo {
      idNode corNode = nullNode;
      idSuperArc corArc = nullSuperArc;
    };

    class Graph : public Allocable {
    private:
      // update operator =
      FTRAtomicVector<std::tuple<idVertex, bool>> leaves_;
      FTRAtomicVector<Node> nodes_;
      FTRAtomicVector<SuperArc> arcs_;

      std::vector<SegmInfo> segmentation_;

#ifdef TTK_ENABLE_FTR_VERT_STATS
      std::vector<idVertex> nbTouch_;
      std::vector<idSuperArc> nbArcActif_;
      idVertex avoided_;
#endif

    public:
      // For openmp capture, we need direct access...
      std::vector<valence> valDown_, valUp_;

      Graph();
      Graph(Graph &&other) noexcept = default;
      Graph(const Graph &other) = delete;
      virtual ~Graph();

      Graph &operator=(Graph &&other) noexcept {
        if(this != &other) {
          leaves_ = std::move(other.leaves_);
          nodes_ = std::move(other.nodes_);
          arcs_ = std::move(other.arcs_);
          segmentation_ = std::move(other.segmentation_);
          valDown_ = std::move(other.valDown_);
          valUp_ = std::move(other.valUp_);
#ifdef TTK_ENABLE_FTR_VERT_STATS
          nbTouch_ = std::move(other.nbTouch_);
          nbArcActif_ = std::move(other.nbArcActif_);
          avoided_ = std::move(other.avoided_);
#endif
        }
        return *this;
      }

      Graph &operator=(Graph &other) = delete;

      // Accessor on structure
      // ---------------------

      idNode getNumberOfNodes(void) const {
        return nodes_.size();
      }

      idSuperArc getNumberOfArcs(void) const {
        return arcs_.size();
      }

      idSuperArc getNumberOfVisibleArcs(void) const {
        idSuperArc res = 0;
        for(const auto &arc : arcs_) {
          if(arc.isVisible())
            ++res;
        }
        return res;
      }

      idNode getNumberOfLeaves(void) const {
        return leaves_.size();
      }

      idVertex getLeaf(const idNode id) const {
        return std::get<0>(leaves_[id]);
      }

      bool isLeafFromMin(const idNode id) const {
        return std::get<1>(leaves_[id]);
      }

      const Node &getNode(const idNode id) const {
        return nodes_[id];
      }

      Node &getNode(const idNode id) {
        return nodes_[id];
      }

      const SuperArc &getArc(const idSuperArc id) const {
        return arcs_[id];
      }

      SuperArc &getArc(const idSuperArc id) {
        return arcs_[id];
      }

      bool isVisited(const idVertex v) const {
        return segmentation_[v].corArc != nullSuperArc
               || segmentation_[v].corNode != nullNode;
      }

      void visit(const idVertex v, const idSegmentation id, bool isArc = true) {
        if(isArc) {
          if(segmentation_[v].corArc == nullSuperArc) {
            segmentation_[v].corArc = id;
          }
        } else {
          if(segmentation_[v].corNode == nullNode) {
            segmentation_[v].corNode = id;
          }
        }
      }

#ifdef TTK_ENABLE_FTR_VERT_STATS
      void incTouch(const idVertex v) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
        nbTouch_[v]++;
      }

      idVertex getNbTouch(const idVertex v) const {
        return nbTouch_[v];
      }

      void setNbArcActive(const idVertex v, const idSuperArc nb) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
        nbArcActif_[v] = nb;
      }

      idSuperArc getNbArcActive(const idVertex v) const {
        return nbArcActif_[v];
      }

      void incAvoid() {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
        avoided_++;
      }

      idVertex getNbAvoided() const {
        return avoided_;
      }

      idVertex getNbMultiTouch() const {
        auto gt1 = [](uint i) { return i > 1; };
        return std::count_if(nbTouch_.cbegin(), nbTouch_.cend(), gt1);
      }
#endif

      void setArc(const idVertex v, const idSegmentation id) {
        segmentation_[v].corArc = id;
      }

      std::tuple<idNode, idSuperArc> visit(const idVertex v) const {
        return {segmentation_[v].corNode, segmentation_[v].corArc};
      }

      bool isArc(const idVertex v) const {
        return segmentation_[v].corArc != nullSuperArc;
      }

      bool isArc(const idVertex v, const idSuperArc id) const {
        return segmentation_[v].corArc == id;
      }

      bool isNode(const idVertex v) const {
        return segmentation_[v].corNode != nullNode;
      }

      bool isNode(const idVertex v, const idNode id) const {
        return segmentation_[v].corNode != id;
      }

      idNode getNodeId(const idVertex v) const {
        return segmentation_[v].corNode;
      }

      idSuperArc getArcId(const idVertex v) const {
        return segmentation_[v].corArc;
      }

      const Node &getDownNode(const idSuperArc a) const {
        return getNode(getArc(a).getDownNodeId());
      }

      // direct access for openmp capture
      const valence &valUp(const idVertex v) const {
        return valUp_[v];
      }

      valence &valUp(const idVertex v) {
        return valUp_[v];
      }

      const valence &valDown(const idVertex v) const {
        return valDown_[v];
      }

      valence &valDown(const idVertex v) {
        return valDown_[v];
      }

      // Build structure
      // ---------------

      void addLeaf(const idVertex v, bool isMax) {
        leaves_.emplace_back(std::make_tuple(v, isMax));
      }

      // return either a new node or the existing one if this vertex
      // already corresponds to a node.
      idNode makeNode(const idVertex v) {
        return std::get<0>(getOrCreateNode(v));
      }

      // create a node in v if no node exists yet, then
      // return the id of the node on v.
      // The bool is true if a node has been created
      // This method is not thread safe
      std::tuple<idNode, bool> getOrCreateNode(const idVertex v) {
        if(isNode(v)) {
          return {getNodeId(v), false};
        }

        const idNode newNode = nodes_.getNext();
        nodes_[newNode].setVerterIdentifier(v);
        visit(v, newNode, false);
        return {newNode, true};
      }

      idSuperArc openArc(const idNode downId, Propagation *p = nullptr) {
        idSuperArc newArc = arcs_.getNext();
        arcs_[newArc].setDownNodeId(downId);
        if(p) {
          arcs_[newArc].setUfProp(p->getId());
        }
#ifndef TTK_ENABLE_KAMIKAZE
        else {
          std::cout << "NULL PROP IN ARC " << static_cast<unsigned>(newArc)
                    << std::endl;
        }
#endif
        return newArc;
      }

      void closeArc(const idSuperArc arc, const idNode upId) {
        arcs_[arc].setUpNodeId(upId);
      }

      void makeArc(const idNode downId, const idNode upId) {
        arcs_.emplace_back(SuperArc{downId, upId});
      }

      idSuperArc makeHiddenArc(Propagation *const lp) {
        idSuperArc newArc = arcs_.getNext();
        arcs_[newArc].hide();
        arcs_[newArc].setUfProp(lp->getId());
        return newArc;
      }

      // Process
      // -------

      // sort leaves vector by scalar value,
      // can be done in parallel
      template <typename ScalarType>
      void sortLeaves(const Scalars<ScalarType> *s,
                      const bool parallel = false) {
        auto compare_fun
          = [&](std::tuple<idVertex, bool> a, std::tuple<idVertex, bool> b) {
              return s->isLower(std::get<0>(a), std::get<0>(b));
            };
        if(parallel) {
          TTK_PSORT(
            this->threadNumber_, leaves_.begin(), leaves_.end(), compare_fun);
        } else {
          std::sort(leaves_.begin(), leaves_.end(), compare_fun);
        }
      }

      void shuffleLeaves() {
        std::shuffle(leaves_.begin(), leaves_.end(), std::random_device());
      }

      // some arc may be pending due to symbolic merge during computation
      // if tasks from both min and max.
      // here we replace them by one consistent arc
      template <typename ScalarType>
      void mergeArcs(const Scalars<ScalarType> *const s);

      // Link nodes to arcs when arc are completely created
      template <typename ScalarType>
      void arcs2nodes(const Scalars<ScalarType> *const s);

      // Build the list of regular vertices of each arc
      template <typename ScalarType>
      void buildArcSegmentation(const Scalars<ScalarType> *const s);

      // Tools
      // -----

      std::string print(const int verbosity) const;

      std::string printArc(const idSuperArc arcId) const;

      std::string printNode(const idNode nodeId) const;

      std::string printVisit(const idVertex v) const;

      // DEBUG function
      std::string printVisit() const;

      // Initialize functions
      // --------------------

      void alloc() override;

      void init() override;

    private:
      // tools

      // ensure that main arc have valid up/down node even if the merge of
      // the two arc occured during the computation, leaving some unfinished
      // arcs.
      template <typename ScalarType>
      void consolidateArc(const idSuperArc mainArc,
                          const idSuperArc mergedArc,
                          const Scalars<ScalarType> *const s);
    };

  } // namespace ftr
} // namespace ttk

#include "Graph_Template.h"

#endif /* end of include guard: GRAPH_H */
