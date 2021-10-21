#ifndef GRAPH_TEMPLATE_H
#define GRAPH_TEMPLATE_H

#ifndef NDEBUG
// #include<sstream>
// #define DEBUG(msg) {std::stringstream s; s << msg << std::endl; std::cout <<
// s.str();}
#define DEBUG(msg)
#else
#define DEBUG(msg)
#endif

#include "Graph.h"

#include <unordered_map>

namespace ttk {
  namespace ftr {
    template <typename ScalarType>
    void Graph::mergeArcs(const Scalars<ScalarType> *const s) {
      std::unordered_map<idSuperArc, idSuperArc> mapArcs;
      std::map<std::pair<idVertex, idVertex>, idSuperArc> masterArcs;

      // int totalArc = getNumberOfArcs();
      int merged = 0;

      // Arc created by increasing and decreasing tasks are reversed in
      // terms of up/down node. We use this property to merge such arcs. in
      // order to distinguish between arcs around a loop, and arc is
      // described by the fisrt and last vertex of its segmentation if any.
      // Otherwise by its up / last node.

      const idSuperArc nbArcs = arcs_.size();
      for(idSuperArc arcId = 0; arcId < nbArcs; ++arcId) {
        const SuperArc &arc = getArc(arcId);
        if(getArc(arcId).isVisible()) {
          idNode upNodeId = getArc(arcId).getUpNodeId();
          const idNode downNodeId = getArc(arcId).getDownNodeId();
          if(upNodeId == nullNode) {
            if(getArc(arcId).getEnd() != nullVertex) {
              upNodeId = getNodeId(getArc(arcId).getEnd());
              if(upNodeId == nullNode) {
                upNodeId = makeNode(getArc(arcId).getEnd());
              }
              getArc(arcId).setUpNodeId(upNodeId);
            } else {
              getArc(arcId).hide();
              continue;
            }
          }
          std::pair<idVertex, idVertex> arcVerts;
          if(getArc(arcId).isEmpty()) {
            arcVerts
              = std::make_pair(getNode(upNodeId).getVertexIdentifier(),
                               getNode(downNodeId).getVertexIdentifier());
          } else {
            arcVerts = std::make_pair(
              getArc(arcId).getFirstReg(), getArc(arcId).getLastReg());
          }
          auto revertArcVerts
            = std::make_pair(std::get<1>(arcVerts), std::get<0>(arcVerts));
          if(masterArcs.count(revertArcVerts)) {
            getArc(arcId).merge(masterArcs[revertArcVerts]);
            DEBUG("Merge " << printArc(arcId) << " in "
                           << printArc(masterArcs[revertArcVerts]));
            DEBUG(" using " << std::get<0>(revertArcVerts) << " "
                            << std::get<1>(revertArcVerts));
          } else {
            masterArcs[arcVerts] = arcId;
          }
        }

        if(arc.merged()) {
          ++merged;
          // std::cout << "arc merged: " << printArc(arcId) << std::endl;
          const idSuperArc target = arc.mergedIn();
          if(mapArcs.count(target) == 0) {
            mapArcs[arcId] = target;
            consolidateArc<ScalarType>(target, arcId, s);
          }
        }
      }

      // std::cout << "Merged: " << merged << " / " << totalArc << std::endl;

      if(!mapArcs.size())
        return;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for(idVertex v = 0; v < nbElmt_; ++v) {
        const idSuperArc vArc = segmentation_[v].corArc;
        if(mapArcs.count(vArc)) {
          segmentation_[v].corArc = mapArcs[vArc];
        }
      }
    }

    template <typename ScalarType>
    void Graph::arcs2nodes(const Scalars<ScalarType> *const s) {
      const idSuperArc nbArcs = getNumberOfArcs();
      const idNode nbNodes = getNumberOfNodes();

      // fix up down
      for(idSuperArc arcId = 0; arcId < nbArcs; ++arcId) {
        if(!getArc(arcId).isVisible())
          continue;

        const idNode upNodeId = getArc(arcId).getUpNodeId();
        const idNode downNodeId = getArc(arcId).getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
        if(upNodeId == nullNode || downNodeId == nullNode) {
          continue;
        }
#endif

        const idVertex upVertId = getNode(upNodeId).getVertexIdentifier();
        const idVertex downVertId = getNode(downNodeId).getVertexIdentifier();

        if(s->isLower(upVertId, downVertId)) {
          getArc(arcId).setUpNodeId(downNodeId);
          getArc(arcId).setDownNodeId(upNodeId);
        } else {
          getArc(arcId).setUpNodeId(upNodeId);
          getArc(arcId).setDownNodeId(downNodeId);
        }
      }

      // reserve good size
      std::vector<valence> upVal(nbNodes, 0), downVal(nbNodes, 0);
      // count
      for(idSuperArc arcId = 0; arcId < nbArcs; ++arcId) {
        if(!getArc(arcId).isVisible())
          continue;

        const idNode upNodeId = getArc(arcId).getUpNodeId();
        const idNode downNodeId = getArc(arcId).getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
        if(upNodeId == nullNode || downNodeId == nullNode) {
          continue;
        }
#endif

        ++downVal[upNodeId];
        ++upVal[downNodeId];
      }
      // alloc
      for(idNode nodeId = 0; nodeId < nbNodes; ++nodeId) {
        getNode(nodeId).reserveUpArc(upVal[nodeId]);
        getNode(nodeId).reserveDownArc(downVal[nodeId]);
      }

      // set the id
      for(idSuperArc arcId = 0; arcId < nbArcs; ++arcId) {
        if(!getArc(arcId).isVisible())
          continue;

        const idNode upNodeId = getArc(arcId).getUpNodeId();
        const idNode downNodeId = getArc(arcId).getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
        if(upNodeId == nullNode || downNodeId == nullNode) {
          continue;
        }
#endif

        getNode(upNodeId).addDownArc(arcId);
        getNode(downNodeId).addUpArc(arcId);
      }
    }

    template <typename ScalarType>
    void Graph::buildArcSegmentation(const Scalars<ScalarType> *const s) {
      const idVertex nbVerts = s->getSize();
      const idSuperArc nbArcs = getNumberOfArcs();
      std::vector<idVertex> arcSizes(getNumberOfArcs(), 0);
      this->printMsg("Building arc segmentation");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for
#endif
      for(idVertex v = 0; v < nbVerts; ++v) {
#ifndef TTK_ENABLE_KAMIKAZE
        if(segmentation_[v].corArc == nullSuperArc) {
          std::cout << "[Graph]: build segmentation, null arc on vertex: " << v
                    << std::endl;
          continue;
        }
#endif
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
        arcSizes[segmentation_[v].corArc]++;
      }

      // alloc
      for(idSuperArc i = 0; i < nbArcs; ++i) {
        getArc(i).segmentation().reserve(arcSizes[i]);
      }

      // sorted list
      for(idVertex v = 0; v < nbVerts; ++v) {
        const idVertex regV = s->getSortedVert(v);
        const idSuperArc corArc = segmentation_[regV].corArc;
        getArc(corArc).segmentation().emplace_back(regV);
      }
    }

    template <typename ScalarType>
    void Graph::consolidateArc(const idSuperArc mainArc,
                               const idSuperArc mergedArc,
                               const Scalars<ScalarType> *const s) {
      const idNode mainDown = getArc(mainArc).getDownNodeId();
      const idNode mergedDown = getArc(mergedArc).getDownNodeId();

#ifndef TTK_ENABLE_KAMIKAZE
      if(mainDown == nullNode || mergedDown == nullNode) {
        std::cout << "[Graph]: consolidate error, arc have a null down node"
                  << std::endl;
      }
#endif

      // std::cout << "arc before " << printArc(mainArc) << std::endl;

      const idVertex vMainDown = getNode(mainDown).getVertexIdentifier();
      const idVertex vMergedDown = getNode(mergedDown).getVertexIdentifier();
      if(s->isLower(vMainDown, vMergedDown)) {
        // getArc(mainArc).setDownNodeId(mainDown); // already the case
        getArc(mainArc).setUpNodeId(mergedDown);
      } else {
        getArc(mainArc).setDownNodeId(mergedDown);
        getArc(mainArc).setUpNodeId(mainDown);
      }
      getArc(mainArc).restore();

      // std::cout << "arc after " << printArc(mainArc) << std::endl;
    }

  } // namespace ftr

} // namespace ttk

#endif /* end of include guard: GRAPH_TEMPLATE_H */
