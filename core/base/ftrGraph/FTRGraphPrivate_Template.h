#pragma once

// local includes
#include "FTRGraph.h"
#include "FTRGraph_Template.h"
#include "FTRTasks.h"

// c++ incldues
#include <unordered_map>

// trick to print a full line in an atomic operation, avoid mixed up redulsts in
// parallel
#ifndef NDEBUG
// #define PRINT(msg) {std::stringstream s; s << msg << std::endl; std::cout <<
// s.str();}
#define PRINT(msg)
#else
#define PRINT(msg)
#endif

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::growthFromSeed(
  const idVertex seed, Propagation *localProp, idSuperArc currentArc) {

  TTK_FORCE_USE(seed);
  PRINT(seed << " :: " << localProp->getNbArcs());
  if(!localProp->getNbArcs()) {
    PRINT(seed << " Stop direct");
#ifdef TTK_ENABLE_FTR_TASK_STATS
    {
      // This only when tasks grows exclusively from min
      // This propagation is dying here
      idVertex curProp;
      float curTime;
#pragma omp atomic capture seq_cst
      {
        curProp = nbProp_;
        --nbProp_;
      }
#pragma omp critical(stats)
      { curTime = sweepStart_.getElapsedTime(); }
      propTimes_[curProp - 1] = curTime;
    }
#endif
    return;
  }

  // topology
  bool isJoinLast = false;
  bool isJoin = false, isSplit = false;

  // containers
  // vitsit
  Star star;
  Comp comp;

  while(!isJoin && !isSplit && !localProp->empty()) {
    localProp->nextVertex();
    const idVertex curVert = localProp->getCurVertex();
    idSuperArc mergeIn = nullSuperArc;

    if(localProp->getNbArcs() == 0) {
#ifdef TTK_ENABLE_FTR_TASK_STATS
      {
        // This only when tasks grows exclusively from min
        // This propagation is dying here
        idVertex curProp;
        float curTime;
#pragma omp atomic capture seq_cst
        {
          curProp = nbProp_;
          --nbProp_;
        }
#pragma omp critical(stats)
        { curTime = sweepStart_.getElapsedTime(); }
        propTimes_[curProp - 1] = curTime;
      }
#endif
      return;
    }

    PRINT("<" << curVert << " " << localProp->goUp() << " a "
              << localProp->getNbArcs());

#ifdef TTK_ENABLE_FTR_VERT_STATS
    graph_.incTouch(curVert);
    graph_.setNbArcActive(curVert, localProp->getNbArcs());
#endif

    visitStar(localProp, star);

    if(propagations_.hasVisitedOpposite(
         curVert, localProp) /*&& !graph_.isNode(curVert)*/) {
      bool ignoreVert = false;
      for(auto edge : star.lower) {
        const idSuperArc tmpLowArc = dynGraph(localProp).getSubtreeArc(edge);
        if(tmpLowArc != nullSuperArc && !graph_.getArc(tmpLowArc).isVisible()
           && graph_.getArc(tmpLowArc).getPropagation() == localProp) {
          PRINT("-" << curVert << " " << graph_.printArc(tmpLowArc));
#ifdef TTK_ENABLE_FTR_VERT_STATS
          graph_.incAvoid();
#endif
          ignoreVert = true;
          break;
        } else if(tmpLowArc != nullSuperArc) {
          PRINT("!" << graph_.printArc(tmpLowArc));
        }
      }
      if(ignoreVert)
        continue;
    }

#ifndef TTK_DISABLE_FTR_LAZY
    if(valences_.lower[curVert] < 2 && valences_.upper[curVert] < 2) {

      // not a local min (for local min, currentArc is already set)
      if(star.lower.size()) {
        // not a min nor a saddle: 1 CC below
        currentArc = dynGraph(localProp).getSubtreeArc(star.lower[0]);
        if(currentArc == nullSuperArc) {
          PRINT("n-" << curVert);
          continue;
        }
        if(valences_.upper[curVert] && valences_.lower[curVert]) {
          // not saddle neither extrema
          graph_.getArc(currentArc).visit(curVert);
        }
      }

      mergeIn = visit(localProp, currentArc);

      lazyUpdatePreimage(localProp, currentArc);

      // ensure we will always recover this arc from the upper neighbors
      for(const idEdge dgNode : star.upper) {
        dynGraph(localProp).setCorArc(dgNode, currentArc);
      }
    } else {
      // locally apply the lazy one the current growing arc
      for(const idEdge e : star.lower) {
        const idSuperArc a = dynGraph(localProp).getNode(e)->findRootArc();
        if(a != nullSuperArc && graph_.getArc(a).isVisible()
           && graph_.getArc(a).getPropagation()->getId()
                == localProp->getId()) {
          lazyApply(localProp, a);
        }
      }
#else
    {
#endif
      comp.lower = lowerComps(star.lower, localProp);

      if(comp.lower.size() > 1) {
        isJoin = true;
        isJoinLast = checkLast(localProp, star.lower);
        break;
      } else {
        if(comp.lower.size()) {
          currentArc = (*comp.lower.begin())->getCorArc();
          if(currentArc == nullSuperArc) {
            PRINT("n--" << curVert);
            continue;
          }
        }
        mergeIn = visit(localProp, currentArc);
      }
      updatePreimage(localProp, currentArc);
      comp.upper = upperComps(star.upper, localProp);
      if(comp.upper.size() > 1) {
        isSplit = true;
      }

      if(!isJoin && !isSplit && comp.upper.size()) {
        // this arc is not empty (not saddle not max)
        graph_.getArc(currentArc).visit(curVert);
      }
    }

    // merge current arc
    if(mergeIn != nullSuperArc && graph_.getArc(currentArc).isVisible()) {
      graph_.getArc(currentArc).merge(mergeIn);
      localProp->lessArc();
      PRINT(">" << graph_.printArc(currentArc) << " "
                << graph_.printArc(mergeIn));
      // Can't stop here
      // try interleaved launch sequential on dragon.vtu
      // if (localProp->getNbArcs() == 0) return;
    }

    // stop on leaves
    if(!star.upper.size()) {
      // We have reached a local extrema (max from this propagation)
      const idNode leafNode = graph_.makeNode(curVert);
      graph_.closeArc(currentArc, leafNode);
      PRINT("^" << graph_.printArc(currentArc));
      if(graph_.getArc(currentArc).isVisible()) {
        // do not decrease on merged arcs
        localProp->lessArc();
      }
      if(localProp->getNbArcs() == 0) {
        // no more active arcs
        PRINT(curVert << " stop");
#ifdef TTK_ENABLE_FTR_TASK_STATS
        {
          // This only when tasks grows exclusively from min
          // This propagation is dying here
          idVertex curProp;
          float curTime;
#pragma omp atomic capture seq_cst
          {
            curProp = nbProp_;
            --nbProp_;
          }
#pragma omp critical(stats)
          { curTime = sweepStart_.getElapsedTime(); }
          propTimes_[curProp - 1] = curTime;
        }
#endif
        return;
      }
    }

    // add upper star for futur visit
    localGrowth(localProp, star.upper);
  } // end propagation while

  // get the corresponging critical point on which
  // the propagation has stopped (join, split, max)
  const idVertex upVert = localProp->getCurVertex();

  if(localProp->getNbArcs() == 0 || localProp->empty()) {
    PRINT(upVert << " stop " << localProp->getNbArcs());
#ifdef TTK_ENABLE_FTR_TASK_STATS
    {
      // This only when tasks grows exclusively from min
      // This propagation is dying here
      idVertex curProp;
      float curTime;
#pragma omp atomic capture seq_cst
      {
        curProp = nbProp_;
        --nbProp_;
      }
#pragma omp critical(stats)
      { curTime = sweepStart_.getElapsedTime(); }
      propTimes_[curProp - 1] = curTime;
    }
#endif
    return;
  }

  PRINT(upVert << " active " << localProp->getNbArcs());

  // reached node id and wether it has been created by this task or already
  // existed
  idNode saddleNode;
  idSuperArc joinParentArc{};
  bool hideFromHere = false; // if true, new arc are hidden to stop propagation.

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
  {
    bool alreadyNode = graph_.isNode(upVert);
    hideFromHere = alreadyNode;
    if(isJoin) {
      PRINT(upVert << " join " << isJoinLast << " h " << hideFromHere);
      if(isJoinLast) {
        graph_.makeNode(upVert);
      }
    }
    if(isSplit) {
      PRINT(upVert << " split " << hideFromHere);
      graph_.makeNode(upVert);
    }
  }

  // arc management
  if(isJoinLast) {
    // ensure we have the good values here, even if other tasks were doing
    // stuff
    {
      // here to solve a 1 over thousands execution bug in parallel
      // TODO Still required ??
      visitStar(localProp, star);
      comp.lower = lowerComps(star.lower, localProp);
    }
    saddleNode = graph_.getNodeId(upVert);
    idSuperArc visibleMerged = mergeAtSaddle(saddleNode, localProp, comp.lower);
    localProp->lessArc(visibleMerged - 1);

    localGrowth(localProp, star.upper);

    joinParentArc = graph_.openArc(saddleNode, localProp);
    visit(localProp, joinParentArc);
    updatePreimage(localProp, joinParentArc);
    comp.upper = upperComps(star.upper, localProp);

    // do not propagate
    if(hideFromHere) {
      if(graph_.getArc(joinParentArc).hide()) {
        localProp->lessArc();
      }
    }

    // split detection required after the merge
    if(comp.upper.size() > 1) {
      // this node is both join and split
      isSplit = true;
      // will be replaced be new arcs of the split
      if(graph_.getArc(joinParentArc).hide()) {
        localProp->lessArc();
      }
    }

    PRINT("+" << graph_.printArc(joinParentArc));
  } else if(!isJoin) {
    // only one arc coming here
    saddleNode = graph_.getNodeId(upVert);
    graph_.closeArc(currentArc, saddleNode);
    if(localProp->getNbArcs() == 0) {
      // no more active arcs
      PRINT(upVert << " stop");
#ifdef TTK_ENABLE_FTR_TASK_STATS
      {
        // This only when tasks grows exclusively from min
        // This propagation is dying here
        idVertex curProp;
        float curTime;
#pragma omp atomic capture seq_cst
        {
          curProp = nbProp_;
          --nbProp_;
        }
#pragma omp critical(stats)
        { curTime = sweepStart_.getElapsedTime(); }
        propTimes_[curProp - 1] = curTime;
      }
#endif
      return;
    }
    if(graph_.getArc(currentArc).isVisible()) {
      // if not visible, already decremented the counter
      localProp->lessArc();
    }
    PRINT("|" << graph_.printArc(currentArc));
  }

  if(isSplit) {
    splitAtSaddle(localProp, comp.upper, hideFromHere);
    if(!hideFromHere) {
      localProp->moreArc(comp.upper.size());
    }
  }

  // starting from the saddle
  if(isSplit && (!isJoin || isJoinLast)) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task OPTIONAL_PRIORITY(PriorityLevel::Low)
#endif
    growthFromSeed(upVert, localProp);

  } else if(isJoinLast) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp task OPTIONAL_PRIORITY(PriorityLevel::Average)
#endif
    growthFromSeed(upVert, localProp, joinParentArc);
  }
#ifdef TTK_ENABLE_FTR_TASK_STATS
  else {
    // This only when tasks grows exclusively from min
    // This propagation is dying here
    idVertex curProp;
    float curTime;
#pragma omp atomic capture seq_cst
    {
      curProp = nbProp_;
      --nbProp_;
    }
#pragma omp critical(stats)
    { curTime = sweepStart_.getElapsedTime(); }
    propTimes_[curProp - 1] = curTime;
  }
#endif
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::growthSequential(
  const idVertex begin, const idVertex stop) {
  Star star;
  Comp comp;

  const bool fromMin = begin < stop;
  Propagation *localProp
    = newPropagation(scalars_.getSortedVert(begin), fromMin);
  const idVertex incr = fromMin ? 1 : -1;
  for(idVertex idv = begin; idv < stop; idv = idv + incr) {
    const idVertex curVert = scalars_.getSortedVert(idv);
    localProp->setCurvert(curVert); // never use Fibo Heap here

    idSuperArc currentArc = nullSuperArc;
    visitStar(localProp, star);

    if(valences_.lower[curVert] == 0) { // min
      const idNode minNode = graph_.makeNode(curVert);
      currentArc = graph_.openArc(minNode, localProp);
    }

#ifndef TTK_DISABLE_FTR_LAZY
    if(valences_.lower[curVert] < 2 && valences_.upper[curVert] < 2) {

      // simple reeb regular, lazyness
      if(star.lower.size()) {
        // not a min nor a saddle: 1 CC below (need findSubtree)
        currentArc = dynGraph(localProp).getSubtreeArc(star.lower[0]);
        if(valences_.upper[curVert] && valences_.lower[curVert]) {
          // not saddle neither extrema
          graph_.getArc(currentArc).visit(curVert);
        }
      }
      if(star.upper.size()) {
        for(const idEdge dgNode : star.upper) {
          dynGraph(localProp).setCorArc(dgNode, currentArc);
        }
      } else { // max
        const idNode maxNode = graph_.makeNode(curVert);
        graph_.closeArc(currentArc, maxNode);
      }

      visit(localProp, currentArc);

      lazyUpdatePreimage(localProp, currentArc);
    } else {

      // locally aply the lazy one the current growing arc
      for(const idEdge e : star.lower) {
        const idSuperArc a = dynGraph(localProp).getNode(e)->findRootArc();
        if(!lazy_.isEmpty(a)) {
          // process lazy
          // sort both list
          // for each:
          // if del < add: delete in real RG
          // if add < del: add in real RG
          // else: drop both, computation avoided
          lazyApply(localProp, a);
        }
      }
#else
    {
#endif
      bool isJoin = false;
      comp.lower = lowerComps(star.lower, localProp);
      if(comp.lower.size() == 1) { // regular
        currentArc = (*comp.lower.begin())->getCorArc();
      } else if(comp.lower.size() > 1) { // join saddle
        const idNode sadNode = graph_.makeNode(curVert);
        currentArc = graph_.openArc(sadNode, localProp);
        mergeAtSaddle(sadNode, comp.lower);
        isJoin = true;
      }

      graph_.visit(curVert, currentArc);
      propagations_.visit(curVert, localProp);
      updatePreimage(localProp, currentArc);

      comp.upper = upperComps(star.upper, localProp);
      if(!comp.upper.size()) { // max
        const idNode maxNode = graph_.makeNode(curVert);
        graph_.closeArc(currentArc, maxNode);
      } else if(comp.upper.size() < 2) {
        if(!isJoin) {
          // this arc is not empty
          graph_.getArc(currentArc).visit(curVert);
        }
      } else {
        if(isJoin) {
          graph_.getArc(currentArc).hide();
        } else { // split
          const idNode splitNode = graph_.makeNode(curVert);
          graph_.closeArc(currentArc, splitNode);
        }
        splitAtSaddle(localProp, comp.upper);
      }
    }
  } // end for each vertex
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::visitStar(
  const Propagation *const localProp, Star &star) const {

  star.lower.clear();
  star.upper.clear();

  const idEdge nbAdjEdges
    = mesh_.getVertexEdgeNumber(localProp->getCurVertex());
  star.lower.reserve(nbAdjEdges);
  star.upper.reserve(nbAdjEdges);

  for(idEdge e = 0; e < nbAdjEdges; ++e) {
    idEdge edgeId;
    mesh_.getVertexEdge(localProp->getCurVertex(), e, edgeId);
    idVertex edgeLowerVert, edgeUpperVert;
    std::tie(edgeLowerVert, edgeUpperVert)
      = mesh_.getOrderedEdge(edgeId, localProp->goUp());
    if(edgeLowerVert == localProp->getCurVertex()) {
      star.upper.emplace_back(edgeId);
    } else {
      star.lower.emplace_back(edgeId);
    }
  }
}

template <typename ScalarType, typename triangulationType>
std::set<ttk::ftr::DynGraphNode<ttk::ftr::idVertex> *>
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::lowerComps(
    const std::vector<idEdge> &finishingEdges,
    const Propagation *const localProp) {
  return dynGraph(localProp).findRoot(finishingEdges);
}

template <typename ScalarType, typename triangulationType>
std::set<ttk::ftr::DynGraphNode<ttk::ftr::idVertex> *>
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::upperComps(
    const std::vector<idEdge> &startingEdges,
    const Propagation *const localProp) {
  return dynGraph(localProp).findRoot(startingEdges);
}

template <typename ScalarType, typename triangulationType>
bool ttk::ftr::FTRGraph<ScalarType, triangulationType>::checkStop(
  const std::vector<DynGraphNode<idVertex> *> &compVect) {
  for(const auto *dgNode : compVect) {
    const idSuperArc arc = dgNode->getCorArc();
    if(arc != nullSuperArc && !graph_.getArc(arc).isVisible()) {
      return true;
    }
  }
  return false;
}

template <typename ScalarType, typename triangulationType>
std::pair<ttk::ftr::valence, ttk::ftr::valence>
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::getLinkNbCC(
    const idVertex curVert,
    LocalForests &localForests,
    const VertCompFN &comp) {
  // traduce edge id in a local id for the forests
  std::unordered_map<idEdge, std::size_t> mapNeighDown, mapNeighUp;
  std::size_t nextId = 0;

  localForests.up.reset();
  localForests.down.reset();
  const idVertex oldUpCC = localForests.up.getNbCC();
  const idVertex oldDownCC = localForests.down.getNbCC();
  const idCell nbTri = mesh_.getVertexTriangleNumber(curVert);
  for(idCell t = 0; t < nbTri; ++t) {
    idCell curTri;
    mesh_.getVertexTriangle(curVert, t, curTri);
    // edges found, max 2 per triangles.
    std::size_t downSide[2], upSide[2];
    unsigned char curDownSide = 0, curUpSide = 0;
    for(idVertex v = 0; v < 3; ++v) {
      idEdge triEdge;
      mesh_.getTriangleEdge(curTri, v, triEdge);
      idVertex v0, v1;
      mesh_.getEdgeVertex(triEdge, 0, v0);
      mesh_.getEdgeVertex(triEdge, 1, v1);

      // only consider the edges in the star
      if(v0 == curVert) {
        if(comp(curVert, v1)) {
          if(mapNeighUp.count(triEdge)) {
            upSide[curUpSide++] = mapNeighUp[triEdge];
          } else {
            upSide[curUpSide++] = nextId;
            mapNeighUp[triEdge] = nextId++;
          }
        } else {
          if(mapNeighDown.count(triEdge)) {
            downSide[curDownSide++] = mapNeighDown[triEdge];
          } else {
            downSide[curDownSide++] = nextId;
            mapNeighDown[triEdge] = nextId++;
          }
        }
      } else if(v1 == curVert) {
        if(comp(curVert, v0)) {
          if(mapNeighUp.count(triEdge)) {
            upSide[curUpSide++] = mapNeighUp[triEdge];
          } else {
            upSide[curUpSide++] = nextId;
            mapNeighUp[triEdge] = nextId++;
          }
        } else {
          if(mapNeighDown.count(triEdge)) {
            downSide[curDownSide++] = mapNeighDown[triEdge];
          } else {
            downSide[curDownSide++] = nextId;
            mapNeighDown[triEdge] = nextId++;
          }
        }
      }
    }

    // if both edges of the triangle were up or down, we add a link btwn
    // them This is how the number of component is reduced
    if(curDownSide == 2) {
      localForests.down.insertEdge(downSide[0], downSide[1], 0, nullSuperArc);
    } else if(curUpSide == 2) {
      localForests.up.insertEdge(upSide[0], upSide[1], 0, nullSuperArc);
    }
  }

  const valence down
    = mapNeighDown.size() - (oldDownCC - localForests.down.getNbCC());
  const valence up = mapNeighUp.size() - (oldUpCC - localForests.up.getNbCC());
  return {down, up};
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updatePreimage(
  const Propagation *const localProp, const idSuperArc curArc) {
  const idCell nbAdjTriangles
    = mesh_.getVertexTriangleNumber(localProp->getCurVertex());

  orderedTriangle oTriangle;

  // TODO  SORT TRANGLES HERE ?

  for(idCell t = 0; t < nbAdjTriangles; ++t) {
    // Classify current cell
    idCell curTriangleid;
    mesh_.getVertexTriangle(localProp->getCurVertex(), t, curTriangleid);

    mesh_.getOrderedTriangle(curTriangleid, localProp->goUp(), oTriangle);
    vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localProp);

    // Update DynGraph
    // We can have an end pos on an unvisited triangle
    // in case of saddle points
    switch(curVertPos) {
      case vertPosInTriangle::Start:
        updatePreimageStartCell(oTriangle, localProp, curArc);
        break;
      case vertPosInTriangle::Middle:
        updatePreimageMiddleCell(oTriangle, localProp, curArc);
        break;
      case vertPosInTriangle::End:
        updatePreimageEndCell(oTriangle, localProp, curArc);
        break;
      default:
        std::cout << "[FTR]: update preimage error, unknown vertPos type"
                  << std::endl;
        break;
    }
  }
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updatePreimageStartCell(
  const orderedTriangle &oTriangle,
  const Propagation *const localProp,
  const idSuperArc curArc) {
  const orderedEdge e0
    = mesh_.getOrderedEdge(std::get<0>(oTriangle), localProp->goUp());
  const orderedEdge e1
    = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
  const idVertex w = getWeight(e0, e1, localProp);

  // this order for history
  dynGraph(localProp).insertEdge(
    std::get<1>(oTriangle), std::get<0>(oTriangle), w, curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::
  updatePreimageMiddleCell(const orderedTriangle &oTriangle,
                           const Propagation *const localProp,
                           const idSuperArc curArc) {
  // Check if exist ?
  // If not, the triangle will be visited again once a merge have occured.
  // So we do not add the edge now
  dynGraph(localProp).removeEdge(
    std::get<0>(oTriangle), std::get<1>(oTriangle));

  // keep history inside the dyngraph structure
  dynGraph(localProp).setCorArc(std::get<0>(oTriangle), curArc);
  // dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);

  const orderedEdge e1
    = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
  const orderedEdge e2
    = mesh_.getOrderedEdge(std::get<2>(oTriangle), localProp->goUp());
  const idVertex w = getWeight(e1, e2, localProp);

  // this order for history
  dynGraph(localProp).insertEdge(
    std::get<1>(oTriangle), std::get<2>(oTriangle), w, curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updatePreimageEndCell(
  const orderedTriangle &oTriangle,
  const Propagation *const localProp,
  const idSuperArc curArc) {
  dynGraph(localProp).removeEdge(
    std::get<1>(oTriangle), std::get<2>(oTriangle));
  dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);
  // dynGraph(localProp).setCorArc(std::get<2>(oTriangle), curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateDynGraphCurArc(
  const idVertex seed,
  const idEdge neigEdge,
  const idSuperArc curArc,
  const Propagation *const localProp) {
  idVertex v0;
  idVertex v1;
  mesh_.getEdgeVertex(neigEdge, 0, v0);
  mesh_.getEdgeVertex(neigEdge, 1, v1);

  const idVertex other = (v0 == seed) ? v1 : v0;

  if(localProp->compare(seed, other)) {
    dynGraph(localProp).setSubtreeArc(neigEdge, curArc);
  }
}

#ifndef TTK_DISABLE_FTR_LAZY
template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::lazyUpdatePreimage(
  Propagation *const localProp, const idSuperArc curArc) {
  const idCell nbAdjTriangles
    = mesh_.getVertexTriangleNumber(localProp->getCurVertex());

  orderedTriangle oTriangle;

  for(idCell t = 0; t < nbAdjTriangles; ++t) {
    // Classify current cell
    idCell curTriangleid;
    mesh_.getVertexTriangle(localProp->getCurVertex(), t, curTriangleid);

    mesh_.getOrderedTriangle(curTriangleid, localProp->goUp(), oTriangle);
    vertPosInTriangle curVertPos = getVertPosInTriangle(oTriangle, localProp);

    // Update DynGraph
    // We can have an end pos on an unvisited triangle
    // in case of saddle points
    switch(curVertPos) {
      case vertPosInTriangle::Start:
        updateLazyStart(oTriangle, localProp, curArc);
        break;
      case vertPosInTriangle::Middle:
        updateLazyMiddle(oTriangle, localProp, curArc);
        break;
      case vertPosInTriangle::End:
        updateLazyEnd(oTriangle, localProp, curArc);
        break;
      default:
        std::cout << "[FTR]: lazy update preimage error, unknown vertPos type"
                  << std::endl;
        break;
    }
  }
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateLazyStart(
  const orderedTriangle &oTriangle,
  Propagation *const ttkNotUsed(localProp),
  const idSuperArc curArc) {
  lazy_.addEmplace(std::get<0>(oTriangle), std::get<1>(oTriangle), curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateLazyMiddle(
  const orderedTriangle &oTriangle,
  Propagation *const localProp,
  const idSuperArc curArc) {
  lazy_.delEmplace(std::get<0>(oTriangle), std::get<1>(oTriangle), curArc);
  dynGraph(localProp).removeEdge(
    std::get<0>(oTriangle), std::get<1>(oTriangle));
  dynGraph(localProp).setCorArc(std::get<0>(oTriangle), curArc);
  dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);
  lazy_.addEmplace(std::get<1>(oTriangle), std::get<2>(oTriangle), curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateLazyEnd(
  const orderedTriangle &oTriangle,
  Propagation *const localProp,
  const idSuperArc curArc) {
  lazy_.delEmplace(std::get<1>(oTriangle), std::get<2>(oTriangle), curArc);
  dynGraph(localProp).removeEdge(
    std::get<1>(oTriangle), std::get<2>(oTriangle));
  dynGraph(localProp).setCorArc(std::get<1>(oTriangle), curArc);
  dynGraph(localProp).setCorArc(std::get<2>(oTriangle), curArc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateLazyAdd(
  const Propagation *const localProp,
  const linkEdge edge,
  const idSuperArc arc) {
  const orderedEdge e0
    = mesh_.getOrderedEdge(std::get<0>(edge), localProp->goUp());
  const orderedEdge e1
    = mesh_.getOrderedEdge(std::get<1>(edge), localProp->goUp());
  const idVertex w = getWeight(e0, e1, localProp);
  dynGraph(localProp).insertEdge(std::get<1>(edge), std::get<0>(edge), w, arc);
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateLazyDel(
  const Propagation *const localProp,
  const linkEdge edge,
  const idSuperArc ttkNotUsed(arc)) {
  dynGraph(localProp).removeEdge(std::get<0>(edge), std::get<1>(edge));
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::lazyApply(
  Propagation *const localProp, const idSuperArc a) {
  auto add = lazy_.addGetNext(a);
  while(add != nullLink) {
    updateLazyAdd(localProp, add, a);
    add = lazy_.addGetNext(a);
  }
}

#endif

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::updateDynGraphCurArc(
  const idVertex seed,
  const idSuperArc curArc,
  const Propagation *const localProp) {
  const idVertex nbEdgesNeigh = mesh_.getVertexEdgeNumber(seed);
  for(idVertex nid = 0; nid < nbEdgesNeigh; ++nid) {
    idEdge edgeId;
    mesh_.getVertexEdge(seed, nid, edgeId);

    updateDynGraphCurArc(seed, edgeId, curArc, localProp);
  }
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::localGrowth(
  Propagation *const localProp, const std::vector<idEdge> &upperEdges) {
  const idVertex curVert = localProp->getCurVertex();
  for(const idEdge e : upperEdges) {
    idVertex v0, v1;
    mesh_.getEdgeVertex(e, 0, v0);
    mesh_.getEdgeVertex(e, 1, v1);
    const idVertex other = (v0 == curVert) ? v1 : v0;
    if(localProp->compare(curVert, other)) {
      if(!propagations_.willVisit(other, localProp)) {
        localProp->addNewVertex(other);
        propagations_.toVisit(other, localProp);
      }
    }
  }
}

template <typename ScalarType, typename triangulationType>
bool ttk::ftr::FTRGraph<ScalarType, triangulationType>::checkLast(
  Propagation *const localProp, const std::vector<idEdge> &starVect) {
  const idVertex curSaddle = localProp->getCurVertex();
  AtomicUF *curId = localProp->getId();
  valence decr = 0;

  // NOTE:
  // Using propagation id allows to decrement by the number of time this
  // propagation has reached the saddle, even if the propagation take care
  // of several of these arcs (after a Hole-split).
  for(idEdge edgeId : starVect) {
    const idSuperArc edgeArc = dynGraph(localProp).getSubtreeArc(edgeId);
    if(edgeArc == nullSuperArc) {
      continue;
    }
    AtomicUF *tmpId = graph_.getArc(edgeArc).getPropagation()->getId();
    if(tmpId == curId) {
      graph_.getArc(edgeArc).setEnd(curSaddle);
      ++decr;
    }
  }

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#endif // __GNUC__

  valence oldVal = 0;
  if(localProp->goUp()) {
    // for gcc 4.8 and old openMP
    valence *const vd = &graph_.valDown_[curSaddle];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
    {
      oldVal = *vd;
      *vd -= decr;
    }

  } else {
    valence *const vu = &graph_.valUp_[curSaddle];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
    {
      oldVal = *vu;
      *vu -= decr;
    }
  }

  if(oldVal == -1) {
    // First task to touch this saddle, compute the valence
    idVertex totalVal = starVect.size();
    valence newVal = 0;
    if(localProp->goUp()) {
      valence *const vd = &graph_.valDown_[curSaddle];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      {
        newVal = *vd;
        *vd += (totalVal + 1);
      }

    } else {
      valence *const vu = &graph_.valUp_[curSaddle];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      {
        newVal = *vu;
        *vu += (totalVal + 1);
      }
    }
    oldVal = decr + newVal + (totalVal + 1);
  }

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif // __GNUC__

  return oldVal == decr;
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idSuperArc
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::mergeAtSaddle(
    const idNode saddleId,
    Propagation *localProp,
    const std::set<DynGraphNode<idVertex> *> &compVect) {

#ifndef TTK_ENABLE_KAMIKAZE
  if(compVect.size() < 2) {
    std::cerr << "[FTR]: merge at saddle with only one lower CC" << std::endl;
  }
#endif

  idSuperArc visibleClosed = 0;
  for(auto *dgNode : compVect) {
    // read in the history (lower comp already contains roots)
    const idSuperArc endingArc = dgNode->getCorArc();
    graph_.closeArc(endingArc, saddleId);
    PRINT("/" << graph_.printArc(endingArc));
    if(graph_.getArc(endingArc).isVisible()) {
      ++visibleClosed;
    }
    Propagation *arcProp = graph_.getArc(endingArc).getPropagation();
    localProp->merge(*arcProp);
  }
  return visibleClosed;
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idSuperArc
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::mergeAtSaddle(
    const idNode saddleId, const std::set<DynGraphNode<idVertex> *> &compVect) {
  // version for the sequential arc growth, do not merge the propagations

#ifndef TTK_ENABLE_KAMIKAZE
  if(compVect.size() < 2) {
    std::cerr << "[FTR]: merge at saddle with only one lower CC" << std::endl;
  }
#endif

  idSuperArc visibleClosed = 0;
  for(auto *dgNode : compVect) {
    const idSuperArc endingArc = dgNode->getCorArc();
    graph_.closeArc(endingArc, saddleId);
    PRINT("/" << graph_.printArc(endingArc));
    if(graph_.getArc(endingArc).isVisible()) {
      ++visibleClosed;
    }
  }
  return visibleClosed;
}

template <typename ScalarType, typename triangulationType>
void ttk::ftr::FTRGraph<ScalarType, triangulationType>::splitAtSaddle(
  Propagation *const localProp,
  const std::set<DynGraphNode<idVertex> *> &compVect,
  const bool hidden) {
  const idVertex curVert = localProp->getCurVertex();
  const idNode curNode = graph_.getNodeId(curVert);

  for(auto *dgNode : compVect) {
    const idSuperArc newArc = graph_.openArc(curNode, localProp);
    dgNode->setRootArc(newArc);
    visit(localProp, newArc);

    if(hidden)
      graph_.getArc(newArc).hide();

    PRINT("v" << graph_.printArc(newArc));
  }
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idSuperArc ttk::ftr::FTRGraph<ScalarType, triangulationType>::visit(
  Propagation *const localProp, const idSuperArc curArc) {
  const idVertex curVert = localProp->getCurVertex();
  Visit opposite;
  idSuperArc retArc;
  // #pragma omp critical
  {
    propagations_.visit(curVert, localProp);
    opposite = propagations_.visitOpposite(curVert, localProp);
    bool done;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
    done = opposite.done;
    if(!done) {
      graph_.visit(curVert, curArc);
      retArc = nullSuperArc;
    } else {
      retArc = graph_.getArcId(curVert);
    }
  }
  return retArc;
}

/// Tools

template <typename ScalarType, typename triangulationType>
ttk::ftr::Propagation *
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::newPropagation(
    const idVertex leaf, const bool fromMin) {
  VertCompFN comp;
  if(fromMin)
    comp = [&](idVertex a, idVertex b) { return scalars_.isHigher(a, b); };
  else
    comp = [&](idVertex a, idVertex b) { return scalars_.isLower(a, b); };
  return propagations_.newPropagation(leaf, comp, fromMin);
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idVertex ttk::ftr::FTRGraph<ScalarType, triangulationType>::getWeight(
  const orderedEdge &e0,
  const orderedEdge &e1,
  const Propagation *const localProp) {
  const idVertex end0 = std::get<1>(e0);
  const idVertex end1 = std::get<1>(e1);

  if(localProp->compare(end0, end1)) {
    if(localProp->goDown()) {
      return -scalars_.getMirror(end0);
    }
    return scalars_.getMirror(end0);
  }

  if(localProp->goDown()) {
    return -scalars_.getMirror(end1);
  }
  return scalars_.getMirror(end1);
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::vertPosInTriangle
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::getVertPosInTriangle(
    const orderedTriangle &oTriangle,
    const Propagation *const localProp) const {
  orderedEdge firstEdge
    = mesh_.getOrderedEdge(std::get<0>(oTriangle), localProp->goUp());
  if(std::get<0>(firstEdge) == localProp->getCurVertex()) {
    return vertPosInTriangle::Start;
  } else if(std::get<1>(firstEdge) == localProp->getCurVertex()) {
    return vertPosInTriangle::Middle;
  } else {
    return vertPosInTriangle::End;
  }
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idVertex
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::getEndVertexInTriangle(
    const orderedTriangle &oTriangle,
    const Propagation *const localProp) const {
  const orderedEdge &higherEdge
    = mesh_.getOrderedEdge(std::get<1>(oTriangle), localProp->goUp());
  return std::get<1>(higherEdge);
}

template <typename ScalarType, typename triangulationType>
ttk::ftr::idEdge
  ttk::ftr::FTRGraph<ScalarType, triangulationType>::getEdgeFromOTri(
    const orderedTriangle oTri, const idVertex v0, const idVertex v1) {
  idVertex edge0Vert, edge1Vert;

  mesh_.getEdgeVertex(std::get<0>(oTri), 0, edge0Vert);
  mesh_.getEdgeVertex(std::get<0>(oTri), 1, edge1Vert);
  if((edge0Vert == v0 && edge1Vert == v1)
     || (edge0Vert == v1 && edge1Vert == v0)) {
    return std::get<0>(oTri);
  }

  mesh_.getEdgeVertex(std::get<1>(oTri), 0, edge0Vert);
  mesh_.getEdgeVertex(std::get<1>(oTri), 1, edge1Vert);
  if((edge0Vert == v0 && edge1Vert == v1)
     || (edge0Vert == v1 && edge1Vert == v0)) {
    return std::get<1>(oTri);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  mesh_.getEdgeVertex(std::get<2>(oTri), 0, edge0Vert);
  mesh_.getEdgeVertex(std::get<2>(oTri), 1, edge1Vert);
  if((edge0Vert == v0 && edge1Vert == v1)
     || (edge0Vert == v1 && edge1Vert == v0)) {
    return std::get<2>(oTri);
  }

  std::cout << "[FTR]: edge not found in triangle " << v0 << " " << v1
            << std::endl;
  return nullEdge;
#else
  return std::get<2>(oTri);
#endif
}
