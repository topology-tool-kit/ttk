/// \ingroup base
/// \class ttk:FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#include <iterator>
#include <string>

#include "FTMTree_CT.h"

using namespace std;
using namespace ttk;

using namespace ftm;

FTMTree_CT::FTMTree_CT(Params *const params,
                       Triangulation *mesh,
                       Scalars *const scalars)
  : FTMTree_MT(params, mesh, scalars, TreeType::Contour),
    jt_(new FTMTree_MT(params, mesh, scalars, TreeType::Join)),
    st_(new FTMTree_MT(params, mesh, scalars, TreeType::Split)) {
}

FTMTree_CT::~FTMTree_CT() {
  if(jt_) {
    delete jt_;
    jt_ = nullptr;
  }
  if(st_) {
    delete st_;
    st_ = nullptr;
  }
}

void FTMTree_CT::build(TreeType tt) {
  DebugTimer mergeTreesTime;

  const bool bothMT = tt == TreeType::Contour || tt == TreeType::Join_Split;

  initComp();

  if(bothMT) {
    // single leaf search for both tree
    // When executed from CT, both minima and maxima are extracted
    DebugTimer precomputeTime;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
    {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
      { leafSearch(); }
    }
    printTime(precomputeTime, "[FTM] leafSearch", -1, 3);
  }

#ifdef TTK_ENABLE_OMP_PRIORITY
  {
    // Set priority
    if(st_->getNumberOfLeaves() < jt_->getNumberOfLeaves())
      st_->setPrior();
    else
      jt_->setPrior();
  }
#endif

  // JT & ST
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
  {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
    {if(tt == TreeType::Join || bothMT){
#ifdef TTK_ENABLE_OPENMP
#pragma omp task untied if(threadNumber_ > 1)
#endif
      jt_->build(tt == TreeType::Contour);
}
if(tt == TreeType::Split || bothMT) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task untied if(threadNumber_ > 1)
#endif
  st_->build(tt == TreeType::Contour);
}
}
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

printTime(mergeTreesTime, "[FTM] merge trees ", -1, 3);

// Combine

if(tt == TreeType::Contour) {

  DebugTimer combineFullTime;
  insertNodes();

  DebugTimer combineTime;
  combine();
  printTime(combineTime, "[FTM] combine trees", -1, 4);
  printTime(combineFullTime, "[FTM] combine full", -1, 3);
}

// Debug

if(debugLevel_ > 3) {
  cout << "- [FTM] final number of nodes :";
  switch(tt) {
    case TreeType::Join:
      cout << jt_->getNumberOfNodes();
      break;
    case TreeType::Split:
      cout << st_->getNumberOfNodes();
      break;
    case TreeType::Join_Split:
      cout << jt_->getNumberOfNodes() + st_->getNumberOfNodes();
      break;
    default:
      cout << getNumberOfNodes();
  }
  cout << endl;
}
}

int FTMTree_CT::combine() {
  DebugTimer stepTime;
  queue<pair<bool, idNode>> growingNodes, remainingNodes;

  const bool DEBUG = false;

  // Reserve
  mt_data_.nodes->reserve(jt_->getNumberOfNodes());
  mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs() + 2);
  mt_data_.leaves->reserve(jt_->getNumberOfLeaves() + st_->getNumberOfLeaves());

  // Add JT & ST Leaves to growingNodes

  // Add leves to growing nodes
  const auto &nbSTLeaves = st_->getNumberOfLeaves();
  if(nbSTLeaves > 1) {
    for(idNode n = 0; n < nbSTLeaves; ++n) {
      const auto &nId = st_->getLeave(n);
      growingNodes.emplace(false, nId);
    }
  } else {
    move(jt_);
    return 0;
  }

  // count how many leaves can be added, if more than one : ok!
  const auto &nbJTLeaves = jt_->getNumberOfLeaves();
  if(nbJTLeaves > 1) {
    for(idNode n = 0; n < nbJTLeaves; ++n) {
      const auto &nId = jt_->getLeave(n);
      growingNodes.emplace(true, nId);
    }
  } // else can't clone, not same up and down

  if(DEBUG) {
    cout << "growingNodes : " << growingNodes.size()
         << " in : " << stepTime.getElapsedTime() << endl;
  }

  // Warning, have a reserve here, can't make it at the begnining, need build
  // output
  mt_data_.leaves->reserve(jt_->getLeaves().size() + st_->getLeaves().size());
  mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs());
  mt_data_.nodes->reserve(jt_->getNumberOfNodes());

  if(growingNodes.empty()) {
    cout << "[FTMTree_CT::combine ] Nothing to combine" << endl;
  }

#ifdef TTK_ENABLE_FTM_TREE_DUAL_QUEUE_COMBINE
  do {
    while(!remainingNodes.empty()) {
      bool isJT;
      idNode currentNodeId;
      FTMTree_MT *xt;

      tie(isJT, currentNodeId) = remainingNodes.front();
      remainingNodes.pop();
      if(isJT) {
        // node come frome jt
        xt = jt_;
      } else {
        // node come from st
        xt = st_;
      }
      if(xt->getNode(currentNodeId)->getNumberOfUpSuperArcs() == 1) {
        growingNodes.emplace(isJT, currentNodeId);
        if(DEBUG) {
          cout << "repush in growing:" << isJT
               << "::" << xt->printNode(currentNodeId) << endl;
        }
      }
    }
#endif

    while(!growingNodes.empty()) {
      idNode currentNodeId;
      bool isJT;

      // INFO QUEUE

      tie(isJT, currentNodeId) = growingNodes.front();
      growingNodes.pop();

      FTMTree_MT *xt = (isJT) ? jt_ : st_;
      FTMTree_MT *yt = (isJT) ? st_ : jt_;

      // INFO JT / ST

      // i <- Get(Q)
      const Node *currentNode = xt->getNode(currentNodeId);

      if(DEBUG) {
        if(xt == jt_)
          cout << endl << "JT ";
        else
          cout << endl << "ST ";
        cout << "node : " << currentNode->getVertexId() << endl;
      }

      // "choose a non-root leaf that is not a split in ST" so we ignore such
      // nodes
      if(currentNode->getNumberOfUpSuperArcs() == 0) {
        if(DEBUG) {
          cout << " ignore already processed" << endl;
        }
        continue;
      }

      idNode correspondingNodeId
        = yt->getCorrespondingNodeId(currentNode->getVertexId());

      if(yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
        if(DEBUG) {
          cout << "put remain:" << isJT << "::" << xt->printNode(currentNodeId)
               << endl;
          cout << " which is in yt : " << yt->printNode(correspondingNodeId)
               << endl;
        }
#ifdef TTK_ENABLE_FTM_TREE_DUAL_QUEUE_COMBINE
        remainingNodes.emplace(isJT, currentNodeId);
#else
      growingNodes.emplace(isJT, currentNodeId);
#endif
        continue;
      }

      // NODES IN CT

      idNode node1, node2;
      SimplexId curVert = currentNode->getVertexId();
      // NODE1
      if(isCorrespondingNode(curVert)) {
        // already a node in the tree
        node1 = getCorrespondingNodeId(curVert);
      } else {
        // create a new node
        node1 = makeNode(currentNode);

        // check if leaf
        if(!currentNode->getNumberOfDownSuperArcs())
          mt_data_.leaves->emplace_back(node1);
        else if(!currentNode->getNumberOfUpSuperArcs())
          mt_data_.leaves->emplace_back(node1);
      }

      // j <- GetAdj(XT, i)
      idSuperArc curUpArc = currentNode->getUpSuperArcId(0);
      idNode parentId = xt->getSuperArc(curUpArc)->getUpNodeId();
      const Node *parentNode = xt->getNode(parentId);

      if(DEBUG) {
        cout << " parent node :" << parentNode->getVertexId() << endl;
      }

      SimplexId parVert = parentNode->getVertexId();
      // NODE2
      if(isCorrespondingNode(parVert)) {
        // already a node in the tree
        node2 = getCorrespondingNodeId(parVert);
      } else {
        // create a new node
        node2 = makeNode(parentNode);
        if(!parentNode->getNumberOfUpSuperArcs())
          mt_data_.leaves->emplace_back(node2);
      }

      // CREATE ARC

      idSuperArc processArc = currentNode->getUpSuperArcId(0);

      // create the arc in in the good direction
      // and add it to crossing if needed
      idSuperArc createdArc;
      if(scalars_->isLower(
           currentNode->getVertexId(),
           parentNode->getVertexId())) { // take care of the order
        createdArc = makeSuperArc(node1, node2);
      } else {
        createdArc = makeSuperArc(node2, node1);
      }

      // Segmentation
      if(params_->segm) {
        createCTArcSegmentation(createdArc, isJT, processArc);
      }

      if(DEBUG) {
        cout << "create arc : " << printArc(createdArc) << endl;
      }

      // DEL NODES

      // DelNode(XT, i)
      {
        if(DEBUG) {
          cout << " delete xt (" << (xt == jt_) << ") ";
          cout << "node :" << xt->printNode(currentNodeId) << endl;
        }

        xt->delNode(currentNodeId);
      }

      // DelNode(YT, i)
      {
        if(DEBUG) {
          cout << " delete yt (" << isJT << ") node :";
          cout << yt->printNode(correspondingNodeId) << endl;
        }

        yt->delNode(correspondingNodeId);
      }

      // PROCESS QUEUE

      if(parentNode->getNumberOfDownSuperArcs() == 0
         && parentNode->getNumberOfUpSuperArcs()) {
        growingNodes.emplace(isJT, parentId);

        if(DEBUG) {
          cout << "will see : " << parentNode->getVertexId() << endl;
        }
      }
    }
#ifdef TTK_ENABLE_FTM_TREE_DUAL_QUEUE_COMBINE
  } while(!remainingNodes.empty());
#endif

  if(DEBUG) {
    printTree2();
  }

  return 0;
}

void FTMTree_CT::createCTArcSegmentation(idSuperArc ctArc,
                                         const bool isJT,
                                         idSuperArc xtArc) {
  const FTMTree_MT *xt = (isJT) ? jt_ : st_;

  /*Here we prefere to create lots of small region, each arc having its own
   * segmentation with no overlap instead of having a same vertice in several
   * arc and using vert2tree to decide because we do not want to maintain
   * vert2tree information during the whole computation*/
  const list<Region> &xtRegions = xt->getSuperArc(xtArc)->getRegions();
  for(const Region &reg : xtRegions) {
    segm_it cur = reg.segmentBegin;
    segm_it end = reg.segmentEnd;
    segm_it tmpBeg = reg.segmentBegin;
    // each element inside this region
    for(; cur != end; ++cur) {
      if(isCorrespondingNull(*cur)) {
        updateCorrespondingArc(*cur, ctArc);
      } else {
        // already set, we finish a region
        if(cur != tmpBeg) {
          getSuperArc(ctArc)->concat(tmpBeg, cur);
        }
        // if several contiguous vertices are discarded
        // cur will be equals to tmpBeg and we will not create empty regions
        tmpBeg = cur + 1;
      }
    }
    // close last region
    if(cur != tmpBeg) {
      getSuperArc(ctArc)->concat(tmpBeg, cur);
    }
  }
}

void FTMTree_CT::finalizeSegmentation(void) {
  DebugTimer finSegmTime;
  const auto &nbArc = getNumberOfSuperArcs();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(idSuperArc i = 0; i < nbArc; i++) {
    getSuperArc(i)->createSegmentation(scalars_);
  }

  printTime(finSegmTime, "[FTM] post-process segm", -1, 4);
}

void FTMTree_CT::insertNodes(void) {
  vector<idNode> sortedJTNodes = jt_->sortedNodes(true);
  vector<idNode> sortedSTNodes = st_->sortedNodes(true);

  for(const idNode &t : sortedSTNodes) {

    SimplexId vertId = st_->getNode(t)->getVertexId();
    if(jt_->isCorrespondingNode(vertId)) {
      continue;
    }
    jt_->insertNode(st_->getNode(t));
  }

  for(const idNode &t : sortedJTNodes) {

    SimplexId vertId = jt_->getNode(t)->getVertexId();
    if(st_->isCorrespondingNode(vertId)) {
      continue;
    }
    st_->insertNode(jt_->getNode(t));
  }
}

int FTMTree_CT::leafSearch() {
  const auto nbScalars = scalars_->size;
  const auto chunkSize = getChunkSize();
  const auto chunkNb = getChunkCount();

  // Extrema extract and launch tasks
  for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId)
#endif
    {
      const SimplexId lowerBound = chunkId * chunkSize;
      const SimplexId upperBound = min(nbScalars, (chunkId + 1) * chunkSize);
      for(SimplexId v = lowerBound; v < upperBound; ++v) {
        const auto &neighNumb = mesh_->getVertexNeighborNumber(v);
        valence upval = 0;
        valence downval = 0;

        for(valence n = 0; n < neighNumb; ++n) {
          SimplexId neigh;
          mesh_->getVertexNeighbor(v, n, neigh);
          if(scalars_->isLower(neigh, v)) {
            ++downval;
          } else {
            ++upval;
          }
        }

        jt_->setValence(v, downval);
        st_->setValence(v, upval);

        if(!downval) {
          jt_->makeNode(v);
        }

        if(!upval) {
          st_->makeNode(v);
        }
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
  return 0;
}
