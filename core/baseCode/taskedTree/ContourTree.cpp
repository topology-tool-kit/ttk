/*
 * file: ContourTree.cpp
 * description: ContourTree processing package.
 * author: Gueunet Charles
 * date: Juin 2015
 */

#include <iterator>
#include <string>

#include "ContourTree.h"

// ------------
// CONSTRUCTOR
// ------------
// {

ContourTree::ContourTree(Params *const params, Triangulation *mesh, Scalars *const scalars)
    : MergeTree(params, mesh, scalars, TreeType::Contour),
      jt_(new MergeTree(params, mesh, scalars, TreeType::Join)),
      st_(new MergeTree(params, mesh, scalars, TreeType::Split))
{
}

ContourTree::~ContourTree()
{
   if (jt_) {
      delete jt_;
      jt_ = nullptr;
   }
   if (st_) {
      delete st_;
      st_ = nullptr;
   }
}

// }
// -------
// PROCESS
// -------
// {
int ContourTree::vertexPrecomputation()
{
   const auto  nbScalars = scalars_->size;
   const auto  chunkSize = getChunkSize();
   const auto  chunkNb   = getChunkCount();

   // --------------------------------
   // Extrema extract and launch tasks
   // --------------------------------
   for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#pragma omp task firstprivate(chunkId)
      {
         const idVertex lowerBound = chunkId * chunkSize;
         const idVertex upperBound = min(nbScalars, (chunkId + 1) * chunkSize);
         for (idVertex v = lowerBound; v < upperBound; ++v) {
            const auto &neighNumb   = mesh_->getVertexNeighborNumber(v);
            valence     upval       = 0;
            valence     downval     = 0;

            for (valence n = 0; n < neighNumb; ++n) {
               idVertex neigh;
               mesh_->getVertexNeighbor(v, n, neigh);
               if(scalars_->isLower(neigh, v)){
                  ++downval;
               } else {
                  ++upval;
               }
            }

            jt_->setValence(v, downval);
            st_->setValence(v, upval);

            if (!downval) {
               jt_->makeNode(v);
            }

            if (!upval) {
               st_->makeNode(v);
            }
         }
      }
   }

#pragma omp taskwait
   return 0;
}

void ContourTree::build(TreeType tt)
{
   DebugTimer mergeTreesTime;

   const bool ct = tt == TreeType::Contour;

   // if(ct){
   //    DebugTimer precomputeTime;
   //    vertexPrecomputation();
   //    printTime(precomputeTime, "2 precompute ");
   // }

   // -------
   // JT & ST
   // -------

#pragma omp parallel num_threads(threadNumber_)
   {
#pragma omp single nowait
       {
          if (tt == TreeType::Join || ct) {
#pragma omp task untied if(ct)
             jt_->build(ct);
          }
          if (tt == TreeType::Split || ct) {
#pragma omp task untied if(ct)
             st_->build(ct);
          }
       }
#pragma omp taskwait
  }

   printTime(mergeTreesTime, "7 merge trees ");

   // -------
   // Combine
   // -------

   if (tt == TreeType::Contour) {
      if (debugLevel_ >= 4) {
         jt_->printTree2();
         st_->printTree2();
      }

      DebugTimer combineFullTime;
      insertNodes();

      DebugTimer combineTime;
      combine();
      printTime(combineTime, "8 combine trees");
      finalizeSegmentation();
      printTime(combineFullTime, "combine full", -1, 3);
   }

   // -----
   // Debug
   // -----

   if (debugLevel_ > 3) {
      switch (tt) {
         case TreeType::Join:
            jt_->printTree2();
            break;
         case TreeType::Split:
            st_->printTree2();
            break;
         default:
            printTree2();
      }
   } else {
      cout << " nodes :";
      switch (tt) {
         case TreeType::Join:
            cout << jt_->getNumberOfNodes();
            break;
         case TreeType::Split:
            cout << st_->getNumberOfNodes();
            break;
         default:
            cout << getNumberOfNodes();
      }
      cout << endl;
   }
}

void ContourTree::insertNodes(void)
{
   vector<idNode> sortedJTNodes = jt_->sortedNodes(true);
   vector<idNode> sortedSTNodes = st_->sortedNodes(true);

   for (const idNode& t : sortedSTNodes) {

      idVertex vertId = st_->getNode(t)->getVertexId();
      if (jt_->isCorrespondingNode(vertId)) {
          continue;
      } else {
         idSuperArc baseArc = jt_->getCorrespondingSuperArcId(vertId);
         if (jt_->getSuperArc(baseArc)->isMerged()) {
            jt_->updateCorrespondingArc(vertId, jt_->getSuperArc(baseArc)->getReplacantArcId());
         }
      }
      jt_->insertNode(st_->getNode(t));
   }

   for (const idNode& t : sortedJTNodes) {

      idVertex vertId = jt_->getNode(t)->getVertexId();
      if (st_->isCorrespondingNode(vertId)) {
          continue;
      } else {
         idSuperArc baseArc = st_->getCorrespondingSuperArcId(vertId);
         if (st_->getSuperArc(baseArc)->isMerged()) {
            st_->updateCorrespondingArc(vertId, st_->getSuperArc(baseArc)->getReplacantArcId());
         }
      }
      st_->insertNode(jt_->getNode(t));
   }
}

int ContourTree::combine()
{
   DebugTimer stepTime;
   queue<pair<bool, idNode>> growingNodes, remainingNodes;

   const bool DEBUG = false;


   // -------
   // Reserve
   // -------
   treeData_.nodes->reserve(jt_->getNumberOfNodes());
   treeData_.superArcs->reserve(jt_->getNumberOfSuperArcs()+2);
   treeData_.leaves->reserve(jt_->getNumberOfLeaves()+st_->getNumberOfLeaves());

   // ----------------------------------
   // Add JT & ST Leaves to growingNodes
   // ----------------------------------

   // Add leves to growing nodes
   const auto& nbSTLeaves = st_->getNumberOfLeaves();
   if (nbSTLeaves > 1) {
      for (idNode n = 0; n < nbSTLeaves; ++n) {
         const auto &nId = st_->getLeave(n);
         growingNodes.emplace(false, nId);
      }
   } else {
      clone(jt_);
      return 0;
   }

   // count how many leaves can be added, if more than one : ok!
   const auto& nbJTLeaves = jt_->getNumberOfLeaves();
   if (nbJTLeaves > 1) {
      for (idNode n = 0; n < nbJTLeaves; ++n) {
         const auto &nId = jt_->getLeave(n);
         growingNodes.emplace(true, nId);
      }
   } // else can't clone, not same up and down

   if (DEBUG) {
      cout << "growingNodes : " << growingNodes.size() << " in : " << stepTime.getElapsedTime()
           << endl;
   }


   // Warning, have a reserve here, can't make it at the begnining, need build output
   treeData_.leaves->reserve(jt_->getLeaves().size() + st_->getLeaves().size());
   treeData_.superArcs->reserve(jt_->getNumberOfSuperArcs());
   treeData_.nodes->reserve(jt_->getNumberOfNodes());

   if (growingNodes.empty()) {
      cout << "[ContourTree::combine ] Nothing to combine" << endl;
   }

#ifdef withDualQueueCombine
   do {
      while (!remainingNodes.empty()) {
         bool       isJT;
         idNode     currentNodeId;
         MergeTree *xt;

         tie(isJT, currentNodeId) = remainingNodes.front();
         remainingNodes.pop();
         if (isJT) {
            // node come frome jt
            xt = jt_;
         } else {
            // node come from st
            xt = st_;
         }
         if (xt->getNode(currentNodeId)->getNumberOfUpSuperArcs() == 1) {
            growingNodes.emplace(isJT, currentNodeId);
            if (DEBUG) {
               cout << "repush in growing:" << isJT << "::" << xt->printNode(currentNodeId) << endl;
            }
         }
      }
#endif

      while (!growingNodes.empty()) {
         idNode currentNodeId;
         bool   isJT;

         // stats[iteration] = growingNodes.size();
         //++iteration;
         // ------------
         // INFO QUEUE
         // ------------
         // {
         tie(isJT, currentNodeId) = growingNodes.front();
         growingNodes.pop();

         MergeTree *xt = (isJT) ? jt_ : st_;
         MergeTree *yt = (isJT) ? st_ : jt_;

         // }

         // ------------
         // INFO JT / ST
         // ------------
         // {
         // i <- Get(Q)
         const Node *currentNode = xt->getNode(currentNodeId);

         if (DEBUG) {
            if (xt == jt_)
               cout << endl << "JT ";
            else
               cout << endl << "ST ";
            cout << "node : " << currentNode->getVertexId() << endl;
         }

         // "choose a non-root leaf that is not a split in ST" so we ignore such nodes
         if (currentNode->getNumberOfUpSuperArcs() == 0) {
            if (DEBUG) {
               cout << " ignore already processed" << endl;
            }
            continue;
         }

         idNode correspondingNodeId = yt->getCorrespondingNodeId(currentNode->getVertexId());

         if (yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
            if (DEBUG) {
               cout << "put remain:" << isJT << "::" << xt->printNode(currentNodeId) << endl;
               cout << " which is in yt : " << yt->printNode(correspondingNodeId) << endl;
            }
#ifdef withDualQueueCombine
            remainingNodes.emplace(isJT, currentNodeId);
#else
            growingNodes.emplace(isJT, currentNodeId);
#endif
            continue;
         }
         // }

         // -----------
         // NODES IN CT
         // -----------
         // {
         idNode   node1, node2;
         idVertex curVert = currentNode->getVertexId();
         // NODE1
         if (isCorrespondingNode(curVert)) {
            // already a node in the tree
            node1 = getCorrespondingNodeId(curVert);
         } else {
            // create a new node
            node1 = makeNode(currentNode);

            // check if leaf
            if (!currentNode->getNumberOfDownSuperArcs())
               treeData_.leaves->emplace_back(node1);
            else if (!currentNode->getNumberOfUpSuperArcs())
               treeData_.leaves->emplace_back(node1);
         }

         // j <- GetAdj(XT, i)
         idSuperArc curUpArc = currentNode->getUpSuperArcId(0);
         if (xt->getSuperArc(curUpArc)->isMerged()) {
            curUpArc = xt->getSuperArc(curUpArc)->getReplacantArcId();
         }
         idNode      parentId   = xt->getSuperArc(curUpArc)->getUpNodeId();
         const Node *parentNode = xt->getNode(parentId);

         if (DEBUG) {
            cout << " parent node :" << parentNode->getVertexId() << endl;
         }

         idVertex parVert = parentNode->getVertexId();
         // NODE2
         if (isCorrespondingNode(parVert)) {
            // already a node in the tree
            node2 = getCorrespondingNodeId(parVert);
         } else {
            // create a new node
            node2 = makeNode(parentNode);
            if (!parentNode->getNumberOfUpSuperArcs())
               treeData_.leaves->emplace_back(node2);
         }
         // }
         // ----------
         // CREATE ARC
         // ----------
         // {
         idSuperArc processArc = currentNode->getUpSuperArcId(0);

         // create the arc in in the good direction
         // and add it to crossing if needed
         idSuperArc createdArc;
         if (scalars_->isLower(currentNode->getVertexId(),
                               parentNode->getVertexId())) {  // take care of the order
            createdArc = makeSuperArc(node1, node2);
         } else {
            createdArc = makeSuperArc(node2, node1);
         }

         // Segmentation
         createCTArcSegmentation(createdArc, isJT, processArc);

         if (DEBUG) {
            cout << "create arc : " << printArc(createdArc) << endl;
         }

         // }

         // ---------
         // DEL NODES
         // ---------
         // {
         // DelNode(XT, i)
         {
            if (DEBUG) {
               cout << " delete xt (" << (xt == jt_) << ") ";
               cout << "node :" << xt->printNode(currentNodeId) << endl;
            }

            xt->delNode(currentNodeId);
         }

         // DelNode(YT, i)
         {
            if (DEBUG) {
               cout << " delete yt (" << isJT << ") node :";
               cout << yt->printNode(correspondingNodeId) << endl;
            }

            yt->delNode(correspondingNodeId);
         }
         // }

         // -------------
         // PROCESS QUEUE
         // -------------
         // {
         if (parentNode->getNumberOfDownSuperArcs() == 0 && parentNode->getNumberOfUpSuperArcs()) {
            growingNodes.emplace(isJT, parentId);

            if (DEBUG) {
               cout << "will see : " << parentNode->getVertexId() << endl;
            }
         }
         // }
      }
#ifdef withDualQueueCombine
   } while (!remainingNodes.empty());
#endif

   return 0;
}

void ContourTree::createCTArcSegmentation(idSuperArc ctArc, const bool isJT, idSuperArc xtArc)
{
   const MergeTree *xt = (isJT) ? jt_ : st_;

   /*Here we prefere to create lots of small region, each arc having its own segmentation with
    * no overlap instead of having a same vertice in several arc and using vert2tree to decide
    * because we do not want to maintain vert2tree information during the whole computation*/
   const list<Region> &xtRegions = xt->getSuperArc(xtArc)->getRegions();
   for (const Region &reg : xtRegions) {
      segm_it cur    = reg.segmentBegin;
      segm_it end    = reg.segmentEnd;
      segm_it tmpBeg = reg.segmentBegin;
      // each element inside this region
      for (; cur != end; ++cur) {
         if (isCorrespondingNull(*cur)) {
            updateCorrespondingArc(*cur, ctArc);
         } else {
            // already set, we finish a region
            if (cur != tmpBeg) {
               getSuperArc(ctArc)->concat(tmpBeg, cur);
            }
            // if several contiguous vertices are discarded
            // cur will be equals to tmpBeg and we will not create empty regions
            tmpBeg = cur + 1;
         }
      }
      // close last region
      if (cur != tmpBeg) {
         getSuperArc(ctArc)->concat(tmpBeg, cur);
      }
   }
}

void ContourTree::finalizeSegmentation(void)
{
   DebugTimer  finSegmTime;
   const auto &nbArc = getNumberOfSuperArcs();

#pragma omp parallel for schedule(dynamic)
   for (idSuperArc i = 0; i < nbArc; i++) {
      getSuperArc(i)->createSegmentation(scalars_);
   }

   printTime(finSegmTime, "Finalize segm", -1, 3);
}

// }
