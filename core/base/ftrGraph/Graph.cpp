#include "Graph.h"

#include<iostream>

using namespace std;
using namespace ttk;
using namespace ftr;

Graph::Graph()
{
}

Graph::~Graph()
{
}

void Graph::print(const int verbosity) const
{
   if (verbosity >= 1) {
      cout << "Graph:" << endl;
      cout << "leaves: " << leaves_.size() << endl;
      cout << "nodes: " << nodes_.size() << endl;
      cout << "arcs: " << arcs_.size() << endl;
   }

   if(verbosity >= 2) {
      cout << "Leaves: " << endl;
      for (const idVertex v : leaves_) {
         cout << v << " ";
      }
      cout << endl;
   }

   if(verbosity >= 3) {
      cout << "Nodes:" << endl;
      const idNode nbn = nodes_.size();
      for(idNode i = 0; i < nbn; ++i) {
         cout << printNode(i) << endl;
      }
      cout << "Arcs:" << endl;
      const idSuperArc nba = arcs_.size();
      for(idSuperArc i = 0; i < nba; ++i) {
         cout << printArc(i) << endl;
      }
   }
}

void Graph::alloc()
{
#ifndef TTK_ENABLE_KAMIKAZE
   if (nbVerts_ == nullVertex) {
      cout << "[FTR Graph]: ERROR, setNumberOfVertices not called before alloc in Graph" << endl;
   }
#endif
   leaves_.reserve(nbVerts_/3);
   nodes_.reserve(nbVerts_/2);
   arcs_.reserve(nbVerts_/2);
   segmentation_.resize(nbVerts_);
   valences_.resize(nbVerts_);
}

void Graph::init()
{
   fillVector<std::forward_list<idSegmentation>>(segmentation_,
                                                 std::forward_list<idSegmentation>{});
   fillVector<valence>(valences_, -1);
}

// private

std::string Graph::printArc(const idSuperArc arcId) const
{
   std::stringstream res;
   res << "a" << arcId << ":(";
   const idNode did = arcs_[arcId].getDownNodeId();
   const idNode uid = arcs_[arcId].getUpNodeId();
   if (did != nullNode){
      res << nodes_[did].getVertexIdentifier();
   } else {
      res << "X";
   }
   res << " - ";
   if (uid != nullNode){
      res << nodes_[uid].getVertexIdentifier();
   } else {
      res << "X";
   }
   res << ")";
   return res.str();
}

std::string Graph::printNode(const idNode nodeId) const
{
   std::stringstream res;
   res << "n" << nodeId << "v" << nodes_[nodeId].getVertexIdentifier() << ":[v";
   const idSuperArc nbd = nodes_[nodeId].getNbDownArcs();
   for(idSuperArc i = 0; i < nbd; ++i) {
      res << " " << nodes_[nodeId].getDownArc(i);
   }
   res << " - ^";
   const idSuperArc nbu = nodes_[nodeId].getNbUpArcs();
   for(idSuperArc i = 0; i < nbu; ++i) {
      res << " " << nodes_[nodeId].getUpArc(i);
   }
   res << "]";
   return res.str();
}


std::string Graph::printVisit(const idVertex v) const
{
   std::stringstream res;
   for (const idSegmentation tmp : segmentation_[v]) {
      res << " " << tmp;
   }
   return res.str();
}

