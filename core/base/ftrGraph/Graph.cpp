#include "Graph.h"

#include <iostream>
#include <unordered_map>

using namespace std;
using namespace ttk;
using namespace ftr;

Graph::Graph() = default;

Graph::~Graph() = default;

std::string Graph::print(const int verbosity) const {
  stringstream res;
  if(verbosity >= 1) {
    res << "Graph:" << std::endl;
    res << "leaves: " << leaves_.size() << std::endl;
    res << "nodes: " << nodes_.size() << std::endl;
    res << "arcs: " << arcs_.size() << std::endl;
  }

  if(verbosity >= 2) {
    res << "visible arcs: " << getNumberOfVisibleArcs() << std::endl;
  }

  if(verbosity >= 3) {
    res << "Leaves: " << std::endl;
    for(const auto &v : leaves_) {
      res << get<0>(v) << " ";
    }
    res << std::endl;
  }

  if(verbosity >= 4) {
    res << "Nodes:" << std::endl;
    const idNode nbn = nodes_.size();
    for(idNode i = 0; i < nbn; ++i) {
      res << printNode(i) << std::endl;
    }
    res << "Arcs:" << std::endl;
    const idSuperArc nba = arcs_.size();
    for(idSuperArc i = 0; i < nba; ++i) {
      res << printArc(i) << std::endl;
    }
  }

  res << "visible " << getNumberOfVisibleArcs() << " / " << getNumberOfArcs()
      << std::endl;

  return res.str();
}

std::string Graph::printArc(const idSuperArc arcId) const {
  std::stringstream res;
  res << "a" << arcId << ":";
  if(!arcs_[arcId].isVisible()) {
    res << " hidden ";
  }
  if(arcs_[arcId].merged()) {
    res << " merged in " << arcs_[arcId].mergedIn();
  }
  res << "(";
  const idNode did = arcs_[arcId].getDownNodeId();
  const idNode uid = arcs_[arcId].getUpNodeId();
  if(did != nullNode) {
    res << nodes_[did].getVertexIdentifier();
  } else {
    res << "X";
  }
  res << " - ";
  if(uid != nullNode) {
    res << nodes_[uid].getVertexIdentifier();
  } else {
    res << "X";
  }
  res << " ";
  if(arcs_[arcId].isEmpty()) {
    res << "--";
  } else {
    res << arcs_[arcId].getFirstReg() << ".." << arcs_[arcId].getLastReg();
  }
  if(arcs_[arcId].getEnd() != nullVertex) {
    res << " > " << arcs_[arcId].getEnd();
  }
  res << ")";
  return res.str();
}

std::string Graph::printNode(const idNode nodeId) const {
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

std::string Graph::printVisit(const idVertex v) const {
  std::stringstream res;
  if(isArc(v))
    res << " a: " << getArcId(v);
  if(isNode(v))
    res << " n: " << getNodeId(v);
  return res.str();
}

// DEBUG function

std::string Graph::printVisit() const {
  std::stringstream res;
  res << "Segmentation: " << std::endl;
  for(idVertex s = 0; s < nbElmt_; ++s) {
    res << s << " : ";
    res << printVisit(s);
    res << std::endl;
  }
  return res.str();
}

// initialization

void Graph::alloc() {
#ifndef TTK_ENABLE_KAMIKAZE
  if(nbElmt_ == nullVertex) {
    this->printErr("setNumberOfElmt not called before alloc in Graph");
  }
#endif
  leaves_.reserve(nbElmt_);
  nodes_.reserve(nbElmt_ * 2);
  arcs_.reserve(nbElmt_ * 2);
  segmentation_.resize(nbElmt_);
  valUp_.resize(nbElmt_);
  valDown_.resize(nbElmt_);

#ifdef TTK_ENABLE_FTR_VERT_STATS
  nbTouch_.resize(nbElmt_);
  nbArcActif_.resize(nbElmt_);
  avoided_ = 0;
#endif
}

void Graph::init() {
  // Consider max valence is 10
  fillVector<SegmInfo>(segmentation_, SegmInfo{});
  fillVector<valence>(valUp_, -1);
  fillVector<valence>(valDown_, -1);
#ifdef TTK_ENABLE_FTR_VERT_STATS
  fillVector<idVertex>(nbTouch_, 0);
  fillVector<idSuperArc>(nbArcActif_, 0);
#endif
}
