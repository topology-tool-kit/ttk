#include <BettiNumbers.h>
#include <set>
#include <map>
#include <vector>
#include <iostream>

ttk::BettiNumbers::BettiNumbers() {
  this->setDebugMsgPrefix("BettiNumbers");
}

int ttk::BettiNumbers::execute() {
  if(!triangulation_) {
    std::cerr << "[BettiNumbers] Erreur : triangulation_ non initialisée." << std::endl;
    return 0;
  }

  Timer timer;

  // --------------------------------------------------------------------
  // 1) Préconditionnement
  // --------------------------------------------------------------------
  triangulation_->preconditionVertexNeighbors();
  triangulation_->preconditionEdges();
  triangulation_->preconditionTriangles();
  triangulation_->preconditionTriangleEdges();
  triangulation_->preconditionTriangleStars();
  triangulation_->preconditionEdgeTriangles();

  // --------------------------------------------------------------------
  // 2) Comptages globaux
  // --------------------------------------------------------------------
  const SimplexId V = triangulation_->getNumberOfVertices();
  const SimplexId E = triangulation_->getNumberOfEdges();
  const SimplexId F = triangulation_->getNumberOfTriangles();

  // --------------------------------------------------------------------
  // 3) B0 (composantes de volume)
  // --------------------------------------------------------------------
  std::vector<ttk::UnionFind> vUF(V);
  for(SimplexId e = 0; e < E; ++e) {
    SimplexId a, b;
    triangulation_->getEdgeVertex(e, 0, a);
    triangulation_->getEdgeVertex(e, 1, b);
    ttk::UnionFind::makeUnion(&vUF[a], &vUF[b]);
  }
  std::set<ttk::UnionFind*> compV;
  for(SimplexId v = 0; v < V; ++v)
    compV.insert(vUF[v].find());
  this->B0_ = static_cast<int>(compV.size());

  // --------------------------------------------------------------------
  // 4) Extraction des triangles de bord
  // --------------------------------------------------------------------
  std::vector<SimplexId> boundaryTris;
  boundaryTris.reserve(F);
  for(SimplexId t = 0; t < F; ++t) {
    if(triangulation_->getTriangleStarNumber(t) == 1)
      boundaryTris.push_back(t);
  }
  const SimplexId BT = static_cast<SimplexId>(boundaryTris.size());

  // --------------------------------------------------------------------
  // 5) B0(boundary) par Union-Find sur ces triangles
  // --------------------------------------------------------------------
  std::vector<ttk::UnionFind> triUF(BT);
  std::map<SimplexId, SimplexId> triIndex;
  for(SimplexId i = 0; i < BT; ++i)
    triIndex[boundaryTris[i]] = i;

  for(SimplexId i = 0; i < BT; ++i) {
    const SimplexId t = boundaryTris[i];
    for(int j = 0; j < 3; ++j) {
      SimplexId e;
      triangulation_->getTriangleEdge(t, j, e);
      const int nt = triangulation_->getEdgeTriangleNumber(e);
      for(int k = 0; k < nt; ++k) {
        SimplexId t2;
        triangulation_->getEdgeTriangle(e, k, t2);
        if(t2 == t) continue;
        auto it = triIndex.find(t2);
        if(it != triIndex.end())
          ttk::UnionFind::makeUnion(&triUF[i], &triUF[it->second]);
      }
    }
  }
  std::set<ttk::UnionFind*> compT;
  for(SimplexId i = 0; i < BT; ++i)
    compT.insert(triUF[i].find());
  const int B0_boundary = static_cast<int>(compT.size());

  // --------------------------------------------------------------------
  // 6) B2 selon Dey & Guha : B2 = B0(boundary) - B0(volume)
  // --------------------------------------------------------------------
  this->B2_ = B0_boundary - this->B0_;

  // --------------------------------------------------------------------
  // 7) B1 = somme des genres des composantes de bord
  // --------------------------------------------------------------------
  // Regrouper triangles par composante
  std::map<ttk::UnionFind*, std::vector<SimplexId>> trisByComp;
  for(SimplexId i = 0; i < BT; ++i)
    trisByComp[triUF[i].find()].push_back(boundaryTris[i]);

  this->B1_ = 0;
  for(auto &kv : trisByComp) {
    auto &tris = kv.second;
    std::set<SimplexId> vs, es;

    // Correction : on insère **toutes** les arêtes de la surface,
    // pas seulement celles à triangle-count < 2.
    for(auto t : tris) {
      // sommets
      for(int j = 0; j < 3; ++j) {
        SimplexId v;
        triangulation_->getTriangleVertex(t, j, v);
        vs.insert(v);
      }
      // arêtes
      for(int j = 0; j < 3; ++j) {
        SimplexId e;
        triangulation_->getTriangleEdge(t, j, e);
        es.insert(e);
      }
    }

    const int Fi = static_cast<int>(tris.size());
    const int Vi = static_cast<int>(vs.size());
    const int Ei = static_cast<int>(es.size());
    const int chi = Vi - Ei + Fi;      // caractéristique d’Euler de la surface
    const int gi  = (2 - chi) / 2;     // genre  (Th. 3.3)
    this->B1_ += gi;
  }

  // --------------------------------------------------------------------
  // 8) Affichage des résultats
  // --------------------------------------------------------------------
  this->printMsg({
    {"#V",  std::to_string(V)},
    {"#E",  std::to_string(E)},
    {"#F",  std::to_string(F)},
    {"#BT", std::to_string(BT)},
    {"B0",  std::to_string(B0_)},
    {"B1",  std::to_string(B1_)},
    {"B2",  std::to_string(B2_)}
  });
  this->printMsg("Complete", timer.getElapsedTime(), this->threadNumber_);

  return 1;
}

