#include <BettiNumbers.h>

ttk::BettiNumbers::BettiNumbers() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("BettiNumbers");
}

int ttk::BettiNumbers::execute() {
    
    const SimplexId V = triangulation_->getNumberOfVertices();
    const SimplexId E = triangulation_->getNumberOfEdges();
    const SimplexId F = triangulation_->getNumberOfTriangles();
   
    std::vector<ttk::UnionFind> vUF(V); 

    for (SimplexId e=0; e<E; e++){
        SimplexId a,b;
        triangulation_->getEdgeVertex(e, 0, a);
        triangulation_->getEdgeVertex(e, 1, b);
        ttk::UnionFind::makeUnion(&vUF[a], &vUF[b]);
    }

    std::set<ttk::UnionFind*> compV;
    
    for (SimplexId v=0; v<V; v++){
        compV.insert(vUF[v].find());
    }

    this->B0 = static_cast<int>(compV.size());

    this->printMsg({"B0 =", std::to_string(B0)});

    return 1;
}

