#pragma once

#include <Debug.h>
#include <PersistenceDiagramUtils.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

namespace ttk {

  class GabowTarjan : public Debug {

    struct Edge {
      int v1{};
      int v2{};
      double weight{};

      Edge(int p1, int p2, double w) : v1{p1}, v2{p2}, weight{w} {
      }

      bool operator<(const Edge &other) const {
        return weight < other.weight;
      }
    };

  public:
    GabowTarjan() {
      this->setDebugMsgPrefix("Gabow-Tarjan");
    }

    double Distance();

    void printCurrentMatching();

    int run(std::vector<MatchingType> &matchings);

    inline void setInput(int rowSize_,
                         int colSize_,
                         std::vector<std::vector<double>> *C_) {
      Cptr = C_;

      Size1 = (unsigned int)rowSize_ - 1;
      Size2 = (unsigned int)colSize_ - 1;
      if(Size1 <= 0 || Size2 <= 0) {
        this->printMsg("One or more empty diagram(s).");
      }

      MaxSize = Size1 + Size2;
      Edges.clear();

      // Connect diagonal points.
      for(unsigned int i = Size1; i < MaxSize; ++i)
        for(unsigned int j = MaxSize + Size2; j < 2 * MaxSize; ++j) {
          const Edge localEdge(i, j, (double)0);
          Edges.emplace_back(localEdge);
        }

      // Connect real points.
      for(unsigned int i = 0; i < Size1; ++i) {
        unsigned int k = MaxSize;
        for(unsigned int j = 0; j < Size2; ++j) {
          auto val = (*C_)[i][j];
          const Edge localEdge(i, k++, val);
          Edges.emplace_back(localEdge);
        }
      }

      // Connect real points with their diagonal.
      for(unsigned int i = 0; i < Size1; ++i) {
        auto val = (*C_)[i][Size2];
        const Edge localEdge(i, MaxSize + Size2 + i, val);
        Edges.emplace_back(localEdge);
      }

      for(unsigned int j = 0, k = MaxSize; j < Size2; ++j, ++k) {
        auto val = (*C_)[Size1][j];
        const Edge localEdge(Size1 + (k - MaxSize), k, val);
        Edges.emplace_back(localEdge);
      }

      std::sort(Edges.begin(), Edges.end());
    }

    inline void clear() {
      MaxSize = 0;
      Size1 = 0;
      Size2 = 0;
      Edges.clear();
      Pair.clear();
      Connections.clear();
      Layers.clear();
    }

  private:
    // Original cost matrix.
    std::vector<std::vector<double>> *Cptr;

    /*
     * Total number of persistencePairs
     */
    unsigned int MaxSize;
    unsigned int Size1;
    unsigned int Size2;

    /*
     * Edges between all nodes given by
     * persistencePairs in both persistence diagrams
     */
    std::vector<Edge> Edges;

    /*
     * Pairing between vertices in persistence diagrams,
     * -1 means the vertex is not paired to any other vertex.
     * Non negative number, index of the vertex to
     * which the given vertex is paired.
     */
    std::vector<int> Pair;

    /*
     * Connection matrix, gives edges between
     * vertices in the first and second diagram
     */
    std::vector<std::vector<int>> Connections;

    /*
     * Layers used in the Hopcroft-Karp algorithm
     * The layer information is required for a NIL vertex
     * and all vertices in PersistencePairs1
     * Hence 0 slot is used for NIL and all the vertices
     * in PersistencePairs1 are shifted by one.
     * So to read layer of the vertex 0 we access Layers[1]
     */
    std::vector<int> Layers;

    bool DFS(int v);
    bool BFS();

    // Hopcroft-Karp algorithm: find a maximal matching
    void HopcroftKarp(unsigned int &matching);
  };

} // namespace ttk
