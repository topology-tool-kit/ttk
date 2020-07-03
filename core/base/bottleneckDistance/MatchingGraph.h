#pragma once

class Edge {

public:
  int v1;
  int v2;
  double weight;

  Edge(int p1, int p2, double w) {
    v1 = p1;
    v2 = p2;
    weight = w;
  }

  bool operator<(const Edge &other) const {
    return weight < other.weight;
  }
};
