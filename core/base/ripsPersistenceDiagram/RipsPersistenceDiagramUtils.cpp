#include <RipsPersistenceDiagramUtils.h>

ttk::rpd::UnionFind::UnionFind(unsigned n) {
  parent_.resize(n);
  rank_.resize(n, 0);
  for(unsigned i = 0; i < n; i++)
    parent_[i] = i;
}

int ttk::rpd::UnionFind::find(int x) {
  if(parent_[x] == x)
    return x;
  return parent_[x] = find(parent_[x]); // path compression
}

void ttk::rpd::UnionFind::merge(int x, int y) {
  const int rootX = find(x);
  const int rootY = find(y);
  if(rootX != rootY) {
    if(rank_[rootX] > rank_[rootY])
      parent_[rootY] = rootX;
    else if(rank_[rootX] < rank_[rootY])
      parent_[rootX] = rootY;
    else {
      parent_[rootY] = rootX;
      rank_[rootX]++;
    }
  }
}

int ttk::rpd::UnionFind::mergeRet(int x, int y) {
  const int rootX = find(x);
  const int rootY = find(y);
  if(rootX != rootY) {
    if(rank_[rootX] > rank_[rootY])
      return parent_[rootY] = rootX;
    else if(rank_[rootX] < rank_[rootY])
      return parent_[rootX] = rootY;
    else {
      rank_[rootX]++;
      return parent_[rootY] = rootX;
    }
  }
  return rootX;
}

bool ttk::rpd::UnionFind::isRoot(int x) const {
  return parent_[x] == x;
}