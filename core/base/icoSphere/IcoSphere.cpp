#include <IcoSphere.h>

typedef tuple<float, float, float> Vertex;

typedef long long IndexType;
typedef tuple<IndexType, IndexType, IndexType> Triangle;

typedef vector<Vertex> VertexList;
typedef vector<Triangle> TriangleList;
typedef map<pair<IndexType, IndexType>, IndexType> EdgeVertexMap;
typedef pair<VertexList, TriangleList> IndexedMesh;

ttk::IcoSphere::IcoSphere() {
}
ttk::IcoSphere::~IcoSphere() {
}

using namespace std;

Vertex normalizeVertex(Vertex &a) {
  Vertex result;

  float length = sqrt(get<0>(a) * get<0>(a) + get<1>(a) * get<1>(a)
                      + get<2>(a) * get<2>(a));

  get<0>(result) = get<0>(a) / length;
  get<1>(result) = get<1>(a) / length;
  get<2>(result) = get<2>(a) / length;

  return result;
}

Vertex addVertices(Vertex &a, Vertex &b) {
  Vertex result;

  get<0>(result) = get<0>(a) + get<0>(b);
  get<1>(result) = get<1>(a) + get<1>(b);
  get<2>(result) = get<2>(a) + get<2>(b);

  return result;
}

IndexType vertex_for_edge(EdgeVertexMap &edgeVertexMap,
                          VertexList &vertices,
                          IndexType first,
                          IndexType second) {
  EdgeVertexMap::key_type key(first, second);
  if(key.first > key.second)
    swap(key.first, key.second);

  auto inserted = edgeVertexMap.insert({key, vertices.size()});
  if(inserted.second) {
    auto &edge0 = vertices[first];
    auto &edge1 = vertices[second];
    auto temp = addVertices(edge0, edge1);
    auto vertex = normalizeVertex(temp);
    vertices.push_back(vertex);
  }

  return inserted.first->second;
}

TriangleList subdivide(VertexList &vertices, TriangleList &triangles) {
  EdgeVertexMap edgeVertexMap;
  TriangleList results;

  size_t maxIndex = triangles.size();

  for(size_t i = 0; i < maxIndex; i++) {
    Triangle &each = triangles[i];

    IndexType mid[3];

    mid[0]
      = vertex_for_edge(edgeVertexMap, vertices, get<0>(each), get<1>(each));
    mid[1]
      = vertex_for_edge(edgeVertexMap, vertices, get<1>(each), get<2>(each));
    mid[2]
      = vertex_for_edge(edgeVertexMap, vertices, get<2>(each), get<0>(each));

    results.push_back(make_tuple(get<0>(each), mid[0], mid[2]));
    results.push_back(make_tuple(get<1>(each), mid[1], mid[0]));
    results.push_back(make_tuple(get<2>(each), mid[2], mid[1]));
    results.push_back(make_tuple(mid[0], mid[1], mid[2]));
  }
  return results;
}

int ttk::IcoSphere::generate(
  // Input
  size_t subdivisions,
  float radius,
  float *center,

  // Output
  vector<tuple<float, float, float>> &vertices,
  vector<tuple<long long, long long, long long>> &triangles) const {

  Timer t;

  const float X = .525731112119133606f;
  const float Z = .850650808352039932f;
  const float N = 0.f;

  vertices
    = {make_tuple(-X, N, Z), make_tuple(X, N, Z),   make_tuple(-X, N, -Z),
       make_tuple(X, N, -Z), make_tuple(N, Z, X),   make_tuple(N, Z, -X),
       make_tuple(N, -Z, X), make_tuple(N, -Z, -X), make_tuple(Z, X, N),
       make_tuple(-Z, X, N), make_tuple(Z, -X, N),  make_tuple(-Z, -X, N)};

  triangles = {make_tuple(0, 4, 1),  make_tuple(0, 9, 4),  make_tuple(9, 5, 4),
               make_tuple(4, 5, 8),  make_tuple(4, 8, 1),  make_tuple(8, 10, 1),
               make_tuple(8, 3, 10), make_tuple(5, 3, 8),  make_tuple(5, 2, 3),
               make_tuple(2, 7, 3),  make_tuple(7, 10, 3), make_tuple(7, 6, 10),
               make_tuple(7, 11, 6), make_tuple(11, 0, 6), make_tuple(0, 1, 6),
               make_tuple(6, 1, 10), make_tuple(9, 0, 11), make_tuple(9, 11, 2),
               make_tuple(9, 2, 5),  make_tuple(7, 2, 11)};

  for(size_t i = 0; i < subdivisions; i++)
    triangles = subdivide(vertices, triangles);

  size_t n = vertices.size();
  for(size_t i = 0; i < n; i++) {
    Vertex &v = vertices[i];
    get<0>(v) = get<0>(v) * radius + center[0];
    get<1>(v) = get<1>(v) * radius + center[1];
    get<2>(v) = get<2>(v) * radius + center[2];
  }

  // Print performance
  {
    stringstream msg;
    msg << "[ttkIcoSphere] Generated in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}