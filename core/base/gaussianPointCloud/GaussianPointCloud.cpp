#include <GaussianPointCloud.h>

typedef tuple<float,float,float> Vertex;

typedef long long IndexType;
typedef tuple<IndexType,IndexType,IndexType> Triangle;

typedef vector<Vertex> VertexList;
typedef vector<Triangle> TriangleList;
typedef map< pair<IndexType,IndexType>, IndexType > EdgeVertexMap;
typedef pair<VertexList, TriangleList> IndexedMesh;

ttk::GaussianPointCloud::GaussianPointCloud(){}
ttk::GaussianPointCloud::~GaussianPointCloud(){}

using namespace std;

int ttk::GaussianPointCloud::generate(
    // Input
    size_t subdivisions,
    float radius,
    float* center,

    // Output
    vector< tuple<float,float,float> >& vertices,
    vector< tuple<long long,long long,long long> >& triangles
) const{

    Timer t;

    const float X=.525731112119133606f;
    const float Z=.850650808352039932f;
    const float N=0.f;

    vertices = {
        make_tuple(-X,N,Z), make_tuple(X,N,Z), make_tuple(-X,N,-Z), make_tuple(X,N,-Z),
        make_tuple(N,Z,X), make_tuple(N,Z,-X), make_tuple(N,-Z,X), make_tuple(N,-Z,-X),
        make_tuple(Z,X,N), make_tuple(-Z,X, N), make_tuple(Z,-X,N), make_tuple(-Z,-X, N)
    };

    triangles = {
        make_tuple(0,4,1),make_tuple(0,9,4),make_tuple(9,5,4),make_tuple(4,5,8),make_tuple(4,8,1),
        make_tuple(8,10,1),make_tuple(8,3,10),make_tuple(5,3,8),make_tuple(5,2,3),make_tuple(2,7,3),
        make_tuple(7,10,3),make_tuple(7,6,10),make_tuple(7,11,6),make_tuple(11,0,6),make_tuple(0,1,6),
        make_tuple(6,1,10),make_tuple(9,0,11),make_tuple(9,11,2),make_tuple(9,2,5),make_tuple(7,2,11)
    };

//     for(size_t i=0; i<subdivisions; i++)
//         triangles = subdivide(vertices, triangles);

    size_t n = vertices.size();
    for(size_t i=0; i<n; i++){
        Vertex& v = vertices[i];
        get<0>(v) = get<0>(v)*radius+center[0];
        get<1>(v) = get<1>(v)*radius+center[1];
        get<2>(v) = get<2>(v)*radius+center[2];
    }

    // Print performance
    {
        stringstream msg;
        msg << "[ttkGaussianPointCloud] Generated in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 0;
}
