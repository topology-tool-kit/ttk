/// \ingroup base
/// \class ttk::SuperlevelSetOverlapTracking
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %superlevelSetOverlapTracking processing package.
///
/// %SuperlevelSetOverlapTracking is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkSuperlevelSetOverlapTracking.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <unordered_set>

using namespace std;

namespace ttk{

    struct Node {
        size_t size;
        float x;
        float y;
        float z;
        size_t branch;
    };

    struct Edge {
        size_t n0;
        size_t n1;
        size_t branch;

        bool operator==(const Edge &other) const {
            return (n0 == other.n0 && n1 == other.n1);
        }
    };

    struct EdgeHasher {
        size_t operator()(const Edge& e) const {
            using std::size_t;
            using std::hash;
            using std::string;

            return ((hash<int>()(e.n0) ^ (hash<int>()(e.n1) << 1)) >> 1);
        }
    };

    class SuperlevelSetOverlapTracking : public Debug{

        public:

        SuperlevelSetOverlapTracking(){
            dim.resize(3,0);
            noLabel = -1;
            emptyLabel = -2;
        };

        ~SuperlevelSetOverlapTracking(){};

        template <class dataType>
            int computeLabels(
                const dataType &level,
                dataType* data,
                float* labels,
                vector<Node>& nodes,
                size_t timeIndex
            ) const;

        template <class dataType>
            int execute(
                float* labels0,
                float* labels1,
                vector<Edge>& edges
            ) const;

        inline int setDim(int* dim){
            this->dim[0] = dim[0];
            this->dim[1] = dim[1];
            this->dim[2] = dim[2];

            this->size2 = dim[0]*dim[1];
            this->size3 = this->size2*dim[2];

            return 0;
        }

        inline size_t xyzToIndex(size_t x, size_t y, size_t z) const {
            return z*size2 + y*dim[0] + x;
        }

        inline void indexToXYZ(size_t& index, size_t& x, size_t& y, size_t& z) const {
            z = index/size2;
            size_t w = index-z*size2;
            y = w/dim[0];
            x = w%dim[0];
        }

        inline void dialateCell(float* labels, size_t& x, size_t& y, size_t& z, const float& label) const {
            labels[ xyzToIndex(x+1,y  ,z  ) ] = label;
            labels[ xyzToIndex(x-1,y  ,z  ) ] = label;
            labels[ xyzToIndex(x  ,y+1,z  ) ] = label;
            labels[ xyzToIndex(x  ,y-1,z  ) ] = label;
            labels[ xyzToIndex(x  ,y  ,z+1) ] = label;
            labels[ xyzToIndex(x  ,y  ,z-1) ] = label;
        }

        inline void dialate(
                float* labels,
                const float& label
            ) const {

            vector<float> temp(size3);
            for(size_t i=0; i<size3; i++)
                temp[i] = labels[i];

            for(size_t x = dim[0]-2; x>1; x--)
                for(size_t y = dim[1]-2; y>1; y--)
                    for(size_t z = dim[2]-2; z>1; z--)
                        if( temp[xyzToIndex(x,y,z)]==label )
                            dialateCell(labels, x, y, z, label);
        }

        inline void processCell(size_t& x, size_t& y, size_t& z, float* labels, vector<size_t>& stack, size_t& stackIndex, float& id) const {
            size_t i = xyzToIndex(x, y, z);
            if(labels[i]==emptyLabel){
                labels[i] = id;
                stack[++stackIndex] = i;
            }
        }

        inline int processNeighbours(float* labels, vector<size_t>& stack, size_t& stackIndex, size_t& vertexIndex, float& id) const {

            size_t x,y,z;
            indexToXYZ(vertexIndex, x, y, z);
            size_t x0 = x==0 ? 0 : x-1;
            size_t y0 = y==0 ? 0 : y-1;
            size_t z0 = z==0 ? 0 : z-1;
            size_t x1 = min(dim[0]-1,x+1);
            size_t y1 = min(dim[1]-1,y+1);
            size_t z1 = min(dim[2]-1,z+1);

            // Top
            processCell(x0, y0, z0, labels, stack, stackIndex, id);
            processCell(x0,  y, z0, labels, stack, stackIndex, id);
            processCell(x0, y1, z0, labels, stack, stackIndex, id);

            processCell( x, y0, z0, labels, stack, stackIndex, id);
            processCell( x,  y, z0, labels, stack, stackIndex, id);
            processCell( x, y1, z0, labels, stack, stackIndex, id);

            processCell(x1, y0, z0, labels, stack, stackIndex, id);
            processCell(x1,  y, z0, labels, stack, stackIndex, id);
            processCell(x1, y1, z0, labels, stack, stackIndex, id);

            // Middle
            processCell(x0, y0,  z, labels, stack, stackIndex, id);
            processCell(x0,  y,  z, labels, stack, stackIndex, id);
            processCell(x0, y1,  z, labels, stack, stackIndex, id);

            processCell( x, y0,  z, labels, stack, stackIndex, id);
            processCell( x, y1,  z, labels, stack, stackIndex, id);

            processCell(x1, y0,  z, labels, stack, stackIndex, id);
            processCell(x1,  y,  z, labels, stack, stackIndex, id);
            processCell(x1, y1,  z, labels, stack, stackIndex, id);

            // Bottom
            processCell(x0, y0, z1, labels, stack, stackIndex, id);
            processCell(x0,  y, z1, labels, stack, stackIndex, id);
            processCell(x0, y1, z1, labels, stack, stackIndex, id);

            processCell( x, y0, z1, labels, stack, stackIndex, id);
            processCell( x,  y, z1, labels, stack, stackIndex, id);
            processCell( x, y1, z1, labels, stack, stackIndex, id);

            processCell(x1, y0, z1, labels, stack, stackIndex, id);
            processCell(x1,  y, z1, labels, stack, stackIndex, id);
            processCell(x1, y1, z1, labels, stack, stackIndex, id);

            return 0;
        }

        inline int floodFill3D(
            float* labels, vector<size_t>& stack, size_t& vertexIndex0, float& id, size_t& size
        ) const {
            size_t stackIndex = 1;
            stack[stackIndex] = vertexIndex0;

            labels[vertexIndex0] = id;

            size = 0;
            while(stackIndex>0){
                size++;
                size_t vertexIndex = stack[stackIndex--];
                processNeighbours(labels, stack, stackIndex, vertexIndex, id);
            }

            return 0;
        }

        inline int computeBranches(
            vector< vector<Node> >& timeNodeMap,
            vector< vector<Edge> >& timeEdgeMap
        ) const {
            const size_t nT = timeNodeMap.size();
            const size_t NO_ID = numeric_limits<size_t>::max();

            // -----------------------------------------------------------------
            // Compute Max Predecessors
            // -----------------------------------------------------------------
            // <ID, size, degree>
            vector< vector<tuple<size_t,size_t,size_t>> > timeNodeMaxPredMap(nT); // One too many but ok
            vector< vector<tuple<size_t,size_t,size_t>> > timeNodeMaxSuccMap(nT); // One too many but ok

            // Prepare Maps
            for(size_t t=0; t<nT; t++){
                size_t n = timeNodeMap[t].size();

                auto& nodeMaxPredMap = timeNodeMaxPredMap[t];
                auto& nodeMaxSuccMap = timeNodeMaxSuccMap[t];

                nodeMaxPredMap.resize( n );
                nodeMaxSuccMap.resize( n );

                for(size_t i=0; i<n; i++){
                    get<0>(nodeMaxPredMap[i]) = NO_ID;
                    get<0>(nodeMaxSuccMap[i]) = NO_ID;
                }
            }

            // Find Max Succ and Pred
            for(size_t t=0; t<nT-1; t++){
                auto& edges = timeEdgeMap[t];
                auto& nodes0 = timeNodeMap[t];
                auto& nodes1 = timeNodeMap[t+1];

                auto& node0MaxSuccMap = timeNodeMaxSuccMap[t];
                auto& node1MaxPredMap = timeNodeMaxPredMap[t+1];

                for(auto& e: edges){
                    auto& n0 = nodes0[ e.n0 ];
                    auto& n1 = nodes1[ e.n1 ];
                    auto& maxPred = node1MaxPredMap[ e.n1 ];
                    get<2>(maxPred)++;

                    auto& maxSucc = node0MaxSuccMap[ e.n0 ];
                    get<2>(maxSucc)++;

                    if( get<1>(maxPred) <= n0.size ){
                        get<0>(maxPred) = e.n0;
                        get<1>(maxPred) = n0.size;
                    }
                    if( get<1>(maxSucc) <= n1.size ){
                        get<0>(maxSucc) = e.n1;
                        get<1>(maxSucc) = n1.size;
                    }
                }
            }

            // Initiallize Branches
            {
                for(size_t t=0; t<nT-1; t++){
                    auto& edges = timeEdgeMap[t];
                    for(auto& e: edges)
                        e.branch = NO_ID;
                    auto& nodes = timeNodeMap[t];
                    for(auto& n: nodes){
                        n.branch = NO_ID;
                    }
                }

                auto& nodes = timeNodeMap[nT-1];
                for(auto& n: nodes)
                    n.branch = NO_ID;
            }

            // Create Initial Labels
            size_t branchCounter = 0;
            for(size_t t=0; t<nT-1; t++){
                auto& edges = timeEdgeMap[t];

                auto& node0MaxPredMap = timeNodeMaxPredMap[t];
                auto& node0MaxSuccMap = timeNodeMaxSuccMap[t];
                auto& node1MaxPredMap = timeNodeMaxPredMap[t+1];

                auto& nodes0 = timeNodeMap[t];
                auto& nodes1 = timeNodeMap[t+1];

                for(auto& e: edges){
                    auto& n0 = nodes0[e.n0];
                    auto& n1 = nodes1[e.n1];

                    size_t inDegree = get<2>(node0MaxPredMap[e.n0]);

                    if(inDegree==0){ // If Birth Node
                        // cout<<"B "<<t<<"_"<<e.n0<<" - "<<t+1<<"_"<<e.n1<<": "<<get<0>(node0MaxSuccMap[e.n0]) << " " << get<0>(node1MaxPredMap[e.n1]) <<endl;

                        if(n0.branch==NO_ID) // Create first label
                            n0.branch = branchCounter++;

                        // if(get<0>(node0MaxSuccMap[e.n0])==e.n1 && get<0>(node1MaxPredMap[e.n1])==e.n0)
                        //     n1.branch = n0.branch;
                        // else if (n1.branch==NO_ID)
                        //     n1.branch = branchCounter++;

                        // if(n1.branch==NO_ID)
                    }

                    size_t outDegree = get<2>(node0MaxSuccMap[e.n0]);
                    if (outDegree>1){ // If Saddle
                        // cout<<"S "<<t<<"_"<<e.n0<<" - "<<t+1<<"_"<<e.n1<<": "<<get<0>(node0MaxSuccMap[e.n0]) << " " << get<0>(node1MaxPredMap[e.n1]) <<endl;
                        // cout<<!(get<0>(node0MaxSuccMap[e.n0])==e.n1 && get<0>(node1MaxPredMap[e.n1])==e.n0)<<endl;
                        if(n1.branch==NO_ID && !(get<0>(node0MaxSuccMap[e.n0])==e.n1 && get<0>(node1MaxPredMap[e.n1])==e.n0) )
                            n1.branch = branchCounter++;
                    }
                }
            }

            // Label all followers with no Label
            for(size_t t=0; t<nT-1; t++){
                auto& edges = timeEdgeMap[t];
                auto& nodes0 = timeNodeMap[t];
                auto& nodes1 = timeNodeMap[t+1];
                auto& node0MaxSuccMap = timeNodeMaxSuccMap[t];
                auto& node1MaxPredMap = timeNodeMaxPredMap[t+1];

                for(auto& e: edges){
                    auto& n0 = nodes0[e.n0];
                    auto& n1 = nodes1[e.n1];

                    if(n1.branch == NO_ID && get<0>(node0MaxSuccMap[e.n0])==e.n1 && get<0>(node1MaxPredMap[e.n1])==e.n0)
                        n1.branch = n0.branch;
                }
            }

            for(size_t t=0; t<nT-1; t++){
                auto& edges = timeEdgeMap[t];
                auto& nodes0 = timeNodeMap[t];
                auto& nodes1 = timeNodeMap[t+1];
                auto& node1MaxPredMap = timeNodeMaxPredMap[t+1];

                for(auto& e: edges){
                    auto& n0 = nodes0[e.n0];
                    auto& n1 = nodes1[e.n1];

                    e.branch = n0.branch==n1.branch
                        ? n0.branch
                        : get<0>(node1MaxPredMap[e.n1])==e.n0
                            ? n1.branch
                            : n0.branch;
                }
            }

            return 0;
        }

        protected:

        vector<size_t>    dim;
        size_t    size2;
        size_t    size3;
        float    noLabel;
        float    emptyLabel;
    };
}

// template functions
template <class dataType> int ttk::SuperlevelSetOverlapTracking::computeLabels(
    const dataType &level,
    dataType* data,
    float* labels,
    vector<Node>& nodes,
    size_t timeIndex
) const{

    // Initial Labels
    for(size_t i=0; i<size3; i++)
        labels[i] = data[i]<level ? noLabel : emptyLabel;
        // labels[i] = data[i]>=level ? noLabel : emptyLabel;

    dialate(labels, noLabel);
    dialate(labels, noLabel);
    dialate(labels, noLabel);

    dialate(labels, emptyLabel);
    dialate(labels, emptyLabel);

    float id = 0;
    vector<size_t> stack(size3);

    for(size_t i=0; i<size3; i++){
        if(labels[i]!=emptyLabel) continue;

        size_t size;
        floodFill3D(labels, stack, i, id, size);

        nodes.push_back( {size, (float) timeIndex, (float) id, 0.f, 0} );

        id++;
    }

    return 0;
}

template <class dataType> int ttk::SuperlevelSetOverlapTracking::execute(
    float* labels0,
    float* labels1,
    vector<Edge>& edges
) const{

    Timer t;

    size_t size3 = dim[0]*dim[1]*dim[2];

    unordered_set<Edge,EdgeHasher> edgesTemp;

    for(size_t i=0; i<size3; i++){
        if( labels0[i]!=noLabel && labels1[i]!=noLabel ){ // if there is overlap
            edgesTemp.insert( {(size_t) labels0[i], (size_t) labels1[i], 0} );
        }
    }

    edges.resize( edgesTemp.size() );
    size_t i=0;
    for(auto& e: edgesTemp){
        edges[i++] = e;
    }

    // {
    //     stringstream msg;
    //     // msg << "[SuperlevelSetOverlapTracking] Data-set (" << vertexNumber
    //     //   << " points) processed in "
    //     //   << t.getElapsedTime() << " s. (" << threadNumber_
    //     //   << " thread(s))."
    //     //   << endl;
    //     dMsg(cout, msg.str(), timeMsg);
    // }

  return 0;
}
