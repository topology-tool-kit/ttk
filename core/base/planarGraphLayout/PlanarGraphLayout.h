/// \ingroup base
/// \class ttk::PlanarGraphLayout
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %planarGraphLayout processing package.
///
/// %PlanarGraphLayout is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkPlanarGraphLayout.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include                  <unordered_set>

#ifdef TTK_ENABLE_GRAPHVIZ
#include <graphviz/cgraph.h>
#include <graphviz/gvc.h>
#endif

using namespace std;

namespace ttk{

    class PlanarGraphLayout : public Debug{

        public:
            PlanarGraphLayout(){};
            ~PlanarGraphLayout(){};

            template <class dataType> int execute(
                Triangulation* triangulation,
                size_t* nodeBranchMap,
                size_t nTime,
                dataType* outputVertices
            ) const;
    };
}

// template functions
template <class dataType> int ttk::PlanarGraphLayout::execute(
    Triangulation* triangulation,
    size_t* nodeBranchMap,
    size_t nTime,
    dataType* outputVertices
) const{

    triangulation->preprocessEdges();

    #ifndef TTK_ENABLE_GRAPHVIZ
    {
        std::stringstream msg;
        msg << "[PlanarGraphLayout] ERROR: This filter requires GraphViz" << endl;
        dMsg(std::cout, msg.str(), timeMsg);
        return 1;
    }
    #endif

    #ifdef TTK_ENABLE_GRAPHVIZ
    Timer t;

    const size_t nVertices = triangulation->getNumberOfVertices();
    const size_t nEdges = triangulation->getNumberOfEdges();

    // // =====================================================================
    // // Compute Dot String
    // // =====================================================================
    string dotString = "digraph g {rankdir=LR;";
    string nodeString = "";
    string edgeString = "";
    string rankString = "";

    // // Set Default Node style
    nodeString += "node[label=\"\",shape=box];";

    auto tl = [](size_t t){return "\"t"+to_string(t)+"\"";};
    auto nl = [](size_t id){return to_string(id);};

    for(size_t t=0; t<nTime-1; t++){
        edgeString += tl(t)+"->"+tl(t+1)+";";
    }

    vector< vector<size_t> > timeNodeMap;
    timeNodeMap.resize( nTime );

    for(size_t i=0; i<nVertices; i++){
        float x,y,z;
        triangulation->getVertexPoint(i,x,y,z);
        timeNodeMap[ (size_t) x ].push_back( i );
    }

    // complete rank string
    for(size_t t=0; t<nTime; t++){
        rankString += "{rank=same "+ tl(t);

        auto& nodes = timeNodeMap[t];
        for(auto& n: nodes){
            rankString += " "+nl(n);
        }

        rankString += "}";
    }

    // edges
    vector<vector<size_t>> nodeEdgeMap( nVertices );
    int n0;
    int n1;
    edgeString += "edge[weight=1];";
    for(size_t i=0; i<nEdges; i++){

        triangulation->getEdgeVertex(i, 0, n0);
        triangulation->getEdgeVertex(i, 1, n1);
        auto& b0 = nodeBranchMap[ n0 ];
        auto& b1 = nodeBranchMap[ n1 ];

        edgeString += nl(n0)+"->"+nl(n1)+ (b0!=b1 ? "[weight=0];":";");
    }

    string dot = dotString + nodeString + edgeString + rankString + "}";
    // cout<<dot<<endl;

    Agraph_t* G = agmemread(dot.c_str());
    GVC_t *gvc = gvContext();
    gvLayout(gvc, G, "dot");

    for(size_t i=0; i<nVertices; i++){
        Agnode_t* n = agnode(G, const_cast<char*>(nl(i).c_str()), 0);
        if(n!=nullptr){
            auto c = ND_coord(n);
            outputVertices[ i*3+1 ] = c.y/100;
        }
    }

    {
        std::stringstream msg;
        msg << "[PlanarGraphLayout] Layout Computed in "
            << t.getElapsedTime() << " s."
            << endl;
        dMsg(std::cout, msg.str(), timeMsg);
    }

    #endif

  return 0;
}
