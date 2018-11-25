/// \ingroup base
/// \class ttk::PlanarGraphLayout
/// \author Wiebke Koepp (wiebke.koepp@gmail.com) and Jonas Lukasczyk (jl@jluk.de)
/// \date 1.11.2018
///
/// \brief TTK %planarGraphLayout processing package.
///
/// %PlanarGraphLayout is a TTK processing package that TODO.
///
/// \sa ttk::Triangulation

#pragma once

#include <map>

// base code includes
#include <Wrapper.h>

#if TTK_ENABLE_GRAPHVIZ
#include <graphviz/cgraph.h>
#include <graphviz/gvc.h>
#endif

using namespace std;

namespace ttk{

    class PlanarGraphLayout : public Debug{

        public:

            PlanarGraphLayout(){};
            ~PlanarGraphLayout(){};

            // Execute the geometry approximation.
            template <class dataType> int execute(
                // Input
                void* levelsP,
                long long* topology,
                size_t nVertices,
                size_t nEdges,

                // Output
                float* layout
            ) const;
    };
}

template <class dataType> int ttk::PlanarGraphLayout::execute(
    // Input
    void* levelsP,
    long long* topology,
    size_t nVertices,
    size_t nEdges,

    // Output
    float* layout
) const{

    auto levels = (dataType*) levelsP;

    Timer t;
    double t0=0;
    // Print Input
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Computing layout for graph with" << endl
            << "[ttkPlanarGraphLayout]     "<< nVertices << " vertices" << endl
            << "[ttkPlanarGraphLayout]     "<< nEdges << " edges" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Get labels and clear layout
    map<double, size_t> levelToLevelIndex;
    {
        for(size_t i=0; i<nVertices; i++){
            levelToLevelIndex[ levels[i] ] = 0;
            layout[i] = 0;
        }
        size_t i=0;
        for(auto& t: levelToLevelIndex)
            t.second = i++;
    }
    size_t nLevels = levelToLevelIndex.size();

    // =========================================================================
    // Compute Dot String
    // =====================================================================
    string dotString = "digraph g {rankdir=LR;";
    string nodeString = "";
    string edgeString = "";
    string rankString = "";

    // Set Default Node style
    nodeString += "node[label=\"\",shape=box];";

    auto tl = [](size_t t){return "\"t"+to_string(t)+"\"";};
    auto nl = [](size_t id){return to_string(id);};

    // Add level chain
    edgeString += tl(0);
    for(size_t t=1; t<nLevels-1; t++)
        edgeString += "->"+tl(t+1);
    edgeString += ";";

    // ranks
    {
        vector< vector<size_t> > levelIndexToNodesIndex( nLevels );
        for(size_t i=0; i<nVertices; i++)
            levelIndexToNodesIndex[
                levelToLevelIndex[ levels[i] ]
            ].push_back( i );

        for(size_t t=0; t<nLevels; t++){
            rankString += "{rank=same "+ tl(t);

            auto& nodes = levelIndexToNodesIndex[t];
            for(auto& i: nodes){
                rankString += " "+nl(i);
            }

            rankString += "}";
        }
    }

    // edges
    {
        edgeString += "edge[weight=1];";

        size_t n = nEdges*3;
        for(size_t i=0; i<n; i+=3){
            long long i0 = topology[i+1];
            long long i1 = topology[i+2];
            edgeString += nl(i0)+"->"+nl(i1)+";";
        }
    }

    string dot = dotString + nodeString + edgeString + rankString + "}";
    // Print Status
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Dot String generated in " << (t.getElapsedTime()-t0) << " s." << endl;
        dMsg(cout, msg.str(), timeMsg);
        t0 = t.getElapsedTime();
    }
    dMsg(cout, dot, advancedInfoMsg);

    // =========================================================================
    // Compute Layout
    // =========================================================================
    #if TTK_ENABLE_GRAPHVIZ
        dMsg(cout, "[ttkPlanarGraphLayout] Computing layout ... ", timeMsg);
        Agraph_t* G = agmemread( dot.c_str() );
        GVC_t *gvc = gvContext();
        gvLayout(gvc, G, "dot");

        for(size_t i=0; i<nVertices; i++){
            Agnode_t* n = agnode(G, const_cast<char*>(nl(i).c_str()), 0);
            if(n!=nullptr){
                auto c = ND_coord(n);
                layout[i] = c.y/100.;
            }
        }

        {
            stringstream msg;
            msg << "done (" << (t.getElapsedTime()-t0) <<  " s)" << endl;
            dMsg(cout, msg.str(), timeMsg);
            t0 = t.getElapsedTime();
        }
    #else
        dMsg(cout, "[ttkPlanarGraphLayout] ERROR: This filter requires GraphViz to compute the layout.\n", fatalMsg);
        return 1;
    #endif

    // Print performance
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Layout computed in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 0;
}
