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
                void* sequencesP,
                float* sizes,
                unsigned long long* branches,
                unsigned int* levels,
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
    void* sequencesP,
    float* sizes,
    unsigned long long* branches,
    unsigned int* levels,
    long long* topology,
    size_t nVertices,
    size_t nEdges,

    // Output
    float* layout
) const{

    // Init input
    auto sequences = (dataType*) sequencesP;
    bool useSizes = sizes!=nullptr;
    bool useBranches = branches!=nullptr;
    bool useLevels = levels!=nullptr;

    Timer t;
    double t0=0;

    // Print Input
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Computing layout for graph with" << endl
            << "[ttkPlanarGraphLayout]  - "<< nVertices << " vertices" << endl
            << "[ttkPlanarGraphLayout]  - "<< nEdges << " edges" << endl;
        if(useSizes)    msg << "[ttkPlanarGraphLayout]  - using sizes" << endl;
        if(useBranches) msg << "[ttkPlanarGraphLayout]  - using branches" << endl;
        if(useLevels)   msg << "[ttkPlanarGraphLayout]  - using levels" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    // Get labels and clear layout
    map<double, size_t> sequenceToSequenceIndex;
    {
        for(size_t i=0; i<nVertices; i++){
            sequenceToSequenceIndex[ sequences[i] ] = 0;
            layout[i] = 0;
        }
        size_t i=0;
        for(auto& t: sequenceToSequenceIndex)
            t.second = i++;
    }
    size_t nSequences = sequenceToSequenceIndex.size();

    // =========================================================================
    // Compute Dot String
    // =====================================================================
    string dotString = "digraph g {rankdir=LR;";
    string nodeString = "";
    string edgeString = "";
    string rankString = "";

    auto tl = [](size_t t){return "\"t"+to_string(t)+"\"";};
    auto nl = [](size_t id){return to_string(id);};

    // Set default node style
    nodeString += "node[label=\"\",shape=box,width=1,height=1];";

    if(useSizes)
        for(size_t i=0; i<nVertices; i++)
            nodeString += nl(i)+"[height="+to_string(sizes[i])+"];";

    // Add sequence chain
    edgeString += tl(0);
    for(size_t t=1; t<nSequences-1; t++)
        edgeString += "->"+tl(t+1);
    edgeString += "[weight=1];";

    // ranks
    {
        vector< vector<size_t> > sequenceIndexToNodesIndex( nSequences );
        for(size_t i=0; i<nVertices; i++)
            sequenceIndexToNodesIndex[
                sequenceToSequenceIndex[ sequences[i] ]
            ].push_back( i );

        for(size_t t=0; t<nSequences; t++){
            rankString += "{rank=same "+ tl(t);

            auto& nodes = sequenceIndexToNodesIndex[t];
            for(auto& i: nodes)
                rankString += " "+nl(i);

            rankString += "}";
        }
    }

    // edges
    {
        size_t n = nEdges*3;
        for(size_t i=0; i<n; i+=3){
            long long i0 = topology[i+1];
            long long i1 = topology[i+2];
            edgeString += nl(i0)+"->"+nl(i1);

            if(useBranches){
                auto b0 = branches[i0];
                auto b1 = branches[i1];
                edgeString += b0==b1
                    ? "[weight=1]"
                    : "[weight=0]";
            }

            edgeString += ";";
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
    dMsg(cout, dot+"\n", advancedInfoMsg);

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
