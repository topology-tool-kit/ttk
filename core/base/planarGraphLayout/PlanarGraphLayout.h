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

            template <typename idType, typename sequenceType> int execute(
                // Input
                const sequenceType* pointSequences,
                const float* pointSizes,
                const idType* pointBranchIds,
                const idType* levels,
                const idType* topology,
                const size_t& nPoints,
                const size_t& nEdges,

                // Output
                float* layout
            ) const;

            template <typename idType> int extractLevel(
                // Input
                const idType& level,
                const idType* levels,
                const idType* topology,
                const size_t& nPoints,
                const size_t& nEdges,

                // Output
                vector<size_t>& nodeIndicies,
                vector<size_t>& edgeIndicies
            ) const;

            template <typename idType> int getParents(
                // Input
                const idType& level,
                const idType* levels,
                const idType* topology,
                const size_t& nPoints,
                const size_t& nEdges,

                // Output
                vector<size_t>& nodeIndicies,
                vector<size_t>& edgeIndicies
            ) const;

            template <typename idType, typename sequenceType> int computeDotString(
                const sequenceType* pointSequences,
                const float* pointSizes,
                const idType* pointBranchIds,
                const map<sequenceType, size_t>& sequenceValueToIndexMap,

                const idType* topology,
                const vector<size_t>& nodeIndicies,
                const vector<size_t>& edgeIndicies,

                string& dotString
            ) const;

            // Compute Dot Layout
            int computeDotLayout(
                const vector<size_t>& nodeIndicies,
                const string& dotString,
                float* layout
            ) const {
                #if TTK_ENABLE_GRAPHVIZ
                    Timer t;

                    dMsg(cout, "[ttkPlanarGraphLayout] Computing layout ... ", timeMsg);

                    auto nl = [](size_t id){return to_string(id);};

                    Agraph_t* G = agmemread( dotString.data() );
                    GVC_t* gvc = gvContext();
                    gvLayout(gvc, G, "dot");

                    for(auto i: nodeIndicies){
                        Agnode_t* n = agnode(G, const_cast<char*>(nl(i).data()), 0);
                        if(n!=nullptr)
                            layout[i] = ND_coord(n).y;
                    }

                    gvFreeLayout(gvc, G);
                    agclose(G);
                    gvFreeContext(gvc);

                    for(auto i: nodeIndicies)
                        layout[i] = layout[i]/72; // points to inches

                    {
                        stringstream msg;
                        msg << "done (" << t.getElapsedTime() <<  " s)" << endl;
                        dMsg(cout, msg.str(), timeMsg);
                    }

                    return 1;
                #else
                    dMsg(cout, "[ttkPlanarGraphLayout] ERROR: This filter requires GraphViz to compute the layout.\n", fatalMsg);
                    return 0;
                #endif
            };
    };
}

// =============================================================================
// Extract Level
// =============================================================================
template <typename idType> int ttk::PlanarGraphLayout::extractLevel(
    // Input
    const idType& level,
    const idType* levels,
    const idType* topology,
    const size_t& nPoints,
    const size_t& nEdges,

    // Output
    vector<size_t>& nodeIndicies,
    vector<size_t>& edgeIndicies
) const {

    // If levels==nullptr then return all points and edges
    if(levels==nullptr){
        nodeIndicies.resize(nPoints);
        for(size_t i=0; i<nPoints; i++)
            nodeIndicies[i] = i;

        edgeIndicies.resize(nEdges);
        for(size_t i=0; i<nEdges; i++)
            edgeIndicies[i] = i;

        return 1;
    }

    // Get number of nodes at level
    for(size_t i=0; i<nPoints; i++)
        if(levels[i]==level)
            nodeIndicies.push_back( i );

    // Get number of edges at level
    size_t nEdges3 = nEdges*3;
    for(size_t i=0; i<nEdges3; i+=3){
        auto n0l = levels[ topology[i+1] ];
        auto n1l = levels[ topology[i+2] ];
        if( n0l==level && n0l==n1l )
            edgeIndicies.push_back( i/3 );
    }

    return 1;
}

// =============================================================================
// Get Parents
// =============================================================================
template <typename idType> int ttk::PlanarGraphLayout::getParents(
    // Input
    const idType& level,
    const idType* levels,
    const idType* topology,
    const size_t& nPoints,
    const size_t& nEdges,

    // Output
    vector<size_t>& nodeIndicies,
    vector<size_t>& edgeIndicies
) const {

    // If levels==nullptr then return all points and edges
    if(levels==nullptr){
        nodeIndicies.resize(nPoints);
        for(size_t i=0; i<nPoints; i++)
            nodeIndicies[i] = i;

        edgeIndicies.resize(nEdges);
        for(size_t i=0; i<nEdges; i++)
            edgeIndicies[i] = i;

        return 1;
    }

    // Get number of nodes at level
    for(size_t i=0; i<nPoints; i++)
        if(levels[i]==level)
            nodeIndicies.push_back( i );

    // Get number of edges at level
    size_t nEdges3 = nEdges*3;
    for(size_t i=0; i<nEdges3; i+=3){
        auto n0l = levels[ topology[i+1] ];
        auto n1l = levels[ topology[i+2] ];
        if( n0l==level && n0l==n1l )
            edgeIndicies.push_back( i/3 );
    }

    return 1;
}

// =============================================================================
// Compute Dot String
// =============================================================================
template <typename idType, typename sequenceType> int ttk::PlanarGraphLayout::computeDotString(
    const sequenceType* pointSequences,
    const float* pointSizes,
    const idType* pointBranchIds,
    const map<sequenceType, size_t>& sequenceValueToIndexMap,

    const idType* topology,
    const vector<size_t>& nodeIndicies,
    const vector<size_t>& edgeIndicies,

    string& dotString
) const {

    Timer t;

    bool usePointSequence = pointSequences!=nullptr;
    bool usePointSize = pointSizes!=nullptr;
    bool usePointBranch = pointBranchIds!=nullptr;

    string headString = "digraph g {rankdir=LR;";
    string nodeString = "";
    string edgeString = "";
    string rankString = "";

    // lambda functions that generate string representations of nodes
    auto tl = [](size_t t){return "\"t"+to_string(t)+"\"";};
    auto nl = [](size_t id){return to_string(id);};

    // Nodes
    {
        // Set default node style
        nodeString += "node[label=\"\",shape=box,width=1,height=1];";

        // If usePointSize then add size via node height
        if(usePointSize)
            for(auto& i: nodeIndicies)
                nodeString += nl(i)+"[height="+to_string(pointSizes[i])+"];";
    }

    // Ranks
    if(usePointSequence){

        size_t nSequenceValues = sequenceValueToIndexMap.size();

        // Sequence Chain
        {
            edgeString += tl(0);
            for(size_t t=1; t<nSequenceValues; t++)
                edgeString += "->"+tl(t);
            edgeString += "[weight=1];";
        }

        // Rank String
        vector< vector<size_t> > sequenceIndexToPointIndexMap( nSequenceValues );
        for(auto& i: nodeIndicies)
            sequenceIndexToPointIndexMap[
                sequenceValueToIndexMap.find(
                    pointSequences[i]
                )->second
            ].push_back( i );

        for(size_t t=0; t<nSequenceValues; t++){
            rankString += "{rank=same "+ tl(t);

            auto& nodes = sequenceIndexToPointIndexMap[t];
            for(auto& i: nodes)
                rankString += " "+nl(i);

            rankString += "}";
        }
    }

    // Edges
    {
        for(auto& edgeIndex: edgeIndicies){
            size_t temp = edgeIndex*3;
            auto& i0 = topology[ temp+1 ];
            auto& i1 = topology[ temp+2 ];
            edgeString += nl(i0)+"->"+nl(i1);

            if(usePointBranch){
                auto b0 = pointBranchIds[i0];
                auto b1 = pointBranchIds[i1];
                edgeString += b0==b1
                    ? "[weight=1]"
                    : "[weight=0]";
            }

            edgeString += ";";
        }
    }

    // Build Dot String
    {
        dotString = headString + nodeString + edgeString + rankString + "}";
    }

    // Print Status
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Dot String generated in " << t.getElapsedTime() << " s." << endl;
        dMsg(cout, msg.str(), timeMsg);

        dMsg(cout, "\n"+dotString+"\n", advancedInfoMsg);
    }

    return 1;
}

// =============================================================================
// Execute
// =============================================================================
template <typename idType, typename sequenceType> int ttk::PlanarGraphLayout::execute(
    // Input
    const sequenceType* pointSequences,
    const float* pointSizes,
    const idType* pointBranchIds,
    const idType* levels,
    const idType* topology,
    const size_t& nPoints,
    const size_t& nEdges,

    // Output
    float* layout
) const{

    Timer t;

    // Init Input
    bool usePointSequence = pointSequences!=nullptr;
    bool usePointSize = pointSizes!=nullptr;
    bool usePointBranch = pointBranchIds!=nullptr;
    bool useLevels = levels!=nullptr;

    // Print Input
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Computing layout for graph with" << endl
            << "[ttkPlanarGraphLayout]  - "<< nPoints << " vertices" << endl
            << "[ttkPlanarGraphLayout]  - "<< nEdges << " edges" << endl;
        if(usePointSequence) msg << "[ttkPlanarGraphLayout]  - using sequences" << endl;
        if(usePointSize)     msg << "[ttkPlanarGraphLayout]  - using sizes" << endl;
        if(usePointBranch)   msg << "[ttkPlanarGraphLayout]  - using branches" << endl;
        if(useLevels)        msg << "[ttkPlanarGraphLayout]  - using levels" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    if(useLevels && !usePointSize){
        dMsg(cout, "[ttkPlanarGraphLayout] ERROR: When 'UseLevels' is enabled then 'UseSizes' must also be enabled.\n", fatalMsg);
        return 0;
    }

    // Global SequenceValue to SequenceIndex map
    map<sequenceType, size_t> sequenceValueToIndexMap;
    if(usePointSequence){
        for(size_t i=0; i<nPoints; i++){
            sequenceValueToIndexMap[ pointSequences[i] ] = 0;
            layout[i] = 0;
        }
        size_t i=0;
        for(auto& t: sequenceValueToIndexMap)
            t.second = i++;
    }

    idType nLevels = 0;
    if(useLevels){
        for(size_t i=0; i<nPoints; i++)
            if(nLevels<levels[i]) nLevels=levels[i];
    }

    // Compute initial layout for each level
    for(idType l=0; l<=nLevels; l++){
        vector<size_t> nodeIndicies;
        vector<size_t> edgeIndicies;

        // Extract nodes and edges at certain level
        {
            int status = this->extractLevel<idType>(
                // Input
                l,
                levels,
                topology,
                nPoints,
                nEdges,

                // Output
                nodeIndicies,
                edgeIndicies
            );
            if(status!=1) return 0;
        }

        // Compute Dot String
        string dotString;
        {
            int status = this->computeDotString<idType, sequenceType>(
                // Input
                pointSequences,
                pointSizes,
                pointBranchIds,
                sequenceValueToIndexMap,

                topology,
                nodeIndicies,
                edgeIndicies,

                // Output
                dotString
            );
            if(status!=1) return 0;
        }

        // Compute Dot Layout
        {
            int status = this->computeDotLayout(
                nodeIndicies,
                dotString,
                layout
            );
            if(status!=1) return 0;
        }
    }

    // If nLevels>0 then compute slots
    if(nLevels>0){

        struct ChildrenComparator {
            const float* layout;

            ChildrenComparator(const float* layout) : layout(layout){};

            inline bool operator() (const size_t& i, const size_t& j) {
                return layout[i]<layout[j];
            }
        };

        auto comparator = ChildrenComparator(layout);

        vector<vector<size_t>> nodeIndexChildrenIndexMap( nPoints );

        size_t nEdges3 = nEdges*3;
        for(size_t i=0; i<nEdges3; i+=3){
            auto n0 = topology[ i+1 ];
            auto n1 = topology[ i+2 ];
            if( (levels[n0]+1) == levels[n1] )
                nodeIndexChildrenIndexMap[ n0 ].push_back(n1);
        }

        // Adjust positions from bottom to top
        for(idType l=0; l<nLevels; l++){
            vector<size_t> nodeIndicies;
            vector<size_t> edgeIndicies;

            // get nodes at current level (parents)
            this->extractLevel<idType>(
                // Input
                l,
                levels,
                topology,
                nPoints,
                nEdges,

                // Output
                nodeIndicies,
                edgeIndicies
            );

            // for each parent adjust position of children
            for(auto& parent: nodeIndicies){
                auto& children = nodeIndexChildrenIndexMap[ parent ];
                size_t nChildren = children.size();
                if(nChildren<1) continue;

                // sort children
                sort(children.begin(), children.end(), comparator);

                // size of parent
                float sizeParent = pointSizes[parent];

                // size of child
                float sizeChildren = 0;
                for(auto& child: children)
                    sizeChildren += pointSizes[child];


                // gap space
                float gap = sizeParent - sizeChildren;
                float gapDelta = (gap/(nChildren+1))/2;

                float y = layout[ parent ] + sizeParent*0.5 - gapDelta;
                for(auto& child: children){
                    float temp = gapDelta + pointSizes[child]/2;
                    layout[child] = y - temp;
                    y -= 2*temp;
                }
            }
        }
    }

    // Print performance
    {
        stringstream msg;
        msg << "[ttkPlanarGraphLayout] Layout computed in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 1;
}
