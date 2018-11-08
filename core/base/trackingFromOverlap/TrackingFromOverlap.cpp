#include <TrackingFromOverlap.h>

#include <algorithm>
#include <unordered_set>

using namespace ttk;
using namespace std;

struct CoordinateComparator {
    float* coords;

    CoordinateComparator(float* coords) : coords(coords){};

    inline bool operator() (const size_t& i, const size_t& j) {
        size_t ic = i*3;
        size_t jc = j*3;
        return coords[ic]==coords[jc]
            ? coords[ic+1]==coords[jc+1]
                ? coords[ic+2]<coords[jc+2]
                : coords[ic+1]<coords[jc+1]
            : coords[ic]<coords[jc];
    }
};

ttk::TrackingFromOverlap::TrackingFromOverlap(){
    this->prevTCD = nullptr;
    this->reset();
}
ttk::TrackingFromOverlap::~TrackingFromOverlap(){}

int ttk::TrackingFromOverlap::reset(){
    this->timeNodeLabelMap = vector< vector<labelType> >();
    this->timeEdgesMap = vector< vector<Edge> >();

    if(this->prevTCD!=nullptr)
        delete this->prevTCD;
    this->prevTCD = nullptr;

    // Print Status
    {
        stringstream msg;
        msg << "[ttkTrackingFromOverlap] Initialized "<<endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 0;
}

vector< vector<labelType> >& ttk::TrackingFromOverlap::getTimeNodeLabelMap(){
    return this->timeNodeLabelMap;
};

vector< vector<Edge> >& ttk::TrackingFromOverlap::getTimeEdgesMap(){
    return this->timeEdgesMap;
};

int ttk::TrackingFromOverlap::processTimestep(
    float*     pointCoords,
    labelType* pointLabels,
    size_t     nPoints
){
    Timer t;
    auto tD = t.getElapsedTime();

    // Print Status
    {
        stringstream msg;
        msg << "[ttkTrackingFromOverlap] Creating nodes ... "<<flush;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Initialize nodes
    this->timeNodeLabelMap.resize( this->timeNodeLabelMap.size()+1 ); // Add a node set
    unordered_map<labelType, size_t> labelIndexMap;
    vector<labelType>& nodeLabels = this->timeNodeLabelMap[ this->timeNodeLabelMap.size()-1 ];
    {
        unordered_set<labelType> labels;
        for(size_t i=0; i<nPoints; i++)
            labels.insert( pointLabels[i] );
        nodeLabels.resize( labels.size() );
        size_t i=0;
        for(auto& x: labels){
            nodeLabels[i] = x;
            labelIndexMap[x] = i;
            i++;
        }
    }

    // Print Status
    {
        auto tTemp = t.getElapsedTime();
        stringstream msg;
        msg << "done ("<<(tTemp-tD)<<" s)."<<endl
            << "[ttkTrackingFromOverlap] Sorting  nodes ... "<<flush;
        tD = tTemp;
        dMsg(cout, msg.str(), timeMsg);
    }

    // Sort indicies for grid comparison
    TrackingComputationData* currTCD = new TrackingComputationData();
    {
        currTCD->nPoints = nPoints;

        currTCD->pointLabelIndicies.resize( nPoints );
        for(size_t i=0; i<nPoints; i++)
            currTCD->pointLabelIndicies[i] = labelIndexMap[ pointLabels[i] ];

        currTCD->pointCoords.resize( nPoints*3 );
        for(size_t i=0; i<nPoints*3; i++)
            currTCD->pointCoords[i] = pointCoords[i];

        currTCD->sortedIndicies.resize( nPoints );
        for(size_t i=0; i<nPoints; i++)
            currTCD->sortedIndicies[i] = i;

        CoordinateComparator c = CoordinateComparator(pointCoords);
        sort(currTCD->sortedIndicies.begin(), currTCD->sortedIndicies.end(), c);
    }

    // Tracking
    auto compare = [](vector<float>& prevPointCoords, vector<float>& currPointCoords, size_t i, size_t j) {
        size_t i3 = i*3;
        size_t j3 = j*3;

        float iX = prevPointCoords[i3++];
        float iY = prevPointCoords[i3++];
        float iZ = prevPointCoords[i3];

        float jX = currPointCoords[j3++];
        float jY = currPointCoords[j3++];
        float jZ = currPointCoords[j3];

        return iX==jX
            ? iY==jY
                ? iZ==jZ
                    ? 0
                    : iZ<jZ
                        ? -1 : 1
                : iY<jY
                    ? -1 : 1
            : iX<jX
                ? -1 : 1;
    };

    // Print Status
    {
        auto tTemp = t.getElapsedTime();
        stringstream msg;
        msg << "done ("<<(tTemp-tD)<<" s)."<<endl;
        dMsg(cout, msg.str(), timeMsg);
        tD=tTemp;
    }

    if(this->prevTCD != nullptr){
        // Print Status
        {
            stringstream msg;
            msg << "[ttkTrackingFromOverlap] Tracking       ... "<<flush;
            dMsg(cout, msg.str(), timeMsg);
        }
        TrackingComputationData* prevTCD = this->prevTCD;

        this->timeEdgesMap.resize( this->timeEdgesMap.size()+1 ); // Add an edge set
        auto& edges = this->timeEdgesMap[ this->timeEdgesMap.size()-1 ];

        unordered_map<size_t, unordered_map<size_t, size_t>> edgesMap0;

        {
            size_t i = 0; // iterator for prev timestep
            size_t j = 0; // iterator for curr timestep

            while(i<prevTCD->nPoints && j<currTCD->nPoints){
                int c = compare(
                    prevTCD->pointCoords,
                    currTCD->pointCoords,
                    prevTCD->sortedIndicies[i],
                    currTCD->sortedIndicies[j]
                );
                if(c == 0){ // Same
                    size_t index0 = prevTCD->sortedIndicies[i];
                    size_t index1 = currTCD->sortedIndicies[j];

                    size_t labelIndex0 = prevTCD->pointLabelIndicies[ index0 ];
                    size_t labelIndex1 = currTCD->pointLabelIndicies[ index1 ];

                    // Find edge and increase overlap counter
                    auto edges0 = edgesMap0.find(labelIndex0); // Edges from node label0 from prev timestep to nodes of current timestep
                    if(edges0 == edgesMap0.end()){
                        edgesMap0[labelIndex0] = unordered_map<size_t, size_t>();
                        edges0 = edgesMap0.find(labelIndex0);
                    }
                    auto edge = edges0->second.find(labelIndex1);
                    if(edge == edges0->second.end()){
                        edges0->second[labelIndex1] = 0;
                        edge = edges0->second.find(labelIndex1);
                    }
                    edge->second++;

                    i++;
                    j++;
                } else if (c>0){ // prev larger than curr
                    j++;
                } else { // curr larger than prev
                    i++;
                }
            }
        }

        size_t nEdges = 0;
        for(auto& i: edgesMap0)
            nEdges += i.second.size();

        edges.resize(nEdges);
        size_t q=0;
        for(auto& i: edgesMap0)
            for(auto& j: i.second){
                auto& edge = edges[q++];
                edge.i = i.first;
                edge.j = j.first;
                edge.overlap = j.second;
            }

        // Print Status
        {
            auto tTemp = t.getElapsedTime();
            stringstream msg;
            msg << "done ("<<(tTemp-tD)<<" s)."<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }


    // Finalizing
    if(this->prevTCD != nullptr)
        delete this->prevTCD;
    this->prevTCD = currTCD;

    return 0;
}
