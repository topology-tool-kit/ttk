#include <OverlapTracking.h>

#include <algorithm>
#include <unordered_set>
#include <unordered_map>

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

ttk::OverlapTracking::OverlapTracking(){
    this->reset();
}
ttk::OverlapTracking::~OverlapTracking(){}

int ttk::OverlapTracking::reset(){
    cout<<"Reset"<<endl;
    this->timeNodesMap.clear();
    this->timeEdgesMap.clear();
    this->firstTimeStep = true;

    return 0;
}

int ttk::OverlapTracking::processTimestep(
    float* pointCoords,
    labelType*   pointLabels,
    size_t nPoints
){
    cout<<"Process "<<nPoints<<endl;

    // Initialize nodes
    unordered_map<size_t, size_t> labelIndexMap;
    vector<Node> nodes;
    {
        unordered_set<size_t> labels;
        for(size_t i=0; i<nPoints; i++)
            labels.insert( (size_t) pointLabels[i] );
        nodes.resize( labels.size() );
        size_t i=0;
        for(auto& x: labels){
            nodes[i].id = x;
            labelIndexMap[x] = i;
            i++;
        }
    }
    cout<<labelIndexMap.size()<<endl;
    cout<<nodes.size()<<endl;
    cout<<nodes[0].id<<endl;
    cout<<"-------------"<<endl;

    // Sort indicies for grid comparison
    TrackingComputationData newTCD;
    {
        newTCD.sortedIndicies.resize( nPoints );
        for(size_t i=0; i<nPoints; i++)
            newTCD.sortedIndicies[i] = i;

        CoordinateComparator c = CoordinateComparator(pointCoords);
        sort(newTCD.sortedIndicies.begin(), newTCD.sortedIndicies.end(), c);

        float oldX = pointCoords[ newTCD.sortedIndicies[0]*3 ];
        for(size_t i=0; i<nPoints; i++){
            float currentX = pointCoords[ newTCD.sortedIndicies[i]*3 ];
            if(oldX!=currentX){
                newTCD.slicesI.push_back( i );
                newTCD.slicesV.push_back( oldX );
                oldX=currentX;
            }
        }
    }
    newTCD.pointCoords = pointCoords;
    newTCD.pointLabels = pointLabels;
    newTCD.nPoints = nPoints;

    // Tracking
    auto compare = [](float* pointCoordsI, float* pointCoordsJ, size_t i, size_t j) {
        size_t i3 = i*3;
        size_t j3 = j*3;

        float iX = pointCoordsI[i3++];
        float iY = pointCoordsI[i3++];
        float iZ = pointCoordsI[i3];

        float jX = pointCoordsJ[j3++];
        float jY = pointCoordsJ[j3++];
        float jZ = pointCoordsJ[j3];

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

    size_t count = 0;
    for(size_t i=0; i<nPoints; i++)
        if(pointLabels[i]==0) count++;

    cout<<" -> "<<count<<endl;

    if(!this->firstTimeStep){
        TrackingComputationData& oldTCD = this->oldTCD;

        // unordered_map<size_t, Edge> newIndexToEdgeMap;

        // for(size_t i=0; i<newTCD.slicesI.size(); i++){
        //     float newX = newTCD.slicesV[i];

        //     // Scan for x in old
        //     size_t j=0;
        //     bool found = false;
        //     for(j=0; j<oldTCD.slicesV.size(); j++){
        //         if( oldTCD.slicesV[j] == newX ){
        //             found = true;
        //             break;
        //         }
        //     }
        //     if(!found) cout<<"not found "<<newX<<endl;

        //     if(found){
        //         cout<<"Track Slice"<<endl;
                // size_t k=newTCD.slicesI[i];
                // size_t l=oldTCD.slicesI[j];
                // size_t kLimit = i==newTCD.slicesI.size()-1 ? newTCD.nPoints : newTCD.slicesI[i+1];
                // size_t lLimit = j==oldTCD.slicesI.size()-1 ? oldTCD.nPoints : oldTCD.slicesI[j+1];

                size_t i = 0;
                size_t j = 0;

                unordered_map<labelType, unordered_map<labelType, size_t>> edges;

                while(i<newTCD.nPoints && j<oldTCD.nPoints){
                    int c = compare(
                        newTCD.pointCoords,
                        oldTCD.pointCoords,
                        newTCD.sortedIndicies[i],
                        oldTCD.sortedIndicies[j]
                    );
                    if(c == 0){ // Same
                        size_t index0 = newTCD.sortedIndicies[i];
                        size_t index1 = oldTCD.sortedIndicies[j];

                        auto label0 = newTCD.pointLabels[ index0 ];
                        auto label1 = oldTCD.pointLabels[ index1 ];

                        auto edgesBack = edges.find(label0);
                        if(edgesBack == edges.end()){
                            edges.insert({ label0, unordered_map<labelType, size_t>() });
                            edgesBack = edges.find(label0);
                        }

                        auto edge = edgesBack->second.find(label1);
                        if(edge == edgesBack->second.end()){
                            edgesBack->second[label1] = 0;
                            edgesBack->second.insert({ label1, 0 });
                            edge = edgesBack->second.find(label1);
                        }

                        edge->second++;

                        // cout
                        //     <<"o " << newTCD.pointLabels[ index0 ]<<" - " << oldTCD.pointLabels[ index1 ] << " ["
                        //     << newTCD.pointCoords[ index0*3 ]
                        //     << " "
                        //     << newTCD.pointCoords[ index0*3+1 ]
                        //     << " "
                        //     << newTCD.pointCoords[ index0*3+2 ]
                        //     << "] ["
                        //     << oldTCD.pointCoords[ index1*3 ]
                        //     << " "
                        //     << oldTCD.pointCoords[ index1*3+1 ]
                        //     << " "
                        //     << oldTCD.pointCoords[ index1*3+2 ]
                        //     << "]" << endl;
                        i++;
                        j++;
                    } else if (c>0){ // new larger than old
                        j++;
                    } else { // old larger than new
                        i++;
                    }
                }

                // for(auto& i: edges){
                //     cout<<i.first<<": "<<endl;
                //     for(auto& j: i.second)
                //         cout<<"    "<<j.first<<" -> "<<j.second<<endl;
                // }

                // 5497
                // 3413

                // 6223
                // 3884
            // }

        //     // float oldX = oldTCD.slicesV[j];
        //     // cout<<newX<< " " <<oldX<<endl;
        //     // while( oldX<newX){
        //     //     j++;
        //     //     oldX = oldTCD.pointCoords[ oldTCD.sortedIndicies[j] ];
        //     // }
        //     // if(oldX == newX){
        //     //     cout<< "Track "<<i<<" "<<j<<endl;
        //     //     j++;
        //     // }
        // }
    }

    // Finalizing
    this->firstTimeStep = false;
    this->oldTCD = newTCD;

    this->timeNodesMap.push_back( nodes );
    // this->timeEdgesMap.push_back( edges );


    return 0;
}
