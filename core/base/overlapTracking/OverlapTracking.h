/// \ingroup base
/// \class ttk::OverlapTracking
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %overlapTracking processing package.
///
/// %OverlapTracking is a TTK processing package that TODO

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

typedef long long int labelType;

struct Node {
    size_t id;

    Node(size_t id=0)
        : id(id){
    }
};

struct Edge {
    size_t i;
    size_t j;
    size_t overlap;

    Edge(size_t i=0, size_t j=0, size_t overlap=0)
        : i(i), j(j), overlap(overlap) {
    }
};

struct TrackingComputationData {
    vector<size_t> sortedIndicies;
    vector<size_t> slicesI;
    vector<float>  slicesV;

    float*     pointCoords;
    labelType* pointLabels;
    size_t     nPoints;
};

namespace ttk{
    class OverlapTracking : public Debug{
        public:
            OverlapTracking();
            ~OverlapTracking();

            int reset();

            int processTimestep(
                float* pointCoords,
                labelType* pointLabels,
                size_t nPoints
            );

        private:
            // Tracking Graph
            vector< vector<Node> > timeNodesMap; // Nodes at time t
            vector< vector<Edge> > timeEdgesMap; // Edges from time t to t+1

            // Previous Timestep
            bool firstTimeStep;
            TrackingComputationData oldTCD;
    };
}