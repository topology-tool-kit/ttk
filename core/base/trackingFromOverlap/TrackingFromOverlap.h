/// \ingroup base
/// \class ttk::TrackingFromOverlap
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %trackingFromOverlap processing package.
///
/// %TrackingFromOverlap is a TTK processing package that TODO

#pragma once

#include <unordered_map>

// base code includes
#include <Wrapper.h>

using namespace std;

typedef long long labelType;

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
    vector<float> pointCoords;
    vector<size_t> pointLabelIndicies;
    size_t nPoints;
};

namespace ttk{
    class TrackingFromOverlap : public Debug{
        public:
            TrackingFromOverlap();
            ~TrackingFromOverlap();

            int reset();

            vector< vector<labelType> >& getTimeNodeLabelMap();
            vector< vector<Edge> >& getTimeEdgesMap();

            int processTimestep(
                float* pointCoords,
                labelType* pointLabels,
                size_t nPoints
            );

        private:
            // Tracking Graph
            vector< vector<labelType> > timeNodeLabelMap; // Nodes at time t
            vector< vector<Edge> > timeEdgesMap; // Edges from time t to t+1

            // Previous Timestep
            TrackingComputationData* prevTCD;
    };
}