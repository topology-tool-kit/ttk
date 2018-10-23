/// \ingroup base
/// \class ttk::OverlapTracking
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK %overlapTracking processing package.
///
/// %OverlapTracking is a TTK processing package that TODO

#pragma once

#include <unordered_map>

// base code includes
#include <Wrapper.h>

using namespace std;

typedef long long labelType;

struct Node {
    labelType label;

    Node(labelType label=0)
        : label(label){
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
    vector<float> pointCoords;
    vector<size_t> pointLabelIndicies;
    size_t nPoints;

    TrackingComputationData(float* pointCoords, labelType* pointLabels, size_t nPoints, unordered_map<labelType, size_t>& labelIndexMap){
        this->pointCoords.resize(nPoints*3);
        this->pointLabelIndicies.resize(nPoints);
        this->nPoints = nPoints;
        for(size_t i=0; i<nPoints; i++){
            size_t j = i*3;
            this->pointCoords[j] = pointCoords[j];
            this->pointCoords[j+1] = pointCoords[j+1];
            this->pointCoords[j+2] = pointCoords[j+2];
            this->pointLabelIndicies[i] = labelIndexMap[ pointLabels[i] ];
        }
    }
};

namespace ttk{
    class OverlapTracking : public Debug{
        public:
            OverlapTracking();
            ~OverlapTracking();

            int reset();

            vector< vector<Node> >& getTimeNodesMap();
            vector< vector<Edge> >& getTimeEdgesMap();

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
            TrackingComputationData* prevTCD;
    };
}