/// \ingroup base
/// \class ttk::MeshGraph
/// \author Wiebke Koepp (wiebke.koepp@gmail.com) and Jonas Lukasczyk (jl@jluk.de)
/// \date 1.11.2018
///
/// \brief TTK %meshGraph processing package.
///
/// %MeshGraph is a TTK processing package that TODO.
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Wrapper.h>

using namespace std;

namespace ttk{

    class MeshGraph : public Debug{

        public:

            MeshGraph(){};
            ~MeshGraph(){};

            inline size_t computeNumberOfOutputPoints(size_t nInputPoints, size_t nInputCells, size_t nSubdivisions) const{
                return nInputPoints*2 + nInputCells*(nSubdivisions*2); // one vertex becomes two + edge subdivisons
            };
            inline size_t computeOutputCellSize(size_t nSubdivisions) const{
                return 5+nSubdivisions*2; // cellDim + 4 corners + 2 for each subdivision
            };
            inline size_t computeOutputTopologySize(size_t nInputCells, size_t nSubdivisions) const{
                return nInputCells*this->computeOutputCellSize( nSubdivisions );
            };

            // Execute the geometry approximation.
            template <class dataType> int execute(
                // Input
                const float* inputPoints,
                const dataType* inputTopology,
                const double* inputPointSizes,
                size_t nInputPoints,
                size_t nInputCells,
                size_t nSubdivisions,

                // Output
                float* outputPoints,
                dataType* outputTopology
            ) const;
    };
}

template <class dataType> int ttk::MeshGraph::execute(
    // Input
    const float* inputPoints,
    const dataType* inputTopology,
    const double* inputPointSizes,
    size_t nInputPoints,
    size_t nInputCells,
    size_t nSubdivisions,

    // Output
    float* outputPoints,
    dataType* outputTopology
) const{

    Timer t;
    double t0;

    // Print Input
    {
        stringstream msg;
        msg << "[ttkMeshGraph] Computing mesh for graph with" << endl
            << "[ttkMeshGraph]  - "<< nInputPoints << " points" << endl
            << "[ttkMeshGraph]  - "<< nInputCells << " edges" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    auto getInputPointData = [](
        size_t pointIndex,
        const float* inputPoints,
        const double* inputPointSizes,
        float data[4]
    ) {
        size_t i= pointIndex*3;
        data[0] = inputPoints[i++];
        data[1] = inputPoints[i++];
        data[2] = inputPoints[i];

        data[3] = (float) inputPointSizes[pointIndex];
    };

    size_t nOutputPoints = this->computeNumberOfOutputPoints(nInputPoints, nInputCells, nSubdivisions); // one vertex becomes two + subdivisons on edges
    // outputPoints: [2*inputPoints in inputPoints order; subdivisionPoints in cell order]

    // compute corners
    {
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t i=0; i<nInputPoints; i++){
            float data[4];
            getInputPointData(i, inputPoints, inputPointSizes,  data);

            size_t q = i*6;
            outputPoints[q++] = data[0];
            outputPoints[q++] = data[1]+data[3]/2;
            outputPoints[q++] = data[2];

            outputPoints[q++] = data[0];
            outputPoints[q++] = data[1]-data[3]/2;
            outputPoints[q++] = data[2];
        }
    }

    // compute subdivisions
    size_t subdivisionOffset = nInputPoints*2;
    size_t nSubdivisionPoints = nSubdivisions*2;
    size_t outputPointsSubdivisonOffset = nSubdivisionPoints*3;
    {
        dataType q = subdivisionOffset*3;
        float nSubdivisionsP1 = nSubdivisions+1;

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t i=0; i<nInputCells; i++){
            size_t temp = i*3+1;
            dataType n0 = inputTopology[ temp++ ];
            dataType n1 = inputTopology[ temp   ];

            size_t no0 = n0*6;
            size_t no1 = n1*6;

            float dx0 = (outputPoints[no1] - outputPoints[no0])/nSubdivisionsP1;
            float dy0 = (outputPoints[no1+1] - outputPoints[no0+1])/nSubdivisionsP1;
            float dx1 = (outputPoints[no1+3] - outputPoints[no0+3])/nSubdivisionsP1;
            float dy1 = (outputPoints[no1+4] - outputPoints[no0+4])/nSubdivisionsP1;

            size_t q2 = q + i*outputPointsSubdivisonOffset;
            for(size_t j=1; j<=nSubdivisions; j++){
                outputPoints[q2++] = outputPoints[no0  ] + j*dx0;
                outputPoints[q2++] = outputPoints[no0+1] + j*dy0;
                outputPoints[q2++] = outputPoints[no0+2];
                outputPoints[q2++] = outputPoints[no0+3] + j*dx1;
                outputPoints[q2++] = outputPoints[no0+4] + j*dy1;
                outputPoints[q2++] = outputPoints[no0+5];
            }
        }
    }

    // Generate Cells
    /*
        c0 - s0u - s1u - s2u ... c3
        |                         |
        c1 - s0d - s1d - s2d ... c2
    */
    size_t cellSize = this->computeOutputCellSize( nSubdivisions );
    size_t cellDim = cellSize-1;

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber_)
    #endif
    for(size_t i=0; i<nInputCells; i++){
        size_t q = i*3+1;
        dataType in0 = inputTopology[q++]*2;
        dataType in1 = inputTopology[q]*2;

        dataType c0 = in0;
        dataType c1 = in0+1;
        dataType c2 = in1+1;
        dataType c3 = in1;

        size_t q2 = cellSize*i;
        outputTopology[q2++] = cellDim;
        outputTopology[q2++] = c0;
        outputTopology[q2++] = c1;

        for(size_t j=0; j<nSubdivisions; j++)
            outputTopology[q2++] = subdivisionOffset + i*nSubdivisionPoints + j*2 + 1;

        outputTopology[q2++] = c2;
        outputTopology[q2++] = c3;

        for(int j=nSubdivisions-1; j>=0; j--)
            outputTopology[q2++] = subdivisionOffset + i*nSubdivisionPoints + j*2;
    }


    // =========================================================================
    // Comment Level 0
    // =========================================================================

    // -------------------------------------------------------------------------
    // Comment Level 1
    // -------------------------------------------------------------------------

    // Comment Level 2

    // Print performance
    {
        stringstream msg;
        msg << "[ttkMeshGraph] Mesh computed in " << t.getElapsedTime() << " s. (" << threadNumber_ << " thread(s))." << endl;
        dMsg(cout, msg.str(), timeMsg);
    }

    return 0;
}
