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
            template <typename idType, typename sizeType> int execute(
                // Input
                const float* inputPoints,
                const idType* inputTopology,
                const sizeType* inputPointSizes,
                const size_t nInputPoints,
                const size_t nInputCells,
                const size_t nSubdivisions,

                const size_t axis1,
                const size_t axis2,

                // Output
                float* outputPoints,
                idType* outputTopology
            ) const;
    };
}

template <typename idType, typename sizeType> int ttk::MeshGraph::execute(
    // Input
    const float* inputPoints,
    const idType* inputTopology,
    const sizeType* inputPointSizes,
    const size_t nInputPoints,
    const size_t nInputCells,
    const size_t nSubdivisions,

    const size_t axis1,
    const size_t axis2,

    // Output
    float* outputPoints,
    idType* outputTopology
) const{

    Timer t;
    double t0=0;

    // Print Input
    {
        stringstream msg;
        msg << "[ttkMeshGraph] Computing mesh for graph with" << endl
            << "[ttkMeshGraph]  - "<< nInputPoints << " points" << endl
            << "[ttkMeshGraph]  - "<< nInputCells << " edges" << endl;
        dMsg(cout, msg.str(), infoMsg);
    }

    size_t axis3 = (axis1==0 && axis2==1) || (axis1==1 && axis2==0)
        ? 2
        : (axis1==1 && axis2==2) || (axis1==2 && axis2==1)
            ? 1
            : (axis1==0 && axis2==2) || (axis1==2 && axis2==0)
                ? 0
                : 4;

    if(axis3==4){
        dMsg(cout, "[ttkMeshGraph] ERROR: Input axes have incorrect configuration.\n", fatalMsg);
        return 0;
    }

    auto getInputPointData = [](
        size_t pointIndex,
        const float* inputPoints,
        const sizeType* inputPointSizes,
        float data[4]
    ) {
        size_t i= pointIndex*3;
        data[0] = inputPoints[i++];
        data[1] = inputPoints[i++];
        data[2] = inputPoints[i];

        data[3] = (float) inputPointSizes[pointIndex];
    };

    size_t subdivisionOffset = nInputPoints*2;
    size_t nSubdivisionPoints = nSubdivisions*2;
    size_t outputPointsSubdivisonOffset = nSubdivisionPoints*3;

    // =========================================================================
    // Compute Output Point Locations
    // =========================================================================
    // outputPoints: [
    //     corners: 2*inputPoints in inputPoints order;
    //     SubPoints: 2 per subdivison in cell order]
    // ]
    {
        dMsg(cout, "[ttkMeshGraph] Computing output points ... ", timeMsg);
        t0 = t.getElapsedTime();

        // -------------------------------------------------------------------------
        // Compute Corners
        // -------------------------------------------------------------------------
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t i=0; i<nInputPoints; i++){
            float data[4];
            getInputPointData(i, inputPoints, inputPointSizes, data);

            size_t q = i*6;
            outputPoints[q  ] = data[0];
            outputPoints[q+1] = data[1];
            outputPoints[q+2] = data[2];

            outputPoints[q+3] = data[0];
            outputPoints[q+4] = data[1];
            outputPoints[q+5] = data[2];

            outputPoints[q+  axis2] += data[3]/2;
            outputPoints[q+3+axis2] -= data[3]/2;
        }

        // -------------------------------------------------------------------------
        // Compute SubPoints
        // -------------------------------------------------------------------------
        size_t q = subdivisionOffset*3;
        float nSubdivisionsP1 = nSubdivisions+1;

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t i=0; i<nInputCells; i++){
            size_t temp = i*3+1;
            size_t n0 = (size_t) inputTopology[ temp++ ];
            size_t n1 = (size_t) inputTopology[ temp   ];

            size_t no0 = n0*6;
            size_t no1 = n1*6;

            float dx0 = (outputPoints[no1+axis1] - outputPoints[no0+axis1])/nSubdivisionsP1;
            float dy0 = (outputPoints[no1+axis2] - outputPoints[no0+axis2])/nSubdivisionsP1;
            float dx1 = (outputPoints[no1+axis1+3] - outputPoints[no0+axis1+3])/nSubdivisionsP1;
            float dy1 = (outputPoints[no1+axis2+3] - outputPoints[no0+axis2+3])/nSubdivisionsP1;

            size_t q2 = q + i*outputPointsSubdivisonOffset;
            for(size_t j=1; j<=nSubdivisions; j++){
                outputPoints[q2+axis1] = outputPoints[no0+axis1] + j*dx0;
                outputPoints[q2+axis2] = outputPoints[no0+axis2] + j*dy0;
                outputPoints[q2+axis3] = outputPoints[no0+axis3];
                outputPoints[q2+axis1+3] = outputPoints[no0+axis1+3] + j*dx1;
                outputPoints[q2+axis2+3] = outputPoints[no0+axis2+3] + j*dy1;
                outputPoints[q2+axis3+3] = outputPoints[no0+axis3+3];

                q2 += 6;
            }
        }

        // Print Status
        {
            stringstream msg;
            msg << "done ("<<(t.getElapsedTime()-t0)<<" s)."<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    // =========================================================================
    // Compute Output Cells
    // =========================================================================
    //  c0 - s0u - s1u - s2u ... c3
    //  |                         |
    //  c1 - s0d - s1d - s2d ... c2
    {
        dMsg(cout, "[ttkMeshGraph] Computing output cells  ... ", timeMsg);
        t0 = t.getElapsedTime();

        size_t cellSize = this->computeOutputCellSize( nSubdivisions );
        idType cellDim = ((idType)cellSize)-1;

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(threadNumber_)
        #endif
        for(size_t i=0; i<nInputCells; i++){
            size_t q = i*3+1;
            idType in0 = inputTopology[q++]*2;
            idType in1 = inputTopology[q]*2;

            idType c0 = in0;
            idType c1 = in0+1;
            idType c2 = in1+1;
            idType c3 = in1;

            size_t q2 = cellSize*i;
            outputTopology[q2++] = cellDim;
            outputTopology[q2++] = c0;
            outputTopology[q2++] = c1;

            size_t temp = subdivisionOffset + i*nSubdivisionPoints;

            for(size_t j=0; j<nSubdivisions; j++)
                outputTopology[q2++] = (idType) (temp + j*2 + 1);

            outputTopology[q2++] = c2;
            outputTopology[q2++] = c3;

            for(int j=nSubdivisions-1; j>=0; j--)
                outputTopology[q2++] = (idType) (temp + j*2);
        }

        // Print Status
        {
            stringstream msg;
            msg << "done ("<<(t.getElapsedTime()-t0)<<" s)."<<endl;
            dMsg(cout, msg.str(), timeMsg);
        }
    }

    return 1;
}
